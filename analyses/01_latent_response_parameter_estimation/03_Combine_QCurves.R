# ============================================================================
# 03_Combine_QCurves.R
# ----------------------------------------------------------------------------
# WHAT THIS SCRIPT DOES
#   Takes the fitted curves from 01 (West Africa) and 02 (DRC conflict /
#   conflict++), puts every scenario on a common 0..730 day grid, builds the two
#   hypothetical West-Africa-with-conflict scenarios from them, and writes one
#   combined CSV in the same schema as the original analysis outputs.
#
# THE COMMON GRID
#   Every scenario is reported on integer days 0..730 (731 rows) with
#   tau = relative_day / 730. Each fitted curve is placed on its REAL outbreak-day
#   axis and held flat (last value repeated) beyond its own support out to day 730
#   -- it is NOT stretched to fill the horizon. (West Africa data reach ~day 357,
#   the DRC conflict curve ~day 457; both plateau from there to 730.)
#
# FIVE scenarios are produced (the DRC no-conflict scenario is intentionally held
# out pending the Model A vs Model B decision; see DRC_no_conflict_checking.R):
#   1. worst_west_africa                     (Model A fit, script 01)
#   2. drc_conflict                          (Model B fit, script 02)
#   3. drc_conflict_plusplus                 (Model B fit, script 02)
#   4. worst_west_africa_conflict            (constructed here from 1 + 2)
#   5. worst_west_africa_conflict_plusplus   (constructed here from 4)
#
# Two construction steps are folded in here (rather than as post-hoc patches):
#   * the West-Africa-with-conflict hybrid Q on real days (scenario 4), and
#   * the conflict++ collapse (scenario 5): force parameters to poor-response
#     endpoints over days 200-300 -- the same window where scenario 4's conflict
#     dip falls -- with q_value forced to 0 there.
# Note: under the ORIGINAL methodology ipc_helper IS the fitted latent_IPC; no
# separate q-scaling of IPC is applied (that is a revised-methodology step).
#
# Inputs : data-processed/WestAfrica_QCurve/WestAfrica_QCurve_Fit.rds,
#          DRC_QCurve/DRC_QCurve_Conflict_Fit.rds, DRC_QCurve/DRC_QCurve_ConflictPlusPlus_Fit.rds
# Output : data-processed/combined_original_methodology_outputs.csv
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(ggplot2)   # display-only grid plot of the combined scenarios at the end
})

source(here::here("analyses", "01_latent_response_parameter_estimation", "helpers.R"))

# Common reporting grid: integer days 0..730, with tau = day / 730.
days    <- 0:HORIZON_DAYS
tau_out <- days / HORIZON_DAYS

# ----------------------------------------------------------------------------
# Generic assembly helpers
# ----------------------------------------------------------------------------

# Interpolate a fitted curve (columns: parameter, relative_day, mean) onto the
# 0..730 grid ON ITS REAL DAY AXIS, holding the last value flat beyond the
# curve's own support out to day 730 (make_interp uses rule = 2). NO stretching:
# a fit that only reaches day 357 plateaus from 357 to 730.
# Returns long format: parameter, relative_day, value.
curve_to_daygrid_long <- function(curve_summ) {
  curve_summ %>%
    group_by(parameter) %>%
    group_modify(function(g, ...) {
      f <- make_interp(g$relative_day, g$mean)                  # interpolant on real outbreak days
      tibble(relative_day = days, value = as.numeric(f(days)))  # eval 0..730; hold flat past support
    }) %>%
    ungroup()
}

# Interpolate a q_value curve (columns: relative_day, q_value) onto the 0..730
# grid on its real day axis, holding flat past its support.
qvalue_to_daygrid <- function(q_curve) {
  f <- make_interp(q_curve$relative_day, q_curve$q_value)
  clip01(f(days))
}

# Turn per-parameter day-grid values (long: parameter, relative_day, value) plus
# a per-day q_value vector into ONE scenario block in the final output schema.
# This is where the internal parameter names (p_hosp, p_ETU, latent_IPC, ...) are
# renamed to the published column names (prob_hosp, prop_etu, ipc_helper, ...).
assemble_scenario <- function(scenario_key, scenario_name, param_long, q_value) {
  wide <- param_long %>%
    pivot_wider(names_from = parameter, values_from = value) %>%   # one column per parameter
    arrange(relative_day)                                          # rows in day order 0..730

  tibble(
    methodology              = "original",
    scenario_key             = scenario_key,
    scenario                 = scenario_name,
    time_index               = days + 1L,        # 1..731
    relative_day             = days,             # 0..730
    tau                      = tau_out,          # day / 730
    prob_hosp                = wide$p_hosp,
    delay_hosp               = wide$delay_hosp,
    prob_unsafe_funeral_comm = wide$p_unsafe_funeral_comm,
    prob_unsafe_funeral_hosp = wide$p_unsafe_funeral_hosp,
    prob_unsafe_funeral_etu  = 0,                       # always 0 in this analysis
    prop_etu                 = wide$p_ETU,
    ipc_helper               = wide$latent_IPC,         # original methodology: fitted IPC
    q_value                  = q_value
  )
}

# Diagnostic q_value from a set of per-parameter Q curves: the row-wise mean
# across parameters, rescaled to [0,1]. (q_value is a summary index only; the
# parameter columns above are what actually drive the downstream model.)
qvalue_from_param_Q <- function(param_Q_wide) {
  rescale_01(rowMeans(select(param_Q_wide, all_of(PARAM_LEVELS)), na.rm = TRUE))
}

# ----------------------------------------------------------------------------
# Load the fitted objects
# ----------------------------------------------------------------------------
wa_fit                    <- readRDS(file.path(DIR_PROCESSED, "WestAfrica_QCurve/WestAfrica_QCurve_Fit.rds"))
drc_conflict_fit          <- readRDS(file.path(DIR_PROCESSED, "DRC_QCurve/DRC_QCurve_Conflict_Fit.rds"))
drc_conflict_plusplus_fit <- readRDS(file.path(DIR_PROCESSED, "DRC_QCurve/DRC_QCurve_ConflictPlusPlus_Fit.rds"))

# ----------------------------------------------------------------------------
# Scenario 1: worst_west_africa
# ----------------------------------------------------------------------------
# Straight read-out of the Model A fit. q_value is the diagnostic mean of the six
# estimated per-parameter Q_j curves (each runs 0 -> 1), rescaled to [0,1].
wa_q_value_curve <- wa_fit$q_summ %>%
  select(parameter, relative_day, mean) %>%
  pivot_wider(names_from = parameter, values_from = mean) %>%
  mutate(q_value = qvalue_from_param_Q(.)) %>%
  select(relative_day, q_value)

scen_wa <- assemble_scenario(
  "worst_west_africa", "Worst_WestAfrica",
  curve_to_daygrid_long(wa_fit$curve_summ),     # the fitted parameter curves
  qvalue_to_daygrid(wa_q_value_curve)
)

# ----------------------------------------------------------------------------
# Scenarios 2 & 3: drc_conflict and drc_conflict_plusplus
# ----------------------------------------------------------------------------
# Straight read-outs of the two Model B fits. For DRC the scenario's q_value is
# the shared empirical SDB curve itself (not a mean of per-parameter curves).
scen_drc_conflict <- assemble_scenario(
  "middle_drc_conflict", "Middle_DRC_ConflictSmoothed",
  curve_to_daygrid_long(drc_conflict_fit$curve_summ),
  qvalue_to_daygrid(drc_conflict_fit$q_grid)
)

scen_drc_conflict_pp <- assemble_scenario(
  "middle_drc_conflict_plusplus", "Middle_DRC_ConflictSmoothed_PlusPlus",
  curve_to_daygrid_long(drc_conflict_plusplus_fit$curve_summ),
  qvalue_to_daygrid(drc_conflict_plusplus_fit$q_grid)
)

# ----------------------------------------------------------------------------
# Scenario 4: worst_west_africa_conflict  (constructed)
# ----------------------------------------------------------------------------
# Idea: keep the West Africa response MAGNITUDES, but disrupt its TIMING with the
# DRC conflict. For each parameter j we build a hybrid response-quality curve on
# REAL outbreak days (no stretching):
#
#     Q_hybrid_j(day) = Q_WA_j(day)  x  Q_DRC_conflict(day)
#
# Both factors are evaluated on the real day axis and held flat beyond their own
# support out to day 730 (West Africa reaches ~357, the DRC conflict curve ~457).
# So the DRC conflict dip lands at its REAL time (~days 200-300) -- which means
# the conflict in this scenario coincides with the conflict++ collapse window
# applied below. The hybrid Q is then mapped onto the WA fitted start/end magnitudes.

# DRC conflict shared Q on the real day axis (held flat past its ~457-day support).
drcQ_f <- make_interp(drc_conflict_fit$q_grid$relative_day, drc_conflict_fit$q_grid$q_value)

# West Africa fitted start/end magnitudes per parameter (the endpoints the hybrid
# Q interpolates between): theta_start at the earliest day, theta_end at the latest.
theta_bounds <- wa_fit$curve_summ %>%
  group_by(parameter) %>%
  summarise(theta_start = mean[which.min(relative_day)],
            theta_end   = mean[which.max(relative_day)], .groups = "drop")

# Build the per-parameter hybrid Q and the resulting native-unit trajectory.
hybrid <- lapply(PARAM_LEVELS, function(p) {
  wq    <- filter(wa_fit$q_summ, parameter == p)
  waQ_f <- make_interp(wq$relative_day, wq$mean)     # WA Q_j on the real day axis

  # Pointwise product of the WA response and the DRC conflict modulator, both on
  # real days (each held flat past its own support).
  q_raw   <- clip01(waQ_f(days)) * clip01(drcQ_f(days))
  q_scale <- max(q_raw, na.rm = TRUE)
  if (!is.finite(q_scale) || q_scale <= 0) stop("Hybrid Q non-positive for ", p)
  q_hybrid <- clip01(q_raw / q_scale)        # per-parameter max-scaling (NO min-subtraction,
                                             # so late dips are preserved, not floored to 0)
  q_hybrid[days == 0] <- 0                    # force the response to start at 0 on day 0

  tibble(parameter = p, relative_day = days, q_hybrid = q_hybrid)
}) %>% bind_rows()

# Map the hybrid Q (0..1) back onto each parameter's WA fitted magnitude range:
#   value = theta_start + (theta_end - theta_start) * Q_hybrid.
wa_conflict_long <- hybrid %>%
  left_join(theta_bounds, by = "parameter") %>%
  mutate(value = theta_start + (theta_end - theta_start) * q_hybrid) %>%
  select(parameter, relative_day, value)

# Diagnostic q_value: mean of the per-parameter hybrid Q, rescaled to [0,1].
wa_conflict_qvalue <- hybrid %>%
  select(parameter, relative_day, q_hybrid) %>%
  pivot_wider(names_from = parameter, values_from = q_hybrid) %>%
  { qvalue_from_param_Q(.) }

scen_wa_conflict <- assemble_scenario(
  "worst_west_africa_conflict", "Worst_WestAfrica_Conflict",
  wa_conflict_long, wa_conflict_qvalue
)

# ----------------------------------------------------------------------------
# Scenario 5: worst_west_africa_conflict_plusplus  (constructed from scenario 4)
# ----------------------------------------------------------------------------
# The "++" scenario adds a hard temporary response collapse over days 200-300:
# every patched parameter is forced to its WORST value seen on the WA-conflict
# trajectory (minimum for good-response params, maximum for adverse params), and
# q_value is forced to 0. prob_unsafe_funeral_etu is left untouched (it is 0).
# (This is the corrected "deterioration, not improvement" collapse direction.)
PLUSPLUS_WINDOW <- c(200, 300)
good_response_params <- c("prob_hosp", "prop_etu", "ipc_helper")             # worst = minimum
adverse_params       <- c("delay_hosp", "prob_unsafe_funeral_comm", "prob_unsafe_funeral_hosp")  # worst = maximum

# Poor-response endpoint per parameter, taken from the extrema of the WA-conflict
# baseline trajectory (min for good-response, max for adverse).
poor_endpoint <- c(
  sapply(good_response_params, function(p) min(scen_wa_conflict[[p]], na.rm = TRUE)),
  sapply(adverse_params,       function(p) max(scen_wa_conflict[[p]], na.rm = TRUE))
)

# Start from a copy of scenario 4, relabel it, then overwrite the collapse window.
scen_wa_conflict_pp <- scen_wa_conflict %>%
  mutate(
    scenario_key = "worst_west_africa_conflict_plusplus",
    scenario     = "Worst_WestAfrica_Conflict_PlusPlus"
  )

# Inside days 200-300: force each patched parameter to its poor endpoint, and
# force q_value to 0. Outside the window, the scenario-4 trajectory is kept.
in_window <- scen_wa_conflict_pp$relative_day >= PLUSPLUS_WINDOW[1] &
             scen_wa_conflict_pp$relative_day <= PLUSPLUS_WINDOW[2]
for (p in names(poor_endpoint)) {
  scen_wa_conflict_pp[[p]][in_window] <- poor_endpoint[[p]]
}
scen_wa_conflict_pp$q_value[in_window] <- 0

# ----------------------------------------------------------------------------
# Combine and write
# ----------------------------------------------------------------------------
combined <- bind_rows(
  scen_wa,
  scen_drc_conflict,
  scen_drc_conflict_pp,
  scen_wa_conflict,
  scen_wa_conflict_pp
)

# Sanity checks before writing: every scenario must have exactly 731 daily rows
# spanning 0..730, and every probability-type column must stay within [0,1].
stopifnot(all(table(combined$scenario_key) == (HORIZON_DAYS + 1L)))
prob_cols <- c("prob_hosp", "prob_unsafe_funeral_comm", "prob_unsafe_funeral_hosp",
               "prob_unsafe_funeral_etu", "prop_etu", "ipc_helper", "q_value")
for (col in prob_cols) {
  v <- combined[[col]]
  if (any(v < -1e-8 | v > 1 + 1e-8, na.rm = TRUE)) {
    stop("Column ", col, " has values outside [0,1].")
  }
}

write_csv(combined, file.path(DIR_PROCESSED, "combined_original_methodology_outputs.csv"))
message("\n03_Combine_QCurves.R complete. Wrote combined_original_methodology_outputs.csv ",
        "to data-processed/ (", nrow(combined), " rows, ",
        length(unique(combined$scenario_key)), " scenarios).")

# ----------------------------------------------------------------------------
# Plot the combined output: parameter (row) x scenario (column)   (display only)
# ----------------------------------------------------------------------------
# A quick visual check of everything just written to the CSV. Each cell is one
# parameter's trajectory over the 0..730 day horizon for one scenario, laid out
# as parameter rows x scenario columns. The bottom row is the diagnostic q_value
# index. NOTE: these are the native output columns of `combined` (probabilities,
# the delay in days, the IPC index) -- the per-parameter normalised Q_j curves
# are not stored in the combined matrix, only the single q_value index is. The
# plot is printed to the graphics device and deliberately NOT saved.

# Scenario columns and parameter rows in a sensible reading order.
scenario_order <- c("worst_west_africa", "middle_drc_conflict", "middle_drc_conflict_plusplus",
                    "worst_west_africa_conflict", "worst_west_africa_conflict_plusplus")
scenario_labels <- c(
  worst_west_africa                    = "WA",
  middle_drc_conflict                  = "DRC conflict",
  middle_drc_conflict_plusplus         = "DRC conflict++",
  worst_west_africa_conflict           = "WA conflict",
  worst_west_africa_conflict_plusplus  = "WA conflict++"
)
quantity_order <- c("prob_hosp", "delay_hosp", "prob_unsafe_funeral_comm",
                    "prob_unsafe_funeral_hosp", "prob_unsafe_funeral_etu",
                    "prop_etu", "ipc_helper", "q_value")

# Reshape the wide combined matrix to one row per (scenario, parameter, day).
combined_long <- combined %>%
  select(scenario_key, relative_day, all_of(quantity_order)) %>%
  pivot_longer(cols = all_of(quantity_order), names_to = "quantity", values_to = "value") %>%
  mutate(scenario_key = factor(scenario_key, levels = scenario_order),
         quantity     = factor(quantity,     levels = quantity_order))

p_grid <- ggplot(combined_long, aes(relative_day, value)) +
  geom_line(linewidth = 0.6, colour = "#1f77b4") +
  # free y per ROW (shared within a row across scenarios, so a parameter is
  # comparable left-to-right; delay_hosp is in days, the rest are in [0,1]).
  facet_grid(quantity ~ scenario_key, scales = "free_y", switch = "y",
             labeller = labeller(scenario_key = scenario_labels)) +
  labs(title = "Combined scenario curves: parameter (row) x scenario (column)",
       x = "Relative outbreak day", y = NULL) +
  theme_bw(base_size = 9) +
  theme(strip.text.y.left = element_text(angle = 0),
        strip.placement = "outside",
        panel.grid.minor = element_blank())

print(p_grid)   # display only; not saved
