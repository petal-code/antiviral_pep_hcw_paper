# ============================================================================
# 03_Combine_QCurves(revisedMethod).R
# ----------------------------------------------------------------------------
# WHAT THIS SCRIPT DOES
#   Takes the fitted/derived curves from 01 (West Africa, Model C) and 02 (DRC
#   conflict / conflict++, deterministic endpoint mapping), puts every scenario
#   on a common 0..730 day grid, and writes one combined CSV in the same schema
#   as the original analysis outputs.
#
# THE COMMON GRID
#   Every scenario is reported on integer days 0..730 (731 rows) with
#   tau = relative_day / 730. Each curve is placed on its REAL outbreak-day axis
#   and held flat (last value repeated) beyond its own support out to day 730 --
#   it is NOT stretched to fill the horizon. (West Africa data reach ~day 357, the
#   DRC conflict curve ~day 457; both plateau from there to 730.)
#
# THREE scenarios are produced (this mirrors the current original-methodology
# pipeline minus the tweak-dependent scenario: the revised methodology has no
# "tweaks", so there is no with/without-tweaks pair, and the hypothetical
# West-Africa-with-conflict scenarios are likewise not built here):
#   1. worst_west_africa            (Model C endpoint-constrained fit, script 01)
#   2. drc_conflict                 (deterministic endpoint mapping, script 02)
#   3. drc_conflict_plusplus        (deterministic endpoint mapping, script 02)
#
# WHERE THE REVISED METHODOLOGY SHOWS UP HERE
#   Nowhere as special-case code -- and that is the point. The two revised-method
#   changes are already baked into the upstream fits:
#     * endpoints are locked to literature extrema in 01/02 (Model C); and
#     * the DRC IPC/PPE index (ipc_helper = fitted latent_IPC) is the q-scaled
#       endpoint mapping produced in 02, so it is simply read through here, just
#       like every other parameter. (Under the ORIGINAL methodology ipc_helper is
#       the Model B fitted latent_IPC and there is no q-scaling; that is the one
#       substantive difference, and it lives in 02, not here.)
#
# Inputs : data-processed/WestAfrica_QCurve_revisedMethod/WestAfrica_QCurve_Fit(revisedMethod).rds,
#          DRC_QCurve_revisedMethod/DRC_QCurve_Conflict_Fit(revisedMethod).rds,
#          DRC_QCurve_revisedMethod/DRC_QCurve_ConflictPlusPlus_Fit(revisedMethod).rds
# Output : data-processed/combined_revised_methodology_outputs.csv
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(ggplot2)   # display-only grid plot of the combined scenarios at the end
})

source(here::here("analyses", "01_latent_response_parameter_estimation_revisedMethodology",
                  "helpers(revisedMethod).R"))

# Common reporting grid: integer days 0..730, with tau = day / 730.
days    <- 0:HORIZON_DAYS
tau_out <- days / HORIZON_DAYS

# ----------------------------------------------------------------------------
# Generic assembly helpers (identical in spirit to the original-methodology 03)
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
    methodology              = "revised",
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
    ipc_helper               = wide$latent_IPC,         # revised: q-scaled (DRC) or Model C fit (WA)
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
wa_fit                    <- readRDS(file.path(DIR_PROCESSED,
  "WestAfrica_QCurve_revisedMethod/WestAfrica_QCurve_Fit(revisedMethod).rds"))
drc_conflict_fit          <- readRDS(file.path(DIR_PROCESSED,
  "DRC_QCurve_revisedMethod/DRC_QCurve_Conflict_Fit(revisedMethod).rds"))
drc_conflict_plusplus_fit <- readRDS(file.path(DIR_PROCESSED,
  "DRC_QCurve_revisedMethod/DRC_QCurve_ConflictPlusPlus_Fit(revisedMethod).rds"))

# ----------------------------------------------------------------------------
# Scenario 1: worst_west_africa
# ----------------------------------------------------------------------------
# Straight read-out of the Model C fit. q_value is the diagnostic mean of the six
# estimated per-parameter Q_j curves (each runs 0 -> 1), rescaled to [0,1].
# NOTE: we pivot on tau (the shared 0..1 prediction grid), not relative_day,
# because in the revised methodology each parameter can sit on its own day axis
# (a terminal-zero parameter reaches Q = 1 before the scenario end). tau is
# identical across parameters, so the per-parameter Q_j align cleanly; the day
# axis is then recovered as tau * max_day.
wa_q_value_curve <- wa_fit$q_summ %>%
  select(parameter, tau, mean) %>%
  pivot_wider(names_from = parameter, values_from = mean) %>%
  mutate(q_value      = qvalue_from_param_Q(.),
         relative_day = tau * wa_fit$max_day) %>%
  select(relative_day, q_value)

scen_wa <- assemble_scenario(
  "worst_west_africa", "Worst_WestAfrica",
  curve_to_daygrid_long(wa_fit$curve_summ),     # the fitted parameter curves
  qvalue_to_daygrid(wa_q_value_curve)
)

# ----------------------------------------------------------------------------
# Scenarios 2 & 3: drc_conflict and drc_conflict_plusplus
# ----------------------------------------------------------------------------
# Straight read-outs of the two deterministic DRC mappings. For DRC the
# scenario's q_value is the shared empirical SDB curve itself (not a mean of
# per-parameter curves).
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
# Combine and write
# ----------------------------------------------------------------------------
combined <- bind_rows(
  scen_wa,
  scen_drc_conflict,
  scen_drc_conflict_pp
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

write_csv(combined, file.path(DIR_PROCESSED, "combined_revised_methodology_outputs.csv"))
message("\n03_Combine_QCurves(revisedMethod).R complete. Wrote combined_revised_methodology_outputs.csv ",
        "to data-processed/ (", nrow(combined), " rows, ",
        length(unique(combined$scenario_key)), " scenarios).")

# ----------------------------------------------------------------------------
# Plot the combined output: parameter (row) x scenario (column)   (display only)
# ----------------------------------------------------------------------------
# A quick visual check of everything just written to the CSV. Each cell is one
# parameter's trajectory over the 0..730 day horizon for one scenario, laid out
# as parameter rows x scenario columns. The bottom row is the diagnostic q_value
# index. The plot is printed to the graphics device and deliberately NOT saved.

# Scenario columns and parameter rows in a sensible reading order.
scenario_order <- c("worst_west_africa", "middle_drc_conflict", "middle_drc_conflict_plusplus")
scenario_labels <- c(
  worst_west_africa            = "WA",
  middle_drc_conflict          = "DRC conflict",
  middle_drc_conflict_plusplus = "DRC conflict++"
)
quantity_order <- c("prob_hosp", "delay_hosp", "prob_unsafe_funeral_comm",
                    "prob_unsafe_funeral_hosp", "prob_unsafe_funeral_etu",
                    "prop_etu", "ipc_helper", "q_value")

# Reshape the wide combined matrix to one row per (scenario, parameter, day).
combined_long <- combined %>%
  select(scenario_key, relative_day, all_of(quantity_order)) %>%
  pivot_longer(cols = all_of(quantity_order), names_to = "quantity", values_to = "value") %>%
  mutate(scenario_key = factor(scenario_key, levels = scenario_order),
         quantity     = factor(quantity,     levels = quantity_order)) %>%
  filter(quantity != "prob_unsafe_funeral_etu",
         quantity != "prob_unsafe_funeral_hosp")

p_grid <- ggplot(combined_long, aes(relative_day, value)) +
  geom_line(linewidth = 0.6, colour = "#1f77b4") +
  # free y per ROW (shared within a row across scenarios, so a parameter is
  # comparable left-to-right; delay_hosp is in days, the rest are in [0,1]).
  facet_grid(quantity ~ scenario_key, scales = "free_y", switch = "y",
             labeller = labeller(scenario_key = scenario_labels)) +
  labs(title = "Combined scenario curves (revised methodology): parameter (row) x scenario (column)",
       x = "Relative outbreak day", y = NULL) +
  theme_bw(base_size = 9) +
  theme(strip.text.y.left = element_text(angle = 0),
        strip.placement = "outside",
        panel.grid.minor = element_blank())

print(p_grid)   # display only; not saved
