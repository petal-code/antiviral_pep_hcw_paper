# ============================================================================
# 03_Combine_QCurves.R
# ----------------------------------------------------------------------------
# Assemble the final scenario matrices on the common 0..730 day grid and write a
# single combined CSV in the same schema as the original outputs.
#
# Inputs (from data-processed/):
#   wa_fit.rds, drc_conflict_fit.rds, drc_conflict_plusplus_fit.rds, drc_durations.csv
# Output:
#   data-processed/combined_original_methodology_outputs.csv
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
#   * the West-Africa-with-conflict hybrid Q + time stretch, and
#   * the conflict++ collapse (force parameters to poor-response endpoints over
#     days 200-300), with q_value forced to 0 there.
# Note: under the ORIGINAL methodology ipc_helper IS the fitted latent_IPC; no
# separate q-scaling of IPC is applied (that is a revised-methodology step).
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
})

source("helpers.R")

# Common reporting grid: integer days 0..730, with tau = day / 730.
days    <- 0:HORIZON_DAYS
tau_out <- days / HORIZON_DAYS

# ----------------------------------------------------------------------------
# Generic assembly helpers
# ----------------------------------------------------------------------------

# Interpolate a fitted curve (columns: parameter, tau, mean) onto the day grid.
# Returns long format: parameter, relative_day, value.
curve_to_daygrid_long <- function(curve_summ) {
  curve_summ %>%
    group_by(parameter) %>%
    group_modify(function(g, ...) {
      f <- make_interp(g$tau, g$mean)
      tibble(relative_day = days, value = as.numeric(f(tau_out)))
    }) %>%
    ungroup()
}

# Interpolate a q_value curve (columns: tau, q_value) onto the day grid.
qvalue_to_daygrid <- function(q_curve) {
  f <- make_interp(q_curve$tau, q_curve$q_value)
  clip01(f(tau_out))
}

# Turn per-parameter day-grid values (long: parameter, relative_day, value) plus
# a per-day q_value vector into one scenario block in the final output schema.
# The internal parameter names are renamed to the published column names here.
assemble_scenario <- function(scenario_key, scenario_name, param_long, q_value) {
  wide <- param_long %>%
    pivot_wider(names_from = parameter, values_from = value) %>%
    arrange(relative_day)

  tibble(
    methodology              = "original",
    scenario_key             = scenario_key,
    scenario                 = scenario_name,
    time_index               = days + 1L,
    relative_day             = days,
    tau                      = tau_out,
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
# across parameters, rescaled to [0,1].
qvalue_from_param_Q <- function(param_Q_wide) {
  rescale_01(rowMeans(select(param_Q_wide, all_of(PARAM_LEVELS)), na.rm = TRUE))
}

# ----------------------------------------------------------------------------
# Load the fitted objects
# ----------------------------------------------------------------------------
wa_fit                    <- readRDS("data-processed/wa_fit.rds")
drc_conflict_fit          <- readRDS("data-processed/drc_conflict_fit.rds")
drc_conflict_plusplus_fit <- readRDS("data-processed/drc_conflict_plusplus_fit.rds")
drc_durations             <- read_csv("data-processed/drc_durations.csv", show_col_types = FALSE)

# ----------------------------------------------------------------------------
# Scenario 1: worst_west_africa
# ----------------------------------------------------------------------------
# q_value is the diagnostic mean of the six estimated per-parameter Q_j curves.
wa_q_value_curve <- wa_fit$q_summ %>%
  select(parameter, tau, mean) %>%
  pivot_wider(names_from = parameter, values_from = mean) %>%
  mutate(q_value = qvalue_from_param_Q(.)) %>%
  select(tau, q_value)

scen_wa <- assemble_scenario(
  "worst_west_africa", "Worst_WestAfrica",
  curve_to_daygrid_long(wa_fit$curve_summ),
  qvalue_to_daygrid(wa_q_value_curve)
)

# ----------------------------------------------------------------------------
# Scenarios 2 & 3: drc_conflict and drc_conflict_plusplus
# ----------------------------------------------------------------------------
# For the DRC scenarios the q_value is the shared empirical Q curve itself.
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
# Idea: keep the West Africa response MAGNITUDES, but reshape the TIMING with the
# DRC conflict disruption and stretch the timeline. Concretely, for each
# parameter j we build a hybrid response-quality curve
#
#     Q_hybrid_j(t) = Q_WA_j(WA clock)  x  Q_DRC_conflict(extended clock)
#
# where the West Africa response matures on its ORIGINAL clock (then plateaus)
# while the DRC conflict modulator runs across the full, stretched window. The
# hybrid Q is then mapped back onto the West Africa fitted start/end magnitudes.

# Time-stretch multiplier = duration(DRC conflict) / duration(DRC no-conflict).
dur <- setNames(drc_durations$max_day, drc_durations$scenario)
duration_multiplier <- dur[["drc_conflict"]] / dur[["drc_no_conflict"]]
if (!is.finite(duration_multiplier) || duration_multiplier < 1) {
  warning("Invalid duration multiplier; falling back to 1.50.")
  duration_multiplier <- 1.50
}
message("West Africa conflict time-stretch multiplier = ", round(duration_multiplier, 3))

# Two clocks on the reporting grid:
tau_wa_progress <- pmin(tau_out * duration_multiplier, 1)  # WA matures then plateaus
tau_drc         <- tau_out                                 # conflict runs full window

# Interpolants: WA per-parameter Q_j, and the DRC conflict shared Q.
drcQ_f <- make_interp(drc_conflict_fit$q_grid$tau, drc_conflict_fit$q_grid$q_value)

# West Africa fitted start/end magnitudes per parameter (the endpoints the hybrid
# Q interpolates between).
theta_bounds <- wa_fit$curve_summ %>%
  group_by(parameter) %>%
  summarise(theta_start = mean[which.min(tau)],
            theta_end   = mean[which.max(tau)], .groups = "drop")

# Build the per-parameter hybrid Q and the resulting native-unit trajectory.
hybrid <- lapply(PARAM_LEVELS, function(p) {
  wq  <- filter(wa_fit$q_summ, parameter == p)
  waQ_f <- make_interp(wq$tau, wq$mean)

  q_raw   <- clip01(waQ_f(tau_wa_progress)) * clip01(drcQ_f(tau_drc))
  q_scale <- max(q_raw, na.rm = TRUE)
  if (!is.finite(q_scale) || q_scale <= 0) stop("Hybrid Q non-positive for ", p)
  q_hybrid <- clip01(q_raw / q_scale)        # per-parameter max-scaling (no min-subtraction)
  q_hybrid[days == 0] <- 0                    # force the response to start at 0

  tibble(parameter = p, relative_day = days, q_hybrid = q_hybrid)
}) %>% bind_rows()

# Native magnitudes: theta_start + (theta_end - theta_start) * Q_hybrid.
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
PLUSPLUS_WINDOW <- c(200, 300)
good_response_params <- c("prob_hosp", "prop_etu", "ipc_helper")
adverse_params       <- c("delay_hosp", "prob_unsafe_funeral_comm", "prob_unsafe_funeral_hosp")

# Poor-response endpoints from extrema of the WA-conflict baseline trajectory.
poor_endpoint <- c(
  sapply(good_response_params, function(p) min(scen_wa_conflict[[p]], na.rm = TRUE)),
  sapply(adverse_params,       function(p) max(scen_wa_conflict[[p]], na.rm = TRUE))
)

scen_wa_conflict_pp <- scen_wa_conflict %>%
  mutate(
    scenario_key = "worst_west_africa_conflict_plusplus",
    scenario     = "Worst_WestAfrica_Conflict_PlusPlus"
  )

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

# Sanity checks: every scenario must have exactly 731 daily rows spanning 0..730,
# and probability-type columns must stay within [0,1].
stopifnot(all(table(combined$scenario_key) == (HORIZON_DAYS + 1L)))
prob_cols <- c("prob_hosp", "prob_unsafe_funeral_comm", "prob_unsafe_funeral_hosp",
               "prob_unsafe_funeral_etu", "prop_etu", "ipc_helper", "q_value")
for (col in prob_cols) {
  v <- combined[[col]]
  if (any(v < -1e-8 | v > 1 + 1e-8, na.rm = TRUE)) {
    stop("Column ", col, " has values outside [0,1].")
  }
}

write_csv(combined, "data-processed/combined_original_methodology_outputs.csv")
message("\n03_Combine_QCurves.R complete. Wrote ",
        "data-processed/combined_original_methodology_outputs.csv (",
        nrow(combined), " rows, ", length(unique(combined$scenario_key)), " scenarios).")
