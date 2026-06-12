# ============================================================================
# 03_Combine_QCurves.R
# ----------------------------------------------------------------------------
# WHAT THIS SCRIPT DOES
#   Takes the fitted curves from 01 (West Africa, with and without tweaks) and 02
#   (DRC conflict / conflict++), puts every scenario on a common 0..730 day grid,
#   and writes one combined CSV in the same schema as the original analysis outputs.
#
# THE COMMON GRID
#   Every scenario is reported on integer days 0..730 (731 rows) with
#   tau = relative_day / 730. Each fitted curve is placed on its REAL outbreak-day
#   axis and held flat (last value repeated) beyond its own support out to day 730
#   -- it is NOT stretched to fill the horizon. (West Africa data reach ~day 357,
#   the DRC conflict curve ~day 457; both plateau from there to 730.)
#
# FOUR scenarios are produced (the West-Africa-with-conflict hypotheticals have
# been removed for now; the DRC no-conflict scenario is also intentionally held out
# pending the Model A vs Model B decision; see DRC_no_conflict_checking.R):
#   1. worst_west_africa                     (Model A fit WITH tweaks, script 01)
#   2. drc_conflict                          (Model B fit, script 02)
#   3. drc_conflict_plusplus                 (Model B fit, script 02)
#   4. worst_west_africa_notweaks            (Model A fit WITHOUT tweaks: the clean
#                                             no-tweak baseline, read straight in)
#
# Note: under the ORIGINAL methodology ipc_helper IS the fitted latent_IPC; no
# separate q-scaling of IPC is applied (that is a revised-methodology step).
#
# Inputs : data-processed/WestAfrica_QCurve/WestAfrica_QCurve_Fit.rds,
#          WestAfrica_QCurve/WestAfrica_QCurve_Fit_NoTweaks.rds,
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
wa_notweaks_fit           <- readRDS(file.path(DIR_PROCESSED, "WestAfrica_QCurve/WestAfrica_QCurve_Fit_NoTweaks.rds"))
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
# Scenario 4: worst_west_africa_notweaks
# ----------------------------------------------------------------------------
# The SAME West Africa Model A fit as scenario 1, but with the targeted tweak
# priors switched OFF -- the clean no-tweak baseline (see west_africa_checking.R
# for the with-vs-without overlay). Read straight in exactly like scenario 1, so
# the two sit side by side in the grid plot and the effect of the tweaks is visible.
wa_nt_q_value_curve <- wa_notweaks_fit$q_summ %>%
  select(parameter, relative_day, mean) %>%
  pivot_wider(names_from = parameter, values_from = mean) %>%
  mutate(q_value = qvalue_from_param_Q(.)) %>%
  select(relative_day, q_value)

scen_wa_notweaks <- assemble_scenario(
  "worst_west_africa_notweaks", "Worst_WestAfrica_NoTweaks",
  curve_to_daygrid_long(wa_notweaks_fit$curve_summ),     # the no-tweak fitted parameter curves
  qvalue_to_daygrid(wa_nt_q_value_curve)
)

# ----------------------------------------------------------------------------
# Combine and write
# ----------------------------------------------------------------------------
combined <- bind_rows(
  scen_wa,
  scen_wa_notweaks,
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
scenario_order <- c("worst_west_africa", "worst_west_africa_notweaks",
                    "middle_drc_conflict", "middle_drc_conflict_plusplus")
scenario_labels <- c(
  worst_west_africa                    = "WA",
  worst_west_africa_notweaks           = "WA (no tweaks)",
  middle_drc_conflict                  = "DRC conflict",
  middle_drc_conflict_plusplus         = "DRC conflict++"
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

p_grid <- ggplot(combined_long, aes(relative_day, value, group = scenario_key)) +
  geom_line(linewidth = 0.6, colour = "#1f77b4") +
  # free y per ROW (shared within a row across scenarios, so a parameter is
  # comparable left-to-right; delay_hosp is in days, the rest are in [0,1]).
  facet_grid(quantity ~ ., scales = "free_y", switch = "y",
             labeller = labeller(scenario_key = scenario_labels)) +
  labs(title = "Combined scenario curves: parameter (row) x scenario (column)",
       x = "Relative outbreak day", y = NULL) +
  theme_bw(base_size = 9) +
  theme(strip.text.y.left = element_text(angle = 0),
        strip.placement = "outside",
        panel.grid.minor = element_blank())

print(p_grid)   # display only; not saved

library(ggh4x)

p_grid +
  facetted_pos_scales(
    y = list(
      quantity == "delay_hosp" ~ scale_y_continuous(limits = c(0, 6)),
      quantity != "prob_hosp" ~ scale_y_continuous(limits = c(0, 1))
    )
  )

# ----------------------------------------------------------------------------
# Fitted curves vs data: West Africa (A) and DRC conflict (B)   (display only)
# ----------------------------------------------------------------------------
# Side-by-side fit-vs-data check for the two headline scenarios, restricted to
# five parameters. Each region is its own 2-column x 3-row panel of the fitted
# curve (blue posterior mean + 90% interval) with the data the model saw on top:
# literature anchors (orange) and, in the DRC community unsafe-funeral panel, the
# SDB community proxy points (grey) that actually drive that curve. This mirrors
# the per-scenario figures in 01/02, placed together as (A) West Africa and
# (B) DRC conflict. Text is kept to the axes only. Display only; not saved.
library(patchwork)

# The five parameters to show, in the requested panel order. The fit objects and
# anchors use the INTERNAL names; the published-output columns are noted alongside
# (p_hosp = prob_hosp, p_ETU = prop_etu, latent_IPC = ipc_helper, ...).
plot_params  <- c("p_hosp", "delay_hosp", "p_unsafe_funeral_comm", "p_ETU", "latent_IPC")
panel_levels <- unname(PANEL_LOOKUP[plot_params])

# The fit objects store only the fitted curves, not the data they were fit to, so
# re-read the prepped inputs purely for the observed points.
wa_prep  <- readRDS(file.path(DIR_PROCESSED, "WestAfrica_QCurve/WestAfrica_QCurve_PreppedData.rds"))
drc_prep <- readRDS(file.path(DIR_PROCESSED, "DRC_QCurve/DRC_QCurve_PreppedData.rds"))

# Restrict a per-parameter table to the five plotted parameters and turn
# `parameter` into an ordered panel factor (so facets read in the requested order).
to_panels <- function(df) {
  df %>%
    filter(parameter %in% plot_params) %>%
    mutate(panel = factor(PANEL_LOOKUP[parameter], levels = panel_levels))
}

# West Africa: scenario-1 fit (with tweaks) + literature anchors.
wa_curve_df  <- to_panels(wa_fit$curve_summ)
wa_anchor_df <- to_panels(wa_prep$anchors)

# DRC conflict: fit + anchors, plus the SDB community proxy behind the community
# unsafe-funeral panel. That panel's mean was deterministically overridden in 02,
# so its q5/q95 are stale -> collapse the ribbon there to the mean (as 02 does).
drc_curve_df <- to_panels(drc_conflict_fit$curve_summ) %>%
  mutate(q5  = if_else(parameter == "p_unsafe_funeral_comm", mean, q5),
         q95 = if_else(parameter == "p_unsafe_funeral_comm", mean, q95))
drc_anchor_df <- to_panels(drc_prep$anchors)
drc_sdb_df <- drc_prep$conflict_qseries %>%
  filter(n_eligible_sum > 0) %>%
  transmute(relative_day, value_used = unsafe_funeral_comm_proxy,
            panel = factor(PANEL_LOOKUP[["p_unsafe_funeral_comm"]], levels = panel_levels))

# One region's 2x3 fit-vs-data panel; no title/subtitle (axes text only).
fit_vs_data_panel <- function(curve_df, anchor_df, extra_points = NULL) {
  p <- ggplot(curve_df, aes(relative_day, mean)) +
    geom_ribbon(aes(ymin = q5, ymax = q95), fill = "#1f77b4", alpha = 0.20) +
    geom_line(colour = "#1f77b4", linewidth = 0.9)
  if (!is.null(extra_points)) {
    p <- p + geom_point(data = extra_points, aes(relative_day, value_used),
                        inherit.aes = FALSE, colour = "grey55", size = 1, alpha = 0.7)
  }
  p +
    geom_point(data = anchor_df, aes(relative_day, value_used),
               inherit.aes = FALSE, colour = "#ff7f0e", size = 2) +
    facet_wrap(~ panel, scales = "free_y", ncol = 2) +
    labs(x = "Relative outbreak day", y = NULL) +
    theme_bw(base_size = 9) +
    theme(strip.text = element_text(size = 7))
}

p_wa  <- fit_vs_data_panel(wa_curve_df,  wa_anchor_df)
p_drc <- fit_vs_data_panel(drc_curve_df, drc_anchor_df, extra_points = drc_sdb_df)

# (A) West Africa and (B) DRC conflict, side by side.   display only; not saved
print((p_wa | p_drc) +
        plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")"))
