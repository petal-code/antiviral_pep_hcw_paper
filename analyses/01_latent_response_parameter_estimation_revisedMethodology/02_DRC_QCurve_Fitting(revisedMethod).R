# ============================================================================
# 02_DRC_QCurve_Fitting(revisedMethod).R
# ----------------------------------------------------------------------------
# WHAT THIS SCRIPT DOES
#   Builds the DRC "conflict" and "conflict++" response-parameter curves under
#   the REVISED methodology, and saves them for the combine step (03).
#
# WHY THERE IS NO STAN HERE (and how this differs from the original 02 / Model B)
#   The original-methodology DRC fit (Model B) supplied the empirical Warsame Q
#   as fixed data and ESTIMATED each parameter's two endpoints with Stan. The
#   revised methodology instead LOCKS the endpoints to early/late literature
#   extrema (the same lock_endpoints() rule used for West Africa) and then maps
#   the fixed empirical Q straight onto them:
#
#       theta_j(t) = theta_start_j + (theta_end_j - theta_start_j) * Q(t)
#
#   With both the shape (the empirical Q) and the endpoints fixed, nothing is left
#   to sample - the curve is a deterministic mapping. (A smooth Stan fit would
#   also erase the jagged, conflict-interrupted Q shape we specifically want to
#   keep, which is the same reason Model B never smoothed Q either.)
#
#   THREE special cases, all mirroring the original methodology:
#     * community unsafe funerals follow the ABSOLUTE Warsame proxy (1 - success)
#       directly, not the relative-Q endpoint mapping;
#     * the IPC/PPE index (latent_IPC) is just an ordinary increasing parameter
#       here, so its endpoint mapping latent_IPC(t) = ipc_low + (ipc_high-ipc_low)*Q
#       IS the revised "q-scaled IPC" rule - no separate patch is needed; and
#     * "conflict++" uses the Q series that already has the success->0 collapse
#       baked in over days 200-300 (done in 00), so every parameter collapses to
#       its worst endpoint there automatically (Q = 0 maps to theta_start).
#
# Input  : data-processed/DRC_QCurve_revisedMethod/DRC_QCurve_PreppedData(revisedMethod).rds
#          (from 00; uses $anchors, $conflict_qseries, $conflict_plusplus_qseries)
# Output : data-processed/DRC_QCurve_revisedMethod/DRC_QCurve_Conflict_Fit(revisedMethod).rds
#          data-processed/DRC_QCurve_revisedMethod/DRC_QCurve_ConflictPlusPlus_Fit(revisedMethod).rds
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(ggplot2)   # display-only plots of the curves at the end
})

source(here::here("analyses", "01_latent_response_parameter_estimation_revisedMethodology",
                  "helpers(revisedMethod).R"))

# ----------------------------------------------------------------------------
# 1. Read the bundled DRC inputs
# ----------------------------------------------------------------------------
drc_prep    <- readRDS(file.path(DIR_PROCESSED,
  "DRC_QCurve_revisedMethod/DRC_QCurve_PreppedData(revisedMethod).rds"))
drc_anchors <- drc_prep$anchors

# ----------------------------------------------------------------------------
# 2. Per-parameter metadata and LOCKED endpoints
# ----------------------------------------------------------------------------
# As in West Africa, the only per-parameter inputs the revised mapping needs are
# the direction and the locked endpoints. The DRC horizon spans the literature
# anchors and the conflict SDB series.
drc_max_day <- max(c(drc_anchors$relative_day,
                     drc_prep$conflict_qseries$relative_day), na.rm = TRUE)  # ~457

param_meta <- drc_anchors %>%
  group_by(parameter) %>%
  summarise(
    direction   = first(direction),
    lower_bound = first(lower_bound),
    upper_bound = first(upper_bound),
    .groups = "drop"
  ) %>%
  mutate(param_id = match(parameter, PARAM_LEVELS)) %>%
  arrange(param_id)

# Lock the endpoints. latent_IPC uses a LATER evidence window (from day 275)
# than the other parameters, matching the revised-methodology DRC specification;
# the unsafe-funeral parameters get the terminal-zero treatment where an explicit
# zero is observed. (See lock_endpoints() in helpers for the full rule.)
endpoint_table <- lock_endpoints(
  drc_anchors,
  scenario_duration_days  = drc_max_day,
  early_window_day        = 50,
  late_start_day          = 325,
  late_start_day_by_param = c(latent_IPC = 275)
)

param_meta <- param_meta %>%
  left_join(endpoint_table, by = c("parameter", "direction")) %>%
  arrange(param_id)

message("DRC locked endpoints (latent_IPC start/end IS the q-scaled IPC range):")
print(param_meta %>% select(parameter, direction, theta_start, theta_end,
                            endpoint_day_for_tau, start_source, end_source))

# ----------------------------------------------------------------------------
# 3. Build one DRC scenario (called twice: conflict, then conflict++)
# ----------------------------------------------------------------------------
# Small guard used below: clamp probability-type parameters to [0,1]; leave
# delay_hosp (in days) untouched. Defined here (not in helpers) because it is
# specific to this script's native-unit columns.
clip01_if_prob <- function(parameter, x) {
  if (parameter == "delay_hosp") x else clip01(x)
}

# Given a Q series (the shared empirical curve) this deterministically maps every
# parameter onto its locked endpoints, overrides community unsafe funerals with
# the absolute Warsame proxy, and returns the per-parameter native-unit curves on
# the same day grid plus the shared Q grid (so 03 can consume it exactly like the
# original-methodology DRC fit object).
build_drc_scenario <- function(qseries, label) {

  message("\n==== Building DRC scenario (revised): ", label, " ====")

  # The shared horizon spans both the literature anchors and the SDB Q series.
  max_day <- max(c(drc_anchors$relative_day, qseries$relative_day), na.rm = TRUE)

  # The grid is the empirical Q support (it already includes the forced day-0
  # start). The mapping is piecewise-linear in Q, so these support points capture
  # the curve exactly; 03 interpolates between them onto the 0..730 day axis.
  grid_day <- qseries$relative_day
  q_value  <- clip01(qseries$q_value)

  # Absolute community unsafe-funeral proxy (1 - success) on the same grid.
  ufc_proxy <- clip01(qseries$unsafe_funeral_comm_proxy)

  # Map each parameter onto its locked endpoints.
  curve_summ <- lapply(seq_len(nrow(param_meta)), function(i) {
    pm <- param_meta[i, ]
    p  <- pm$parameter

    if (p == "p_unsafe_funeral_comm") {
      # Special case: follow the absolute Warsame proxy directly (its floor is
      # 1 - max(success), not 0), exactly as in the original methodology.
      mean_vec <- ufc_proxy
    } else {
      # Deterministic endpoint mapping. Direction is already baked into the
      # endpoints, so Q = 0 -> worst (theta_start), Q = 1 -> best (theta_end).
      # For latent_IPC this IS the revised q-scaled IPC rule.
      mean_vec <- pm$theta_start + (pm$theta_end - pm$theta_start) * q_value

      # Terminal-zero plateau: if the end endpoint is a forced zero reached before
      # the horizon, hold it at zero from that day on.
      if (!is.na(pm$end_source) && pm$end_source == "terminal_zero_anchor" &&
          is.finite(pm$endpoint_day_for_tau)) {
        mean_vec[grid_day >= pm$endpoint_day_for_tau] <- pm$theta_end
      }
    }

    tibble(
      parameter    = p,
      param_id     = pm$param_id,
      tau          = grid_day / max_day,
      relative_day = grid_day,
      mean         = clip01_if_prob(p, mean_vec),
      # The empirical Q is treated as FIXED, so there is no posterior spread to
      # report here; q5/q95 equal the mean (03 uses only `mean` anyway).
      q5           = clip01_if_prob(p, mean_vec),
      q95          = clip01_if_prob(p, mean_vec)
    )
  }) %>% bind_rows() %>% arrange(param_id, relative_day)

  # The shared Q on the grid; 03 uses this as the scenario's q_value.
  q_grid <- tibble(relative_day = grid_day, tau = grid_day / max_day, q_value = q_value)

  list(curve_summ = curve_summ, q_grid = q_grid,
       param_meta = param_meta, endpoint_table = endpoint_table, max_day = max_day)
}

# ----------------------------------------------------------------------------
# 4. Build both conflict scenarios and save
# ----------------------------------------------------------------------------
# The two Q series differ only in that the "++" one has the forced collapse baked
# in (done in 00); the SAME deterministic mapping is applied to each.
drc_conflict_fit          <- build_drc_scenario(drc_prep$conflict_qseries,          "drc_conflict")
drc_conflict_plusplus_fit <- build_drc_scenario(drc_prep$conflict_plusplus_qseries, "drc_conflict_plusplus")

saveRDS(drc_conflict_fit, file.path(DIR_PROCESSED,
        "DRC_QCurve_revisedMethod/DRC_QCurve_Conflict_Fit(revisedMethod).rds"))
saveRDS(drc_conflict_plusplus_fit, file.path(DIR_PROCESSED,
        "DRC_QCurve_revisedMethod/DRC_QCurve_ConflictPlusPlus_Fit(revisedMethod).rds"))

# ----------------------------------------------------------------------------
# 5. Plot each scenario with its data on top  (display only)
# ----------------------------------------------------------------------------
# One facet per parameter: the deterministic mean curve (blue), the literature
# anchors (orange), and - in the community unsafe-funeral panel only - the SDB
# data points (1 - success) that drive that curve (grey). The plots are printed
# to the active graphics device and deliberately NOT saved.
plot_drc_scenario <- function(fit, qseries, label) {

  curve_plot_df <- fit$curve_summ %>%
    mutate(panel = factor(PANEL_LOOKUP[parameter], levels = unname(PANEL_LOOKUP)))

  anchor_plot_df <- drc_anchors %>%
    mutate(panel = factor(PANEL_LOOKUP[parameter], levels = unname(PANEL_LOOKUP)))

  sdb_plot_df <- qseries %>%
    filter(n_eligible_sum > 0) %>%
    transmute(relative_day,
              value = unsafe_funeral_comm_proxy,
              panel = factor(PANEL_LOOKUP[["p_unsafe_funeral_comm"]], levels = unname(PANEL_LOOKUP)))

  ggplot(curve_plot_df, aes(relative_day, mean)) +
    geom_line(colour = "#1f77b4", linewidth = 0.9) +
    geom_point(data = sdb_plot_df, aes(relative_day, value),
               inherit.aes = FALSE, colour = "grey55", size = 1, alpha = 0.7) +
    geom_point(data = anchor_plot_df, aes(relative_day, value_used),
               inherit.aes = FALSE, colour = "#ff7f0e", size = 2) +
    facet_wrap(~ panel, scales = "free_y", ncol = 2) +
    labs(title = paste0("DRC ", label, " (revised): endpoint-mapped parameter curves"),
         subtitle = "Blue = empirical-Q mapped onto locked endpoints; orange = literature anchors; grey = SDB community proxy",
         x = "Relative outbreak day", y = NULL) +
    theme_bw(base_size = 11) +
    theme(strip.text = element_text(face = "bold"))
}

print(plot_drc_scenario(drc_conflict_fit,          drc_prep$conflict_qseries,          "conflict"))     # display only
print(plot_drc_scenario(drc_conflict_plusplus_fit, drc_prep$conflict_plusplus_qseries, "conflict++"))   # display only
