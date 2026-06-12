# ============================================================
# Hypothetical West Africa-with-conflict++++ scenario
# Extended-duration version
#
# Purpose
# -------
# This script constructs a fifth scenario for the BP input matrix:
#   West Africa magnitudes + DRC-like conflict disruption of response timing.
#
# It is NOT a new empirical Stan fit, because there are no new data anchors for
# this hypothetical counterfactual. Instead it performs a posterior-preserving
# transform:
#   1. use the existing West Africa posterior Q_j(t) curves and fitted
#      parameter start/end magnitudes;
#   2. use the final DRC conflict Q(t) as a conflict-disruption modulator;
#   3. extend the West Africa response timeline using the empirical ratio
#      duration(DRC conflict) / duration(DRC no-conflict), if available;
#   4. construct a new Q_j(t) over the extended West Africa timeline;
#   5. map that Q_j(t) onto the original West Africa fitted magnitudes;
#   6. export plots, CSVs, and optionally add a fifth sheet to the Excel
#      scenario-matrix workbook.
#
# Key interpretation
# ------------------
# Q is a relative response-maturity index. Parameter magnitudes are inherited
# from the original West Africa fit. Conflict changes timing/shape/duration,
# not the underlying West Africa magnitude bounds.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(stringr)
  library(tibble)
})

# ------------------------------------------------------------
# User settings
# ------------------------------------------------------------
input_dir <- "."

# Latest West Africa outputs from:
# West_Africa_USE_partial_pooling_independent_bounds_UFaCD_tweaked.R
wa_out_prefix <- "west_africa_partial_pool_normalisedQ_estimated_bounds_tweaked_ufc"
wa_q_path     <- file.path(input_dir, paste0(wa_out_prefix, "_Q_summaries.csv"))
wa_curve_path <- file.path(input_dir, paste0(wa_out_prefix, "_curve_summaries.csv"))

# Final DRC conflict outputs from:
# drc_conflict_warsame_weekly_Q_relative_absolute_UFC_rolling4_sparse_greydots.R
# Prefer Q_summaries, fall back to shared Q points if present.
drc_conflict_prefix <- "drc_conflict_warsame_weekly_Q_relative_absolute_UFC_rolling4_sparse_greydots"
drc_q_summaries_path <- file.path(input_dir, paste0(drc_conflict_prefix, "_Q_summaries.csv"))
drc_q_points_path    <- file.path(input_dir, paste0(drc_conflict_prefix, "_shared_conflict_Q_points.csv"))

# Optional BP matrices used only to estimate the duration-extension ratio.
# If these are not found, duration_multiplier_manual is used instead.
drc_conflict_bp_path <- file.path(input_dir, paste0(drc_conflict_prefix, "_bp_input_matrix.csv"))
drc_no_conflict_bp_path <- file.path(
  input_dir,
  "drc_non_conflict_warsame_only_Q_minus_conflict_time_original_plateau_bp_input_matrix.csv"
)

use_drc_duration_ratio <- TRUE
# Used if DRC BP matrices are not available.
duration_multiplier_manual <- 1.50

# Number of output time points for the extended WA-with-conflict matrix.
# If NULL, uses at least as many points as the original West Africa Q grid.
n_pred_extended <- NULL

# Existing scenario-matrix workbook to update.
update_workbook <- FALSE
matrix_workbook_path <- file.path(
  input_dir,
  "final_four_scenario_matrices_with 0_1_Qcurves.xlsx"
)
output_workbook_path <- file.path(
  input_dir,
  "final_five_scenario_matrices_with_extended_WA_conflict.xlsx"
)

new_sheet_name    <- "Worst_WestAfrica_Conflict_PlusPlus"
new_scenario_name <- "Worst_WestAfrica_Conflict_PlusPlus"

existing_scenario_sheets <- c(
  "Worst_WestAfrica",
  "Middle_DRC_NoConflict",
  "Middle_DRC_ConflictSmoothed",
  "Best_Composite_EastAfrica"
)
existing_scenario_names <- c(
  "worst_west_africa",
  "middle_drc_no_conflict",
  "middle_drc_conflict",
  "best_composite_east_africa"
)

out_prefix <- "west_africa_conflict_plusplus_original_method_directCollapse"

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
clip01 <- function(x) pmin(pmax(x, 0), 1)

rescale_01 <- function(x) {
  r <- range(x, na.rm = TRUE)
  if (!all(is.finite(r))) stop("Cannot rescale: non-finite range.")
  if (diff(r) <= 0) stop("Cannot rescale: zero-width range.")
  (x - r[1]) / diff(r)
}

check_file <- function(path, label) {
  if (!file.exists(path)) {
    stop(
      label, " not found at: ", path,
      "\nRun the upstream script first or update the User settings block."
    )
  }
  invisible(path)
}

normalise_q_support <- function(df, tau_col, q_col, label) {
  out <- df %>%
    dplyr::transmute(
      tau = as.numeric(.data[[tau_col]]),
      q   = as.numeric(.data[[q_col]])
    ) %>%
    dplyr::filter(is.finite(tau), is.finite(q)) %>%
    dplyr::arrange(tau) %>%
    dplyr::group_by(tau) %>%
    dplyr::summarise(q = mean(q, na.rm = TRUE), .groups = "drop")

  if (nrow(out) < 2) stop(label, " has fewer than two usable Q support points.")

  tau_range <- range(out$tau, na.rm = TRUE)
  if (tau_range[1] < -1e-8 || tau_range[2] > 1 + 1e-8) {
    out <- out %>% dplyr::mutate(tau = rescale_01(tau))
  }

  # Ensure interpolation endpoints are present.
  if (min(out$tau, na.rm = TRUE) > 1e-8) {
    out <- dplyr::bind_rows(tibble(tau = 0, q = 0), out)
  }
  if (max(out$tau, na.rm = TRUE) < 1 - 1e-8) {
    last_q <- out$q[which.max(out$tau)]
    out <- dplyr::bind_rows(out, tibble(tau = 1, q = last_q))
  }

  out <- out %>%
    dplyr::arrange(tau) %>%
    dplyr::group_by(tau) %>%
    dplyr::summarise(q = mean(q, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(q = clip01(q))

  # Max-scale only. Do not subtract the minimum: this preserves late partial
  # deterioration rather than redefining it as baseline failure.
  q_max <- max(out$q, na.rm = TRUE)
  if (!is.finite(q_max) || q_max <= 0) {
    stop(label, " cannot be normalised: maximum Q is not positive.")
  }

  out %>% dplyr::mutate(q = clip01(q / q_max))
}

make_interp <- function(x, y) {
  stats::approxfun(
    x = as.numeric(x),
    y = as.numeric(y),
    method = "linear",
    rule = 2,
    ties = "ordered"
  )
}

# ------------------------------------------------------------
# Read West Africa Q and fitted parameter magnitudes
# ------------------------------------------------------------
check_file(wa_q_path, "West Africa Q summaries")
check_file(wa_curve_path, "West Africa curve summaries")

wa_q <- readr::read_csv(wa_q_path, show_col_types = FALSE) %>%
  dplyr::mutate(
    parameter    = as.character(parameter),
    param_id     = as.integer(param_id),
    grid_id      = as.integer(grid_id),
    tau          = as.numeric(tau),
    relative_day = as.numeric(relative_day),
    mean         = as.numeric(mean),
    median       = as.numeric(median),
    q5           = as.numeric(q5),
    q95          = as.numeric(q95),
    panel_title  = as.character(panel_title)
  ) %>%
  dplyr::arrange(param_id, tau)

if (any(!is.finite(wa_q$mean))) stop("West Africa Q mean contains non-finite values.")
if (any(wa_q$mean < -1e-6 | wa_q$mean > 1 + 1e-6)) {
  stop("West Africa Q mean contains values outside [0,1].")
}

wa_curve <- readr::read_csv(wa_curve_path, show_col_types = FALSE) %>%
  dplyr::mutate(
    parameter    = as.character(parameter),
    param_id     = as.integer(param_id),
    grid_id      = as.integer(grid_id),
    tau          = as.numeric(tau),
    relative_day = as.numeric(relative_day),
    mean         = as.numeric(mean),
    median       = as.numeric(median),
    q5           = as.numeric(q5),
    q95          = as.numeric(q95),
    panel_title  = as.character(panel_title)
  )

wa_original_duration <- max(wa_q$relative_day, na.rm = TRUE)
if (!is.finite(wa_original_duration) || wa_original_duration <= 0) {
  stop("Could not determine original West Africa duration from Q summaries.")
}

# ------------------------------------------------------------
# Read DRC conflict Q modulator
# ------------------------------------------------------------
if (file.exists(drc_q_summaries_path)) {
  drc_q_raw <- readr::read_csv(drc_q_summaries_path, show_col_types = FALSE)
  if (!all(c("tau", "mean") %in% names(drc_q_raw))) {
    stop("DRC Q summaries file exists but lacks tau and mean columns: ", drc_q_summaries_path)
  }
  drc_q <- normalise_q_support(drc_q_raw, "tau", "mean", "DRC conflict Q summaries")
} else if (file.exists(drc_q_points_path)) {
  drc_q_raw <- readr::read_csv(drc_q_points_path, show_col_types = FALSE)
  if (all(c("tau_q", "q_conflict_shape") %in% names(drc_q_raw))) {
    drc_q <- normalise_q_support(drc_q_raw, "tau_q", "q_conflict_shape", "DRC conflict Q points")
  } else if (all(c("tau", "q_conflict_shape") %in% names(drc_q_raw))) {
    drc_q <- normalise_q_support(drc_q_raw, "tau", "q_conflict_shape", "DRC conflict Q points")
  } else {
    stop("DRC Q points file exists, but expected tau_q/q_conflict_shape or tau/q_conflict_shape columns.")
  }
} else {
  stop(
    "Could not find DRC conflict Q input. Expected one of:\n",
    " - ", drc_q_summaries_path, "\n",
    " - ", drc_q_points_path
  )
}

drc_interp <- make_interp(drc_q$tau, drc_q$q)

# ------------------------------------------------------------
# Determine duration extension
# ------------------------------------------------------------
duration_multiplier <- duration_multiplier_manual
duration_source <- "manual"

if (isTRUE(use_drc_duration_ratio) && file.exists(drc_conflict_bp_path) && file.exists(drc_no_conflict_bp_path)) {
  drc_conflict_bp <- readr::read_csv(drc_conflict_bp_path, show_col_types = FALSE)
  drc_no_conflict_bp <- readr::read_csv(drc_no_conflict_bp_path, show_col_types = FALSE)

  if ("relative_day" %in% names(drc_conflict_bp) && "relative_day" %in% names(drc_no_conflict_bp)) {
    drc_conflict_duration <- max(as.numeric(drc_conflict_bp$relative_day), na.rm = TRUE)
    drc_no_conflict_duration <- max(as.numeric(drc_no_conflict_bp$relative_day), na.rm = TRUE)

    if (is.finite(drc_conflict_duration) && is.finite(drc_no_conflict_duration) && drc_no_conflict_duration > 0) {
      duration_multiplier <- drc_conflict_duration / drc_no_conflict_duration
      duration_source <- paste0(
        "DRC conflict/no-conflict duration ratio = ",
        signif(drc_conflict_duration, 4), " / ", signif(drc_no_conflict_duration, 4)
      )
    }
  }
}

if (!is.finite(duration_multiplier) || duration_multiplier < 1) {
  warning("Invalid duration multiplier; resetting to manual value.")
  duration_multiplier <- duration_multiplier_manual
  duration_source <- "manual fallback"
}

wa_conflict_duration <- wa_original_duration * duration_multiplier

cat("\nDuration settings:\n")
cat(" - Original West Africa duration: ", signif(wa_original_duration, 5), " days\n", sep = "")
cat(" - Duration multiplier: ", signif(duration_multiplier, 5), " (", duration_source, ")\n", sep = "")
cat(" - Extended WA-with-conflict duration: ", signif(wa_conflict_duration, 5), " days\n", sep = "")

n_grid_wa <- length(unique(wa_q$tau))
if (is.null(n_pred_extended)) {
  n_pred_extended <- max(150L, n_grid_wa)
}

extended_grid <- tibble::tibble(
  grid_id = seq_len(n_pred_extended),
  relative_day = seq(0, wa_conflict_duration, length.out = n_pred_extended),
  tau = relative_day / wa_conflict_duration,
  # The intrinsic West Africa response clock progresses on the original
  # West Africa timescale and then plateaus.
  tau_wa_progress = pmin(relative_day / wa_original_duration, 1),
  # The conflict modulator evolves across the extended conflict-affected window.
  tau_drc_conflict = tau
)

# ------------------------------------------------------------
# Construct extended WA-with-conflict Q_j curves
# ------------------------------------------------------------
param_info <- wa_q %>%
  dplyr::distinct(parameter, param_id, panel_title) %>%
  dplyr::arrange(param_id)

q_hybrid_list <- lapply(seq_len(nrow(param_info)), function(i) {
  p <- param_info$parameter[[i]]
  p_df <- wa_q %>% dplyr::filter(parameter == p) %>% dplyr::arrange(tau)

  f_mean <- make_interp(p_df$tau, clip01(p_df$mean))
  f_med  <- make_interp(p_df$tau, clip01(p_df$median))
  f_q5   <- make_interp(p_df$tau, clip01(p_df$q5))
  f_q95  <- make_interp(p_df$tau, clip01(p_df$q95))

  out <- extended_grid %>%
    dplyr::mutate(
      parameter = p,
      param_id = param_info$param_id[[i]],
      panel_title = param_info$panel_title[[i]],
      q_west_africa_mean   = clip01(f_mean(tau_wa_progress)),
      q_west_africa_median = clip01(f_med(tau_wa_progress)),
      q_west_africa_q5     = clip01(f_q5(tau_wa_progress)),
      q_west_africa_q95    = clip01(f_q95(tau_wa_progress)),
      q_drc_conflict_modulator = clip01(drc_interp(tau_drc_conflict)),
      q_raw_mean   = q_west_africa_mean   * q_drc_conflict_modulator,
      q_raw_median = q_west_africa_median * q_drc_conflict_modulator,
      q_raw_q5     = q_west_africa_q5     * q_drc_conflict_modulator,
      q_raw_q95    = q_west_africa_q95    * q_drc_conflict_modulator
    )

  # Scale by the maximum of the mean hybrid curve only. This keeps Q relative
  # and preserves late deterioration if the conflict modulator falls again.
  q_scale <- max(out$q_raw_mean, na.rm = TRUE)
  if (!is.finite(q_scale) || q_scale <= 0) stop("Hybrid Q has non-positive maximum for parameter: ", p)

  out %>%
    dplyr::mutate(
      q_hybrid_scale = q_scale,
      mean   = clip01(q_raw_mean / q_scale),
      median = clip01(q_raw_median / q_scale),
      q5     = clip01(q_raw_q5 / q_scale),
      q95    = clip01(q_raw_q95 / q_scale),
      mean   = dplyr::if_else(grid_id == min(grid_id), 0, mean),
      median = dplyr::if_else(grid_id == min(grid_id), 0, median),
      q5     = dplyr::if_else(grid_id == min(grid_id), 0, q5),
      q95    = dplyr::if_else(grid_id == min(grid_id), 0, q95)
    )
})

q_hybrid <- dplyr::bind_rows(q_hybrid_list)

# Defensive repair: q_hybrid_scale is only a diagnostic scale factor. If the
# column is absent for any reason, reconstruct it from each parameter's raw
# hybrid mean curve before downstream selects.
if (!("q_hybrid_scale" %in% names(q_hybrid))) {
  q_hybrid <- q_hybrid %>%
    dplyr::group_by(parameter) %>%
    dplyr::mutate(q_hybrid_scale = max(q_raw_mean, na.rm = TRUE)) %>%
    dplyr::ungroup()
}

# Fail early with a clearer error if the diagnostic scaling column is still absent
# or unusable before any downstream select() calls.
if (!("q_hybrid_scale" %in% names(q_hybrid))) {
  stop("Internal error: q_hybrid_scale was not created. Check q_hybrid construction.")
}
if (any(!is.finite(q_hybrid$q_hybrid_scale)) || any(q_hybrid$q_hybrid_scale <= 0)) {
  stop("Internal error: q_hybrid_scale contains non-finite or non-positive values.")
}

q_hybrid <- q_hybrid %>%
  dplyr::mutate(
    q5 = pmin(q5, q95),
    q95 = pmax(q5, q95)
  )

q_range_check <- q_hybrid %>%
  dplyr::group_by(parameter) %>%
  dplyr::summarise(
    q_min = min(mean, na.rm = TRUE),
    q_max = max(mean, na.rm = TRUE),
    q_final = mean[which.max(relative_day)],
    .groups = "drop"
  )
print(q_range_check)

if (any(abs(q_range_check$q_min - 0) > 1e-8) || any(abs(q_range_check$q_max - 1) > 1e-8)) {
  stop("Hybrid Q curves should span 0-1. Check construction.")
}

q_hybrid_out <- q_hybrid %>%
  dplyr::mutate(
    # Store the hybrid-Q rescaling factor in sd purely as a diagnostic placeholder
    # so this output has the same column layout as *_Q_summaries.csv.
    # Keep q_hybrid_scale as its own explicit column too.
    sd = q_hybrid_scale
  ) %>%
  dplyr::select(
    parameter, param_id, grid_id, tau, relative_day,
    mean, median, sd, q5, q95, panel_title,
    tau_wa_progress, tau_drc_conflict,
    q_west_africa_mean, q_drc_conflict_modulator,
    q_raw_mean, q_hybrid_scale
  )

# ------------------------------------------------------------
# Map new Q_j curves onto original West Africa fitted magnitudes
# ------------------------------------------------------------
theta_bounds <- wa_curve %>%
  dplyr::group_by(parameter) %>%
  dplyr::summarise(
    theta_start_mean   = mean[which.min(tau)],
    theta_end_mean     = mean[which.max(tau)],
    theta_start_median = median[which.min(tau)],
    theta_end_median   = median[which.max(tau)],
    theta_start_q5     = q5[which.min(tau)],
    theta_end_q5       = q5[which.max(tau)],
    theta_start_q95    = q95[which.min(tau)],
    theta_end_q95      = q95[which.max(tau)],
    .groups = "drop"
  )

hybrid_curve <- q_hybrid_out %>%
  dplyr::left_join(theta_bounds, by = "parameter") %>%
  dplyr::mutate(
    mean = theta_start_mean + (theta_end_mean - theta_start_mean) * mean,
    median = theta_start_median + (theta_end_median - theta_start_median) * median,
    curve_q5_raw = theta_start_q5 + (theta_end_q5 - theta_start_q5) * q5,
    curve_q95_raw = theta_start_q95 + (theta_end_q95 - theta_start_q95) * q95,
    q5 = pmin(curve_q5_raw, curve_q95_raw),
    q95 = pmax(curve_q5_raw, curve_q95_raw)
  ) %>%
  dplyr::select(
    parameter, param_id, grid_id, tau, relative_day,
    mean, median, sd, q5, q95, panel_title,
    tau_wa_progress, tau_drc_conflict,
    q_west_africa_mean, q_drc_conflict_modulator,
    q_raw_mean, q_hybrid_scale,
    dplyr::starts_with("theta_start"), dplyr::starts_with("theta_end")
  )



# ------------------------------------------------------------
# West Africa conflict++ direct-collapse sensitivity
# ------------------------------------------------------------
# Hypothetical counterpart of DRC conflict++: during the same day-200 to
# day-300 window, force the response index to the poor-response state and force
# parameter curves back to their starting poor-response values. Community unsafe
# funerals are treated as complete SDB collapse and set to 1.
wa_plusplus_window_day <- c(200, 300)
collapse_lookup <- hybrid_curve %>%
  arrange(relative_day) %>%
  group_by(parameter) %>%
  summarise(
    collapse_mean = dplyr::first(mean),
    collapse_median = dplyr::first(median),
    collapse_q5 = dplyr::first(q5),
    collapse_q95 = dplyr::first(q95),
    .groups = "drop"
  ) %>%
  mutate(
    collapse_mean = if_else(parameter == "p_unsafe_funeral_comm", 1, collapse_mean),
    collapse_median = if_else(parameter == "p_unsafe_funeral_comm", 1, collapse_median),
    collapse_q5 = if_else(parameter == "p_unsafe_funeral_comm", 1, collapse_q5),
    collapse_q95 = if_else(parameter == "p_unsafe_funeral_comm", 1, collapse_q95)
  )

hybrid_curve <- hybrid_curve %>%
  left_join(collapse_lookup, by = "parameter") %>%
  mutate(
    is_west_africa_plusplus_collapse_period = is.finite(relative_day) &
      relative_day >= wa_plusplus_window_day[1] & relative_day <= wa_plusplus_window_day[2],
    mean = if_else(is_west_africa_plusplus_collapse_period, collapse_mean, mean),
    median = if_else(is_west_africa_plusplus_collapse_period, collapse_median, median),
    q5 = if_else(is_west_africa_plusplus_collapse_period, collapse_q5, q5),
    q95 = if_else(is_west_africa_plusplus_collapse_period, collapse_q95, q95)
  ) %>%
  select(-starts_with("collapse_"))

q_hybrid_out <- q_hybrid_out %>%
  mutate(
    is_west_africa_plusplus_collapse_period = is.finite(relative_day) &
      relative_day >= wa_plusplus_window_day[1] & relative_day <= wa_plusplus_window_day[2],
    mean = if_else(is_west_africa_plusplus_collapse_period, 0, mean),
    median = if_else(is_west_africa_plusplus_collapse_period, 0, median),
    q5 = if_else(is_west_africa_plusplus_collapse_period, 0, q5),
    q95 = if_else(is_west_africa_plusplus_collapse_period, 0, q95)
  )

# ------------------------------------------------------------
# BP-ready matrix
# ------------------------------------------------------------
q_wide <- q_hybrid_out %>%
  dplyr::select(parameter, relative_day, tau, q = mean) %>%
  tidyr::pivot_wider(names_from = parameter, values_from = q, names_prefix = "q_") %>%
  dplyr::mutate(
    q_value = rowMeans(dplyr::select(., dplyr::starts_with("q_")), na.rm = TRUE),
    q_value = rescale_01(q_value)
  ) %>%
  dplyr::select(relative_day, tau, q_value)

bp_matrix <- hybrid_curve %>%
  dplyr::select(parameter, relative_day, tau, mean) %>%
  tidyr::pivot_wider(names_from = parameter, values_from = mean) %>%
  dplyr::mutate(
    prob_hosp = p_hosp,
    delay_hosp = delay_hosp,
    prob_unsafe_funeral_comm = p_unsafe_funeral_comm,
    prob_unsafe_funeral_hosp = p_unsafe_funeral_hosp,
    prob_unsafe_funeral_etu = 0,
    prop_etu = p_ETU,
    ipc_helper = latent_IPC
  ) %>%
  dplyr::select(
    relative_day,
    tau,
    prob_hosp,
    delay_hosp,
    prob_unsafe_funeral_comm,
    prob_unsafe_funeral_hosp,
    prob_unsafe_funeral_etu,
    prop_etu,
    ipc_helper
  ) %>%
  dplyr::left_join(q_wide, by = c("relative_day", "tau"))

range_check <- bp_matrix %>%
  dplyr::summarise(
    prob_hosp_min = min(prob_hosp, na.rm = TRUE),
    prob_hosp_max = max(prob_hosp, na.rm = TRUE),
    delay_hosp_min = min(delay_hosp, na.rm = TRUE),
    delay_hosp_max = max(delay_hosp, na.rm = TRUE),
    p_ufc_min = min(prob_unsafe_funeral_comm, na.rm = TRUE),
    p_ufc_max = max(prob_unsafe_funeral_comm, na.rm = TRUE),
    p_ufh_min = min(prob_unsafe_funeral_hosp, na.rm = TRUE),
    p_ufh_max = max(prob_unsafe_funeral_hosp, na.rm = TRUE),
    prop_etu_min = min(prop_etu, na.rm = TRUE),
    prop_etu_max = max(prop_etu, na.rm = TRUE),
    ipc_helper_min = min(ipc_helper, na.rm = TRUE),
    ipc_helper_max = max(ipc_helper, na.rm = TRUE),
    q_value_min = min(q_value, na.rm = TRUE),
    q_value_max = max(q_value, na.rm = TRUE)
  )
print(range_check)

prob_cols <- c(
  "prob_hosp", "prob_unsafe_funeral_comm", "prob_unsafe_funeral_hosp",
  "prob_unsafe_funeral_etu", "prop_etu", "ipc_helper", "q_value"
)
if (any(bp_matrix[prob_cols] < -1e-8 | bp_matrix[prob_cols] > 1 + 1e-8, na.rm = TRUE)) {
  stop("At least one probability/index column falls outside [0,1].")
}

# ------------------------------------------------------------
# Plots
# ------------------------------------------------------------
plot_q_df <- q_hybrid_out %>%
  dplyr::select(
    parameter, panel_title, relative_day,
    q_west_africa_mean, q_drc_conflict_modulator, mean
  ) %>%
  tidyr::pivot_longer(
    cols = c(q_west_africa_mean, q_drc_conflict_modulator, mean),
    names_to = "curve",
    values_to = "q"
  ) %>%
  dplyr::mutate(
    curve = dplyr::recode(
      curve,
      q_west_africa_mean = "Original West Africa Q_j(t), stretched/plateaued",
      q_drc_conflict_modulator = "DRC conflict modulator Q_conflict(t)",
      mean = "Extended West Africa-with-conflict++ Q_j(t)"
    ),
    curve = factor(
      curve,
      levels = c(
        "Original West Africa Q_j(t), stretched/plateaued",
        "DRC conflict modulator Q_conflict(t)",
        "Extended West Africa-with-conflict++ Q_j(t)"
      )
    )
  )

p_q <- ggplot(plot_q_df, aes(x = relative_day, y = q, linetype = curve, linewidth = curve)) +
  geom_line(colour = "grey25", alpha = 0.85) +
  geom_line(
    data = plot_q_df %>% dplyr::filter(curve == "Extended West Africa-with-conflict++ Q_j(t)"),
    aes(x = relative_day, y = q),
    inherit.aes = FALSE,
    linewidth = 1.0,
    colour = "#1f77b4"
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_linewidth_manual(values = c(0.55, 0.65, 1.0), guide = "none") +
  facet_wrap(~ panel_title, ncol = 2, scales = "fixed") +
  theme_bw(base_size = 11) +
  labs(
    title = "Extended hypothetical West Africa scenario with DRC conflict modulation",
    subtitle = paste0(
      "Duration multiplier = ", signif(duration_multiplier, 3),
      "; Q_hybrid = Q_WA × Q_DRC_conflict, then max-scaled per parameter"
    ),
    x = "Extended West Africa relative outbreak day",
    y = "Q(t)",
    linetype = NULL
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

print(p_q)

ggsave(
  filename = paste0(out_prefix, "_Qplot.png"),
  plot = p_q,
  width = 12,
  height = 10,
  dpi = 180
)

p_theta <- ggplot(hybrid_curve, aes(x = relative_day, y = mean)) +
  geom_ribbon(aes(ymin = q5, ymax = q95), fill = "#c7dcec", alpha = 0.9) +
  geom_line(linewidth = 0.9, colour = "#1f77b4") +
  facet_wrap(~ panel_title, scales = "free_y", ncol = 2) +
  theme_bw(base_size = 11) +
  labs(
    title = "Extended hypothetical West Africa conflict++ parameter trajectories",
    subtitle = "Original West Africa start/end magnitudes retained; DRC conflict-modulated Q shape plus day-200 to day-300 conflict++ collapse",
    x = "Extended West Africa relative outbreak day",
    y = NULL
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.minor = element_blank()
  )

print(p_theta)

ggsave(
  filename = paste0(out_prefix, "_parameter_curves.png"),
  plot = p_theta,
  width = 12,
  height = 10,
  dpi = 180
)

# ------------------------------------------------------------
# Save CSV outputs
# ------------------------------------------------------------
readr::write_csv(q_hybrid_out, paste0(out_prefix, "_Q_summaries.csv"))
readr::write_csv(drc_q, paste0(out_prefix, "_DRC_conflict_modulator.csv"))
readr::write_csv(q_range_check, paste0(out_prefix, "_Q_range_check.csv"))
readr::write_csv(hybrid_curve, paste0(out_prefix, "_curve_summaries.csv"))
readr::write_csv(bp_matrix, paste0(out_prefix, "_bp_input_matrix.csv"))
readr::write_csv(range_check, paste0(out_prefix, "_bp_range_check.csv"))

# ------------------------------------------------------------
# Optional: update scenario matrix workbook with fifth scenario
# ------------------------------------------------------------
if (isTRUE(update_workbook)) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package 'openxlsx' is required to update the Excel workbook. Install it or set update_workbook <- FALSE.")
  }

  check_file(matrix_workbook_path, "Scenario matrix workbook")

  # IMPORTANT:
  # Do not load the existing workbook and remove/add sheets in place. In some
  # openxlsx versions that can leave stale workbook XML and produce errors such
  # as:
  #   Error in gsub(...): invalid 'replacement' argument
  # when saveWorkbook() is called.
  #
  # Instead, build a clean workbook from scratch: read the existing sheets we
  # want to keep, add/replace the fifth scenario, rebuild Combined_Long, and
  # then save to a new file.
  old_sheet_names <- openxlsx::getSheetNames(matrix_workbook_path)

  sheets_to_replace <- c(new_sheet_name, "Combined_Long", "Notes_WA_Conflict")
  sheets_to_keep <- setdiff(old_sheet_names, sheets_to_replace)

  old_sheet_data <- list()
  for (s in sheets_to_keep) {
    old_sheet_data[[s]] <- openxlsx::readWorkbook(matrix_workbook_path, sheet = s)
  }

  add_df_sheet <- function(wb, sheet_name, df) {
    # Excel sheet names have a hard 31-character limit. Our names are already
    # valid, but this keeps the helper robust if labels are changed later.
    safe_name <- substr(sheet_name, 1, 31)
    openxlsx::addWorksheet(wb, safe_name)
    openxlsx::writeData(wb, sheet = safe_name, x = df)
    openxlsx::freezePane(wb, sheet = safe_name, firstActiveRow = 2)
    openxlsx::setColWidths(wb, sheet = safe_name, cols = seq_len(ncol(df)), widths = "auto")
    invisible(wb)
  }

  # Rebuild Combined_Long from the scenario sheets in the original workbook plus
  # the new fifth scenario.
  combined_list <- list()
  for (i in seq_along(existing_scenario_sheets)) {
    s <- existing_scenario_sheets[[i]]
    scenario_nm <- existing_scenario_names[[i]]

    if (s %in% names(old_sheet_data)) {
      df <- old_sheet_data[[s]]
    } else if (s %in% old_sheet_names) {
      df <- openxlsx::readWorkbook(matrix_workbook_path, sheet = s)
    } else {
      warning("Expected scenario sheet not found in workbook and will be skipped: ", s)
      next
    }

    df <- df %>%
      dplyr::select(-dplyr::any_of("scenario")) %>%
      dplyr::mutate(scenario = scenario_nm, .before = 1)

    combined_list[[scenario_nm]] <- df
  }

  combined_list[[new_scenario_name]] <- bp_matrix %>%
    dplyr::mutate(scenario = new_scenario_name, .before = 1)

  combined_long <- dplyr::bind_rows(combined_list)

  notes <- tibble::tibble(
    item = c(
      "New scenario",
      "Interpretation",
      "Duration extension",
      "Q construction",
      "Parameter magnitudes",
      "Empirical points",
      "q_value column",
      "Generated by script"
    ),
    detail = c(
      paste0(new_sheet_name, " / ", new_scenario_name),
      "Hypothetical West Africa outbreak with DRC-like conflict disruption imposed on response maturation.",
      paste0(
        "West Africa response window extended by multiplier ", signif(duration_multiplier, 4),
        " (", duration_source, ")."
      ),
      "Q_hybrid_j(t) = Q_WA_j(t_on_original_WA_clock) * Q_DRC_conflict(t_on_extended_clock), then each parameter-specific Q_hybrid_j is max-scaled to span 0-1.",
      "Original West Africa fitted start/end magnitudes are retained; the new Q curves alter timing/shape/duration, not magnitude endpoints.",
      "No new empirical points are added because this is a hypothetical counterfactual scenario.",
      "q_value is diagnostic only: row-wise mean of parameter-specific Q_hybrid_j values, rescaled to 0-1. BP inputs are the parameter columns.",
      "west_africa_with_conflict_extended_duration_update_matrix_FIX3.R"
    )
  )

  # Create a clean workbook.
  wb_out <- openxlsx::createWorkbook()

  # Preserve existing non-replaced sheets in original order.
  for (s in sheets_to_keep) {
    add_df_sheet(wb_out, s, old_sheet_data[[s]])
  }

  # Add/replace the fifth scenario and rebuilt diagnostic sheets.
  add_df_sheet(wb_out, new_sheet_name, bp_matrix)
  add_df_sheet(wb_out, "Combined_Long", combined_long)
  add_df_sheet(wb_out, "Notes_WA_Conflict", notes)

  openxlsx::saveWorkbook(wb_out, output_workbook_path, overwrite = TRUE)
}

cat("\nDone. Files written:\n")
cat(" - ", paste0(out_prefix, "_Qplot.png"), "\n", sep = "")
cat(" - ", paste0(out_prefix, "_parameter_curves.png"), "\n", sep = "")
cat(" - ", paste0(out_prefix, "_Q_summaries.csv"), "\n", sep = "")
cat(" - ", paste0(out_prefix, "_curve_summaries.csv"), "\n", sep = "")
cat(" - ", paste0(out_prefix, "_bp_input_matrix.csv"), "\n", sep = "")
cat(" - ", paste0(out_prefix, "_bp_range_check.csv"), "\n", sep = "")
if (isTRUE(update_workbook)) {
  cat(" - ", output_workbook_path, "\n", sep = "")
}
