# ============================================================
# West Africa with conflict++: original hybrid-Q logic, updated to
# latest endpoint-constrained West Africa magnitudes
# ============================================================
#
# Purpose
# -------
# Construct the hypothetical West Africa-with-conflict++ scenario using the
# original counterfactual logic that gave the sensible conflict-modulated shape:
#
#   Q_hybrid_j(t) =
#       Q_WA_j(t on the original West Africa response clock)
#       × Q_DRC_conflict(t on the extended conflict clock)
#   followed by max-scaling within each parameter so Q_hybrid_j spans 0–1.
#
# This is the key correction relative to the previous "latest WA" attempt:
#   - DO NOT evaluate Q_WA_j on the extended 0–1 timeline.
#   - Instead evaluate Q_WA_j on the original West Africa response clock:
#         tau_wa_progress = min(relative_day_extended / WA_original_duration, 1)
#   - This lets the West Africa response mature on its original timescale,
#     while the DRC conflict modulator creates interruptions across the
#     extended counterfactual duration.
#
# Interpretation
# --------------
#   West Africa supplies: endpoint-constrained start/end magnitudes and
#                         baseline response shape.
#   DRC conflict supplies: disruption/modulation and duration extension.
#
# This scenario has no independent empirical anchors, so the early/late
# literature-extrema endpoint logic is NOT applied here directly.
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

# Latest West Africa outputs from script 01
wa_out_prefix <- "worst_west_africa_endpoint_constrained_zero_plateau"
wa_q_path     <- file.path(input_dir, paste0(wa_out_prefix, "_Q_summaries.csv"))
wa_curve_path <- file.path(input_dir, paste0(wa_out_prefix, "_curve_summaries.csv"))
wa_matrix_path <- file.path(input_dir, paste0(wa_out_prefix, "_matrix.csv"))

# Preferred DRC conflict modulator sources.
# The old-style Q summaries/shared-Q files reproduce the original method most
# closely. The current endpoint-constrained conflict matrix is allowed as a
# fallback.
drc_conflict_prefix_old <- "drc_conflict_warsame_weekly_Q_relative_absolute_UFC_rolling4_sparse_greydots"
drc_q_summaries_path_old <- file.path(input_dir, paste0(drc_conflict_prefix_old, "_Q_summaries.csv"))
drc_q_points_path_old    <- file.path(input_dir, paste0(drc_conflict_prefix_old, "_shared_conflict_Q_points.csv"))

drc_conflict_matrix_candidates <- c(
  file.path(input_dir, "drc_conflict_endpoint_constrained_preserveQ_greyWarsame_CI_trimmedUFC_matrix.csv"),
  file.path(input_dir, "drc_conflict_endpoint_constrained_preserveQ_greyWarsame_CI_matrix.csv"),
  file.path(input_dir, "drc_conflict_endpoint_constrained_zero_plateau_matrix.csv"),
  file.path(input_dir, "drc_conflict_endpoint_constrained_preserveQ_directUFC_greyPoints_matrix.csv")
)

# DRC no-conflict output from script 02, used only to estimate duration multiplier.
drc_no_conflict_matrix_candidates <- c(
  file.path(input_dir, "drc_no_conflict_endpoint_constrained_zero_plateau_matrix.csv"),
  file.path(input_dir, "drc_no_conflict_endpoint_constrained_short_horizon_filtered_matrix.csv"),
  file.path(input_dir, "drc_non_conflict_allparams_minus_conflict_time_corrected_warsame_bp_input_matrix.csv")
)

use_drc_duration_ratio <- TRUE
duration_multiplier_manual <- 2.37

# Number of output time points. NULL = preserve at least original WA Q resolution,
# scaled by duration multiplier.
n_pred_extended <- NULL

out_prefix <- "west_africa_with_conflict_plusplus_originalHybridQ_latestWA_completeCI"
scenario_name <- "Worst_WestAfrica_Conflict_PlusPlus"

# Preserve endpoint-constrained West Africa zero plateaus, scaled in time.
preserve_scaled_terminal_zero_plateaus <- TRUE
zero_tolerance <- 1e-10

# If the selected DRC conflict modulator source does not contain q5/q95-style
# uncertainty columns, use this approximate uncertainty band around Q_DRC.
# This is only a plotting/uncertainty-propagation fallback; the mean curve and
# BP matrix values are unchanged.
drc_modulator_ci_half_width_fallback <- 0.12

# Optional workbook update is deliberately off by default. This avoids accidental
# workbook corruption while iterating. The CSV is the important BP input.
update_workbook <- FALSE

matrix_workbook_path <- file.path(input_dir, "final_four_scenario_matrices_with 0_1_Qcurves.xlsx")
output_workbook_path <- file.path(input_dir, "final_five_scenario_matrices_with_WA_conflict_originalHybridQ_latestWA.xlsx")
new_sheet_name <- "Worst_WestAfrica_Conflict_PlusPlus"
new_scenario_name <- "worst_west_africa_conflict_plusplus"

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
      "\nRun upstream scripts first or update the User settings block.",
      call. = FALSE
    )
  }
  invisible(path)
}

first_existing_file <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

normalise_q_support <- function(df, tau_col, q_col, label) {
  out <- df %>%
    transmute(
      tau = as.numeric(.data[[tau_col]]),
      q   = as.numeric(.data[[q_col]])
    ) %>%
    filter(is.finite(tau), is.finite(q)) %>%
    arrange(tau) %>%
    group_by(tau) %>%
    summarise(q = mean(q, na.rm = TRUE), .groups = "drop")

  if (nrow(out) < 2) stop(label, " has fewer than two usable Q support points.")

  tau_range <- range(out$tau, na.rm = TRUE)
  if (tau_range[1] < -1e-8 || tau_range[2] > 1 + 1e-8) {
    out <- out %>% mutate(tau = rescale_01(tau))
  }

  # Ensure interpolation endpoints exist.
  if (min(out$tau, na.rm = TRUE) > 1e-8) {
    out <- bind_rows(tibble(tau = 0, q = 0), out)
  }
  if (max(out$tau, na.rm = TRUE) < 1 - 1e-8) {
    out <- bind_rows(out, tibble(tau = 1, q = out$q[which.max(out$tau)]))
  }

  out <- out %>%
    arrange(tau) %>%
    group_by(tau) %>%
    summarise(q = mean(q, na.rm = TRUE), .groups = "drop") %>%
    mutate(q = clip01(q))

  # Original DRC logic used max-scaling only, not min-max scaling.
  q_max <- max(out$q, na.rm = TRUE)
  if (!is.finite(q_max) || q_max <= 0) {
    stop(label, " cannot be normalised: maximum Q is not positive.")
  }

  out %>% mutate(q = clip01(q / q_max))
}


first_matching_col <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0) return(NULL)
  hit[[1]]
}

normalise_q_support_complete <- function(df, tau_col, mean_col, label,
                                         low_col = NULL,
                                         high_col = NULL,
                                         fallback_half_width = 0.12) {
  if (!tau_col %in% names(df)) stop(label, " lacks tau column: ", tau_col)
  if (!mean_col %in% names(df)) stop(label, " lacks mean/Q column: ", mean_col)

  q_low_expr <- if (!is.null(low_col) && low_col %in% names(df)) {
    as.numeric(df[[low_col]])
  } else {
    rep(NA_real_, nrow(df))
  }

  q_high_expr <- if (!is.null(high_col) && high_col %in% names(df)) {
    as.numeric(df[[high_col]])
  } else {
    rep(NA_real_, nrow(df))
  }

  out <- tibble(
    tau = as.numeric(df[[tau_col]]),
    q_mean = as.numeric(df[[mean_col]]),
    q_low = q_low_expr,
    q_high = q_high_expr
  ) %>%
    filter(is.finite(tau), is.finite(q_mean)) %>%
    arrange(tau) %>%
    group_by(tau) %>%
    summarise(
      q_mean = mean(q_mean, na.rm = TRUE),
      q_low = if (all(is.na(q_low))) NA_real_ else mean(q_low, na.rm = TRUE),
      q_high = if (all(is.na(q_high))) NA_real_ else mean(q_high, na.rm = TRUE),
      .groups = "drop"
    )

  if (nrow(out) < 2) stop(label, " has fewer than two usable Q support points.")

  tau_range <- range(out$tau, na.rm = TRUE)
  if (tau_range[1] < -1e-8 || tau_range[2] > 1 + 1e-8) {
    out <- out %>% mutate(tau = rescale_01(tau))
  }

  # Ensure interpolation endpoints exist. For low/high, endpoint rows inherit
  # the mean unless proper uncertainty is already available.
  if (min(out$tau, na.rm = TRUE) > 1e-8) {
    out <- bind_rows(tibble(tau = 0, q_mean = 0, q_low = 0, q_high = 0), out)
  }
  if (max(out$tau, na.rm = TRUE) < 1 - 1e-8) {
    last_i <- which.max(out$tau)
    out <- bind_rows(
      out,
      tibble(
        tau = 1,
        q_mean = out$q_mean[last_i],
        q_low = out$q_low[last_i],
        q_high = out$q_high[last_i]
      )
    )
  }

  out <- out %>%
    arrange(tau) %>%
    group_by(tau) %>%
    summarise(
      q_mean = mean(q_mean, na.rm = TRUE),
      q_low = if (all(is.na(q_low))) NA_real_ else mean(q_low, na.rm = TRUE),
      q_high = if (all(is.na(q_high))) NA_real_ else mean(q_high, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      q_mean = clip01(q_mean),
      q_low = clip01(q_low),
      q_high = clip01(q_high)
    )

  q_max <- max(out$q_mean, na.rm = TRUE)
  if (!is.finite(q_max) || q_max <= 0) {
    stop(label, " cannot be normalised: maximum Q is not positive.")
  }

  out <- out %>%
    mutate(
      q_mean = clip01(q_mean / q_max),
      q_low = clip01(q_low / q_max),
      q_high = clip01(q_high / q_max)
    )

  # If the source had no usable uncertainty columns, create a full-length
  # fallback band around the DRC modulator so late conflict wiggles are not
  # plotted as artificially certain.
  if (all(is.na(out$q_low)) || all(is.na(out$q_high))) {
    out <- out %>%
      mutate(
        q_low = clip01(q_mean - fallback_half_width),
        q_high = clip01(q_mean + fallback_half_width)
      )
  }

  out %>%
    mutate(
      q_low_tmp = pmin(q_low, q_high, q_mean, na.rm = TRUE),
      q_high_tmp = pmax(q_low, q_high, q_mean, na.rm = TRUE),
      q_low = q_low_tmp,
      q_high = q_high_tmp,
      q = q_mean
    ) %>%
    select(tau, q, q_low, q_high)
}


make_interp <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- as.numeric(x[ok])
  y <- as.numeric(y[ok])
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]

  df <- tibble(x = x, y = y) %>%
    group_by(x) %>%
    summarise(y = mean(y, na.rm = TRUE), .groups = "drop")

  approxfun(df$x, df$y, method = "linear", rule = 2, ties = "ordered")
}

find_zero_plateau_day <- function(df, parameter, zero_tol = 1e-10) {
  if (!parameter %in% names(df) || !"relative_day" %in% names(df)) return(NA_real_)
  vals <- as.numeric(df[[parameter]])
  days <- as.numeric(df$relative_day)
  idx <- which(is.finite(vals) & is.finite(days) & days > 0 & abs(vals) <= zero_tol)
  if (length(idx) == 0) return(NA_real_)
  min(days[idx], na.rm = TRUE)
}

# ------------------------------------------------------------
# Read West Africa Q, curve summaries, and matrix
# ------------------------------------------------------------
check_file(wa_q_path, "Latest West Africa Q summaries")
check_file(wa_curve_path, "Latest West Africa curve summaries")
check_file(wa_matrix_path, "Latest West Africa matrix")

wa_q <- read_csv(wa_q_path, show_col_types = FALSE) %>%
  mutate(
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
  arrange(param_id, tau)

wa_curve <- read_csv(wa_curve_path, show_col_types = FALSE) %>%
  mutate(
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

wa_matrix <- read_csv(wa_matrix_path, show_col_types = FALSE)

if (any(!is.finite(wa_q$mean))) stop("West Africa Q mean contains non-finite values.")
if (any(wa_q$mean < -1e-6 | wa_q$mean > 1 + 1e-6)) {
  stop("West Africa Q mean contains values outside [0,1].")
}

wa_original_duration <- max(wa_q$relative_day, na.rm = TRUE)
if (!is.finite(wa_original_duration) || wa_original_duration <= 0) {
  stop("Could not determine original West Africa duration from Q summaries.")
}

# ------------------------------------------------------------
# Read DRC conflict Q modulator, including uncertainty where available
# ------------------------------------------------------------
drc_q_source <- NA_character_

if (file.exists(drc_q_summaries_path_old)) {
  drc_q_raw <- read_csv(drc_q_summaries_path_old, show_col_types = FALSE)
  if (!all(c("tau", "mean") %in% names(drc_q_raw))) {
    stop("Old DRC Q summaries file exists but lacks tau and mean columns: ", drc_q_summaries_path_old)
  }

  low_col <- first_matching_col(drc_q_raw, c("q5", "q025", "lower", "low", "lwr"))
  high_col <- first_matching_col(drc_q_raw, c("q95", "q975", "upper", "high", "upr"))

  drc_q <- normalise_q_support_complete(
    drc_q_raw, "tau", "mean", "Old DRC conflict Q summaries",
    low_col = low_col,
    high_col = high_col,
    fallback_half_width = drc_modulator_ci_half_width_fallback
  )
  drc_q_source <- drc_q_summaries_path_old

} else if (file.exists(drc_q_points_path_old)) {
  drc_q_raw <- read_csv(drc_q_points_path_old, show_col_types = FALSE)

  if (all(c("tau_q", "q_conflict_shape") %in% names(drc_q_raw))) {
    tau_col <- "tau_q"
    mean_col <- "q_conflict_shape"
  } else if (all(c("tau", "q_conflict_shape") %in% names(drc_q_raw))) {
    tau_col <- "tau"
    mean_col <- "q_conflict_shape"
  } else {
    stop("DRC Q points file exists, but expected tau_q/q_conflict_shape or tau/q_conflict_shape columns.")
  }

  low_col <- first_matching_col(drc_q_raw, c("q_conflict_low", "q_low", "q5", "lower", "low", "lwr"))
  high_col <- first_matching_col(drc_q_raw, c("q_conflict_high", "q_high", "q95", "upper", "high", "upr"))

  drc_q <- normalise_q_support_complete(
    drc_q_raw, tau_col, mean_col, "Old DRC conflict Q points",
    low_col = low_col,
    high_col = high_col,
    fallback_half_width = drc_modulator_ci_half_width_fallback
  )
  drc_q_source <- drc_q_points_path_old

} else {
  # Fallback to latest endpoint-constrained conflict matrix q_value.
  drc_conflict_matrix_path <- first_existing_file(drc_conflict_matrix_candidates)
  if (is.na(drc_conflict_matrix_path)) {
    stop(
      "Could not find any DRC conflict Q input. Tried old Q files and latest conflict matrices.",
      call. = FALSE
    )
  }

  drc_conflict_matrix <- read_csv(drc_conflict_matrix_path, show_col_types = FALSE)
  if (!all(c("tau", "q_value") %in% names(drc_conflict_matrix))) {
    stop("DRC conflict matrix must contain tau and q_value: ", drc_conflict_matrix_path)
  }

  low_col <- first_matching_col(drc_conflict_matrix, c("q_value_low", "q_low", "q5", "lower", "low", "lwr"))
  high_col <- first_matching_col(drc_conflict_matrix, c("q_value_high", "q_high", "q95", "upper", "high", "upr"))

  drc_q <- normalise_q_support_complete(
    drc_conflict_matrix, "tau", "q_value", "Latest DRC conflict matrix q_value",
    low_col = low_col,
    high_col = high_col,
    fallback_half_width = drc_modulator_ci_half_width_fallback
  )
  drc_q_source <- drc_conflict_matrix_path
}

drc_interp <- make_interp(drc_q$tau, drc_q$q)
drc_low_interp <- make_interp(drc_q$tau, drc_q$q_low)
drc_high_interp <- make_interp(drc_q$tau, drc_q$q_high)

# ------------------------------------------------------------
# Determine duration extension
# ------------------------------------------------------------
duration_multiplier <- duration_multiplier_manual
duration_source <- "manual"

if (isTRUE(use_drc_duration_ratio)) {
  # Prefer latest conflict matrix duration if available.
  drc_conflict_duration <- NA_real_
  drc_no_conflict_duration <- NA_real_

  drc_conflict_matrix_path_for_duration <- first_existing_file(drc_conflict_matrix_candidates)
  if (!is.na(drc_conflict_matrix_path_for_duration)) {
    dfc <- read_csv(drc_conflict_matrix_path_for_duration, show_col_types = FALSE)
    if ("relative_day" %in% names(dfc)) {
      drc_conflict_duration <- max(as.numeric(dfc$relative_day), na.rm = TRUE)
    }
  }

  drc_no_conflict_matrix_path <- first_existing_file(drc_no_conflict_matrix_candidates)
  if (!is.na(drc_no_conflict_matrix_path)) {
    dfn <- read_csv(drc_no_conflict_matrix_path, show_col_types = FALSE)
    if ("relative_day" %in% names(dfn)) {
      drc_no_conflict_duration <- max(as.numeric(dfn$relative_day), na.rm = TRUE)
    }
  }

  if (is.finite(drc_conflict_duration) && is.finite(drc_no_conflict_duration) && drc_no_conflict_duration > 0) {
    duration_multiplier <- drc_conflict_duration / drc_no_conflict_duration
    duration_source <- paste0(
      "DRC conflict/no-conflict duration ratio = ",
      signif(drc_conflict_duration, 4), " / ", signif(drc_no_conflict_duration, 4)
    )
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
cat(" - DRC Q source: ", drc_q_source, "\n", sep = "")

n_grid_wa <- length(unique(wa_q$tau))
if (is.null(n_pred_extended)) {
  n_pred_extended <- max(150L, ceiling(n_grid_wa * duration_multiplier))
}

extended_grid <- tibble(
  grid_id = seq_len(n_pred_extended),
  relative_day = seq(0, wa_conflict_duration, length.out = n_pred_extended),
  tau = relative_day / wa_conflict_duration,
  # IMPORTANT: original WA response clock, not stretched over conflict duration.
  tau_wa_progress = pmin(relative_day / wa_original_duration, 1),
  # DRC conflict evolves across the extended conflict-affected window.
  tau_drc_conflict = tau
)

# ------------------------------------------------------------
# Construct original-style extended WA-with-conflict Q_j curves
# ------------------------------------------------------------
param_info <- wa_q %>%
  distinct(parameter, param_id, panel_title) %>%
  arrange(param_id)

q_hybrid_list <- lapply(seq_len(nrow(param_info)), function(i) {
  p <- param_info$parameter[[i]]
  p_df <- wa_q %>% filter(parameter == p) %>% arrange(tau)

  f_mean <- make_interp(p_df$tau, clip01(p_df$mean))
  f_med  <- make_interp(p_df$tau, clip01(p_df$median))
  f_q5   <- make_interp(p_df$tau, clip01(p_df$q5))
  f_q95  <- make_interp(p_df$tau, clip01(p_df$q95))

  out <- extended_grid %>%
    mutate(
      parameter = p,
      param_id = param_info$param_id[[i]],
      panel_title = param_info$panel_title[[i]],

      q_west_africa_mean   = clip01(f_mean(tau_wa_progress)),
      q_west_africa_median = clip01(f_med(tau_wa_progress)),
      q_west_africa_q5     = clip01(f_q5(tau_wa_progress)),
      q_west_africa_q95    = clip01(f_q95(tau_wa_progress)),

      q_drc_conflict_modulator = clip01(drc_interp(tau_drc_conflict)),
      q_drc_conflict_low = clip01(drc_low_interp(tau_drc_conflict)),
      q_drc_conflict_high = clip01(drc_high_interp(tau_drc_conflict)),

      # Propagate both West Africa fit uncertainty and DRC conflict-modulator
      # uncertainty. The DRC uncertainty is what keeps the ribbon visible over
      # the full extended conflict period.
      q_raw_mean   = q_west_africa_mean   * q_drc_conflict_modulator,
      q_raw_median = q_west_africa_median * q_drc_conflict_modulator,

      q_raw_lo_1 = q_west_africa_q5  * q_drc_conflict_low,
      q_raw_lo_2 = q_west_africa_q5  * q_drc_conflict_high,
      q_raw_hi_1 = q_west_africa_q95 * q_drc_conflict_low,
      q_raw_hi_2 = q_west_africa_q95 * q_drc_conflict_high,

      q_raw_q5 = pmin(q_raw_lo_1, q_raw_lo_2, q_raw_hi_1, q_raw_hi_2, na.rm = TRUE),
      q_raw_q95 = pmax(q_raw_lo_1, q_raw_lo_2, q_raw_hi_1, q_raw_hi_2, na.rm = TRUE)
    )

  # Scale by max of mean hybrid curve per parameter, as in the original script.
  q_scale <- max(out$q_raw_mean, na.rm = TRUE)
  if (!is.finite(q_scale) || q_scale <= 0) {
    stop("Hybrid Q has non-positive maximum for parameter: ", p)
  }

  out %>%
    mutate(
      q_hybrid_scale = q_scale,
      mean   = clip01(q_raw_mean / q_scale),
      median = clip01(q_raw_median / q_scale),
      q5     = clip01(q_raw_q5 / q_scale),
      q95    = clip01(q_raw_q95 / q_scale),
      mean   = if_else(grid_id == min(grid_id), 0, mean),
      median = if_else(grid_id == min(grid_id), 0, median),
      q5     = if_else(grid_id == min(grid_id), 0, q5),
      q95    = if_else(grid_id == min(grid_id), 0, q95)
    )
})

q_hybrid <- bind_rows(q_hybrid_list) %>%
  mutate(
    q5_tmp = pmin(q5, q95),
    q95_tmp = pmax(q5, q95),
    q5 = q5_tmp,
    q95 = q95_tmp
  ) %>%
  select(-q5_tmp, -q95_tmp)

q_range_check <- q_hybrid %>%
  group_by(parameter) %>%
  summarise(
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
  mutate(sd = q_hybrid_scale) %>%
  select(
    parameter, param_id, grid_id, tau, relative_day,
    mean, median, sd, q5, q95, panel_title,
    tau_wa_progress, tau_drc_conflict,
    q_west_africa_mean, q_drc_conflict_modulator,
    q_drc_conflict_low, q_drc_conflict_high,
    q_raw_mean, q_hybrid_scale
  )

# ------------------------------------------------------------
# Map new Q_j curves onto latest endpoint-constrained West Africa magnitudes
# ------------------------------------------------------------
theta_bounds <- wa_curve %>%
  group_by(parameter) %>%
  summarise(
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
  left_join(theta_bounds, by = "parameter") %>%
  mutate(
    mean = theta_start_mean + (theta_end_mean - theta_start_mean) * mean,
    median = theta_start_median + (theta_end_median - theta_start_median) * median,
    curve_q5_raw = theta_start_q5 + (theta_end_q5 - theta_start_q5) * q5,
    curve_q95_raw = theta_start_q95 + (theta_end_q95 - theta_start_q95) * q95,
    q5 = pmin(curve_q5_raw, curve_q95_raw),
    q95 = pmax(curve_q5_raw, curve_q95_raw)
  )

# Preserve terminal-zero plateaus from endpoint-constrained West Africa.
if (isTRUE(preserve_scaled_terminal_zero_plateaus)) {
  zero_params <- unique(hybrid_curve$parameter[str_detect(hybrid_curve$parameter, "unsafe_funeral")])

  for (p in zero_params) {
    zero_day_wa <- find_zero_plateau_day(wa_matrix, p, zero_tolerance)
    if (is.finite(zero_day_wa)) {
      zero_day_ext <- zero_day_wa * duration_multiplier

      hybrid_curve <- hybrid_curve %>%
        mutate(
          mean = if_else(parameter == p & relative_day >= zero_day_ext, 0, mean),
          median = if_else(parameter == p & relative_day >= zero_day_ext, 0, median),
          q5 = if_else(parameter == p & relative_day >= zero_day_ext, 0, q5),
          q95 = if_else(parameter == p & relative_day >= zero_day_ext, 0, q95)
        )

      q_hybrid_out <- q_hybrid_out %>%
        mutate(
          mean = if_else(parameter == p & relative_day >= zero_day_ext, 1, mean),
          median = if_else(parameter == p & relative_day >= zero_day_ext, 1, median),
          q5 = if_else(parameter == p & relative_day >= zero_day_ext, 1, q5),
          q95 = if_else(parameter == p & relative_day >= zero_day_ext, 1, q95)
        )
    }
  }
}



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

hybrid_curve <- hybrid_curve %>%
  select(
    parameter, param_id, grid_id, tau, relative_day,
    mean, median, sd, q5, q95, panel_title,
    tau_wa_progress, tau_drc_conflict,
    q_west_africa_mean, q_drc_conflict_modulator,
    q_drc_conflict_low, q_drc_conflict_high,
    q_raw_mean, q_hybrid_scale,
    starts_with("theta_start"), starts_with("theta_end")
  )

# ------------------------------------------------------------
# BP-ready matrix
# ------------------------------------------------------------
q_wide <- q_hybrid_out %>%
  select(parameter, relative_day, tau, q = mean) %>%
  pivot_wider(names_from = parameter, values_from = q, names_prefix = "q_") %>%
  mutate(
    q_value = rowMeans(select(., starts_with("q_")), na.rm = TRUE),
    q_value = rescale_01(q_value)
  ) %>%
  select(relative_day, tau, q_value)

bp_matrix <- hybrid_curve %>%
  select(parameter, relative_day, tau, mean) %>%
  pivot_wider(names_from = parameter, values_from = mean) %>%
  mutate(
    prob_hosp = if ("p_hosp" %in% names(.)) p_hosp else prob_hosp,
    delay_hosp = delay_hosp,
    prob_unsafe_funeral_comm = if ("p_unsafe_funeral_comm" %in% names(.)) p_unsafe_funeral_comm else prob_unsafe_funeral_comm,
    prob_unsafe_funeral_hosp = if ("p_unsafe_funeral_hosp" %in% names(.)) p_unsafe_funeral_hosp else prob_unsafe_funeral_hosp,
    prob_unsafe_funeral_etu = 0,
    prop_etu = if ("p_ETU" %in% names(.)) p_ETU else prop_etu,
    ipc_helper = if ("latent_IPC" %in% names(.)) latent_IPC else ipc_helper,
    scenario = scenario_name,
    .before = 1
  ) %>%
  select(
    scenario,
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
  left_join(q_wide, by = c("relative_day", "tau"))

# Defensive clipping for probabilities/indices.
prob_cols <- c(
  "prob_hosp", "prob_unsafe_funeral_comm", "prob_unsafe_funeral_hosp",
  "prob_unsafe_funeral_etu", "prop_etu", "ipc_helper", "q_value"
)
bp_matrix <- bp_matrix %>%
  mutate(across(all_of(intersect(prob_cols, names(.))), clip01))

range_check <- bp_matrix %>%
  summarise(
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

if (any(bp_matrix[prob_cols] < -1e-8 | bp_matrix[prob_cols] > 1 + 1e-8, na.rm = TRUE)) {
  stop("At least one probability/index column falls outside [0,1].")
}

# ------------------------------------------------------------
# Plots
# ------------------------------------------------------------
plot_q_df <- q_hybrid_out %>%
  select(
    parameter, panel_title, relative_day,
    q_west_africa_mean, q_drc_conflict_modulator, mean
  ) %>%
  pivot_longer(
    cols = c(q_west_africa_mean, q_drc_conflict_modulator, mean),
    names_to = "curve",
    values_to = "q"
  ) %>%
  mutate(
    curve = recode(
      curve,
      q_west_africa_mean = "Latest endpoint-constrained West Africa Q_j(t) on original clock",
      q_drc_conflict_modulator = "DRC conflict modulator Q_conflict(t)",
      mean = "Extended West Africa-with-conflict++ Q_j(t)"
    ),
    curve = factor(
      curve,
      levels = c(
        "Latest endpoint-constrained West Africa Q_j(t) on original clock",
        "DRC conflict modulator Q_conflict(t)",
        "Extended West Africa-with-conflict++ Q_j(t)"
      )
    )
  )

p_q <- ggplot(plot_q_df, aes(x = relative_day, y = q, linetype = curve, linewidth = curve)) +
  geom_line(colour = "grey25", alpha = 0.85) +
  geom_line(
    data = plot_q_df %>% filter(curve == "Extended West Africa-with-conflict++ Q_j(t)"),
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
    title = "West Africa-with-conflict++ Q construction",
    subtitle = paste0(
      "Q_hybrid = Q_WA(original clock) × Q_DRC_conflict(extended clock), max-scaled; duration multiplier = ",
      signif(duration_multiplier, 3)
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

p_theta <- ggplot(hybrid_curve, aes(x = relative_day, y = mean)) +
  geom_ribbon(aes(ymin = q5, ymax = q95), fill = "#c7dcec", alpha = 0.9) +
  geom_line(linewidth = 0.9, colour = "#1f77b4") +
  facet_wrap(~ panel_title, scales = "free_y", ncol = 2) +
  theme_bw(base_size = 11) +
  labs(
    title = "Extended hypothetical West Africa conflict++ parameter trajectories, complete CI",
    subtitle = "Latest West Africa endpoints retained; CI propagates WA fit + DRC conflict-modulator uncertainty",
    x = "Extended West Africa relative outbreak day",
    y = NULL
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.minor = element_blank()
  )

print(p_q)
print(p_theta)

# ------------------------------------------------------------
# Save CSV/plot outputs
# ------------------------------------------------------------
write_csv(q_hybrid_out, paste0(out_prefix, "_Q_summaries.csv"))
write_csv(drc_q, paste0(out_prefix, "_DRC_conflict_modulator.csv"))
write_csv(q_range_check, paste0(out_prefix, "_Q_range_check.csv"))
write_csv(hybrid_curve, paste0(out_prefix, "_curve_summaries.csv"))
write_csv(bp_matrix, paste0(out_prefix, "_bp_input_matrix.csv"))
write_csv(range_check, paste0(out_prefix, "_bp_range_check.csv"))

ggsave(paste0(out_prefix, "_Qplot.png"), p_q, width = 12, height = 10, dpi = 180)
ggsave(paste0(out_prefix, "_parameter_curves.png"), p_theta, width = 12, height = 10, dpi = 180)
ggsave(paste0(out_prefix, "_Qplot.pdf"), p_q, width = 12, height = 10)
ggsave(paste0(out_prefix, "_parameter_curves.pdf"), p_theta, width = 12, height = 10)

# ------------------------------------------------------------
# Optional workbook update
# ------------------------------------------------------------
if (isTRUE(update_workbook)) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package 'openxlsx' required to update workbook. Install it or set update_workbook <- FALSE.")
  }
  check_file(matrix_workbook_path, "Scenario matrix workbook")

  old_sheet_names <- openxlsx::getSheetNames(matrix_workbook_path)
  sheets_to_replace <- c(new_sheet_name, "Combined_Long", "Notes_WA_Conflict")
  sheets_to_keep <- setdiff(old_sheet_names, sheets_to_replace)

  old_sheet_data <- list()
  for (s in sheets_to_keep) {
    old_sheet_data[[s]] <- openxlsx::readWorkbook(matrix_workbook_path, sheet = s)
  }

  add_df_sheet <- function(wb, sheet_name, df) {
    safe_name <- substr(sheet_name, 1, 31)
    openxlsx::addWorksheet(wb, safe_name)
    openxlsx::writeData(wb, sheet = safe_name, x = df)
    openxlsx::freezePane(wb, sheet = safe_name, firstActiveRow = 2)
    openxlsx::setColWidths(wb, sheet = safe_name, cols = seq_len(ncol(df)), widths = "auto")
    invisible(wb)
  }

  combined_list <- list()
  for (i in seq_along(existing_scenario_sheets)) {
    s <- existing_scenario_sheets[[i]]
    scenario_nm <- existing_scenario_names[[i]]

    if (s %in% names(old_sheet_data)) {
      df <- old_sheet_data[[s]]
    } else if (s %in% old_sheet_names) {
      df <- openxlsx::readWorkbook(matrix_workbook_path, sheet = s)
    } else {
      warning("Expected scenario sheet not found and skipped: ", s)
      next
    }

    combined_list[[scenario_nm]] <- df %>%
      select(-any_of("scenario")) %>%
      mutate(scenario = scenario_nm, .before = 1)
  }

  combined_list[[new_scenario_name]] <- bp_matrix %>%
    mutate(scenario = new_scenario_name, .before = 1)

  combined_long <- bind_rows(combined_list)

  notes <- tibble(
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
      paste0("West Africa response window extended by multiplier ", signif(duration_multiplier, 4), " (", duration_source, ")."),
      "Q_hybrid_j(t) = Q_WA_j(t on original West Africa response clock) * Q_DRC_conflict(t on extended conflict clock), then each parameter-specific Q_hybrid_j is max-scaled to span 0-1.",
      "Latest endpoint-constrained West Africa fitted start/end magnitudes are retained; Q curves alter timing/shape/duration, not magnitude endpoints.",
      "No new empirical points are added because this is a hypothetical counterfactual scenario.",
      "q_value is diagnostic only: row-wise mean of parameter-specific Q_hybrid_j values, rescaled to 0-1. BP inputs are the parameter columns.",
      basename(out_prefix)
    )
  )

  wb_out <- openxlsx::createWorkbook()
  for (s in sheets_to_keep) {
    add_df_sheet(wb_out, s, old_sheet_data[[s]])
  }
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
