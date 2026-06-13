# =============================================================================
# helper_functions_figure_1to4.R
#
# Shared helper functions for figures 1 to 4.
# Source this at the top of each figure-generation script.
#
# Simulation design:
#   - All runs are with OBV enabled.
#   - Each RDS file contains out$tdf and out$prevented_completed.
#
# Scenario reconstruction from a single arm's results:
#   "Without OBV" : tdf + prevented_completed  (use_prevented = TRUE)
#   "With OBV X%" : tdf only                   (use_prevented = FALSE)
#                   averted = deaths in prevented_completed
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(here)
library(patchwork)

# =============================================================================
# Constants
# =============================================================================
SCENARIO_LABELS <- c(
  WestAfrica = "West Africa (Worst)",
  DRC        = "DRC (Middle, PlusPlus)"
)
SCENARIO_COLORS <- c(
  WestAfrica = "#d95f02",
  DRC        = "#1b9e77"
)

OBV_EFFICACY_LEVELS <- c("obv_10", "obv_20", "obv_30", "obv_40", "obv_50",
                         "obv_60", "obv_70", "obv_80", "obv_90")
OBV_EFFICACY_LABELS <- c("10%", "20%", "30%", "40%", "50%",
                         "60%", "70%", "80%", "90%")
OBV_EFFICACY_VALUES <- c(obv_10 = 0.10, obv_20 = 0.20, obv_30 = 0.30,
                         obv_40 = 0.40, obv_50 = 0.50, obv_60 = 0.60,
                         obv_70 = 0.70, obv_80 = 0.80, obv_90 = 0.90)

COVERAGE_LEVELS <- c("full", "ramp_high", "ramp_low")
COVERAGE_LABELS <- c("Full (100%)", "Ramp high (0%->80%)", "Ramp low (0%->50%)")
COVERAGE_COLORS <- c(full = "#1a9641", ramp_high = "#fdae61", ramp_low = "#d7191c")

COVERAGE_LABELS <- c("Full (100%)", "Ramp high (0%->80%)", "Ramp low (0%->50%)")
COVERAGE_COLORS <- c(full = "#1a9641", ramp_high = "#fdae61", ramp_low = "#d7191c")

# Coverage curves
.make_clamped_spline <- function(t_knots, y_knots,
                                 deriv_start = 0, deriv_end = 0) {
  n  <- length(t_knots)
  h  <- diff(t_knots)
  
  rhs <- numeric(n)
  rhs[1] <- 3 * (y_knots[2] - y_knots[1]) / h[1] - 3 * deriv_start
  rhs[n] <- 3 * deriv_end - 3 * (y_knots[n] - y_knots[n - 1]) / h[n - 1]
  for (i in 2:(n - 1)) {
    rhs[i] <- 3 * ((y_knots[i + 1] - y_knots[i]) / h[i] -
                     (y_knots[i]     - y_knots[i - 1]) / h[i - 1])
  }
  
  diag_main <- c(2 * h[1], rep(0, n - 2), 2 * h[n - 1])
  diag_main[1] <- 2 * h[1]; diag_main[n] <- 2 * h[n - 1]
  for (i in 2:(n - 1)) diag_main[i] <- 2 * (h[i - 1] + h[i])
  
  c_vec <- h; d_vec <- rhs
  c_vec[1] <- c_vec[1] / diag_main[1]
  d_vec[1] <- d_vec[1] / diag_main[1]
  for (i in 2:n) {
    denom <- diag_main[i] - (if (i > 1) h[i - 1] else 0) * c_vec[i - 1]
    if (i < n) c_vec[i] <- h[i] / denom
    d_vec[i] <- (d_vec[i] - (if (i > 1) h[i - 1] else 0) * d_vec[i - 1]) / denom
  }
  
  M <- numeric(n); M[n] <- d_vec[n]
  for (i in (n - 1):1) M[i] <- d_vec[i] - c_vec[i] * M[i + 1]
  
  function(t_eval) {
    vapply(t_eval, function(t) {
      if (t <= t_knots[1]) return(y_knots[1])
      if (t >= t_knots[n]) return(y_knots[n])
      i  <- min(max(findInterval(t, t_knots, rightmost.closed = TRUE), 1L), n - 1L)
      hi <- h[i]; a <- (t_knots[i + 1] - t) / hi; b <- (t - t_knots[i]) / hi
      a * y_knots[i] + b * y_knots[i + 1] +
        ((a^3 - a) * M[i] + (b^3 - b) * M[i + 1]) * hi^2 / 6
    }, numeric(1))
  }
}

COVERAGE_SPECS <- list(
  full = list(fn = function(t) rep(1.0, length(t))),
  ramp_high = list(fn = local({
    spline_part <- .make_clamped_spline(
      t_knots = c(0, 90, 180), y_knots = c(0.0, 0.40, 0.80)
    )
    function(t) ifelse(t <= 180, spline_part(t), 0.80)
  })),
  ramp_low = list(fn = .make_clamped_spline(
    t_knots = c(0, 75, 365), y_knots = c(0.0, 0.0, 0.50)
  ))
)



# old version
# COVERAGE_SPECS <- list(
#   full = list(fn = function(t) rep(1.0, length(t))),
#   ramp_high = list(fn = .make_clamped_spline(
#     t_knots = c(0, 26 * 7, 52 * 7), y_knots = c(0.20, 0.50, 0.80)
#   )),
#   ramp_low = list(fn = .make_clamped_spline(
#     t_knots = c(0, 26 * 7, 52 * 7), y_knots = c(0.00, 0.25, 0.50)
#   ))
# )

COVERAGE_SPECS <- list(
  full = list(fn = function(t) rep(1.0, length(t))),
  ramp_high = list(fn = local({
    spline_part <- .make_clamped_spline(
      t_knots = c(0, 90, 180), y_knots = c(0.0, 0.40, 0.80)
    )
    function(t) ifelse(t <= 180, spline_part(t), 0.80)
  })),
  ramp_low = list(fn = .make_clamped_spline(
    t_knots = c(0, 75, 365), y_knots = c(0.0, 0.0, 0.50)
  ))
)

SCENARIO_X_MAX_DAYS <- c(WestAfrica = 60 * 7, DRC = 80 * 7)
DEFAULT_ARM <- "full_obv80"

# =============================================================================
# Null-coalescing operator
# =============================================================================
`%||%` <- function(a, b) if (!is.null(a)) a else b

# =============================================================================
# coverage_at_time
# =============================================================================
coverage_at_time <- function(t, coverage_spec) {
  pmin(pmax(coverage_spec$fn(t), 0), 1)
}

# =============================================================================
# compute_hcw_days_lost
#
# HCW days lost = time_symptom_onset_absolute -> duration for all infected HCW.
# =============================================================================
compute_hcw_days_lost <- function(cases, duration) {
  hcw <- cases[!is.na(cases$class) & cases$class == "HCW", ]
  if (nrow(hcw) == 0) return(0)
  days_lost <- pmax(duration - hcw$time_symptom_onset_absolute, 0)
  sum(days_lost, na.rm = TRUE)
}

# =============================================================================
# make_particle_df
# =============================================================================
make_particle_df <- function(run_df) {
  run_df %>%
    group_by(scenario, particle_id, arm) %>%
    summarise(
      n_hcw_deaths       = mean(n_hcw_deaths),
      hcw_days_lost      = sum(hcw_days_lost),
      prevented_hcw      = sum(prevented_hcw),
      counterfactual_hcw = sum(counterfactual_hcw),
      baseline_days_lost = sum(baseline_days_lost),
      .groups = "drop"
    ) %>%
    mutate(
      pct_hcw_deaths_averted = ifelse(
        counterfactual_hcw > 0,
        100 * prevented_hcw / counterfactual_hcw,
        NA_real_
      ),
      pct_days_lost_averted = ifelse(
        !is.na(baseline_days_lost) & baseline_days_lost > 0,
        100 * (baseline_days_lost - hcw_days_lost) / baseline_days_lost,
        NA_real_
      )
    )
}

# =============================================================================
# theme_fig
# =============================================================================
theme_fig <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(legend.position = "top")
}

# =============================================================================
# extract_weekly_ts  (Figure 1)
# =============================================================================
extract_weekly_ts <- function(arm_dir,
                              bin_width = 7,
                              n_workers = 10L,
                              base_dir  = here("outputs", "simulation")) {
  library(future)
  library(future.apply)
  
  dir_path <- file.path(base_dir, arm_dir)
  files    <- list.files(dir_path, pattern = "\\.rds$", full.names = TRUE)
  if (length(files) == 0) stop("No RDS files found in: ", dir_path)
  message(sprintf("Extracting weekly ts from %d files in arm '%s'...",
                  length(files), arm_dir))
  
  plan(multisession, workers = min(n_workers, future::availableCores()))
  on.exit(plan(sequential), add = TRUE)
  
  results <- future_lapply(seq_along(files), function(i) {
    if (i %% 10 == 0)
      message(sprintf("  Processing file %d / %d...", i, length(files)))
    
    f <- files[[i]]
    x <- readRDS(f)
    
    fname <- tools::file_path_sans_ext(basename(f))
    parts <- regmatches(fname, regexec("^(.+)_p(\\d+)_r(\\d+)$", fname))[[1]]
    sc    <- if (length(parts) == 4) parts[2] else NA_character_
    pid   <- if (length(parts) == 4) as.integer(parts[3]) else NA_integer_
    rep   <- if (length(parts) == 4) as.integer(parts[4]) else NA_integer_
    
    cap    <- SCENARIO_X_MAX_DAYS[sc]
    breaks <- seq(0, ceiling(cap / bin_width) * bin_width, by = bin_width)
    mids   <- breaks[-length(breaks)] + bin_width / 2
    
    .bin <- function(times) {
      t_clip <- times[!is.na(times) & times >= 0 & times <= max(breaks)]
      hist(t_clip, breaks = breaks, plot = FALSE)$counts
    }
    
    .extract <- function(cases, arm_label) {
      is_hcw <- !is.na(cases$class) & cases$class == "HCW"
      died   <- !is.na(cases$outcome) & cases$outcome
      metrics <- list(
        deaths               = cases$time_outcome_absolute[died],
        infections           = cases$time_infection_absolute,
        hcw_deaths_incidence = cases$time_outcome_absolute[died & is_hcw],
        hcw_deaths           = cases$time_outcome_absolute[died & is_hcw]
      )
      # All returned as incidence; cumsum for hcw_deaths applied after aggregation
      do.call(rbind, lapply(names(metrics), function(m) {
        data.frame(scenario = sc, particle_id = pid, rep = rep,
                   arm = arm_label, week = mids, metric = m,
                   value = .bin(metrics[[m]]), stringsAsFactors = FALSE)
      }))
    }
    
    tdf       <- x$tdf
    prevented <- x$prevented_completed
    if (!is.null(prevented) && nrow(prevented) > 0) {
      missing_cols <- setdiff(names(tdf), names(prevented))
      for (col in missing_cols) prevented[[col]] <- NA
      prevented  <- prevented[, names(tdf), drop = FALSE]
      cases_base <- rbind(tdf, prevented)
    } else {
      cases_base <- tdf
    }
    
    rbind(.extract(cases_base, "baseline"), .extract(tdf, "obv"))
  }, future.packages = c("here"), future.seed = TRUE)
  
  do.call(rbind, results)
}

# =============================================================================
# extract_run_summary  (Figures 2, 3, 4)
# =============================================================================
extract_run_summary <- function(arm_dir,
                                arm_label = arm_dir,
                                n_workers = 10L,
                                obv_return = FALSE,
                                base_dir  = here("outputs", "simulation")) {
  library(future)
  library(future.apply)
  
  dir_path <- file.path(base_dir, arm_dir)
  files    <- list.files(dir_path, pattern = "\\.rds$", full.names = TRUE)
  if (length(files) == 0) stop("No RDS files found in: ", dir_path)
  message(sprintf("Extracting run summaries from %d files in arm '%s'...",
                  length(files), arm_dir))
  
  plan(multisession, workers = min(n_workers, future::availableCores()))
  on.exit(plan(sequential), add = TRUE)
  
  results <- future_lapply(seq_along(files), function(i) {
    if (i %% 10 == 0)
      message(sprintf("  Processing file %d / %d...", i, length(files)))
    
    f <- files[[i]]
    x <- readRDS(f)
    
    fname <- tools::file_path_sans_ext(basename(f))
    parts <- regmatches(fname, regexec("^(.+)_p(\\d+)_r(\\d+)$", fname))[[1]]
    sc    <- if (length(parts) == 4) parts[2] else NA_character_
    pid   <- if (length(parts) == 4) as.integer(parts[3]) else NA_integer_
    rep   <- if (length(parts) == 4) as.integer(parts[4]) else NA_integer_
    
    tdf       <- x$tdf
    prevented <- x$prevented_completed
    duration  <- max(tdf$time_outcome_absolute, na.rm = TRUE)
    
    if (!is.null(prevented) && nrow(prevented) > 0) {
      missing_cols <- setdiff(names(tdf), names(prevented))
      for (col in missing_cols) prevented[[col]] <- NA
      prevented  <- prevented[, names(tdf), drop = FALSE]
      cases_base <- rbind(tdf, prevented)
    } else {
      cases_base <- tdf
    }
    
    .hcw_deaths <- function(cases) {
      is_hcw <- !is.na(cases$class) & cases$class == "HCW"
      died   <- !is.na(cases$outcome) & cases$outcome
      sum(died & is_hcw)
    }
    .days_lost <- function(cases, is_prevented = FALSE) {
      hcw <- cases[!is.na(cases$class) & cases$class == "HCW", ]
      if (nrow(hcw) == 0) return(0)
      if (obv_return) {
        # OBV recipients who survived return to work at time_outcome_absolute.
        # prevented_completed HCW are all OBV recipients who survived.
        # tdf HCW with obv_pep_received=TRUE and survived also return early.
        died         <- !is.na(hcw$outcome) & hcw$outcome
        obv_recv     <- is_prevented | (!is.na(hcw$obv_pep_received) & hcw$obv_pep_received)
        early_return <- obv_recv & !died
        absence_end  <- ifelse(early_return, hcw$time_outcome_absolute, duration)
      } else {
        # Default: all infected HCW absent until simulation end
        absence_end <- rep(duration, nrow(hcw))
      }
      sum(pmax(absence_end - hcw$time_symptom_onset_absolute, 0), na.rm = TRUE)
    }
    
    n_base    <- .hcw_deaths(cases_base)
    n_obv     <- .hcw_deaths(tdf)
    # Baseline: prevented HCW are NOT OBV recipients (counterfactual world)
    days_base <- .days_lost(cases_base, is_prevented = FALSE)
    # OBV arm: tdf as-is; prevented HCW counted separately as early returners
    days_obv  <- .days_lost(tdf, is_prevented = FALSE) +
      if (obv_return && !is.null(x$prevented_completed) &&
          nrow(x$prevented_completed) > 0) {
        prev_hcw <- x$prevented_completed[
          !is.na(x$prevented_completed$class) &
            x$prevented_completed$class == "HCW", ]
        if (nrow(prev_hcw) > 0)
          sum(pmax(prev_hcw$time_outcome_absolute -
                     prev_hcw$time_symptom_onset_absolute, 0), na.rm = TRUE)
        else 0
      } else 0
    
    # Store baseline and OBV values in a single row per run.
    # This avoids duplicate baseline rows when multiple arms are combined.
    data.frame(
      scenario           = sc,
      particle_id        = pid,
      rep                = rep,
      arm                = arm_label,
      n_hcw_deaths       = n_obv,
      hcw_days_lost      = days_obv,
      counterfactual_hcw = n_base,
      prevented_hcw      = n_base - n_obv,
      baseline_days_lost = days_base,
      stringsAsFactors   = FALSE
    )
  }, future.packages = c("here"), future.seed = TRUE)
  
  do.call(rbind, results)
}

# =============================================================================
# extract_dose_summary  (Figure 5)
# =============================================================================
extract_dose_summary <- function(arm_dir,
                                 eff_name,
                                 n_workers = 10L,
                                 base_dir  = here("outputs", "simulation")) {
  library(future)
  library(future.apply)
  
  efficacy <- OBV_EFFICACY_VALUES[[eff_name]]
  dir_path <- file.path(base_dir, arm_dir)
  files    <- list.files(dir_path, pattern = "\\.rds$", full.names = TRUE)
  if (length(files) == 0) stop("No RDS files found in: ", dir_path)
  message(sprintf("Extracting dose summaries from %d files in arm '%s'...",
                  length(files), arm_dir))
  
  plan(multisession, workers = min(n_workers, future::availableCores()))
  on.exit(plan(sequential), add = TRUE)
  
  results <- future_lapply(seq_along(files), function(i) {
    if (i %% 10 == 0)
      message(sprintf("  Processing file %d / %d...", i, length(files)))
    
    f <- files[[i]]
    x <- readRDS(f)
    
    fname <- tools::file_path_sans_ext(basename(f))
    parts <- regmatches(fname, regexec("^(.+)_p(\\d+)_r(\\d+)$", fname))[[1]]
    sc    <- if (length(parts) == 4) parts[2] else NA_character_
    pid   <- if (length(parts) == 4) as.integer(parts[3]) else NA_integer_
    rep   <- if (length(parts) == 4) as.integer(parts[4]) else NA_integer_
    
    tdf       <- x$tdf
    prevented <- x$prevented_completed
    
    is_hcw_tdf  <- !is.na(tdf$class) & tdf$class == "HCW"
    is_hcw_prev <- if (!is.null(prevented) && nrow(prevented) > 0)
      !is.na(prevented$class) & prevented$class == "HCW"
    else logical(0)
    
    data.frame(
      scenario        = sc,
      particle_id     = pid,
      rep             = rep,
      eff_name        = eff_name,
      efficacy        = efficacy,
      doses_B         = sum(is_hcw_tdf) + sum(is_hcw_prev),
      n_prevented_hcw = sum(is_hcw_prev),
      stringsAsFactors = FALSE
    )
  }, future.packages = c("here"), future.seed = TRUE)
  
  do.call(rbind, results)
}

# =============================================================================
# save_figure_data
# =============================================================================
save_figure_data <- function(df, filename,
                             out_dir = here("output_figgen")) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(out_dir, filename)
  write.csv(df, path, row.names = FALSE)
  message(sprintf("Saved figure data: %s", path))
  invisible(path)
}