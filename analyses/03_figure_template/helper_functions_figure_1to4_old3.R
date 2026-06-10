# =============================================================================
# helper_functions_figure_1to4.R
#
# Shared helper functions for figures 1 to 4.
# Source this at the top of each figure-generation script.
#
# Key design: OBV efficacy and coverage are applied post-hoc to baseline
# simulation output (tdf). No separate OBV simulation runs are needed.
#
# Takeoff filter: only runs with n_deaths >= TAKEOFF_DEATH_THRESHOLD are used,
# consistent with the ABC fitting scheme (aggregate_decoupled / run_abc_particle_decoupled).
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

OBV_EFFICACY_LEVELS <- c("obv_50", "obv_60", "obv_70", "obv_80", "obv_90")
OBV_EFFICACY_LABELS <- c("50%", "60%", "70%", "80%", "90%")
OBV_EFFICACY_VALUES <- c(obv_50 = 0.50, obv_60 = 0.60, obv_70 = 0.70,
                         obv_80 = 0.80, obv_90 = 0.90)

COVERAGE_LEVELS <- c("full", "ramp_high", "ramp_low")
COVERAGE_LABELS <- c("Full (100%)", "Ramp high (20%->100%)", "Ramp low (20%->50%)")
COVERAGE_COLORS <- c(full = "#1a9641", ramp_high = "#fdae61", ramp_low = "#d7191c")

# Piecewise-linear coverage curves: breakpoints in days, values in [0, 1]
COVERAGE_SPECS <- list(
  full      = list(times = c(0, 90),           values = c(1.00, 1.00)),
  ramp_high = list(times = c(0, 30, 60, 90),   values = c(0.20, 0.47, 0.73, 1.00)),
  ramp_low  = list(times = c(0, 30, 60, 90),   values = c(0.20, 0.30, 0.40, 0.50))
)

# Days post-recovery before an OBV-treated HCW can return to work
OBV_RETURN_TO_WORK_DAYS <- 7

# Minimum deaths for a run to be considered "taken off" --
# must match the ABC fitting scheme (TAKEOFF_DEATH_THRESHOLD in calibration scripts)
TAKEOFF_DEATH_THRESHOLD <- 100L

# x-axis caps per scenario (in days; divide by 7 for weeks in plotting scripts)
SCENARIO_X_MAX_DAYS <- c(WestAfrica = 60 * 7, DRC = 80 * 7)

# =============================================================================
# load_results
#
# Loads all RDS files from the baseline simulation output directory and
# filters to taken-off runs only (n_deaths >= TAKEOFF_DEATH_THRESHOLD),
# consistent with the ABC fitting scheme.
# =============================================================================
load_results <- function(base_dir = here("outputs", "simulation_baseline"),
                         takeoff_threshold = TAKEOFF_DEATH_THRESHOLD) {
  files <- list.files(base_dir, pattern = "\\.rds$", full.names = TRUE)
  if (length(files) == 0) stop("No RDS files found in: ", base_dir)
  message(sprintf("Loading %d baseline files...", length(files)))
  results <- lapply(files, readRDS)
  
  # Filter to taken-off runs only (consistent with ABC fitting)
  n_before  <- length(results)
  results   <- Filter(function(x) x$n_deaths >= takeoff_threshold, results)
  n_after   <- length(results)
  message(sprintf("  Takeoff filter (n_deaths >= %d): %d / %d runs retained (%.0f%%)",
                  takeoff_threshold, n_after, n_before,
                  100 * n_after / max(n_before, 1)))
  results
}

# =============================================================================
# coverage_at_time
# =============================================================================
coverage_at_time <- function(t, times, values) {
  if (t <= times[1])              return(values[1])
  if (t >= times[length(times)])  return(values[length(values)])
  i    <- findInterval(t, times)
  frac <- (t - times[i]) / (times[i + 1L] - times[i])
  values[i] + frac * (values[i + 1L] - values[i])
}

# =============================================================================
# apply_obv_posthoc
# =============================================================================
apply_obv_posthoc <- function(cases, efficacy, coverage_spec, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  hcw <- cases[!is.na(cases$class) & cases$class == "HCW", ]
  n   <- nrow(hcw)
  
  if (n == 0) {
    return(list(prevented_hcw = 0L, counterfactual_hcw = 0L,
                obv_received = logical(0), prevented_flag = logical(0)))
  }
  
  died <- !is.na(hcw$outcome) & hcw$outcome
  
  cov_prob <- vapply(hcw$time_infection_absolute, function(t) {
    coverage_at_time(t, coverage_spec$times, coverage_spec$values)
  }, numeric(1))
  
  obv_received   <- runif(n) < cov_prob
  prevented_flag <- died & obv_received & (runif(n) < efficacy)
  
  list(
    prevented_hcw      = sum(prevented_flag),
    counterfactual_hcw = sum(died),
    obv_received       = obv_received,
    prevented_flag     = prevented_flag
  )
}

# =============================================================================
# compute_hcw_days_lost
# =============================================================================
compute_hcw_days_lost <- function(cases, duration,
                                  obv_received    = NULL,
                                  prevented       = NULL,
                                  obv_return_days = OBV_RETURN_TO_WORK_DAYS) {
  hcw <- cases[!is.na(cases$class) & cases$class == "HCW", ]
  if (nrow(hcw) == 0) return(0)
  
  absence_start <- ifelse(
    !is.na(hcw$time_hospitalisation_absolute),
    hcw$time_hospitalisation_absolute,
    hcw$time_symptom_onset_absolute
  )
  
  died            <- !is.na(hcw$outcome) & hcw$outcome
  received        <- if (!is.null(obv_received)) obv_received else rep(FALSE, nrow(hcw))
  prevented_death <- if (!is.null(prevented))    prevented    else rep(FALSE, nrow(hcw))
  
  died_after_obv <- died & !prevented_death
  absence_end    <- rep(NA_real_, nrow(hcw))
  
  absence_end[died_after_obv] <- duration
  
  surv_no_obv <- !died_after_obv & !received
  absence_end[surv_no_obv] <- pmin(
    hcw$time_outcome_absolute[surv_no_obv] + 6 * 30,
    duration
  )
  
  surv_obv <- !died_after_obv & received
  absence_end[surv_obv] <- hcw$time_outcome_absolute[surv_obv] + obv_return_days
  
  days_lost <- pmax(absence_end - absence_start, 0)
  sum(days_lost, na.rm = TRUE)
}

# =============================================================================
# build_run_df_obv
# =============================================================================
build_run_df_obv <- function(results,
                             efficacy_name = "baseline",
                             coverage_name = NULL,
                             seed_offset   = 0L) {
  is_baseline <- efficacy_name == "baseline"
  efficacy    <- if (is_baseline) NA_real_ else OBV_EFFICACY_VALUES[[efficacy_name]]
  cov_spec    <- if (is_baseline || is.null(coverage_name)) NULL
  else COVERAGE_SPECS[[coverage_name]]
  
  do.call(rbind, lapply(results, function(x) {
    if (is_baseline) {
      days_lost     <- compute_hcw_days_lost(x$tdf, x$duration)
      prevented_hcw <- 0L
    } else {
      run_seed  <- seed_offset + x$particle_id * 1000L + x$rep
      obv       <- apply_obv_posthoc(x$tdf, efficacy, cov_spec, seed = run_seed)
      days_lost <- compute_hcw_days_lost(x$tdf, x$duration,
                                         obv_received = obv$obv_received,
                                         prevented    = obv$prevented_flag)
      prevented_hcw <- obv$prevented_hcw
    }
    
    data.frame(
      scenario           = x$scenario,
      particle_id        = x$particle_id,
      rep                = x$rep,
      arm                = efficacy_name,
      coverage_scenario  = if (is_baseline) "baseline" else coverage_name,
      n_infections       = x$n_infections,
      n_hcw_deaths       = x$n_hcw_deaths - prevented_hcw,
      counterfactual_hcw = x$n_hcw_deaths,
      prevented_hcw      = prevented_hcw,
      hcw_days_lost      = days_lost,
      stringsAsFactors   = FALSE
    )
  }))
}

# =============================================================================
# make_particle_df
# =============================================================================
make_particle_df <- function(run_df) {
  pdf <- run_df %>%
    group_by(scenario, particle_id, arm, coverage_scenario) %>%
    summarise(
      n_infections       = mean(n_infections),
      n_hcw_deaths       = mean(n_hcw_deaths),
      hcw_days_lost      = mean(hcw_days_lost),
      prevented_hcw      = sum(prevented_hcw),
      counterfactual_hcw = sum(counterfactual_hcw),
      .groups = "drop"
    ) %>%
    mutate(
      pct_hcw_deaths_averted = ifelse(
        arm != "baseline" & counterfactual_hcw > 0,
        100 * prevented_hcw / counterfactual_hcw,
        NA_real_
      )
    )
  
  base_days <- pdf %>%
    filter(arm == "baseline") %>%
    select(scenario, particle_id, baseline_days_lost = hcw_days_lost)
  
  pdf %>%
    left_join(base_days, by = c("scenario", "particle_id")) %>%
    mutate(
      pct_days_lost_averted = ifelse(
        arm != "baseline" & !is.na(baseline_days_lost) & baseline_days_lost > 0,
        100 * (baseline_days_lost - hcw_days_lost) / baseline_days_lost,
        NA_real_
      )
    )
}

# =============================================================================
# build_weekly_ts
#
# Builds a time series of quantiles across posterior particles.
#
# Key changes from previous version:
#   - Only taken-off runs are included (filtered in load_results()).
#   - max_day is capped at SCENARIO_X_MAX_DAYS per scenario to avoid the
#     sparse-trailing-zero problem from short-lived runs stretching the grid.
#   - Quantiles are computed across particle means (mean over reps first),
#     consistent with ABC averaging over replicates.
#
# metric options:
#   "hcw_deaths"          - cumulative HCW deaths
#   "hcw_deaths_incidence"- incident HCW deaths per bin
#   "hcw_infections"      - cumulative HCW infections
#   "infections"          - incident infections (all population)
#   "deaths"              - incident deaths (all population)
# =============================================================================
build_weekly_ts <- function(results,
                            metric        = c("hcw_deaths", "hcw_infections",
                                              "infections", "deaths",
                                              "hcw_deaths_incidence"),
                            bin_width     = 28,
                            efficacy_name = "baseline",
                            coverage_name = NULL,
                            seed_offset   = 0L) {
  metric <- match.arg(metric)
  
  is_baseline <- efficacy_name == "baseline"
  efficacy    <- if (is_baseline) NA_real_ else OBV_EFFICACY_VALUES[[efficacy_name]]
  cov_spec    <- if (is_baseline || is.null(coverage_name)) NULL
  else COVERAGE_SPECS[[coverage_name]]
  
  # Build per-scenario breaks capped at SCENARIO_X_MAX_DAYS
  scenarios_present <- unique(sapply(results, function(x) x$scenario))
  breaks_by_sc <- lapply(setNames(scenarios_present, scenarios_present), function(sc) {
    cap <- SCENARIO_X_MAX_DAYS[sc]
    seq(0, ceiling(cap / bin_width) * bin_width, by = bin_width)
  })
  
  rows <- do.call(rbind, lapply(results, function(x) {
    cases  <- x$tdf
    is_hcw <- cases$class == "HCW"
    breaks <- breaks_by_sc[[x$scenario]]
    died   <- !is.na(cases$outcome) & cases$outcome
    
    if (!is_baseline) {
      run_seed       <- seed_offset + x$particle_id * 1000L + x$rep
      obv            <- apply_obv_posthoc(cases, efficacy, cov_spec, seed = run_seed)
      prevented_full <- logical(nrow(cases))
      prevented_full[which(is_hcw)] <- obv$prevented_flag
      died <- died & !prevented_full
    }
    
    mids <- breaks[-length(breaks)] + bin_width / 2
    
    times <- switch(metric,
                    hcw_deaths           = cases$time_outcome_absolute[died & is_hcw],
                    hcw_deaths_incidence = cases$time_outcome_absolute[died & is_hcw],
                    hcw_infections       = cases$time_infection_absolute[is_hcw],
                    deaths               = cases$time_outcome_absolute[died],
                    infections           = cases$time_infection_absolute
    )
    
    # Clip times to the scenario's break range before binning
    t_clipped <- times[!is.na(times) & times >= 0 & times <= max(breaks)]
    counts    <- hist(t_clipped, breaks = breaks, plot = FALSE)$counts
    
    data.frame(
      scenario    = x$scenario,
      arm         = efficacy_name,
      particle_id = x$particle_id,
      rep         = x$rep,
      week        = mids,
      value       = counts,
      stringsAsFactors = FALSE
    )
  }))
  
  is_incidence <- metric %in% c("infections", "deaths", "hcw_deaths_incidence")
  
  # Mean over reps per particle, then quantiles over particles
  rows %>%
    group_by(scenario, arm, particle_id, week) %>%
    summarise(value = mean(value), .groups = "drop") %>%
    arrange(scenario, arm, particle_id, week) %>%
    group_by(scenario, arm, particle_id) %>%
    mutate(value = if (is_incidence) value else cumsum(value)) %>%
    ungroup() %>%
    group_by(scenario, arm, week) %>%
    summarise(
      q025 = quantile(value, 0.025),
      q25  = quantile(value, 0.25),
      q50  = quantile(value, 0.50),
      q75  = quantile(value, 0.75),
      q975 = quantile(value, 0.975),
      .groups = "drop"
    )
}

# =============================================================================
# theme_fig
# =============================================================================
theme_fig <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      legend.position = "top"
    )
}