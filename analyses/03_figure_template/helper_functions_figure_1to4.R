# =============================================================================
# helper_functions_figure_1to4.R
#
# Shared helper functions for figures 1 to 4.
# Source this at the top of each figure-generation script.
#
# Key design: OBV efficacy and coverage are applied post-hoc to baseline
# simulation output (tdf). No separate OBV simulation runs are needed.
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

# =============================================================================
# load_results
#
# Loads all RDS files from the baseline simulation output directory.
# Each file is a single simulation run (one particle x rep combination).
# =============================================================================
load_results <- function(base_dir = here("outputs", "simulation_baseline")) {
  files <- list.files(base_dir, pattern = "\\.rds$", full.names = TRUE)
  if (length(files) == 0) stop("No RDS files found in: ", base_dir)
  message(sprintf("Loading %d baseline files...", length(files)))
  lapply(files, readRDS)
}

# =============================================================================
# coverage_at_time
#
# Interpolates a piecewise-linear coverage curve at time t (in days).
# Values outside the defined range are clamped to the boundary values.
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
#
# Applies OBV intervention post-hoc to a single run's HCW case data.
# OBV is modelled as affecting individual outcomes only; transmission is
# unchanged (so the counterfactual infection count equals baseline).
#
# For each infected HCW at time t:
#   1. Coverage probability is read from the coverage curve at t.
#   2. Bernoulli draw: does this HCW receive OBV?
#   3. Among HCW who both died and received OBV: Bernoulli draw with
#      probability = efficacy determines whether the death is prevented.
#
# Arguments:
#   cases         - data frame of individual cases (x$tdf from a simulation run)
#   efficacy      - scalar in [0, 1]; probability of preventing death given OBV
#   coverage_spec - list(times, values) defining the coverage curve
#   seed          - optional RNG seed for reproducibility
#
# Returns a list:
#   prevented_hcw      - number of HCW deaths prevented
#   counterfactual_hcw - total HCW deaths without OBV (equals baseline)
#   obv_received       - logical vector (length = n HCW rows) for OBV receipt
#   prevented_flag     - logical vector (length = n HCW rows) for prevented
#                        deaths; pass this directly to compute_hcw_days_lost
#                        to keep prevented individuals consistent with the count
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
    # Do NOT re-draw runif for prevented_flag elsewhere; re-drawing produces a
    # different set of prevented individuals, making days-lost inconsistent with
    # the prevented count.
    prevented_flag     = prevented_flag
  )
}

# =============================================================================
# compute_hcw_days_lost
#
# Computes total HCW work-days lost for a single run under three cases:
#   - Died (after OBV if applicable): absent from absence_start to simulation end
#   - Survived, no OBV: absent until min(recovery + 6 months, simulation end)
#   - Survived, received OBV: earlier return, absent until recovery + OBV_RETURN_TO_WORK_DAYS
#
# Absence start is hospitalisation time if recorded, otherwise symptom onset.
#
# Arguments:
#   cases           - data frame of individual cases (x$tdf)
#   duration        - simulation duration in days (x$duration)
#   obv_received    - logical vector aligned to HCW rows (from apply_obv_posthoc),
#                     or NULL for baseline (treated as all FALSE)
#   prevented       - logical vector of prevented deaths (from apply_obv_posthoc),
#                     or NULL for baseline (treated as all FALSE)
#   obv_return_days - days post-recovery before OBV-treated HCW return to work
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

  absence_end <- rep(NA_real_, nrow(hcw))

  absence_end[died_after_obv] <- duration

  surv_no_obv <- !died_after_obv & !received
  absence_end[surv_no_obv] <- pmin(
    hcw$time_outcome_absolute[surv_no_obv] + 6 * 30,
    duration
  )

  # Covers both: (a) HCW who survived naturally and received OBV, and
  # (b) HCW who would have died but were saved by OBV.
  surv_obv <- !died_after_obv & received
  absence_end[surv_obv] <- hcw$time_outcome_absolute[surv_obv] + obv_return_days

  days_lost <- pmax(absence_end - absence_start, 0)
  sum(days_lost, na.rm = TRUE)
}

# =============================================================================
# build_run_df_obv
#
# Builds a per-run summary data frame for one OBV arm. For baseline,
# uses raw simulation counts. For OBV arms, applies apply_obv_posthoc
# with a deterministic per-run seed for reproducibility.
#
# Arguments:
#   results       - list of simulation run objects from load_results()
#   efficacy_name - one of OBV_EFFICACY_LEVELS, or "baseline"
#   coverage_name - one of COVERAGE_LEVELS, or NULL (required unless baseline)
#   seed_offset   - integer added to the per-run seed (use to vary across calls)
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
#
# Aggregates run-level data to the particle level and computes % averted.
#
# Uses a burden-weighted approach: sum prevented and counterfactual deaths
# across reps before dividing, so that particles with more HCW deaths
# contribute proportionally more to the estimate. This guarantees
# pct_hcw_deaths_averted is in [0, 100] with no negative effectiveness.
#
# pct_days_lost_averted is computed relative to the baseline days-lost for
# the same particle.
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
# Builds a time series of quantiles across posterior particles for a given
# metric, binned by bin_width days.
#
# metric options:
#   "hcw_deaths"          - cumulative HCW deaths over time
#   "hcw_deaths_incidence"- incident HCW deaths per bin (not cumulative)
#   "hcw_infections"      - cumulative HCW infections over time
#   "infections"          - incident infections (entire population, not cumulative)
#   "deaths"              - incident deaths (entire population, not cumulative)
#
# For OBV arms, prevented deaths are excluded from the death time series.
# The output mid-point times are in the same unit as bin_width (days by default;
# divide by 7 in the calling script when weekly x-axes are needed).
#
# Returns a data frame with columns:
#   scenario, arm, week (bin midpoint), q025, q25, q50, q75, q975
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

  max_day <- max(sapply(results, function(x) x$duration), na.rm = TRUE)
  breaks  <- seq(0, ceiling(max_day / bin_width) * bin_width, by = bin_width)
  mids    <- breaks[-length(breaks)] + bin_width / 2

  rows <- do.call(rbind, lapply(results, function(x) {
    cases  <- x$tdf
    is_hcw <- cases$class == "HCW"
    died   <- !is.na(cases$outcome) & cases$outcome

    if (!is_baseline) {
      run_seed       <- seed_offset + x$particle_id * 1000L + x$rep
      obv            <- apply_obv_posthoc(cases, efficacy, cov_spec, seed = run_seed)
      prevented_full <- logical(nrow(cases))
      prevented_full[which(is_hcw)] <- obv$prevented_flag
      died <- died & !prevented_full
    }

    times <- switch(metric,
      hcw_deaths           = cases$time_outcome_absolute[died & is_hcw],
      hcw_deaths_incidence = cases$time_outcome_absolute[died & is_hcw],
      hcw_infections       = cases$time_infection_absolute[is_hcw],
      deaths               = cases$time_outcome_absolute[died],
      infections           = cases$time_infection_absolute
    )
    counts <- hist(times[!is.na(times)], breaks = breaks, plot = FALSE)$counts

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
#
# Shared ggplot theme for all figures.
# =============================================================================
theme_fig <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      legend.position = "top"
    )
}
