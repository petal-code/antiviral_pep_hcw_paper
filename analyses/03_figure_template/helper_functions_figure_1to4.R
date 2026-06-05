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

# Coverage curve specs: piecewise-linear breakpoints (time in days, value in [0,1])
COVERAGE_SPECS <- list(
  full      = list(times = c(0, 90),           values = c(1.00, 1.00)),
  ramp_high = list(times = c(0, 30, 60, 90),   values = c(0.20, 0.47, 0.73, 1.00)),
  ramp_low  = list(times = c(0, 30, 60, 90),   values = c(0.20, 0.30, 0.40, 0.50))
)

# Arm colors — baseline uses scenario color, OBV arms use fixed gradient
ARM_COLORS_OBV <- c(
  obv_50 = "#fee090",
  obv_60 = "#fdae61",
  obv_70 = "#f46d43",
  obv_80 = "#d73027",
  obv_90 = "#a50026"
)
ARM_LABELS <- c(
  baseline = "Baseline",
  obv_50   = "OBV 50%",
  obv_60   = "OBV 60%",
  obv_70   = "OBV 70%",
  obv_80   = "OBV 80%",
  obv_90   = "OBV 90%"
)

get_arm_colors <- function(sc, arms) {
  colors <- c(baseline = unname(SCENARIO_COLORS[sc]), ARM_COLORS_OBV)
  colors[arms]
}

# Days post-recovery before OBV-treated HCW can return to work
OBV_RETURN_TO_WORK_DAYS <- 7

# =============================================================================
# Load all RDS files from the baseline simulation directory
# =============================================================================
load_results <- function(base_dir = here("outputs", "simulation_baseline")) {
  files <- list.files(base_dir, pattern = "\\.rds$", full.names = TRUE)
  if (length(files) == 0) stop("No RDS files found in: ", base_dir)
  message(sprintf("Loading %d baseline files...", length(files)))
  lapply(files, readRDS)
}

# =============================================================================
# Interpolate a piecewise-linear coverage curve at a given time point
# =============================================================================
coverage_at_time <- function(t, times, values) {
  # Clamp to the defined range
  if (t <= times[1])  return(values[1])
  if (t >= times[length(times)]) return(values[length(values)])
  # Find the enclosing interval and linearly interpolate
  i <- findInterval(t, times)
  frac <- (t - times[i]) / (times[i + 1L] - times[i])
  values[i] + frac * (values[i + 1L] - values[i])
}

# =============================================================================
# Apply OBV post-hoc to a single run's HCW cases
#
# For each HCW who was infected at time t:
#   1. Draw coverage probability from the coverage curve at t
#   2. Bernoulli sample: does this HCW receive OBV?
#   3. Among those who receive OBV AND died: Bernoulli sample with efficacy
#      to determine whether the death was prevented
#
# Returns a list:
#   prevented_hcw  - number of deaths prevented by OBV in this run
#   counterfactual_hcw - total HCW deaths that would have occurred without OBV
#                        (= original deaths; same as baseline since OBV does
#                         not change transmission, only individual outcomes)
#   obv_received   - logical vector over HCW rows indicating OBV receipt
# =============================================================================
apply_obv_posthoc <- function(cases, efficacy, coverage_spec, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  hcw  <- cases[!is.na(cases$class) & cases$class == "HCW", ]
  n    <- nrow(hcw)
  
  if (n == 0) {
    return(list(prevented_hcw = 0L, counterfactual_hcw = 0L,
                obv_received = logical(0)))
  }
  
  died <- !is.na(hcw$outcome) & hcw$outcome
  
  # Coverage probability at infection time for each HCW
  cov_prob <- vapply(hcw$time_infection_absolute, function(t) {
    coverage_at_time(t, coverage_spec$times, coverage_spec$values)
  }, numeric(1))
  
  # Each HCW independently receives OBV with probability = coverage
  obv_received <- runif(n) < cov_prob
  
  # Among HCW who died AND received OBV: each death is prevented with prob = efficacy
  prevented_flag <- died & obv_received & (runif(n) < efficacy)
  
  list(
    prevented_hcw      = sum(prevented_flag),
    counterfactual_hcw = sum(died),   # unchanged from baseline
    obv_received       = obv_received,
    prevented_flag     = prevented_flag   # full logical vector; use this for days-lost calc
  )
}

# =============================================================================
# Compute HCW days lost for a single run, with optional post-hoc OBV receipt
#
# obv_received: logical vector aligned to HCW rows of cases (from apply_obv_posthoc),
#               or NULL for no-OBV (baseline) computation.
# prevented:    logical vector -- which HCW deaths were prevented by OBV.
#               If NULL, no deaths are prevented (baseline).
# =============================================================================
compute_hcw_days_lost <- function(cases, duration,
                                  obv_received = NULL,
                                  prevented    = NULL,
                                  obv_return_days = OBV_RETURN_TO_WORK_DAYS) {
  hcw <- cases[!is.na(cases$class) & cases$class == "HCW", ]
  if (nrow(hcw) == 0) return(0)
  
  # Start of absence: hospitalisation time if available, else symptom onset
  absence_start <- ifelse(
    !is.na(hcw$time_hospitalisation_absolute),
    hcw$time_hospitalisation_absolute,
    hcw$time_symptom_onset_absolute
  )
  
  died <- !is.na(hcw$outcome) & hcw$outcome
  
  # OBV receipt and prevention default to FALSE if not supplied (baseline case)
  received  <- if (!is.null(obv_received)) obv_received else rep(FALSE, nrow(hcw))
  prevented_death <- if (!is.null(prevented)) prevented else rep(FALSE, nrow(hcw))
  
  # Effective alive/dead status after OBV intervention
  died_after_obv <- died & !prevented_death
  
  absence_end <- rep(NA_real_, nrow(hcw))
  
  # Died: absent until simulation end
  absence_end[died_after_obv] <- duration
  
  # Survived without OBV: absent until min(recovery + 6 months, simulation end)
  surv_no_obv <- !died_after_obv & !received
  absence_end[surv_no_obv] <- pmin(
    hcw$time_outcome_absolute[surv_no_obv] + 6 * 30,
    duration
  )
  
  # Survived with OBV (either would have died but prevented, or survived naturally):
  # earlier return to work
  surv_obv <- !died_after_obv & received
  absence_end[surv_obv] <- hcw$time_outcome_absolute[surv_obv] + obv_return_days
  
  days_lost <- pmax(absence_end - absence_start, 0)
  sum(days_lost, na.rm = TRUE)
}

# =============================================================================
# Build per-run summary data frame for a given OBV scenario
#
# efficacy_name: one of OBV_EFFICACY_LEVELS, or "baseline"
# coverage_name: one of COVERAGE_LEVELS, or NULL for baseline
# seed_offset:   integer added to particle/rep index for reproducible sampling
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
      # No OBV -- use raw counts from simulation
      days_lost     <- compute_hcw_days_lost(x$tdf, x$duration)
      prevented_hcw <- 0L
    } else {
      # Deterministic seed per run so post-hoc sampling is reproducible
      run_seed <- seed_offset + x$particle_id * 1000L + x$rep
      obv      <- apply_obv_posthoc(x$tdf, efficacy, cov_spec, seed = run_seed)
      
      # Use the prevented_flag returned directly -- do NOT re-draw runif here,
      # as that would produce a different set of prevented individuals than
      # the count in obv$prevented_hcw, causing days-lost to be inconsistent.
      days_lost     <- compute_hcw_days_lost(x$tdf, x$duration,
                                             obv_received = obv$obv_received,
                                             prevented    = obv$prevented_flag)
      prevented_hcw <- obv$prevented_hcw
    }
    
    data.frame(
      scenario          = x$scenario,
      particle_id       = x$particle_id,
      rep               = x$rep,
      arm               = efficacy_name,
      coverage_scenario = if (is_baseline) "baseline" else coverage_name,
      n_infections      = x$n_infections,
      n_hcw_deaths      = x$n_hcw_deaths - prevented_hcw,
      counterfactual_hcw = x$n_hcw_deaths,
      prevented_hcw     = prevented_hcw,
      hcw_days_lost     = days_lost,
      stringsAsFactors  = FALSE
    )
  }))
}

# =============================================================================
# Aggregate: mean over reps per particle, then attach counterfactual for % averted
#
# Follows the burden-weighted approach from the reference implementation:
#   pct_averted = 100 * sum(prevented) / sum(counterfactual)  per particle
# This guarantees values in [0, 100] with no negative effectiveness.
# =============================================================================
make_particle_df <- function(run_df) {
  # Sum prevented and counterfactual over reps (burden-weighted denominator)
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
  
  # Baseline days-lost for % days-lost-averted calculation
  base_days <- pdf %>%
    filter(arm == "baseline") %>%
    select(scenario, particle_id,
           baseline_days_lost = hcw_days_lost)
  
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
# Build particle_df for all efficacy x coverage combinations (Figure 3 / 4)
# =============================================================================
make_particle_df_grid <- function(results,
                                  efficacy_levels = OBV_EFFICACY_LEVELS,
                                  coverage_levels = COVERAGE_LEVELS,
                                  seed_offset     = 0L) {
  rows <- vector("list",
                 length(efficacy_levels) * length(coverage_levels) + 1L)
  idx <- 1L
  
  # Baseline row
  rows[[idx]] <- build_run_df_obv(results, "baseline")
  idx <- idx + 1L
  
  for (eff in efficacy_levels) {
    for (cov in coverage_levels) {
      rows[[idx]] <- build_run_df_obv(results, eff, cov,
                                      seed_offset = seed_offset)
      idx <- idx + 1L
    }
  }
  
  make_particle_df(do.call(rbind, rows))
}

# =============================================================================
# Weekly cumulative time series for a metric (HCW deaths or infections)
# =============================================================================
build_weekly_ts <- function(results, metric = c("hcw_deaths", "hcw_infections"),
                            bin_width = 28,
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
      run_seed  <- seed_offset + x$particle_id * 1000L + x$rep
      obv       <- apply_obv_posthoc(cases, efficacy, cov_spec, seed = run_seed)
      hcw_rows  <- which(is_hcw)
      # Mark prevented deaths so they are excluded from the time series
      prevented_idx <- hcw_rows[obv$prevented_hcw > 0]   # approximate: use count
      # Full prevention vector for time-series exclusion
      hcw_died_idx  <- which(is_hcw & died)
      set.seed(run_seed + 1L)
      n_prevent <- obv$prevented_hcw
      if (n_prevent > 0 && length(hcw_died_idx) >= n_prevent) {
        prevented_full <- logical(nrow(cases))
        prevented_full[sample(hcw_died_idx, n_prevent)] <- TRUE
        died <- died & !prevented_full
      }
    }
    
    times <- if (metric == "hcw_deaths") {
      cases$time_outcome_absolute[died & is_hcw]
    } else {
      cases$time_infection_absolute[is_hcw]
    }
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
  
  # Mean over reps, then quantiles over particles, then cumsum
  rows %>%
    group_by(scenario, arm, particle_id, week) %>%
    summarise(value = mean(value), .groups = "drop") %>%
    arrange(scenario, arm, particle_id, week) %>%
    group_by(scenario, arm, particle_id) %>%
    mutate(value = cumsum(value)) %>%
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
# Shared ggplot theme
# =============================================================================
theme_fig <- function(base_size = 12) {
  # theme_bw(base_size = base_size) +
  theme_classic(base_size = base_size) +
    theme(
      plot.title       = element_text(face = "bold", size = base_size + 1),
      plot.subtitle    = element_text(color = "grey40", size = base_size - 2),
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position  = "top"
    )
}