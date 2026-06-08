# =============================================================================
#
# Shared helper functions for fig 1 to 4
# Source this at the top of each figgen script.
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(here)

# =============================================================================
# Constants
# =============================================================================
SCENARIO_LABELS <- c(
  WestAfrica = "West Africa (Worst)",
  DRC        = "DRC (Middle, PlusPlus)"
)
SCENARIO_COLORS <- c(
  WestAfrica = "#d73097",
  DRC        = "#2166ac"
)
OBV_EFFICACY_LEVELS <- c("obv_50","obv_60","obv_70","obv_80","obv_90")
OBV_EFFICACY_LABELS <- c("50%","60%","70%","80%","90%")
COVERAGE_LEVELS     <- c("full","ramp_high","ramp_low")
COVERAGE_LABELS     <- c("Full (100%)","Ramp high (20%->100%)","Ramp low (20%->50%)")
COVERAGE_COLORS     <- c(full = "#1a9641", ramp_high = "#fdae61", ramp_low = "#d7191c")
# Arm colors for time series plots (baseline vs OBV arms)
# Arm colors — baseline uses scenario color, OBV arms use fixed gradient
# Use get_arm_colors(sc) in figgen scripts instead of ARM_COLORS directly
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
# Load all RDS files from a coverage scenario directory
# Returns a flat list of result objects
# =============================================================================
load_results <- function(coverage_scenario,
                         base_dir = here("outputs", "simulation_fig1to3")) {
  dir_path <- file.path(base_dir, coverage_scenario)
  files    <- list.files(dir_path, pattern = "\\.rds$", full.names = TRUE)
  if (length(files) == 0) stop("No RDS files found in: ", dir_path)
  message(sprintf("Loading %d files from %s...", length(files), coverage_scenario))
  lapply(files, readRDS)
}

# =============================================================================
# Flatten results list to summary data frame (scalar metrics)
# =============================================================================
flatten_to_df <- function(results) {
  do.call(rbind, lapply(results, function(x) {
    data.frame(
      scenario          = x$scenario,
      particle_id       = x$particle_id,
      arm               = x$arm,
      rep               = x$rep,
      coverage_scenario = x$coverage_scenario,
      n_infections      = x$n_infections,
      n_hcw_infections  = x$n_hcw_infections,
      n_deaths          = x$n_deaths,
      n_hcw_deaths      = x$n_hcw_deaths,
      duration          = x$duration,
      stringsAsFactors  = FALSE
    )
  }))
}

# =============================================================================
# Compute HCW days lost from a single tdf (cases data frame)
#
# Rules:
#   No OBV (obv_pep_received == FALSE or NA):
#     Died    : time_hospitalisation_absolute (or time_symptom_onset_absolute
#               if not hospitalised) -> simulation end (duration)
#     Survived: same start -> min(time_outcome_absolute + 6*30, duration)
#
#   OBV received (obv_pep_received == TRUE):
#     Died    : same as no-OBV died
#     Survived: start -> time_outcome_absolute + OBV_RETURN_TO_WORK_DAYS
#               (earlier return to work assumption)
#
# Returns total HCW days lost for this run.
# =============================================================================
compute_hcw_days_lost <- function(cases, duration,
                                  obv_return_days = OBV_RETURN_TO_WORK_DAYS) {
  hcw <- cases[!is.na(cases$class) & cases$class == "HCW", ]
  if (nrow(hcw) == 0) return(0)
  
  # Start of absence: hospitalisation time if hospitalised, else symptom onset
  absence_start <- ifelse(
    !is.na(hcw$time_hospitalisation_absolute),
    hcw$time_hospitalisation_absolute,
    hcw$time_symptom_onset_absolute
  )
  
  died     <- !is.na(hcw$outcome) & hcw$outcome
  received <- !is.na(hcw$obv_pep_received) & hcw$obv_pep_received
  
  # Absence end
  absence_end <- rep(NA_real_, nrow(hcw))
  
  # Died (OBV or not): absence until simulation end
  absence_end[died]  <- duration
  
  # Survived, no OBV: absence until min(recovery + 6 months, simulation end)
  surv_no_obv <- !died & !received
  absence_end[surv_no_obv] <- pmin(
    hcw$time_outcome_absolute[surv_no_obv] + 6 * 30,
    duration
  )
  
  # Survived, OBV: absence until recovery + return-to-work days
  surv_obv <- !died & received
  absence_end[surv_obv] <- hcw$time_outcome_absolute[surv_obv] + obv_return_days
  
  days_lost <- pmax(absence_end - absence_start, 0)
  sum(days_lost, na.rm = TRUE)
}

# =============================================================================
# Build per-run data frame including HCW days lost
# =============================================================================
build_run_df <- function(results, obv_return_days = OBV_RETURN_TO_WORK_DAYS) {
  do.call(rbind, lapply(results, function(x) {
    days_lost <- compute_hcw_days_lost(x$tdf, x$duration,
                                       obv_return_days = obv_return_days)
    data.frame(
      scenario          = x$scenario,
      particle_id       = x$particle_id,
      arm               = x$arm,
      rep               = x$rep,
      coverage_scenario = x$coverage_scenario,
      n_infections      = x$n_infections,
      n_hcw_deaths      = x$n_hcw_deaths,
      hcw_days_lost     = days_lost,
      stringsAsFactors  = FALSE
    )
  }))
}

# =============================================================================
# Aggregate: mean over reps per particle, then attach baseline for % averted
# =============================================================================
make_particle_df <- function(run_df) {
  # Mean over reps
  pdf <- run_df %>%
    group_by(scenario, particle_id, arm, coverage_scenario) %>%
    summarise(
      n_infections  = mean(n_infections),
      n_hcw_deaths  = mean(n_hcw_deaths),
      hcw_days_lost = mean(hcw_days_lost),
      .groups = "drop"
    )
  
  # Join baseline values
  base <- pdf %>%
    filter(arm == "baseline") %>%
    select(scenario, particle_id,
           baseline_hcw_deaths  = n_hcw_deaths,
           baseline_days_lost   = hcw_days_lost)
  
  pdf %>%
    left_join(base, by = c("scenario", "particle_id")) %>%
    mutate(
      pct_hcw_deaths_averted = ifelse(
        arm != "baseline" & baseline_hcw_deaths > 0,
        100 * (baseline_hcw_deaths - n_hcw_deaths) / baseline_hcw_deaths,
        NA_real_
      ),
      pct_days_lost_averted = ifelse(
        arm != "baseline" & baseline_days_lost > 0,
        100 * (baseline_days_lost - hcw_days_lost) / baseline_days_lost,
        NA_real_
      )
    )
}

# =============================================================================
# Weekly time series for a metric from tdf list
# =============================================================================
build_weekly_ts <- function(results, metric = c("hcw_deaths", "hcw_infections"),
                            bin_width = 28, arms = NULL) {
  metric <- match.arg(metric)
  if (!is.null(arms)) results <- Filter(function(x) x$arm %in% arms, results)
  
  max_day <- max(sapply(results, function(x) x$duration), na.rm = TRUE)
  breaks  <- seq(0, ceiling(max_day / bin_width) * bin_width, by = bin_width)
  mids    <- breaks[-length(breaks)] + bin_width / 2
  
  rows <- do.call(rbind, lapply(results, function(x) {
    cases  <- x$tdf
    is_hcw <- cases$class == "HCW"
    died   <- !is.na(cases$outcome) & cases$outcome
    
    times <- if (metric == "hcw_deaths") {
      cases$time_outcome_absolute[died & is_hcw]
    } else {
      cases$time_infection_absolute[is_hcw]
    }
    counts <- hist(times[!is.na(times)], breaks = breaks, plot = FALSE)$counts
    
    data.frame(
      scenario          = x$scenario,
      arm               = x$arm,
      coverage_scenario = x$coverage_scenario,
      particle_id       = x$particle_id,
      rep               = x$rep,
      week              = mids,
      value             = counts,
      stringsAsFactors  = FALSE
    )
  }))
  
  # Mean over reps, then quantiles over particles, then cumsum
  rows %>%
    group_by(scenario, arm, coverage_scenario, particle_id, week) %>%
    summarise(value = mean(value), .groups = "drop") %>%
    arrange(scenario, arm, coverage_scenario, particle_id, week) %>%
    group_by(scenario, arm, coverage_scenario, particle_id) %>%
    mutate(value = cumsum(value)) %>%
    ungroup() %>%
    group_by(scenario, arm, coverage_scenario, week) %>%
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
  theme_bw(base_size = base_size) +
    theme(
      plot.title       = element_text(face = "bold", size = base_size + 1),
      plot.subtitle    = element_text(color = "grey40", size = base_size - 2),
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position  = "top"
    )
}
