# =============================================================================
# 02_extract_figure3new.R
# Extract weekly HCW death time series, particle-level summaries, and
# PEP uptake/DPC summaries for the four coverage/DPC scenarios x three
# efficacy arms (figure3_new).
#
# Outputs:
#   figure_3new_weekly_ts.csv             : weekly mean +/- 95% CI HCW deaths
#   figure_3new_particle_summary.csv      : particle-level deaths averted % per arm
#   figure_3new_pep_uptake_summary.csv    : % HCW receiving PEP and DPC distribution per run
#
# All arms are run under the DRC_conflict trajectory. The scenario column
# parsed from filenames (DRC_conflict) is normalised to "DRC" for downstream
# grouping consistency.
# =============================================================================
source(here::here("analyses", "03_figure_template", "helper_functions_figure_1to4.R"))
library(tidyr)

ARM_NAMES_3NEW <- c(
  "with_conflict_mid", "with_conflict_lo", "with_conflict_hi",
  "cov_conflict_mid",  "cov_conflict_lo",  "cov_conflict_hi",
  "dpc_conflict_mid",  "dpc_conflict_lo",  "dpc_conflict_hi",
  "optimistic_mid",    "optimistic_lo",    "optimistic_hi",
  "no_pep"
)

# =============================================================================
# extract_weekly_ts_safe
#
# Parallel-safe replacement for extract_weekly_ts. Skips files where
# max(time_outcome_absolute) is non-finite rather than erroring out.
# =============================================================================
extract_weekly_ts_safe <- function(arm_dir, bin_width = 7,
                                   n_workers = 10L,
                                   base_dir = here("outputs", "simulation")) {
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
    if (i %% 100 == 0) message(sprintf("  file %d / %d", i, length(files)))
    f <- files[[i]]
    
    tryCatch({
      x <- readRDS(f)
      
      fname <- tools::file_path_sans_ext(basename(f))
      parts <- regmatches(fname, regexec("^(.+)_p(\\d+)_r(\\d+)$", fname))[[1]]
      sc    <- if (length(parts) == 4) parts[2] else NA_character_
      pid   <- if (length(parts) == 4) as.integer(parts[3]) else NA_integer_
      rep   <- if (length(parts) == 4) as.integer(parts[4]) else NA_integer_
      
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
      
      cap <- max(tdf$time_outcome_absolute, na.rm = TRUE)
      if (!is.finite(cap)) {
        message(sprintf("  SKIP (non-finite cap): %s", basename(f)))
        return(NULL)
      }
      
      cap    <- ceiling(cap / bin_width) * bin_width
      breaks <- seq(0, cap, by = bin_width)
      mids   <- head(breaks, -1) + bin_width / 2
      
      .extract <- function(cases, arm_label) {
        is_hcw <- !is.na(cases$class) & cases$class == "HCW"
        died   <- !is.na(cases$outcome) & cases$outcome
        metrics <- list(
          deaths               = cases$time_outcome_absolute[died],
          infections           = cases$time_infection_absolute,
          hcw_deaths_incidence = cases$time_outcome_absolute[died & is_hcw],
          hcw_deaths           = cases$time_outcome_absolute[died & is_hcw]
        )
        do.call(rbind, lapply(names(metrics), function(m) {
          vals   <- metrics[[m]]
          vals   <- vals[!is.na(vals) & vals <= cap]
          counts <- hist(vals, breaks = breaks, plot = FALSE)$counts
          data.frame(
            scenario = sc, particle_id = pid, rep = rep,
            arm = arm_label, week = mids, metric = m,
            value = counts, stringsAsFactors = FALSE
          )
        }))
      }
      
      rbind(.extract(cases_base, "baseline"), .extract(tdf, "obv"))
      
    }, error = function(e) {
      message(sprintf("  ERROR in file %s: %s", basename(f), conditionMessage(e)))
      NULL
    })
  }, future.packages = c("here"), future.seed = TRUE)
  
  results <- results[!sapply(results, is.null)]
  if (length(results) == 0) return(NULL)
  do.call(rbind, results)
}

# =============================================================================
# extract_pep_uptake_summary
#
# For each run, among HCW exposures eligible for antiviral PEP:
#   - % who received PEP (tdf recipients + all prevented cases)
#   - DPC distribution among infected HCW in tdf who received PEP
# =============================================================================
extract_pep_uptake_summary <- function(arm_dir,
                                       arm_label = arm_dir,
                                       n_workers = 10L,
                                       base_dir  = here("outputs", "simulation")) {
  library(future)
  library(future.apply)
  
  dir_path <- file.path(base_dir, arm_dir)
  files    <- list.files(dir_path, pattern = "\\.rds$", full.names = TRUE)
  if (length(files) == 0) stop("No RDS files found in: ", dir_path)
  message(sprintf("Extracting PEP uptake summaries from %d files in arm '%s'...",
                  length(files), arm_dir))
  
  plan(multisession, workers = min(n_workers, future::availableCores()))
  on.exit(plan(sequential), add = TRUE)
  
  results <- future_lapply(seq_along(files), function(i) {
    f <- files[[i]]
    
    tryCatch({
      x <- readRDS(f)
      
      fname <- tools::file_path_sans_ext(basename(f))
      parts <- regmatches(fname, regexec("^(.+)_p(\\d+)_r(\\d+)$", fname))[[1]]
      sc    <- if (length(parts) == 4) parts[2] else NA_character_
      pid   <- if (length(parts) == 4) as.integer(parts[3]) else NA_integer_
      rep   <- if (length(parts) == 4) as.integer(parts[4]) else NA_integer_
      
      tdf       <- x$tdf
      prevented <- x$prevented_completed
      
      hcw_tdf <- tdf[!is.na(tdf$class) & tdf$class == "HCW", ]
      
      # Use data.frame() instead of prevented[FALSE, ] to avoid length-0
      # vector vs scalar mismatch in subsequent data.frame() calls
      hcw_prevented <- if (!is.null(prevented) && nrow(prevented) > 0) {
        prevented[!is.na(prevented$class) & prevented$class == "HCW", ]
      } else {
        data.frame()
      }
      
      eligible_tdf   <- hcw_tdf[!is.na(hcw_tdf$obv_pep_eligible) & hcw_tdf$obv_pep_eligible, ]
      n_eligible     <- nrow(eligible_tdf) + nrow(hcw_prevented)
      n_received_tdf <- sum(!is.na(eligible_tdf$obv_pep_received) & eligible_tdf$obv_pep_received)
      n_received     <- n_received_tdf + nrow(hcw_prevented)
      pct_received   <- ifelse(n_eligible > 0, 100 * n_received / n_eligible, NA_real_)
      
      dpc_vals <- eligible_tdf$obv_pep_dpc[
        !is.na(eligible_tdf$obv_pep_received) & eligible_tdf$obv_pep_received
      ]
      dpc_vals <- dpc_vals[!is.na(dpc_vals)]
      
      data.frame(
        scenario     = sc,
        particle_id  = pid,
        rep          = rep,
        arm          = arm_label,
        n_eligible   = n_eligible,
        n_received   = n_received,
        pct_received = pct_received,
        dpc_mean     = if (length(dpc_vals) > 0) mean(dpc_vals)   else NA_real_,
        dpc_median   = if (length(dpc_vals) > 0) median(dpc_vals) else NA_real_,
        dpc_min      = if (length(dpc_vals) > 0) min(dpc_vals)    else NA_real_,
        dpc_max      = if (length(dpc_vals) > 0) max(dpc_vals)    else NA_real_,
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      message(sprintf("ERROR in file: %s -- %s", f, conditionMessage(e)))
      NULL
    })
  }, future.packages = c("here"), future.seed = TRUE)
  
  results <- results[!sapply(results, is.null)]
  do.call(rbind, results)
}

# =============================================================================
# 1. Weekly time series (incident + cumulative HCW deaths)
# =============================================================================
message("Extracting weekly time series for figure3_new arms...")

raw_ts_list <- lapply(ARM_NAMES_3NEW, function(arm_name) {
  message(sprintf("  Processing arm: %s", arm_name))
  df <- extract_weekly_ts_safe(
    arm_dir   = file.path("conflict_dpc", arm_name),
    bin_width = 7,
    n_workers = 10L
  )
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df <- df[df$arm == "obv", ]
  df$arm      <- arm_name
  df$scenario <- "DRC"  # normalise: DRC_conflict -> DRC
  df
})

raw_ts_3new <- do.call(rbind, raw_ts_list[!sapply(raw_ts_list, is.null)])

# Determine the global maximum week across ALL arms and ALL particles.
# Padding to global max ensures n_val stays constant across the full time
# range, preventing the cumulative mean from dropping at the tail.
global_max_week <- max(raw_ts_3new$week, na.rm = TRUE)

ts_quantiles_3new <- raw_ts_3new %>%
  group_by(scenario, arm, particle_id, week, metric) %>%
  summarise(value = mean(value), .groups = "drop") %>%
  # Pad each particle to global_max_week with 0 incidence
  group_by(scenario, arm, metric) %>%
  group_modify(~ {
    all_weeks <- seq(min(.x$week), global_max_week, by = 7)
    tidyr::complete(.x,
                    particle_id = unique(.x$particle_id),
                    week        = all_weeks,
                    fill        = list(value = 0))
  }) %>%
  ungroup() %>%
  arrange(scenario, arm, particle_id, metric, week) %>%
  # Apply cumsum for cumulative metric only; incidence stays as-is
  group_by(scenario, arm, particle_id, metric) %>%
  mutate(value = if (unique(metric) == "hcw_deaths") cumsum(value) else value) %>%
  ungroup() %>%
  group_by(scenario, arm, week, metric) %>%
  summarise(
    mean_val = mean(value, na.rm = TRUE),
    sd_val   = sd(value,  na.rm = TRUE),
    n_val    = sum(!is.na(value)),
    .groups  = "drop"
  ) %>%
  mutate(
    se_val = sd_val / sqrt(n_val),
    ci_lo  = mean_val - 1.96 * se_val,
    ci_hi  = mean_val + 1.96 * se_val,
    week   = week / 7  # convert day midpoints to weeks for plot script
  )

save_figure_data(ts_quantiles_3new, "figure_3new_weekly_ts.csv")
message("Weekly time series saved.")

# =============================================================================
# 2. Particle-level run summary (deaths averted % per arm)
# =============================================================================
message("Extracting particle-level run summaries for figure3_new arms...")

run_df_3new <- do.call(rbind, lapply(ARM_NAMES_3NEW, function(arm_name) {
  df <- extract_run_summary(
    arm_dir   = file.path("conflict_dpc", arm_name),
    arm_label = arm_name,
    n_workers = 10L
  )
  df$scenario <- "DRC"  # normalise: DRC_conflict -> DRC
  df
}))

# Use conflict-dpc variant: all arms share no_pep as common denominator
particle_df_3new <- make_particle_df(run_df_3new)
save_figure_data(particle_df_3new, "figure_3new_particle_summary.csv")
message("Particle-level summary saved.")

# =============================================================================
# 3. PEP uptake and DPC distribution summary (excludes no_pep arm)
# =============================================================================
message("Extracting PEP uptake and DPC summaries for figure3_new arms...")

pep_uptake_3new <- do.call(rbind, lapply(ARM_NAMES_3NEW, function(arm_name) {
  if (arm_name == "no_pep") return(NULL)
  df <- extract_pep_uptake_summary(
    arm_dir   = file.path("conflict_dpc", arm_name),
    arm_label = arm_name,
    n_workers = 10L
  )
  df$scenario <- "DRC"  # normalise: DRC_conflict -> DRC
  df
}))

save_figure_data(pep_uptake_3new, "figure_3new_pep_uptake_summary.csv")
message("PEP uptake summary saved.")

message("Figure 3new extraction complete.")