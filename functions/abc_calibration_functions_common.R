# abc_calibration_functions_common.R
# -----------------------------------------------------------------------------
# Generic ABC-SMC helpers shared by BOTH fitting approaches (HCW-risk and
# NPI-efficacy). Nothing here depends on which parameters are fitted, so this
# file is sourced by abc_calibration_functions_hcwRisk.R and
# abc_calibration_functions_npi.R (and by their run scripts and PSOCK workers).
#
# Contents:
#   AVAILABLE_ABC_METRICS    : per-replicate metrics abc_summarise() can compute
#   DEFAULT_SUMMARY_STATS    : default FITTED summary set (takeoff + 3 means)
#   weekly_incidence_counts() : event day-indices -> fixed-width (weekly) counts
#   peak_stats_from_death_days(): weekly death curve -> (time_to_peak, peak_height)
#   abc_summarise()          : per-replicate summary of a configurable metric set
#   run_one_abc_replicate()  : do.call(branching_process_main) + abc_summarise()
#   aggregate_abc_reps()     : per-replicate metric matrix -> fitted summary vector
#   run_abc_particle()       : n_reps replicates -> aggregated fitted summary vector
#   save_abc_config()        : serialise worker config; advertise via FIBER_ABC_CONFIG
#   make_abc_output_dir()    : fresh timestamped per-run output directory
#   with_abc_output_dir()    : run an expression with cwd set to that directory
#   abc_progress()           : print progress from output_step*/tolerance_step* files
#   abc_compare_steps()      : weighted summary across completed steps
#   reconstruct_abc_result() : rebuild an ABC_sequential()-style result from disk
#
# WHICH SUMMARIES ARE FITTED is configurable, not hard-coded: a scheme passes a
# `summary_stats` vector (any subset of "takeoff" + AVAILABLE_ABC_METRICS) and
# run_abc_particle() returns exactly those, in that order, so the same functions
# fit the historical four-number summary OR an extended set (e.g. adding the
# weekly-death-curve features time_to_peak / peak_height) with no forked code.
#
# The disk-inspection trio take `param_names` (default = the HCW-risk triple);
# callers with a different parameterisation pass their own, e.g.
# c("R0", "prop_funeral", "npi_scaler"). They also take `stat_names` (default =
# the four canonical summaries takeoff/n_deaths/n_hcw_deaths/duration); a scheme
# that fits extra summaries (e.g. the NPI-peak scheme's time_to_peak/peak_height)
# passes the full vector so the on-disk columns are labelled correctly. Requires
# fiber to be loaded.
# -----------------------------------------------------------------------------

# Per-replicate metrics abc_summarise() can compute from one
# branching_process_main() run, and the DEFAULT set of FITTED summary statistics
# a model returns / ABC matches. "takeoff" is the fraction of replicates that
# exceeded the death threshold (a property of the replicate ensemble, not of a
# single run), so it is a fittable summary but NOT a per-replicate metric.
AVAILABLE_ABC_METRICS    <- c("n_cases", "n_deaths", "n_hcw_deaths", "duration",
                              "time_to_peak", "peak_height")
DEFAULT_SUMMARY_STATS    <- c("takeoff", "n_deaths", "n_hcw_deaths", "duration")
DEFAULT_PEAK_BIN_WIDTH   <- 7L            # days per incidence bin (weekly)
DEFAULT_PEAK_TIME_ORIGIN <- "first_death" # "first_death" or "outbreak_start"


# Bin a vector of (0-based) event day-indices into fixed-width bins; returns the
# integer count in each bin from the first up to the bin holding the latest
# event. Pure base R. Non-finite days are dropped; negatives clamped to 0.
weekly_incidence_counts <- function(day_index, bin_width = DEFAULT_PEAK_BIN_WIDTH) {
  day_index <- day_index[is.finite(day_index)]
  if (length(day_index) == 0L) return(integer(0))
  day_index <- as.integer(day_index)
  day_index[day_index < 0L] <- 0L
  bins <- (day_index %/% bin_width) + 1L          # 1-based bin index
  tabulate(bins, nbins = max(bins))
}


# Peak of the (weekly) death-incidence curve for ONE replicate. `death_days` are
# the floored, 0-based, absolute day-indices of every death; `first_death_day`
# anchors the "first_death" origin. Returns c(time_to_peak, peak_height) (both 0
# with no deaths). The peak week is the FIRST week attaining the maximum.
#   time_origin = "first_death"    -> whole weeks from the first-death week to the
#                                     peak week (>= 0; anchored like `duration`).
#   time_origin = "outbreak_start" -> midpoint day of the peak week, from t = 0.
peak_stats_from_death_days <- function(death_days,
                                       first_death_day,
                                       bin_width   = DEFAULT_PEAK_BIN_WIDTH,
                                       time_origin = DEFAULT_PEAK_TIME_ORIGIN) {
  weekly <- weekly_incidence_counts(death_days, bin_width)
  if (length(weekly) == 0L || all(weekly == 0L)) {
    return(c(time_to_peak = 0, peak_height = 0))
  }
  peak_bin    <- which.max(weekly)
  peak_height <- as.numeric(weekly[peak_bin])
  time_to_peak <- if (identical(time_origin, "outbreak_start")) {
    (peak_bin - 0.5) * bin_width
  } else {
    fd_bin <- (as.integer(floor(first_death_day)) %/% bin_width) + 1L
    max((peak_bin - fd_bin) * bin_width, 0)
  }
  c(time_to_peak = as.numeric(time_to_peak), peak_height = peak_height)
}


# Per-replicate summary. Computes the requested `metrics` (any subset of
# AVAILABLE_ABC_METRICS) from one branching_process_main() run and returns them
# as a named numeric vector in `metrics` order. The default reproduces the
# historical four-number summary exactly. The weekly-death-curve features
# (time_to_peak / peak_height) are only computed when requested.
abc_summarise <- function(out,
                          metrics     = c("n_cases", "n_deaths", "n_hcw_deaths", "duration"),
                          bin_width   = DEFAULT_PEAK_BIN_WIDTH,
                          time_origin = DEFAULT_PEAK_TIME_ORIGIN) {
  vals <- c(n_cases = 0, n_deaths = 0, n_hcw_deaths = 0, duration = 0,
            time_to_peak = 0, peak_height = 0)

  tdf <- out$tdf
  tdf <- tdf[!is.na(tdf$time_infection_absolute), , drop = FALSE]
  if (nrow(tdf) == 0L) return(vals[metrics])

  deaths <- !is.na(tdf$outcome) & tdf$outcome
  hcw    <- !is.na(tdf$class) & tdf$class == "HCW"

  vals["n_cases"]      <- nrow(tdf)
  vals["n_deaths"]     <- sum(deaths)
  vals["n_hcw_deaths"] <- sum(deaths & hcw)

  if (vals[["n_deaths"]] > 0L) {
    death_t       <- tdf$time_outcome_absolute[deaths]
    death_t       <- death_t[!is.na(death_t)]
    # Duration is measured from the first death to the last event, so it tracks
    # the part of the outbreak that would be observable in surveillance data.
    t_first_death <- min(death_t)
    t_last_event  <- max(tdf$time_outcome_absolute, na.rm = TRUE)
    vals["duration"] <- max(t_last_event - t_first_death, 0)

    if (any(c("time_to_peak", "peak_height") %in% metrics)) {
      peak <- peak_stats_from_death_days(
        death_days      = as.integer(floor(death_t)),
        first_death_day = t_first_death,
        bin_width       = bin_width,
        time_origin     = time_origin
      )
      vals["time_to_peak"] <- peak[["time_to_peak"]]
      vals["peak_height"]  <- peak[["peak_height"]]
    }
  }

  vals[metrics]
}


run_one_abc_replicate <- function(args, ...) {
  out <- do.call(branching_process_main, args)
  abc_summarise(out, ...)
}


# Aggregate a [metric x replicate] matrix into the FITTED summary vector.
# `took_off` flags the replicates that exceeded the death threshold;
# `summary_stats` names the statistics to return (and their order). "takeoff"
# maps to the take-off fraction; every other name maps to the mean over taken-off
# replicates (0 when none took off). `reps` must have a row for each non-"takeoff"
# entry of `summary_stats`.
aggregate_abc_reps <- function(reps, took_off, summary_stats = DEFAULT_SUMMARY_STATS) {
  any_off <- any(took_off)
  out <- vapply(summary_stats, function(s) {
    if (identical(s, "takeoff")) return(mean(took_off))
    if (!any_off) return(0)
    mean(reps[s, took_off])
  }, numeric(1))
  names(out) <- summary_stats
  out
}


# Run ONE ABC particle: `n_reps` stochastic replicates of `args`, summarise each
# on the metrics needed for `summary_stats`, then aggregate to the fitted summary
# vector. Shared by the scheme-specific parallel/serial model functions, which
# differ only in how they build `args`. NULL config values fall back to defaults.
run_abc_particle <- function(args,
                             n_reps,
                             takeoff_death_threshold,
                             summary_stats = NULL,
                             bin_width     = NULL,
                             time_origin   = NULL) {
  if (is.null(summary_stats)) summary_stats <- DEFAULT_SUMMARY_STATS
  if (is.null(bin_width))     bin_width     <- DEFAULT_PEAK_BIN_WIDTH
  if (is.null(time_origin))   time_origin   <- DEFAULT_PEAK_TIME_ORIGIN

  unknown <- setdiff(summary_stats, c("takeoff", AVAILABLE_ABC_METRICS))
  if (length(unknown) > 0L) {
    stop("Unknown summary_stats: ", paste(unknown, collapse = ", "),
         ". Available: takeoff, ", paste(AVAILABLE_ABC_METRICS, collapse = ", "),
         ".", call. = FALSE)
  }

  # Per-replicate metrics needed: the fitted ones (minus the take-off fraction)
  # plus n_deaths, which is always required to score the take-off threshold.
  rep_metrics <- union("n_deaths", setdiff(summary_stats, "takeoff"))

  reps <- vapply(seq_len(n_reps), function(i) {
    out <- do.call(branching_process_main, args)
    abc_summarise(out, metrics = rep_metrics,
                  bin_width = bin_width, time_origin = time_origin)
  }, numeric(length(rep_metrics)))
  # vapply returns a bare vector when rep_metrics has length 1; make it a matrix
  # so reps["n_deaths", ] and aggregate_abc_reps() behave uniformly.
  if (is.null(dim(reps))) {
    reps <- matrix(reps, nrow = 1L, dimnames = list(rep_metrics, NULL))
  } else {
    rownames(reps) <- rep_metrics
  }

  took_off <- reps["n_deaths", ] >= takeoff_death_threshold
  aggregate_abc_reps(reps, took_off, summary_stats)
}


save_abc_config <- function(config, file = tempfile(fileext = ".rds")) {
  required <- c("setup_path", "functions_path", "r0_path",
                "scenario_csv", "scenario_id")
  missing_keys <- setdiff(required, names(config))
  if (length(missing_keys) > 0L) {
    stop("save_abc_config(): config is missing required key(s): ",
         paste(missing_keys, collapse = ", "), call. = FALSE)
  }
  saveRDS(config, file = file)
  Sys.setenv(FIBER_ABC_CONFIG = file)
  invisible(file)
}


make_abc_output_dir <- function(base_dir,
                                scenario_id,
                                label = NULL,
                                suffix = NULL,
                                subdir = "abc_outputs",
                                timestamp = TRUE) {
  if (missing(base_dir) || is.null(base_dir) || !nzchar(base_dir)) {
    stop("`base_dir` is required.", call. = FALSE)
  }
  if (missing(scenario_id) || is.null(scenario_id) || !nzchar(scenario_id)) {
    stop("`scenario_id` is required.", call. = FALSE)
  }

  parts <- scenario_id
  if (isTRUE(timestamp)) {
    parts <- paste(parts, format(Sys.time(), "%Y%m%d_%H%M%S"), sep = "_")
  }
  if (!is.null(label) && nzchar(label)) {
    parts <- paste(parts, label, sep = "_")
  }
  # `suffix` is appended last so any run-identifying tag (e.g. the
  # NBREPS_X_NBSIMUL_Y settings tag) lands at the very END of the directory
  # name, after the timestamp and label.
  if (!is.null(suffix) && nzchar(suffix)) {
    parts <- paste(parts, suffix, sep = "_")
  }

  out <- file.path(base_dir, subdir, parts)
  dedup_suffix <- 0L
  candidate <- out
  while (dir.exists(candidate)) {
    dedup_suffix <- dedup_suffix + 1L
    candidate <- paste0(out, "_", sprintf("%02d", dedup_suffix))
  }
  out <- candidate

  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  out
}


with_abc_output_dir <- function(output_dir, expr) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  old <- setwd(output_dir)
  on.exit(setwd(old), add = TRUE)
  force(expr)
}


abc_progress <- function(dir = getwd(),
                         tolerance_target = 1.0,
                         param_names = c("R0", "prop_funeral", "hcw_risk_scalar"),
                         stat_names  = c("takeoff", "n_deaths", "n_hcw_deaths", "duration")) {

  step_num     <- function(f) as.integer(sub(".*_step([0-9]+)$", "\\1", f))
  sort_by_step <- function(f) f[order(step_num(f))]

  out_files <- sort_by_step(list.files(dir, pattern = "^output_step[0-9]+$",       full.names = TRUE))
  tol_files <- sort_by_step(list.files(dir, pattern = "^tolerance_step[0-9]+$",    full.names = TRUE))
  sim_files <- sort_by_step(list.files(dir, pattern = "^n_simul_tot_step[0-9]+$",  full.names = TRUE))

  if (length(out_files) == 0L) {
    cat("No steps completed yet.\n")
    return(invisible(NULL))
  }

  cat(sprintf("Steps completed: %d\n", length(out_files)))

  if (length(sim_files) > 0L) {
    cum_sims <- as.numeric(readLines(sim_files[length(sim_files)]))
    cat(sprintf("Total simulations so far: %s\n", format(cum_sims, big.mark = ",")))
  }

  if (length(tol_files) > 0L) {
    tols <- vapply(tol_files, function(f) as.numeric(readLines(f)), numeric(1))
    cat("\nTolerance trajectory:\n")
    for (i in seq_along(tols)) {
      cat(sprintf("  After step %d: %.3f\n", step_num(tol_files)[i], tols[i]))
    }
    cat(sprintf("\nTarget tolerance: %.3f\n", tolerance_target))

    current <- tols[length(tols)]
    if (current <= tolerance_target) {
      cat("Reached target -- algorithm should terminate after this step.\n")
    } else if (length(tols) >= 2L) {
      ratio <- tols[length(tols)] / tols[length(tols) - 1L]
      if (ratio < 1) {
        remaining <- ceiling(log(tolerance_target / current) / log(ratio))
        cat(sprintf(
          "Projected remaining steps: ~%d (assuming current shrinkage rate continues)\n",
          remaining
        ))
      } else {
        cat("Tolerance is not decreasing -- the algorithm may be stuck.\n")
      }
    }
  }

  last_out <- read.table(out_files[length(out_files)], header = FALSE)
  colnames(last_out) <- c("weight", param_names, stat_names)
  cat(sprintf("\nLatest particle cloud (step %d):\n",
              step_num(out_files)[length(out_files)]))
  for (nm in param_names) {
    cat(sprintf("  %-14s median = %.3f, 95%% CI = [%.3f, %.3f]\n",
                paste0(nm, ":"),
                median(last_out[[nm]]),
                quantile(last_out[[nm]], 0.025),
                quantile(last_out[[nm]], 0.975)))
  }

  invisible(last_out)
}


abc_compare_steps <- function(dir = getwd(),
                              param_names = c("R0", "prop_funeral", "hcw_risk_scalar"),
                              stat_names  = c("takeoff", "n_deaths", "n_hcw_deaths", "duration")) {

  step_num     <- function(f) as.integer(sub(".*_step([0-9]+)$", "\\1", f))
  sort_by_step <- function(f) f[order(step_num(f))]

  out_files <- sort_by_step(list.files(dir, pattern = "^output_step[0-9]+$",      full.names = TRUE))
  tol_files <- sort_by_step(list.files(dir, pattern = "^tolerance_step[0-9]+$",   full.names = TRUE))
  sim_files <- sort_by_step(list.files(dir, pattern = "^n_simul_tot_step[0-9]+$", full.names = TRUE))

  if (length(out_files) == 0L) {
    cat("No steps completed yet.\n")
    return(invisible(NULL))
  }

  wq <- function(x, w, p) {
    ord <- order(x); xs <- x[ord]; ws <- w[ord] / sum(w)
    approx(cumsum(ws), xs, p, rule = 2, ties = "ordered")$y
  }

  tol_lookup <- if (length(tol_files) > 0L) {
    setNames(vapply(tol_files, function(f) as.numeric(readLines(f)), numeric(1)),
             as.character(step_num(tol_files)))
  } else c()
  cum_lookup <- if (length(sim_files) > 0L) {
    setNames(vapply(sim_files, function(f) as.numeric(readLines(f)), numeric(1)),
             as.character(step_num(sim_files)))
  } else c()

  rows <- list()
  prev_cum <- 0
  for (f in out_files) {
    s  <- step_num(f)
    df <- read.table(f, header = FALSE)
    colnames(df) <- c("weight", param_names, stat_names)
    w <- df$weight / sum(df$weight)

    cum_s     <- if (as.character(s) %in% names(cum_lookup)) cum_lookup[[as.character(s)]] else NA_real_
    sims_step <- if (!is.na(cum_s)) cum_s - prev_cum else NA_real_
    if (!is.na(cum_s)) prev_cum <- cum_s

    row <- data.frame(
      step           = s,
      tol            = if (as.character(s) %in% names(tol_lookup)) {
        round(tol_lookup[[as.character(s)]], 3)
      } else NA_real_,
      sims_this_step = sims_step,
      cum_sims       = cum_s,
      ESS            = round(1 / sum(w^2), 1),
      stringsAsFactors = FALSE
    )

    # Weighted mean of EVERY fitted summary statistic, generic over stat_names so
    # whatever was fitted is reported (the canonical four, or an extended set with
    # e.g. time_to_peak / peak_height). The canonical four keep their established
    # column names + rounding (takeoff to 3 dp, counts/duration to integer); any
    # other stat gets a mean_<stat> column rounded by magnitude.
    stat_label <- list(takeoff = list(col = "mean_takeoff",  dp = 3L),
                       n_deaths = list(col = "mean_deaths",   dp = 0L),
                       n_hcw_deaths = list(col = "mean_hcw",  dp = 0L),
                       duration = list(col = "mean_duration", dp = 0L))
    for (nm in stat_names) {
      m   <- sum(w * df[[nm]])
      fmt <- stat_label[[nm]]
      col <- if (!is.null(fmt)) fmt$col else paste0("mean_", nm)
      dp  <- if (!is.null(fmt)) fmt$dp  else if (abs(m) < 10) 3L else 0L
      row[[col]] <- round(m, dp)
    }

    for (nm in param_names) {
      row[[paste0(nm, "_med")]] <- round(wq(df[[nm]], w, 0.500), 3)
      row[[paste0(nm, "_lo")]]  <- round(wq(df[[nm]], w, 0.025), 3)
      row[[paste0(nm, "_hi")]]  <- round(wq(df[[nm]], w, 0.975), 3)
    }

    rows[[length(rows) + 1L]] <- row
  }

  do.call(rbind, rows)
}


reconstruct_abc_result <- function(dir = getwd(),
                                   step = NULL,
                                   param_names = c("R0", "prop_funeral", "hcw_risk_scalar"),
                                   stat_names  = c("takeoff", "n_deaths", "n_hcw_deaths", "duration")) {

  step_num <- function(f) as.integer(sub(".*_step([0-9]+)$", "\\1", f))

  out_files <- list.files(dir, pattern = "^output_step[0-9]+$", full.names = TRUE)
  out_files <- out_files[order(step_num(out_files))]
  if (length(out_files) == 0L) stop("No output_step files found in ", dir)

  if (is.null(step)) {
    target_file <- out_files[length(out_files)]
    step <- step_num(target_file)
  } else {
    target_file <- file.path(dir, paste0("output_step", step))
    if (!file.exists(target_file)) {
      stop("output_step", step, " not found in ", dir, call. = FALSE)
    }
  }

  message("Reconstructing result from step ", step)

  df <- read.table(target_file, header = FALSE)
  colnames(df) <- c("weight", param_names, stat_names)

  tol_file <- file.path(dir, paste0("tolerance_step", step))
  epsilon  <- if (file.exists(tol_file)) as.numeric(readLines(tol_file)) else NA_real_

  sim_file <- file.path(dir, paste0("n_simul_tot_step", step))
  nsim     <- if (file.exists(sim_file)) as.numeric(readLines(sim_file)) else NA_real_

  # Normalisation SDs are computed from step 1 by ABC_sequential. Recompute the
  # same way so the reconstructed object behaves like a live one.
  step1_file <- file.path(dir, "output_step1")
  if (file.exists(step1_file)) {
    step1 <- read.table(step1_file, header = FALSE)
    stat_cols <- seq(2L + length(param_names), 1L + length(param_names) + length(stat_names))
    stats_norm <- apply(step1[, stat_cols], 2, sd)
  } else {
    stats_norm <- apply(df[, stat_names], 2, sd)
  }

  list(
    param                   = as.matrix(df[, param_names]),
    stats                   = as.matrix(df[, stat_names]),
    weights                 = df$weight / sum(df$weight),
    stats_normalization     = as.numeric(stats_norm),
    epsilon                 = epsilon,
    nsim                    = nsim,
    computime               = NA_real_,
    step_reconstructed_from = step
  )
}


