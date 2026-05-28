# helper_functions.R  (analyses/04_obeldesivir_impact)
# =============================================================================
# Analysis-specific helpers for the obeldesivir (OBV PEP) impact analysis.
#
# This analysis takes the posterior from the DRC ABC-SMC calibration
# (analyses/02_model_fits), converts the 3 fitted parameters into the 4 fiber
# model parameters, downsamples to a manageable posterior sample, and runs the
# fiber branching-process model with and without obeldesivir post-exposure
# prophylaxis (PEP).
#
# Contents
#   --- Path / IO ---
#     `%||%`                       : null-coalescing operator.
#     find_repo_root()             : walk up from a directory to the repo root.
#     find_latest_abc_run_dir()    : newest timestamped ABC run dir for a scenario.
#     read_abc_posterior_step()    : read an output_step<k> particle cloud.
#
#   --- Posterior handling ---
#     downsample_posterior()       : weighted resample of the posterior particles.
#     derive_model_parameters()    : 3 fitted params -> 4 fiber model parameters,
#                                    using the SAME mapping as the calibration's
#                                    build_abc_model_args().
#
#   --- Simulation (runs on future workers) ---
#     simulate_one()               : run ONE fiber replicate (one parameter set,
#                                    one stochastic rep, one arm) and return
#                                    compact per-replicate outputs, including the
#                                    day-indices of every death and HCW death so
#                                    epidemic curves can be reconstructed.
#
#   --- Aggregation / summary ---
#     bin_counts()                 : bin a vector of event day-indices into
#                                    fixed-width time bins.
#     q_summary()                  : median + chosen quantiles of a vector.
#
# Requires fiber to be loaded (library(fiber)) and, for derive_model_parameters()
# / building model args, the calibration helpers in
# analyses/02_model_fits/helper_functions to be sourced (so build_abc_model_args()
# and solve_offspring_means_for_R0() are available).
# =============================================================================


# -----------------------------------------------------------------------------
# Path / IO
# -----------------------------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

# Walk up the directory tree from `start` until a directory containing `marker`
# (the .Rproj file) is found. Returns NULL if not found.
find_repo_root <- function(start = getwd(), marker = "obv_hcw_paper.Rproj") {
  d <- normalizePath(start, winslash = "/", mustWork = FALSE)
  repeat {
    if (file.exists(file.path(d, marker))) return(d)
    parent <- dirname(d)
    if (identical(parent, d)) return(NULL)
    d <- parent
  }
}

# Find the most recent timestamped ABC run directory for a scenario. ABC runs
# are written to <abc_outputs>/<scenario_id>_YYYYMMDD_HHMMSS[...]; the directory
# names sort chronologically, so the last one is the newest.
find_latest_abc_run_dir <- function(abc_outputs_dir, scenario_id) {
  if (!dir.exists(abc_outputs_dir)) {
    stop("ABC outputs directory not found: ", abc_outputs_dir, call. = FALSE)
  }
  dirs <- list.dirs(abc_outputs_dir, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[grepl(scenario_id, basename(dirs), fixed = TRUE)]
  if (length(dirs) == 0L) {
    stop("No ABC run directories for scenario '", scenario_id, "' under ",
         abc_outputs_dir, call. = FALSE)
  }
  dirs[order(basename(dirs))][length(dirs)]
}

# Read one output_step<k> file written by EasyABC::ABC_sequential. If `step` is
# NULL the latest completed step is used. Returns a data.frame with the standard
# column names plus a "step" attribute.
read_abc_posterior_step <- function(run_dir,
                                     step = NULL,
                                     param_names = c("R0", "prop_funeral", "hcw_risk_scalar"),
                                     stat_names  = c("takeoff", "n_deaths", "n_hcw_deaths", "duration")) {
  step_of <- function(f) as.integer(sub(".*_step([0-9]+)$", "\\1", f))
  files   <- list.files(run_dir, pattern = "^output_step[0-9]+$", full.names = TRUE)
  if (length(files) == 0L) stop("No output_step files found in ", run_dir, call. = FALSE)
  files <- files[order(step_of(files))]

  target <- if (is.null(step)) {
    files[length(files)]
  } else {
    f <- file.path(run_dir, paste0("output_step", step))
    if (!file.exists(f)) stop("output_step", step, " not found in ", run_dir, call. = FALSE)
    f
  }

  df <- utils::read.table(target, header = FALSE)
  colnames(df) <- c("weight", param_names, stat_names)
  attr(df, "step")      <- step_of(target)
  attr(df, "step_file") <- target
  df
}


# -----------------------------------------------------------------------------
# Posterior handling
# -----------------------------------------------------------------------------

# Weighted resample of an ABC particle cloud (posterior). Particles are sampled
# WITH replacement with probability proportional to their ABC weights, which is
# the standard way to turn a weighted ABC-SMC population into an unweighted
# posterior sample. Returns the 3 fitted parameters for the drawn particles.
downsample_posterior <- function(posterior,
                                 n_sets,
                                 seed = 1,
                                 param_names = c("R0", "prop_funeral", "hcw_risk_scalar")) {
  if (!"weight" %in% names(posterior)) stop("`posterior` must contain a 'weight' column.", call. = FALSE)
  w <- posterior$weight / sum(posterior$weight)
  set.seed(seed)
  idx <- sample(seq_len(nrow(posterior)), size = n_sets, replace = TRUE, prob = w)
  out <- posterior[idx, param_names, drop = FALSE]
  out <- cbind(set_id = seq_len(n_sets), particle = idx, out)
  rownames(out) <- NULL
  out
}

# Convert the 3 fitted ABC parameters (R0, prop_funeral, hcw_risk_scalar) into
# the 4 fiber model parameters, using EXACTLY the mapping the calibration used
# (build_abc_model_args() in analyses/02_model_fits/helper_functions):
#
#   mn_offspring_genPop           = (1 - prop_funeral) * R0 / D
#   mn_offspring_funeral          =      prop_funeral  * R0 / F
#   prob_hcw_cond_genPop_hospital = min(hcw_base_prob * hcw_risk_scalar, 1)
#   prob_hcw_cond_hcw_hospital    = min(hcw_base_prob * hcw_risk_scalar, 1)
#
# D and F are the scenario-level direct / funeral R0 multipliers returned by
# solve_offspring_means_for_R0(); pass the values you computed once for the
# scenario. This function is just for a transparent, saved record of the derived
# parameters; the actual model args fed to the simulator are built with
# build_abc_model_args() so the two can never drift apart.
derive_model_parameters <- function(theta, D, F_fun, hcw_base_prob = 0.25) {
  pf  <- theta$prop_funeral
  R0  <- theta$R0
  hrs <- theta$hcw_risk_scalar
  data.frame(
    set_id                        = theta$set_id,
    particle                      = theta$particle,
    R0                            = R0,
    prop_funeral                  = pf,
    hcw_risk_scalar               = hrs,
    mn_offspring_genPop           = (1 - pf) * R0 / D,
    mn_offspring_funeral          =      pf  * R0 / F_fun,
    prob_hcw_cond_genPop_hospital = pmin(hcw_base_prob * hrs, 1),
    prob_hcw_cond_hcw_hospital    = pmin(hcw_base_prob * hrs, 1),
    stringsAsFactors              = FALSE
  )
}


# -----------------------------------------------------------------------------
# Simulation (executed on future workers)
# -----------------------------------------------------------------------------
# Run a single fiber replicate. `job` is a one-row list describing what to run;
# `args_list` is the list of pre-built fiber argument lists (one per parameter
# set, produced in the main session with build_abc_model_args()); `obv_cfg`
# carries the obeldesivir settings for the OBV arm.
#
# Only fiber + base R are used here, so the function is safe to ship to
# future/PSOCK workers (load fiber on the worker via future.packages = "fiber").
#
# Reproducibility / pairing: the same `seed` is used for a (set, rep) pair in
# both arms. Each arm calls set.seed(seed) inside branching_process_main(), so
# the two trajectories share their early stochastic history and only diverge
# once the first OBV draw fires (see the OBV reproducibility caveat in
# branching_process_main()). This is a deliberate variance-reduction choice for
# the paired with/without comparison; we still summarise over many reps.
#
# Returns a compact list: scalar summaries plus the integer day-index of every
# death (and every HCW death), from which epidemic curves are reconstructed.
simulate_one <- function(job, args_list, obv_cfg) {
  a <- args_list[[job$set_id]]
  a$seed <- job$seed

  if (identical(job$arm, "obv")) {
    a$obv_pep_enabled          <- TRUE
    a$obv_pep_coverage         <- obv_cfg$coverage
    a$obv_pep_adherence        <- obv_cfg$adherence
    a$obv_pep_efficacy         <- obv_cfg$efficacy
    a$obv_pep_dpc              <- obv_cfg$dpc            %||% 1
    a$obv_pep_target_class     <- obv_cfg$target_class
    a$obv_pep_target_locations <- obv_cfg$target_locations
  }

  out <- do.call(fiber::branching_process_main, a)

  tdf <- out$tdf
  tdf <- tdf[!is.na(tdf$time_infection_absolute), , drop = FALSE]

  deaths <- !is.na(tdf$outcome) & tdf$outcome
  hcw    <- !is.na(tdf$class)   & tdf$class == "HCW"

  death_time     <- tdf$time_outcome_absolute[deaths]
  hcw_death_time <- tdf$time_outcome_absolute[deaths & hcw]

  list(
    set_id         = job$set_id,
    rep_id         = job$rep_id,
    arm            = job$arm,
    seed           = job$seed,
    n_cases        = nrow(tdf),
    n_deaths       = sum(deaths),
    n_hcw_deaths   = sum(deaths & hcw),
    # day-index (0-based) of each death; used to rebuild incidence time series.
    death_days     = as.integer(floor(death_time[!is.na(death_time)])),
    hcw_death_days = as.integer(floor(hcw_death_time[!is.na(hcw_death_time)]))
  )
}


# -----------------------------------------------------------------------------
# Aggregation / summary
# -----------------------------------------------------------------------------

# Bin a vector of (0-based) event day-indices into `n_bins` bins of width
# `bin_width` days. Returns an integer vector of counts of length n_bins. Any
# event beyond the last bin is clipped into it (should not happen when n_bins is
# derived from the global maximum day).
bin_counts <- function(day_index, bin_width, n_bins) {
  if (length(day_index) == 0L) return(integer(n_bins))
  b <- (day_index %/% bin_width) + 1L
  b[b > n_bins] <- n_bins
  b[b < 1L]     <- 1L
  tabulate(b, nbins = n_bins)
}

# Median plus an arbitrary set of quantiles of a numeric vector, returned as a
# named numeric vector. NA values are dropped.
q_summary <- function(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(setNames(rep(NA_real_, length(probs) + 1L), c("mean", paste0("q", probs * 100))))
  c(mean = mean(x), setNames(stats::quantile(x, probs = probs, names = FALSE), paste0("q", probs * 100)))
}
