# abc_calibration_functions_npi_peak.R
# -----------------------------------------------------------------------------
# NPI-EFFICACY ABC approach, EXTENDED with two epidemic-CURVE summary statistics
# on top of the four scalar summaries the base NPI scheme fits to. This is an
# EXPERIMENTAL scheme: it fits the SAME three parameters as
# abc_calibration_functions_npi.R -- (R0, prop_funeral, npi_scaler) -- and reuses
# that file's per-particle machinery (build_abc_model_args(),
# npi_efficacy_from_scaler(), the metadata helpers), but it asks ABC to match TWO
# extra features of the weekly death curve:
#
#   time_to_peak : time (in days) to the peak week of DEATH incidence.
#   peak_height  : deaths in that peak week (i.e. peak weekly death incidence).
#
# So the fitted summary vector grows from
#   (takeoff, n_deaths, n_hcw_deaths, duration)                         [4 stats]
# to
#   (takeoff, n_deaths, n_hcw_deaths, duration, time_to_peak, peak_height) [6].
#
# WHY a separate file (not edits to abc_calibration_functions_npi.R): the base
# scheme is used by the production DRC / West-Africa fits and its model function
# returns a 4-vector. Adding stats there would change every existing fit's
# summary length. Keeping the extended scheme alongside lets the two run side by
# side (exactly the "two fitting schemes without duplicated plumbing" split the
# README describes), and a downstream consumer can pick whichever it wants.
#
# DESIGN NOTES
#  * Weekly binning. Deaths are binned into fixed-width (default 7-day) bins on
#    ABSOLUTE outbreak time (day 0 = seeding), matching the trajectory binning in
#    the run scripts (bin_counts() / BIN_WIDTH). peak_height = the largest weekly
#    count; the peak week is the FIRST week attaining that maximum (which.max),
#    a deterministic tie-break.
#  * time_to_peak origin. Controlled by `time_origin`:
#      "first_death"   (default) -> whole weeks from the FIRST-death week to the
#                       peak week. Anchored on the first death, like `duration`
#                       in abc_summarise(), so it tracks the surveillance-
#                       observable curve and is robust to the takeoff/seeding lag.
#                       Always >= 0 (deaths cannot peak before the first death).
#      "outbreak_start" -> midpoint day of the peak week measured from t = 0.
#  * Aggregation across replicates mirrors the base scheme EXACTLY: each
#    replicate is summarised on its OWN curve, then the per-replicate summaries
#    are averaged across the replicates that took off (NOT a peak of the averaged
#    curve, which would be biased low by between-replicate misalignment).
#  * Scale. time_to_peak (~hundreds of days) and peak_height (~tens) sit on very
#    different scales from n_deaths (~thousands); ABC_sequential() normalises each
#    summary by its step-1 SD, so no manual weighting is needed here.
#
# Contents (peak-scheme specific; everything else is reused from the NPI file):
#   DEFAULT_PEAK_BIN_WIDTH / DEFAULT_PEAK_TIME_ORIGIN : binning defaults
#   weekly_incidence_counts()      : vector of event day-indices -> weekly counts
#   peak_stats_from_death_days()    : weekly death days -> (time_to_peak, peak_height)
#   abc_summarise_peak()            : per-replicate 6-stat summary
#   run_one_abc_replicate_peak()    : do.call(branching_process_main) + summary
#   fiber_abc_model_peak()          : in-process single-core ABC model (6 stats)
#   bootstrap_abc_worker_peak()     : per-worker setup (also sources the NPI file)
#   fiber_abc_model_parallel_peak() : the function ABC_sequential() calls
#   prior_predictive_check_peak()   : draw from priors, run the non-parallel model
#
# Requires setup_model_parameters.R, calculate_model_approx_r0.R,
# abc_calibration_functions_common.R AND abc_calibration_functions_npi.R sourced
# first (for build_abc_model_args(), npi_efficacy_from_scaler(), the make/attach/
# write_npi_run_metadata() helpers, and DEFAULT_NPI_SPEC), plus library(fiber).
# The PSOCK-worker self-bootstrap below sources the NPI file for you.
# -----------------------------------------------------------------------------

# Binning defaults. Weekly bins, anchored on the first death (see header).
DEFAULT_PEAK_BIN_WIDTH   <- 7L
DEFAULT_PEAK_TIME_ORIGIN <- "first_death"   # or "outbreak_start"


# Bin a vector of (0-based) event day-indices into fixed-width bins and return
# the integer count in each bin from the first bin up to the bin holding the
# latest event. Pure base R (safe to ship to PSOCK workers without
# simulation_helpers.R). Non-finite / negative days are dropped / clamped.
weekly_incidence_counts <- function(day_index, bin_width = DEFAULT_PEAK_BIN_WIDTH) {
  day_index <- day_index[is.finite(day_index)]
  if (length(day_index) == 0L) return(integer(0))
  day_index <- as.integer(day_index)
  day_index[day_index < 0L] <- 0L
  bins   <- (day_index %/% bin_width) + 1L        # 1-based bin index
  n_bins <- max(bins)
  tabulate(bins, nbins = n_bins)
}


# Peak-of-death-curve summary for ONE replicate. `death_days` are the (floored,
# 0-based, absolute) day-indices of every death; `first_death_day` anchors the
# "first_death" origin. Returns c(time_to_peak, peak_height); both 0 when there
# are no deaths.
peak_stats_from_death_days <- function(death_days,
                                       first_death_day,
                                       bin_width   = DEFAULT_PEAK_BIN_WIDTH,
                                       time_origin = DEFAULT_PEAK_TIME_ORIGIN) {
  weekly <- weekly_incidence_counts(death_days, bin_width)
  if (length(weekly) == 0L || all(weekly == 0L)) {
    return(c(time_to_peak = 0, peak_height = 0))
  }

  peak_bin    <- which.max(weekly)               # FIRST week attaining the max
  peak_height <- as.numeric(weekly[peak_bin])

  time_to_peak <- if (identical(time_origin, "outbreak_start")) {
    (peak_bin - 0.5) * bin_width                 # midpoint of peak week, from t=0
  } else {
    # whole weeks from the first-death week to the peak week (>= 0).
    fd_bin <- (as.integer(floor(first_death_day)) %/% bin_width) + 1L
    max((peak_bin - fd_bin) * bin_width, 0)
  }

  c(time_to_peak = as.numeric(time_to_peak), peak_height = peak_height)
}


# Per-replicate summary: the four base summaries (see abc_summarise() in the
# common file) PLUS time_to_peak and peak_height from the weekly death curve.
# Returns a named numeric(6) in the canonical order. Mirrors abc_summarise()'s
# definitions of n_cases / n_deaths / n_hcw_deaths / duration exactly so the two
# schemes agree on the shared four.
abc_summarise_peak <- function(out,
                               bin_width   = DEFAULT_PEAK_BIN_WIDTH,
                               time_origin = DEFAULT_PEAK_TIME_ORIGIN) {
  empty <- c(n_cases = 0, n_deaths = 0, n_hcw_deaths = 0, duration = 0,
             time_to_peak = 0, peak_height = 0)

  tdf <- out$tdf
  tdf <- tdf[!is.na(tdf$time_infection_absolute), , drop = FALSE]
  if (nrow(tdf) == 0L) return(empty)

  deaths <- !is.na(tdf$outcome) & tdf$outcome
  hcw    <- !is.na(tdf$class)   & tdf$class == "HCW"

  n_cases      <- nrow(tdf)
  n_deaths     <- sum(deaths)
  n_hcw_deaths <- sum(deaths & hcw)

  if (n_deaths == 0L) {
    empty["n_cases"]      <- n_cases
    empty["n_hcw_deaths"] <- n_hcw_deaths
    return(empty)
  }

  death_t       <- tdf$time_outcome_absolute[deaths]
  death_t       <- death_t[!is.na(death_t)]
  t_first_death <- min(death_t)
  t_last_event  <- max(tdf$time_outcome_absolute, na.rm = TRUE)
  duration      <- max(t_last_event - t_first_death, 0)

  peak <- peak_stats_from_death_days(
    death_days      = as.integer(floor(death_t)),
    first_death_day = t_first_death,
    bin_width       = bin_width,
    time_origin     = time_origin
  )

  c(n_cases      = n_cases,
    n_deaths     = n_deaths,
    n_hcw_deaths = n_hcw_deaths,
    duration     = duration,
    time_to_peak = unname(peak[["time_to_peak"]]),
    peak_height  = unname(peak[["peak_height"]]))
}


run_one_abc_replicate_peak <- function(args,
                                       bin_width   = DEFAULT_PEAK_BIN_WIDTH,
                                       time_origin = DEFAULT_PEAK_TIME_ORIGIN) {
  out <- do.call(branching_process_main, args)
  abc_summarise_peak(out, bin_width = bin_width, time_origin = time_origin)
}


# In-process (single-core) ABC model for the extended scheme. Same contract as
# fiber_abc_model() in the NPI file but returns the 6-stat vector. Used by
# prior_predictive_check_peak(); the parallel worker path uses
# fiber_abc_model_parallel_peak().
fiber_abc_model_peak <- function(theta,
                                 base,
                                 tv,
                                 invariants,
                                 npi_spec = DEFAULT_NPI_SPEC,
                                 general_hospital_quarantine_efficacy,
                                 safe_funeral_efficacy,
                                 n_replicates = 30,
                                 seeding_cases = 25,
                                 takeoff_death_threshold = 100,
                                 bin_width   = DEFAULT_PEAK_BIN_WIDTH,
                                 time_origin = DEFAULT_PEAK_TIME_ORIGIN,
                                 include_n_cases = FALSE) {

  R0           <- theta[1]
  prop_funeral <- theta[2]
  npi_scaler   <- theta[3]

  args <- build_abc_model_args(
    R0 = R0, prop_funeral = prop_funeral, npi_scaler = npi_scaler,
    base = base, tv = tv, invariants = invariants, npi_spec = npi_spec,
    general_hospital_quarantine_efficacy = general_hospital_quarantine_efficacy,
    safe_funeral_efficacy = safe_funeral_efficacy,
    seeding_cases = seeding_cases
  )

  reps <- vapply(
    seq_len(n_replicates),
    function(i) run_one_abc_replicate_peak(args, bin_width = bin_width,
                                           time_origin = time_origin),
    numeric(6)
  )
  rownames(reps) <- c("n_cases", "n_deaths", "n_hcw_deaths", "duration",
                      "time_to_peak", "peak_height")

  took_off     <- reps["n_deaths", ] >= takeoff_death_threshold
  takeoff_frac <- mean(took_off)

  if (!any(took_off)) {
    out <- c(takeoff = 0, n_deaths = 0, n_hcw_deaths = 0, duration = 0,
             time_to_peak = 0, peak_height = 0)
    if (include_n_cases) out <- c(out, n_cases = 0)
    return(out)
  }

  out <- c(takeoff      = takeoff_frac,
           n_deaths     = mean(reps["n_deaths",     took_off]),
           n_hcw_deaths = mean(reps["n_hcw_deaths", took_off]),
           duration     = mean(reps["duration",     took_off]),
           time_to_peak = mean(reps["time_to_peak", took_off]),
           peak_height  = mean(reps["peak_height",  took_off]))

  if (include_n_cases) {
    out <- c(out, n_cases = mean(reps["n_cases", took_off]))
  }
  out
}


# Per-worker setup for the extended scheme. Identical to bootstrap_abc_worker()
# in the NPI file EXCEPT it (a) also sources the NPI approach file that lives
# alongside this one -- the peak scheme reuses build_abc_model_args() /
# npi_efficacy_from_scaler() / DEFAULT_NPI_SPEC from there -- and (b) carries the
# weekly-binning settings (peak_bin_width, peak_time_origin) in the config.
bootstrap_abc_worker_peak <- function(setup_path,
                                      functions_path,
                                      r0_path,
                                      scenario_csv,
                                      scenario_id,
                                      abc_config = list(),
                                      model_overrides = list()) {
  # Source helpers FIRST so DEFAULT_SCALAR_INPUTS / DEFAULT_NPI_SPEC exist before
  # default_config (below) references them. functions_path is THIS (peak) file;
  # the NPI file sits next to it (dirname()), as does the common generics file.
  source(setup_path)
  source(file.path(dirname(functions_path), "abc_calibration_functions_npi.R"))
  source(functions_path)
  source(r0_path)
  source(file.path(dirname(functions_path), "abc_calibration_functions_common.R"))

  if (!"fiber" %in% loadedNamespaces()) {
    library(fiber)
  }

  default_config <- list(
    check_final_size        = 30000,
    takeoff_death_threshold = 100,
    n_reps                  = 30,
    seeding_cases           = 25,
    setup_R0_n              = 100000L,
    setup_R0_seed           = 42L,
    npi_spec                = DEFAULT_NPI_SPEC,
    general_hospital_quarantine_efficacy =
      DEFAULT_SCALAR_INPUTS$general_hospital_quarantine_efficacy,
    safe_funeral_efficacy   = DEFAULT_SCALAR_INPUTS$safe_funeral_efficacy,
    peak_bin_width          = DEFAULT_PEAK_BIN_WIDTH,
    peak_time_origin        = DEFAULT_PEAK_TIME_ORIGIN
  )
  abc_config <- utils::modifyList(default_config, abc_config)

  scenario_matrix <- read_scenario_matrix(scenario_csv)

  scalar_overrides <- model_overrides
  if (!"check_final_size" %in% names(scalar_overrides)) {
    scalar_overrides$check_final_size <- abc_config$check_final_size
  }

  mp <- make_model_parameters(
    scenario_id = scenario_id,
    scenario_matrix = scenario_matrix,
    overrides = scalar_overrides
  )

  invariants <- compute_R0_invariants(
    args = mp$args,
    n    = abc_config$setup_R0_n,
    seed = abc_config$setup_R0_seed
  )

  assign("base_args",      mp$base_args,                          envir = globalenv())
  assign("tv_args_model",  mp$tv_args,                            envir = globalenv())
  assign("R0_invariants",  invariants,                            envir = globalenv())
  assign(".npi_spec",      abc_config$npi_spec,                   envir = globalenv())
  assign(".fixed_general_hosp",
         abc_config$general_hospital_quarantine_efficacy,         envir = globalenv())
  assign(".fixed_safe_funeral",
         abc_config$safe_funeral_efficacy,                        envir = globalenv())
  assign(".abc_config",    abc_config,                            envir = globalenv())
  assign(".fiber_abc_ready", TRUE,                                envir = globalenv())

  invisible(list(
    invariants = invariants,
    abc_config = abc_config,
    scenario_label = mp$scenario_label
  ))
}


# The function ABC_sequential() calls on each PSOCK worker. Mirrors
# fiber_abc_model_parallel() in the NPI file but summarises with
# abc_summarise_peak() and returns the 6-stat vector. Self-bootstraps on first
# call (reads FIBER_ABC_CONFIG, sources this file, then
# bootstrap_abc_worker_peak()).
fiber_abc_model_parallel_peak <- function(theta_with_seed) {
  if (!isTRUE(get0(".fiber_abc_ready", envir = globalenv()))) {
    cfg_path <- Sys.getenv("FIBER_ABC_CONFIG")
    if (cfg_path == "" || !file.exists(cfg_path)) {
      stop("FIBER_ABC_CONFIG env var not set or file missing. ",
           "Call save_abc_config(<config>) in the main process before ",
           "running ABC_sequential().", call. = FALSE)
    }
    cfg <- readRDS(cfg_path)
    source(cfg$functions_path)

    abc_config_arg <- if (is.null(cfg$abc_config)) list() else cfg$abc_config
    overrides_arg  <- if (is.null(cfg$model_overrides)) list() else cfg$model_overrides

    bootstrap_abc_worker_peak(
      setup_path      = cfg$setup_path,
      functions_path  = cfg$functions_path,
      r0_path         = cfg$r0_path,
      scenario_csv    = cfg$scenario_csv,
      scenario_id     = cfg$scenario_id,
      abc_config      = abc_config_arg,
      model_overrides = overrides_arg
    )
  }

  cfg_run <- get(".abc_config", envir = globalenv())

  set.seed(theta_with_seed[1])

  args <- build_abc_model_args(
    R0              = theta_with_seed[2],
    prop_funeral    = theta_with_seed[3],
    npi_scaler      = theta_with_seed[4],
    base            = get("base_args",     envir = globalenv()),
    tv              = get("tv_args_model", envir = globalenv()),
    invariants      = get("R0_invariants", envir = globalenv()),
    npi_spec        = get(".npi_spec",     envir = globalenv()),
    general_hospital_quarantine_efficacy = get(".fixed_general_hosp", envir = globalenv()),
    safe_funeral_efficacy                = get(".fixed_safe_funeral", envir = globalenv()),
    seeding_cases   = cfg_run$seeding_cases
  )

  bin_width <- cfg_run$peak_bin_width
  if (is.null(bin_width)) bin_width <- DEFAULT_PEAK_BIN_WIDTH
  time_origin <- cfg_run$peak_time_origin
  if (is.null(time_origin)) time_origin <- DEFAULT_PEAK_TIME_ORIGIN

  reps <- vapply(seq_len(cfg_run$n_reps), function(i) {
    out <- do.call(branching_process_main, args)
    abc_summarise_peak(out, bin_width = bin_width, time_origin = time_origin)
  }, numeric(6))
  rownames(reps) <- c("n_cases", "n_deaths", "n_hcw_deaths", "duration",
                      "time_to_peak", "peak_height")

  took_off <- reps["n_deaths", ] >= cfg_run$takeoff_death_threshold
  if (!any(took_off)) {
    return(c(takeoff = 0, n_deaths = 0, n_hcw_deaths = 0, duration = 0,
             time_to_peak = 0, peak_height = 0))
  }
  c(takeoff      = mean(took_off),
    n_deaths     = mean(reps["n_deaths",     took_off]),
    n_hcw_deaths = mean(reps["n_hcw_deaths", took_off]),
    duration     = mean(reps["duration",     took_off]),
    time_to_peak = mean(reps["time_to_peak", took_off]),
    peak_height  = mean(reps["peak_height",  took_off]))
}


# Prior predictive check for the extended scheme. Same as prior_predictive_check()
# in the NPI file but runs fiber_abc_model_peak() (6 stats + optional n_cases).
prior_predictive_check_peak <- function(n_draws,
                                        prior_list,
                                        base,
                                        tv,
                                        invariants,
                                        npi_spec = DEFAULT_NPI_SPEC,
                                        general_hospital_quarantine_efficacy,
                                        safe_funeral_efficacy,
                                        param_names = c("R0", "prop_funeral", "npi_scaler"),
                                        parallel = FALSE,
                                        n_replicates = 30,
                                        seeding_cases = 25,
                                        takeoff_death_threshold = 100,
                                        bin_width   = DEFAULT_PEAK_BIN_WIDTH,
                                        time_origin = DEFAULT_PEAK_TIME_ORIGIN,
                                        include_n_cases = TRUE) {

  sample_one <- function(spec, n) {
    dist   <- spec[1]
    params <- as.numeric(spec[-1])
    switch(dist,
           "unif"        = runif(n,  min = params[1], max = params[2]),
           "normal"      = rnorm(n,  mean = params[1], sd = params[2]),
           "lognormal"   = rlnorm(n, meanlog = params[1], sdlog = params[2]),
           "exponential" = rexp(n,   rate = params[1]),
           stop("Unsupported prior distribution: ", dist)
    )
  }

  draws <- as.data.frame(lapply(prior_list, sample_one, n = n_draws))
  names(draws) <- param_names

  theta_list <- lapply(seq_len(nrow(draws)), function(i) as.numeric(draws[i, ]))

  run_one_theta <- function(theta) {
    fiber_abc_model_peak(
      theta = theta,
      base = base, tv = tv, invariants = invariants, npi_spec = npi_spec,
      general_hospital_quarantine_efficacy = general_hospital_quarantine_efficacy,
      safe_funeral_efficacy = safe_funeral_efficacy,
      n_replicates = n_replicates,
      seeding_cases = seeding_cases,
      takeoff_death_threshold = takeoff_death_threshold,
      bin_width = bin_width, time_origin = time_origin,
      include_n_cases = include_n_cases
    )
  }

  sims_list <- if (parallel) {
    if (!requireNamespace("future.apply", quietly = TRUE)) {
      stop("Package 'future.apply' is required when parallel = TRUE.",
           call. = FALSE)
    }
    if (requireNamespace("progressr", quietly = TRUE)) {
      progressr::with_progress({
        p <- progressr::progressor(steps = length(theta_list))
        fn <- function(theta) { out <- run_one_theta(theta); p(); out }
        future.apply::future_lapply(theta_list, fn, future.seed = TRUE)
      })
    } else {
      future.apply::future_lapply(theta_list, run_one_theta, future.seed = TRUE)
    }
  } else {
    if (requireNamespace("progressr", quietly = TRUE)) {
      progressr::with_progress({
        p <- progressr::progressor(steps = length(theta_list))
        lapply(theta_list, function(theta) { out <- run_one_theta(theta); p(); out })
      })
    } else {
      lapply(theta_list, run_one_theta)
    }
  }

  sims <- do.call(rbind, sims_list)
  cbind(draws, as.data.frame(sims))
}
