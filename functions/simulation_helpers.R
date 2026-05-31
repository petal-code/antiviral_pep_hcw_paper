# simulation_helpers.R
# -----------------------------------------------------------------------------
# Run individual fiber replicates downstream of a fit and summarise their
# outputs. Used by analyses that take pre-built model argument lists (e.g. from
# build_abc_model_args()) and run many stochastic replicates, optionally across
# intervention arms.
#
# Contents
#   `%||%`         : null-coalescing operator.
#   simulate_one() : run ONE fiber replicate (one parameter set, one stochastic
#                    rep, one arm) and return compact per-replicate outputs,
#                    including the day-indices of every death and HCW death so
#                    epidemic curves can be reconstructed.
#   bin_counts()   : bin a vector of event day-indices into fixed-width bins.
#   q_summary()    : mean + chosen quantiles of a vector.
#
# simulate_one() only uses fiber + base R, so it is safe to ship to
# future/PSOCK workers (load fiber on the worker via future.packages = "fiber").
# -----------------------------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a


# Run a single fiber replicate. `job` is a one-row list describing what to run;
# `args_list` is the list of pre-built fiber argument lists (one per parameter
# set, produced in the main session with build_abc_model_args()); `obv_cfg`
# carries the obeldesivir settings for the OBV arm.
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

# Mean plus an arbitrary set of quantiles of a numeric vector, returned as a
# named numeric vector. Non-finite values are dropped.
q_summary <- function(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(setNames(rep(NA_real_, length(probs) + 1L), c("mean", paste0("q", probs * 100))))
  c(mean = mean(x), setNames(stats::quantile(x, probs = probs, names = FALSE), paste0("q", probs * 100)))
}
