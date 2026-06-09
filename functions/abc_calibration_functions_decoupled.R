# abc_calibration_functions_decoupled.R
# -----------------------------------------------------------------------------
# DECOUPLED-EFFICACY ABC approach. A self-contained scheme that:
#   (1) fits etu_efficacy and ppe_efficacy DIRECTLY and INDEPENDENTLY (no shared
#       npi_scaler), so chasing HCW deaths (low PPE) no longer drags ETU -- and
#       hence overall epidemic size -- down with it;
#   (2) adds hcw_risk_scalar (scales the symmetric prob_hcw_cond_*_hospital base)
#       as the HCW-exposure lever;
#   (3) fits a TRANSFORMED summary set: log(n_deaths), log(n_hcw_deaths),
#       hcw_fraction (ratio of means), d_p05_p95 (5-95% span of death dates) and
#       log(peak_height). The logs put the heavy-tailed counts on a relative-error
#       scale; hcw_fraction up-weights the HCW share; d_p05_p95 is a low-noise,
#       tail-robust timing summary that still carries the prop_funeral signal.
#
# Canonical fittable parameters (DECOUPLED_PARAM_NAMES), in this order:
#   R0              : baseline reproduction number for a genPop seeding case (t=0).
#   prop_funeral    : share of R0 attributable to funeral transmission at t = 0.
#   etu_efficacy    : ETU hospital-quarantine efficacy (enters D / offspring means).
#   ppe_efficacy    : PPE efficacy protecting HCWs (paired with ppe_coverage_hcw(t)).
#   hcw_risk_scalar : prob_hcw_cond_*_hospital = min(hcw_base_prob * scalar, 1).
# Any SUBSET may be fitted; the rest are held at fixed values (per-run FIXED_PARAMS,
# falling back to DECOUPLED_PARAM_DEFAULTS).
#
# IDENTIFIABILITY NOTE: with only aggregate (not time-resolved) HCW data, PPE and
# hcw_risk_scalar both act on the HCW-share dimension and form a (prior-bounded)
# ridge -- fit both, but report (ppe, hcw_risk) JOINTLY, and bound hcw_risk via
# its prior. See pairs() of the posterior in the run scripts.
#
# REUSES (does NOT modify) the shared infrastructure:
#   * abc_summarise() + the disk-inspection trio (abc_progress / abc_compare_steps
#     / reconstruct_abc_result) + make_abc_output_dir() / save_abc_config() /
#     with_abc_output_dir() from abc_calibration_functions_common.R
#   * D_from_invariants(), F_from_invariants(), solve_offspring_means(),
#     compute_R0_invariants() from calculate_model_approx_r0.R
#   * make_model_parameters(), read_scenario_matrix(), DEFAULT_SCALAR_INPUTS from
#     setup_model_parameters.R
# The summary aggregation (logs / hcw_fraction / d_p05_p95) is NEW and lives here
# (run_abc_particle_decoupled / aggregate_decoupled), so the shared common file is
# left untouched.
#
# Approach-specific functions ONLY. Required source order (run script AND worker
# bootstrap):
#   setup_model_parameters.R -> abc_calibration_functions_common.R ->
#   THIS FILE -> calculate_model_approx_r0.R   (and library(fiber)).
# bootstrap_abc_worker_decoupled() derives the common path from
# dirname(functions_path), exactly like the other schemes.
# -----------------------------------------------------------------------------


# Canonical fittable parameters and the order ABC stores them.
DECOUPLED_PARAM_NAMES <- c("R0", "prop_funeral", "etu_efficacy",
                           "ppe_efficacy", "hcw_risk_scalar")

# Values used for any canonical parameter NOT being fitted (override per run via
# FIXED_PARAMS). hcw_risk_scalar = 1 keeps prob_hcw at the base; the efficacies
# default to plausible mid-range values.
DECOUPLED_PARAM_DEFAULTS <- list(
  R0              = 1.45,
  prop_funeral    = 0.30,
  etu_efficacy    = 0.85,
  ppe_efficacy    = 0.65,
  hcw_risk_scalar = 1.5
)

# Summaries this scheme can return (any subset may be fitted). The default set is
# the one agreed for the decoupled fit; raw / extra metrics are available too.
DECOUPLED_AVAILABLE_SUMMARIES <- c(
  "takeoff",
  "n_deaths", "log_n_deaths",
  "n_hcw_deaths", "log_n_hcw_deaths",
  "peak_height", "log_peak_height",
  "hcw_fraction",
  "d_p05_p95",
  "duration", "time_to_peak"
)
DEFAULT_DECOUPLED_SUMMARIES <- c("log_n_deaths", "log_n_hcw_deaths",
                                 "hcw_fraction", "d_p05_p95", "log_peak_height")

DECOUPLED_PEAK_BIN_WIDTH   <- 7L
DECOUPLED_PEAK_TIME_ORIGIN <- "first_death"


# -----------------------------------------------------------------------------
# Per-replicate metrics + aggregation (NEW: logs / hcw_fraction / d_p05_p95).
# -----------------------------------------------------------------------------

# Per-replicate raw metrics for ONE branching_process_main() run: the standard
# counts/curve features via the shared abc_summarise(), plus d_p05_p95 (the day
# span between the 5th and 95th percentile of death dates -- a tail-robust,
# low-noise "duration"). Returns a named numeric vector of length 6.
per_rep_metrics_decoupled <- function(out,
                                      bin_width   = DECOUPLED_PEAK_BIN_WIDTH,
                                      time_origin = DECOUPLED_PEAK_TIME_ORIGIN,
                                      lo = 0.05, hi = 0.95) {
  b <- abc_summarise(
    out,
    metrics     = c("n_deaths", "n_hcw_deaths", "peak_height", "duration", "time_to_peak"),
    bin_width   = bin_width,
    time_origin = time_origin
  )
  tdf <- out$tdf
  tdf <- tdf[!is.na(tdf$time_infection_absolute), , drop = FALSE]
  dt  <- sort(tdf$time_outcome_absolute[!is.na(tdf$outcome) & tdf$outcome])
  n   <- length(dt)
  d_p05_p95 <- if (n >= 1L) dt[ceiling(hi * n)] - dt[ceiling(lo * n)] else 0
  c(b, d_p05_p95 = d_p05_p95)
}


# Aggregate a [metric x replicate] matrix into the requested fitted summaries.
# Means are over TAKEN-OFF replicates (the `took_off` mask = n_deaths >=
# takeoff_death_threshold); hcw_fraction is the ratio of the two means (stable);
# "log_*" summaries log the corresponding mean. "takeoff" is the fraction of
# replicates that took off, scored on its OWN threshold
# (takeoff_fraction_threshold, default 1 death) so it stays DECOUPLED from the
# conditioning `took_off` mask -- it reports a low-bar "did it seed at all?"
# fraction while the means remain conditioned on real (>= takeoff_death_threshold)
# outbreaks. Returns a named vector in `summary_stats` order.
aggregate_decoupled <- function(reps, took_off,
                                summary_stats = DEFAULT_DECOUPLED_SUMMARIES,
                                takeoff_fraction_threshold = 1) {
  any_off <- any(took_off)
  m  <- function(metric) if (any_off) mean(reps[metric, took_off]) else 0
  md <- m("n_deaths"); mh <- m("n_hcw_deaths"); mp <- m("peak_height")
  slog <- function(x) log(pmax(x, 1e-6))

  val <- function(s) switch(
    s,
    takeoff          = mean(reps["n_deaths", ] >= takeoff_fraction_threshold),
    n_deaths         = md,
    log_n_deaths     = slog(md),
    n_hcw_deaths     = mh,
    log_n_hcw_deaths = slog(mh),
    peak_height      = mp,
    log_peak_height  = slog(mp),
    hcw_fraction     = if (md > 0) mh / md else 0,
    d_p05_p95        = m("d_p05_p95"),
    duration         = m("duration"),
    time_to_peak     = m("time_to_peak"),
    stop("Unknown summary statistic: ", s, ". Available: ",
         paste(DECOUPLED_AVAILABLE_SUMMARIES, collapse = ", "), ".", call. = FALSE)
  )

  out <- vapply(summary_stats, val, numeric(1))
  names(out) <- summary_stats
  out
}


# Run ONE ABC particle: n_reps stochastic replicates of `args`, summarise each,
# then aggregate to the fitted (transformed) summary vector. The decoupled analogue
# of the common run_abc_particle(); kept here so the shared file is untouched.
run_abc_particle_decoupled <- function(args,
                                       n_reps,
                                       takeoff_death_threshold,
                                       summary_stats = DEFAULT_DECOUPLED_SUMMARIES,
                                       bin_width     = DECOUPLED_PEAK_BIN_WIDTH,
                                       time_origin   = DECOUPLED_PEAK_TIME_ORIGIN,
                                       takeoff_fraction_threshold = 1) {
  unknown <- setdiff(summary_stats, DECOUPLED_AVAILABLE_SUMMARIES)
  if (length(unknown) > 0L) {
    stop("Unknown summary_stats: ", paste(unknown, collapse = ", "),
         ". Available: ", paste(DECOUPLED_AVAILABLE_SUMMARIES, collapse = ", "),
         ".", call. = FALSE)
  }

  reps <- vapply(seq_len(n_reps), function(i) {
    out <- do.call(branching_process_main, args)
    per_rep_metrics_decoupled(out, bin_width = bin_width, time_origin = time_origin)
  }, numeric(6))
  rownames(reps) <- c("n_deaths", "n_hcw_deaths", "peak_height",
                      "duration", "time_to_peak", "d_p05_p95")

  took_off <- reps["n_deaths", ] >= takeoff_death_threshold
  aggregate_decoupled(reps, took_off, summary_stats,
                      takeoff_fraction_threshold = takeoff_fraction_threshold)
}


# -----------------------------------------------------------------------------
# Parameter assembly + model-args construction (decoupled efficacies).
# -----------------------------------------------------------------------------

# Splice the fitted subset `theta_fitted` (in `fit_params` order) onto the
# complete `fixed_values` list (all DECOUPLED_PARAM_NAMES).
assemble_decoupled_theta <- function(theta_fitted, fit_params, fixed_values) {
  if (length(theta_fitted) != length(fit_params)) {
    stop("length(theta_fitted) (", length(theta_fitted), ") must equal ",
         "length(fit_params) (", length(fit_params), ").", call. = FALSE)
  }
  full <- fixed_values[DECOUPLED_PARAM_NAMES]
  names(full) <- DECOUPLED_PARAM_NAMES
  for (i in seq_along(fit_params)) full[[fit_params[i]]] <- theta_fitted[[i]]
  full
}


# (R0, prop_funeral, etu_efficacy, ppe_efficacy, hcw_risk_scalar) -> fiber args.
# etu_efficacy enters D (per-particle, via the cached invariants); ppe_efficacy is
# applied directly (it only thins HCW recipients, so it does NOT enter D in the
# genPop-dominant R0 approximation); hcw_risk_scalar scales hospital HCW exposure.
# Efficacies are clamped to [0, 1] so informative (e.g. normal) priors whose tails
# stray out of range never break the model.
build_abc_model_args_decoupled <- function(R0,
                                           prop_funeral,
                                           etu_efficacy,
                                           ppe_efficacy,
                                           hcw_risk_scalar,
                                           base,
                                           tv,
                                           invariants,
                                           general_hospital_quarantine_efficacy,
                                           safe_funeral_efficacy,
                                           hcw_base_prob = 0.25,
                                           seeding_cases = 25) {
  etu_efficacy <- min(max(etu_efficacy, 0), 1)
  ppe_efficacy <- min(max(ppe_efficacy, 0), 1)

  D <- D_from_invariants(
    invariants,
    etu_efficacy = etu_efficacy,
    general_hospital_quarantine_efficacy = general_hospital_quarantine_efficacy
  )
  F_fun <- F_from_invariants(invariants, safe_funeral_efficacy)
  means <- solve_offspring_means(R0, prop_funeral, D, F_fun)

  hcw_hospital <- pmin(hcw_base_prob * hcw_risk_scalar, 1.0)

  args <- c(base, tv)
  args$mn_offspring_genPop  <- means$mn_genPop
  args$mn_offspring_funeral <- means$mn_funeral

  args$ppe_efficacy <- ppe_efficacy
  args$etu_efficacy <- etu_efficacy
  args$general_hospital_quarantine_efficacy <- general_hospital_quarantine_efficacy
  args$safe_funeral_efficacy <- safe_funeral_efficacy

  args$prob_hcw_cond_genPop_hospital <- hcw_hospital
  args$prob_hcw_cond_hcw_hospital    <- hcw_hospital

  args$seed          <- NULL
  args$seeding_cases <- seeding_cases
  args
}


# In-process single-core ABC model (prior_predictive_check_decoupled() + testing).
fiber_abc_model_decoupled <- function(theta,
                                      fit_params,
                                      fixed_values,
                                      base,
                                      tv,
                                      invariants,
                                      general_hospital_quarantine_efficacy,
                                      safe_funeral_efficacy,
                                      hcw_base_prob = 0.25,
                                      n_replicates = 30,
                                      seeding_cases = 25,
                                      takeoff_death_threshold = 100,
                                      summary_stats = DEFAULT_DECOUPLED_SUMMARIES,
                                      bin_width     = DECOUPLED_PEAK_BIN_WIDTH,
                                      time_origin   = DECOUPLED_PEAK_TIME_ORIGIN,
                                      takeoff_fraction_threshold = 1) {
  full <- assemble_decoupled_theta(theta, fit_params, fixed_values)
  args <- build_abc_model_args_decoupled(
    R0 = full$R0, prop_funeral = full$prop_funeral,
    etu_efficacy = full$etu_efficacy, ppe_efficacy = full$ppe_efficacy,
    hcw_risk_scalar = full$hcw_risk_scalar,
    base = base, tv = tv, invariants = invariants,
    general_hospital_quarantine_efficacy = general_hospital_quarantine_efficacy,
    safe_funeral_efficacy = safe_funeral_efficacy,
    hcw_base_prob = hcw_base_prob, seeding_cases = seeding_cases
  )
  run_abc_particle_decoupled(
    args, n_reps = n_replicates, takeoff_death_threshold = takeoff_death_threshold,
    summary_stats = summary_stats, bin_width = bin_width, time_origin = time_origin,
    takeoff_fraction_threshold = takeoff_fraction_threshold
  )
}


# -----------------------------------------------------------------------------
# Per-worker bootstrap + the function ABC_sequential() calls.
# -----------------------------------------------------------------------------
bootstrap_abc_worker_decoupled <- function(setup_path,
                                           functions_path,
                                           r0_path,
                                           scenario_csv,
                                           scenario_id,
                                           abc_config = list(),
                                           model_overrides = list()) {
  fn_dir <- dirname(functions_path)
  source(setup_path)
  source(functions_path)
  source(r0_path)
  source(file.path(fn_dir, "abc_calibration_functions_common.R"))

  if (!"fiber" %in% loadedNamespaces()) {
    library(fiber)
  }

  default_config <- list(
    check_final_size           = 30000,
    takeoff_death_threshold    = 100,
    takeoff_fraction_threshold = 1,
    n_reps                  = 30,
    seeding_cases           = 25,
    setup_R0_n              = 100000L,
    setup_R0_seed           = 42L,
    hcw_base_prob           = 0.25,
    general_hospital_quarantine_efficacy =
      DEFAULT_SCALAR_INPUTS$general_hospital_quarantine_efficacy,
    safe_funeral_efficacy   = DEFAULT_SCALAR_INPUTS$safe_funeral_efficacy,
    fit_params              = DECOUPLED_PARAM_NAMES,
    fixed_values            = DECOUPLED_PARAM_DEFAULTS,
    summary_stats           = DEFAULT_DECOUPLED_SUMMARIES,
    peak_bin_width          = DECOUPLED_PEAK_BIN_WIDTH,
    peak_time_origin        = DECOUPLED_PEAK_TIME_ORIGIN
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

  assign("base_args",      mp$base_args,            envir = globalenv())
  assign("tv_args_model",  mp$tv_args,              envir = globalenv())
  assign("R0_invariants",  invariants,              envir = globalenv())
  assign(".hcw_base_prob", abc_config$hcw_base_prob, envir = globalenv())
  assign(".fit_params",    abc_config$fit_params,   envir = globalenv())
  assign(".fixed_values",  abc_config$fixed_values, envir = globalenv())
  assign(".fixed_general_hosp",
         abc_config$general_hospital_quarantine_efficacy, envir = globalenv())
  assign(".fixed_safe_funeral",
         abc_config$safe_funeral_efficacy,          envir = globalenv())
  assign(".abc_config",    abc_config,              envir = globalenv())
  assign(".fiber_abc_ready", TRUE,                  envir = globalenv())

  invisible(list(
    invariants = invariants,
    abc_config = abc_config,
    scenario_label = mp$scenario_label
  ))
}


fiber_abc_model_parallel_decoupled <- function(theta_with_seed) {
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

    bootstrap_abc_worker_decoupled(
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
  theta_fitted <- theta_with_seed[-1]

  full <- assemble_decoupled_theta(
    theta_fitted, cfg_run$fit_params, cfg_run$fixed_values
  )

  args <- build_abc_model_args_decoupled(
    R0              = full$R0,
    prop_funeral    = full$prop_funeral,
    etu_efficacy    = full$etu_efficacy,
    ppe_efficacy    = full$ppe_efficacy,
    hcw_risk_scalar = full$hcw_risk_scalar,
    base            = get("base_args",     envir = globalenv()),
    tv              = get("tv_args_model", envir = globalenv()),
    invariants      = get("R0_invariants", envir = globalenv()),
    general_hospital_quarantine_efficacy = get(".fixed_general_hosp", envir = globalenv()),
    safe_funeral_efficacy                = get(".fixed_safe_funeral", envir = globalenv()),
    hcw_base_prob   = get(".hcw_base_prob", envir = globalenv()),
    seeding_cases   = cfg_run$seeding_cases
  )

  run_abc_particle_decoupled(
    args,
    n_reps                  = cfg_run$n_reps,
    takeoff_death_threshold = cfg_run$takeoff_death_threshold,
    summary_stats           = cfg_run$summary_stats,
    bin_width               = cfg_run$peak_bin_width,
    time_origin             = cfg_run$peak_time_origin,
    takeoff_fraction_threshold =
      if (is.null(cfg_run$takeoff_fraction_threshold)) 1 else cfg_run$takeoff_fraction_threshold
  )
}


# -----------------------------------------------------------------------------
# Run-level validation + order-safe assembly (mirrors prepare_combined_run()).
# -----------------------------------------------------------------------------
prepare_decoupled_run <- function(fit_params, priors_named, fixed_params,
                                  summary_stats, observed_named) {
  # --- parameters ---
  if (length(fit_params) == 0L) stop("`fit_params` must name >= 1 parameter.", call. = FALSE)
  if (anyDuplicated(fit_params)) {
    stop("`fit_params` has duplicates: ",
         paste(fit_params[duplicated(fit_params)], collapse = ", "), ".", call. = FALSE)
  }
  unknown_p <- setdiff(fit_params, DECOUPLED_PARAM_NAMES)
  if (length(unknown_p)) {
    stop("Unknown fit_params: ", paste(unknown_p, collapse = ", "),
         ". Available: ", paste(DECOUPLED_PARAM_NAMES, collapse = ", "), ".", call. = FALSE)
  }
  miss_prior <- setdiff(fit_params, names(priors_named))
  if (length(miss_prior)) {
    stop("Missing prior(s) for fitted parameter(s): ",
         paste(miss_prior, collapse = ", "), ".", call. = FALSE)
  }
  priors <- unname(priors_named[fit_params])

  user_fixed   <- if (is.null(fixed_params)) list() else fixed_params
  fixed_values <- utils::modifyList(DECOUPLED_PARAM_DEFAULTS, user_fixed)
  fixed_values <- fixed_values[DECOUPLED_PARAM_NAMES]

  # --- summaries ---
  if (length(summary_stats) == 0L) stop("`summary_stats` must name >= 1 metric.", call. = FALSE)
  if (anyDuplicated(summary_stats)) {
    stop("`summary_stats` has duplicates: ",
         paste(summary_stats[duplicated(summary_stats)], collapse = ", "), ".", call. = FALSE)
  }
  unknown_s <- setdiff(summary_stats, DECOUPLED_AVAILABLE_SUMMARIES)
  if (length(unknown_s)) {
    stop("Unknown summary_stats: ", paste(unknown_s, collapse = ", "),
         ". Available: ", paste(DECOUPLED_AVAILABLE_SUMMARIES, collapse = ", "), ".", call. = FALSE)
  }
  miss_obs <- setdiff(summary_stats, names(observed_named))
  if (length(miss_obs)) {
    stop("Missing observed target(s) for summary(ies): ",
         paste(miss_obs, collapse = ", "), ".", call. = FALSE)
  }
  observed_summaries <- observed_named[summary_stats]

  list(
    priors             = priors,
    observed_summaries = observed_summaries,
    fixed_values       = fixed_values,
    fit_params         = fit_params,
    summary_stats      = summary_stats
  )
}


# -----------------------------------------------------------------------------
# Observed d_p05_p95 target from a weekly death series.
# -----------------------------------------------------------------------------
# The model's d_p05_p95 is the day span between the 5th and 95th percentile of
# individual death DATES. The observed analogue from a weekly death curve is the
# day span between the points where cumulative deaths cross `lo` and `hi`,
# obtained by interpolating the cumulative-death curve (weekly resolution).
#   weekly_deaths : deaths per week, in chronological order.
#   week_day      : optional day index for each week (default = week midpoints).
observed_d_p05_p95 <- function(weekly_deaths, lo = 0.05, hi = 0.95,
                               bin_width = 7, week_day = NULL) {
  w <- as.numeric(weekly_deaths); w[!is.finite(w)] <- 0
  if (sum(w) <= 0) stop("`weekly_deaths` must sum to a positive number.", call. = FALSE)
  if (is.null(week_day)) week_day <- (seq_along(w) - 0.5) * bin_width
  cum   <- cumsum(w) / sum(w)
  dayq  <- function(q) approx(cum, week_day, q, rule = 2, ties = "ordered")$y
  unname(dayq(hi) - dayq(lo))
}


# -----------------------------------------------------------------------------
# Run-provenance metadata (self-describing .rds + source-able sidecar).
# -----------------------------------------------------------------------------
make_decoupled_run_metadata <- function(scenario_id = NULL,
                                        fit_params = NULL,
                                        priors = NULL,
                                        fixed_values = NULL,
                                        summary_stats = NULL,
                                        observed_summaries = NULL,
                                        hcw_base_prob = 0.25,
                                        general_hospital_quarantine_efficacy = NULL,
                                        safe_funeral_efficacy = NULL,
                                        peak_settings = NULL,
                                        extra = list()) {
  c(list(
    approach           = "decoupled",
    scenario_id        = scenario_id,
    fit_params         = fit_params,
    priors             = priors,
    fixed_values       = fixed_values,
    summary_stats      = summary_stats,
    observed_summaries = observed_summaries,
    hcw_base_prob      = hcw_base_prob,
    fixed_efficacies   = list(
      general_hospital_quarantine_efficacy = general_hospital_quarantine_efficacy,
      safe_funeral_efficacy                = safe_funeral_efficacy
    ),
    peak_settings      = peak_settings,
    # d_p05_p95 is the 5-95% span of death DATES; reconstruct observed via
    # observed_d_p05_p95() from a weekly death series (NOT the old "duration").
    duration_definition = "d_p05_p95 = day span between 5th and 95th percentile of death dates",
    created_at         = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")
  ), extra)
}

attach_decoupled_run_metadata <- function(result, metadata) {
  result$run_metadata <- metadata
  result
}

write_decoupled_run_metadata <- function(metadata, rds_path) {
  stopifnot(is.character(rds_path), length(rds_path) == 1L)
  stem     <- sub("\\.rds$", "", rds_path, ignore.case = TRUE)
  txt_path <- paste0(stem, "_metadata.txt")
  con <- file(txt_path, open = "wt"); on.exit(close(con), add = TRUE)
  writeLines(c(
    "# Auto-generated DECOUPLED-efficacy ABC run metadata.",
    "# source() this file to load the list `decoupled_run_metadata`.",
    "# NOTE: d_p05_p95 is the 5-95% span of death DATES -- recompute the observed",
    "#       target with observed_d_p05_p95() from a weekly death series.",
    ""
  ), con)
  cat("decoupled_run_metadata <- ", file = con)
  dput(metadata, file = con)
  paths <- txt_path
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    json_path <- paste0(stem, "_metadata.json")
    jsonlite::write_json(metadata, json_path, auto_unbox = TRUE, pretty = TRUE, null = "null")
    paths <- c(paths, json_path)
  }
  invisible(paths)
}


# -----------------------------------------------------------------------------
# Prior predictive check (optional; draw from priors, run the in-process model).
# -----------------------------------------------------------------------------
prior_predictive_check_decoupled <- function(n_draws,
                                             prior_list,
                                             fit_params,
                                             fixed_values,
                                             base,
                                             tv,
                                             invariants,
                                             general_hospital_quarantine_efficacy,
                                             safe_funeral_efficacy,
                                             hcw_base_prob = 0.25,
                                             parallel = FALSE,
                                             n_replicates = 20,
                                             seeding_cases = 25,
                                             takeoff_death_threshold = 100,
                                             summary_stats = DEFAULT_DECOUPLED_SUMMARIES,
                                             bin_width = DECOUPLED_PEAK_BIN_WIDTH,
                                             time_origin = DECOUPLED_PEAK_TIME_ORIGIN,
                                             takeoff_fraction_threshold = 1) {
  if (length(prior_list) != length(fit_params)) {
    stop("length(prior_list) must equal length(fit_params).", call. = FALSE)
  }
  sample_one <- function(spec, n) {
    params <- as.numeric(spec[-1])
    switch(spec[1],
           "unif"        = runif(n,  params[1], params[2]),
           "normal"      = rnorm(n,  params[1], params[2]),
           "lognormal"   = rlnorm(n, params[1], params[2]),
           "exponential" = rexp(n,   params[1]),
           stop("Unsupported prior distribution: ", spec[1]))
  }
  draws <- as.data.frame(lapply(prior_list, sample_one, n = n_draws))
  names(draws) <- fit_params
  theta_list <- lapply(seq_len(nrow(draws)), function(i) as.numeric(draws[i, ]))

  run_one <- function(theta) {
    fiber_abc_model_decoupled(
      theta = theta, fit_params = fit_params, fixed_values = fixed_values,
      base = base, tv = tv, invariants = invariants,
      general_hospital_quarantine_efficacy = general_hospital_quarantine_efficacy,
      safe_funeral_efficacy = safe_funeral_efficacy, hcw_base_prob = hcw_base_prob,
      n_replicates = n_replicates, seeding_cases = seeding_cases,
      takeoff_death_threshold = takeoff_death_threshold,
      summary_stats = summary_stats, bin_width = bin_width, time_origin = time_origin,
      takeoff_fraction_threshold = takeoff_fraction_threshold
    )
  }
  sims_list <- if (parallel && requireNamespace("future.apply", quietly = TRUE)) {
    future.apply::future_lapply(theta_list, run_one, future.seed = TRUE)
  } else {
    lapply(theta_list, run_one)
  }
  cbind(draws, as.data.frame(do.call(rbind, sims_list)))
}
