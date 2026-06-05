# abc_calibration_functions_combined.R
# -----------------------------------------------------------------------------
# COMBINED ABC approach: ONE parameterisation that unions the NPI-efficacy and
# HCW-risk schemes and lets the run script choose, per run, WHICH parameters are
# fitted and WHICH summary statistics are matched. Nothing here is hard-coded to
# a fixed parameter triple or summary set.
#
# Canonical fittable parameters (COMBINED_PARAM_NAMES), in this order:
#   R0              : baseline reproduction number for a typical genPop seeding
#                     case at t = 0 (under this particle's efficacies).
#   prop_funeral    : share of R0 attributable to funeral transmission at t = 0.
#   npi_scaler      : s in [-1, 1] positioning ppe_efficacy + etu_efficacy on
#                     their NPI_SPEC [min, max] intervals (s = 0 -> midpoints).
#                     Mapped by npi_efficacy_from_scaler() (NPI-efficacy scheme).
#   hcw_risk_scalar : multiplies the symmetric prob_hcw_cond_*_hospital base
#                     (hcw_base_prob, default 0.25): prob = min(base * scalar, 1).
#                     (HCW-risk scheme; mirrors map_abc_params_to_model() in
#                     abc_calibration_functions_hcwRisk.R.)
#
# Any SUBSET of these four may be fitted; the rest are held at fixed values
# (COMBINED_PARAM_DEFAULTS, overridable per run via the run script's
# FIXED_PARAMS). Two useful reductions:
#   * fix hcw_risk_scalar = 1 (-> prob_hcw at the 0.25 base) and fit
#     (R0, prop_funeral, npi_scaler) == the pure NPI-efficacy scheme.
#   * fix npi_scaler at the value whose NPI_SPEC midpoints give your chosen
#     efficacies and fit (R0, prop_funeral, hcw_risk_scalar) == the HCW-risk
#     scheme (parameterised through NPI_SPEC rather than make_base_args()).
#
# Because etu_efficacy may vary per particle (when npi_scaler is fitted), D and F
# are recomputed per particle from cached compute_R0_invariants() -- exactly as
# in the NPI-efficacy scheme; when npi_scaler is fixed they simply come out
# constant. The replicate loop + summary selection are the SHARED
# run_abc_particle() from the common file, so WHICH summaries are fitted is the
# usual abc_config$summary_stats knob (any subset of "takeoff" +
# AVAILABLE_ABC_METRICS).
#
# Reuses (does NOT redefine):
#   * npi_efficacy_from_scaler(), DEFAULT_NPI_SPEC and the run-metadata helpers
#     make/attach/write_npi_run_metadata() from abc_calibration_functions_npi.R
#   * D_from_invariants(), F_from_invariants(), solve_offspring_means(),
#     compute_R0_invariants() from calculate_model_approx_r0.R
#   * run_abc_particle(), abc_summarise(), AVAILABLE_ABC_METRICS, the disk-
#     inspection trio and make_abc_output_dir()/save_abc_config() from the common
#     file
#   * make_model_parameters(), read_scenario_matrix(), DEFAULT_SCALAR_INPUTS from
#     setup_model_parameters.R
#
# Approach-specific functions ONLY. Required source order (run script AND worker
# bootstrap) -- npi BEFORE this file so its helpers exist:
#   setup_model_parameters.R -> abc_calibration_functions_common.R ->
#   abc_calibration_functions_npi.R -> THIS FILE -> calculate_model_approx_r0.R
# and library(fiber). bootstrap_abc_worker_combined() derives the npi/common
# paths from dirname(functions_path), the same way the other workers do.
# -----------------------------------------------------------------------------


# Canonical fittable parameters and the order in which the assembled parameter
# list is laid out (build_abc_model_args_combined() consumes these by name).
COMBINED_PARAM_NAMES <- c("R0", "prop_funeral", "npi_scaler", "hcw_risk_scalar")

# Values used for any canonical parameter NOT being fitted (override per run via
# the run script's FIXED_PARAMS). hcw_risk_scalar = 1 keeps prob_hcw_cond_*_hospital
# at the symmetric base (0.25); npi_scaler = 0 gives the central (midpoint)
# efficacies of NPI_SPEC.
COMBINED_PARAM_DEFAULTS <- list(
  R0              = 1.45,
  prop_funeral    = 0.35,
  npi_scaler      = 0,
  hcw_risk_scalar = 1.0
)


# Assemble the full named (R0, prop_funeral, npi_scaler, hcw_risk_scalar) list
# from the fitted subset `theta_fitted` (values in `fit_params` order) spliced
# onto the complete `fixed_values` list (which must hold all COMBINED_PARAM_NAMES;
# prepare_combined_run() guarantees this). The fixed_values entries for fitted
# parameters are placeholders -- they are overwritten here.
assemble_combined_theta <- function(theta_fitted, fit_params, fixed_values) {
  if (length(theta_fitted) != length(fit_params)) {
    stop("length(theta_fitted) (", length(theta_fitted), ") must equal ",
         "length(fit_params) (", length(fit_params), ").", call. = FALSE)
  }
  full <- fixed_values[COMBINED_PARAM_NAMES]
  names(full) <- COMBINED_PARAM_NAMES
  for (i in seq_along(fit_params)) full[[fit_params[i]]] <- theta_fitted[[i]]
  full
}


# (R0, prop_funeral, npi_scaler, hcw_risk_scalar) -> a ready-to-run fiber args
# list. Unions the two schemes: the NPI mapping sets the efficacies and (via the
# per-particle D / F) the offspring means; the HCW-risk mapping sets the hospital
# HCW-exposure probabilities.
build_abc_model_args_combined <- function(R0,
                                          prop_funeral,
                                          npi_scaler,
                                          hcw_risk_scalar,
                                          base,
                                          tv,
                                          invariants,
                                          npi_spec = DEFAULT_NPI_SPEC,
                                          general_hospital_quarantine_efficacy,
                                          safe_funeral_efficacy,
                                          hcw_base_prob = 0.25,
                                          seeding_cases = 25) {
  # --- NPI-efficacy mapping: scaler -> (ppe, etu), then per-particle D, F, means.
  effs <- npi_efficacy_from_scaler(npi_scaler, npi_spec)
  D <- D_from_invariants(
    invariants,
    etu_efficacy = effs$etu_efficacy,
    general_hospital_quarantine_efficacy = general_hospital_quarantine_efficacy
  )
  F_fun <- F_from_invariants(invariants, safe_funeral_efficacy)
  means <- solve_offspring_means(R0, prop_funeral, D, F_fun)

  # --- HCW-risk mapping: scale the symmetric hospital HCW-exposure base. Mirrors
  # map_abc_params_to_model() in abc_calibration_functions_hcwRisk.R.
  hcw_hospital <- pmin(hcw_base_prob * hcw_risk_scalar, 1.0)

  args <- c(base, tv)
  args$mn_offspring_genPop  <- means$mn_genPop
  args$mn_offspring_funeral <- means$mn_funeral

  # Fitted + fixed efficacies for this particle.
  args$ppe_efficacy <- effs$ppe_efficacy
  args$etu_efficacy <- effs$etu_efficacy
  args$general_hospital_quarantine_efficacy <- general_hospital_quarantine_efficacy
  args$safe_funeral_efficacy <- safe_funeral_efficacy

  # HCW per-contact exposure in hospital settings (both classes share the base).
  args$prob_hcw_cond_genPop_hospital <- hcw_hospital
  args$prob_hcw_cond_hcw_hospital    <- hcw_hospital

  args$seed          <- NULL
  args$seeding_cases <- seeding_cases
  args
}


# In-process single-core ABC model (used by prior_predictive_check_combined() and
# for ad-hoc testing). `theta` holds the fitted parameters in `fit_params` order.
fiber_abc_model_combined <- function(theta,
                                     fit_params,
                                     fixed_values,
                                     base,
                                     tv,
                                     invariants,
                                     npi_spec = DEFAULT_NPI_SPEC,
                                     general_hospital_quarantine_efficacy,
                                     safe_funeral_efficacy,
                                     hcw_base_prob = 0.25,
                                     n_replicates = 30,
                                     seeding_cases = 25,
                                     takeoff_death_threshold = 100,
                                     summary_stats = DEFAULT_SUMMARY_STATS,
                                     bin_width     = DEFAULT_PEAK_BIN_WIDTH,
                                     time_origin   = DEFAULT_PEAK_TIME_ORIGIN) {

  full <- assemble_combined_theta(theta, fit_params, fixed_values)
  args <- build_abc_model_args_combined(
    R0 = full$R0, prop_funeral = full$prop_funeral,
    npi_scaler = full$npi_scaler, hcw_risk_scalar = full$hcw_risk_scalar,
    base = base, tv = tv, invariants = invariants, npi_spec = npi_spec,
    general_hospital_quarantine_efficacy = general_hospital_quarantine_efficacy,
    safe_funeral_efficacy = safe_funeral_efficacy,
    hcw_base_prob = hcw_base_prob, seeding_cases = seeding_cases
  )

  # Aggregation over replicates + selection of the fitted summaries is the SHARED
  # run_abc_particle() (functions/abc_calibration_functions_common.R).
  run_abc_particle(
    args,
    n_reps                  = n_replicates,
    takeoff_death_threshold = takeoff_death_threshold,
    summary_stats           = summary_stats,
    bin_width               = bin_width,
    time_origin             = time_origin
  )
}


# Per-worker setup for the parallel path. Mirrors bootstrap_abc_worker() in the
# NPI file, but also stashes the combined-scheme choices (fit_params,
# fixed_values, hcw_base_prob) in globalenv for fiber_abc_model_parallel_combined().
# npi/common are sourced via dirname(functions_path) -- no extra config keys.
bootstrap_abc_worker_combined <- function(setup_path,
                                          functions_path,
                                          r0_path,
                                          scenario_csv,
                                          scenario_id,
                                          abc_config = list(),
                                          model_overrides = list()) {
  fn_dir <- dirname(functions_path)
  # Source helpers FIRST (npi before this file so npi_efficacy_from_scaler /
  # DEFAULT_NPI_SPEC exist before default_config references them, and common so
  # AVAILABLE_ABC_METRICS / run_abc_particle() exist).
  source(setup_path)
  source(file.path(fn_dir, "abc_calibration_functions_npi.R"))
  source(functions_path)
  source(r0_path)
  source(file.path(fn_dir, "abc_calibration_functions_common.R"))

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
    hcw_base_prob           = 0.25,
    general_hospital_quarantine_efficacy =
      DEFAULT_SCALAR_INPUTS$general_hospital_quarantine_efficacy,
    safe_funeral_efficacy   = DEFAULT_SCALAR_INPUTS$safe_funeral_efficacy,
    fit_params              = COMBINED_PARAM_NAMES,
    fixed_values            = COMBINED_PARAM_DEFAULTS,
    summary_stats           = DEFAULT_SUMMARY_STATS,
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
  assign(".hcw_base_prob", abc_config$hcw_base_prob,              envir = globalenv())
  assign(".fit_params",    abc_config$fit_params,                 envir = globalenv())
  assign(".fixed_values",  abc_config$fixed_values,               envir = globalenv())
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


# The function ABC_sequential() calls on each worker. Self-bootstraps on a cold
# PSOCK worker from the FIBER_ABC_CONFIG tempfile, then maps the fitted-parameter
# vector (theta_with_seed[-1], in .fit_params order) onto model args.
fiber_abc_model_parallel_combined <- function(theta_with_seed) {
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

    bootstrap_abc_worker_combined(
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

  full <- assemble_combined_theta(
    theta_fitted, cfg_run$fit_params, cfg_run$fixed_values
  )

  args <- build_abc_model_args_combined(
    R0              = full$R0,
    prop_funeral    = full$prop_funeral,
    npi_scaler      = full$npi_scaler,
    hcw_risk_scalar = full$hcw_risk_scalar,
    base            = get("base_args",     envir = globalenv()),
    tv              = get("tv_args_model", envir = globalenv()),
    invariants      = get("R0_invariants", envir = globalenv()),
    npi_spec        = get(".npi_spec",     envir = globalenv()),
    general_hospital_quarantine_efficacy = get(".fixed_general_hosp", envir = globalenv()),
    safe_funeral_efficacy                = get(".fixed_safe_funeral", envir = globalenv()),
    hcw_base_prob   = get(".hcw_base_prob", envir = globalenv()),
    seeding_cases   = cfg_run$seeding_cases
  )

  run_abc_particle(
    args,
    n_reps                  = cfg_run$n_reps,
    takeoff_death_threshold = cfg_run$takeoff_death_threshold,
    summary_stats           = cfg_run$summary_stats,
    bin_width               = cfg_run$peak_bin_width,
    time_origin             = cfg_run$peak_time_origin
  )
}


# -----------------------------------------------------------------------------
# Run-level validation + order-safe assembly.
# -----------------------------------------------------------------------------
# The run script declares its choices as NAMED lists keyed by parameter / summary
# name, so it never has to keep two positional lists in sync. This function
# validates them against the canonical sets and returns everything ABC needs, in
# the right order:
#   fit_params    : character vector, subset of COMBINED_PARAM_NAMES (fit order).
#   priors_named  : named list, one EasyABC prior spec per fitted parameter.
#   fixed_params  : named list of fixed values for the NON-fitted parameters
#                   (entries for fitted params are ignored; any non-fitted param
#                   absent here falls back to COMBINED_PARAM_DEFAULTS).
#   summary_stats : character vector, subset of c("takeoff", AVAILABLE_ABC_METRICS).
#   observed_named: named numeric, one observed target per fitted summary.
# Returns list(priors, observed_summaries, fixed_values, fit_params, summary_stats).
prepare_combined_run <- function(fit_params, priors_named, fixed_params,
                                 summary_stats, observed_named) {
  # --- parameters ---
  if (length(fit_params) == 0L) {
    stop("`fit_params` must name at least one parameter.", call. = FALSE)
  }
  if (anyDuplicated(fit_params)) {
    stop("`fit_params` contains duplicates: ",
         paste(fit_params[duplicated(fit_params)], collapse = ", "), ".", call. = FALSE)
  }
  unknown_p <- setdiff(fit_params, COMBINED_PARAM_NAMES)
  if (length(unknown_p)) {
    stop("Unknown fit_params: ", paste(unknown_p, collapse = ", "),
         ". Available: ", paste(COMBINED_PARAM_NAMES, collapse = ", "), ".",
         call. = FALSE)
  }
  miss_prior <- setdiff(fit_params, names(priors_named))
  if (length(miss_prior)) {
    stop("Missing prior(s) for fitted parameter(s): ",
         paste(miss_prior, collapse = ", "), ".", call. = FALSE)
  }
  priors <- unname(priors_named[fit_params])

  # Complete fixed-value list: canonical defaults overlaid with the user's
  # fixed_params. Fitted entries are placeholders (overwritten per particle).
  user_fixed   <- if (is.null(fixed_params)) list() else fixed_params
  fixed_values <- utils::modifyList(COMBINED_PARAM_DEFAULTS, user_fixed)
  fixed_values <- fixed_values[COMBINED_PARAM_NAMES]   # canonical names + order

  # --- summaries ---
  if (length(summary_stats) == 0L) {
    stop("`summary_stats` must name at least one metric.", call. = FALSE)
  }
  if (anyDuplicated(summary_stats)) {
    stop("`summary_stats` contains duplicates: ",
         paste(summary_stats[duplicated(summary_stats)], collapse = ", "), ".",
         call. = FALSE)
  }
  unknown_s <- setdiff(summary_stats, c("takeoff", AVAILABLE_ABC_METRICS))
  if (length(unknown_s)) {
    stop("Unknown summary_stats: ", paste(unknown_s, collapse = ", "),
         ". Available: takeoff, ", paste(AVAILABLE_ABC_METRICS, collapse = ", "),
         ".", call. = FALSE)
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
# Prior predictive check (optional; draw from priors, run the non-parallel model).
# -----------------------------------------------------------------------------
# Returns a data frame: one row per draw, the sampled fitted parameters followed
# by the simulated summaries. Handy for sizing priors and -- with summary_stats
# including peak_height -- for sanity-checking the observed peak target.
prior_predictive_check_combined <- function(n_draws,
                                            prior_list,
                                            fit_params,
                                            fixed_values,
                                            base,
                                            tv,
                                            invariants,
                                            npi_spec = DEFAULT_NPI_SPEC,
                                            general_hospital_quarantine_efficacy,
                                            safe_funeral_efficacy,
                                            hcw_base_prob = 0.25,
                                            parallel = FALSE,
                                            n_replicates = 30,
                                            seeding_cases = 25,
                                            takeoff_death_threshold = 100,
                                            summary_stats = c(DEFAULT_SUMMARY_STATS, "n_cases"),
                                            bin_width = DEFAULT_PEAK_BIN_WIDTH,
                                            time_origin = DEFAULT_PEAK_TIME_ORIGIN) {

  if (length(prior_list) != length(fit_params)) {
    stop("length(prior_list) (", length(prior_list), ") must equal ",
         "length(fit_params) (", length(fit_params), ").", call. = FALSE)
  }

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
  names(draws) <- fit_params

  theta_list <- lapply(seq_len(nrow(draws)), function(i) as.numeric(draws[i, ]))

  run_one_theta <- function(theta) {
    fiber_abc_model_combined(
      theta = theta, fit_params = fit_params, fixed_values = fixed_values,
      base = base, tv = tv, invariants = invariants, npi_spec = npi_spec,
      general_hospital_quarantine_efficacy = general_hospital_quarantine_efficacy,
      safe_funeral_efficacy = safe_funeral_efficacy, hcw_base_prob = hcw_base_prob,
      n_replicates = n_replicates, seeding_cases = seeding_cases,
      takeoff_death_threshold = takeoff_death_threshold,
      summary_stats = summary_stats, bin_width = bin_width, time_origin = time_origin
    )
  }

  sims_list <- if (parallel) {
    if (!requireNamespace("future.apply", quietly = TRUE)) {
      stop("Package 'future.apply' is required when parallel = TRUE.", call. = FALSE)
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
