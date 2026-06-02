# abc_calibration_functions_npi.R
# -----------------------------------------------------------------------------
# NPI-EFFICACY ABC approach. Fits (R0, prop_funeral, npi_scaler); the single
# scaler positions ppe_efficacy and etu_efficacy on explicit [min, max]
# intervals (see DEFAULT_NPI_SPEC / npi_efficacy_from_scaler). general-hospital
# and safe-funeral efficacies are FIXED; prob_hcw_cond_*_hospital fixed at 0.25
# (NO HCW-risk scalar). Because etu_efficacy is fitted, D and F are recomputed
# per particle from cached compute_R0_invariants().
#
# WHICH SUMMARIES ARE FITTED is configurable via abc_config$summary_stats (any
# subset of "takeoff" + AVAILABLE_ABC_METRICS, default = the historical four).
# The replicate loop + summary selection are shared with the in-process model via
# run_abc_particle() in the common file, so the same code fits the four-number
# summary OR an extended set (e.g. adding time_to_peak / peak_height) -- set
# summary_stats (and, optionally, peak_bin_width / peak_time_origin) in the run
# script's ABC_CONFIG and match observed_summaries' length/order.
#
# Approach-specific functions ONLY. The generic helpers live in
# abc_calibration_functions_common.R.
#
# Contents:
#   DEFAULT_NPI_SPEC                 : default [min,max] for the two fitted effs
#   npi_efficacy_from_scaler()       : scaler s in [-1,1] -> (ppe, etu) efficacies
#   npi_efficacy_from_scaler_maxscale(): alternative one-sided scaler
#   build_abc_model_args()           : scaler -> effs, per-particle D/F, model args
#   fiber_abc_model()                : in-process single-core ABC model
#   bootstrap_abc_worker()           : per-worker setup (base/tv/invariants/effs)
#   fiber_abc_model_parallel()       : the function ABC_sequential() calls
#   prior_predictive_check()         : draw from priors, run the non-parallel model
#   make_npi_run_metadata()          : serialisable record of npi_spec + settings
#   attach_npi_run_metadata()        : embed that record in the ABC result object
#   write_npi_run_metadata()         : write a .txt (source-able) + .json sidecar
#
# Requires setup_model_parameters.R, calculate_model_approx_r0.R, and
# abc_calibration_functions_common.R sourced first, and library(fiber).
# -----------------------------------------------------------------------------

DEFAULT_NPI_SPEC <- list(
  ppe_efficacy = list(min = 0.20, max = 0.90),   # PLACEHOLDER
  etu_efficacy = list(min = 0.50, max = 0.95)    # PLACEHOLDER
)


npi_efficacy_from_scaler <- function(s, npi_spec = DEFAULT_NPI_SPEC) {
  if (!is.numeric(s) || length(s) != 1L || is.na(s) || s < -1 || s > 1) {
    stop("`npi_scaler` (s) must be a single number in [-1, 1].", call. = FALSE)
  }
  out <- list()
  for (nm in names(npi_spec)) {
    lo <- npi_spec[[nm]]$min
    hi <- npi_spec[[nm]]$max
    if (is.null(lo) || is.null(hi) || hi < lo) {
      stop("npi_spec[['", nm, "']] must have min <= max.", call. = FALSE)
    }
    mid  <- (lo + hi) / 2
    half <- (hi - lo) / 2
    out[[nm]] <- mid + s * half
  }
  out
}


npi_efficacy_from_scaler_maxscale <- function(q, npi_spec = DEFAULT_NPI_SPEC) {
  if (!is.numeric(q) || length(q) != 1L || is.na(q) || q <= 0 || q > 1) {
    stop("`npi_scaler` (q) must be a single number in (0, 1].", call. = FALSE)
  }
  out <- list()
  for (nm in names(npi_spec)) out[[nm]] <- q * npi_spec[[nm]]$max
  out
}


build_abc_model_args <- function(R0,
                                 prop_funeral,
                                 npi_scaler,
                                 base,
                                 tv,
                                 invariants,
                                 npi_spec = DEFAULT_NPI_SPEC,
                                 general_hospital_quarantine_efficacy,
                                 safe_funeral_efficacy,
                                 seeding_cases = 25) {
  effs <- npi_efficacy_from_scaler(npi_scaler, npi_spec)

  D <- D_from_invariants(
    invariants,
    etu_efficacy = effs$etu_efficacy,
    general_hospital_quarantine_efficacy = general_hospital_quarantine_efficacy
  )
  F_fun <- F_from_invariants(invariants, safe_funeral_efficacy)
  means <- solve_offspring_means(R0, prop_funeral, D, F_fun)

  args <- c(base, tv)
  args$mn_offspring_genPop  <- means$mn_genPop
  args$mn_offspring_funeral <- means$mn_funeral

  # Fitted + fixed efficacies for this particle.
  args$ppe_efficacy <- effs$ppe_efficacy
  args$etu_efficacy <- effs$etu_efficacy
  args$general_hospital_quarantine_efficacy <- general_hospital_quarantine_efficacy
  args$safe_funeral_efficacy <- safe_funeral_efficacy

  args$seed          <- NULL
  args$seeding_cases <- seeding_cases
  args
}


fiber_abc_model <- function(theta,
                            base,
                            tv,
                            invariants,
                            npi_spec = DEFAULT_NPI_SPEC,
                            general_hospital_quarantine_efficacy,
                            safe_funeral_efficacy,
                            n_replicates = 30,
                            seeding_cases = 25,
                            takeoff_death_threshold = 100,
                            summary_stats = DEFAULT_SUMMARY_STATS,
                            bin_width     = DEFAULT_PEAK_BIN_WIDTH,
                            time_origin   = DEFAULT_PEAK_TIME_ORIGIN) {

  args <- build_abc_model_args(
    R0 = theta[1], prop_funeral = theta[2], npi_scaler = theta[3],
    base = base, tv = tv, invariants = invariants, npi_spec = npi_spec,
    general_hospital_quarantine_efficacy = general_hospital_quarantine_efficacy,
    safe_funeral_efficacy = safe_funeral_efficacy,
    seeding_cases = seeding_cases
  )

  # Aggregation over replicates + selection of the fitted summaries is shared
  # with the parallel path via run_abc_particle() (functions/.../common.R).
  run_abc_particle(
    args,
    n_reps                  = n_replicates,
    takeoff_death_threshold = takeoff_death_threshold,
    summary_stats           = summary_stats,
    bin_width               = bin_width,
    time_origin             = time_origin
  )
}


bootstrap_abc_worker <- function(setup_path,
                                 functions_path,
                                 r0_path,
                                 scenario_csv,
                                 scenario_id,
                                 abc_config = list(),
                                 model_overrides = list()) {
  # Source helpers FIRST so DEFAULT_SCALAR_INPUTS / DEFAULT_NPI_SPEC exist
  # before default_config (below) references them. This matters on a fresh
  # self-bootstrapping PSOCK worker, whose globalenv initially holds only the
  # functions from cfg$functions_path.
  source(setup_path)
  source(functions_path)
  source(r0_path)
  ## Generic ABC helpers (abc_summarise, etc.) live alongside the approach file.
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
    # Which summaries to fit (default = historical four) + weekly-curve binning
    # for the optional time_to_peak / peak_height summaries.
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


fiber_abc_model_parallel <- function(theta_with_seed) {
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

    bootstrap_abc_worker(
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

  run_abc_particle(
    args,
    n_reps                  = cfg_run$n_reps,
    takeoff_death_threshold = cfg_run$takeoff_death_threshold,
    summary_stats           = cfg_run$summary_stats,
    bin_width               = cfg_run$peak_bin_width,
    time_origin             = cfg_run$peak_time_origin
  )
}


prior_predictive_check <- function(n_draws,
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
                                   summary_stats = c(DEFAULT_SUMMARY_STATS, "n_cases"),
                                   bin_width = DEFAULT_PEAK_BIN_WIDTH,
                                   time_origin = DEFAULT_PEAK_TIME_ORIGIN) {

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
    fiber_abc_model(
      theta = theta,
      base = base, tv = tv, invariants = invariants, npi_spec = npi_spec,
      general_hospital_quarantine_efficacy = general_hospital_quarantine_efficacy,
      safe_funeral_efficacy = safe_funeral_efficacy,
      n_replicates = n_replicates,
      seeding_cases = seeding_cases,
      takeoff_death_threshold = takeoff_death_threshold,
      summary_stats = summary_stats,
      bin_width = bin_width,
      time_origin = time_origin
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


# -----------------------------------------------------------------------------
# RUN PROVENANCE: store the npi_spec (and other settings) WITH the fit output.
# -----------------------------------------------------------------------------
# A fitted npi_scaler is meaningless without the [min, max] interval it indexes
# (npi_efficacy_from_scaler() needs npi_spec to recover ppe/etu efficacies). That
# spec lived only in the (mutable) run script and an ephemeral save_abc_config()
# tempfile, so a finished .rds could not be interpreted on its own. These three
# helpers fix that: build a serialisable record, attach it to the ABC result so
# the .rds is self-describing, and write a human-readable + source()-able sidecar
# next to it.

# Build a structured, serialisable record of everything needed to interpret a
# fitted npi_scaler downstream -- chiefly npi_spec, plus the fixed efficacies,
# priors, observed targets and scenario id for completeness. `extra` lets callers
# bolt on anything else (e.g. result_filename, git sha).
make_npi_run_metadata <- function(npi_spec,
                                  scenario_id = NULL,
                                  priors = NULL,
                                  observed_summaries = NULL,
                                  general_hospital_quarantine_efficacy = NULL,
                                  safe_funeral_efficacy = NULL,
                                  param_names = c("R0", "prop_funeral", "npi_scaler"),
                                  extra = list()) {
  c(list(
    npi_spec           = npi_spec,
    scenario_id        = scenario_id,
    param_names        = param_names,
    priors             = priors,
    observed_summaries = observed_summaries,
    fixed_efficacies   = list(
      general_hospital_quarantine_efficacy = general_hospital_quarantine_efficacy,
      safe_funeral_efficacy                = safe_funeral_efficacy
    ),
    created_at         = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")
  ), extra)
}

# Attach the metadata to an ABC_sequential() result so the .rds is self-describing.
# Stores a top-level $npi_spec (the cheap read path downstream code uses) plus the
# full $run_metadata record. Adding named elements does not disturb $param /
# $weights / $stats consumers.
attach_npi_run_metadata <- function(result, metadata) {
  result$npi_spec     <- metadata$npi_spec
  result$run_metadata <- metadata
  result
}

# Write a sidecar next to a result .rds. Always writes "<stem>_metadata.txt",
# which is both human-readable AND source()-able (it assigns `npi_run_metadata`),
# so a fit whose .rds predates metadata embedding can still be documented and
# read back. Also writes a .json copy when jsonlite is available. Returns the
# path(s) written.
write_npi_run_metadata <- function(metadata, rds_path) {
  stopifnot(is.character(rds_path), length(rds_path) == 1L)
  stem     <- sub("\\.rds$", "", rds_path, ignore.case = TRUE)
  txt_path <- paste0(stem, "_metadata.txt")

  con <- file(txt_path, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(c(
    "# Auto-generated NPI-efficacy ABC run metadata.",
    "# source() this file to load the list `npi_run_metadata`.",
    "# npi_spec maps npi_scaler in [-1, 1] -> efficacy:",
    "#   efficacy = midpoint + npi_scaler * (max - min) / 2  (per fitted efficacy).",
    "# Required by any downstream analysis that reconstructs efficacies from the fit.",
    ""
  ), con)
  cat("npi_run_metadata <- ", file = con)
  dput(metadata, file = con)

  paths <- txt_path
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    json_path <- paste0(stem, "_metadata.json")
    jsonlite::write_json(metadata, json_path, auto_unbox = TRUE,
                         pretty = TRUE, null = "null")
    paths <- c(paths, json_path)
  }
  invisible(paths)
}


