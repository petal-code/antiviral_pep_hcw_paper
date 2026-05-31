# abc_calibration_functions_hcwRisk.R
# -----------------------------------------------------------------------------
# HCW-RISK ABC approach. Fits (R0, prop_funeral, hcw_risk_scalar); the scalar
# scales both prob_hcw_cond_*_hospital from a symmetric base. Conditional
# efficacies (etu/ppe/general/safe_funeral) are FIXED inputs from make_base_args();
# D and F are computed once (efficacies do not vary across particles).
#
# Approach-specific functions ONLY. The generic helpers live in
# abc_calibration_functions_common.R.
#
# Contents:
#   map_abc_params_to_model() : (R0, prop_funeral, hcw_risk_scalar) -> 4 model pars
#   build_abc_model_args()    : splice mapped pars onto base + tv args
#   fiber_abc_model()         : in-process single-core ABC model
#   bootstrap_abc_worker()    : per-worker setup (base/tv/D/F into globalenv)
#   fiber_abc_model_parallel(): the function ABC_sequential() calls on workers
#   prior_predictive_check()  : draw from priors, run the non-parallel model
#
# Requires setup_model_parameters.R, calculate_model_approx_r0.R, and
# abc_calibration_functions_common.R sourced first, and library(fiber).
# -----------------------------------------------------------------------------

map_abc_params_to_model <- function(R0, prop_funeral, hcw_risk_scalar,
                                    D, F_fun, hcw_base_prob = 0.25) {
  hcw_hospital <- pmin(hcw_base_prob * hcw_risk_scalar, 1.0)
  list(
    mn_offspring_genPop           = (1 - prop_funeral) * R0 / D,
    mn_offspring_funeral          =      prop_funeral  * R0 / F_fun,
    prob_hcw_cond_genPop_hospital = hcw_hospital,
    prob_hcw_cond_hcw_hospital    = hcw_hospital
  )
}


build_abc_model_args <- function(R0,
                                 prop_funeral,
                                 hcw_risk_scalar,
                                 base,
                                 tv,
                                 D,
                                 F_fun,
                                 seeding_cases = 25,
                                 hcw_base_prob = 0.25) {
  model_pars <- map_abc_params_to_model(
    R0 = R0, prop_funeral = prop_funeral, hcw_risk_scalar = hcw_risk_scalar,
    D = D, F_fun = F_fun, hcw_base_prob = hcw_base_prob
  )

  args <- c(base, tv)
  args$mn_offspring_genPop           <- model_pars$mn_offspring_genPop
  args$mn_offspring_funeral          <- model_pars$mn_offspring_funeral
  args$prob_hcw_cond_genPop_hospital <- model_pars$prob_hcw_cond_genPop_hospital
  args$prob_hcw_cond_hcw_hospital    <- model_pars$prob_hcw_cond_hcw_hospital
  args$seed          <- NULL
  args$seeding_cases <- seeding_cases
  args
}


fiber_abc_model <- function(theta,
                            base,
                            tv,
                            D,
                            F_fun,
                            n_replicates = 30,
                            seeding_cases = 25,
                            takeoff_death_threshold = 100,
                            hcw_base_prob = 0.25,
                            include_n_cases = FALSE) {

  R0              <- theta[1]
  prop_funeral    <- theta[2]
  hcw_risk_scalar <- theta[3]

  args <- build_abc_model_args(
    R0 = R0, prop_funeral = prop_funeral, hcw_risk_scalar = hcw_risk_scalar,
    base = base, tv = tv, D = D, F_fun = F_fun,
    seeding_cases = seeding_cases, hcw_base_prob = hcw_base_prob
  )

  reps <- vapply(
    seq_len(n_replicates),
    function(i) run_one_abc_replicate(args),
    numeric(4)
  )
  rownames(reps) <- c("n_cases", "n_deaths", "n_hcw_deaths", "duration")

  took_off     <- reps["n_deaths", ] >= takeoff_death_threshold
  takeoff_frac <- mean(took_off)

  if (!any(took_off)) {
    out <- c(takeoff = 0, n_deaths = 0, n_hcw_deaths = 0, duration = 0)
    if (include_n_cases) out <- c(out, n_cases = 0)
    return(out)
  }

  out <- c(takeoff      = takeoff_frac,
           n_deaths     = mean(reps["n_deaths",     took_off]),
           n_hcw_deaths = mean(reps["n_hcw_deaths", took_off]),
           duration     = mean(reps["duration",     took_off]))

  if (include_n_cases) {
    out <- c(out, n_cases = mean(reps["n_cases", took_off]))
  }
  out
}


bootstrap_abc_worker <- function(setup_path,
                                 functions_path,
                                 r0_path,
                                 scenario_csv,
                                 scenario_id,
                                 abc_config = list(),
                                 model_overrides = list()) {
  default_config <- list(
    check_final_size        = 30000,
    takeoff_death_threshold = 100,
    n_reps                  = 30,
    seeding_cases           = 25,
    hcw_base_prob           = 0.25,
    setup_R0_n              = 100000L,
    setup_R0_seed           = 42L,
    setup_funeral_share     = 0.5
  )
  abc_config <- utils::modifyList(default_config, abc_config)

  source(setup_path)
  source(functions_path)
  source(r0_path)
  ## Generic ABC helpers (abc_summarise, etc.) live alongside the approach file.
  source(file.path(dirname(functions_path), "abc_calibration_functions_common.R"))

  if (!"fiber" %in% loadedNamespaces()) {
    library(fiber)
  }

  scenario_matrix <- read_scenario_matrix(scenario_csv)

  # check_final_size lives in DEFAULT_SCALAR_INPUTS, so any caller override
  # for it should be merged on top of the abc_config default.
  scalar_overrides <- model_overrides
  if (!"check_final_size" %in% names(scalar_overrides)) {
    scalar_overrides$check_final_size <- abc_config$check_final_size
  }

  mp <- make_model_parameters(
    scenario_id = scenario_id,
    scenario_matrix = scenario_matrix,
    overrides = scalar_overrides
  )

  setup_solve <- solve_offspring_means_for_R0(
    R0   = 1.0,
    args = mp$args,
    proportion_transmission_from_funerals = abc_config$setup_funeral_share,
    n    = abc_config$setup_R0_n,
    seed = abc_config$setup_R0_seed
  )

  assign("base_args",           mp$base_args,                       envir = globalenv())
  assign("tv_args_model",       mp$tv_args,                         envir = globalenv())
  assign("D_direct_multiplier", setup_solve$D_direct_multiplier,    envir = globalenv())
  assign("F_funeral_multiplier", setup_solve$F_funeral_multiplier,  envir = globalenv())
  assign(".abc_config",         abc_config,                         envir = globalenv())
  assign(".fiber_abc_ready",    TRUE,                               envir = globalenv())

  invisible(list(
    setup_solve = setup_solve,
    abc_config = abc_config,
    scenario_label = mp$scenario_label
  ))
}


fiber_abc_model_parallel <- function(theta_with_seed) {
  if (!isTRUE(get0(".fiber_abc_ready", envir = globalenv()))) {
    # The worker's globalenv only contains this function (exported by
    # ABC_sequential's PSOCK setup). Custom helpers like bootstrap_abc_worker()
    # are NOT yet defined here, so the self-bootstrap path can only use
    # base R until we have sourced the helper files.
    cfg_path <- Sys.getenv("FIBER_ABC_CONFIG")
    if (cfg_path == "" || !file.exists(cfg_path)) {
      stop("FIBER_ABC_CONFIG env var not set or file missing. ",
           "Call save_abc_config(<config>) in the main process before ",
           "running ABC_sequential().", call. = FALSE)
    }
    cfg <- readRDS(cfg_path)

    # We only need bootstrap_abc_worker right now; that function will
    # source setup / r0 itself before doing anything that depends on them.
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
    hcw_risk_scalar = theta_with_seed[4],
    base            = get("base_args",            envir = globalenv()),
    tv              = get("tv_args_model",        envir = globalenv()),
    D               = get("D_direct_multiplier",  envir = globalenv()),
    F_fun           = get("F_funeral_multiplier", envir = globalenv()),
    seeding_cases   = cfg_run$seeding_cases,
    hcw_base_prob   = cfg_run$hcw_base_prob
  )

  reps <- vapply(seq_len(cfg_run$n_reps), function(i) {
    out <- do.call(branching_process_main, args)
    abc_summarise(out)
  }, numeric(4))
  rownames(reps) <- c("n_cases", "n_deaths", "n_hcw_deaths", "duration")

  took_off <- reps["n_deaths", ] >= cfg_run$takeoff_death_threshold
  if (!any(took_off)) {
    return(c(takeoff = 0, n_deaths = 0, n_hcw_deaths = 0, duration = 0))
  }
  c(takeoff      = mean(took_off),
    n_deaths     = mean(reps["n_deaths",     took_off]),
    n_hcw_deaths = mean(reps["n_hcw_deaths", took_off]),
    duration     = mean(reps["duration",     took_off]))
}


prior_predictive_check <- function(n_draws,
                                   prior_list,
                                   base,
                                   tv,
                                   D,
                                   F_fun,
                                   param_names = c("R0", "prop_funeral", "hcw_risk_scalar"),
                                   parallel = FALSE,
                                   n_replicates = 30,
                                   seeding_cases = 25,
                                   takeoff_death_threshold = 100,
                                   hcw_base_prob = 0.25,
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
    fiber_abc_model(
      theta = theta,
      base = base, tv = tv, D = D, F_fun = F_fun,
      n_replicates = n_replicates,
      seeding_cases = seeding_cases,
      takeoff_death_threshold = takeoff_death_threshold,
      hcw_base_prob = hcw_base_prob,
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
        fn <- function(theta) {
          out <- run_one_theta(theta)
          p()
          out
        }
        future.apply::future_lapply(theta_list, fn, future.seed = TRUE)
      })
    } else {
      future.apply::future_lapply(theta_list, run_one_theta, future.seed = TRUE)
    }
  } else {
    if (requireNamespace("progressr", quietly = TRUE)) {
      progressr::with_progress({
        p <- progressr::progressor(steps = length(theta_list))
        lapply(theta_list, function(theta) {
          out <- run_one_theta(theta)
          p()
          out
        })
      })
    } else {
      lapply(theta_list, run_one_theta)
    }
  }

  sims <- do.call(rbind, sims_list)
  cbind(draws, as.data.frame(sims))
}


