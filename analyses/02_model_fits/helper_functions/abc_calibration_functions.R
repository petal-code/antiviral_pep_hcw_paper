# abc_calibration_functions.R
# -----------------------------------------------------------------------------
# Helper functions for ABC-SMC calibration of fiber's branching-process model.
#
# Contents (grouped by purpose):
#
#   --- Summarising a single model output ---
#     abc_summarise()                  : single-simulation summary used as the
#                                        per-replicate input to ABC.
#
#   --- Turning ABC parameters into model args ---
#     build_abc_model_args()           : maps (R0, prop_funeral,
#                                        hcw_risk_scalar) onto fiber inputs
#                                        using pre-computed D and F
#                                        multipliers from
#                                        solve_offspring_means_for_R0().
#     run_one_abc_replicate()          : do.call(branching_process_main) +
#                                        abc_summarise() for one replicate.
#
#   --- ABC model wrappers ---
#     fiber_abc_model()                : in-process model (single core).
#     save_abc_config()                : serialises the everything-workers-
#                                        need-to-know config to disk and
#                                        advertises it via the
#                                        FIBER_ABC_CONFIG env var. PSOCK
#                                        workers spawned by ABC_sequential
#                                        inherit the env var and read it
#                                        back in fiber_abc_model_parallel.
#     bootstrap_abc_worker()           : per-worker setup that builds
#                                        base_args, tv_args_model, D, F and
#                                        stashes them in globalenv so that
#                                        fiber_abc_model_parallel() can be
#                                        called by EasyABC::ABC_sequential().
#     fiber_abc_model_parallel()       : the function ABC_sequential calls
#                                        on each worker. On first call it
#                                        self-bootstraps from the config
#                                        file pointed at by FIBER_ABC_CONFIG;
#                                        thereafter it reads precomputed
#                                        state from globalenv.
#
#   --- Prior-predictive checking ---
#     prior_predictive_check()         : draws from priors and runs the
#                                        non-parallel model to inspect what
#                                        the priors imply about the
#                                        summaries.
#
#   --- Output-directory helpers ---
#     make_abc_output_dir()            : returns a freshly-created, per-run
#                                        subdirectory so ABC_sequential's
#                                        output_step* / tolerance_step* /
#                                        n_simul_tot_step* files don't
#                                        pollute the repo root.
#     with_abc_output_dir()            : runs an expression with the working
#                                        directory temporarily set to that
#                                        subdirectory; restores cwd even on
#                                        error.
#
#   --- Inspecting / reconstructing in-progress or completed ABC runs ---
#     abc_progress()                   : prints progress from output_step*
#                                        / tolerance_step* / n_simul_tot_step*
#                                        files written by ABC_sequential.
#     abc_compare_steps()              : weighted summary across all
#                                        completed steps.
#     reconstruct_abc_result()         : rebuilds an ABC_sequential()-style
#                                        result object from disk.
#
# Requires setup_model_parameters.R and calculate_model_approx_r0.R to be
# sourced first, and fiber to be loaded (library(fiber)).
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Per-replicate summary
# -----------------------------------------------------------------------------

abc_summarise <- function(out) {
  tdf <- out$tdf
  tdf <- tdf[!is.na(tdf$time_infection_absolute), , drop = FALSE]

  if (nrow(tdf) == 0L) {
    return(c(n_cases = 0, n_deaths = 0, n_hcw_deaths = 0, duration = 0))
  }

  deaths <- !is.na(tdf$outcome) & tdf$outcome
  hcw    <- !is.na(tdf$class) & tdf$class == "HCW"

  n_cases      <- nrow(tdf)
  n_deaths     <- sum(deaths)
  n_hcw_deaths <- sum(deaths & hcw)

  if (n_deaths == 0L) {
    duration <- 0
  } else {
    # Duration is measured from the first death to the last event, so it tracks
    # the part of the outbreak that would be observable in surveillance data.
    t_first_death <- min(tdf$time_outcome_absolute[deaths], na.rm = TRUE)
    t_last_event  <- max(tdf$time_outcome_absolute,         na.rm = TRUE)
    duration      <- max(t_last_event - t_first_death, 0)
  }

  c(n_cases      = n_cases,
    n_deaths     = n_deaths,
    n_hcw_deaths = n_hcw_deaths,
    duration     = duration)
}


# -----------------------------------------------------------------------------
# Parameter -> model args
# -----------------------------------------------------------------------------
# Single-type R0 inversion:
#   R0_direct  = mn_offspring_genPop  * D
#   R0_funeral = mn_offspring_funeral * F
# so for a target R0 split (1 - prop_funeral) / prop_funeral:
#   mn_offspring_genPop  = (1 - prop_funeral) * R0 / D
#   mn_offspring_funeral =      prop_funeral  * R0 / F
# D and F come from solve_offspring_means_for_R0() and depend only on the
# fixed scenario inputs at t = 0, so they can be computed once per scenario.
#
# hcw_risk_scalar scales BOTH HCW-given-hospital probabilities by the same
# multiplier (capped at 1), starting from a symmetric base `hcw_base_prob`:
#   prob_hcw_cond_genPop_hospital <- pmin(hcw_base_prob * hcw_risk_scalar, 1)
#   prob_hcw_cond_hcw_hospital    <- pmin(hcw_base_prob * hcw_risk_scalar, 1)
# A symmetric base is a more honest reflection of prior uncertainty about
# the relative magnitudes than the older asymmetric defaults (0.12, 0.20).

build_abc_model_args <- function(R0,
                                 prop_funeral,
                                 hcw_risk_scalar,
                                 base,
                                 tv,
                                 D,
                                 F_fun,
                                 seeding_cases = 25,
                                 hcw_base_prob = 0.25) {
  mn_genPop  <- (1 - prop_funeral) * R0 / D
  mn_funeral <-      prop_funeral  * R0 / F_fun

  args <- c(base, tv)
  args$mn_offspring_genPop           <- mn_genPop
  args$mn_offspring_funeral          <- mn_funeral
  args$prob_hcw_cond_genPop_hospital <- pmin(hcw_base_prob * hcw_risk_scalar, 1.0)
  args$prob_hcw_cond_hcw_hospital    <- pmin(hcw_base_prob * hcw_risk_scalar, 1.0)
  args$seed          <- NULL
  args$seeding_cases <- seeding_cases
  args
}

run_one_abc_replicate <- function(args) {
  out <- do.call(branching_process_main, args)
  abc_summarise(out)
}


# -----------------------------------------------------------------------------
# In-process ABC model (single-core)
# -----------------------------------------------------------------------------
# A replicate "took off" if it produced >= takeoff_death_threshold deaths.
# Reports P(takeoff) plus the mean of each conditional summary among
# replicates that took off; if no replicate takes off, returns zeros so the
# distance to observed is large.

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


# -----------------------------------------------------------------------------
# Worker config
# -----------------------------------------------------------------------------
# ABC_sequential() spawns its own PSOCK cluster (via the n_cluster argument)
# and calls the model function with just a theta_with_seed vector. The per-
# worker context (paths, ABC tuning, model overrides) therefore has to be
# discoverable by the worker itself.
#
# save_abc_config() writes the config to an RDS file in the main process and
# sets the FIBER_ABC_CONFIG environment variable to the file path. PSOCK
# workers inherit env vars from the parent R session; fiber_abc_model_parallel
# inlines the readRDS() on first call so the bootstrap can run before any
# of these helper functions are sourced on the worker.

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


# -----------------------------------------------------------------------------
# Output-directory helpers
# -----------------------------------------------------------------------------
# ABC_sequential writes its intermediate files (output_step*, tolerance_step*,
# n_simul_tot_step*) to the current working directory. To stop these polluting
# the repo, make a per-run subdirectory and chdir into it just for the ABC
# call. abc_progress() / abc_compare_steps() / reconstruct_abc_result() all
# already accept a `dir` argument, so they can read back from the same place.
#
# make_abc_output_dir() returns the path of a freshly-created subdirectory of
# the form <base_dir>/<subdir>/<scenario_id>[_YYYYMMDD_HHMMSS][_<label>], and
# disambiguates with a numeric suffix if the path somehow already exists.

make_abc_output_dir <- function(base_dir,
                                scenario_id,
                                label = NULL,
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

  out <- file.path(base_dir, subdir, parts)
  suffix <- 0L
  candidate <- out
  while (dir.exists(candidate)) {
    suffix <- suffix + 1L
    candidate <- paste0(out, "_", sprintf("%02d", suffix))
  }
  out <- candidate

  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  out
}

# with_abc_output_dir() runs an arbitrary expression with the working
# directory temporarily set to `output_dir`. `expr` is evaluated lazily, so
# anything written to cwd by ABC_sequential lands in `output_dir`; the
# original cwd is restored even if the expression errors.

with_abc_output_dir <- function(output_dir, expr) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  old <- setwd(output_dir)
  on.exit(setwd(old), add = TRUE)
  force(expr)
}


# -----------------------------------------------------------------------------
# Worker bootstrap
# -----------------------------------------------------------------------------
# Reads paths + ABC tuning + model overrides, sources the helper files,
# loads the fiber package, precomputes base_args, tv_args_model, D, F, and
# stashes everything in globalenv. The .fiber_abc_ready sentinel prevents
# re-running on subsequent calls.

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


# -----------------------------------------------------------------------------
# Parallel ABC model
# -----------------------------------------------------------------------------
# Called by EasyABC::ABC_sequential() with use_seed = TRUE. theta_with_seed
# is c(seed, R0, prop_funeral, hcw_risk_scalar). Returns the 4-element
# summary vector matched against observed_summaries.
#
# On its FIRST call inside a given worker, the function reads the config
# advertised by FIBER_ABC_CONFIG, sources the helper files, and calls
# bootstrap_abc_worker() to precompute base_args / tv_args_model / D / F.
# All subsequent calls reuse the values from globalenv.

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


# -----------------------------------------------------------------------------
# Prior predictive check
# -----------------------------------------------------------------------------

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


# -----------------------------------------------------------------------------
# Inspecting ABC progress on disk
# -----------------------------------------------------------------------------
# EasyABC::ABC_sequential writes per-step files to the working directory:
#   output_step<k>          : the particle population (weight + params + stats)
#   tolerance_step<k>       : the tolerance for step k (k >= 2)
#   n_simul_tot_step<k>     : the cumulative number of simulations after step k
# The functions below read those files so that progress can be inspected
# while a run is ongoing.

abc_progress <- function(dir = getwd(),
                         tolerance_target = 1.0,
                         param_names = c("R0", "prop_funeral", "hcw_risk_scalar")) {

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
  colnames(last_out) <- c("weight", param_names,
                          "takeoff", "n_deaths", "n_hcw_deaths", "duration")
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
                              param_names = c("R0", "prop_funeral", "hcw_risk_scalar")) {

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
    colnames(df) <- c("weight", param_names,
                      "takeoff", "n_deaths", "n_hcw_deaths", "duration")
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
      mean_takeoff   = round(sum(w * df$takeoff), 3),
      mean_deaths    = round(sum(w * df$n_deaths)),
      mean_hcw       = round(sum(w * df$n_hcw_deaths)),
      mean_duration  = round(sum(w * df$duration)),
      stringsAsFactors = FALSE
    )

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
