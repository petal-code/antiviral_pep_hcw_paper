# abc_calibration_functions_new_approach.R
# -----------------------------------------------------------------------------
# ABC-SMC calibration helpers for the *revamped* NPI parameterisation.
#
# Fitted parameters (3):
#   R0           : baseline reproduction number for a genPop seeding case at
#                  t = 0 (under this particle's efficacies).
#   prop_funeral : share of R0 attributable to funeral transmission at t = 0.
#   npi_scaler   : a SINGLE shared scaler that moves the two fitted conditional
#                  efficacies (ppe_efficacy, etu_efficacy) together. See
#                  npi_efficacy_from_scaler().
#
# general_hospital_quarantine_efficacy and safe_funeral_efficacy are held FIXED
# (passed through the worker config). prob_hcw_cond_*_hospital are fixed at the
# honest 0.25 base — there is no HCW-risk scalar in this approach.
#
# Because etu_efficacy is fitted, the direct multiplier D depends on the
# particle. We therefore cache the efficacy-INDEPENDENT R0 invariants once per
# scenario (compute_R0_invariants()) and recompute the cheap closed-form D and
# F per particle (D_from_invariants() / F_from_invariants()).
#
# Requires setup_model_parameters_new_approach.R and
# calculate_model_approx_r0_new_approach.R sourced first, and library(fiber).
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# NPI efficacy scaler
# -----------------------------------------------------------------------------
# >>> PLACEHOLDER bounds — UPDATE WITH REAL NUMBERS. <<<
# Approach B (symmetric, two-sided around a central value). Each fitted efficacy
# lives on an explicit [min, max]; a single scaler s in [-1, 1] positions BOTH
# of them on their respective intervals:
#
#     eff_i(s) = mid_i + s * half_range_i,
#       mid_i = (min_i + max_i) / 2,  half_range_i = (max_i - min_i) / 2
#
#   s = -1 -> min_i,   s = 0 -> central (midpoint),   s = +1 -> max_i.
#
# The RELATIVE widths (max_i - min_i) set how strongly the scaler leans on PPE
# vs ETU (and hence the HCW-skew-vs-burden mix that makes the scaler
# identifiable) — choose them deliberately, not just the endpoints.

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

# --- Approach A (alternative; one-sided downward scaling from a best-case max).
# To swap: point build_abc_model_args() at this function instead, and change the
# npi_scaler prior to c("unif", q_lo, 1.0). Kept here for easy comparison.
#   eff_i(q) = q * max_i,   q in (0, 1].
npi_efficacy_from_scaler_maxscale <- function(q, npi_spec = DEFAULT_NPI_SPEC) {
  if (!is.numeric(q) || length(q) != 1L || is.na(q) || q <= 0 || q > 1) {
    stop("`npi_scaler` (q) must be a single number in (0, 1].", call. = FALSE)
  }
  out <- list()
  for (nm in names(npi_spec)) out[[nm]] <- q * npi_spec[[nm]]$max
  out
}


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
# Maps (R0, prop_funeral, npi_scaler) onto fiber inputs:
#   1. npi_scaler -> (ppe_efficacy, etu_efficacy) via npi_efficacy_from_scaler()
#   2. D = D_from_invariants(inv, etu_efficacy, general_hosp)   [per particle]
#      F = F_from_invariants(inv, safe_funeral_efficacy)        [fixed sfe]
#   3. mn_offspring_* = solve_offspring_means(R0, prop_funeral, D, F)
# The fixed efficacies (general_hosp, safe_funeral) are passed in.

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

run_one_abc_replicate <- function(args) {
  out <- do.call(branching_process_main, args)
  abc_summarise(out)
}


# -----------------------------------------------------------------------------
# In-process ABC model (single-core)
# -----------------------------------------------------------------------------

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
# Sources the new-approach helpers, loads fiber, builds base_args / tv_args,
# computes the efficacy-INDEPENDENT R0 invariants ONCE, and stashes everything
# (plus npi_spec and the fixed efficacies) in globalenv.

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
    safe_funeral_efficacy   = DEFAULT_SCALAR_INPUTS$safe_funeral_efficacy
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


# -----------------------------------------------------------------------------
# Parallel ABC model
# -----------------------------------------------------------------------------
# Called by EasyABC::ABC_sequential() with use_seed = TRUE. theta_with_seed is
# c(seed, R0, prop_funeral, npi_scaler). Returns the 4-element summary vector.

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
                                   invariants,
                                   npi_spec = DEFAULT_NPI_SPEC,
                                   general_hospital_quarantine_efficacy,
                                   safe_funeral_efficacy,
                                   param_names = c("R0", "prop_funeral", "npi_scaler"),
                                   parallel = FALSE,
                                   n_replicates = 30,
                                   seeding_cases = 25,
                                   takeoff_death_threshold = 100,
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
      base = base, tv = tv, invariants = invariants, npi_spec = npi_spec,
      general_hospital_quarantine_efficacy = general_hospital_quarantine_efficacy,
      safe_funeral_efficacy = safe_funeral_efficacy,
      n_replicates = n_replicates,
      seeding_cases = seeding_cases,
      takeoff_death_threshold = takeoff_death_threshold,
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


# -----------------------------------------------------------------------------
# Inspecting ABC progress on disk
# -----------------------------------------------------------------------------

abc_progress <- function(dir = getwd(),
                         tolerance_target = 1.0,
                         param_names = c("R0", "prop_funeral", "npi_scaler")) {

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
                              param_names = c("R0", "prop_funeral", "npi_scaler")) {

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
                                   param_names = c("R0", "prop_funeral", "npi_scaler"),
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
