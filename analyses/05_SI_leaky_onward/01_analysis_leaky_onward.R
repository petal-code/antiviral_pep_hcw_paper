# =============================================================================
# 01_analysis_leaky_onward.R
#
# Reviewer-response sensitivity analysis: "leaky obeldesivir".
#
# The main model assumes obeldesivir (OBV), when effective, fully prevents a
# health-worker (HCW) infection -- the treated person leaves the transmission
# chain entirely. A reviewer notes this is optimistic (a real antiviral might
# only REDUCE onward transmission). This script quantifies what that does to
# OBV's estimated DEATH impact -- the paper's main outcome, focused on HCW
# deaths -- for the DRC and West Africa archetypes, using the same posterior
# draws and the same 80%-coverage / 80%-efficacy programme as the main figures.
#
# METHOD (see analyses/05_SI_leaky_onward/leaky_onward_helpers.R and
#         fiber::estimate_leaky_onward for the full rationale)
#   For each posterior draw (200, weighted) x stochastic replicate (10) we run
#   the WITH-OBV base model once (Figure-1 arm: cov80_obv80), take its
#   `prevented_completed` set of averted infections, and forward-simulate them:
#     * A1  no-OBV counterfactual (OBV off, full transmission) -> total deaths
#           OBV averts, index + downstream (the number the current index-only
#           reporting leaves uncounted).
#     * A2  leaky-OBV across residual transmissibility r = 0..1 -> as OBV is made
#           progressively leakier, how many deaths still occur, keeping OBV's
#           protection against DEATH on everyone it treats (index and downstream).
#   Deaths that occur under leak = untreated would-be deaths only; deaths averted
#   at leakiness r = A1 - (leaked deaths). We median the 10 replicates within
#   each draw, then aggregate across draws (done in 02_plot_leaky_onward.R).
#
# This is heavy (200 x 10 x 2 base runs, plus 12 forward sims each; the forward
# sims near r = 1 grow large -- capped at MAX_TREE_SIZE with a warning). Set
# QUICK_TEST <- TRUE for a tiny end-to-end smoke run first.
#
# Output: a compact per-run metrics table at
#   analyses/05_SI_leaky_onward/_intermediate/leaky_onward_per_run.rds
# (gitignored, regenerable), read by 02_plot_leaky_onward.R.
# =============================================================================

library(here)
library(future)
library(future.apply)

# ---- Install the exact fiber branch this analysis requires ------------------
# estimate_leaky_onward() and out$sim_info$params live on this branch only.
INSTALL_FIBER <- TRUE
FIBER_REF     <- "claude/antiviral-efficacy-obeldesivir-e4tewv"
if (INSTALL_FIBER) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  remotes::install_github("petal-code/fiber", ref = FIBER_REF, upgrade = "never")
}
library(fiber)
if (!"estimate_leaky_onward" %in% getNamespaceExports("fiber")) {
  stop("Installed fiber lacks estimate_leaky_onward(); set INSTALL_FIBER <- TRUE ",
       "and install ref '", FIBER_REF, "'.", call. = FALSE)
}

source(here("functions", "setup_model_parameters.R"))
source(here("functions", "calculate_model_approx_r0.R"))
source(here("functions", "abc_calibration_functions_common.R"))
source(here("functions", "abc_calibration_functions_decoupled.R"))
source(here("functions", "abc_posterior.R"))
source(here("analyses", "05_SI_leaky_onward", "leaky_onward_helpers.R"))

env_decoupled <- new.env(parent = globalenv())
sys.source(here("functions", "abc_calibration_functions_decoupled.R"), envir = env_decoupled)

# =============================================================================
# Configuration
# =============================================================================
QUICK_TEST <- FALSE                       # TRUE -> tiny smoke run (see below)

N_WORKERS     <- 100L
N_PARTICLES   <- 200L
N_REPS        <- 10L
SEEDING_CASES <- 25L
RESAMPLE_SEED <- 42L                       # SAME draws as the main figures
SEED_BASE     <- 20260701L                 # SAME base-run seeds as Figure 1
HCW_BASE_PROB <- 0.25

TAKEOFF_DEATH_THRESHOLD <- 100L            # per-rep take-off, matches the fit
MAX_RETRIES             <- 50L

# Leaky sweep + forward-sim controls
R_GRID        <- seq(0, 1, by = 0.1)       # residual transmissibility
LEAKY_SEED    <- 77000L                    # base seed for the forward sims

# OBV programme = the Figure-1 arm: 80% coverage, 80% efficacy, dpc 0, HCWs in
# every setting. (Coverage is the reviewer-requested 80%.)
OBV <- list(
  coverage         = function(t) rep(0.8, length(t)),
  adherence        = 1.0,
  dpc              = 0,
  efficacy         = 0.80,
  target_class     = "HCW",
  target_locations = c("hospital", "community", "funeral")
)

if (QUICK_TEST) {
  N_WORKERS   <- 2L
  N_PARTICLES <- 4L
  N_REPS      <- 2L
  R_GRID      <- c(0, 0.5, 1)
}

stopifnot(N_PARTICLES %% N_WORKERS == 0L)
PARTICLES_PER_WORKER <- N_PARTICLES %/% N_WORKERS

INT_DIR <- here("analyses", "05_SI_leaky_onward", "_intermediate")
dir.create(INT_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Scenarios (identical fits / sizes to 01_analysis_figure1.R). MAX_TREE_SIZE is
# set to the base run's check_final_size, so each forward sim is bounded by the
# same cap as the base epidemic (warns rather than truncating silently).
# =============================================================================
SCENARIOS <- list(
  DRC = list(
    id           = "Middle_DRC_ConflictSmoothed_PlusPlus",
    rds          = here("outputs", "02_ABC_model_fits_Final",
                        "fiber_ABC_SMC_Middle_DRC_ConflictSmoothed_PlusPlus_Decoupled_20260607_215621_check_NP5_NS4_NBREPS_30_NBSIMUL_590.RDS"),
    scenario_csv = here("data-processed", "final_six_scenario_values_original_approach.csv"),
    param_names  = c("R0", "prop_funeral", "etu_efficacy", "ppe_efficacy", "hcw_risk_scalar"),
    check_final_size = 10000
  ),
  WestAfrica = list(
    id           = "Worst_WestAfrica",
    rds          = here("outputs", "02_ABC_model_fits_Final",
                        "fiber_ABC_SMC_Worst_WestAfrica_Decoupled_20260608_162044_check_NP5_NS4_NBREPS_30_NBSIMUL_472.RDS"),
    scenario_csv = here("data-processed", "final_six_scenario_values_original_approach.csv"),
    param_names  = c("R0", "prop_funeral", "etu_efficacy", "ppe_efficacy", "hcw_risk_scalar"),
    check_final_size = 40000
  )
)

# Robustly resolve the fit RDS path (filenames may differ in case on Linux).
resolve_rds <- function(path) {
  if (file.exists(path)) return(path)
  hit <- list.files(dirname(path), pattern = paste0("^", basename(tools::file_path_sans_ext(path)), "\\.rds$"),
                    ignore.case = TRUE, full.names = TRUE)
  if (length(hit) >= 1L) return(hit[1])
  stop("Fit RDS not found near: ", path, call. = FALSE)
}

tryCatch(check_model_function_version(), error = function(e) NULL)

# =============================================================================
# Pre-compute per-particle base fiber args per scenario (mirrors Figure 1)
# =============================================================================
message("Pre-computing particle args...")
scenario_setups <- lapply(names(SCENARIOS), function(sc_name) {
  sc  <- SCENARIOS[[sc_name]]
  res <- readRDS(resolve_rds(sc$rds))

  posterior <- data.frame(weight = res$weights)
  for (j in seq_along(sc$param_names)) posterior[[sc$param_names[j]]] <- res$param[, j]
  theta <- downsample_posterior(posterior, n_sets = N_PARTICLES,
                                seed = RESAMPLE_SEED, param_names = sc$param_names)

  scenario_matrix <- read_scenario_matrix(sc$scenario_csv)
  mp <- make_model_parameters(
    scenario_id     = sc$id,
    scenario_matrix = scenario_matrix,
    overrides       = list(seeding_cases = SEEDING_CASES, check_final_size = sc$check_final_size)
  )
  inv <- compute_R0_invariants(args = mp$args, n = 50000, seed = 42L)

  build_decoupled <- env_decoupled$build_abc_model_args_decoupled
  base_particle_args <- lapply(seq_len(N_PARTICLES), function(p) {
    args <- build_decoupled(
      R0              = theta$R0[p],
      prop_funeral    = theta$prop_funeral[p],
      etu_efficacy    = theta$etu_efficacy[p],
      ppe_efficacy    = theta$ppe_efficacy[p],
      hcw_risk_scalar = theta$hcw_risk_scalar[p],
      base            = mp$base_args,
      tv              = mp$tv_args,
      invariants      = inv,
      general_hospital_quarantine_efficacy = DEFAULT_SCALAR_INPUTS$general_hospital_quarantine_efficacy,
      safe_funeral_efficacy                = DEFAULT_SCALAR_INPUTS$safe_funeral_efficacy,
      hcw_base_prob   = HCW_BASE_PROB,
      seeding_cases   = SEEDING_CASES
    )
    args$check_final_size <- sc$check_final_size
    args
  })

  list(sc_name = sc_name, sc_idx = match(sc_name, names(SCENARIOS)),
       base_particle_args = base_particle_args, check_final_size = sc$check_final_size)
})
names(scenario_setups) <- names(SCENARIOS)

# =============================================================================
# Worker: for one particle, run all reps x scenarios, computing the leaky-onward
# metrics for each. Returns a long metrics data.frame. Resumable: one compact rds
# per (scenario, particle); skip if present.
# =============================================================================
run_particle <- function(job) {
  out_rows <- list()
  for (p in seq(job$p_from, job$p_to)) {
    for (sc_name in names(scenario_setups)) {
      setup   <- scenario_setups[[sc_name]]
      sc_idx  <- setup$sc_idx
      cfs     <- setup$check_final_size
      pfile   <- file.path(INT_DIR, sprintf("%s_p%03d.rds", sc_name, p))
      if (file.exists(pfile)) { out_rows[[length(out_rows) + 1L]] <- readRDS(pfile); next }

      per_p <- list()
      for (rep in seq_len(N_REPS)) {
        # Base-run seed: same construction as Figure 1 so draws/reps are matched.
        base_seed <- SEED_BASE +
          (sc_idx - 1L) * N_PARTICLES * N_REPS +
          (p      - 1L) * N_REPS +
          (rep    - 1L)

        args <- setup$base_particle_args[[p]]
        args$obv_pep_enabled          <- TRUE
        args$obv_pep_coverage         <- OBV$coverage
        args$obv_pep_adherence        <- OBV$adherence
        args$obv_pep_dpc              <- OBV$dpc
        args$obv_pep_efficacy         <- OBV$efficacy
        args$obv_pep_target_class     <- OBV$target_class
        args$obv_pep_target_locations <- OBV$target_locations

        # Take-off retry (matches the fit): resample seed until the outbreak
        # establishes (>= TAKEOFF_DEATH_THRESHOLD deaths), max MAX_RETRIES.
        out <- NULL; retry <- 0L; cur <- base_seed
        repeat {
          args$seed <- cur
          out <- do.call(branching_process_main, args)
          tdf <- out$tdf[!is.na(out$tdf$time_infection_absolute), ]
          if (sum(!is.na(tdf$outcome) & tdf$outcome) >= TAKEOFF_DEATH_THRESHOLD) break
          retry <- retry + 1L; cur <- cur + 1L
          if (retry >= MAX_RETRIES) break
        }

        # Leaky-onward metrics for this base run. Forward-sim seed derived from
        # the (matched) base seed so the whole analysis is reproducible.
        m <- compute_leaky_onward_metrics(
          out           = out,
          seed          = LEAKY_SEED + base_seed,
          r_grid        = R_GRID,
          max_tree_size = cfs
        )
        m$scenario <- sc_name; m$particle <- p; m$rep <- rep
        per_p[[length(per_p) + 1L]] <- m
      }
      per_p_df <- do.call(rbind, per_p)
      saveRDS(per_p_df, pfile)
      out_rows[[length(out_rows) + 1L]] <- per_p_df
    }
  }
  do.call(rbind, out_rows)
}

jobs <- lapply(seq_len(N_WORKERS), function(w) {
  list(worker_id = w,
       p_from = (w - 1L) * PARTICLES_PER_WORKER + 1L,
       p_to   =  w       * PARTICLES_PER_WORKER)
})

message(sprintf("%d workers | %d particles each | %d scenarios x %d reps x %d r-values",
                N_WORKERS, PARTICLES_PER_WORKER, length(SCENARIOS), N_REPS, length(R_GRID)))

n_workers <- min(N_WORKERS, future::availableCores(), length(jobs))
plan(multisession, workers = n_workers)
t0 <- proc.time()
results <- future_lapply(
  jobs, run_particle,
  future.globals = list(
    scenario_setups = scenario_setups, INT_DIR = INT_DIR, N_PARTICLES = N_PARTICLES,
    N_REPS = N_REPS, SEED_BASE = SEED_BASE, LEAKY_SEED = LEAKY_SEED, R_GRID = R_GRID,
    OBV = OBV, TAKEOFF_DEATH_THRESHOLD = TAKEOFF_DEATH_THRESHOLD, MAX_RETRIES = MAX_RETRIES,
    compute_leaky_onward_metrics = compute_leaky_onward_metrics,
    count_tree_deaths = count_tree_deaths, count_tdf_deaths = count_tdf_deaths
  ),
  future.packages = c("fiber", "here"), future.seed = TRUE
)
plan(sequential)
message(sprintf("Done in %.1f min.", (proc.time() - t0)[3] / 60))

per_run <- do.call(rbind, results)
saveRDS(per_run, file.path(INT_DIR, "leaky_onward_per_run.rds"))
message("Saved per-run metrics: ", file.path(INT_DIR, "leaky_onward_per_run.rds"))
