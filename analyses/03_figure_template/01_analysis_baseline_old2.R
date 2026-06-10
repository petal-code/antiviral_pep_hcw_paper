# =============================================================================
# 01_analysis_baseline.R
#
# Runs baseline (no OBV) simulations only.
# OBV efficacy and coverage effects are applied post-hoc in helper functions.
#
# Both scenarios use the DECOUPLED fit:
#   (R0, prop_funeral, etu_efficacy, ppe_efficacy, hcw_risk_scalar)
#
# Jobs are split by posterior particle -- each worker handles a slice of
# particles and runs all (scenario x rep) combinations for them.
#
# N_WORKERS must be a divisor of N_PARTICLES (e.g. 1,2,4,5,10,20,25,50,100).
# =============================================================================

library(here)
library(future)
library(future.apply)
library(fiber)

source(here("functions", "setup_model_parameters.R"))
source(here("functions", "calculate_model_approx_r0.R"))
source(here("functions", "abc_calibration_functions_common.R"))
source(here("functions", "abc_calibration_functions_decoupled.R"))
source(here("functions", "abc_posterior.R"))

env_decoupled <- new.env(parent = globalenv())
sys.source(here("functions", "abc_calibration_functions_decoupled.R"), envir = env_decoupled)

# =============================================================================
# Configuration
# =============================================================================

# >>> SET THIS: must be a divisor of N_PARTICLES (1,2,4,5,10,20,25,50,100)
N_WORKERS        <- 10L

N_PARTICLES      <- 200
N_REPS           <- 10
SEEDING_CASES    <- 25
CHECK_FINAL_SIZE <- 50000
RESAMPLE_SEED    <- 42L
SEED_BASE        <- 20260601L

HCW_BASE_PROB <- 0.25

stopifnot(N_PARTICLES %% N_WORKERS == 0L)
PARTICLES_PER_WORKER <- N_PARTICLES %/% N_WORKERS

SCENARIOS <- list(
  WestAfrica = list(
    id           = "Worst_WestAfrica",
    rds          = here("outputs", "02_ABC_model_fits_Final",
                        "fiber_ABC_SMC_Worst_WestAfrica_Decoupled_20260608_162044_check_NP5_NS4_NBREPS_30_NBSIMUL_472.RDS"),
    scenario_csv = here("data-processed", "final_six_scenario_values_original_approach.csv"),
    param_names  = c("R0", "prop_funeral", "etu_efficacy", "ppe_efficacy", "hcw_risk_scalar"),
    check_final_size = 40000
  ),
  DRC = list(
    id           = "Middle_DRC_ConflictSmoothed_PlusPlus",
    rds          = here("outputs", "02_ABC_model_fits_Final",
                        "fiber_ABC_SMC_Middle_DRC_ConflictSmoothed_PlusPlus_Decoupled_20260607_215621_check_NP5_NS4_NBREPS_30_NBSIMUL_590.RDS"),
    scenario_csv = here("data-processed", "final_six_scenario_values_original_approach.csv"),
    param_names  = c("R0", "prop_funeral", "etu_efficacy", "ppe_efficacy", "hcw_risk_scalar"),
    check_final_size = 10000
  )
)

OUT_DIR <- here("outputs", "simulation_baseline")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

check_model_function_version()

# =============================================================================
# Pre-compute particle args in main session
# =============================================================================
message("Pre-computing particle args...")

scenario_setups <- lapply(names(SCENARIOS), function(sc_name) {
  sc <- SCENARIOS[[sc_name]]
  message(sprintf("  %s (decoupled fit)...", sc_name))
  
  res       <- readRDS(sc$rds)
  posterior <- data.frame(weight = res$weights)
  for (j in seq_along(sc$param_names)) {
    posterior[[sc$param_names[j]]] <- res$param[, j]
  }
  theta <- downsample_posterior(posterior, n_sets = N_PARTICLES,
                                seed = RESAMPLE_SEED, param_names = sc$param_names)
  message(sprintf("    Posterior: %d particles -> resampled %d", nrow(posterior), N_PARTICLES))
  
  scenario_matrix <- read_scenario_matrix(sc$scenario_csv)
  mp <- make_model_parameters(
    scenario_id     = sc$id,
    scenario_matrix = scenario_matrix,
    overrides       = list(seeding_cases    = SEEDING_CASES,
                           check_final_size = sc$check_final_size)
  )
  
  inv <- compute_R0_invariants(args = mp$args, n = 50000, seed = 42L)
  
  build_decoupled <- env_decoupled$build_abc_model_args_decoupled
  
  particle_args <- lapply(seq_len(N_PARTICLES), function(p) {
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
    args$check_final_size <- CHECK_FINAL_SIZE
    args
  })
  
  list(
    sc_name       = sc_name,
    sc_idx        = match(sc_name, names(SCENARIOS)),
    theta         = theta,
    particle_args = particle_args
  )
})
names(scenario_setups) <- names(SCENARIOS)

# =============================================================================
# Build job list: one job = one particle slice
# =============================================================================
jobs <- lapply(seq_len(N_WORKERS), function(w) {
  p_from <- (w - 1L) * PARTICLES_PER_WORKER + 1L
  p_to   <-  w       * PARTICLES_PER_WORKER
  list(worker_id = w, p_from = p_from, p_to = p_to)
})

message(sprintf(
  "%d workers | %d particles each | %d scenarios x %d reps",
  N_WORKERS, PARTICLES_PER_WORKER, length(SCENARIOS), N_REPS
))

# =============================================================================
# Run in parallel — baseline only, no OBV args set
# =============================================================================
n_workers <- min(N_WORKERS, future::availableCores(), length(jobs))
plan(multisession, workers = n_workers)
t_start <- proc.time()

future_lapply(jobs, function(job) {
  for (p in seq(job$p_from, job$p_to)) {
    for (sc_name in names(scenario_setups)) {
      setup  <- scenario_setups[[sc_name]]
      sc_idx <- setup$sc_idx
      
      for (r in seq_len(N_REPS)) {
        fname    <- sprintf("%s_p%03d_r%d.rds", sc_name, p, r)
        out_path <- file.path(OUT_DIR, fname)
        if (file.exists(out_path)) next
        
        seed <- SEED_BASE +
          (sc_idx - 1L) * N_PARTICLES * N_REPS +
          (p      - 1L) * N_REPS +
          (r      - 1L)
        
        args      <- setup$particle_args[[p]]
        args$seed <- seed
        
        out   <- do.call(fiber::branching_process_main, args)
        tdf   <- out$tdf
        cases <- tdf[!is.na(tdf$time_infection_absolute), ]
        
        is_hcw <- cases$class == "HCW"
        died   <- !is.na(cases$outcome) & cases$outcome
        
        result <- list(
          scenario         = sc_name,
          particle_id      = p,
          rep              = r,
          R0               = setup$theta$R0[p],
          prop_funeral     = setup$theta$prop_funeral[p],
          etu_efficacy     = setup$theta$etu_efficacy[p],
          ppe_efficacy     = setup$theta$ppe_efficacy[p],
          hcw_risk_scalar  = setup$theta$hcw_risk_scalar[p],
          tdf              = cases,
          n_infections     = nrow(cases),
          n_hcw_infections = sum(is_hcw),
          n_deaths         = sum(died),
          n_hcw_deaths     = sum(died & is_hcw),
          duration         = max(cases$time_outcome_absolute, na.rm = TRUE)
        )
        
        saveRDS(result, out_path)
        message(sprintf("  [w%02d | p%03d | %s | r%d] -> %d inf, %d HCW deaths",
                        job$worker_id, p, sc_name, r,
                        result$n_infections, result$n_hcw_deaths))
      }
    }
  }
  
  message(sprintf("Worker %d done (particles %d-%d).",
                  job$worker_id, job$p_from, job$p_to))
  NULL
},
future.globals = list(
  scenario_setups = scenario_setups,
  OUT_DIR         = OUT_DIR,
  N_PARTICLES     = N_PARTICLES,
  N_REPS          = N_REPS,
  SEED_BASE       = SEED_BASE
),
future.packages = c("fiber", "here"),
future.seed     = TRUE)

plan(sequential)

elapsed <- proc.time() - t_start
n_files <- length(list.files(OUT_DIR, pattern = "\\.rds$"))
message(sprintf("Done in %.1f minutes. Total files: %d",
                elapsed["elapsed"] / 60, n_files))