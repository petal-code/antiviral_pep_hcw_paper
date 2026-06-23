# =============================================================================
# 01_analysis_figure4_baseline.R
#
# Runs baseline (no antiviral) simulations for Figure 4.
# West Africa archetype only.
#
# Policy A / B / stockpile / efficacy sweeps are all applied post-hoc
# in the extract script, so this script runs a single baseline arm.
#
# Output: one RDS per particle x rep under
#   outputs/simulation/figure4_baseline/WestAfrica/
#
# N_WORKERS must be a divisor of N_PARTICLES.
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
N_WORKERS     <- 100L
N_PARTICLES   <- 200L
N_REPS        <- 10L
SEEDING_CASES <- 25L
RESAMPLE_SEED <- 42L
SEED_BASE     <- 20260704L
HCW_BASE_PROB <- 0.25

TAKEOFF_DEATH_THRESHOLD <- 100L
MAX_RETRIES             <- 50L

stopifnot(N_PARTICLES %% N_WORKERS == 0L)
PARTICLES_PER_WORKER <- N_PARTICLES %/% N_WORKERS

# =============================================================================
# 1. Scenario: West Africa only
# =============================================================================
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

OUT_BASE <- here("outputs", "simulation", "figure4_baseline")
for (sc_name in names(SCENARIOS)) {
  dir.create(file.path(OUT_BASE, sc_name), recursive = TRUE, showWarnings = FALSE)
}

check_model_function_version()

# =============================================================================
# 2. Pre-compute particle args
# =============================================================================
message("Pre-computing particle args...")

scenario_setups <- lapply(names(SCENARIOS), function(sc_name) {
  sc <- SCENARIOS[[sc_name]]
  message(sprintf("  %s...", sc_name))
  
  res       <- readRDS(sc$rds)
  posterior <- data.frame(weight = res$weights)
  for (j in seq_along(sc$param_names)) {
    posterior[[sc$param_names[j]]] <- res$param[, j]
  }
  theta <- downsample_posterior(posterior, n_sets = N_PARTICLES,
                                seed = RESAMPLE_SEED, param_names = sc$param_names)
  
  scenario_matrix <- read_scenario_matrix(sc$scenario_csv)
  mp <- make_model_parameters(
    scenario_id     = sc$id,
    scenario_matrix = scenario_matrix,
    overrides       = list(seeding_cases    = SEEDING_CASES,
                           check_final_size = sc$check_final_size)
  )
  inv <- compute_R0_invariants(args = mp$args, n = 50000, seed = 42L)
  
  build_decoupled <- env_decoupled$build_abc_model_args_decoupled
  
  # Build baseline args with antiviral disabled
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
    args$check_final_size  <- sc$check_final_size
    args$obv_pep_enabled   <- FALSE   # no antiviral; policy applied post-hoc
    args
  })
  
  # Save ppe_efficacy per particle for Policy A post-hoc calculation in extract script
  ppe_lookup <- data.frame(
    scenario     = sc_name,
    particle_id  = seq_len(N_PARTICLES),
    ppe_efficacy = theta$ppe_efficacy
  )
  saveRDS(ppe_lookup,
          file.path(OUT_BASE, sc_name,
                    sprintf("ppe_efficacy_lookup_%s.rds", sc_name)))
  
  list(
    sc_name            = sc_name,
    sc_idx             = match(sc_name, names(SCENARIOS)),
    theta              = theta,
    base_particle_args = base_particle_args
  )
})
names(scenario_setups) <- names(SCENARIOS)

# =============================================================================
# 3. Build worker job list
# =============================================================================
jobs <- lapply(seq_len(N_WORKERS), function(w) {
  p_from <- (w - 1L) * PARTICLES_PER_WORKER + 1L
  p_to   <-  w       * PARTICLES_PER_WORKER
  list(worker_id = w, p_from = p_from, p_to = p_to)
})

message(sprintf(
  "%d workers | %d particles each | %d scenarios x %d reps = %d total runs",
  N_WORKERS, PARTICLES_PER_WORKER,
  length(SCENARIOS), N_REPS,
  length(SCENARIOS) * N_PARTICLES * N_REPS
))

# =============================================================================
# 4. Run in parallel
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
        fname    <- sprintf("%s_p%03d_r%02d.rds", sc_name, p, r)
        out_path <- file.path(OUT_BASE, sc_name, fname)
        if (file.exists(out_path)) next
        
        seed <- SEED_BASE +
          (sc_idx - 1L) * N_PARTICLES * N_REPS +
          (p      - 1L) * N_REPS +
          (r      - 1L)
        
        args <- setup$base_particle_args[[p]]
        
        out          <- NULL
        retry        <- 0L
        current_seed <- seed
        repeat {
          args$seed <- current_seed
          out       <- do.call(fiber::branching_process_main, args)
          tdf_all   <- out$tdf[!is.na(out$tdf$time_infection_absolute), ]
          n_deaths  <- sum(!is.na(tdf_all$outcome) & tdf_all$outcome)
          
          if (n_deaths >= TAKEOFF_DEATH_THRESHOLD) break
          
          retry        <- retry + 1L
          current_seed <- current_seed + 1L
          if (retry >= MAX_RETRIES) {
            message(sprintf(
              "  [w%02d | %s | p%03d | r%02d] WARNING: max retries reached (n_deaths=%d)",
              job$worker_id, sc_name, p, r, n_deaths
            ))
            break
          }
        }
        
        saveRDS(out, out_path)
      }
    }
  }
  
  message(sprintf("Worker %d done (particles %d-%d).", job$worker_id, job$p_from, job$p_to))
  NULL
},
future.globals = list(
  scenario_setups         = scenario_setups,
  OUT_BASE                = OUT_BASE,
  N_PARTICLES             = N_PARTICLES,
  N_REPS                  = N_REPS,
  SEED_BASE               = SEED_BASE,
  TAKEOFF_DEATH_THRESHOLD = TAKEOFF_DEATH_THRESHOLD,
  MAX_RETRIES             = MAX_RETRIES
),
future.packages = c("fiber", "here"),
future.seed     = TRUE)

plan(sequential)

elapsed <- proc.time() - t_start
n_files <- length(list.files(file.path(OUT_BASE, "WestAfrica"), pattern = "\\.rds$"))
message(sprintf("Done in %.1f minutes. Total files: %d / %d expected.",
                elapsed["elapsed"] / 60, n_files, N_PARTICLES * N_REPS))