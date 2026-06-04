# =============================================================================
#
# Jobs are split by posterior particle — each worker handles a slice of
# particles and runs all (coverage_scenario x arm x rep) combinations for them.
#
# N_WORKERS must be a divisor of N_PARTICLES (e.g. 1,2,4,5,10,20,25,50,100).
# Each worker runs: (N_PARTICLES / N_WORKERS) x 3 coverage x arms x 5 reps
#
# No manual settings needed — just set N_WORKERS and source().
# =============================================================================

library(here)
library(future)
library(future.apply)
library(fiber)

source(here("functions", "setup_model_parameters.R"))
source(here("functions", "calculate_model_approx_r0.R"))
source(here("functions", "abc_calibration_functions_common.R"))
source(here("functions", "abc_posterior.R"))

env_hcw <- new.env(parent = globalenv())
env_npi <- new.env(parent = globalenv())
sys.source(here("functions", "abc_calibration_functions_hcwRisk.R"), envir = env_hcw)
sys.source(here("functions", "abc_calibration_functions_npi.R"),     envir = env_npi)

# =============================================================================
# Configuration
# =============================================================================

# >>> SET THIS: must be a divisor of N_PARTICLES (1,2,4,5,10,20,25,50,100)
N_WORKERS        <- 10L

N_PARTICLES      <- 10
N_REPS           <- 2
SEEDING_CASES    <- 25
CHECK_FINAL_SIZE <- 200
RESAMPLE_SEED    <- 42L
SEED_BASE        <- 20260601L

stopifnot(N_PARTICLES %% N_WORKERS == 0L)
PARTICLES_PER_WORKER <- N_PARTICLES %/% N_WORKERS

SCENARIOS <- list(
  WestAfrica = list(
    id           = "Worst_WestAfrica",
    fit_type     = "hcw",
    rds          = here("inst", "fiber_ABC_SMC_Worst_WestAfrica_2026-05-26.rds"),
    scenario_csv = here("data-processed", "final_six_scenario_values_original_approach.csv"),
    param_names  = c("R0", "prop_funeral", "hcw_risk_scalar"),
    npi_spec     = NULL
  ),
  DRC = list(
    id           = "Middle_DRC_ConflictSmoothed_PlusPlus",
    fit_type     = "npi",
    rds          = here("outputs", "02_ABC_model_fits_NPI_Eff",
                        "fiber_ABC_SMC_Middle_DRC_ConflictSmoothed_PlusPlus_NPIeff_20260601_154658.rds"),
    scenario_csv = here("data-processed", "final_six_scenario_values_original_approach.csv"),
    param_names  = c("R0", "prop_funeral", "npi_scaler"),
    npi_spec     = list(
      ppe_efficacy = list(min = 0.30, max = 0.90),
      etu_efficacy = list(min = 0.60, max = 0.95)
    )
  )
)

OBV_EFFICACIES <- c(
  baseline = NA,
  obv_50   = 0.50,
  obv_60   = 0.60,
  obv_70   = 0.70,
  obv_80   = 0.80,
  obv_90   = 0.90
)

OBV_ADHERENCE        <- 1.0
OBV_DPC              <- 1
OBV_TARGET_CLASS     <- "HCW"
OBV_TARGET_LOCATIONS <- "hospital"

COVERAGE_SCENARIOS <- c("full", "ramp_high", "ramp_low")

COVERAGE_SPECS <- list(
  full      = list(type = "scalar", value = 1.0),
  ramp_high = list(type = "ramp", times  = c(0, 30, 60, 90),
                   values = c(0.20, 0.47, 0.73, 1.00)),
  ramp_low  = list(type = "ramp", times  = c(0, 30, 60, 90),
                   values = c(0.20, 0.30, 0.40, 0.50))
)

for (cs in COVERAGE_SCENARIOS) {
  dir.create(here("outputs", "simulation_fig1to3", cs),
             recursive = TRUE, showWarnings = FALSE)
}

check_model_function_version()

# =============================================================================
# Pre-compute particle args in main session
# =============================================================================
message("Pre-computing particle args...")

scenario_setups <- lapply(names(SCENARIOS), function(sc_name) {
  sc <- SCENARIOS[[sc_name]]
  message(sprintf("  %s (%s fit)...", sc_name, sc$fit_type))
  
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
                           check_final_size = CHECK_FINAL_SIZE)
  )
  
  inv <- compute_R0_invariants(args = mp$args, n = 50000, seed = 42L)
  
  build_args_hcw <- env_hcw$build_abc_model_args
  build_args_npi <- env_npi$build_abc_model_args
  
  particle_args <- lapply(seq_len(N_PARTICLES), function(p) {
    if (sc$fit_type == "hcw") {
      D     <- D_from_invariants(
        inv,
        etu_efficacy                         = mp$args$etu_efficacy,
        general_hospital_quarantine_efficacy = mp$args$general_hospital_quarantine_efficacy
      )
      F_fun <- F_from_invariants(inv, safe_funeral_efficacy = mp$args$safe_funeral_efficacy)
      args  <- build_args_hcw(
        R0              = theta$R0[p],
        prop_funeral    = theta$prop_funeral[p],
        hcw_risk_scalar = theta$hcw_risk_scalar[p],
        base            = mp$base_args,
        tv              = mp$tv_args,
        D               = D,
        F_fun           = F_fun,
        seeding_cases   = SEEDING_CASES
      )
    } else {
      args <- build_args_npi(
        R0           = theta$R0[p],
        prop_funeral = theta$prop_funeral[p],
        npi_scaler   = theta$npi_scaler[p],
        base         = mp$base_args,
        tv           = mp$tv_args,
        invariants   = inv,
        npi_spec     = sc$npi_spec,
        general_hospital_quarantine_efficacy = DEFAULT_SCALAR_INPUTS$general_hospital_quarantine_efficacy,
        safe_funeral_efficacy                = DEFAULT_SCALAR_INPUTS$safe_funeral_efficacy,
        seeding_cases = SEEDING_CASES
      )
    }
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
# N_WORKERS jobs, each handles PARTICLES_PER_WORKER particles
# x 2 scenarios x 3 coverage x arms x 5 reps
# =============================================================================
jobs <- lapply(seq_len(N_WORKERS), function(w) {
  p_from <- (w - 1L) * PARTICLES_PER_WORKER + 1L
  p_to   <-  w       * PARTICLES_PER_WORKER
  list(worker_id = w, p_from = p_from, p_to = p_to)
})

message(sprintf(
  "%d workers | %d particles each | %d scenarios x %d coverage x arms x %d reps",
  N_WORKERS, PARTICLES_PER_WORKER, length(SCENARIOS),
  length(COVERAGE_SCENARIOS), N_REPS
))

# =============================================================================
# Run in parallel
# =============================================================================
n_workers <- min(N_WORKERS, future::availableCores(), length(jobs))
plan(multisession, workers = n_workers)
t_start <- proc.time()

future_lapply(jobs, function(job) {
  for (p in seq(job$p_from, job$p_to)) {
    for (sc_name in names(scenario_setups)) {
      setup  <- scenario_setups[[sc_name]]
      sc_idx <- setup$sc_idx
      
      for (cs in COVERAGE_SCENARIOS) {
        spec    <- COVERAGE_SPECS[[cs]]
        out_dir <- here::here("outputs", "simulation_fig1to3", cs)
        
        obv_coverage <- if (spec$type == "scalar") {
          spec$value
        } else {
          make_time_varying(times = spec$times, values = spec$values)
        }
        
        arms_to_run <- if (cs == "full") {
          names(OBV_EFFICACIES)
        } else {
          names(OBV_EFFICACIES)[names(OBV_EFFICACIES) != "baseline"]
        }
        
        for (arm_name in arms_to_run) {
          arm_idx <- match(arm_name, names(OBV_EFFICACIES))
          
          for (r in seq_len(N_REPS)) {
            fname    <- sprintf("%s_p%03d_%s_r%d.rds", sc_name, p, arm_name, r)
            out_path <- file.path(out_dir, fname)
            if (file.exists(out_path)) next
            
            seed <- SEED_BASE +
              (sc_idx  - 1L) * N_PARTICLES * length(OBV_EFFICACIES) * N_REPS +
              (p       - 1L) * length(OBV_EFFICACIES) * N_REPS +
              (arm_idx - 1L) * N_REPS +
              (r       - 1L)
            
            args      <- setup$particle_args[[p]]
            args$seed <- seed
            
            if (!is.na(OBV_EFFICACIES[arm_name])) {
              args$obv_pep_enabled          <- TRUE
              args$obv_pep_coverage         <- obv_coverage
              args$obv_pep_adherence        <- OBV_ADHERENCE
              args$obv_pep_efficacy         <- OBV_EFFICACIES[arm_name]
              args$obv_pep_dpc              <- OBV_DPC
              args$obv_pep_target_class     <- OBV_TARGET_CLASS
              args$obv_pep_target_locations <- OBV_TARGET_LOCATIONS
            }
            
            out   <- do.call(fiber::branching_process_main, args)
            tdf   <- out$tdf
            cases <- tdf[!is.na(tdf$time_infection_absolute), ]
            
            is_hcw <- cases$class == "HCW"
            died   <- !is.na(cases$outcome) & cases$outcome
            
            result <- list(
              scenario          = sc_name,
              particle_id       = p,
              arm               = arm_name,
              rep               = r,
              coverage_scenario = cs,
              R0                = setup$theta$R0[p],
              prop_funeral      = setup$theta$prop_funeral[p],
              third_param       = setup$theta[[names(setup$theta)[3]]][p],
              tdf               = cases,
              n_infections      = nrow(cases),
              n_hcw_infections  = sum(is_hcw),
              n_deaths          = sum(died),
              n_hcw_deaths      = sum(died & is_hcw),
              duration          = max(cases$time_outcome_absolute, na.rm = TRUE),
              num_treated       = out$sim_info$obv_pep_num_treated
            )
            
            saveRDS(result, out_path)
          }
        }
      }
    }
  }
  
  message(sprintf("Worker %d done (particles %d-%d).",
                  job$worker_id, job$p_from, job$p_to))
  NULL
},
future.globals = list(
  scenario_setups      = scenario_setups,
  COVERAGE_SCENARIOS   = COVERAGE_SCENARIOS,
  COVERAGE_SPECS       = COVERAGE_SPECS,
  N_PARTICLES          = N_PARTICLES,
  N_REPS               = N_REPS,
  SEED_BASE            = SEED_BASE,
  OBV_EFFICACIES       = OBV_EFFICACIES,
  OBV_ADHERENCE        = OBV_ADHERENCE,
  OBV_DPC              = OBV_DPC,
  OBV_TARGET_CLASS     = OBV_TARGET_CLASS,
  OBV_TARGET_LOCATIONS = OBV_TARGET_LOCATIONS
),
future.packages = c("fiber", "here"),
future.seed     = TRUE)

plan(sequential)

elapsed <- proc.time() - t_start
n_files <- sum(sapply(COVERAGE_SCENARIOS, function(cs)
  length(list.files(here("outputs", "simulation_fig1to3", cs),
                    pattern = "\\.rds$"))))
message(sprintf("Done in %.1f minutes. Total files: %d",
                elapsed["elapsed"] / 60, n_files))