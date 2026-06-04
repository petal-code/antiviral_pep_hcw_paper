# =============================================================================
# 01_analysis_figure_4.R
#
#   Efficacy : 50, 60, 70, 80, 90%  (5 levels)
#   Coverage : 10, 30, 50, 70, 90%  (5 levels, fixed scalar)
#   => 25 combinations x 100 particles x 2 scenarios x 5 reps = 25,000 runs
#
# N_WORKERS must be a divisor of N_PARTICLES (1,2,4,5,10,20,25,50,100).
# Each worker handles a particle slice x all grid combos x all reps.
#
# Output: outputs/simulation_fig4/
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
SEED_BASE        <- 20260701L   # different from fig1to3 to avoid seed collision

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

OBV_EFFICACY_GRID <- seq(0.50, 0.90, by = 0.10)
OBV_COVERAGE_GRID <- seq(0.10, 0.90, by = 0.20)

OBV_ADHERENCE        <- 1.0
OBV_DPC              <- 1
OBV_TARGET_CLASS     <- "HCW"
OBV_TARGET_LOCATIONS <- "hospital"

GRID_DF <- expand.grid(
  obv_efficacy = OBV_EFFICACY_GRID,
  obv_coverage = OBV_COVERAGE_GRID,
  KEEP.OUT.ATTRS = FALSE
)
GRID_DF$grid_idx <- seq_len(nrow(GRID_DF))

OUT_DIR <- here("outputs", "simulation_fig4")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

message(sprintf("Grid: %d efficacy x %d coverage = %d combinations",
                length(OBV_EFFICACY_GRID), length(OBV_COVERAGE_GRID), nrow(GRID_DF)))
message(sprintf("Total simulations: %d",
                nrow(GRID_DF) * N_PARTICLES * length(SCENARIOS) * N_REPS))

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
# =============================================================================
jobs <- lapply(seq_len(N_WORKERS), function(w) {
  p_from <- (w - 1L) * PARTICLES_PER_WORKER + 1L
  p_to   <-  w       * PARTICLES_PER_WORKER
  list(worker_id = w, p_from = p_from, p_to = p_to)
})

message(sprintf(
  "%d workers | %d particles each | %d scenarios x %d grid combos x %d reps",
  N_WORKERS, PARTICLES_PER_WORKER, length(SCENARIOS), nrow(GRID_DF), N_REPS
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
      
      for (g in seq_len(nrow(GRID_DF))) {
        eff     <- GRID_DF$obv_efficacy[g]
        cov     <- GRID_DF$obv_coverage[g]
        eff_pct <- round(eff * 100)
        cov_pct <- round(cov * 100)
        
        for (r in seq_len(N_REPS)) {
          fname    <- sprintf("%s_p%03d_eff%02d_cov%02d_r%d.rds",
                              sc_name, p, eff_pct, cov_pct, r)
          out_path <- file.path(OUT_DIR, fname)
          if (file.exists(out_path)) next
          
          seed <- SEED_BASE +
            (sc_idx - 1L) * N_PARTICLES * nrow(GRID_DF) * N_REPS +
            (p      - 1L) * nrow(GRID_DF) * N_REPS +
            (g      - 1L) * N_REPS +
            (r      - 1L)
          
          args                          <- setup$particle_args[[p]]
          args$seed                     <- seed
          args$obv_pep_enabled          <- TRUE
          args$obv_pep_coverage         <- cov
          args$obv_pep_adherence        <- OBV_ADHERENCE
          args$obv_pep_efficacy         <- eff
          args$obv_pep_dpc              <- OBV_DPC
          args$obv_pep_target_class     <- OBV_TARGET_CLASS
          args$obv_pep_target_locations <- OBV_TARGET_LOCATIONS
          
          out   <- do.call(fiber::branching_process_main, args)
          tdf   <- out$tdf
          cases <- tdf[!is.na(tdf$time_infection_absolute), ]
          
          is_hcw <- cases$class == "HCW"
          died   <- !is.na(cases$outcome) & cases$outcome
          
          result <- list(
            scenario         = sc_name,
            particle_id      = p,
            obv_efficacy     = eff,
            obv_coverage     = cov,
            rep              = r,
            R0               = setup$theta$R0[p],
            prop_funeral     = setup$theta$prop_funeral[p],
            third_param      = setup$theta[[names(setup$theta)[3]]][p],
            tdf              = cases,
            n_infections     = nrow(cases),
            n_hcw_infections = sum(is_hcw),
            n_deaths         = sum(died),
            n_hcw_deaths     = sum(died & is_hcw),
            duration         = max(cases$time_outcome_absolute, na.rm = TRUE),
            num_treated      = out$sim_info$obv_pep_num_treated
          )
          
          saveRDS(result, out_path)
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
  GRID_DF              = GRID_DF,
  OUT_DIR              = OUT_DIR,
  N_PARTICLES          = N_PARTICLES,
  N_REPS               = N_REPS,
  SEED_BASE            = SEED_BASE,
  OBV_ADHERENCE        = OBV_ADHERENCE,
  OBV_DPC              = OBV_DPC,
  OBV_TARGET_CLASS     = OBV_TARGET_CLASS,
  OBV_TARGET_LOCATIONS = OBV_TARGET_LOCATIONS
),
future.packages = c("fiber", "here"),
future.seed     = TRUE)

plan(sequential)

elapsed <- proc.time() - t_start
n_files <- length(list.files(OUT_DIR, pattern = "\\.rds$"))
message(sprintf("Done in %.1f minutes. Total files: %d",
                elapsed["elapsed"] / 60, n_files))