# =============================================================================
# 01_analysis_figure3.R
#
# Runs OBV simulations exclusively for ramp coverage combinations.
# Each arm is saved to a separate subdirectory under outputs/simulation/.
#
# Arms (10 total):
#   ramp_high_obv50 ~ ramp_high_obv90  (5) : ramp-high coverage, 50-90%   -- Figure 3
#   ramp_low_obv50  ~ ramp_low_obv90   (5) : ramp-low coverage, 50-90%    -- Figure 3
#
# Coverage splines (consistent with helper COVERAGE_SPECS):
#   ramp_high : clamped cubic spline (0wk,0.20)->(26wk,0.50)->(52wk,0.80)
#   ramp_low  : clamped cubic spline (0wk,0.00)->(26wk,0.25)->(52wk,0.50)
#
# Each RDS file contains the full fiber output (out$tdf + out$prevented_completed
# + out$sim_info). "Without OBV" is reconstructed post-hoc by combining
# tdf + prevented_completed from any arm.
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
N_WORKERS    <- 50L
N_PARTICLES  <- 200
N_REPS       <- 10
SEEDING_CASES <- 25L
RESAMPLE_SEED <- 42L
SEED_BASE     <- 20260701L
HCW_BASE_PROB <- 0.25

TAKEOFF_DEATH_THRESHOLD <- 100L
MAX_RETRIES             <- 50L

stopifnot(N_PARTICLES %% N_WORKERS == 0L)
PARTICLES_PER_WORKER <- N_PARTICLES %/% N_WORKERS

# =============================================================================
# Coverage spline builder (consistent with helper COVERAGE_SPECS)
# =============================================================================
make_clamped_spline_fn <- function(t_knots, y_knots,
                                   deriv_start = 0, deriv_end = 0) {
  n <- length(t_knots)
  h <- diff(t_knots)
  
  rhs <- numeric(n)
  rhs[1] <- 3 * (y_knots[2] - y_knots[1]) / h[1] - 3 * deriv_start
  rhs[n] <- 3 * deriv_end - 3 * (y_knots[n] - y_knots[n - 1]) / h[n - 1]
  for (i in 2:(n - 1)) {
    rhs[i] <- 3 * ((y_knots[i + 1] - y_knots[i]) / h[i] -
                     (y_knots[i]     - y_knots[i - 1]) / h[i - 1])
  }
  
  diag_main <- numeric(n)
  diag_main[1] <- 2 * h[1]
  diag_main[n] <- 2 * h[n - 1]
  for (i in 2:(n - 1)) diag_main[i] <- 2 * (h[i - 1] + h[i])
  
  c_vec <- h; d_vec <- rhs
  c_vec[1] <- c_vec[1] / diag_main[1]
  d_vec[1] <- d_vec[1] / diag_main[1]
  for (i in 2:n) {
    denom <- diag_main[i] - (if (i > 1) h[i - 1] else 0) * c_vec[i - 1]
    if (i < n) c_vec[i] <- h[i] / denom
    d_vec[i] <- (d_vec[i] - (if (i > 1) h[i - 1] else 0) * d_vec[i - 1]) / denom
  }
  
  M <- numeric(n)
  M[n] <- d_vec[n]
  for (i in (n - 1):1) M[i] <- d_vec[i] - c_vec[i] * M[i + 1]
  
  function(t) {
    pmin(pmax(vapply(t, function(tv) {
      if (tv <= t_knots[1])  return(y_knots[1])
      if (tv >= t_knots[n])  return(y_knots[n])
      i  <- min(max(findInterval(tv, t_knots, rightmost.closed = TRUE), 1L), n - 1L)
      hi <- h[i]
      a  <- (t_knots[i + 1] - tv) / hi
      b  <- (tv - t_knots[i])     / hi
      a * y_knots[i] + b * y_knots[i + 1] +
        ((a^3 - a) * M[i] + (b^3 - b) * M[i + 1]) * hi^2 / 6
    }, numeric(1)), 0), 1)
  }
}

# Coverage curve functions (days as input)
COVERAGE_FNS <- list(
  ramp_high = make_clamped_spline_fn(
    t_knots = c(0, 26 * 7, 52 * 7),
    y_knots = c(0.20, 0.50, 0.80)
  ),
  ramp_low  = make_clamped_spline_fn(
    t_knots = c(0, 26 * 7, 52 * 7),
    y_knots = c(0.00, 0.25, 0.50)
  )
)

# =============================================================================
# Arm definitions
#
# ramp_high_obv50~obv90       (5 arms) : ramp-high coverage, 50-90%      -- Figure 3
# ramp_low_obv50~obv90        (5 arms) : ramp-low coverage, 50-90%       -- Figure 3
# =============================================================================
RAMP_EFFICACIES  <- c(0.50, 0.60, 0.70, 0.80, 0.90)

ARMS <- do.call(rbind, lapply(c("ramp_high", "ramp_low"), function(cov) {
  do.call(rbind, lapply(RAMP_EFFICACIES, function(eff) {
    data.frame(arm_name  = sprintf("%s_obv%02d", cov, round(eff * 100)),
               coverage  = cov, efficacy = eff, const_cov = NA_real_,
               stringsAsFactors = FALSE)
  }))
}))

message(sprintf("Total arms: %d (ramp: %d)", nrow(ARMS), nrow(ARMS)))
print(ARMS$arm_name)

# =============================================================================
# Scenarios
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

# Create output directories
OUT_BASE <- here("outputs", "simulation")
for (arm_name in ARMS$arm_name) {
  dir.create(file.path(OUT_BASE, arm_name), recursive = TRUE, showWarnings = FALSE)
}

check_model_function_version()

# =============================================================================
# Pre-compute particle args per scenario
# =============================================================================
message("Pre-computing particle args...")

scenario_setups <- lapply(names(SCENARIOS), function(sc_name) {
  sc <- SCENARIOS[[sc_name]]
  message(sprintf("  %s...", sc_name))
  
  res        <- readRDS(sc$rds)
  posterior  <- data.frame(weight = res$weights)
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
  
  # Base args per particle (no OBV settings yet)
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
  
  list(
    sc_name           = sc_name,
    sc_idx            = match(sc_name, names(SCENARIOS)),
    theta             = theta,
    base_particle_args = base_particle_args
  )
})
names(scenario_setups) <- names(SCENARIOS)

# =============================================================================
# Build job list
# =============================================================================
jobs <- lapply(seq_len(N_WORKERS), function(w) {
  p_from <- (w - 1L) * PARTICLES_PER_WORKER + 1L
  p_to   <-  w       * PARTICLES_PER_WORKER
  list(worker_id = w, p_from = p_from, p_to = p_to)
})

message(sprintf(
  "%d workers | %d particles each | %d scenarios x %d arms x %d reps = %d total runs",
  N_WORKERS, PARTICLES_PER_WORKER,
  length(SCENARIOS), nrow(ARMS), N_REPS,
  length(SCENARIOS) * nrow(ARMS) * N_PARTICLES * N_REPS
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
      
      for (arm_idx in seq_len(nrow(ARMS))) {
        arm      <- ARMS[arm_idx, ]
        efficacy <- arm$efficacy
        arm_name <- arm$arm_name
        
        # Resolve coverage function (Only Ramp Coverage in this script)
        cov_fn <- COVERAGE_FNS[[arm$coverage]]
        
        for (r in seq_len(N_REPS)) {
          fname    <- sprintf("%s_p%03d_r%02d.rds", sc_name, p, r)
          out_path <- file.path(OUT_BASE, arm_name, fname)
          if (file.exists(out_path)) next
          
          seed <- SEED_BASE +
            (sc_idx  - 1L) * N_PARTICLES * nrow(ARMS) * N_REPS +
            (arm_idx - 1L) * N_PARTICLES * N_REPS +
            (p       - 1L) * N_REPS +
            (r       - 1L)
          
          args <- setup$base_particle_args[[p]]
          
          # OBV settings
          args$obv_pep_enabled          <- TRUE
          args$obv_pep_coverage         <- cov_fn   # function(t)
          args$obv_pep_adherence        <- 1.0
          args$obv_pep_dpc              <- 0
          args$obv_pep_efficacy         <- efficacy
          args$obv_pep_target_class     <- "HCW"
          args$obv_pep_target_locations <- c("hospital", "community", "funeral")
          
          # Retry until takeoff
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
                "  [w%02d | %s | %s | p%03d | r%02d] WARNING: max retries (n_deaths=%d)",
                job$worker_id, sc_name, arm_name, p, r, n_deaths
              ))
              break
            }
          }
          
          saveRDS(out, out_path)
          message(sprintf(
            "  [w%02d | %s | %s | p%03d | r%02d] seed=%d retries=%d -> %d deaths, %d prevented",
            job$worker_id, sc_name, arm_name, p, r,
            current_seed, retry,
            sum(!is.na(out$tdf$outcome) & out$tdf$outcome),
            if (!is.null(out$sim_info$obv_pep_num_treated))
              out$sim_info$obv_pep_num_treated$prevented else NA
          ))
        }
      }
    }
  }
  
  message(sprintf("Worker %d done (particles %d-%d).", job$worker_id, job$p_from, job$p_to))
  NULL
},
future.globals = list(
  scenario_setups         = scenario_setups,
  OUT_BASE                = OUT_BASE,
  ARMS                    = ARMS,
  COVERAGE_FNS            = COVERAGE_FNS,
  make_clamped_spline_fn  = make_clamped_spline_fn,
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
n_files <- sum(sapply(ARMS$arm_name, function(a)
  length(list.files(file.path(OUT_BASE, a), pattern = "\\.rds$"))))
message(sprintf("Done in %.1f minutes. Total files: %d", elapsed["elapsed"] / 60, n_files))