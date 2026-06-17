# =============================================================================
# 01_analysis_figure3new_conflict_dpc.R
#
# Runs antiviral simulations under four coverage/DPC scenarios, each combined
# with three efficacy curves derived from the DPC-efficacy lookup table
# (data-processed/DPC_fixed_efficacy_varied_d50.rds).
#
# Coverage/DPC scenarios:
#   with_conflict    : coverage and DPC derived from sdb$value (conflict-disrupted)
#   without_conflict : same up to peak day, then held flat (no disruption)
#   optimistic       : coverage = 100%, DPC = 0 throughout
#   no_pep           : coverage = 0% throughout (no antiviral)
#
# Efficacy arms (from curve_d50_dat, interpolated against the DPC(t) curve):
#   mid : efficacy
#   lo  : eighty_efficacy_lo
#   hi  : eighty_efficacy_hi
#
# Arms (10 total):
#   with_conflict_mid, with_conflict_lo, with_conflict_hi
#   without_conflict_mid, without_conflict_lo, without_conflict_hi
#   optimistic_mid, optimistic_lo, optimistic_hi
#   no_pep
#
# obv_pep_dpc is passed the same DPC(t) function used to derive obv_pep_efficacy,
# so the per-individual treatment delay and the DPC-dependent efficacy are both
# driven by the same underlying coverage/DPC scenario curve.
#
# Each RDS file contains the full fiber output (out$tdf + out$prevented_completed
# + out$sim_info).
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
SEED_BASE     <- 20260703L
HCW_BASE_PROB <- 0.25

TAKEOFF_DEATH_THRESHOLD <- 100L
MAX_RETRIES             <- 50L

stopifnot(N_PARTICLES %% N_WORKERS == 0L)
PARTICLES_PER_WORKER <- N_PARTICLES %/% N_WORKERS

# =============================================================================
# 1. Build coverage(t) and dpc(t) functions for each scenario
# =============================================================================
sdb <- readRDS(here("data-processed", "SDB_communityDeath_blended.rds"))

# With conflict: direct transform of sdb$value
sdb$coverage_conflict <- sdb$value * 80 / max(sdb$value)
sdb$dpc_conflict       <- 1 + 4 * (1 - sdb$value / max(sdb$value))

# Find peak coverage day restricted to day < 200
sub <- sdb[sdb$day < 200, ]
peak_row <- sub[which.max(sub$coverage_conflict), ]
peak_day <- peak_row$day

# Hold DPC at 1 up until peak_day (with conflict)
sdb$dpc_conflict[sdb$day <= peak_day] <- 1

# Without conflict: coverage same up to peak_day, then flat; DPC held at 1 throughout
sdb$coverage_noconflict <- sdb$coverage_conflict
sdb$coverage_noconflict[sdb$day > peak_day] <- peak_row$coverage_conflict
sdb$dpc_noconflict <- 1

message(sprintf("peak_day = %d, peak_coverage = %.2f", peak_day, peak_row$coverage_conflict))

# Build interpolating step functions from the day-indexed sdb table.
# Coverage values in sdb are on a 0-100 scale; obv_pep_coverage expects a
# 0-1 proportion, so divide by 100 here.
make_step_fn <- function(days, values) {
  function(t) {
    approx(x = days, y = values, xout = t, rule = 2)$y
  }
}

cov_with_conflict_fn    <- make_step_fn(sdb$day, sdb$coverage_conflict / 100)
cov_without_conflict_fn <- make_step_fn(sdb$day, sdb$coverage_noconflict / 100)
dpc_with_conflict_fn    <- make_step_fn(sdb$day, sdb$dpc_conflict)
dpc_without_conflict_fn <- make_step_fn(sdb$day, sdb$dpc_noconflict)

cov_optimistic_fn <- function(t) rep(1.0, length(t))
dpc_optimistic_fn <- function(t) rep(0,   length(t))

cov_no_pep_fn <- function(t) rep(0.0, length(t))
dpc_no_pep_fn <- function(t) rep(0,   length(t))   # irrelevant when coverage = 0

# =============================================================================
# 2. Build efficacy(t) functions by interpolating DPC(t) against the lookup table
# =============================================================================
curve_d50_dat <- readRDS(here("data-processed", "DPC_fixed_efficacy_varied_d50.rds"))

make_efficacy_fn_from_dpc <- function(dpc_t_fn, efficacy_col) {
  function(t) {
    dpc_vals <- dpc_t_fn(t)
    approx(x = curve_d50_dat$dpc, y = curve_d50_dat[[efficacy_col]],
           xout = dpc_vals, rule = 2)$y
  }
}

eff_with_conflict_mid <- make_efficacy_fn_from_dpc(dpc_with_conflict_fn, "efficacy")
eff_with_conflict_lo  <- make_efficacy_fn_from_dpc(dpc_with_conflict_fn, "eighty_efficacy_lo")
eff_with_conflict_hi  <- make_efficacy_fn_from_dpc(dpc_with_conflict_fn, "eighty_efficacy_hi")

eff_without_conflict_mid <- make_efficacy_fn_from_dpc(dpc_without_conflict_fn, "efficacy")
eff_without_conflict_lo  <- make_efficacy_fn_from_dpc(dpc_without_conflict_fn, "eighty_efficacy_lo")
eff_without_conflict_hi  <- make_efficacy_fn_from_dpc(dpc_without_conflict_fn, "eighty_efficacy_hi")

eff_optimistic_mid <- make_efficacy_fn_from_dpc(dpc_optimistic_fn, "efficacy")
eff_optimistic_lo  <- make_efficacy_fn_from_dpc(dpc_optimistic_fn, "eighty_efficacy_lo")
eff_optimistic_hi  <- make_efficacy_fn_from_dpc(dpc_optimistic_fn, "eighty_efficacy_hi")

# =============================================================================
# 3. Arm definitions: coverage_fn, dpc_fn, efficacy_fn triples
# =============================================================================
ARM_DEFS <- list(
  with_conflict_mid    = list(cov_fn = cov_with_conflict_fn,    dpc_fn = dpc_with_conflict_fn,    eff_fn = eff_with_conflict_mid),
  with_conflict_lo     = list(cov_fn = cov_with_conflict_fn,    dpc_fn = dpc_with_conflict_fn,    eff_fn = eff_with_conflict_lo),
  with_conflict_hi     = list(cov_fn = cov_with_conflict_fn,    dpc_fn = dpc_with_conflict_fn,    eff_fn = eff_with_conflict_hi),
  without_conflict_mid = list(cov_fn = cov_without_conflict_fn, dpc_fn = dpc_without_conflict_fn, eff_fn = eff_without_conflict_mid),
  without_conflict_lo  = list(cov_fn = cov_without_conflict_fn, dpc_fn = dpc_without_conflict_fn, eff_fn = eff_without_conflict_lo),
  without_conflict_hi  = list(cov_fn = cov_without_conflict_fn, dpc_fn = dpc_without_conflict_fn, eff_fn = eff_without_conflict_hi),
  optimistic_mid       = list(cov_fn = cov_optimistic_fn,       dpc_fn = dpc_optimistic_fn,       eff_fn = eff_optimistic_mid),
  optimistic_lo        = list(cov_fn = cov_optimistic_fn,       dpc_fn = dpc_optimistic_fn,       eff_fn = eff_optimistic_lo),
  optimistic_hi        = list(cov_fn = cov_optimistic_fn,       dpc_fn = dpc_optimistic_fn,       eff_fn = eff_optimistic_hi),
  no_pep               = list(cov_fn = cov_no_pep_fn,           dpc_fn = dpc_no_pep_fn,           eff_fn = function(t) rep(0, length(t)))
)

ARM_NAMES <- names(ARM_DEFS)
message(sprintf("Total arms: %d", length(ARM_NAMES)))
print(ARM_NAMES)

# =============================================================================
# 4. Scenarios
# =============================================================================
SCENARIOS <- list(
  DRC = list(
    id           = "Middle_DRC_ConflictSmoothed_PlusPlus",
    rds          = here("outputs", "02_ABC_model_fits_Final",
                        "fiber_ABC_SMC_Middle_DRC_ConflictSmoothed_PlusPlus_Decoupled_20260607_215621_check_NP5_NS4_NBREPS_30_NBSIMUL_590.RDS"),
    scenario_csv = here("data-processed", "final_six_scenario_values_original_approach.csv"),
    param_names  = c("R0", "prop_funeral", "etu_efficacy", "ppe_efficacy", "hcw_risk_scalar"),
    check_final_size = 10000
  )
)

OUT_BASE <- here("outputs", "simulation", "conflict_dpc")
for (arm_name in ARM_NAMES) {
  dir.create(file.path(OUT_BASE, arm_name), recursive = TRUE, showWarnings = FALSE)
}

check_model_function_version()

# =============================================================================
# 5. Pre-compute particle args per scenario
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
    sc_name            = sc_name,
    sc_idx             = match(sc_name, names(SCENARIOS)),
    theta              = theta,
    base_particle_args = base_particle_args
  )
})
names(scenario_setups) <- names(SCENARIOS)

# =============================================================================
# 6. Build job list
# =============================================================================
jobs <- lapply(seq_len(N_WORKERS), function(w) {
  p_from <- (w - 1L) * PARTICLES_PER_WORKER + 1L
  p_to   <-  w       * PARTICLES_PER_WORKER
  list(worker_id = w, p_from = p_from, p_to = p_to)
})

message(sprintf(
  "%d workers | %d particles each | %d scenarios x %d arms x %d reps = %d total runs",
  N_WORKERS, PARTICLES_PER_WORKER,
  length(SCENARIOS), length(ARM_NAMES), N_REPS,
  length(SCENARIOS) * length(ARM_NAMES) * N_PARTICLES * N_REPS
))

# =============================================================================
# 7. Run in parallel
# =============================================================================
n_workers <- min(N_WORKERS, future::availableCores(), length(jobs))
plan(multisession, workers = n_workers)
t_start <- proc.time()

future_lapply(jobs, function(job) {
  for (p in seq(job$p_from, job$p_to)) {
    for (sc_name in names(scenario_setups)) {
      setup  <- scenario_setups[[sc_name]]
      sc_idx <- setup$sc_idx
      
      for (arm_idx in seq_along(ARM_NAMES)) {
        arm_name <- ARM_NAMES[arm_idx]
        arm_def  <- ARM_DEFS[[arm_name]]
        
        for (r in seq_len(N_REPS)) {
          fname    <- sprintf("%s_p%03d_r%02d.rds", sc_name, p, r)
          out_path <- file.path(OUT_BASE, arm_name, fname)
          if (file.exists(out_path)) next
          
          seed <- SEED_BASE +
            (sc_idx  - 1L) * N_PARTICLES * length(ARM_NAMES) * N_REPS +
            (arm_idx - 1L) * N_PARTICLES * N_REPS +
            (p       - 1L) * N_REPS +
            (r       - 1L)
          
          args <- setup$base_particle_args[[p]]
          
          if (arm_name == "no_pep") {
            args$obv_pep_enabled <- FALSE
          } else {
            args$obv_pep_enabled          <- TRUE
            args$obv_pep_coverage         <- arm_def$cov_fn
            args$obv_pep_adherence        <- 1.0
            args$obv_pep_dpc              <- arm_def$dpc_fn
            args$obv_pep_efficacy         <- arm_def$eff_fn
            args$obv_pep_target_class     <- "HCW"
            args$obv_pep_target_locations <- c("hospital", "community", "funeral")
          }
          
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
  ARM_NAMES               = ARM_NAMES,
  ARM_DEFS                = ARM_DEFS,
  sdb                     = sdb,
  curve_d50_dat           = curve_d50_dat,
  make_step_fn              = make_step_fn,
  make_efficacy_fn_from_dpc = make_efficacy_fn_from_dpc,
  cov_with_conflict_fn      = cov_with_conflict_fn,
  cov_without_conflict_fn   = cov_without_conflict_fn,
  dpc_with_conflict_fn      = dpc_with_conflict_fn,
  dpc_without_conflict_fn   = dpc_without_conflict_fn,
  cov_optimistic_fn         = cov_optimistic_fn,
  dpc_optimistic_fn         = dpc_optimistic_fn,
  cov_no_pep_fn              = cov_no_pep_fn,
  dpc_no_pep_fn              = dpc_no_pep_fn,
  eff_with_conflict_mid      = eff_with_conflict_mid,
  eff_with_conflict_lo       = eff_with_conflict_lo,
  eff_with_conflict_hi       = eff_with_conflict_hi,
  eff_without_conflict_mid   = eff_without_conflict_mid,
  eff_without_conflict_lo    = eff_without_conflict_lo,
  eff_without_conflict_hi    = eff_without_conflict_hi,
  eff_optimistic_mid         = eff_optimistic_mid,
  eff_optimistic_lo          = eff_optimistic_lo,
  eff_optimistic_hi          = eff_optimistic_hi,
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
n_files <- sum(sapply(ARM_NAMES, function(a)
  length(list.files(file.path(OUT_BASE, a), pattern = "\\.rds$"))))
message(sprintf("Done in %.1f minutes. Total files: %d", elapsed["elapsed"] / 60, n_files))