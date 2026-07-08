# =============================================================================
#
# Sensitivity analysis: what if antiviral efficacy stays high for LONGER
# after exposure before it starts to decline (i.e. a more optimistic
# efficacy-vs-DPC curve)?
#
# curve_d50_dat gives efficacy as a function of DPC (days post-exposure to
# treatment). We build shifted versions of this curve where
#   efficacy_shifted(dpc) = efficacy_original(max(dpc - shift, 0))
# i.e. the same efficacy value that originally applied at a given DPC now
# applies `shift` days later -- the drop-off is delayed by `shift` days.
# We sweep shift = 2, 4, 6, 8, 10 days.
#
# The conflict's impact on COVERAGE and DPC (delivery delay) itself is left
# at the baseline intensity (coverage max 80%, DPC max +4 days) -- only the
# efficacy-vs-DPC curve changes here. Everything else (scenario, arms,
# particles, reps, seeds) is identical to the baseline script, so results
# are directly comparable across shift levels.
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
SEED_BASE     <- 20260703L   # kept identical across shifts -> matched-seed comparisons
HCW_BASE_PROB <- 0.25

TAKEOFF_DEATH_THRESHOLD <- 100L
MAX_RETRIES             <- 50L

stopifnot(N_PARTICLES %% N_WORKERS == 0L)
PARTICLES_PER_WORKER <- N_PARTICLES %/% N_WORKERS

# =============================================================================
# 1. Build coverage(t) and dpc(t) functions -- baseline conflict intensity,
#    unchanged from 01_analysis_figure3new_conflict_dpc.R
# =============================================================================
sdb <- readRDS(here("data-processed", "SDB_communityDeath_blended.rds"))

sdb_tweaked <- sdb$value
idx_150_325 <- sdb$day >= 150 & sdb$day <= 325
idx_325_400 <- sdb$day >  325 & sdb$day <= 400

rescale_sdb_segment <- function(day_out, orig_from, orig_to) {
  t        <- (day_out - day_out[1]) / (day_out[length(day_out)] - day_out[1])
  day_orig <- orig_from + t * (orig_to - orig_from)
  approx(sdb$day, sdb$value, xout = day_orig, rule = 2)$y
}

sdb_tweaked[idx_150_325] <- rescale_sdb_segment(sdb$day[idx_150_325], 150, 200)
sdb_tweaked[idx_325_400] <- rescale_sdb_segment(sdb$day[idx_325_400], 200, 350)
sdb$value_tweaked <- sdb_tweaked

sdb$coverage_conflict <- sdb$value_tweaked * 80 / max(sdb$value_tweaked)
sdb$dpc_conflict      <- 1 + 4 * (1 - (sdb$value_tweaked / max(sdb$value_tweaked)))

sub      <- sdb[sdb$day < 200, ]
peak_row <- sub[which.max(sub$coverage_conflict), ]
peak_day <- peak_row$day

sdb$dpc_conflict[sdb$day <= peak_day] <- 1

sdb$coverage_flat <- sdb$coverage_conflict
sdb$coverage_flat[sdb$day > peak_day] <- peak_row$coverage_conflict

message(sprintf("peak_day = %d, peak_coverage = %.2f", peak_day, peak_row$coverage_conflict))

make_step_fn <- function(days, values) {
  function(t) approx(x = days, y = values, xout = t, rule = 2)$y
}

cov_with_conflict_fn <- make_step_fn(sdb$day, sdb$coverage_conflict / 100)
cov_flat_fn          <- make_step_fn(sdb$day, sdb$coverage_flat / 100)
cov_perfect_fn        <- function(t) rep(1.0, length(t))
cov_no_pep_fn         <- function(t) rep(0.0, length(t))

dpc_with_conflict_fn <- make_step_fn(sdb$day, sdb$dpc_conflict)
dpc_flat_fn          <- make_step_fn(sdb$day, rep(1, nrow(sdb)))
dpc_perfect_fn        <- function(t) rep(0, length(t))
dpc_no_pep_fn         <- function(t) rep(0, length(t))

# =============================================================================
# 2. Efficacy curve + shifted variants
#
# obv_pep_efficacy receives individual DPC directly from fiber, so shifting
# the lookup here delays the point at which efficacy starts declining by
# `shift` days, without needing to touch coverage/DPC-delay curves at all.
# =============================================================================
curve_d50_dat <- readRDS(here("data-processed", "DPC_fixed_efficacy_varied_d50.rds"))

make_efficacy_fn_shifted <- function(efficacy_col, shift) {
  force(efficacy_col)
  force(shift)
  function(dpc) {
    approx(x = curve_d50_dat$dpc, y = curve_d50_dat[[efficacy_col]],
           xout = pmax(dpc - shift, 0), rule = 2)$y
  }
}

SHIFTS <- c(2, 4, 6, 8, 10)

# =============================================================================
# 3. Arm definitions, parameterised by the shift-dependent efficacy functions
# =============================================================================
build_arm_defs <- function(eff_mid, eff_lo, eff_hi) {
  list(
    no_pep = list(cov_fn = cov_no_pep_fn, dpc_fn = dpc_no_pep_fn,
                  eff_fn = function(dpc) rep(0, length(dpc))),
    
    with_conflict_mid = list(cov_fn = cov_with_conflict_fn, dpc_fn = dpc_with_conflict_fn, eff_fn = eff_mid),
    with_conflict_lo  = list(cov_fn = cov_with_conflict_fn, dpc_fn = dpc_with_conflict_fn, eff_fn = eff_lo),
    with_conflict_hi  = list(cov_fn = cov_with_conflict_fn, dpc_fn = dpc_with_conflict_fn, eff_fn = eff_hi),
    
    cov_conflict_mid = list(cov_fn = cov_with_conflict_fn, dpc_fn = dpc_flat_fn, eff_fn = eff_mid),
    cov_conflict_lo  = list(cov_fn = cov_with_conflict_fn, dpc_fn = dpc_flat_fn, eff_fn = eff_lo),
    cov_conflict_hi  = list(cov_fn = cov_with_conflict_fn, dpc_fn = dpc_flat_fn, eff_fn = eff_hi),
    
    dpc_conflict_mid = list(cov_fn = cov_flat_fn, dpc_fn = dpc_with_conflict_fn, eff_fn = eff_mid),
    dpc_conflict_lo  = list(cov_fn = cov_flat_fn, dpc_fn = dpc_with_conflict_fn, eff_fn = eff_lo),
    dpc_conflict_hi  = list(cov_fn = cov_flat_fn, dpc_fn = dpc_with_conflict_fn, eff_fn = eff_hi),
    
    optimistic_mid = list(cov_fn = cov_perfect_fn, dpc_fn = dpc_perfect_fn, eff_fn = eff_mid),
    optimistic_lo  = list(cov_fn = cov_perfect_fn, dpc_fn = dpc_perfect_fn, eff_fn = eff_lo),
    optimistic_hi  = list(cov_fn = cov_perfect_fn, dpc_fn = dpc_perfect_fn, eff_fn = eff_hi)
  )
}

ARM_NAMES <- c(
  "no_pep",
  "with_conflict_mid", "with_conflict_lo", "with_conflict_hi",
  "cov_conflict_mid",  "cov_conflict_lo",  "cov_conflict_hi",
  "dpc_conflict_mid",  "dpc_conflict_lo",  "dpc_conflict_hi",
  "optimistic_mid",    "optimistic_lo",    "optimistic_hi"
)

# =============================================================================
# 4. Scenario + particle-arg precompute (identical to baseline, done ONCE
#    and reused across all shift levels -- this is the expensive step, and
#    it doesn't depend on the efficacy curve shift at all)
# =============================================================================
SCENARIO_CSV <- here("data-processed", "final_six_scenario_values_original_approach.csv")

SCENARIOS <- list(
  DRC_conflict = list(
    id               = "Middle_DRC_ConflictSmoothed_PlusPlus",
    rds              = here("outputs", "02_ABC_model_fits_Final",
                            "fiber_ABC_SMC_Middle_DRC_ConflictSmoothed_PlusPlus_Decoupled_20260607_215621_check_NP5_NS4_NBREPS_30_NBSIMUL_590.RDS"),
    scenario_csv     = SCENARIO_CSV,
    param_names      = c("R0", "prop_funeral", "etu_efficacy", "ppe_efficacy", "hcw_risk_scalar"),
    check_final_size = 10000,
    arms             = ARM_NAMES
  )
)

check_model_function_version()

message("Pre-computing particle args (shared across shift levels)...")

scenario_setups <- lapply(names(SCENARIOS), function(sc_name) {
  sc <- SCENARIOS[[sc_name]]
  message(sprintf("  %s (scenario_id: %s)...", sc_name, sc$id))
  
  res       <- readRDS(sc$rds)
  posterior <- data.frame(weight = res$weights)
  for (j in seq_along(sc$param_names)) {
    posterior[[sc$param_names[j]]] <- res$param[, j]
  }
  theta <- downsample_posterior(posterior, n_sets = N_PARTICLES,
                                seed = RESAMPLE_SEED, param_names = sc$param_names)
  
  scenario_matrix <- read_scenario_matrix(sc$scenario_csv)
  
  mp  <- make_model_parameters(
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
    arms               = sc$arms,
    theta              = theta,
    base_particle_args = base_particle_args
  )
})
names(scenario_setups) <- names(SCENARIOS)

total_runs <- length(ARM_NAMES) * N_PARTICLES * N_REPS
message(sprintf(
  "%d workers | %d particles each | %d total runs PER shift level",
  N_WORKERS, PARTICLES_PER_WORKER, total_runs
))

# =============================================================================
# 5. Parallel run for one shift level
# =============================================================================
run_condition <- function(OUT_BASE, ARM_DEFS, ARM_NAMES, scenario_setups) {
  for (arm_name in ARM_NAMES) {
    dir.create(file.path(OUT_BASE, arm_name), recursive = TRUE, showWarnings = FALSE)
  }
  
  jobs <- lapply(seq_len(N_WORKERS), function(w) {
    p_from <- (w - 1L) * PARTICLES_PER_WORKER + 1L
    p_to   <-  w       * PARTICLES_PER_WORKER
    list(worker_id = w, p_from = p_from, p_to = p_to)
  })
  
  n_workers <- min(N_WORKERS, future::availableCores(), length(jobs))
  plan(multisession, workers = n_workers)
  t_start <- proc.time()
  
  future_lapply(jobs, function(job) {
    for (p in seq(job$p_from, job$p_to)) {
      for (sc_name in names(scenario_setups)) {
        setup <- scenario_setups[[sc_name]]
        
        for (arm_name in setup$arms) {
          arm_idx <- which(ARM_NAMES == arm_name)
          arm_def <- ARM_DEFS[[arm_name]]
          
          for (r in seq_len(N_REPS)) {
            fname    <- sprintf("%s_p%03d_r%02d.rds", sc_name, p, r)
            out_path <- file.path(OUT_BASE, arm_name, fname)
            if (file.exists(out_path)) next
            
            seed <- SEED_BASE +
              (which(names(scenario_setups) == sc_name) - 1L) * N_PARTICLES * length(ARM_NAMES) * N_REPS +
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
}

# =============================================================================
# 6. Sweep all shift levels
# =============================================================================
for (shift in SHIFTS) {
  message(sprintf("\n===== Efficacy delay shift: +%d days =====", shift))
  
  eff_mid_s <- make_efficacy_fn_shifted("efficacy",            shift)
  eff_lo_s  <- make_efficacy_fn_shifted("eighty_efficacy_lo",  shift)
  eff_hi_s  <- make_efficacy_fn_shifted("eighty_efficacy_hi",  shift)
  
  ARM_DEFS <- build_arm_defs(eff_mid_s, eff_lo_s, eff_hi_s)
  OUT_BASE <- here("outputs", "simulation", sprintf("conflict_dpc_sens_delay%02d", shift))
  
  run_condition(OUT_BASE, ARM_DEFS, ARM_NAMES, scenario_setups)
}

message("\nEfficacy-delay sensitivity sweep complete.")