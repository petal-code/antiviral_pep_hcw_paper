# =============================================================================
# 01_analysis_figure3new_conflict_dpc.R
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

# ---- Tweak sdb$value in the 150-400 day window only ----
# Shape-preserving x-axis rescale: curve form maintained, trough shifted
# from day 200 to day 325.
#   150-200 (original descent) stretched onto 150-325
#   200-350 (original recovery) squeezed onto 325-400
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

# ---- Coverage and DPC curves derived from tweaked sdb ----
sdb$coverage_conflict <- sdb$value_tweaked * 80 / max(sdb$value_tweaked)
# sdb$dpc_conflict      <- 1 + 9 * (1 - (sdb$value_tweaked / max(sdb$value_tweaked))^2)
sdb$dpc_conflict      <- 1 + 2 * (1 - (sdb$value_tweaked / max(sdb$value_tweaked)))
OUT_BASE <- here("outputs", "simulation", "conflict_dpc_max3")
# Find peak coverage day restricted to day < 200
sub      <- sdb[sdb$day < 200, ]
peak_row <- sub[which.max(sub$coverage_conflict), ]
peak_day <- peak_row$day

# Hold DPC at 1 up until peak_day (before conflict disruption begins)
sdb$dpc_conflict[sdb$day <= peak_day] <- 1

message(sprintf("peak_day = %d, peak_coverage = %.2f", peak_day, peak_row$coverage_conflict))

# ---- Flat coverage: same ramp up to peak_day, held constant thereafter ----
sdb$coverage_flat <- sdb$coverage_conflict
sdb$coverage_flat[sdb$day > peak_day] <- peak_row$coverage_conflict

# ---- Step function builders ----
make_step_fn <- function(days, values) {
  function(t) approx(x = days, y = values, xout = t, rule = 2)$y
}

# Coverage functions (take absolute time t, return coverage proportion)
cov_with_conflict_fn <- make_step_fn(sdb$day, sdb$coverage_conflict / 100)
cov_flat_fn          <- make_step_fn(sdb$day, sdb$coverage_flat / 100)
cov_perfect_fn       <- function(t) rep(1.0, length(t))
cov_no_pep_fn        <- function(t) rep(0.0, length(t))

# DPC functions (take absolute time t, return DPC in days)
dpc_with_conflict_fn <- make_step_fn(sdb$day, sdb$dpc_conflict)
dpc_flat_fn          <- make_step_fn(sdb$day, rep(1, nrow(sdb)))
dpc_perfect_fn       <- function(t) rep(0, length(t))
dpc_no_pep_fn        <- function(t) rep(0, length(t))

# =============================================================================
# 2. Build efficacy functions
#
# obv_pep_efficacy receives individual DPC (days post-exposure) directly
# from fiber — NOT absolute time. All arms therefore share the same efficacy
# functions; differences between arms are captured solely by obv_pep_dpc,
# which maps absolute time t -> DPC.
# =============================================================================
curve_d50_dat <- readRDS(here("data-processed", "DPC_fixed_efficacy_varied_d50.rds"))

make_efficacy_fn_direct <- function(efficacy_col) {
  function(dpc) {
    approx(x = curve_d50_dat$dpc, y = curve_d50_dat[[efficacy_col]],
           xout = dpc, rule = 2)$y
  }
}

eff_mid <- make_efficacy_fn_direct("efficacy")
eff_lo  <- make_efficacy_fn_direct("eighty_efficacy_lo")
eff_hi  <- make_efficacy_fn_direct("eighty_efficacy_hi")

# =============================================================================
# 3. Arm definitions
#
# All arms run under DRC_conflict trajectory.
#   no_pep        : no PEP
#   with_conflict : both coverage and DPC conflict-impacted
#   cov_conflict  : coverage conflict-impacted, DPC flat (=1)
#   dpc_conflict  : DPC conflict-impacted, coverage flat
#   optimistic    : perfect coverage (100%) and DPC = 0
# =============================================================================
ARM_DEFS <- list(
  no_pep           = list(cov_fn = cov_no_pep_fn,        dpc_fn = dpc_no_pep_fn,        eff_fn = function(dpc) rep(0, length(dpc))),
  
  with_conflict_mid = list(cov_fn = cov_with_conflict_fn, dpc_fn = dpc_with_conflict_fn, eff_fn = eff_mid),
  with_conflict_lo  = list(cov_fn = cov_with_conflict_fn, dpc_fn = dpc_with_conflict_fn, eff_fn = eff_lo),
  with_conflict_hi  = list(cov_fn = cov_with_conflict_fn, dpc_fn = dpc_with_conflict_fn, eff_fn = eff_hi),
  
  cov_conflict_mid = list(cov_fn = cov_with_conflict_fn, dpc_fn = dpc_flat_fn,           eff_fn = eff_mid),
  cov_conflict_lo  = list(cov_fn = cov_with_conflict_fn, dpc_fn = dpc_flat_fn,           eff_fn = eff_lo),
  cov_conflict_hi  = list(cov_fn = cov_with_conflict_fn, dpc_fn = dpc_flat_fn,           eff_fn = eff_hi),
  
  dpc_conflict_mid = list(cov_fn = cov_flat_fn,          dpc_fn = dpc_with_conflict_fn,  eff_fn = eff_mid),
  dpc_conflict_lo  = list(cov_fn = cov_flat_fn,          dpc_fn = dpc_with_conflict_fn,  eff_fn = eff_lo),
  dpc_conflict_hi  = list(cov_fn = cov_flat_fn,          dpc_fn = dpc_with_conflict_fn,  eff_fn = eff_hi),
  
  optimistic_mid   = list(cov_fn = cov_perfect_fn,       dpc_fn = dpc_perfect_fn,        eff_fn = eff_mid),
  optimistic_lo    = list(cov_fn = cov_perfect_fn,       dpc_fn = dpc_perfect_fn,        eff_fn = eff_lo),
  optimistic_hi    = list(cov_fn = cov_perfect_fn,       dpc_fn = dpc_perfect_fn,        eff_fn = eff_hi)
)

ARM_NAMES <- names(ARM_DEFS)
message(sprintf("Total arms: %d", length(ARM_NAMES)))
print(ARM_NAMES)

# =============================================================================
# 4. Scenarios: all arms run under DRC_conflict trajectory
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

# =============================================================================
# 5. Pre-compute particle args per scenario
# =============================================================================

for (arm_name in ARM_NAMES) {
  dir.create(file.path(OUT_BASE, arm_name), recursive = TRUE, showWarnings = FALSE)
}

check_model_function_version()

message("Pre-computing particle args...")

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

# =============================================================================
# 6. Build job list
# =============================================================================
jobs <- lapply(seq_len(N_WORKERS), function(w) {
  p_from <- (w - 1L) * PARTICLES_PER_WORKER + 1L
  p_to   <-  w       * PARTICLES_PER_WORKER
  list(worker_id = w, p_from = p_from, p_to = p_to)
})

total_runs <- sum(sapply(names(SCENARIOS), function(sc_name) {
  length(SCENARIOS[[sc_name]]$arms) * N_PARTICLES * N_REPS
}))

message(sprintf(
  "%d workers | %d particles each | %d total runs",
  N_WORKERS, PARTICLES_PER_WORKER, total_runs
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
  make_step_fn            = make_step_fn,
  rescale_sdb_segment     = rescale_sdb_segment,
  make_efficacy_fn_direct = make_efficacy_fn_direct,
  cov_with_conflict_fn    = cov_with_conflict_fn,
  cov_flat_fn             = cov_flat_fn,
  cov_perfect_fn          = cov_perfect_fn,
  cov_no_pep_fn           = cov_no_pep_fn,
  dpc_with_conflict_fn    = dpc_with_conflict_fn,
  dpc_flat_fn             = dpc_flat_fn,
  dpc_perfect_fn          = dpc_perfect_fn,
  dpc_no_pep_fn           = dpc_no_pep_fn,
  eff_mid                 = eff_mid,
  eff_lo                  = eff_lo,
  eff_hi                  = eff_hi,
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