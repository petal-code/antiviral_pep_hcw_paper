# =============================================================================
# 01_analysis_figure3new_tv_param_sensitivity.R
#
# Sensitivity analysis: instead of perturbing the conflict's impact on drug
# DISTRIBUTION (coverage/DPC -- see the conflict_dpc_sensitivity scripts),
# here we perturb the underlying EPIDEMIOLOGICAL time-varying parameters
# themselves, holding drug distribution (coverage/DPC) at baseline.
#
# Five time-varying parameters come from the scenario CSV
# (final_six_scenario_values_original_approach.csv), for the
# "Middle_DRC_ConflictSmoothed_PlusPlus" scenario rows only:
#   prob_hosp                : higher = better  (earlier hospitalisation
#                               suppresses community transmission and unsafe
#                               funerals)
#   delay_hosp                : lower  = better  (time to hospitalisation)
#   prob_unsafe_funeral_comm  : lower  = better  (community unsafe funerals)
#   prob_unsafe_funeral_hosp  : lower  = better  (hospital unsafe funerals)
#   prop_etu                  : higher = better  (ETU treatment proportion)
#
# Four conditions, moving ALL FIVE parameters together in the same direction:
#   goodgood : strongest good-direction perturbation (50%)
#   good     : mild good-direction perturbation      (25%)
#   bad      : mild bad-direction perturbation        (25%)
#   badbad   : strongest bad-direction perturbation   (50%)
#
# Multiplier is 1.25/1.5 for an increase, 0.75/0.5 for a decrease, chosen
# per-parameter so it points in the "good" or "bad" direction as above.
# Probability-valued parameters (prob_hosp, prob_unsafe_funeral_comm,
# prob_unsafe_funeral_hosp, prop_etu) are capped at 1 after scaling;
# delay_hosp (a time, not a probability) is not capped.
#
# Drug distribution (coverage/DPC) is held at baseline intensity throughout
# (coverage max 80%, DPC max +4 days -- see
# 01_analysis_figure3new_conflict_dpc.R). Only "with_conflict" (both
# impacted) arms are run, at all three efficacy levels (mid/lo/hi); no
# separate no_pep arm is simulated -- averted outcomes are read directly off
# each run's own prevented_completed data, as elsewhere in this pipeline.
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
SEED_BASE     <- 20260703L   # kept identical across TV conditions -> matched-seed comparisons
HCW_BASE_PROB <- 0.25

TAKEOFF_DEATH_THRESHOLD <- 100L
MAX_RETRIES             <- 50L

stopifnot(N_PARTICLES %% N_WORKERS == 0L)
PARTICLES_PER_WORKER <- N_PARTICLES %/% N_WORKERS

SCENARIO_ID  <- "Middle_DRC_ConflictSmoothed_PlusPlus"
SCENARIO_CSV <- here("data-processed", "final_six_scenario_values_original_approach.csv")

# =============================================================================
# 1. Baseline coverage(t) and dpc(t) functions (fixed across all TV
#    conditions -- unchanged from 01_analysis_figure3new_conflict_dpc.R)
# =============================================================================
sdb <- readRDS(here("data-processed", "SDB_communityDeath_blended.rds"))

rescale_sdb_segment <- function(sdb_ref, day_out, orig_from, orig_to) {
  t        <- (day_out - day_out[1]) / (day_out[length(day_out)] - day_out[1])
  day_orig <- orig_from + t * (orig_to - orig_from)
  approx(sdb_ref$day, sdb_ref$value, xout = day_orig, rule = 2)$y
}

idx_150_325 <- sdb$day >= 150 & sdb$day <= 325
idx_325_400 <- sdb$day >  325 & sdb$day <= 400

sdb_tweaked <- sdb$value
sdb_tweaked[idx_150_325] <- rescale_sdb_segment(sdb, sdb$day[idx_150_325], 150, 200)
sdb_tweaked[idx_325_400] <- rescale_sdb_segment(sdb, sdb$day[idx_325_400], 200, 350)
sdb$value_tweaked <- sdb_tweaked

sdb$coverage_conflict <- sdb$value_tweaked * 80 / max(sdb$value_tweaked)
sdb$dpc_conflict      <- 1 + 4 * (1 - (sdb$value_tweaked / max(sdb$value_tweaked)))

sub      <- sdb[sdb$day < 200, ]
peak_row <- sub[which.max(sub$coverage_conflict), ]
peak_day <- peak_row$day
sdb$dpc_conflict[sdb$day <= peak_day] <- 1

message(sprintf("peak_day = %d, peak_coverage = %.2f", peak_day, peak_row$coverage_conflict))

make_step_fn <- function(days, values) {
  force(days); force(values)
  function(t) approx(x = days, y = values, xout = t, rule = 2)$y
}

cov_with_conflict_fn <- make_step_fn(sdb$day, sdb$coverage_conflict / 100)
dpc_with_conflict_fn <- make_step_fn(sdb$day, sdb$dpc_conflict)

# =============================================================================
# 2. Efficacy functions (unaffected by TV condition -- identical to baseline)
# =============================================================================
curve_d50_dat <- readRDS(here("data-processed", "DPC_fixed_efficacy_varied_d50.rds"))

make_efficacy_fn_direct <- function(efficacy_col) {
  force(efficacy_col)
  function(dpc) {
    approx(x = curve_d50_dat$dpc, y = curve_d50_dat[[efficacy_col]],
           xout = dpc, rule = 2)$y
  }
}

eff_mid <- make_efficacy_fn_direct("efficacy")
eff_lo  <- make_efficacy_fn_direct("eighty_efficacy_lo")
eff_hi  <- make_efficacy_fn_direct("eighty_efficacy_hi")

ARM_DEFS <- list(
  with_conflict_mid = list(cov_fn = cov_with_conflict_fn, dpc_fn = dpc_with_conflict_fn, eff_fn = eff_mid),
  with_conflict_lo  = list(cov_fn = cov_with_conflict_fn, dpc_fn = dpc_with_conflict_fn, eff_fn = eff_lo),
  with_conflict_hi  = list(cov_fn = cov_with_conflict_fn, dpc_fn = dpc_with_conflict_fn, eff_fn = eff_hi)
)
ARM_NAMES <- names(ARM_DEFS)

# =============================================================================
# 3. TV-parameter perturbation logic
# =============================================================================
DIRECTION_MAP <- c(
  prob_hosp                = "good_high",
  delay_hosp                = "good_low",
  prob_unsafe_funeral_comm = "good_low",
  prob_unsafe_funeral_hosp = "good_low",
  prop_etu                  = "good_high"
)

CAPPED_PARAMS <- c("prob_hosp", "prob_unsafe_funeral_comm", "prob_unsafe_funeral_hosp", "prop_etu")

CONDITIONS <- list(
  goodgood = list(type = "good", tier = 2),
  good     = list(type = "good", tier = 1),
  bad      = list(type = "bad",  tier = 1),
  badbad   = list(type = "bad",  tier = 2)
)

get_multiplier <- function(direction, type, tier) {
  increase <- if (tier == 2) 1.5 else 1.25
  decrease <- if (tier == 2) 0.5 else 0.75
  if (direction == "good_high") {
    if (type == "good") increase else decrease
  } else {
    if (type == "good") decrease else increase
  }
}

apply_tv_condition <- function(scenario_matrix, scenario_id, type, tier) {
  sm   <- scenario_matrix
  rows <- sm$scenario == scenario_id
  for (param in names(DIRECTION_MAP)) {
    mult    <- get_multiplier(DIRECTION_MAP[[param]], type, tier)
    new_val <- sm[[param]][rows] * mult
    if (param %in% CAPPED_PARAMS) new_val <- pmin(new_val, 1)
    sm[[param]][rows] <- new_val
  }
  sm
}

base_scenario_matrix <- read_scenario_matrix(SCENARIO_CSV)

# =============================================================================
# 4. Posterior resample (identical across TV conditions -- matched-seed)
# =============================================================================
check_model_function_version()

RDS_PATH    <- here("outputs", "02_ABC_model_fits_Final",
                    "fiber_ABC_SMC_Middle_DRC_ConflictSmoothed_PlusPlus_Decoupled_20260607_215621_check_NP5_NS4_NBREPS_30_NBSIMUL_590.RDS")
PARAM_NAMES <- c("R0", "prop_funeral", "etu_efficacy", "ppe_efficacy", "hcw_risk_scalar")
CHECK_FINAL_SIZE <- 10000

res       <- readRDS(RDS_PATH)
posterior <- data.frame(weight = res$weights)
for (j in seq_along(PARAM_NAMES)) {
  posterior[[PARAM_NAMES[j]]] <- res$param[, j]
}
theta <- downsample_posterior(posterior, n_sets = N_PARTICLES,
                              seed = RESAMPLE_SEED, param_names = PARAM_NAMES)

build_decoupled <- env_decoupled$build_abc_model_args_decoupled

# =============================================================================
# 5. Build particle args for one TV condition (expensive -- mp/inv depend on
#    the perturbed scenario_matrix, so this must be redone per condition)
# =============================================================================
build_particle_args_for_condition <- function(cond_name, cond) {
  message(sprintf("  Building particle args for condition '%s' (type=%s, tier=%d)...",
                  cond_name, cond$type, cond$tier))
  
  scenario_matrix_cond <- apply_tv_condition(base_scenario_matrix, SCENARIO_ID, cond$type, cond$tier)
  
  mp <- make_model_parameters(
    scenario_id     = SCENARIO_ID,
    scenario_matrix = scenario_matrix_cond,
    overrides       = list(seeding_cases    = SEEDING_CASES,
                           check_final_size = CHECK_FINAL_SIZE)
  )
  inv <- compute_R0_invariants(args = mp$args, n = 50000, seed = 42L)
  
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
    args$check_final_size <- CHECK_FINAL_SIZE
    args
  })
  
  base_particle_args
}

# =============================================================================
# 6. Parallel run for one TV condition
# =============================================================================
run_condition <- function(OUT_BASE, base_particle_args) {
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
      for (arm_name in ARM_NAMES) {
        arm_idx <- which(ARM_NAMES == arm_name)
        arm_def <- ARM_DEFS[[arm_name]]
        
        for (r in seq_len(N_REPS)) {
          fname    <- sprintf("DRC_conflict_p%03d_r%02d.rds", p, r)
          out_path <- file.path(OUT_BASE, arm_name, fname)
          if (file.exists(out_path)) next
          
          seed <- SEED_BASE +
            (arm_idx - 1L) * N_PARTICLES * N_REPS +
            (p       - 1L) * N_REPS +
            (r       - 1L)
          
          args <- base_particle_args[[p]]
          args$obv_pep_enabled          <- TRUE
          args$obv_pep_coverage         <- arm_def$cov_fn
          args$obv_pep_adherence        <- 1.0
          args$obv_pep_dpc              <- arm_def$dpc_fn
          args$obv_pep_efficacy         <- arm_def$eff_fn
          args$obv_pep_target_class     <- "HCW"
          args$obv_pep_target_locations <- c("hospital", "community", "funeral")
          
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
                "  [w%02d | %s | p%03d | r%02d] WARNING: max retries (n_deaths=%d)",
                job$worker_id, arm_name, p, r, n_deaths
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
    OUT_BASE                = OUT_BASE,
    ARM_NAMES               = ARM_NAMES,
    ARM_DEFS                = ARM_DEFS,
    base_particle_args      = base_particle_args,
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
# 7. Sweep all four TV conditions
# =============================================================================
for (cond_name in names(CONDITIONS)) {
  cond <- CONDITIONS[[cond_name]]
  message(sprintf("\n===== TV condition: %s (type=%s, tier=%d) =====",
                  cond_name, cond$type, cond$tier))
  
  base_particle_args <- build_particle_args_for_condition(cond_name, cond)
  OUT_BASE <- here("outputs", "simulation", sprintf("tv_sens_%s", cond_name))
  
  run_condition(OUT_BASE, base_particle_args)
}

message("\nTV-parameter sensitivity sweep complete.")