# =============================================================================
# check_prevented_completed.R
#
# Minimal local sanity check: run ONE DRC-like OBV/PEP simulation and confirm
# out$prevented_completed comes back populated (so the no-PEP baseline can be
# reconstructed as tdf + prevented_completed). Build path is identical to
# 01_analysis_figure2.R / 01_run_SI_sensitivity_DRC.R.
#
#   source(here::here("analyses", "04_SI_sensitivity_DRC", "check_prevented_completed.R"))
# =============================================================================
library(here)
library(fiber)
source(here("functions", "setup_model_parameters.R"))
source(here("functions", "calculate_model_approx_r0.R"))
source(here("functions", "abc_calibration_functions_common.R"))
source(here("functions", "abc_calibration_functions_decoupled.R"))
source(here("functions", "abc_posterior.R"))

env_decoupled <- new.env(parent = globalenv())
sys.source(here("functions", "abc_calibration_functions_decoupled.R"), envir = env_decoupled)
build_decoupled <- env_decoupled$build_abc_model_args_decoupled

PARTICLE <- 1L          # which downsampled particle to use
EFFICACY <- 0.80        # PEP efficacy
SEED     <- 20260801L

# ---- load posterior + build args for one particle (reference, no scaling) ----
rds <- here("outputs", "02_ABC_model_fits_Final",
            "fiber_ABC_SMC_Middle_DRC_ConflictSmoothed_PlusPlus_Decoupled_20260607_215621_check_NP5_NS4_NBREPS_30_NBSIMUL_590.rds")
param_names <- c("R0", "prop_funeral", "etu_efficacy", "ppe_efficacy", "hcw_risk_scalar")
res <- readRDS(rds)
posterior <- as.data.frame(res$param); colnames(posterior) <- param_names
posterior$weight <- res$weights
theta <- downsample_posterior(posterior, n_sets = 200, seed = 42L, param_names = param_names)

scenario_matrix <- read_scenario_matrix(here("data-processed", "final_six_scenario_values_original_approach.csv"))
mp  <- make_model_parameters("Middle_DRC_ConflictSmoothed_PlusPlus", scenario_matrix,
                             overrides = list(seeding_cases = 25L, check_final_size = 15000))
inv <- compute_R0_invariants(args = mp$args, n = 50000, seed = 42L)

args <- build_decoupled(
  R0 = theta$R0[PARTICLE], prop_funeral = theta$prop_funeral[PARTICLE],
  etu_efficacy = theta$etu_efficacy[PARTICLE], ppe_efficacy = theta$ppe_efficacy[PARTICLE],
  hcw_risk_scalar = theta$hcw_risk_scalar[PARTICLE],
  base = mp$base_args, tv = mp$tv_args, invariants = inv,
  general_hospital_quarantine_efficacy = DEFAULT_SCALAR_INPUTS$general_hospital_quarantine_efficacy,
  safe_funeral_efficacy = DEFAULT_SCALAR_INPUTS$safe_funeral_efficacy,
  hcw_base_prob = 0.25, seeding_cases = 25L)
args$check_final_size <- 15000

# ---- OBV / PEP settings (identical to the figure scripts) -------------------
args$obv_pep_enabled          <- TRUE
args$obv_pep_coverage         <- function(t) rep(1.0, length(t))
args$obv_pep_adherence        <- 1.0
args$obv_pep_dpc              <- 0
args$obv_pep_efficacy         <- EFFICACY
args$obv_pep_target_class     <- "HCW"
args$obv_pep_target_locations <- c("hospital", "community", "funeral")

# ---- run one sim, retrying the seed until it takes off ----------------------
seed <- SEED
repeat {
  args$seed <- seed
  out <- do.call(fiber::branching_process_main, args)
  tdf <- out$tdf[!is.na(out$tdf$time_infection_absolute), ]
  if (sum(!is.na(tdf$outcome) & tdf$outcome) >= 100) break
  seed <- seed + 1L
  if (seed - SEED > 50) { message("did not take off in 50 tries"); break }
}

# ---- report -----------------------------------------------------------------
pc  <- out$prevented_completed
hcwd <- function(x) sum(!is.na(x$class) & x$class == "HCW" & !is.na(x$outcome) & x$outcome)
n_obv  <- hcwd(tdf)
n_base <- n_obv + (if (!is.null(pc) && nrow(pc) > 0) hcwd(pc) else 0L)

cat("\n================ prevented_completed check ================\n")
cat(sprintf("particle %d: R0=%.3f  hcw_risk_scalar=%.3f  (seed=%d)\n",
            PARTICLE, theta$R0[PARTICLE], theta$hcw_risk_scalar[PARTICLE], seed))
cat(sprintf("is.null(prevented_completed)      : %s\n", is.null(pc)))
cat(sprintf("nrow(prevented_completed)         : %s\n",
            if (is.null(pc)) "NULL" else nrow(pc)))
cat(sprintf("obv_pep_num_treated$prevented     : %s\n",
            tryCatch(out$sim_info$obv_pep_num_treated$prevented, error = function(e) NA)))
cat(sprintf("prevented_completed columns       : %s\n",
            if (is.null(pc)) "-" else paste(names(pc), collapse = ", ")))
cat(sprintf("total deaths (tdf)                : %d\n", sum(!is.na(tdf$outcome) & tdf$outcome)))
cat(sprintf("HCW deaths WITH PEP (tdf)         : %d\n", n_obv))
cat(sprintf("HCW deaths no-PEP (tdf+prevented) : %d\n", n_base))
cat(sprintf("HCW deaths averted                : %d  (%.0f%%)\n",
            n_base - n_obv, if (n_base > 0) 100 * (n_base - n_obv) / n_base else NA))
if (is.null(pc) || nrow(pc) == 0)
  cat("\n>>> prevented_completed is EMPTY -- OBV prevention was not recorded. <<<\n")
cat("==========================================================\n")
