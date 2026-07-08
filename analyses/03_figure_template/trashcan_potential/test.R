library(fiber)
library(here)
# Single run to inspect tdf and prevented_completed
source(here("functions", "setup_model_parameters.R"))
source(here("functions", "calculate_model_approx_r0.R"))
source(here("functions", "abc_calibration_functions_common.R"))
source(here("functions", "abc_calibration_functions_decoupled.R"))
source(here("functions", "abc_posterior.R"))

env_decoupled <- new.env(parent = globalenv())
sys.source(here("functions", "abc_calibration_functions_decoupled.R"), envir = env_decoupled)

sc <- list(
  id           = "Worst_WestAfrica",
  rds          = here("outputs", "02_ABC_model_fits_Final",
                      "fiber_ABC_SMC_Worst_WestAfrica_Decoupled_20260608_162044_check_NP5_NS4_NBREPS_30_NBSIMUL_472.RDS"),
  scenario_csv = here("data-processed", "final_six_scenario_values_original_approach.csv"),
  param_names  = c("R0", "prop_funeral", "etu_efficacy", "ppe_efficacy", "hcw_risk_scalar"),
  check_final_size = 40000
)

res       <- readRDS(sc$rds)
posterior <- data.frame(weight = res$weights)
for (j in seq_along(sc$param_names)) {
  posterior[[sc$param_names[j]]] <- res$param[, j]
}
theta <- downsample_posterior(posterior, n_sets = 1L, seed = 42L,
                              param_names = sc$param_names)

scenario_matrix <- read_scenario_matrix(sc$scenario_csv)
mp  <- make_model_parameters(
  scenario_id     = sc$id,
  scenario_matrix = scenario_matrix,
  overrides       = list(seeding_cases = 25L, check_final_size = sc$check_final_size)
)
inv <- compute_R0_invariants(args = mp$args, n = 50000, seed = 42L)

args <- env_decoupled$build_abc_model_args_decoupled(
  R0              = theta$R0[1],
  prop_funeral    = theta$prop_funeral[1],
  etu_efficacy    = theta$etu_efficacy[1],
  ppe_efficacy    = theta$ppe_efficacy[1],
  hcw_risk_scalar = theta$hcw_risk_scalar[1],
  base            = mp$base_args,
  tv              = mp$tv_args,
  invariants      = inv,
  general_hospital_quarantine_efficacy = DEFAULT_SCALAR_INPUTS$general_hospital_quarantine_efficacy,
  safe_funeral_efficacy                = DEFAULT_SCALAR_INPUTS$safe_funeral_efficacy,
  hcw_base_prob   = 0.25,
  seeding_cases   = 25L
)

# OBV settings: efficacy = 0 to capture all potential preventions in tdf
args$obv_pep_enabled          <- TRUE
args$obv_pep_coverage         <- 1.0
args$obv_pep_adherence        <- 1.0
args$obv_pep_dpc              <- 0
args$obv_pep_efficacy         <- 0.8       # 0 so all HCW appear in tdf with obv_pep_received flag
args$obv_pep_target_class     <- "HCW"
args$obv_pep_target_locations <- c("hospital", "community", "funeral")
args$check_final_size         <- sc$check_final_size
args$seed                     <- 42L

out <- do.call(fiber::branching_process_main, args)

tdf       <- out$tdf[!is.na(out$tdf$time_infection_absolute), ]
prevented <- out$prevented_completed

# Inspect
cat("tdf rows:", nrow(tdf), "\n")
cat("tdf class distribution:\n"); print(table(tdf$class))
cat("obv_pep_received distribution:\n"); print(table(tdf$obv_pep_received))
cat("\nprevented_completed rows:", nrow(prevented), "\n")
cat("prevented class distribution:\n"); print(table(prevented$class))