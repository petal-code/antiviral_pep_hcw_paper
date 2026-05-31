# DRC_run_abc_calibration.R  (NPI-efficacy approach)
# =============================================================================
# ABC-SMC calibration of the revamped fiber branching-process model against a
# DRC Ebola outbreak (Middle_DRC_ConflictSmoothed scenario).
#
# Fitted parameters (3):
#   R0           : baseline reproduction number for a typical genPop seeding
#                  case at t = 0 (under this particle's efficacies).
#   prop_funeral : share of R0 attributable to funeral transmission at t = 0.
#   npi_scaler   : SINGLE shared scaler in [-1, 1] that positions the two fitted
#                  conditional efficacies (ppe_efficacy, etu_efficacy) on their
#                  [min, max] intervals (Approach B).
#                    s = -1 -> min,  s = 0 -> central,  s = +1 -> max.
#
# FIXED (not fitted): general_hospital_quarantine_efficacy, safe_funeral_efficacy.
#   prob_hcw_cond_*_hospital fixed at 0.25 (no HCW-risk scalar).
#
# >>> PLACEHOLDERS <<<  NPI_SPEC bounds and FIXED_EFFICACIES are PLACEHOLDER
#     values. UPDATE WITH REAL NUMBERS before any production run.
# =============================================================================


# -----------------------------------------------------------------------------
# 1. CONFIGURATION
# -----------------------------------------------------------------------------
ANALYSIS_DIR <- switch(
  Sys.info()[["user"]],
  "cwhittaker" = "C:/Users/cwhittaker/Documents/Research Projects/obv_hcw_paper/analyses/02_ABC_model_fits_NPI_Eff",
  "PETAL_WS_1" = "C:/Users/PETAL_WS_1/Documents/obv_hcw_paper/analyses/02_ABC_model_fits_NPI_Eff",
  "PETAL_WS_2" = "C:/Users/PETAL_WS_2/Documents/obv_hcw_paper/analyses/02_ABC_model_fits_NPI_Eff",
  getwd()
)

FUNCTIONS_DIR  <- normalizePath(file.path(ANALYSIS_DIR, "..", "..", "functions"),
                                mustWork = FALSE)
SETUP_PATH     <- file.path(FUNCTIONS_DIR, "setup_model_parameters_new_approach.R")
FUNCTIONS_PATH <- file.path(FUNCTIONS_DIR, "abc_calibration_functions_new_approach.R")
R0_PATH        <- file.path(FUNCTIONS_DIR, "calculate_model_approx_r0_new_approach.R")
SCENARIO_CSV   <- file.path(ANALYSIS_DIR, "final_four_scenario_values.csv")
SCENARIO_ID    <- "Middle_DRC_ConflictSmoothed"

# >>> PLACEHOLDER: fixed conditional efficacies (NOT fitted). UPDATE. <<<
FIXED_EFFICACIES <- list(
  general_hospital_quarantine_efficacy = 0.30,   # PLACEHOLDER
  safe_funeral_efficacy                = 0.95    # PLACEHOLDER
)

# >>> PLACEHOLDER: bounds for the two fitted efficacies (Approach B). UPDATE. <<<
NPI_SPEC <- list(
  ppe_efficacy = list(min = 0.20, max = 0.90),   # PLACEHOLDER
  etu_efficacy = list(min = 0.50, max = 0.95)    # PLACEHOLDER
)

MODEL_OVERRIDES <- c(
  list(check_final_size = 10000),
  FIXED_EFFICACIES
)

ABC_OUTPUT_BASE  <- ANALYSIS_DIR
ABC_OUTPUT_LABEL <- "NPIeff"

FINAL_OUTPUTS_DIR <- normalizePath(
  file.path(ANALYSIS_DIR, "..", "..", "outputs"),
  mustWork = FALSE
)
if (!dir.exists(FINAL_OUTPUTS_DIR)) {
  dir.create(FINAL_OUTPUTS_DIR, recursive = TRUE, showWarnings = FALSE)
}

ABC_CONFIG <- list(
  takeoff_death_threshold = 100,
  n_reps                  = 150,
  seeding_cases           = 25,
  setup_R0_n              = 100000L,
  setup_R0_seed           = 42L,
  npi_spec                = NPI_SPEC,
  general_hospital_quarantine_efficacy =
    FIXED_EFFICACIES$general_hospital_quarantine_efficacy,
  safe_funeral_efficacy   = FIXED_EFFICACIES$safe_funeral_efficacy
)

ABC_SETTINGS <- list(
  method              = "Delmoral",
  nb_simul            = 690,
  alpha               = 0.5,
  tolerance_target    = 0.2,
  M                   = 1,
  use_seed            = TRUE,
  verbose             = TRUE
)

N_CLUSTER <- if (grepl("PETAL", Sys.info()[["user"]], ignore.case = TRUE)) {
  min(120, parallel::detectCores() - 10)
} else {
  min(10, parallel::detectCores() - 4)
}


# -----------------------------------------------------------------------------
# 2. LIBRARIES + SOURCES
# -----------------------------------------------------------------------------
library(EasyABC)
library(future)
library(future.apply)
library(parallel)
library(progressr)
library(fiber)
handlers("progress")

source(SETUP_PATH)
source(FUNCTIONS_PATH)
source(R0_PATH)

check_model_function_version()


# -----------------------------------------------------------------------------
# 3. BUILD BASE + TIME-VARYING ARGS; COMPUTE EFFICACY-INDEPENDENT R0 INVARIANTS
# -----------------------------------------------------------------------------
scenario_matrix <- read_scenario_matrix(SCENARIO_CSV)

mp <- make_model_parameters(
  scenario_id     = SCENARIO_ID,
  scenario_matrix = scenario_matrix,
  overrides       = MODEL_OVERRIDES
)
base_args     <- mp$base_args
tv_args_model <- mp$tv_args

R0_invariants <- compute_R0_invariants(
  args = mp$args,
  n    = ABC_CONFIG$setup_R0_n,
  seed = ABC_CONFIG$setup_R0_seed
)

central_effs <- npi_efficacy_from_scaler(0, NPI_SPEC)
cat("Central-efficacy check (s = 0):\n")
cat(sprintf("  ppe_efficacy = %.3f, etu_efficacy = %.3f\n",
            central_effs$ppe_efficacy, central_effs$etu_efficacy))
cat(sprintf("  D = %.4f, F = %.4f, Q_g = %.4f\n",
            D_from_invariants(R0_invariants, central_effs$etu_efficacy,
                              FIXED_EFFICACIES$general_hospital_quarantine_efficacy),
            F_from_invariants(R0_invariants, FIXED_EFFICACIES$safe_funeral_efficacy),
            R0_invariants$Q_g))


# -----------------------------------------------------------------------------
# 4. OBSERVED TARGETS AND PRIORS
# -----------------------------------------------------------------------------
observed_summaries <- c(
  takeoff      = 1.0,
  n_deaths     = 2299,
  n_hcw_deaths = 79,     # https://afenet-journal.org/10-37432-jieph-d-25-00072/
  duration     = 450     # ~ Aug 2018 - Nov 2019 main phase
)

priors <- list(
  c("unif", 1.25, 1.65),   # R0
  c("unif", 0.10, 0.40),   # prop_funeral
  c("unif", -1.0, 1.0)     # npi_scaler (s); s=0 -> central efficacies
)


# -----------------------------------------------------------------------------
# 5. PRIOR PREDICTIVE CHECK (optional, slow)
# -----------------------------------------------------------------------------
# set.seed(1)
# pp <- prior_predictive_check(
#   n_draws       = 50,
#   prior_list    = priors,
#   base          = base_args,
#   tv            = tv_args_model,
#   invariants    = R0_invariants,
#   npi_spec      = NPI_SPEC,
#   general_hospital_quarantine_efficacy =
#     FIXED_EFFICACIES$general_hospital_quarantine_efficacy,
#   safe_funeral_efficacy = FIXED_EFFICACIES$safe_funeral_efficacy,
#   parallel      = FALSE,
#   n_replicates  = 5,
#   seeding_cases = ABC_CONFIG$seeding_cases,
#   takeoff_death_threshold = ABC_CONFIG$takeoff_death_threshold
# )
# print(pp); summary(pp)


# -----------------------------------------------------------------------------
# 6. PER-RUN OUTPUT DIRECTORY + WORKER CONFIG
# -----------------------------------------------------------------------------
ABC_OUTPUT_DIR <- make_abc_output_dir(
  base_dir    = ABC_OUTPUT_BASE,
  scenario_id = SCENARIO_ID,
  label       = ABC_OUTPUT_LABEL
)
message("ABC outputs will be written to: ", ABC_OUTPUT_DIR)

save_abc_config(list(
  setup_path      = SETUP_PATH,
  functions_path  = FUNCTIONS_PATH,
  r0_path         = R0_PATH,
  scenario_csv    = SCENARIO_CSV,
  scenario_id     = SCENARIO_ID,
  abc_config      = ABC_CONFIG,
  model_overrides = MODEL_OVERRIDES
))


# -----------------------------------------------------------------------------
# 7. RUN ABC_SEQUENTIAL (Del Moral et al. 2012 adaptive SMC)
# -----------------------------------------------------------------------------
start_time <- Sys.time()
result <- with_abc_output_dir(
  ABC_OUTPUT_DIR,
  ABC_sequential(
    method              = ABC_SETTINGS$method,
    model               = fiber_abc_model_parallel,
    prior               = priors,
    nb_simul            = ABC_SETTINGS$nb_simul,
    summary_stat_target = observed_summaries,
    alpha               = ABC_SETTINGS$alpha,
    tolerance_target    = ABC_SETTINGS$tolerance_target,
    M                   = ABC_SETTINGS$M,
    use_seed            = ABC_SETTINGS$use_seed,
    verbose             = ABC_SETTINGS$verbose,
    n_cluster           = N_CLUSTER
  )
)
end_time <- Sys.time()
print(end_time - start_time)

result_date <- format(start_time, "%Y-%m-%d")
result_filename <- paste0(
  "fiber_ABC_SMC_", SCENARIO_ID,
  if (nzchar(ABC_OUTPUT_LABEL)) paste0("_", ABC_OUTPUT_LABEL) else "",
  "_", result_date, ".rds"
)
saveRDS(result, file = file.path(ABC_OUTPUT_DIR, result_filename))
saveRDS(result, file = file.path(FINAL_OUTPUTS_DIR, result_filename))


# -----------------------------------------------------------------------------
# 8. POSTERIOR INSPECTION
# -----------------------------------------------------------------------------
posterior <- as.data.frame(result$param)
colnames(posterior) <- c("R0", "prop_funeral", "npi_scaler")

print(apply(posterior, 2, quantile, probs = c(0.025, 0.5, 0.975)))

implied_ppe <- vapply(posterior$npi_scaler,
                      function(s) npi_efficacy_from_scaler(s, NPI_SPEC)$ppe_efficacy, numeric(1))
implied_etu <- vapply(posterior$npi_scaler,
                      function(s) npi_efficacy_from_scaler(s, NPI_SPEC)$etu_efficacy, numeric(1))
cat("\nImplied ppe_efficacy: ", paste(round(quantile(implied_ppe, c(0.025,0.5,0.975)), 3), collapse = " / "), "\n")
cat("Implied etu_efficacy: ", paste(round(quantile(implied_etu, c(0.025,0.5,0.975)), 3), collapse = " / "), "\n")

par(mfrow = c(1, 3))
for (j in seq_len(ncol(posterior))) {
  hist(posterior[, j], breaks = 10, main = colnames(posterior)[j],
       xlab = colnames(posterior)[j])
  abline(v = quantile(posterior[, j], c(0.025, 0.5, 0.975)),
         lty = c(2, 1, 2), col = "red")
}
par(mfrow = c(1, 1))


# -----------------------------------------------------------------------------
# 9. PROGRESS / RECONSTRUCTION FROM DISK
# -----------------------------------------------------------------------------
# abc_progress(ABC_OUTPUT_DIR, tolerance_target = ABC_SETTINGS$tolerance_target)
# print(abc_compare_steps(ABC_OUTPUT_DIR))
# result <- reconstruct_abc_result(ABC_OUTPUT_DIR)


# -----------------------------------------------------------------------------
# 10. POSTERIOR PREDICTIVE CHECKS
# -----------------------------------------------------------------------------
sim_stats <- as.data.frame(result$stats)
colnames(sim_stats) <- c("takeoff", "n_deaths", "n_hcw_deaths", "duration")

set.seed(1)
idx <- sample(seq_len(nrow(sim_stats)), size = 10000, replace = TRUE,
              prob = result$weights)
sim_stats_post <- sim_stats[idx, ]

par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (s in names(observed_summaries)) {
  x  <- sim_stats_post[[s]]
  qs <- quantile(x, probs = c(0.025, 0.5, 0.975))
  hist(x, breaks = 10, main = paste0("Posterior predictive: ", s), xlab = s,
       col = adjustcolor("steelblue", alpha = 0.6), border = "white")
  abline(v = qs, col = "darkblue", lty = c(2, 1, 2), lwd = c(1, 2, 1))
  abline(v = observed_summaries[s], col = "red", lwd = 2.5)
}
par(mfrow = c(1, 1))
