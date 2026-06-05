# westAfrica_run_abc_calibration_combined.R  (COMBINED approach: NPI-efficacy + HCW-risk)
# =============================================================================
# ABC-SMC calibration of the revamped fiber branching-process model against the
# 2014-16 West Africa Ebola outbreak (Worst_WestAfrica), using the COMBINED
# parameterisation (functions/abc_calibration_functions_combined.R) that unions
# the NPI-efficacy and HCW-risk schemes. This is the West Africa twin of
# DRC_run_abc_calibration_combined.R.
#
# The COMBINED scheme adds, on top of the pure NPI-efficacy fit:
#   (1) `hcw_risk_scalar` as a FOURTH fittable parameter (it scales the symmetric
#       prob_hcw_cond_*_hospital base), so HCW deaths get a lever that is near-
#       orthogonal to overall size; and
#   (2) the curve-shape summaries `time_to_peak` / `peak_height`, which (when
#       supplied) break the R0 <-> etu_efficacy size-confound that pure n_deaths
#       leaves.
#
# >>> PEAK TARGETS LEFT BLANK FOR YOU TO FILL IN <<<
#     The two curve-shape summaries, `time_to_peak` and `peak_height`, ARE part of
#     the fitted set (SUMMARY_STATS), matching DRC_run_abc_calibration_combined.R --
#     but their OBSERVED values in OBSERVED_NAMED are deliberately left blank
#     (NA_real_) as a reminder to fill them in. Read them off the West Africa weekly
#     DEATH surveillance series:
#         peak_height  = peak weekly DEATHS
#         time_to_peak = days from first death to that peak week
#     A guard just before the ABC run (section 7) STOPS with a clear message while
#     either target is still NA, so an unfilled value can't silently propagate into
#     ABC_sequential(). (Together these two break the R0 <-> etu_efficacy size-
#     confound that pure n_deaths leaves.)
#
# >>> FULLY CONFIGURABLE <<<  Choose WHICH parameters to fit (FIT_PARAMS) and
#     WHICH summaries to match (SUMMARY_STATS) by editing the two blocks in
#     section 1. Non-fitted parameters are held at FIXED_PARAMS (falling back to
#     COMBINED_PARAM_DEFAULTS). priors/observed are keyed BY NAME, so you only
#     edit one place and prepare_combined_run() keeps the ABC inputs in order.
#
#     Reductions:
#       * fit (R0, prop_funeral, npi_scaler) + fix hcw_risk_scalar = 1.0
#         -> the pure NPI-efficacy fit (== westAfrica_run_abc_calibration.R).
#       * fit (R0, prop_funeral, hcw_risk_scalar) + fix npi_scaler
#         -> an HCW-risk-style fit parameterised through NPI_SPEC.
#
# >>> PLACEHOLDERS <<<  NPI_SPEC bounds, the hcw_risk_scalar prior, and the
#     still-blank (NA) observed peak targets are PLACEHOLDER values. UPDATE WITH
#     REAL NUMBERS (peak weekly deaths from the West Africa weekly surveillance
#     series; literature bounds for the efficacies/HCW exposure) before any
#     production run.
# =============================================================================


# -----------------------------------------------------------------------------
# 1. CONFIGURATION
# -----------------------------------------------------------------------------
# Paths resolved from the repo root with here::here() (locates obv_hcw_paper.Rproj),
# so the script runs the same regardless of working directory / machine.
ANALYSIS_DIR   <- here::here("analyses", "02_ABC_model_fits_NPI_Eff")
FUNCTIONS_DIR  <- here::here("functions")
SETUP_PATH     <- file.path(FUNCTIONS_DIR, "setup_model_parameters.R")
COMMON_PATH    <- file.path(FUNCTIONS_DIR, "abc_calibration_functions_common.R")
NPI_PATH       <- file.path(FUNCTIONS_DIR, "abc_calibration_functions_npi.R")       # reused helpers
FUNCTIONS_PATH <- file.path(FUNCTIONS_DIR, "abc_calibration_functions_combined.R")  # this scheme
R0_PATH        <- file.path(FUNCTIONS_DIR, "calculate_model_approx_r0.R")
SCENARIO_CSV   <- here::here("data-processed", "final_six_scenario_values_original_approach.csv")
SCENARIO_ID    <- "Worst_WestAfrica"

# ---- WHICH PARAMETERS TO FIT (edit me) --------------------------------------
# Any subset of c("R0", "prop_funeral", "npi_scaler", "hcw_risk_scalar"); this is
# also the order ABC stores them. Non-fitted ones are held at FIXED_PARAMS.
FIT_PARAMS <- c("R0", "prop_funeral", "npi_scaler", "hcw_risk_scalar")

# Priors for the fitted parameters, keyed BY NAME (only entries for FIT_PARAMS
# are used; order here does not matter -- prepare_combined_run() reorders).
# R0 / prop_funeral inherited from the West Africa NPI-efficacy fit
# (westAfrica_run_abc_calibration.R).
PRIORS_NAMED <- list(
  R0              = c("unif", 1.15, 1.65),
  prop_funeral    = c("unif", 0.10, 0.40),
  npi_scaler      = c("unif", -1.0, 1.0),   
  hcw_risk_scalar = c("unif",  1.0, 1.8)
)

# Fixed values for any parameter NOT in FIT_PARAMS (missing ones fall back to
# COMBINED_PARAM_DEFAULTS). hcw_risk_scalar = 1.0 -> prob_hcw at the 0.25 base;
# npi_scaler = 0 -> central efficacies.
FIXED_PARAMS <- list(
  R0              = 1.40,
  prop_funeral    = 0.35,
  npi_scaler      = 0.0,
  hcw_risk_scalar = 1.0
)

# ---- WHICH SUMMARIES TO FIT (edit me) ---------------------------------------
# Any subset of c("takeoff", "n_cases", "n_deaths", "n_hcw_deaths", "duration",
# "time_to_peak", "peak_height"). Mirrors the DRC combined fit: the two curve-shape
# summaries are INCLUDED, but their observed targets in OBSERVED_NAMED are left
# blank (NA) for you to fill in (the section-7 guard stops the run until you do).
SUMMARY_STATS <- c("n_deaths", "n_hcw_deaths", "duration", "time_to_peak", "peak_height")

# Observed targets, keyed BY NAME (only entries for SUMMARY_STATS are used).
OBSERVED_NAMED <- c(
  # takeoff      = 1.0,
  n_deaths     = 11325,
  n_hcw_deaths = 513,
  duration     = 365,
  time_to_peak = 164,
  peak_height  = 599
)
# time to peak and peak height from: https://en.wikipedia.org/wiki/West_African_Ebola_virus_epidemic_timeline_of_reported_cases_and_deaths
# 12th October deaths peak, starting from 1st May "start" (approx) yields 164 days in (again approx, as that's the end of the reporting period for those deaths)

# HCW per-contact exposure base that hcw_risk_scalar multiplies (prob =
# min(base * scalar, 1)). 0.25 reproduces DEFAULT_SCALAR_INPUTS.
HCW_BASE_PROB <- 0.25

# Weekly DEATH-incidence binning for the curve summaries (time_to_peak/peak_height).
# Inert while the peak summaries are not fitted, but kept so the run is self-
# describing and the settings are ready when you switch them on.
PEAK_BIN_WIDTH   <- 7L
PEAK_TIME_ORIGIN <- "first_death"   # "first_death" or "outbreak_start"

# >>> PLACEHOLDER: bounds for the two NPI efficacies (positioned by npi_scaler). <<<
NPI_SPEC <- list(
  ppe_efficacy = list(min = 0.30, max = 0.90),   # PLACEHOLDER
  etu_efficacy = list(min = 0.60, max = 0.95)    # PLACEHOLDER
)

# West Africa was a far larger outbreak than the DRC; keep the larger final-size
# guard from the West Africa NPI-efficacy fit.
MODEL_OVERRIDES <- list(check_final_size = 35000)

ABC_OUTPUT_BASE  <- ANALYSIS_DIR
ABC_OUTPUT_LABEL <- "Combined"   # distinct from "NPIeff" / "NPIeffPeak"

FINAL_OUTPUTS_DIR <- here::here("outputs", "02_ABC_model_fits_NPI_Eff")
if (!dir.exists(FINAL_OUTPUTS_DIR)) {
  dir.create(FINAL_OUTPUTS_DIR, recursive = TRUE, showWarnings = FALSE)
}

ABC_SETTINGS <- list(
  method              = "Delmoral",
  nb_simul            = 460,    # 4 params -> consider 600-800; quick test ~200
  alpha               = 0.5,
  tolerance_target    = 0.5,
  M                   = 1,
  use_seed            = TRUE,
  verbose             = TRUE
)

# Replicate count / take-off threshold etc. travel to each worker via ABC_CONFIG
# (assembled in section 3, after the function files are sourced).
N_REPS                  <- 80L   # per-particle stochastic replicates (quick test ~30)
SEEDING_CASES           <- 25L
TAKEOFF_DEATH_THRESHOLD <- 100L
SETUP_R0_N              <- 100000L
SETUP_R0_SEED           <- 42L

N_CLUSTER <- if (parallel::detectCores() > 120) {
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
source(COMMON_PATH)
source(NPI_PATH)        # npi_efficacy_from_scaler(), DEFAULT_NPI_SPEC, metadata helpers
source(FUNCTIONS_PATH)  # the COMBINED scheme
source(R0_PATH)

check_model_function_version()


# -----------------------------------------------------------------------------
# 3. VALIDATE CHOICES + ASSEMBLE ORDER-SAFE ABC INPUTS
# -----------------------------------------------------------------------------
# prepare_combined_run() checks FIT_PARAMS / SUMMARY_STATS against the canonical
# sets, confirms every fitted parameter has a prior and every fitted summary an
# observed target, and returns the positional `priors` + `observed_summaries`
# (in fit/summary order) plus the complete `fixed_values` list. Assembled here
# (not in section 1) because it needs the just-sourced function files.
prep <- prepare_combined_run(
  fit_params     = FIT_PARAMS,
  priors_named   = PRIORS_NAMED,
  fixed_params   = FIXED_PARAMS,
  summary_stats  = SUMMARY_STATS,
  observed_named = OBSERVED_NAMED
)
priors             <- prep$priors
observed_summaries <- prep$observed_summaries
fixed_values       <- prep$fixed_values
PARAM_NAMES        <- prep$fit_params   # alias used by the disk-inspection helpers

# Everything the workers need travels via the saved config (section 6).
ABC_CONFIG <- list(
  takeoff_death_threshold = TAKEOFF_DEATH_THRESHOLD,
  n_reps                  = N_REPS,
  seeding_cases           = SEEDING_CASES,
  setup_R0_n              = SETUP_R0_N,
  setup_R0_seed           = SETUP_R0_SEED,
  npi_spec                = NPI_SPEC,
  hcw_base_prob           = HCW_BASE_PROB,
  fit_params              = prep$fit_params,
  fixed_values            = fixed_values,
  summary_stats           = prep$summary_stats,
  peak_bin_width          = PEAK_BIN_WIDTH,
  peak_time_origin        = PEAK_TIME_ORIGIN
)

cat("Fitting parameters : ", paste(PARAM_NAMES, collapse = ", "), "\n", sep = "")
cat("Fixed parameters   : ",
    paste(sprintf("%s = %s", setdiff(COMBINED_PARAM_NAMES, PARAM_NAMES),
                  unlist(fixed_values[setdiff(COMBINED_PARAM_NAMES, PARAM_NAMES)])),
          collapse = ", "), "\n", sep = "")
cat("Fitting summaries  : ", paste(prep$summary_stats, collapse = ", "), "\n", sep = "")
cat("Observed targets   : ",
    paste(sprintf("%s = %g", names(observed_summaries), observed_summaries),
          collapse = ", "), "\n", sep = "")


# -----------------------------------------------------------------------------
# 4. BUILD BASE + TIME-VARYING ARGS; COMPUTE EFFICACY-INDEPENDENT R0 INVARIANTS
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
cat("\nCentral-efficacy check (npi_scaler s = 0):\n")
cat(sprintf("  ppe_efficacy = %.3f, etu_efficacy = %.3f\n",
            central_effs$ppe_efficacy, central_effs$etu_efficacy))
cat(sprintf("  D = %.4f, F = %.4f, Q_g = %.4f\n",
            D_from_invariants(R0_invariants, central_effs$etu_efficacy,
                              DEFAULT_SCALAR_INPUTS$general_hospital_quarantine_efficacy),
            F_from_invariants(R0_invariants, DEFAULT_SCALAR_INPUTS$safe_funeral_efficacy),
            R0_invariants$Q_g))
cat(sprintf("  prob_hcw_cond_*_hospital at hcw_risk_scalar = 1: %.3f\n",
            min(HCW_BASE_PROB * 1.0, 1.0)))


# -----------------------------------------------------------------------------
# 5. PRIOR PREDICTIVE CHECK (optional, slow) -- useful to size the peak target
# -----------------------------------------------------------------------------
# Especially handy here: prep$summary_stats already includes "peak_height" /
# "time_to_peak", so this prints the prior-implied peak weekly deaths -- which tells
# you the right ballpark for the still-blank West Africa peak targets. This block
# runs BEFORE the section-7 guard, so you can size the targets even while NA.
# set.seed(1)
# pp <- prior_predictive_check_combined(
#   n_draws       = 50,
#   prior_list    = priors,
#   fit_params    = PARAM_NAMES,
#   fixed_values  = fixed_values,
#   base          = base_args,
#   tv            = tv_args_model,
#   invariants    = R0_invariants,
#   npi_spec      = NPI_SPEC,
#   general_hospital_quarantine_efficacy =
#     DEFAULT_SCALAR_INPUTS$general_hospital_quarantine_efficacy,
#   safe_funeral_efficacy = DEFAULT_SCALAR_INPUTS$safe_funeral_efficacy,
#   hcw_base_prob = HCW_BASE_PROB,
#   parallel      = FALSE,
#   n_replicates  = 5,
#   seeding_cases = ABC_CONFIG$seeding_cases,
#   takeoff_death_threshold = ABC_CONFIG$takeoff_death_threshold,
#   summary_stats = prep$summary_stats,
#   bin_width     = PEAK_BIN_WIDTH,
#   time_origin   = PEAK_TIME_ORIGIN
# )
# print(pp); summary(pp)


# -----------------------------------------------------------------------------
# 6. PER-RUN OUTPUT DIRECTORY + WORKER CONFIG
# -----------------------------------------------------------------------------
# Tag each per-run output with the settings that drive cost/resolution.
RUN_TAG <- sprintf("NP%d_NS%d_NBREPS_%d_NBSIMUL_%d",
                   length(PARAM_NAMES), length(prep$summary_stats),
                   ABC_CONFIG$n_reps, ABC_SETTINGS$nb_simul)

ABC_OUTPUT_DIR <- make_abc_output_dir(
  base_dir    = ABC_OUTPUT_BASE,
  scenario_id = SCENARIO_ID,
  label       = ABC_OUTPUT_LABEL,
  suffix      = RUN_TAG
)
message("ABC outputs will be written to: ", ABC_OUTPUT_DIR)

# ABC_CONFIG (incl. fit_params, fixed_values, summary_stats, peak settings)
# travels to each worker via the saved config; bootstrap_abc_worker_combined()
# merges it over its defaults. functions_path is the COMBINED file; the worker
# derives the npi/common paths from its directory.
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
# Fill-in reminder: time_to_peak / peak_height are intentionally left blank (NA)
# in OBSERVED_NAMED (section 1). Stop here -- with a clear message -- while any
# fitted target is still NA, so an unfilled value can't propagate into
# ABC_sequential(). (Sections 4-5 above still run, so you can size the targets with
# the prior predictive check first.)
na_targets <- names(observed_summaries)[is.na(observed_summaries)]
if (length(na_targets) > 0L) {
  stop("Observed target(s) still blank (NA): ", paste(na_targets, collapse = ", "),
       ".\n  Fill them in in OBSERVED_NAMED (section 1) before running this script.",
       call. = FALSE)
}

start_time <- Sys.time()
result <- with_abc_output_dir(
  ABC_OUTPUT_DIR,
  ABC_sequential(
    method              = ABC_SETTINGS$method,
    model               = fiber_abc_model_parallel_combined,  # returns summary_stats in order
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

result_stamp <- format(start_time, "%Y%m%d_%H%M%S")
result_filename <- paste0(
  "fiber_ABC_SMC_", SCENARIO_ID,
  if (nzchar(ABC_OUTPUT_LABEL)) paste0("_", ABC_OUTPUT_LABEL) else "",
  "_", result_stamp, "_", RUN_TAG, ".rds"
)

# Self-describing provenance. Reuses the NPI metadata helpers (npi_spec is needed
# to read npi_scaler back); the combined-scheme specifics ride in `extra`.
run_metadata <- make_npi_run_metadata(
  npi_spec           = NPI_SPEC,
  scenario_id        = SCENARIO_ID,
  priors             = priors,
  observed_summaries = observed_summaries,
  param_names        = PARAM_NAMES,
  general_hospital_quarantine_efficacy =
    DEFAULT_SCALAR_INPUTS$general_hospital_quarantine_efficacy,
  safe_funeral_efficacy = DEFAULT_SCALAR_INPUTS$safe_funeral_efficacy,
  extra = list(
    approach        = "combined",
    result_filename = result_filename,
    fit_params      = PARAM_NAMES,
    fixed_values    = fixed_values,
    hcw_base_prob   = HCW_BASE_PROB,
    summary_stats   = prep$summary_stats,
    peak_settings   = list(bin_width = PEAK_BIN_WIDTH, time_origin = PEAK_TIME_ORIGIN),
    abc_settings    = ABC_SETTINGS,
    abc_config      = ABC_CONFIG
  )
)
result <- attach_npi_run_metadata(result, run_metadata)

saveRDS(result, file = file.path(ABC_OUTPUT_DIR, result_filename))
saveRDS(result, file = file.path(FINAL_OUTPUTS_DIR, result_filename))
write_npi_run_metadata(run_metadata, file.path(ABC_OUTPUT_DIR, result_filename))
write_npi_run_metadata(run_metadata, file.path(FINAL_OUTPUTS_DIR, result_filename))


# -----------------------------------------------------------------------------
# 8. POSTERIOR INSPECTION
# -----------------------------------------------------------------------------
posterior <- as.data.frame(result$param)
colnames(posterior) <- PARAM_NAMES

print(apply(posterior, 2, quantile, probs = c(0.025, 0.5, 0.975)))

# Per-draw values of npi_scaler / hcw_risk_scalar whether or not they were fitted
# (fitted -> posterior column; fixed -> constant), so the implied efficacies and
# HCW exposure print correctly for any FIT_PARAMS choice.
draw_values <- function(name) {
  if (name %in% PARAM_NAMES) posterior[[name]] else rep(fixed_values[[name]], nrow(posterior))
}
scaler_draws <- draw_values("npi_scaler")
hcwrs_draws  <- draw_values("hcw_risk_scalar")

implied_ppe      <- vapply(scaler_draws, function(s) npi_efficacy_from_scaler(s, NPI_SPEC)$ppe_efficacy, numeric(1))
implied_etu      <- vapply(scaler_draws, function(s) npi_efficacy_from_scaler(s, NPI_SPEC)$etu_efficacy, numeric(1))
implied_prob_hcw <- pmin(HCW_BASE_PROB * hcwrs_draws, 1.0)

qfmt <- function(x) paste(round(quantile(x, c(0.025, 0.5, 0.975)), 3), collapse = " / ")
cat("\nImplied ppe_efficacy (2.5/50/97.5%):  ", qfmt(implied_ppe), "\n")
cat("Implied etu_efficacy (2.5/50/97.5%):  ", qfmt(implied_etu), "\n")
cat("Implied prob_hcw_hospital (2.5/50/97.5%): ", qfmt(implied_prob_hcw), "\n")

par(mfrow = c(1, length(PARAM_NAMES)))
for (j in seq_len(ncol(posterior))) {
  hist(posterior[, j], breaks = 10, main = colnames(posterior)[j],
       xlab = colnames(posterior)[j])
  abline(v = quantile(posterior[, j], c(0.025, 0.5, 0.975)),
         lty = c(2, 1, 2), col = "red")
}
par(mfrow = c(1, 1))

# Pairwise scatter to eyeball identifiability ridges (e.g. R0 <-> npi_scaler):
# a tight diagonal band = a confounded pair. peak_height/time_to_peak are in the
# fitted set precisely to break the R0 <-> etu_efficacy size-confound.
if (length(PARAM_NAMES) >= 2L) {
  pairs(posterior, main = "Posterior pairs (look for ridges)",
        col = adjustcolor("steelblue", alpha = 0.4), pch = 16, cex = 0.5)
}


# -----------------------------------------------------------------------------
# 9. PROGRESS / RECONSTRUCTION FROM DISK
# -----------------------------------------------------------------------------
# Sections 10-12 below run off the in-memory `result` from section 7. To re-inspect
# (or reload) a COMPLETED run from disk later, point ABC_OUTPUT_DIR at that run's
# directory and uncomment. Pass the fitted PARAM_NAMES + SUMMARY_STATS so the disk-
# inspection helpers label (and summarise) the right columns.
# ABC_OUTPUT_DIR <- file.path(ABC_OUTPUT_BASE, "abc_outputs", "<completed run directory>")
# abc_progress(ABC_OUTPUT_DIR, tolerance_target = ABC_SETTINGS$tolerance_target,
#              param_names = PARAM_NAMES, stat_names = prep$summary_stats)
# print(abc_compare_steps(ABC_OUTPUT_DIR, param_names = PARAM_NAMES,
#                         stat_names = prep$summary_stats))
# result <- reconstruct_abc_result(ABC_OUTPUT_DIR, param_names = PARAM_NAMES,
#                                  stat_names = prep$summary_stats)


# -----------------------------------------------------------------------------
# 10. POSTERIOR PREDICTIVE CHECKS (all fitted summaries)
# -----------------------------------------------------------------------------
sim_stats <- as.data.frame(result$stats)
colnames(sim_stats) <- prep$summary_stats

set.seed(1)
idx <- sample(seq_len(nrow(sim_stats)), size = 10000, replace = TRUE,
              prob = result$weights)
sim_stats_post <- sim_stats[idx, ]

n_stat <- length(observed_summaries)
par(mfrow = c(ceiling(n_stat / 2), 2), mar = c(4, 4, 3, 1))
for (s in names(observed_summaries)) {
  x  <- sim_stats_post[[s]]
  qs <- quantile(x, probs = c(0.025, 0.5, 0.975))
  hist(x, breaks = 10, main = paste0("Posterior predictive: ", s), xlab = s,
       col = adjustcolor("steelblue", alpha = 0.6), border = "white")
  abline(v = qs, col = "darkblue", lty = c(2, 1, 2), lwd = c(1, 2, 1))
  abline(v = observed_summaries[s], col = "red", lwd = 2.5)
}
par(mfrow = c(1, 1))

# Single-panel summary: simulated / observed for every fitted summary statistic.
plot(NA, xlim = c(0.5, n_stat + 0.5), ylim = c(0, 1.5),
     xaxt = "n", xlab = "", ylab = "Simulated / Observed",
     main = "Posterior-predictive fit ratio")
axis(1, at = seq_len(n_stat), labels = names(observed_summaries), las = 2, cex.axis = 0.8)
abline(h = 1, lty = 2, col = "red")
for (i in seq_len(n_stat)) {
  s  <- names(observed_summaries)[i]
  x  <- sim_stats_post[[s]] / observed_summaries[s]
  qs <- quantile(x, c(0.025, 0.5, 0.975))
  segments(i, qs[1], i, qs[3], lwd = 2, col = "darkblue")
  points(i, qs[2], pch = 16, cex = 1.5, col = "darkblue")
}


# -----------------------------------------------------------------------------
# 11. POSTERIOR TRAJECTORY CHECKS
# -----------------------------------------------------------------------------
# Draw N_TRAJ parameter sets from the weighted posterior, simulate one stochastic
# outbreak each, bin weekly incidence of (a) new infections and (b) deaths, and
# plot every trajectory plus the across-draw median + 95% CrI. Each draw's full
# 4-vector is reconstructed with assemble_combined_theta() (fitted from the
# posterior, the rest from fixed_values), so this works for any FIT_PARAMS choice.
library(ggplot2)
library(future.apply)
source(file.path(FUNCTIONS_DIR, "simulation_helpers.R"))   # bin_counts()

N_TRAJ    <- 200L
BIN_WIDTH <- PEAK_BIN_WIDTH
TRAJ_SEED <- 100L

set.seed(TRAJ_SEED)
traj_idx <- sample(seq_len(nrow(posterior)), size = N_TRAJ, replace = TRUE,
                   prob = result$weights)

plan(multisession, workers = min(N_CLUSTER, future::availableCores()))
traj_runs <- future_lapply(seq_along(traj_idx), function(k) {
  i    <- traj_idx[k]
  full <- assemble_combined_theta(as.numeric(posterior[i, PARAM_NAMES]),
                                  PARAM_NAMES, fixed_values)
  a <- build_abc_model_args_combined(
    R0 = full$R0, prop_funeral = full$prop_funeral,
    npi_scaler = full$npi_scaler, hcw_risk_scalar = full$hcw_risk_scalar,
    base = base_args, tv = tv_args_model, invariants = R0_invariants,
    npi_spec = NPI_SPEC,
    general_hospital_quarantine_efficacy = DEFAULT_SCALAR_INPUTS$general_hospital_quarantine_efficacy,
    safe_funeral_efficacy = DEFAULT_SCALAR_INPUTS$safe_funeral_efficacy,
    hcw_base_prob = HCW_BASE_PROB, seeding_cases = ABC_CONFIG$seeding_cases)
  a$seed <- TRAJ_SEED + k
  tdf <- do.call(branching_process_main, a)$tdf
  tdf <- tdf[!is.na(tdf$time_infection_absolute), , drop = FALSE]
  died <- !is.na(tdf$outcome) & tdf$outcome
  list(case_days  = as.integer(floor(tdf$time_infection_absolute)),
       death_days = as.integer(floor(tdf$time_outcome_absolute[died])))
}, future.packages = "fiber", future.seed = TRUE)
plan(sequential)

# Common weekly grid spanning the longest trajectory, then bin every run onto it.
max_day <- max(unlist(lapply(traj_runs, function(r) c(r$case_days, r$death_days))), 0)
n_bins  <- max(1L, ceiling((max_day + 1L) / BIN_WIDTH))
week    <- ((seq_len(n_bins) - 1L) + 0.5) * BIN_WIDTH

traj_long <- do.call(rbind, lapply(seq_along(traj_runs), function(k) {
  r <- traj_runs[[k]]
  rbind(
    data.frame(draw = k, week = week, metric = "Cases",
               incidence = bin_counts(r$case_days,  BIN_WIDTH, n_bins)),
    data.frame(draw = k, week = week, metric = "Deaths",
               incidence = bin_counts(r$death_days, BIN_WIDTH, n_bins))
  )
}))

# (1) Individual trajectories (one faint line per posterior draw), faceted by metric.
print(
  ggplot(traj_long, aes(week, incidence, group = draw)) +
    geom_line(alpha = 0.15, colour = "steelblue", linewidth = 0.3) +
    facet_wrap(~ metric, scales = "free_y") +
    labs(x = "Time since outbreak start (days)", y = "Count per week",
         title = sprintf("Posterior trajectories: %s", SCENARIO_ID),
         subtitle = sprintf("%d posterior draws, one stochastic outbreak each", N_TRAJ)) +
    theme_bw()
)

# (2) Across-draw median + 95% credible interval per week, faceted by metric.
band <- do.call(rbind, lapply(split(traj_long, list(traj_long$metric, traj_long$week), drop = TRUE),
  function(d) data.frame(metric = d$metric[1], week = d$week[1],
                         lo  = quantile(d$incidence, 0.025, names = FALSE),
                         med = quantile(d$incidence, 0.500, names = FALSE),
                         hi  = quantile(d$incidence, 0.975, names = FALSE))))
print(
  ggplot(band, aes(week, med)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "steelblue", alpha = 0.3) +
    geom_line(colour = "steelblue", linewidth = 0.8) +
    facet_wrap(~ metric, scales = "free_y") +
    labs(x = "Time since outbreak start (days)", y = "Count per week",
         title = sprintf("Posterior median and 95%% CrI: %s", SCENARIO_ID),
         subtitle = sprintf("Across %d posterior draws (parameter + stochastic uncertainty)", N_TRAJ)) +
    theme_bw()
)


# -----------------------------------------------------------------------------
# 12. PEAK-METRIC POSTERIOR PREDICTIVE (curve-shape summaries, if fitted)
# -----------------------------------------------------------------------------
# Recompute time_to_peak / peak_height for each posterior-draw trajectory with the
# SAME helper the fit used (peak_stats_from_death_days(), common file) and compare
# to the observed targets. Only the curve summaries actually in SUMMARY_STATS are
# shown -- here that's time_to_peak and peak_height, the fitted curve summaries
# (compared against the targets you fill into OBSERVED_NAMED).
peak_metrics_fitted <- intersect(c("time_to_peak", "peak_height"), prep$summary_stats)
if (length(peak_metrics_fitted) > 0L) {
  peak_post <- do.call(rbind, lapply(traj_runs, function(r) {
    if (length(r$death_days) == 0L) {
      return(data.frame(time_to_peak = NA_real_, peak_height = 0))
    }
    ps <- peak_stats_from_death_days(
      death_days      = r$death_days,
      first_death_day = min(r$death_days),
      bin_width       = PEAK_BIN_WIDTH,
      time_origin     = PEAK_TIME_ORIGIN
    )
    data.frame(time_to_peak = ps[["time_to_peak"]], peak_height = ps[["peak_height"]])
  }))

  cat("\nPeak-metric posterior predictive (median / 95% CrI vs observed):\n")
  for (s in peak_metrics_fitted) {
    qs <- quantile(peak_post[[s]], c(0.025, 0.5, 0.975), na.rm = TRUE)
    cat(sprintf("  %-13s %.1f / [%.1f, %.1f]   (observed %.1f)\n",
                paste0(s, ":"), qs[2], qs[1], qs[3], observed_summaries[s]))
  }

  par(mfrow = c(1, length(peak_metrics_fitted)), mar = c(4, 4, 3, 1))
  for (s in peak_metrics_fitted) {
    x  <- peak_post[[s]][is.finite(peak_post[[s]])]
    qs <- quantile(x, probs = c(0.025, 0.5, 0.975))
    hist(x, breaks = 12, main = paste0("Trajectory ", s), xlab = s,
         col = adjustcolor("seagreen", alpha = 0.6), border = "white")
    abline(v = qs, col = "darkgreen", lty = c(2, 1, 2), lwd = c(1, 2, 1))
    abline(v = observed_summaries[s], col = "red", lwd = 2.5)
  }
  par(mfrow = c(1, 1))
}
