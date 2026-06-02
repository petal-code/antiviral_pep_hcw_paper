# DRC_run_abc_calibration_peak.R  (NPI-efficacy approach, EXTENDED summaries)
# =============================================================================
# EXPERIMENTAL variant of DRC_run_abc_calibration.R. Same model, same scenario,
# same THREE fitted parameters (R0, prop_funeral, npi_scaler) and same priors --
# the ONLY change is that ABC is asked to match two extra features of the weekly
# DEATH-incidence curve on top of the usual four scalar summaries:
#
#   Existing (4):  takeoff, n_deaths, n_hcw_deaths, duration
#   Added    (2):  time_to_peak  -- time to the peak week of death incidence
#                  peak_height   -- deaths in that peak week (deaths / week)
#
# The extended summary/model machinery lives in
#   functions/abc_calibration_functions_npi_peak.R
# (abc_summarise_peak(), fiber_abc_model_parallel_peak(), ...), which reuses the
# base NPI scheme's per-particle mapping. Run this script exactly like the base
# one; outputs are tagged "NPIeffPeak" so they sit next to -- and never clobber --
# the 4-summary fits.
#
# Fitted parameters (3, unchanged):
#   R0           : baseline reproduction number for a typical genPop seeding
#                  case at t = 0 (under this particle's efficacies).
#   prop_funeral : share of R0 attributable to funeral transmission at t = 0.
#   npi_scaler   : SINGLE shared scaler in [-1, 1] positioning the two fitted
#                  conditional efficacies (ppe_efficacy, etu_efficacy) on their
#                  [min, max] intervals. s = -1 -> min, 0 -> central, +1 -> max.
#
# FIXED (not fitted): general_hospital_quarantine_efficacy, safe_funeral_efficacy
#   default to DEFAULT_SCALAR_INPUTS in functions/setup_model_parameters.R.
#   prob_hcw_cond_*_hospital fixed at 0.25 (no HCW-risk scalar).
#
# >>> PLACEHOLDERS <<<  NPI_SPEC bounds AND the two new observed targets
#   (time_to_peak, peak_height) are PLACEHOLDER values -- see notes at each.
#   UPDATE WITH REAL NUMBERS (peak week + peak weekly deaths from the DRC weekly
#   surveillance series) before any production run.
# =============================================================================


# -----------------------------------------------------------------------------
# 1. CONFIGURATION
# -----------------------------------------------------------------------------
# Paths resolved from the repo root with here::here() (locates obv_hcw_paper.Rproj),
# so the script runs the same regardless of working directory / machine.
ANALYSIS_DIR   <- here::here("analyses", "02_ABC_model_fits_NPI_Eff")
FUNCTIONS_DIR  <- here::here("functions")
# Shared model code + the NPI-specific ABC helpers + the EXTENDED (peak) helpers.
# The worker config below points `functions_path` at the PEAK file; its
# self-bootstrap sources the base NPI file (for build_abc_model_args()) and the
# generic common file alongside it (see bootstrap_abc_worker_peak()).
SETUP_PATH     <- file.path(FUNCTIONS_DIR, "setup_model_parameters.R")
COMMON_PATH    <- file.path(FUNCTIONS_DIR, "abc_calibration_functions_common.R")
NPI_PATH       <- file.path(FUNCTIONS_DIR, "abc_calibration_functions_npi.R")
PEAK_PATH      <- file.path(FUNCTIONS_DIR, "abc_calibration_functions_npi_peak.R")
R0_PATH        <- file.path(FUNCTIONS_DIR, "calculate_model_approx_r0.R")
SCENARIO_CSV   <- here::here("data-processed", "final_six_scenario_values_original_approach.csv")
SCENARIO_ID    <- "Middle_DRC_ConflictSmoothed_PlusPlus" # either "Middle_DRC_ConflictSmoothed" or "Middle_DRC_ConflictSmoothed_PlusPlus"

# Names of the 3 fitted parameters and the 6 fitted summary statistics, in the
# order the worker returns them. Reused throughout the post-processing so the
# disk-inspection helpers and column labels stay consistent.
PARAM_NAMES <- c("R0", "prop_funeral", "npi_scaler")
STAT_NAMES  <- c("takeoff", "n_deaths", "n_hcw_deaths", "duration",
                 "time_to_peak", "peak_height")

# >>> PLACEHOLDER: bounds for the two fitted efficacies (Approach B). UPDATE. <<<
NPI_SPEC <- list(
  ppe_efficacy = list(min = 0.30, max = 0.90),   # PLACEHOLDER
  etu_efficacy = list(min = 0.60, max = 0.95)    # PLACEHOLDER
)

MODEL_OVERRIDES <- list(check_final_size = 10000)

ABC_OUTPUT_BASE  <- ANALYSIS_DIR
ABC_OUTPUT_LABEL <- "NPIeffPeak"   # distinct from the base scheme's "NPIeff"

FINAL_OUTPUTS_DIR <- here::here("outputs", "02_ABC_model_fits_NPI_Eff")
if (!dir.exists(FINAL_OUTPUTS_DIR)) {
  dir.create(FINAL_OUTPUTS_DIR, recursive = TRUE, showWarnings = FALSE)
}

# Weekly DEATH-incidence binning for the two new summaries.
#   peak_bin_width   : days per bin (7 = weekly, matching the trajectory plots).
#   peak_time_origin : "first_death" (whole weeks from the first-death week to the
#                      peak week; anchored like `duration`, robust to takeoff lag)
#                      or "outbreak_start" (peak-week midpoint day from t = 0).
PEAK_BIN_WIDTH   <- 7L
PEAK_TIME_ORIGIN <- "first_death"

# general_hospital_quarantine_efficacy and safe_funeral_efficacy are omitted, so
# bootstrap_abc_worker_peak() fills them from DEFAULT_SCALAR_INPUTS.
ABC_CONFIG <- list(
  takeoff_death_threshold = 100,
  n_reps                  = 50, # 100
  seeding_cases           = 25,
  setup_R0_n              = 100000L,
  setup_R0_seed           = 42L,
  npi_spec                = NPI_SPEC,
  peak_bin_width          = PEAK_BIN_WIDTH,
  peak_time_origin        = PEAK_TIME_ORIGIN
)

ABC_SETTINGS <- list(
  method              = "Delmoral",
  nb_simul            = 220, #660
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
source(COMMON_PATH)
source(NPI_PATH)      # base NPI scheme: build_abc_model_args(), metadata helpers
source(PEAK_PATH)     # extended scheme: abc_summarise_peak(), parallel model, ...
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
                              DEFAULT_SCALAR_INPUTS$general_hospital_quarantine_efficacy),
            F_from_invariants(R0_invariants, DEFAULT_SCALAR_INPUTS$safe_funeral_efficacy),
            R0_invariants$Q_g))


# -----------------------------------------------------------------------------
# 4. OBSERVED TARGETS AND PRIORS
# -----------------------------------------------------------------------------
# Order MUST match STAT_NAMES and the worker's return vector.
observed_summaries <- c(
  takeoff      = 1.0,
  n_deaths     = 2299,
  n_hcw_deaths = 79,     # https://afenet-journal.org/10-37432-jieph-d-25-00072/
  duration     = 450,    # ~ Aug 2018 - Nov 2019 main phase
  # >>> PLACEHOLDERS <<< derive from the DRC WEEKLY death series:
  time_to_peak = 250,    # PLACEHOLDER days to the peak DEATH week (origin = PEAK_TIME_ORIGIN)
  peak_height  = 75      # PLACEHOLDER peak weekly deaths (deaths / week)
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
# pp <- prior_predictive_check_peak(
#   n_draws       = 50,
#   prior_list    = priors,
#   base          = base_args,
#   tv            = tv_args_model,
#   invariants    = R0_invariants,
#   npi_spec      = NPI_SPEC,
#   general_hospital_quarantine_efficacy =
#     DEFAULT_SCALAR_INPUTS$general_hospital_quarantine_efficacy,
#   safe_funeral_efficacy = DEFAULT_SCALAR_INPUTS$safe_funeral_efficacy,
#   parallel      = FALSE,
#   n_replicates  = 5,
#   seeding_cases = ABC_CONFIG$seeding_cases,
#   takeoff_death_threshold = ABC_CONFIG$takeoff_death_threshold,
#   bin_width     = PEAK_BIN_WIDTH,
#   time_origin   = PEAK_TIME_ORIGIN
# )
# print(pp); summary(pp)


# -----------------------------------------------------------------------------
# 6. PER-RUN OUTPUT DIRECTORY + WORKER CONFIG
# -----------------------------------------------------------------------------
# Tag every per-run output with the two settings that drive a run's cost and
# resolution -- n_reps and nb_simul -- so runs are self-describing on disk.
RUN_TAG <- sprintf("NBREPS_%d_NBSIMUL_%d", ABC_CONFIG$n_reps, ABC_SETTINGS$nb_simul)

ABC_OUTPUT_DIR <- make_abc_output_dir(
  base_dir    = ABC_OUTPUT_BASE,
  scenario_id = SCENARIO_ID,
  label       = ABC_OUTPUT_LABEL,
  suffix      = RUN_TAG
)
message("ABC outputs will be written to: ", ABC_OUTPUT_DIR)

# functions_path points at the PEAK file so each worker sources the extended
# scheme (and, via its bootstrap, the base NPI file alongside it).
save_abc_config(list(
  setup_path      = SETUP_PATH,
  functions_path  = PEAK_PATH,
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
    model               = fiber_abc_model_parallel_peak,   # 6-stat model
    prior               = priors,
    nb_simul            = ABC_SETTINGS$nb_simul,
    summary_stat_target = observed_summaries,              # length 6
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
# RUN_TAG goes AFTER the timestamp so the name still sorts chronologically.
result_filename <- paste0(
  "fiber_ABC_SMC_", SCENARIO_ID,
  if (nzchar(ABC_OUTPUT_LABEL)) paste0("_", ABC_OUTPUT_LABEL) else "",
  "_", result_stamp, "_", RUN_TAG, ".rds"
)
# Self-describing provenance: store the npi_spec (needed to read back npi_scaler),
# priors, fixed efficacies, targets, AND the peak-binning settings with the fit.
run_metadata <- make_npi_run_metadata(
  npi_spec           = NPI_SPEC,
  scenario_id        = SCENARIO_ID,
  priors             = priors,
  observed_summaries = observed_summaries,
  general_hospital_quarantine_efficacy =
    DEFAULT_SCALAR_INPUTS$general_hospital_quarantine_efficacy,
  safe_funeral_efficacy = DEFAULT_SCALAR_INPUTS$safe_funeral_efficacy,
  extra = list(result_filename = result_filename,
               abc_settings = ABC_SETTINGS, abc_config = ABC_CONFIG,
               stat_names = STAT_NAMES,
               peak_settings = list(bin_width = PEAK_BIN_WIDTH,
                                    time_origin = PEAK_TIME_ORIGIN))
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
# Pass the 6 STAT_NAMES so the disk-inspection helpers label the extra columns.
abc_progress(ABC_OUTPUT_DIR, tolerance_target = ABC_SETTINGS$tolerance_target,
             param_names = PARAM_NAMES, stat_names = STAT_NAMES)
print(abc_compare_steps(ABC_OUTPUT_DIR, param_names = PARAM_NAMES, stat_names = STAT_NAMES))
result <- reconstruct_abc_result(ABC_OUTPUT_DIR, param_names = PARAM_NAMES,
                                 stat_names = STAT_NAMES)


# -----------------------------------------------------------------------------
# 10. POSTERIOR PREDICTIVE CHECKS (all 6 summaries)
# -----------------------------------------------------------------------------
sim_stats <- as.data.frame(result$stats)
colnames(sim_stats) <- STAT_NAMES

set.seed(1)
idx <- sample(seq_len(nrow(sim_stats)), size = 10000, replace = TRUE,
              prob = result$weights)
sim_stats_post <- sim_stats[idx, ]

par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
for (s in names(observed_summaries)) {
  x  <- sim_stats_post[[s]]
  qs <- quantile(x, probs = c(0.025, 0.5, 0.975))
  hist(x, breaks = 10, main = paste0("Posterior predictive: ", s), xlab = s,
       col = adjustcolor("steelblue", alpha = 0.6), border = "white")
  abline(v = qs, col = "darkblue", lty = c(2, 1, 2), lwd = c(1, 2, 1))
  abline(v = observed_summaries[s], col = "red", lwd = 2.5)
}
par(mfrow = c(1, 1))

# Single-panel summary: simulated / observed for every summary statistic.
n_stat <- length(observed_summaries)
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
# plot every trajectory plus the across-draw median + 95% CrI. Folds in BOTH
# parameter and stochastic uncertainty.
library(ggplot2)
library(future.apply)
source(file.path(FUNCTIONS_DIR, "simulation_helpers.R"))   # bin_counts()

N_TRAJ    <- 200L  # posterior draws to simulate (one stochastic outbreak each)
BIN_WIDTH <- PEAK_BIN_WIDTH   # weekly, matching the fitted peak summaries
TRAJ_SEED <- 100L  # base seed for reproducible per-draw simulations

set.seed(TRAJ_SEED)
traj_idx <- sample(seq_len(nrow(posterior)), size = N_TRAJ, replace = TRUE,
                   prob = result$weights)

plan(multisession, workers = min(N_CLUSTER, future::availableCores()))
traj_runs <- future_lapply(seq_along(traj_idx), function(k) {
  i <- traj_idx[k]
  a <- build_abc_model_args(
    R0 = posterior$R0[i], prop_funeral = posterior$prop_funeral[i],
    npi_scaler = posterior$npi_scaler[i],
    base = base_args, tv = tv_args_model, invariants = R0_invariants,
    npi_spec = NPI_SPEC,
    general_hospital_quarantine_efficacy = DEFAULT_SCALAR_INPUTS$general_hospital_quarantine_efficacy,
    safe_funeral_efficacy = DEFAULT_SCALAR_INPUTS$safe_funeral_efficacy,
    seeding_cases = ABC_CONFIG$seeding_cases)
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
week    <- ((seq_len(n_bins) - 1L) + 0.5) * BIN_WIDTH          # bin midpoints (days)

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
# 12. PEAK-METRIC POSTERIOR PREDICTIVE (the two NEW summaries)
# -----------------------------------------------------------------------------
# Recompute time_to_peak / peak_height for each posterior-draw trajectory with the
# SAME helper the fit used (peak_stats_from_death_days()), and compare to the
# observed targets. This is the curve-level analogue of section 10's histograms
# and is the most direct check that the two added summaries are being matched.
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
for (s in c("time_to_peak", "peak_height")) {
  qs <- quantile(peak_post[[s]], c(0.025, 0.5, 0.975), na.rm = TRUE)
  cat(sprintf("  %-13s %.1f / [%.1f, %.1f]   (observed %.1f)\n",
              paste0(s, ":"), qs[2], qs[1], qs[3], observed_summaries[s]))
}

par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
for (s in c("time_to_peak", "peak_height")) {
  x  <- peak_post[[s]][is.finite(peak_post[[s]])]
  qs <- quantile(x, probs = c(0.025, 0.5, 0.975))
  hist(x, breaks = 12, main = paste0("Trajectory ", s), xlab = s,
       col = adjustcolor("seagreen", alpha = 0.6), border = "white")
  abline(v = qs, col = "darkgreen", lty = c(2, 1, 2), lwd = c(1, 2, 1))
  abline(v = observed_summaries[s], col = "red", lwd = 2.5)
}
par(mfrow = c(1, 1))
