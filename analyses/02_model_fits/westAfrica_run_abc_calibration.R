# westAfrica_run_abc_calibration.R
# =============================================================================
# Phase 4 ABC-SMC calibration of the fiber branching-process model against
# the 2014-16 West Africa Ebola outbreak (Worst_WestAfrica scenario).
#
# Fitted parameters (3):
#   R0               : baseline reproduction number for a typical genPop
#                      seeding case (evaluated at t = 0 using the scenario's
#                      unsafe-funeral and hospital-quarantine settings).
#   prop_funeral     : share of R0 attributable to funeral transmission at
#                      t = 0.
#   hcw_risk_scalar  : multiplier (capped at 1.0) applied to a symmetric
#                      hcw_base_prob = 0.25 to produce both
#                      prob_hcw_cond_*_hospital probabilities.
#
# Observed summaries (3, plus a "took off" indicator to handle bimodality):
#   takeoff      : 1.0 (observed outbreak did take off)
#   n_deaths     : total deaths
#   n_hcw_deaths : HCW deaths
#   duration     : first death -> last outcome, in days
#
# Phase 4 production targets (Worst_WestAfrica):
#   nb_simul         = 220
#   n_reps           = 60
#   tolerance_target = 0.27  (~0.04 above the N=60 noise floor of ~0.23)
#   hcw_base_prob    = 0.25  (symmetric base for both HCW-hospital probs)
#   n_cluster        = 110   (PETAL workstation)
# Expected ~9 SMC steps, ~46 h wall time on PETAL.
#
# Sections:
#   1. Configuration (paths, scenario, ABC tuning)
#   2. Libraries + sources
#   3. Build base + time-varying args; compute D, F multipliers
#   4. Observed targets and priors
#   5. (Optional) Prior predictive check
#   6. Per-run output directory + worker config (FIBER_ABC_CONFIG)
#   7. Run ABC_sequential (Del Moral et al. 2012)
#   8. Save result; inspect posterior
#   9. (Optional) Reconstruct from disk + monitor progress mid-run
#  10. Posterior predictive checks
# =============================================================================


# -----------------------------------------------------------------------------
# 1. CONFIGURATION
# -----------------------------------------------------------------------------
# ANALYSIS_DIR is this script's containing directory (analyses/02_model_fits
# inside the obv_hcw_paper repo). The script assumes you've setwd() to this
# directory before sourcing; the machine-specific switch below baked in
# absolute defaults for the most common workstations.

ANALYSIS_DIR <- switch(
  Sys.info()[["user"]],
  "cwhittaker" = "C:/Users/cwhittaker/Documents/Research Projects/obv_hcw_paper/analyses/02_model_fits",
  "PETAL_WS_2" = "C:/Users/PETAL_WS_2/Documents/obv_hcw_paper/analyses/02_model_fits",
  getwd()
)
HELPER_DIR     <- file.path(ANALYSIS_DIR, "helper_functions")
SETUP_PATH     <- file.path(HELPER_DIR,   "setup_model_parameters.R")
FUNCTIONS_PATH <- file.path(HELPER_DIR,   "abc_calibration_functions.R")
R0_PATH        <- file.path(HELPER_DIR,   "calculate_model_approx_r0.R")
SCENARIO_CSV   <- file.path(ANALYSIS_DIR, "final_four_scenario_values.csv")
SCENARIO_ID    <- "Worst_WestAfrica"

# Any scalar-parameter overrides to layer on top of DEFAULT_SCALAR_INPUTS.
# Pass anything you want to differ from the literature-informed defaults
# here as a named list; new parameters that get added to the model in the
# future automatically become overridable without further code changes.
MODEL_OVERRIDES <- list(
  check_final_size = 30000
)

# Where ABC_sequential's intermediate files (output_step*, tolerance_step*,
# n_simul_tot_step*) and the final result RDS get written. Each call to
# make_abc_output_dir() in section 6 creates a fresh timestamped
# subdirectory of <ABC_OUTPUT_BASE>/abc_outputs/ tagged with
# ABC_OUTPUT_LABEL, so successive phases / runs don't overwrite each other.
ABC_OUTPUT_BASE  <- ANALYSIS_DIR
ABC_OUTPUT_LABEL <- "final"

# Repo-level outputs/ folder. The final-result RDS is copied here in addition
# to being written under ABC_OUTPUT_DIR, so manuscript-ready artefacts live
# in one canonical place.
FINAL_OUTPUTS_DIR <- normalizePath(
  file.path(ANALYSIS_DIR, "..", "..", "outputs"),
  mustWork = FALSE
)
if (!dir.exists(FINAL_OUTPUTS_DIR)) {
  dir.create(FINAL_OUTPUTS_DIR, recursive = TRUE, showWarnings = FALSE)
}

# Symmetric base for both prob_hcw_cond_*_hospital probabilities. The fitted
# hcw_risk_scalar multiplies this for both, capped at 1.0 — see
# build_abc_model_args() in abc_calibration_functions.R.
HCW_BASE_PROB <- 0.25

# ABC tuning. These travel into each worker via bootstrap_abc_worker().
ABC_CONFIG <- list(
  check_final_size        = 30000,
  takeoff_death_threshold = 100,         # >= K deaths counts as a take-off
  n_reps                  = 60,          # replicates per particle (per theta)
  seeding_cases           = 25,
  hcw_base_prob           = HCW_BASE_PROB,
  setup_R0_n              = 100000L,
  setup_R0_seed           = 42L,
  setup_funeral_share     = 0.5          # only used to choose D and F units
)

# EasyABC::ABC_sequential settings.
ABC_SETTINGS <- list(
  method              = "Delmoral",
  nb_simul            = 220,
  alpha               = 0.5,
  tolerance_target    = 0.27,            # ~0.04 above the N=60 noise floor
  M                   = 1,
  use_seed            = TRUE,
  verbose             = TRUE
)

# Worker count. Aggressive on the PETAL box (Phase 4 target 110), modest
# on dev workstations.
N_CLUSTER <- if (grepl("PETAL", Sys.info()[["user"]], ignore.case = TRUE)) {
  min(110, parallel::detectCores() - 10)
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
# 3. BUILD BASE + TIME-VARYING ARGS, SOLVE FOR D / F MULTIPLIERS
# -----------------------------------------------------------------------------
# Computed ONCE for the main process. Workers will rebuild their own copies
# from the configuration above in section 6.

scenario_matrix <- read_scenario_matrix(SCENARIO_CSV)

mp <- make_model_parameters(
  scenario_id     = SCENARIO_ID,
  scenario_matrix = scenario_matrix,
  overrides       = MODEL_OVERRIDES
)
base_args      <- mp$base_args
tv_args_model  <- mp$tv_args

# Sanity glance at the prob_hosp curve.
plot(
  scenario_matrix$relative_day[scenario_matrix$scenario == SCENARIO_ID],
  scenario_matrix$prob_hosp[scenario_matrix$scenario == SCENARIO_ID],
  xlab = "Relative day", ylab = "P(hospitalised)",
  main = paste0(SCENARIO_ID, ": prob_hosp(t)")
)

setup_solve <- solve_offspring_means_for_R0(
  R0   = 1.0,
  args = mp$args,
  proportion_transmission_from_funerals = ABC_CONFIG$setup_funeral_share,
  n    = ABC_CONFIG$setup_R0_n,
  seed = ABC_CONFIG$setup_R0_seed
)
D_direct_multiplier  <- setup_solve$D_direct_multiplier
F_funeral_multiplier <- setup_solve$F_funeral_multiplier


# -----------------------------------------------------------------------------
# 4. OBSERVED TARGETS AND PRIORS
# -----------------------------------------------------------------------------

observed_summaries <- c(
  takeoff      = 1.0,
  n_deaths     = 11325,
  n_hcw_deaths = 513,
  duration     = 365    # spatial heterogeneity sustained it; main outbreak ~ a year
)

# Phase 4 priors. R0 and prop_funeral inherited unchanged from Phase 3;
# hcw_risk_scalar widened/shifted to (0.50, 4.00) to reflect the new
# symmetric HCW_BASE_PROB = 0.25 base.
priors <- list(
  c("unif", 1.35, 1.55),   # R0
  c("unif", 0.10, 0.40),   # prop_funeral
  c("unif", 0.50, 4.00)    # hcw_risk_scalar
)


# -----------------------------------------------------------------------------
# 5. PRIOR PREDICTIVE CHECK (optional, slow)
# -----------------------------------------------------------------------------
# Run a small prior predictive check before launching the full ABC. This
# uses the in-process (sequential) model wrapper and is intentionally cheap.
#
# To parallelise, set up a future plan first; the function will then use
# future_lapply when parallel = TRUE.
#
# set.seed(1)
# pp <- prior_predictive_check(
#   n_draws       = 20,
#   prior_list    = priors,
#   base          = base_args,
#   tv            = tv_args_model,
#   D             = D_direct_multiplier,
#   F_fun         = F_funeral_multiplier,
#   parallel      = FALSE,
#   n_replicates  = 5,
#   seeding_cases = ABC_CONFIG$seeding_cases,
#   takeoff_death_threshold = ABC_CONFIG$takeoff_death_threshold,
#   hcw_base_prob = HCW_BASE_PROB
# )
# print(pp)
# summary(pp)


# -----------------------------------------------------------------------------
# 6. PER-RUN OUTPUT DIRECTORY + WORKER CONFIG
# -----------------------------------------------------------------------------
# Make a fresh per-run subdirectory for everything ABC_sequential writes.
# ABC_sequential() spawns its own PSOCK cluster, so we cannot rely on
# clusterExport. Instead, serialise the everything-workers-need-to-know
# config to disk and advertise it via the FIBER_ABC_CONFIG env var, which
# PSOCK workers inherit from this R session. On its first call each worker
# reads the config and self-bootstraps.

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
# Wrap the call so cwd is temporarily ABC_OUTPUT_DIR while step files are
# written; the workers' cwd is independent.

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

result_filename <- paste0(
  "fiber_abc_smc_result_", SCENARIO_ID,
  if (nzchar(ABC_OUTPUT_LABEL)) paste0("_", ABC_OUTPUT_LABEL) else "",
  ".rds"
)
saveRDS(result, file = file.path(ABC_OUTPUT_DIR, result_filename))
saveRDS(result, file = file.path(FINAL_OUTPUTS_DIR, result_filename))


# -----------------------------------------------------------------------------
# 8. POSTERIOR INSPECTION
# -----------------------------------------------------------------------------

posterior <- as.data.frame(result$param)
colnames(posterior) <- c("R0", "prop_funeral", "hcw_risk_scalar")

print(apply(posterior, 2, quantile, probs = c(0.025, 0.5, 0.975)))

par(mfrow = c(1, 3))
for (j in seq_len(ncol(posterior))) {
  hist(posterior[, j], breaks = 15,
       main = colnames(posterior)[j],
       xlab = colnames(posterior)[j])
  abline(v = quantile(posterior[, j], c(0.025, 0.5, 0.975)),
         lty = c(2, 1, 2), col = "red")
}
par(mfrow = c(1, 1))


# -----------------------------------------------------------------------------
# 9. PROGRESS / RECONSTRUCTION FROM DISK
# -----------------------------------------------------------------------------
# These functions work mid-run (peek at progress files) or after the fact
# (rebuild an ABC_sequential()-style result object from output_step*). Point
# them at ABC_OUTPUT_DIR for the current run, or at any previous run's
# subdirectory under <ABC_OUTPUT_BASE>/abc_outputs/.
#
# abc_progress(ABC_OUTPUT_DIR, tolerance_target = ABC_SETTINGS$tolerance_target)
# print(abc_compare_steps(ABC_OUTPUT_DIR))
#
# # Inspect the final step's particle cloud directly:
# last_step_file <- tail(list.files(ABC_OUTPUT_DIR, pattern = "^output_step[0-9]+$",
#                                   full.names = TRUE), 1L)
# step_last <- read.table(last_step_file, header = FALSE)
# colnames(step_last) <- c("weight", "R0", "prop_funeral", "hcw_risk_scalar",
#                          "takeoff", "n_deaths", "n_hcw_deaths", "duration")
# summary(step_last)
#
# # Reconstruct an ABC result object from the latest completed step:
# # result <- reconstruct_abc_result(ABC_OUTPUT_DIR)
# # Or from a specific step:
# # result <- reconstruct_abc_result(ABC_OUTPUT_DIR, step = 3)


# -----------------------------------------------------------------------------
# 10. POSTERIOR PREDICTIVE CHECKS
# -----------------------------------------------------------------------------
# Resample particles by their weights so the histograms reflect the
# weighted posterior.

sim_stats <- as.data.frame(result$stats)
colnames(sim_stats) <- c("takeoff", "n_deaths", "n_hcw_deaths", "duration")

set.seed(1)
idx <- sample(seq_len(nrow(sim_stats)),
              size    = 10000,
              replace = TRUE,
              prob    = result$weights)
sim_stats_post <- sim_stats[idx, ]

par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (s in names(observed_summaries)) {
  x  <- sim_stats_post[[s]]
  qs <- quantile(x, probs = c(0.025, 0.5, 0.975))

  hist(x,
       breaks = 10,
       main   = paste0("Posterior predictive: ", s),
       xlab   = s,
       col    = adjustcolor("steelblue", alpha = 0.6),
       border = "white")
  abline(v = qs, col = "darkblue", lty = c(2, 1, 2), lwd = c(1, 2, 1))
  abline(v = observed_summaries[s], col = "red", lwd = 2.5)
  legend("topleft",
         legend = c(paste0("Observed: ", signif(observed_summaries[s], 3)),
                    paste0("Median: ",   signif(qs[2], 3)),
                    paste0("95% CI: [",  signif(qs[1], 3), ", ",
                           signif(qs[3], 3), "]")),
         col = c("red", "darkblue", NA),
         lty = c(1, 1, NA),
         lwd = c(2.5, 2, NA),
         bty = "n",
         cex = 0.75)
}
par(mfrow = c(1, 1))

# Single-panel summary: simulated vs observed.
plot(NA, xlim = c(0.5, 4.5), ylim = c(0, 1.5),
     xaxt = "n", xlab = "", ylab = "Simulated / Observed",
     main = "Posterior-predictive fit ratio")
axis(1, at = 1:4, labels = names(observed_summaries))
abline(h = 1, lty = 2, col = "red")
for (i in seq_along(observed_summaries)) {
  x  <- sim_stats_post[[names(observed_summaries)[i]]] /
        observed_summaries[names(observed_summaries)[i]]
  qs <- quantile(x, c(0.025, 0.5, 0.975))
  segments(i, qs[1], i, qs[3], lwd = 2, col = "darkblue")
  points(i, qs[2], pch = 16, cex = 1.5, col = "darkblue")
}
