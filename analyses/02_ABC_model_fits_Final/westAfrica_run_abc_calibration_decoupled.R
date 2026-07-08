# westAfrica_run_abc_calibration_decoupled.R  (DECOUPLED-efficacy approach)
# =============================================================================
# ABC-SMC calibration of the revamped fiber branching-process model against the
# Worst_WestAfrica Ebola scenario, using the DECOUPLED scheme
# (functions/abc_calibration_functions_decoupled.R):
#
#   FITTED PARAMETERS (5): R0, prop_funeral, etu_efficacy, ppe_efficacy,
#                          hcw_risk_scalar
#     * etu_efficacy and ppe_efficacy are fit DIRECTLY and INDEPENDENTLY (no
#       shared npi_scaler), so raising HCW deaths via low PPE no longer drags ETU
#       -- and overall epidemic size -- down.
#     * hcw_risk_scalar (bounded prior) is the HCW-exposure lever.
#
#   FITTED SUMMARIES: takeoff, log(n_deaths), log(n_hcw_deaths), hcw_fraction,
#                     d_p05_p95, log(peak_height)   (any subset; see SUMMARY_STATS)
#     * takeoff is the fraction of trajectories that took off, scored on the SAME
#       TAKEOFF_DEATH_THRESHOLD (>= 100 deaths) that conditions the means -- so
#       "took off" means one thing throughout, and the fraction is what conditions
#       the means. Observed target = 1.0 (the real outbreak did take off);
#     * logs put the heavy-tailed counts on a relative-error scale;
#     * hcw_fraction up-weights the HCW share (deliberately redundant with the
#       two count logs);
#     * d_p05_p95 (5-95% span of death dates) is a low-noise, tail-robust timing
#       summary that replaces the old "duration".
#
# >>> DURATION INPUT MUST CHANGE <<<  The observed timing target is NO LONGER the
#     365-day "duration". d_p05_p95 is the day span between the 5th and 95th
#     percentile of DEATH DATES. RECOMPUTE the observed value from the weekly
#     death series with observed_d_p05_p95() (see section 1 + the commented call
#     in section 3). The value below is a PLACEHOLDER.
#
# >>> IDENTIFIABILITY <<<  PPE and hcw_risk_scalar both act on the HCW share and
#     form a prior-bounded ridge (no time-resolved HCW data). Report them JOINTLY;
#     PPE is substantially prior-informed (informative normal prior below). The
#     section-9 pairs() plot shows the ridge.
#
# RUN_PROFILE selects smoke / check / production settings (section 1).
# =============================================================================


# -----------------------------------------------------------------------------
# 1. CONFIGURATION
# -----------------------------------------------------------------------------
ANALYSIS_DIR   <- here::here("analyses", "02_ABC_model_fits_Final")
FUNCTIONS_DIR  <- here::here("functions")
SETUP_PATH     <- file.path(FUNCTIONS_DIR, "setup_model_parameters.R")
COMMON_PATH    <- file.path(FUNCTIONS_DIR, "abc_calibration_functions_common.R")
FUNCTIONS_PATH <- file.path(FUNCTIONS_DIR, "abc_calibration_functions_decoupled.R")  # this scheme
R0_PATH        <- file.path(FUNCTIONS_DIR, "calculate_model_approx_r0.R")
SCENARIO_CSV   <- here::here("data-processed", "final_six_scenario_values_original_approach.csv")
SCENARIO_ID    <- "Worst_WestAfrica"

# ---- RUN PROFILE: "smoke" | "quickcheck" | "check" | "production" -----------
#   smoke      : minutes; confirms the pipeline runs end-to-end (NOT a real fit).
#   quickcheck : <1.5 h on ~118-120 cores; validates a new fiber build end-to-end
#                AND re-times it (nb_simul = 118 = one full wave; WA ~1 h, DRC mins).
#   check      : ~1 hr-ish; a rough posterior to sanity-check shapes/targets.
#   production : the real fit (n_reps=40 from the noise check; stop on plateau).
RUN_PROFILE <- "quickcheck"
.PROFILES <- list(
  smoke      = list(n_reps =  5L, nb_simul =  60L, tolerance_target = 5.00, n_traj =  20L),
  quickcheck = list(n_reps =  8L, nb_simul = 118L, tolerance_target = 1.20, n_traj =  30L),
  check      = list(n_reps = 30L, nb_simul = 472L, tolerance_target = 1, n_traj = 200L),
  production = list(n_reps = 40L, nb_simul = 944L, tolerance_target = 0.35, n_traj = 200L)
)
stopifnot(RUN_PROFILE %in% names(.PROFILES))
.prof <- .PROFILES[[RUN_PROFILE]]

# ---- WHICH PARAMETERS TO FIT ------------------------------------------------
# Any subset of c("R0","prop_funeral","etu_efficacy","ppe_efficacy","hcw_risk_scalar").
FIT_PARAMS <- c("R0", "prop_funeral", "etu_efficacy", "ppe_efficacy", "hcw_risk_scalar")

# Priors keyed BY NAME (only FIT_PARAMS entries are used; order here is irrelevant).
PRIORS_NAMED <- list(
  R0              = c("unif",   1.15, 1.65),
  prop_funeral    = c("unif",   0.10, 0.40),
  etu_efficacy    = c("unif",   0.60, 0.95),   # wide: the size/peak lever
  ppe_efficacy    = c("unif", 0.30, 0.90),   # INFORMATIVE; clamp to [0,1] in build. Confirm centre/sd vs literature.
  hcw_risk_scalar = c("unif",   1.00, 3.00)    # bounded: prob_hcw 0.25..0.75 (can't run away)
)

# Fixed values for any parameter NOT in FIT_PARAMS (fall back to DECOUPLED_PARAM_DEFAULTS).
FIXED_PARAMS <- list(
  R0 = 1.45, prop_funeral = 0.30, etu_efficacy = 0.85, ppe_efficacy = 0.65, hcw_risk_scalar = 1.5
)

# ---- WHICH SUMMARIES TO FIT -------------------------------------------------
# Any subset of DECOUPLED_AVAILABLE_SUMMARIES; comment a line out of BOTH this
# vector and OBSERVED_NAMED to drop a summary from the fit.
SUMMARY_STATS <- c("takeoff", "log_n_deaths", "log_n_hcw_deaths", "hcw_fraction", "log_peak_height", "d_p05_p95")

# Observed targets, ON THE FITTED SCALE (log the counts), keyed BY NAME.
#   raw WA targets: n_deaths = 11325, n_hcw_deaths = 513, peak_height = 599.
OBSERVED_NAMED <- c(
  takeoff          = 1.0,           # the real outbreak took off (>= TAKEOFF_DEATH_THRESHOLD deaths)
  log_n_deaths     = log(11325),
  log_n_hcw_deaths = log(513),
  hcw_fraction     = 513 / 11325,   # = 0.0453
  d_p05_p95        = 274, # linear interpolation from ## info from here: https://en.wikipedia.org/wiki/West_African_Ebola_virus_epidemic_timeline_of_reported_cases_and_deaths 
                          # gives the period 23rd July 2014 - 23 April 2015
  log_peak_height  = log(599) ## info from here: https://en.wikipedia.org/wiki/West_African_Ebola_virus_epidemic_timeline_of_reported_cases_and_deaths 
)

HCW_BASE_PROB    <- 0.25            # prob_hcw_cond_*_hospital = min(base * hcw_risk_scalar, 1)
PEAK_BIN_WIDTH   <- 7L
PEAK_TIME_ORIGIN <- "first_death"

# WA epidemic is ~28k cases; keep check_final_size well above it so summaries are
# not truncated (correctness AND it bounds runaway draws' cost).
MODEL_OVERRIDES <- list(check_final_size = 40000)

ABC_OUTPUT_BASE   <- ANALYSIS_DIR
ABC_OUTPUT_LABEL  <- "Decoupled"
FINAL_OUTPUTS_DIR <- here::here("outputs", "02_ABC_model_fits_Final")
if (!dir.exists(FINAL_OUTPUTS_DIR)) dir.create(FINAL_OUTPUTS_DIR, recursive = TRUE, showWarnings = FALSE)

ABC_SETTINGS <- list(
  method           = "Delmoral",
  nb_simul         = .prof$nb_simul,
  alpha            = 0.5,
  tolerance_target = .prof$tolerance_target,  # don't chase below ~0.25-0.30 (noise floor); stop on plateau
  M                = 1,
  use_seed         = TRUE,
  verbose          = TRUE
)

N_REPS                  <- .prof$n_reps
SEEDING_CASES           <- 25L
TAKEOFF_DEATH_THRESHOLD <- 100L   # >= K deaths = "took off": conditions the means AND the takeoff fraction
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
source(FUNCTIONS_PATH)   # the DECOUPLED scheme
source(R0_PATH)

check_model_function_version()


# -----------------------------------------------------------------------------
# 3. VALIDATE CHOICES + ASSEMBLE ORDER-SAFE ABC INPUTS
# -----------------------------------------------------------------------------
prep <- prepare_decoupled_run(
  fit_params     = FIT_PARAMS,
  priors_named   = PRIORS_NAMED,
  fixed_params   = FIXED_PARAMS,
  summary_stats  = SUMMARY_STATS,
  observed_named = OBSERVED_NAMED
)
priors             <- prep$priors
observed_summaries <- prep$observed_summaries
fixed_values       <- prep$fixed_values
PARAM_NAMES        <- prep$fit_params

ABC_CONFIG <- list(
  takeoff_death_threshold = TAKEOFF_DEATH_THRESHOLD,
  n_reps                  = N_REPS,
  seeding_cases           = SEEDING_CASES,
  setup_R0_n              = SETUP_R0_N,
  setup_R0_seed           = SETUP_R0_SEED,
  hcw_base_prob           = HCW_BASE_PROB,
  fit_params              = prep$fit_params,
  fixed_values            = fixed_values,
  summary_stats           = prep$summary_stats,
  peak_bin_width          = PEAK_BIN_WIDTH,
  peak_time_origin        = PEAK_TIME_ORIGIN
)

cat(sprintf("RUN_PROFILE = %s  |  n_reps = %d, nb_simul = %d, tol = %.2f\n",
            RUN_PROFILE, N_REPS, ABC_SETTINGS$nb_simul, ABC_SETTINGS$tolerance_target))
cat("Fitting parameters : ", paste(PARAM_NAMES, collapse = ", "), "\n", sep = "")
cat("Fitting summaries  : ", paste(prep$summary_stats, collapse = ", "), "\n", sep = "")
cat("Observed targets   : ",
    paste(sprintf("%s=%.4g", names(observed_summaries), observed_summaries), collapse = ", "),
    "\n", sep = "")


# -----------------------------------------------------------------------------
# 4. BUILD BASE + TIME-VARYING ARGS; COMPUTE EFFICACY-INDEPENDENT R0 INVARIANTS
# -----------------------------------------------------------------------------
scenario_matrix <- read_scenario_matrix(SCENARIO_CSV)
mp <- make_model_parameters(scenario_id = SCENARIO_ID, scenario_matrix = scenario_matrix,
                            overrides = MODEL_OVERRIDES)
base_args     <- mp$base_args
tv_args_model <- mp$tv_args
R0_invariants <- compute_R0_invariants(args = mp$args, n = ABC_CONFIG$setup_R0_n,
                                       seed = ABC_CONFIG$setup_R0_seed)


# -----------------------------------------------------------------------------
# 5. PRIOR PREDICTIVE CHECK (optional; sizes the d_p05_p95 target + checks regime)
# -----------------------------------------------------------------------------
# set.seed(1)
# pp <- prior_predictive_check_decoupled(
#   n_draws = 50, prior_list = priors, fit_params = PARAM_NAMES, fixed_values = fixed_values,
#   base = base_args, tv = tv_args_model, invariants = R0_invariants,
#   general_hospital_quarantine_efficacy = DEFAULT_SCALAR_INPUTS$general_hospital_quarantine_efficacy,
#   safe_funeral_efficacy = DEFAULT_SCALAR_INPUTS$safe_funeral_efficacy,
#   hcw_base_prob = HCW_BASE_PROB, n_replicates = 5, seeding_cases = ABC_CONFIG$seeding_cases,
#   takeoff_death_threshold = ABC_CONFIG$takeoff_death_threshold,
#   summary_stats = prep$summary_stats, bin_width = PEAK_BIN_WIDTH, time_origin = PEAK_TIME_ORIGIN)
# print(summary(pp))


# -----------------------------------------------------------------------------
# 6. PER-RUN OUTPUT DIRECTORY + WORKER CONFIG + PRE-RUN PROVENANCE
# -----------------------------------------------------------------------------
RUN_TAG <- sprintf("%s_NP%d_NS%d_NBREPS_%d_NBSIMUL_%d",
                   RUN_PROFILE, length(PARAM_NAMES), length(prep$summary_stats),
                   N_REPS, ABC_SETTINGS$nb_simul)

ABC_OUTPUT_DIR <- make_abc_output_dir(base_dir = ABC_OUTPUT_BASE, scenario_id = SCENARIO_ID,
                                      label = ABC_OUTPUT_LABEL, suffix = RUN_TAG)
message("ABC outputs will be written to: ", ABC_OUTPUT_DIR)

save_abc_config(list(
  setup_path = SETUP_PATH, functions_path = FUNCTIONS_PATH, r0_path = R0_PATH,
  scenario_csv = SCENARIO_CSV, scenario_id = SCENARIO_ID,
  abc_config = ABC_CONFIG, model_overrides = MODEL_OVERRIDES
))

# Build + write the run metadata to ABC_OUTPUT_DIR *before* ABC_sequential() starts.
# A run interrupted mid-fit leaves only the raw output_step*/tolerance_step* files
# (no result.rds); writing the metadata up front means those partial runs can still
# be traced back to the exact parameter combination, priors, fixed values and
# observed targets that produced them. start_time / result_stamp / result_filename
# are fixed here so the pre-run sidecar matches the post-run result.rds copy exactly.
start_time      <- Sys.time()
result_stamp    <- format(start_time, "%Y%m%d_%H%M%S")
result_filename <- paste0("fiber_ABC_SMC_", SCENARIO_ID, "_", ABC_OUTPUT_LABEL,
                          "_", result_stamp, "_", RUN_TAG, ".rds")

run_metadata <- make_decoupled_run_metadata(
  scenario_id = SCENARIO_ID, fit_params = PARAM_NAMES, priors = priors,
  fixed_values = fixed_values, summary_stats = prep$summary_stats,
  observed_summaries = observed_summaries, hcw_base_prob = HCW_BASE_PROB,
  general_hospital_quarantine_efficacy = DEFAULT_SCALAR_INPUTS$general_hospital_quarantine_efficacy,
  safe_funeral_efficacy = DEFAULT_SCALAR_INPUTS$safe_funeral_efficacy,
  peak_settings = list(bin_width = PEAK_BIN_WIDTH, time_origin = PEAK_TIME_ORIGIN),
  extra = list(result_filename = result_filename, run_profile = RUN_PROFILE,
               abc_settings = ABC_SETTINGS, abc_output_dir = ABC_OUTPUT_DIR)
)
# ABC_OUTPUT_DIR's name already encodes scenario/label/timestamp/tag, so the .rds +
# metadata inside it use a SHORT stem ("result") to stay under the Windows 260-char
# MAX_PATH limit (a deep timestamped dir + the long descriptive filename can exceed
# it -- which fails the write even though the short output_step* files write fine).
# The descriptive name is kept for the shallow FINAL_OUTPUTS_DIR copy.
write_decoupled_run_metadata(run_metadata, file.path(ABC_OUTPUT_DIR, "result.rds"))


# -----------------------------------------------------------------------------
# 7. RUN ABC_SEQUENTIAL (Del Moral et al. 2012 adaptive SMC)
# -----------------------------------------------------------------------------
result <- with_abc_output_dir(
  ABC_OUTPUT_DIR,
  ABC_sequential(
    method              = ABC_SETTINGS$method,
    model               = fiber_abc_model_parallel_decoupled,
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
print(Sys.time() - start_time)

# Metadata sidecars were written to ABC_OUTPUT_DIR pre-run (section 6); here we
# attach the same metadata to the completed object and persist the .rds plus the
# descriptive-name copy (and its sidecar) to FINAL_OUTPUTS_DIR.
result <- attach_decoupled_run_metadata(result, run_metadata)

saveRDS(result, file = file.path(ABC_OUTPUT_DIR, "result.rds"))
saveRDS(result, file = file.path(FINAL_OUTPUTS_DIR, result_filename))
write_decoupled_run_metadata(run_metadata, file.path(FINAL_OUTPUTS_DIR, result_filename))


# -----------------------------------------------------------------------------
# 8. POSTERIOR INSPECTION (+ implied prob_hcw, + ridge diagnostic)
# -----------------------------------------------------------------------------
posterior <- as.data.frame(result$param)
colnames(posterior) <- PARAM_NAMES
print(apply(posterior, 2, quantile, probs = c(0.025, 0.5, 0.975)))

draw_values <- function(name) {
  if (name %in% PARAM_NAMES) posterior[[name]] else rep(fixed_values[[name]], nrow(posterior))
}
implied_prob_hcw <- pmin(HCW_BASE_PROB * draw_values("hcw_risk_scalar"), 1.0)
cat("\nImplied prob_hcw_hospital (2.5/50/97.5%): ",
    paste(round(quantile(implied_prob_hcw, c(0.025, 0.5, 0.975)), 3), collapse = " / "), "\n")

par(mfrow = c(1, length(PARAM_NAMES)))
for (j in seq_len(ncol(posterior))) {
  hist(posterior[, j], breaks = 12, main = colnames(posterior)[j], xlab = colnames(posterior)[j])
  abline(v = quantile(posterior[, j], c(0.025, 0.5, 0.975)), lty = c(2, 1, 2), col = "red")
}
par(mfrow = c(1, 1))

# Ridge check: a tight anti-diagonal in (ppe_efficacy, hcw_risk_scalar) = the
# expected PPE<->exposure ridge. Report those two JOINTLY.
if (all(c("ppe_efficacy", "hcw_risk_scalar") %in% PARAM_NAMES)) {
  plot(posterior$ppe_efficacy, posterior$hcw_risk_scalar,
       xlab = "ppe_efficacy", ylab = "hcw_risk_scalar", pch = 16,
       col = adjustcolor("steelblue", 0.4), main = "PPE <-> hcw_risk ridge")
}


# -----------------------------------------------------------------------------
# 9. PROGRESS / RECONSTRUCTION FROM DISK
# -----------------------------------------------------------------------------
abc_progress(ABC_OUTPUT_DIR, tolerance_target = ABC_SETTINGS$tolerance_target,
             param_names = PARAM_NAMES, stat_names = prep$summary_stats)
print(abc_compare_steps(ABC_OUTPUT_DIR, param_names = PARAM_NAMES, stat_names = prep$summary_stats))
result <- reconstruct_abc_result(ABC_OUTPUT_DIR, param_names = PARAM_NAMES,
                                 stat_names = prep$summary_stats)


# -----------------------------------------------------------------------------
# 10. POSTERIOR PREDICTIVE CHECKS (fitted summaries; targets on the FITTED scale)
# -----------------------------------------------------------------------------
sim_stats <- as.data.frame(result$stats)
colnames(sim_stats) <- prep$summary_stats

set.seed(1)
idx <- sample(seq_len(nrow(sim_stats)), size = 10000, replace = TRUE, prob = result$weights)
sim_stats_post <- sim_stats[idx, ]

n_stat <- length(observed_summaries)
par(mfrow = c(ceiling(n_stat / 2), 2), mar = c(4, 4, 3, 1))
for (s in names(observed_summaries)) {
  x <- sim_stats_post[[s]]
  hist(x, breaks = 12, main = paste0("Posterior predictive: ", s), xlab = s,
       col = adjustcolor("steelblue", 0.6), border = "white")
  abline(v = quantile(x, c(0.025, 0.5, 0.975)), col = "darkblue", lty = c(2, 1, 2), lwd = c(1, 2, 1))
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
library(ggplot2)
source(file.path(FUNCTIONS_DIR, "simulation_helpers.R"))   # bin_counts()

N_TRAJ    <- .prof$n_traj
BIN_WIDTH <- PEAK_BIN_WIDTH
TRAJ_SEED <- 100L

set.seed(TRAJ_SEED)
traj_idx <- sample(seq_len(nrow(posterior)), size = N_TRAJ, replace = TRUE, prob = result$weights)

plan(multisession, workers = min(N_CLUSTER, future::availableCores()))
traj_runs <- future_lapply(seq_along(traj_idx), function(k) {
  i    <- traj_idx[k]
  full <- assemble_decoupled_theta(as.numeric(posterior[i, PARAM_NAMES]), PARAM_NAMES, fixed_values)
  a <- build_abc_model_args_decoupled(
    R0 = full$R0, prop_funeral = full$prop_funeral, etu_efficacy = full$etu_efficacy,
    ppe_efficacy = full$ppe_efficacy, hcw_risk_scalar = full$hcw_risk_scalar,
    base = base_args, tv = tv_args_model, invariants = R0_invariants,
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

max_day <- max(unlist(lapply(traj_runs, function(r) c(r$case_days, r$death_days))), 0)
n_bins  <- max(1L, ceiling((max_day + 1L) / BIN_WIDTH))
week    <- ((seq_len(n_bins) - 1L) + 0.5) * BIN_WIDTH

traj_long <- do.call(rbind, lapply(seq_along(traj_runs), function(k) {
  r <- traj_runs[[k]]
  rbind(data.frame(draw = k, week = week, metric = "Cases",
                   incidence = bin_counts(r$case_days,  BIN_WIDTH, n_bins)),
        data.frame(draw = k, week = week, metric = "Deaths",
                   incidence = bin_counts(r$death_days, BIN_WIDTH, n_bins)))
}))

print(ggplot(traj_long, aes(week, incidence, group = draw)) +
        geom_line(alpha = 0.15, colour = "steelblue", linewidth = 0.3) +
        facet_wrap(~ metric, scales = "free_y") +
        labs(x = "Time since outbreak start (days)", y = "Count per week",
             title = sprintf("Posterior trajectories: %s (%s)", SCENARIO_ID, RUN_PROFILE)) +
        theme_bw())

band <- do.call(rbind, lapply(split(traj_long, list(traj_long$metric, traj_long$week), drop = TRUE),
  function(d) data.frame(metric = d$metric[1], week = d$week[1],
                         lo = quantile(d$incidence, 0.025, names = FALSE),
                         med = quantile(d$incidence, 0.500, names = FALSE),
                         hi = quantile(d$incidence, 0.975, names = FALSE))))
print(ggplot(band, aes(week, med)) +
        geom_ribbon(aes(ymin = lo, ymax = hi), fill = "steelblue", alpha = 0.3) +
        geom_line(colour = "steelblue", linewidth = 0.8) +
        facet_wrap(~ metric, scales = "free_y") +
        labs(x = "Time since outbreak start (days)", y = "Count per week",
             title = sprintf("Posterior median + 95%% CrI: %s", SCENARIO_ID)) +
        theme_bw())

# -----------------------------------------------------------------------------
# 12. OUTBREAK DURATION FROM SIMULATED TRAJECTORIES
# -----------------------------------------------------------------------------
# Per-trajectory outbreak duration, recomputed from the SAME posterior-draw
# trajectories as section 11 (`traj_runs`) -- no new model runs. Measured from
# DEATH dates (the scale this scheme's timing summary is defined on; swap
# r$death_days -> r$case_days for a case-based span). Two complementary measures:
#   * d_p05_p95 : day span between the 5th and 95th percentile of death dates --
#                 the tail-robust timing summary this scheme fits. Directly
#                 comparable to the observed target (~378 d; Kivu 4 Oct 2018 -
#                 17 Oct 2019). Recompute the observed value via observed_d_p05_p95().
#   * dur_full  : first-death to last-death span -- the classic, tail-sensitive
#                 "duration", reported for context (no single observed target).
# Trajectories with no deaths give NA; a single death gives a span of 0.
dur_post <- do.call(rbind, lapply(traj_runs, function(r) {
  d <- r$death_days[is.finite(r$death_days)]
  if (length(d) == 0L) {
    return(data.frame(n_deaths = 0L, dur_full = NA_real_, d_p05_p95 = NA_real_))
  }
  qs <- quantile(d, c(0.05, 0.95), names = FALSE)
  data.frame(n_deaths  = length(d),
             dur_full  = max(d) - min(d),
             d_p05_p95 = qs[2] - qs[1])
}))

obs_dur <- c(
  dur_full  = NA_real_,    # no single observed "duration" target in this scheme
  d_p05_p95 = if ("d_p05_p95" %in% names(observed_summaries)) {
    observed_summaries[["d_p05_p95"]]
  } else 378    # placeholder: Kivu 4 Oct 2018 - 17 Oct 2019
)

cat("\nOutbreak duration posterior predictive (median / 95% CrI, days):\n")
for (s in c("dur_full", "d_p05_p95")) {
  qs <- quantile(dur_post[[s]], c(0.025, 0.5, 0.975), na.rm = TRUE)
  obs_txt <- if (is.finite(obs_dur[s])) sprintf("   (observed %.1f)", obs_dur[s]) else ""
  cat(sprintf("  %-11s %.1f / [%.1f, %.1f]%s\n",
              paste0(s, ":"), qs[2], qs[1], qs[3], obs_txt))
}

par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
for (s in c("dur_full", "d_p05_p95")) {
  x  <- dur_post[[s]][is.finite(dur_post[[s]])]
  qs <- quantile(x, probs = c(0.025, 0.5, 0.975))
  hist(x, breaks = 12, main = paste0("Trajectory ", s), xlab = paste0(s, " (days)"),
       col = adjustcolor("steelblue", 0.6), border = "white")
  abline(v = qs, col = "darkblue", lty = c(2, 1, 2), lwd = c(1, 2, 1))
  if (is.finite(obs_dur[s])) abline(v = obs_dur[s], col = "red", lwd = 2.5)
}
par(mfrow = c(1, 1))

