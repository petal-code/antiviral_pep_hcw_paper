# DRC_posterior_trajectory_bands_decoupled.R
# =============================================================================
# Posterior trajectory "spaghetti + bands" plot for the DECOUPLED-efficacy DRC
# fit (Middle_DRC_ConflictSmoothed_PlusPlus).
#
# What it does, exactly as requested:
#   1. Reads the ABC-SMC result.rds from the decoupled DRC run.
#   2. Draws N_POST = 250 posterior parameter samples (weighted, with
#      replacement, by the ABC importance weights).
#   3. For EACH posterior sample, runs N_REPS = 10 stochastic replicates of the
#      fiber branching-process model and bins each replicate into a weekly
#      incidence curve, then takes the element-wise MEDIAN across the 10
#      replicates -> ONE "median trajectory" per posterior sample (250 of them).
#   4. Plots:
#        * every per-sample median trajectory          (thin faint lines)
#        * the MEDIAN of the 250 median trajectories    (bold central line)
#        * the inter-quartile range  (25-75%)           (darker band)
#        * the 95% credible interval (2.5-97.5%)        (lighter band)
#      ...computed point-wise (per weekly bin) ACROSS the 250 median
#      trajectories. Curves are shown for both weekly Cases and weekly Deaths
#      (Deaths is the fitted quantity; Cases is included for context), faceted.
#
# Parallelism: stochastic replicates are run with future / future.apply using a
# `multisession` plan, which is the Windows-safe (PSOCK) backend -- no forking.
# Each replicate is made reproducible by an explicit per-(sample, replicate)
# seed, so results are identical regardless of the number of workers.
#
# This mirrors the simulation pipeline in section 11 of
# DRC_run_abc_calibration_decoupled.R (same args-builder, same trajectory
# extraction, same weekly binning); it only changes the *statistical* layout
# (10 replicates per draw -> per-draw median -> cross-draw median/IQR/95% CrI).
#
# Run from anywhere inside the project (paths resolved with here::here()):
#   Rscript analyses/02_ABC_model_fits_Final/DRC_posterior_trajectory_bands_decoupled.R
# or source() it interactively.
# =============================================================================


# -----------------------------------------------------------------------------
# 1. CONFIGURATION
# -----------------------------------------------------------------------------
FUNCTIONS_DIR <- here::here("functions")
ANALYSIS_DIR  <- here::here("analyses", "02_ABC_model_fits_Final")

# The decoupled DRC ABC-SMC result. Override RDS_PATH if you want a different run.
RDS_PATH <- file.path(
  ANALYSIS_DIR, "abc_outputs",
  "Middle_DRC_ConflictSmoothed_PlusPlus_20260607_215621_Decoupled_check_NP5_NS4_NBREPS_30_NBSIMUL_590",
  "result.rds"
)

# Scenario inputs (must match the calibration run that produced RDS_PATH).
SCENARIO_CSV     <- here::here("data-processed", "final_six_scenario_values_original_approach.csv")
SCENARIO_ID_DFLT <- "Middle_DRC_ConflictSmoothed_PlusPlus"

# ---- Sampling / replication ----
N_POST    <- 250L     # posterior parameter samples
N_REPS    <- 10L      # stochastic replicates per posterior sample
BIN_WIDTH <- 7L       # days per incidence bin (weekly); overridden by run metadata if present
METRICS   <- c("Deaths", "Cases")   # facet order (Deaths first = the fitted quantity)

# ---- Seeds (reproducibility) ----
SAMPLE_SEED <- 2026L  # selects which 250 posterior particles are drawn
BASE_SEED   <- 10000L # per-(sample, replicate) model seeds start here

# ---- Model-build constants (must match the calibration run) ----
SEEDING_CASES           <- 25L
SETUP_R0_N              <- 100000L   # n for compute_R0_invariants() in the fit
SETUP_R0_SEED           <- 42L       # seed   for compute_R0_invariants() in the fit
MODEL_OVERRIDES         <- list(check_final_size = 10000)
TAKEOFF_DEATH_THRESHOLD <- 100L      # used only for the take-off diagnostic below

# If TRUE, the per-sample median trajectory is taken only over replicates that
# "took off" (>= TAKEOFF_DEATH_THRESHOLD deaths); samples with no take-off are
# dropped. Default FALSE = median over all N_REPS replicates, exactly as asked.
CONDITION_ON_TAKEOFF <- TRUE

# ---- Parallel workers (Windows-safe multisession / PSOCK) ----
# NULL = auto. Mirrors the calibration script's heuristic; cap as you like.
N_WORKERS_OVERRIDE <- NULL

# ---- Outputs ----
FIG_DIR <- here::here("figures")
OUT_DIR <- here::here("outputs", "02_ABC_model_fits_Final")
FIG_STEM <- "DRC_decoupled_posterior_trajectory_bands"

# ---- Plot aesthetics ----
BAND_FILL       <- "#3182BD"
CI95_ALPHA      <- 0.18    # lighter band  (95% CrI)
IQR_ALPHA       <- 0.40    # darker band   (IQR)
SPAGHETTI_COL   <- "grey25"
SPAGHETTI_ALPHA <- 0.08
SPAGHETTI_LWD   <- 0.20
MEDIAN_COL      <- "#08306B"
MEDIAN_LWD      <- 1.10


# -----------------------------------------------------------------------------
# 2. LIBRARIES + SOURCES  (order required by the decoupled scheme)
# -----------------------------------------------------------------------------
library(fiber)
library(future)
library(future.apply)
library(ggplot2)

source(file.path(FUNCTIONS_DIR, "setup_model_parameters.R"))              # make_model_parameters(), read_scenario_matrix(), DEFAULT_SCALAR_INPUTS
source(file.path(FUNCTIONS_DIR, "abc_calibration_functions_common.R"))    # abc_summarise(), check_model_function_version()
source(file.path(FUNCTIONS_DIR, "abc_calibration_functions_decoupled.R")) # assemble/build args, DECOUPLED_PARAM_NAMES/DEFAULTS
source(file.path(FUNCTIONS_DIR, "calculate_model_approx_r0.R"))           # compute_R0_invariants(), D/F_from_invariants(), solve_offspring_means()
source(file.path(FUNCTIONS_DIR, "simulation_helpers.R"))                  # bin_counts()

# Soft version check (warn, don't stop, if the helper is unavailable/mismatched).
tryCatch(check_model_function_version(),
         error = function(e) warning("check_model_function_version() failed: ", conditionMessage(e)))


# -----------------------------------------------------------------------------
# 3. LOAD POSTERIOR + RUN METADATA
# -----------------------------------------------------------------------------
if (!file.exists(RDS_PATH)) {
  stop("Could not find the decoupled DRC result at:\n  ", RDS_PATH,
       "\nEdit RDS_PATH at the top of this script to point at your result.rds.", call. = FALSE)
}
result <- readRDS(RDS_PATH)
meta   <- result$run_metadata  # NULL-safe accessors below

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

PARAM_NAMES  <- meta$fit_params    %||% DECOUPLED_PARAM_NAMES
fixed_values <- meta$fixed_values  %||% DECOUPLED_PARAM_DEFAULTS
SCENARIO_ID  <- meta$scenario_id   %||% SCENARIO_ID_DFLT
HCW_BASE_PROB <- meta$hcw_base_prob %||% 0.25
GEN_HOSP_EFF <- meta$fixed_efficacies$general_hospital_quarantine_efficacy %||%
                  DEFAULT_SCALAR_INPUTS$general_hospital_quarantine_efficacy
SAFE_FUN_EFF <- meta$fixed_efficacies$safe_funeral_efficacy %||%
                  DEFAULT_SCALAR_INPUTS$safe_funeral_efficacy
BIN_WIDTH    <- as.integer(meta$peak_settings$bin_width %||% BIN_WIDTH)

posterior <- as.data.frame(result$param)
colnames(posterior) <- PARAM_NAMES
weights <- result$weights
if (is.null(weights) || length(weights) != nrow(posterior)) {
  warning("No usable result$weights found; sampling posterior particles uniformly.")
  weights <- rep(1, nrow(posterior))
}

cat(sprintf("Loaded %d posterior particles from:\n  %s\n", nrow(posterior), RDS_PATH))
cat("Scenario        : ", SCENARIO_ID, "\n", sep = "")
cat("Fitted params   : ", paste(PARAM_NAMES, collapse = ", "), "\n", sep = "")
cat(sprintf("Sampling %d draws x %d replicates = %d simulations; bin width = %d days.\n",
            N_POST, N_REPS, N_POST * N_REPS, BIN_WIDTH))


# -----------------------------------------------------------------------------
# 4. BUILD BASE / TIME-VARYING ARGS + R0 INVARIANTS (must match the fit)
# -----------------------------------------------------------------------------
scenario_matrix <- read_scenario_matrix(SCENARIO_CSV)
mp <- make_model_parameters(scenario_id = SCENARIO_ID, scenario_matrix = scenario_matrix,
                            overrides = MODEL_OVERRIDES)
base_args     <- mp$base_args
tv_args_model <- mp$tv_args
R0_invariants <- compute_R0_invariants(args = mp$args, n = SETUP_R0_N, seed = SETUP_R0_SEED)


# -----------------------------------------------------------------------------
# 5. DRAW 250 POSTERIOR SAMPLES + BUILD THE (sample, replicate) JOB LIST
# -----------------------------------------------------------------------------
set.seed(SAMPLE_SEED)
draw_idx <- sample(seq_len(nrow(posterior)), size = N_POST, replace = TRUE, prob = weights)

# One job per (posterior sample k, replicate r) with a unique, reproducible seed.
jobs <- vector("list", N_POST * N_REPS)
jj <- 1L
for (k in seq_len(N_POST)) {
  for (r in seq_len(N_REPS)) {
    jobs[[jj]] <- list(draw     = k,                 # 1..N_POST
                       draw_row = draw_idx[k],       # row in `posterior`
                       rep      = r,                 # 1..N_REPS
                       seed     = BASE_SEED + (k - 1L) * N_REPS + r)
    jj <- jj + 1L
  }
}


# -----------------------------------------------------------------------------
# 6. RUN ALL REPLICATES IN PARALLEL (future / multisession; Windows-safe)
# -----------------------------------------------------------------------------
# Worker builds the full fiber arg list from the posterior draw (same pipeline as
# section 11 of the calibration script), sets the explicit per-replicate seed,
# runs ONE branching-process replicate, and returns the floored absolute day of
# every infection (cases) and every death. fiber is loaded on each worker via
# future.packages; all other objects/functions are exported automatically.
run_traj_job <- function(job) {
  full <- assemble_decoupled_theta(as.numeric(posterior[job$draw_row, PARAM_NAMES]),
                                   PARAM_NAMES, fixed_values)
  a <- build_abc_model_args_decoupled(
    R0 = full$R0, prop_funeral = full$prop_funeral, etu_efficacy = full$etu_efficacy,
    ppe_efficacy = full$ppe_efficacy, hcw_risk_scalar = full$hcw_risk_scalar,
    base = base_args, tv = tv_args_model, invariants = R0_invariants,
    general_hospital_quarantine_efficacy = GEN_HOSP_EFF,
    safe_funeral_efficacy = SAFE_FUN_EFF,
    hcw_base_prob = HCW_BASE_PROB, seeding_cases = SEEDING_CASES)
  a$seed <- job$seed

  tdf <- do.call(branching_process_main, a)$tdf
  tdf <- tdf[!is.na(tdf$time_infection_absolute), , drop = FALSE]
  died <- !is.na(tdf$outcome) & tdf$outcome

  list(draw       = job$draw,
       rep        = job$rep,
       case_days  = as.integer(floor(tdf$time_infection_absolute)),
       death_days = as.integer(floor(tdf$time_outcome_absolute[died])))
}

N_CLUSTER <- if (parallel::detectCores() > 120) {
  min(120L, parallel::detectCores() - 10L)
} else {
  min(12L, parallel::detectCores() - 4L)
}
N_WORKERS <- N_WORKERS_OVERRIDE %||% max(1L, min(N_CLUSTER, future::availableCores()))

cat(sprintf("Running %d simulations on %d workers (multisession)...\n",
            length(jobs), N_WORKERS))
start_time <- Sys.time()

plan(multisession, workers = N_WORKERS)
results <- future_lapply(jobs, run_traj_job,
                         future.packages = "fiber", future.seed = TRUE)
plan(sequential)

cat("Simulation time: "); print(Sys.time() - start_time)


# -----------------------------------------------------------------------------
# 7. BIN TO WEEKLY INCIDENCE -> PER-SAMPLE MEDIAN TRAJECTORIES
# -----------------------------------------------------------------------------
# Common weekly grid across ALL replicates (absolute model time, from day 0).
max_day <- max(c(unlist(lapply(results, function(r) c(r$case_days, r$death_days))), 0L))
n_bins  <- max(1L, ceiling((max_day + 1L) / BIN_WIDTH))
week    <- ((seq_len(n_bins) - 1L) + 0.5) * BIN_WIDTH

# inc[[metric]] is a [N_POST x N_REPS x n_bins] array of weekly counts.
inc <- list(
  Cases  = array(0L, dim = c(N_POST, N_REPS, n_bins)),
  Deaths = array(0L, dim = c(N_POST, N_REPS, n_bins))
)
# took[draw, rep] = did that replicate reach the take-off death threshold?
took <- matrix(FALSE, N_POST, N_REPS)
for (r in results) {
  inc$Cases [r$draw, r$rep, ] <- bin_counts(r$case_days,  BIN_WIDTH, n_bins)
  inc$Deaths[r$draw, r$rep, ] <- bin_counts(r$death_days, BIN_WIDTH, n_bins)
  took[r$draw, r$rep] <- length(r$death_days) >= TAKEOFF_DEATH_THRESHOLD
}

# Take-off diagnostic (informational; does not change the plot unless
# CONDITION_ON_TAKEOFF = TRUE).
takeoff_frac_per_draw <- rowMeans(took)
cat(sprintf("Take-off (>= %d deaths): %.1f%% of all replicates; per-sample take-off fraction median = %.2f.\n",
            TAKEOFF_DEATH_THRESHOLD, 100 * mean(took), median(takeoff_frac_per_draw)))

# Per-sample median trajectory: element-wise median across the N_REPS replicates.
median_traj <- lapply(METRICS, function(m) {
  arr <- inc[[m]]                          # [N_POST, N_REPS, n_bins]
  if (CONDITION_ON_TAKEOFF) arr[!took] <- NA_integer_   # recycles over the bin dim
  mt <- apply(arr, c(1, 3), median, na.rm = CONDITION_ON_TAKEOFF)  # [N_POST, n_bins]
  mt[!is.finite(mt)] <- NA_real_           # samples with no take-off -> NA rows
  mt
})
names(median_traj) <- METRICS


# -----------------------------------------------------------------------------
# 8. SUMMARISE ACROSS THE 250 MEDIAN TRAJECTORIES (median / IQR / 95% CrI)
# -----------------------------------------------------------------------------
qfun <- function(x, p) stats::quantile(x, p, names = FALSE, na.rm = TRUE)

band <- do.call(rbind, lapply(METRICS, function(m) {
  mt <- median_traj[[m]]                   # [N_POST, n_bins]
  data.frame(
    metric = m, week = week,
    med  = apply(mt, 2, median, na.rm = TRUE),
    q25  = apply(mt, 2, qfun, 0.25),
    q75  = apply(mt, 2, qfun, 0.75),
    lo95 = apply(mt, 2, qfun, 0.025),
    hi95 = apply(mt, 2, qfun, 0.975)
  )
}))

spaghetti <- do.call(rbind, lapply(METRICS, function(m) {
  mt <- median_traj[[m]]                   # [N_POST, n_bins], column-major -> draw fastest
  data.frame(
    metric = m,
    draw   = rep(seq_len(N_POST), times = n_bins),
    week   = rep(week, each = N_POST),
    value  = as.vector(mt)
  )
}))
spaghetti <- spaghetti[is.finite(spaghetti$value), , drop = FALSE]

# Stable facet order (Deaths first).
band$metric      <- factor(band$metric,      levels = METRICS)
spaghetti$metric <- factor(spaghetti$metric, levels = METRICS)


# -----------------------------------------------------------------------------
# 9. PLOT: spaghetti + median-of-medians + IQR (dark) + 95% CrI (light)
# -----------------------------------------------------------------------------
p <- ggplot() +
  # bands first (drawn underneath): lighter 95% CrI, then darker IQR on top
  geom_ribbon(data = band, aes(week, ymin = lo95, ymax = hi95),
              fill = BAND_FILL, alpha = CI95_ALPHA) +
  geom_ribbon(data = band, aes(week, ymin = q25, ymax = q75),
              fill = BAND_FILL, alpha = IQR_ALPHA) +
  # every per-sample median trajectory
  geom_line(data = spaghetti, aes(week, value, group = draw),
            colour = SPAGHETTI_COL, alpha = SPAGHETTI_ALPHA, linewidth = SPAGHETTI_LWD) +
  # median of the 250 median trajectories
  geom_line(data = band, aes(week, med), colour = MEDIAN_COL, linewidth = MEDIAN_LWD) +
  facet_wrap(~ metric, ncol = 1, scales = "free_y") +
  labs(
    x = "Days since outbreak seeding",
    y = sprintf("Incidence per %d-day bin", BIN_WIDTH),
    title = sprintf("DRC decoupled fit: posterior trajectories (%d draws x %d replicates)",
                    N_POST, N_REPS),
    subtitle = paste0("Thin grey lines: per-sample median trajectories.  ",
                      "Bold line: median of medians.  ",
                      "Dark band: IQR (25-75%).  Light band: 95% CrI (2.5-97.5%).")
  ) +
  theme_bw(base_size = 12) +
  theme(plot.subtitle = element_text(size = 8))

if (interactive()) print(p)


# -----------------------------------------------------------------------------
# 10. SAVE FIGURE + UNDERLYING DATA
# -----------------------------------------------------------------------------
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(FIG_DIR, paste0(FIG_STEM, ".pdf")), p, width = 8, height = 8)
ggsave(file.path(FIG_DIR, paste0(FIG_STEM, ".png")), p, width = 8, height = 8, dpi = 300)

saveRDS(list(
  band        = band,
  spaghetti   = spaghetti,
  median_traj = median_traj,
  week        = week,
  n_bins      = n_bins,
  bin_width   = BIN_WIDTH,
  metrics     = METRICS,
  draw_idx    = draw_idx,
  takeoff_frac_per_draw = takeoff_frac_per_draw,
  config      = list(N_POST = N_POST, N_REPS = N_REPS, sample_seed = SAMPLE_SEED,
                     base_seed = BASE_SEED, condition_on_takeoff = CONDITION_ON_TAKEOFF,
                     rds_source = RDS_PATH, scenario_id = SCENARIO_ID)
), file.path(OUT_DIR, paste0(FIG_STEM, ".rds")))

cat("Saved:\n  ", file.path(FIG_DIR, paste0(FIG_STEM, ".pdf")),
    "\n  ", file.path(FIG_DIR, paste0(FIG_STEM, ".png")),
    "\n  ", file.path(OUT_DIR, paste0(FIG_STEM, ".rds")), "\n", sep = "")
