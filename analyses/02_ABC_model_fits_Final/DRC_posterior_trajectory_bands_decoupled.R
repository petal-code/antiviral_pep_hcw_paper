# DRC_posterior_trajectory_bands_decoupled.R
# =============================================================================
# Obeldesivir (OBV) effect on death trajectories for the DECOUPLED-efficacy DRC
# fit (Middle_DRC_ConflictSmoothed_PlusPlus): WITH vs WITHOUT OBV, as a paired
# two-run (total-effect) comparison.
#
# What it does:
#   1. Reads the ABC-SMC result.rds from the decoupled DRC run.
#   2. Draws N_POST = 250 posterior parameter samples (weighted, with
#      replacement, by the ABC importance weights).
#   3. For EACH posterior sample, runs N_REPS = 10 stochastic replicates TWICE:
#        * Without OBV: obv_pep_enabled = FALSE  (the calibrated epidemic).
#        * With OBV:    obv_pep_enabled = TRUE    (PEP 80% efficacy, coverage 1,
#                       adherence 1, delivered to HCWs exposed in hospital).
#      Realised deaths are read directly from each run's tdf (time_outcome_absolute
#      where outcome == TRUE), split into HCW vs total. This captures the TOTAL
#      OBV effect -- prevented index infections AND their averted onward chains --
#      unlike a single-run "prevented_completed" reconstruction, which only adds
#      back index deaths (~0 here, since prevented hospital-exposed HCWs mostly
#      survive) and so shows almost no difference.
#   4. Bins each replicate into a weekly death curve (and running cumulative),
#      then for EACH posterior draw takes the per-rep CENTRAL across its 10
#      replicates -> one "with" and one "without" line per draw, per arm. Across
#      the 250 draws that is 250 lines per arm.
#   5. Produces TWO figures, each with two stacked facets (weekly incidence on
#      top, cumulative on the bottom), all 250 per-draw lines per arm (faint) with
#      the cross-draw central bold on top:
#        * Figure 1: HCW deaths     (figures/DRC_decoupled_obv_HCW_deaths)
#        * Figure 2: Total deaths   (figures/DRC_decoupled_obv_total_deaths)
#
#   CENTRAL TENDENCY: defaults to the MEAN over taken-off replicates, matching how
#   the ABC fit was scored (aggregate_decoupled() logs the MEAN over reps with
#   n_deaths >= take-off threshold). So the WITHOUT-OBV arm reproduces the fitted
#   scale (~79 HCW deaths / ~2,299 total for DRC). Epidemic sizes are heavily
#   right-skewed under a near-critical branching process, so the MEDIAN sits ~4x
#   lower and HCW weekly deaths are 0 in most weeks; set CENTRAL <- "median"
#   (Section 8) only for the typical-outbreak view.
#
# Each arm is conditioned on its OWN taken-off replicates (n_deaths >= threshold),
# mirroring the fit. OBV-off and OBV-on share a seed per (draw, rep) but diverge
# once the first OBV draw fires (see branching_process_main reproducibility note),
# so the arms are compared as distributions (means), not as paired trajectories.
#
# Parallelism: replicates run with future / future.apply on a `multisession`
# (Windows-safe PSOCK) plan, one explicit seed per (sample, replicate) so results
# are reproducible regardless of worker count.
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
N_REPS    <- 10L      # stochastic replicates per posterior sample PER ARM
BIN_WIDTH <- 7L       # days per incidence bin (weekly); overridden by run metadata if present

# Per-rep -> per-draw (and cross-draw) central tendency. "mean" matches the ABC
# fit (which scored the mean over taken-off reps) and reproduces the ~79 HCW /
# ~2,299 total scale; "median" gives the typical (right-skewed) outbreak instead.
CENTRAL <- "mean"

# ---- Obeldesivir (OBV) PEP settings (the WITH-OBV arm) ----
OBV <- list(efficacy = 0.80, coverage = 1.0, adherence = 1.0, dpc = 1L,
            target_class = "HCW", target_locations = "hospital")

# ---- Seeds (reproducibility) ----
SAMPLE_SEED <- 2026L  # selects which 250 posterior particles are drawn
BASE_SEED   <- 10000L # per-(sample, replicate) model seeds start here

# ---- Model-build constants (must match the calibration run) ----
SEEDING_CASES           <- 25L
SETUP_R0_N              <- 100000L   # n for compute_R0_invariants() in the fit
SETUP_R0_SEED           <- 42L       # seed   for compute_R0_invariants() in the fit
MODEL_OVERRIDES         <- list(check_final_size = 10000)
TAKEOFF_DEATH_THRESHOLD <- 100L      # per-rep take-off (n_deaths >= this), per arm; matches the fit

# If TRUE, each draw's per-rep central uses only that arm's taken-off replicates
# (fizzles dropped), matching the fit. Set FALSE to average over all replicates.
CONDITION_ON_TAKEOFF <- TRUE

# ---- Parallel workers (Windows-safe multisession / PSOCK) ----
N_WORKERS_OVERRIDE <- NULL

# ---- Outputs ----
FIG_DIR  <- here::here("figures")
OUT_DIR  <- here::here("outputs", "02_ABC_model_fits_Final")
FIG_STEM_HCW <- "DRC_decoupled_obv_HCW_deaths"
FIG_STEM_TOT <- "DRC_decoupled_obv_total_deaths"
DATA_STEM    <- "DRC_decoupled_obv_trajectory_lines"

# ---- Plot aesthetics ----
ARM_COLS     <- c("Without OBV" = "#D55E00", "With OBV" = "#0072B2")  # Okabe-Ito (colourblind-safe)
RIBBON_ALPHA <- 0.22   # cross-draw IQR ribbon opacity, per arm
MEDIAN_LWD   <- 1.30   # bold cross-draw central line width
XLIM         <- 500    # x-axis (days) display cap; zoom only, data not dropped


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
meta   <- result$run_metadata

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

PARAM_NAMES   <- meta$fit_params    %||% DECOUPLED_PARAM_NAMES
fixed_values  <- meta$fixed_values  %||% DECOUPLED_PARAM_DEFAULTS
SCENARIO_ID   <- meta$scenario_id   %||% SCENARIO_ID_DFLT
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
cat(sprintf("OBV (with-arm)  : efficacy %.0f%%, coverage %.1f, adherence %.1f, target %s/%s\n",
            100 * OBV$efficacy, OBV$coverage, OBV$adherence, OBV$target_class, OBV$target_locations))
cat(sprintf("Sampling %d draws x %d reps x 2 arms = %d simulations; bin width = %d days.\n",
            N_POST, N_REPS, N_POST * N_REPS * 2L, BIN_WIDTH))
if (!is.null(meta$observed_summaries)) {
  os <- meta$observed_summaries
  cat(sprintf("ABC observed targets (mean over taken-off reps): total deaths ~ %.0f, HCW deaths ~ %.0f.\n",
              exp(os[["log_n_deaths"]]), exp(os[["log_n_hcw_deaths"]])))
}


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
# 5. DRAW 250 POSTERIOR SAMPLES + BUILD THE (sample, replicate, arm) JOB LIST
# -----------------------------------------------------------------------------
set.seed(SAMPLE_SEED)
draw_idx <- sample(seq_len(nrow(posterior)), size = N_POST, replace = TRUE, prob = weights)

ARM_ON     <- c(off = FALSE, on = TRUE)                       # internal key -> obv_pep_enabled
ARM_LABELS <- c(off = "Without OBV", on = "With OBV")

jobs <- vector("list", N_POST * N_REPS * length(ARM_ON))
jj <- 1L
for (k in seq_len(N_POST)) {
  for (r in seq_len(N_REPS)) {
    seed_kr <- BASE_SEED + (k - 1L) * N_REPS + r            # shared by both arms (they diverge once OBV fires)
    for (ak in names(ARM_ON)) {
      jobs[[jj]] <- list(draw = k, draw_row = draw_idx[k], rep = r,
                         arm = ak, obv_on = ARM_ON[[ak]], seed = seed_kr)
      jj <- jj + 1L
    }
  }
}


# -----------------------------------------------------------------------------
# 6. RUN ALL REPLICATES IN PARALLEL (future / multisession; Windows-safe)
# -----------------------------------------------------------------------------
# Worker builds the fiber arg list from the posterior draw, sets OBV on/off for
# this arm, runs ONE replicate, and returns the floored day of every realised
# death (total and HCW) plus the total death count (for the take-off filter).
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

  if (isTRUE(job$obv_on)) {
    a$obv_pep_enabled          <- TRUE
    a$obv_pep_coverage         <- OBV$coverage
    a$obv_pep_adherence        <- OBV$adherence
    a$obv_pep_efficacy         <- OBV$efficacy
    a$obv_pep_dpc              <- OBV$dpc
    a$obv_pep_target_class     <- OBV$target_class
    a$obv_pep_target_locations <- OBV$target_locations
  } else {
    a$obv_pep_enabled <- FALSE
  }
  a$seed <- job$seed

  out <- do.call(branching_process_main, a)

  fdays <- function(x) { x <- x[is.finite(x)]; as.integer(floor(x)) }
  tdf  <- out$tdf
  tdf  <- tdf[!is.na(tdf$time_infection_absolute), , drop = FALSE]
  died <- !is.na(tdf$outcome) & tdf$outcome
  hcw  <- !is.na(tdf$class) & tdf$class == "HCW"

  total_days <- fdays(tdf$time_outcome_absolute[died])
  hcw_days   <- fdays(tdf$time_outcome_absolute[died & hcw])
  list(draw = job$draw, rep = job$rep, arm = job$arm,
       total_days = total_days, hcw_days = hcw_days, n_total = length(total_days))
}

N_CLUSTER <- if (parallel::detectCores() > 120) {
  min(120L, parallel::detectCores() - 10L)
} else {
  min(12L, parallel::detectCores() - 4L)
}
N_WORKERS <- N_WORKERS_OVERRIDE %||% max(1L, min(N_CLUSTER, future::availableCores()))

cat(sprintf("Running %d simulations on %d workers (multisession)...\n", length(jobs), N_WORKERS))
start_time <- Sys.time()

plan(multisession, workers = N_WORKERS)
results <- future_lapply(jobs, run_traj_job,
                         future.packages = "fiber", future.seed = TRUE)
plan(sequential)

cat("Simulation time: "); print(Sys.time() - start_time)


# -----------------------------------------------------------------------------
# 7. BIN TO WEEKLY DEATHS (+ CUMULATIVE) PER ARM
# -----------------------------------------------------------------------------
max_day <- max(c(unlist(lapply(results, function(r) r$total_days)), 0L))
n_bins  <- max(1L, ceiling((max_day + 1L) / BIN_WIDTH))
week    <- ((seq_len(n_bins) - 1L) + 0.5) * BIN_WIDTH

METRIC_KEYS <- c(hcw = "HCW deaths", tot = "Total deaths")
MEASURES    <- c(inc = "Weekly incidence", cum = "Cumulative")

# A[[arm]][[metric]]$inc / $cum : [N_POST x N_REPS x n_bins]; A[[arm]]$took : [N_POST x N_REPS]
mk_arr <- function() array(0L, dim = c(N_POST, N_REPS, n_bins))
A <- setNames(lapply(names(ARM_LABELS), function(ak)
  list(hcw  = list(inc = mk_arr(), cum = mk_arr()),
       tot  = list(inc = mk_arr(), cum = mk_arr()),
       took = matrix(FALSE, N_POST, N_REPS))), names(ARM_LABELS))

for (res in results) {
  ak <- res$arm; d <- res$draw; rp <- res$rep
  bh <- bin_counts(res$hcw_days,   BIN_WIDTH, n_bins)
  bt <- bin_counts(res$total_days, BIN_WIDTH, n_bins)
  A[[ak]]$hcw$inc[d, rp, ] <- bh; A[[ak]]$hcw$cum[d, rp, ] <- cumsum(bh)
  A[[ak]]$tot$inc[d, rp, ] <- bt; A[[ak]]$tot$cum[d, rp, ] <- cumsum(bt)
  A[[ak]]$took[d, rp] <- res$n_total >= TAKEOFF_DEATH_THRESHOLD
}

# Magnitude check vs the ABC targets: mean total/HCW deaths over taken-off reps,
# per arm. WITHOUT OBV should land near the fitted ~2,299 total / ~79 HCW.
for (ak in names(ARM_LABELS)) {
  tk <- A[[ak]]$took
  ft <- A[[ak]]$tot$cum[, , n_bins]; fh <- A[[ak]]$hcw$cum[, , n_bins]
  cat(sprintf("%-12s: take-off %.1f%% of reps; mean over taken-off reps -> total = %.0f, HCW = %.1f deaths.\n",
              ARM_LABELS[[ak]], 100 * mean(tk),
              if (any(tk)) mean(ft[tk]) else NA_real_,
              if (any(tk)) mean(fh[tk]) else NA_real_))
}


# -----------------------------------------------------------------------------
# 8. PER-DRAW LINES (per-rep CENTRAL within each draw, conditioned on take-off)
# -----------------------------------------------------------------------------
agg_fun <- if (identical(CENTRAL, "mean")) {
  function(x) mean(x, na.rm = CONDITION_ON_TAKEOFF)
} else {
  function(x) median(x, na.rm = CONDITION_ON_TAKEOFF)
}

# arr3d [N_POST x N_REPS x n_bins] + that arm's took [N_POST x N_REPS]
#   -> mt [N_POST x n_bins] : per-draw central across (taken-off) reps.
per_draw_central <- function(arr3d, took_mat) {
  if (CONDITION_ON_TAKEOFF) arr3d[!took_mat] <- NA_integer_   # drop fizzled reps first (matches the fit)
  mt <- apply(arr3d, c(1, 3), agg_fun)
  mt[!is.finite(mt)] <- NA_real_                              # draws with no take-off at all -> NA line
  mt
}

qfun <- function(x, p) stats::quantile(x, p, names = FALSE, na.rm = TRUE)

draws_list <- list(); band_list <- list(); li <- 1L
for (ak in names(ARM_LABELS)) for (mk in names(METRIC_KEYS)) for (meas in names(MEASURES)) {
  mt <- per_draw_central(A[[ak]][[mk]][[meas]], A[[ak]]$took)   # [N_POST x n_bins]
  draws_list[[li]] <- data.frame(
    draw    = rep(seq_len(N_POST), times = n_bins),
    week    = rep(week, each = N_POST),
    value   = as.vector(mt),
    arm     = ARM_LABELS[[ak]], metric = METRIC_KEYS[[mk]], measure = MEASURES[[meas]],
    stringsAsFactors = FALSE
  )
  band_list[[li]] <- data.frame(
    week    = week,
    med     = apply(mt, 2, agg_fun),       # cross-draw central of the per-draw centrals
    lo      = apply(mt, 2, qfun, 0.25),    # cross-draw IQR of the per-draw centrals
    hi      = apply(mt, 2, qfun, 0.75),
    arm     = ARM_LABELS[[ak]], metric = METRIC_KEYS[[mk]], measure = MEASURES[[meas]],
    stringsAsFactors = FALSE
  )
  li <- li + 1L
}
draws_long <- do.call(rbind, draws_list)
band       <- do.call(rbind, band_list)
for (df_nm in c("draws_long", "band")) {
  assign(df_nm, within(get(df_nm), {
    arm     <- factor(arm,     levels = ARM_LABELS)
    measure <- factor(measure, levels = MEASURES)              # Weekly incidence facet on top
  }))
}

# Headline: cross-draw central cumulative deaths at the final week, with vs
# without OBV, and the implied total % averted.
for (mlab in METRIC_KEYS) {
  d  <- band[band$metric == mlab & band$measure == "Cumulative", ]
  fw <- max(d$week)
  no <- d$med[d$arm == "Without OBV" & d$week == fw]
  ob <- d$med[d$arm == "With OBV"    & d$week == fw]
  cat(sprintf("%-12s cumulative deaths (cross-draw %s): without OBV = %.1f, with OBV = %.1f (averted %.1f%%)\n",
              mlab, CENTRAL, no, ob, if (length(no) && no > 0) 100 * (no - ob) / no else NA_real_))
}


# -----------------------------------------------------------------------------
# 9. PLOT: per-arm cross-draw central line + IQR ribbon, distinct colours
#          two facets: weekly incidence (top), cumulative (bottom)
# -----------------------------------------------------------------------------
make_fig <- function(metric_label) {
  d <- band[band$metric == metric_label, , drop = FALSE]
  ggplot(d, aes(week, med, colour = arm, fill = arm)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = RIBBON_ALPHA, colour = NA) +
    geom_line(linewidth = MEDIAN_LWD, na.rm = TRUE) +
    facet_wrap(~ measure, ncol = 1, scales = "free_y") +
    scale_colour_manual(values = ARM_COLS) +
    scale_fill_manual(values = ARM_COLS) +
    coord_cartesian(xlim = c(0, XLIM)) +
    labs(
      x = "Days since outbreak seeding", y = metric_label,
      colour = NULL, fill = NULL,
      title = sprintf("DRC decoupled fit: %s, obeldesivir effect", tolower(metric_label)),
      subtitle = sprintf(paste0("Two runs per draw x rep: OBV off vs on (PEP %.0f%% efficacy, coverage %.0f, adherence %.0f, HCW in hospital).  ",
                                "Line: cross-draw %s of the per-draw %ss across %d reps (%d draws); band: cross-draw IQR.  ",
                                "'Without OBV' reproduces the ABC fit scale."),
                         100 * OBV$efficacy, OBV$coverage, OBV$adherence, CENTRAL, CENTRAL, N_REPS, N_POST)
    ) +
    theme_bw(base_size = 12) +
    theme(plot.subtitle = element_text(size = 7), legend.position = "top")
}

fig_hcw <- make_fig("HCW deaths")
fig_tot <- make_fig("Total deaths")
if (interactive()) { print(fig_hcw); print(fig_tot) }


# -----------------------------------------------------------------------------
# 10. SAVE FIGURES + UNDERLYING DATA
# -----------------------------------------------------------------------------
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(FIG_DIR, paste0(FIG_STEM_HCW, ".pdf")), fig_hcw, width = 8, height = 8)
ggsave(file.path(FIG_DIR, paste0(FIG_STEM_HCW, ".png")), fig_hcw, width = 8, height = 8, dpi = 300)
ggsave(file.path(FIG_DIR, paste0(FIG_STEM_TOT, ".pdf")), fig_tot, width = 8, height = 8)
ggsave(file.path(FIG_DIR, paste0(FIG_STEM_TOT, ".png")), fig_tot, width = 8, height = 8, dpi = 300)

saveRDS(list(
  draws_long = draws_long,   # N_POST per-draw lines per arm/metric/measure (for spaghetti if wanted)
  band       = band,         # cross-draw central + IQR of those per-draw lines (plotted)
  week       = week,
  n_bins     = n_bins,
  bin_width  = BIN_WIDTH,
  central_fun = CENTRAL,
  obv        = OBV,
  effect     = "total (two runs: OBV off vs on; includes averted onward chains)",
  draw_idx   = draw_idx,
  config     = list(N_POST = N_POST, N_REPS = N_REPS, sample_seed = SAMPLE_SEED,
                    base_seed = BASE_SEED, condition_on_takeoff = CONDITION_ON_TAKEOFF,
                    takeoff_death_threshold = TAKEOFF_DEATH_THRESHOLD,
                    rds_source = RDS_PATH, scenario_id = SCENARIO_ID)
), file.path(OUT_DIR, paste0(DATA_STEM, ".rds")))

cat("Saved:\n  ",
    file.path(FIG_DIR, paste0(FIG_STEM_HCW, ".pdf")), " (+ .png)\n  ",
    file.path(FIG_DIR, paste0(FIG_STEM_TOT, ".pdf")), " (+ .png)\n  ",
    file.path(OUT_DIR, paste0(DATA_STEM, ".rds")), "\n", sep = "")
