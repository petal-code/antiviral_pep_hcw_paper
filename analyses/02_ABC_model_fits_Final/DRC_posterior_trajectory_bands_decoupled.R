# DRC_posterior_trajectory_bands_decoupled.R
# =============================================================================
# Direct obeldesivir (OBV) effect on death trajectories for the DECOUPLED-efficacy
# DRC fit (Middle_DRC_ConflictSmoothed_PlusPlus): WITH vs WITHOUT OBV, from a
# SINGLE OBV run per replicate (no separate baseline arm).
#
# What it does:
#   1. Reads the ABC-SMC result.rds from the decoupled DRC run.
#   2. Draws N_POST = 250 posterior parameter samples (weighted, with
#      replacement, by the ABC importance weights).
#   3. For EACH posterior sample, runs N_REPS = 10 stochastic replicates with
#      obeldesivir PEP ON (80% efficacy, coverage 1, adherence 1, delivered to
#      HCWs exposed in hospital). Each run yields BOTH curves:
#        * WITH OBV     = the realised deaths in the simulated tree (out$tdf).
#        * WITHOUT OBV  = realised deaths PLUS the would-be deaths OBV directly
#          (direct)       prevented -- fiber's out$prevented_completed (the
#                         averted index infections replayed through the same
#                         outcome model), counted at their would-be death time
#                         (time_outcome_absolute where outcome == TRUE).
#      This is the DIRECT effect only: it adds back the prevented index deaths,
#      NOT their averted onward transmission chains (that is the total effect).
#   4. Bins each replicate into a weekly death-incidence curve (and its running
#      cumulative), then for EACH posterior draw takes the median across its 10
#      replicates -> ONE "with OBV" and ONE "without OBV" median line per draw.
#      Across the 250 draws that is 250 lines per arm.
#   5. Produces TWO figures, each with two stacked facets (weekly incidence on
#      top, cumulative on the bottom). Each plots all 250 per-draw median lines
#      per arm (faint), with the cross-draw median drawn bold on top:
#        * Figure 1: HCW deaths     (figures/DRC_decoupled_obv_direct_HCW_deaths)
#        * Figure 2: Total deaths   (figures/DRC_decoupled_obv_direct_total_deaths)
#
#      NOTE on HCW weekly incidence: HCW deaths are rare (a handful per epidemic),
#      so a given week usually has 0 in most replicates and the per-draw MEDIAN
#      across replicates is 0 for most weeks -- the incidence facet is therefore
#      sparse/spiky for HCWs by construction. The cumulative facet is unaffected.
#      Flip the per-rep aggregation to mean (per_draw_medians) if a smooth HCW
#      incidence curve is wanted instead.
#
# Because both curves come from the SAME run, the comparison is exactly paired:
# the realised epidemic is identical and "without OBV" simply adds the prevented
# deaths on top.
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
N_REPS    <- 10L      # stochastic replicates per posterior sample
BIN_WIDTH <- 7L       # days per incidence bin (weekly); overridden by run metadata if present

# ---- Obeldesivir (OBV) PEP settings ----
# Standard delivery from the codebase (OBV_BASE): full coverage/adherence,
# modelled as post-exposure prophylaxis for HCWs exposed in hospital.
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
TAKEOFF_DEATH_THRESHOLD <- 100L      # take-off defined on the realised (with-OBV) total deaths

# If TRUE, per-sample medians use only replicates whose epidemic took off
# (>= TAKEOFF_DEATH_THRESHOLD realised deaths); fizzles are dropped. Same mask
# for both curves (they share the run).
CONDITION_ON_TAKEOFF <- TRUE

# ---- Parallel workers (Windows-safe multisession / PSOCK) ----
# NULL = auto. Mirrors the calibration script's heuristic; cap as you like.
N_WORKERS_OVERRIDE <- NULL

# ---- Outputs ----
FIG_DIR  <- here::here("figures")
OUT_DIR  <- here::here("outputs", "02_ABC_model_fits_Final")
FIG_STEM_HCW <- "DRC_decoupled_obv_direct_HCW_deaths"
FIG_STEM_TOT <- "DRC_decoupled_obv_direct_total_deaths"
DATA_STEM    <- "DRC_decoupled_obv_direct_trajectory_bands"

# ---- Plot aesthetics ----
ARM_COLS        <- c("Without OBV" = "#D55E00", "With OBV" = "#0072B2")  # Okabe-Ito (colourblind-safe)
SPAGHETTI_ALPHA <- 0.06   # opacity of the N_POST faint per-draw median lines
MEDIAN_LWD      <- 1.20   # bold cross-draw median line width
XLIM            <- 500    # x-axis (days) display cap; zoom only, data not dropped


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
cat(sprintf("OBV (direct)    : efficacy %.0f%%, coverage %.1f, adherence %.1f, target %s/%s\n",
            100 * OBV$efficacy, OBV$coverage, OBV$adherence, OBV$target_class, OBV$target_locations))
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

jobs <- vector("list", N_POST * N_REPS)
jj <- 1L
for (k in seq_len(N_POST)) {
  for (r in seq_len(N_REPS)) {
    jobs[[jj]] <- list(draw = k, draw_row = draw_idx[k], rep = r,
                       seed = BASE_SEED + (k - 1L) * N_REPS + r)
    jj <- jj + 1L
  }
}


# -----------------------------------------------------------------------------
# 6. RUN ALL REPLICATES IN PARALLEL (future / multisession; Windows-safe)
# -----------------------------------------------------------------------------
# Worker builds the fiber arg list from the posterior draw, switches OBV PEP on,
# runs ONE replicate, and returns: the floored day of every realised death (total
# and HCW), and the floored would-be death day of every OBV-prevented index
# infection (total and HCW), read from out$prevented_completed.
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

  a$obv_pep_enabled          <- TRUE
  a$obv_pep_coverage         <- OBV$coverage
  a$obv_pep_adherence        <- OBV$adherence
  a$obv_pep_efficacy         <- OBV$efficacy
  a$obv_pep_dpc              <- OBV$dpc
  a$obv_pep_target_class     <- OBV$target_class
  a$obv_pep_target_locations <- OBV$target_locations
  a$seed <- job$seed

  out <- do.call(branching_process_main, a)

  fdays <- function(x) { x <- x[is.finite(x)]; as.integer(floor(x)) }

  # Realised deaths WITH OBV (from the simulated tree).
  tdf  <- out$tdf
  tdf  <- tdf[!is.na(tdf$time_infection_absolute), , drop = FALSE]
  died <- !is.na(tdf$outcome) & tdf$outcome
  hcw  <- !is.na(tdf$class) & tdf$class == "HCW"

  # Would-be deaths OBV directly prevented (averted index infections only), at
  # their counterfactual death time.
  pc <- out$prevented_completed
  if (!is.null(pc) && nrow(pc) > 0) {
    pdied <- !is.na(pc$outcome) & pc$outcome
    phcw  <- !is.na(pc$class) & pc$class == "HCW"
    prevented_total_days <- fdays(pc$time_outcome_absolute[pdied])
    prevented_hcw_days   <- fdays(pc$time_outcome_absolute[pdied & phcw])
  } else {
    prevented_total_days <- integer(0)
    prevented_hcw_days   <- integer(0)
  }

  list(draw = job$draw, rep = job$rep,
       realised_total_days  = fdays(tdf$time_outcome_absolute[died]),
       realised_hcw_days    = fdays(tdf$time_outcome_absolute[died & hcw]),
       prevented_total_days = prevented_total_days,
       prevented_hcw_days   = prevented_hcw_days)
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
# 7. BIN TO WEEKLY INCIDENCE (+ CUMULATIVE) -> PER-SAMPLE MEDIAN TRAJECTORIES
# -----------------------------------------------------------------------------
# Common weekly grid across ALL replicates (absolute model time, from day 0).
max_day <- max(c(unlist(lapply(results, function(r) c(r$realised_total_days, r$prevented_total_days))), 0L))
n_bins  <- max(1L, ceiling((max_day + 1L) / BIN_WIDTH))
week    <- ((seq_len(n_bins) - 1L) + 0.5) * BIN_WIDTH

ARM_LABELS  <- c(noobv = "Without OBV", obv = "With OBV")
METRIC_KEYS <- c(hcw = "HCW deaths", tot = "Total deaths")
MEASURES    <- c(inc = "Weekly incidence", cum = "Cumulative")

# A[[arm]][[metric]]$inc / $cum : [N_POST x N_REPS x n_bins] weekly counts.
mk_arr <- function() array(0L, dim = c(N_POST, N_REPS, n_bins))
A <- setNames(lapply(names(ARM_LABELS), function(ak)
  setNames(lapply(names(METRIC_KEYS), function(mk) list(inc = mk_arr(), cum = mk_arr())),
           names(METRIC_KEYS))), names(ARM_LABELS))

# WITH OBV   = realised deaths; WITHOUT OBV = realised + directly-prevented.
took <- matrix(FALSE, N_POST, N_REPS)
for (res in results) {
  d <- res$draw; rp <- res$rep
  days <- list(
    obv   = list(hcw = res$realised_hcw_days,
                 tot = res$realised_total_days),
    noobv = list(hcw = c(res$realised_hcw_days,   res$prevented_hcw_days),
                 tot = c(res$realised_total_days, res$prevented_total_days))
  )
  for (ak in names(days)) for (mk in names(days[[ak]])) {
    b <- bin_counts(days[[ak]][[mk]], BIN_WIDTH, n_bins)
    A[[ak]][[mk]]$inc[d, rp, ] <- b
    A[[ak]][[mk]]$cum[d, rp, ] <- cumsum(b)
  }
  took[d, rp] <- length(res$realised_total_days) >= TAKEOFF_DEATH_THRESHOLD
}

takeoff_frac_per_draw <- rowMeans(took)
cat(sprintf("Take-off (>= %d realised deaths): %.1f%% of replicates; per-sample take-off fraction median = %.2f.\n",
            TAKEOFF_DEATH_THRESHOLD, 100 * mean(took), median(takeoff_frac_per_draw)))

# Magnitude check (sums over ALL sims): how big is the DIRECT OBV effect on
# deaths? If "OBV-prevented HCW deaths" is ~0 the with/without lines will sit on
# top of each other -- that is the direct effect being genuinely small (most of
# OBV's benefit is indirect, via averted onward chains, which this design does
# NOT count), not a plotting bug. Compare against realised to gauge the ratio.
tot_real_hcw <- sum(vapply(results, function(r) length(r$realised_hcw_days),  0L))
tot_real_all <- sum(vapply(results, function(r) length(r$realised_total_days), 0L))
tot_prev_hcw <- sum(vapply(results, function(r) length(r$prevented_hcw_days),  0L))
tot_prev_all <- sum(vapply(results, function(r) length(r$prevented_total_days), 0L))
cat(sprintf(paste0("Death counts summed over all %d sims:\n",
                   "  realised:  HCW = %d, total = %d\n",
                   "  prevented: HCW = %d, total = %d  (direct effect; %.1f%% of realised HCW deaths)\n"),
            length(results), tot_real_hcw, tot_real_all, tot_prev_hcw, tot_prev_all,
            if (tot_real_hcw > 0) 100 * tot_prev_hcw / tot_real_hcw else NA_real_))


# -----------------------------------------------------------------------------
# 8. PER-DRAW MEDIAN TRAJECTORIES (median across the N_REPS replicates per draw)
# -----------------------------------------------------------------------------
# For each posterior draw: take the median across its N_REPS replicate weekly
# trajectories -> ONE "With OBV" median line and ONE "Without OBV" median line
# per draw. Repeated over all N_POST draws this gives N_POST lines per arm,
# which we plot directly (Section 9) rather than collapsing to a single band.
qfun <- function(x, p) stats::quantile(x, p, names = FALSE, na.rm = TRUE)

# arr3d [N_POST x N_REPS x n_bins] -> mt [N_POST x n_bins]: per-draw median over reps.
per_draw_medians <- function(arr3d) {
  if (CONDITION_ON_TAKEOFF) arr3d[!took] <- NA_integer_  # drop fizzled reps before the per-draw median
  mt <- apply(arr3d, c(1, 3), median, na.rm = CONDITION_ON_TAKEOFF)
  mt[!is.finite(mt)] <- NA_real_                         # draws with no take-off at all -> NA line
  mt
}

# Long frame of all N_POST per-draw median lines (one row per draw x week), for
# every arm x metric x measure; plus the bold cross-draw median for readability.
draws_list <- list(); central_list <- list(); li <- 1L
for (ak in names(ARM_LABELS)) for (mk in names(METRIC_KEYS)) for (meas in names(MEASURES)) {
  mt <- per_draw_medians(A[[ak]][[mk]][[meas]])               # [N_POST x n_bins]
  draws_list[[li]] <- data.frame(
    draw    = rep(seq_len(N_POST), times = n_bins),
    week    = rep(week, each = N_POST),
    value   = as.vector(mt),
    arm     = ARM_LABELS[[ak]], metric = METRIC_KEYS[[mk]], measure = MEASURES[[meas]],
    stringsAsFactors = FALSE
  )
  central_list[[li]] <- data.frame(
    week    = week,
    med     = apply(mt, 2, median, na.rm = TRUE),
    arm     = ARM_LABELS[[ak]], metric = METRIC_KEYS[[mk]], measure = MEASURES[[meas]],
    stringsAsFactors = FALSE
  )
  li <- li + 1L
}
draws_long <- do.call(rbind, draws_list)
central    <- do.call(rbind, central_list)
for (df_nm in c("draws_long", "central")) {
  assign(df_nm, within(get(df_nm), {
    arm     <- factor(arm,     levels = ARM_LABELS)
    measure <- factor(measure, levels = MEASURES)            # Weekly incidence facet on top
  }))
}

# Headline diagnostic: cross-draw median of the per-draw cumulative-death lines
# at the final week, with vs without OBV, and the implied direct % averted.
for (mlab in METRIC_KEYS) {
  d  <- central[central$metric == mlab & central$measure == "Cumulative", ]
  fw <- max(d$week)
  no <- d$med[d$arm == "Without OBV" & d$week == fw]
  ob <- d$med[d$arm == "With OBV"    & d$week == fw]
  cat(sprintf("%-12s cumulative deaths (cross-draw median): without OBV = %.1f, with OBV = %.1f (direct averted %.1f%%)\n",
              mlab, no, ob, if (length(no) && no > 0) 100 * (no - ob) / no else NA_real_))
}


# -----------------------------------------------------------------------------
# 9. PLOT: N_POST per-draw median lines per arm (faint), bold cross-draw median
#          two facets: weekly incidence (top), cumulative (bottom)
# -----------------------------------------------------------------------------
make_fig <- function(metric_label) {
  dl <- draws_long[draws_long$metric == metric_label, , drop = FALSE]
  ce <- central[central$metric == metric_label, , drop = FALSE]
  ggplot() +
    # one faint line per posterior draw (its median across the N_REPS replicates)
    geom_line(data = dl,
              aes(week, value, group = interaction(arm, draw), colour = arm),
              alpha = SPAGHETTI_ALPHA, linewidth = 0.3, na.rm = TRUE) +
    # bold line: median across the N_POST per-draw medians
    geom_line(data = ce, aes(week, med, colour = arm),
              linewidth = MEDIAN_LWD, na.rm = TRUE) +
    facet_wrap(~ measure, ncol = 1, scales = "free_y") +
    scale_colour_manual(values = ARM_COLS) +
    coord_cartesian(xlim = c(0, XLIM)) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, linewidth = 1.2))) +
    labs(
      x = "Days since outbreak seeding", y = metric_label,
      colour = NULL,
      title = sprintf("DRC decoupled fit: %s, direct obeldesivir effect", tolower(metric_label)),
      subtitle = sprintf(paste0("Single OBV run (PEP %.0f%% efficacy, coverage %.0f, adherence %.0f, HCW in hospital).  ",
                                "'Without OBV' = realised + directly-prevented index deaths (excludes averted chains).  ",
                                "Faint lines: each draw's median across %d reps (%d draws); bold: cross-draw median."),
                         100 * OBV$efficacy, OBV$coverage, OBV$adherence, N_REPS, N_POST)
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
  draws_long = draws_long,   # N_POST per-draw median lines per arm/metric/measure
  central    = central,      # cross-draw median of those per-draw medians
  week       = week,
  n_bins     = n_bins,
  bin_width  = BIN_WIDTH,
  obv        = OBV,
  effect     = "direct (realised + prevented index deaths; excludes averted onward chains)",
  draw_idx   = draw_idx,
  takeoff_frac_per_draw = takeoff_frac_per_draw,
  config     = list(N_POST = N_POST, N_REPS = N_REPS, sample_seed = SAMPLE_SEED,
                    base_seed = BASE_SEED, condition_on_takeoff = CONDITION_ON_TAKEOFF,
                    rds_source = RDS_PATH, scenario_id = SCENARIO_ID)
), file.path(OUT_DIR, paste0(DATA_STEM, ".rds")))

cat("Saved:\n  ",
    file.path(FIG_DIR, paste0(FIG_STEM_HCW, ".pdf")), " (+ .png)\n  ",
    file.path(FIG_DIR, paste0(FIG_STEM_TOT, ".pdf")), " (+ .png)\n  ",
    file.path(OUT_DIR, paste0(DATA_STEM, ".rds")), "\n", sep = "")
