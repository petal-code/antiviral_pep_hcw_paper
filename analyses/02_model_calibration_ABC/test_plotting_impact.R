# DRC_posterior_trajectory_bands_decoupled.R
# =============================================================================
# Obeldesivir (OBV) effect on death trajectories for the DECOUPLED-efficacy fits
# (DRC + West Africa), computed from the WITHIN-RUN counterfactual of a SINGLE
# OBV-enabled run per (draw, rep) -- no separate no-OBV arm.
#
# Each branching_process_main() call with obv_pep_enabled = TRUE returns, after
# the tree is finalised (the counterfactual draws happen AFTER the trajectory and
# consume no RNG inside it, so the simulated epidemic is byte-for-byte unchanged):
#   * tdf                  -- the realised (WITH OBV) epidemic; deaths are
#                             time_outcome_absolute where outcome == TRUE.
#   * prevented_completed  -- the infections the OBV gate prevented (the averted
#                             INDEX infections, replayed through the same outcome
#                             model), each carrying its would-be `class`,
#                             `outcome`, and `time_outcome_absolute`.
# So, per run:
#   WITH OBV    deaths = realised (tdf).
#   WITHOUT OBV deaths = realised + prevented index   (direct counterfactual).
#   OBV averted        = prevented index deaths.
# Both curves come from ONE run, perfectly paired -- no OBV RNG desync, no
# cross-arm take-off mismatch, and the comparison can never go negative.
#
# ESTIMAND: prevented_completed is the DIRECT effect -- the would-be deaths of
# the averted INDEX infections only; it excludes the onward chains those
# infections would have seeded. (The gap to the TOTAL effect is exactly those
# averted chains, which a single run cannot contain.) The "Without OBV" band is
# therefore the conservative, direct within-run counterfactual.
#
# It bins each run's two curves into weekly incidence + running cumulative, takes
# the per-rep CENTRAL within each draw, summarises ACROSS draws (median line +
# IQR + 95% CrI), and PRINTS, in the raw-ggplot style:
#   * per scenario: HCW-deaths and total-deaths band figures (incidence facet on
#     top, cumulative below);                                    -> 4 figures
#   * one final 2x2 cowplot panel of HCW deaths -- weekly incidence (top row),
#     cumulative (bottom row), one column per scenario, rel_heights = c(1, 2).
#                                                                 -> 1 figure
# NOTHING is written to disk: every figure is print()ed to the active device.
#
# Run from anywhere inside the project (paths resolved with here::here()):
#   source("analyses/02_model_calibration_ABC/DRC_posterior_trajectory_bands_decoupled.R")
# =============================================================================


# -----------------------------------------------------------------------------
# 1. CONFIGURATION
# -----------------------------------------------------------------------------
FUNCTIONS_DIR <- here::here("functions")
ANALYSIS_DIR  <- here::here("analyses", "02_model_calibration_ABC")
SCENARIO_CSV  <- here::here("data-processed", "final_six_scenario_values_original_approach.csv")

# Both decoupled ABC-SMC results. `rds_subdir` is the folder under
# <ANALYSIS_DIR>/abc_outputs/ that holds result.rds; `scenario_id_default` is
# only a fallback used if the run metadata does not carry a scenario id.
SCENARIOS <- list(
  list(label = "DRC",
       rds_subdir = "Middle_DRC_ConflictSmoothed_PlusPlus_20260607_215621_Decoupled_check_NP5_NS4_NBREPS_30_NBSIMUL_590",
       scenario_id_default = "Middle_DRC_ConflictSmoothed_PlusPlus"),
  list(label = "West Africa",
       rds_subdir = "Worst_WestAfrica_20260608_162044_Decoupled_check_NP5_NS4_NBREPS_30_NBSIMUL_472",
       scenario_id_default = "Worst_WestAfrica")
)
SCEN_LEVELS <- vapply(SCENARIOS, function(s) s$label, character(1))  # column order: DRC, West Africa

# ---- Sampling / replication ----
N_POST    <- 250L     # posterior parameter samples per scenario
N_REPS    <- 10L      # stochastic replicates per posterior sample (ONE OBV run each)
BIN_WIDTH <- 7L       # days per incidence bin (weekly); overridden by run metadata if present

# Per-rep -> per-draw central tendency (WITHIN a draw, across reps). "mean"
# matches the ABC fit (scored on the mean over taken-off reps); "median" gives
# the typical (right-skewed) outbreak.
CENTRAL <- "mean"

# Cross-draw central tendency for the plotted LINE (ACROSS the N_POST draws).
# This is what the figures you liked used (median); the bands are always the
# cross-draw IQR + 95% CrI of the per-draw lines.
BAND_CENTRAL <- "median"   # "median" or "mean"

# ---- Obeldesivir (OBV) PEP settings (always enabled; counterfactual is internal) ----
OBV <- list(efficacy = 0.80, coverage = 1.0, adherence = 1.0, dpc = 1L,
            target_class = "HCW", target_locations = "hospital")

# ---- Seeds (reproducibility) ----
SAMPLE_SEED <- 2026L  # selects which N_POST posterior particles are drawn (per scenario)
BASE_SEED   <- 10000L # per-(sample, replicate) model seeds start here

# ---- Model-build constants (must match the calibration runs) ----
SEEDING_CASES           <- 25L
SETUP_R0_N              <- 100000L
SETUP_R0_SEED           <- 42L
MODEL_OVERRIDES         <- list(check_final_size = 10000)
TAKEOFF_DEATH_THRESHOLD <- 100L
CONDITION_ON_TAKEOFF    <- TRUE

# ---- Parallel workers (Windows-safe multisession / PSOCK) ----
N_WORKERS_OVERRIDE <- NULL

# ---- Plot aesthetics ----
# Internal arm keys -> labels: "off" is the (counterfactual) no-OBV curve,
# "on" is the realised OBV curve. Both are reconstructed from the SAME run.
ARM_LABELS     <- c(off = "Without OBV", on = "With OBV")
ARM_LEVELS     <- unname(ARM_LABELS)                                # factor / colour order
ARM_COLS       <- c("Without OBV" = "#D55E00", "With OBV" = "#0072B2")  # Okabe-Ito (colourblind-safe)
MEASURE_LEVELS <- c("Weekly incidence", "Cumulative")               # facet / row order (incidence first)
XLIM           <- 500       # x-axis (days) display cap; zoom only, data not dropped
CRI_ALPHA      <- 0.15      # 95% CrI ribbon opacity (lighter)
IQR_ALPHA      <- 0.30      # IQR ribbon opacity (darker)
BAND_LWD       <- 1.05      # cross-draw central line width


# -----------------------------------------------------------------------------
# 2. LIBRARIES + SOURCES  (order required by the decoupled scheme)
# -----------------------------------------------------------------------------
library(fiber)
library(future)
library(future.apply)
library(ggplot2)
library(cowplot)   # plot_grid(rel_heights = ...) for the final 2x2 panel

source(file.path(FUNCTIONS_DIR, "setup_model_parameters.R"))              # make_model_parameters(), read_scenario_matrix(), DEFAULT_SCALAR_INPUTS
source(file.path(FUNCTIONS_DIR, "abc_calibration_functions_common.R"))    # abc_summarise(), check_model_function_version()
source(file.path(FUNCTIONS_DIR, "abc_calibration_functions_decoupled.R")) # assemble/build args, DECOUPLED_PARAM_NAMES/DEFAULTS
source(file.path(FUNCTIONS_DIR, "calculate_model_approx_r0.R"))           # compute_R0_invariants(), D/F_from_invariants(), solve_offspring_means()
source(file.path(FUNCTIONS_DIR, "simulation_helpers.R"))                  # bin_counts()

tryCatch(check_model_function_version(),
         error = function(e) warning("check_model_function_version() failed: ", conditionMessage(e)))


# -----------------------------------------------------------------------------
# 3. SHARED HELPERS + SCENARIO INPUTS
# -----------------------------------------------------------------------------
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

# The scenario CSV holds every scenario; the per-scenario id selects rows inside
# make_model_parameters(), so we only need to read it once.
scenario_matrix <- read_scenario_matrix(SCENARIO_CSV)

METRIC_KEYS <- c(hcw = "HCW deaths", tot = "Total deaths")
MEASURES_K  <- c(inc = "Weekly incidence", cum = "Cumulative")

# Per-rep -> per-draw aggregator (WITHIN a draw, across reps).
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

# Cross-draw central for the plotted line, and a quantile helper for the bands.
central_fun <- if (identical(BAND_CENTRAL, "mean")) function(x) mean(x, na.rm = TRUE) else
                                                    function(x) median(x, na.rm = TRUE)
qf <- function(x, p) stats::quantile(x, p, names = FALSE, na.rm = TRUE)


# -----------------------------------------------------------------------------
# 4. PER-SCENARIO ENGINE: posterior draws -> single OBV runs -> per-draw lines
# -----------------------------------------------------------------------------
# Runs ONE OBV-enabled replicate per (draw, rep) for a scenario, reconstructs the
# WITH-OBV (realised) and WITHOUT-OBV (realised + prevented_completed) death
# curves from that single run, and returns the per-draw central lines in long
# form (one row per draw x week x arm x metric x measure), tagged with `scenario`.
compute_scenario_draws <- function(scn) {

  rds_path <- file.path(ANALYSIS_DIR, "abc_outputs", scn$rds_subdir, "result.rds")
  if (!file.exists(rds_path)) {
    stop("Could not find the decoupled result for '", scn$label, "' at:\n  ", rds_path,
         "\nFix `rds_subdir` for this scenario in the SCENARIOS list.", call. = FALSE)
  }

  ## ---- load posterior + run metadata ----
  result <- readRDS(rds_path)
  meta   <- result$run_metadata

  PARAM_NAMES   <- meta$fit_params    %||% DECOUPLED_PARAM_NAMES
  fixed_values  <- meta$fixed_values  %||% DECOUPLED_PARAM_DEFAULTS
  SCENARIO_ID   <- meta$scenario_id   %||% scn$scenario_id_default
  HCW_BASE_PROB <- meta$hcw_base_prob %||% 0.25
  GEN_HOSP_EFF  <- meta$fixed_efficacies$general_hospital_quarantine_efficacy %||%
                     DEFAULT_SCALAR_INPUTS$general_hospital_quarantine_efficacy
  SAFE_FUN_EFF  <- meta$fixed_efficacies$safe_funeral_efficacy %||%
                     DEFAULT_SCALAR_INPUTS$safe_funeral_efficacy
  bin_width     <- as.integer(meta$peak_settings$bin_width %||% BIN_WIDTH)

  posterior <- as.data.frame(result$param)
  colnames(posterior) <- PARAM_NAMES
  weights <- result$weights
  if (is.null(weights) || length(weights) != nrow(posterior)) {
    warning("[", scn$label, "] No usable result$weights; sampling particles uniformly.")
    weights <- rep(1, nrow(posterior))
  }

  cat(sprintf("\n[%s] Loaded %d posterior particles from:\n  %s\n",
              scn$label, nrow(posterior), rds_path))
  cat(sprintf("[%s] Scenario id: %s | fitted params: %s\n",
              scn$label, SCENARIO_ID, paste(PARAM_NAMES, collapse = ", ")))
  if (!is.null(meta$observed_summaries)) {
    os <- meta$observed_summaries
    cat(sprintf("[%s] ABC targets (mean over taken-off reps): total ~ %.0f, HCW ~ %.0f deaths.\n",
                scn$label, exp(os[["log_n_deaths"]]), exp(os[["log_n_hcw_deaths"]])))
  }

  ## ---- base / time-varying args + R0 invariants (must match the fit) ----
  mp <- make_model_parameters(scenario_id = SCENARIO_ID, scenario_matrix = scenario_matrix,
                              overrides = MODEL_OVERRIDES)
  base_args     <- mp$base_args
  tv_args_model <- mp$tv_args
  R0_invariants <- compute_R0_invariants(args = mp$args, n = SETUP_R0_N, seed = SETUP_R0_SEED)

  ## ---- draw posterior samples + build the (draw, rep) job list ----
  set.seed(SAMPLE_SEED)
  draw_idx <- sample(seq_len(nrow(posterior)), size = N_POST, replace = TRUE, prob = weights)

  jobs <- vector("list", N_POST * N_REPS)
  jj <- 1L
  for (k in seq_len(N_POST)) {
    for (r in seq_len(N_REPS)) {
      seed_kr <- BASE_SEED + (k - 1L) * N_REPS + r
      jobs[[jj]] <- list(draw = k, draw_row = draw_idx[k], rep = r, seed = seed_kr)
      jj <- jj + 1L
    }
  }

  ## ---- worker: ONE OBV run; returns realised + prevented death day-indices ----
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

    ## OBV PEP is ALWAYS on: one run yields the realised epidemic AND its
    ## internal counterfactual (out$prevented_completed).
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

    ## WITH OBV: realised deaths (time_outcome_absolute where outcome == TRUE).
    tdf  <- out$tdf
    tdf  <- tdf[!is.na(tdf$time_infection_absolute), , drop = FALSE]
    rdied <- !is.na(tdf$outcome) & tdf$outcome
    rhcw  <- !is.na(tdf$class)  & tdf$class == "HCW"
    with_total <- fdays(tdf$time_outcome_absolute[rdied])
    with_hcw   <- fdays(tdf$time_outcome_absolute[rdied & rhcw])

    ## Deaths OBV PREVENTED: the averted INDEX infections, replayed through the
    ## same outcome model (out$prevented_completed; NULL if nothing prevented).
    ## Each carries its would-be class + death time, so we can add them back.
    pc <- out$prevented_completed
    if (is.null(pc) || nrow(pc) == 0L) {
      prev_total <- integer(0); prev_hcw <- integer(0)
    } else {
      pdied <- !is.na(pc$outcome) & pc$outcome
      phcw  <- !is.na(pc$class)  & pc$class == "HCW"
      prev_total <- fdays(pc$time_outcome_absolute[pdied])
      prev_hcw   <- fdays(pc$time_outcome_absolute[pdied & phcw])
    }

    ## WITHOUT OBV (direct within-run counterfactual) = realised + prevented index.
    without_total <- c(with_total, prev_total)
    without_hcw   <- c(with_hcw,   prev_hcw)

    list(draw = job$draw, rep = job$rep,
         with_total = with_total, with_hcw = with_hcw,
         without_total = without_total, without_hcw = without_hcw,
         n_total = length(without_total))   # take-off on the counterfactual (no-OBV) scale
  }

  ## ---- run all replicates in parallel (uses the active multisession plan) ----
  cat(sprintf("[%s] Running %d OBV simulations (counterfactual from prevented_completed)...\n",
              scn$label, length(jobs)))
  st <- Sys.time()
  results <- future_lapply(jobs, run_traj_job,
                           future.packages = "fiber", future.seed = TRUE)
  cat(sprintf("[%s] Simulation time: ", scn$label)); print(Sys.time() - st)

  ## ---- bin to weekly deaths (+ cumulative); fill BOTH arms from each run ----
  max_day <- max(c(unlist(lapply(results, function(r) r$without_total)), 0L))
  n_bins  <- max(1L, ceiling((max_day + 1L) / bin_width))
  week    <- ((seq_len(n_bins) - 1L) + 0.5) * bin_width

  mk_arr <- function() array(0L, dim = c(N_POST, N_REPS, n_bins))
  A <- setNames(lapply(names(ARM_LABELS), function(ak)
    list(hcw  = list(inc = mk_arr(), cum = mk_arr()),
         tot  = list(inc = mk_arr(), cum = mk_arr()),
         took = matrix(FALSE, N_POST, N_REPS))), names(ARM_LABELS))

  for (res in results) {
    d <- res$draw; rp <- res$rep
    took <- res$n_total >= TAKEOFF_DEATH_THRESHOLD   # one run -> both arms share take-off

    ## WITH OBV ("on") = realised
    bh <- bin_counts(res$with_hcw,   bin_width, n_bins)
    bt <- bin_counts(res$with_total, bin_width, n_bins)
    A[["on"]]$hcw$inc[d, rp, ] <- bh; A[["on"]]$hcw$cum[d, rp, ] <- cumsum(bh)
    A[["on"]]$tot$inc[d, rp, ] <- bt; A[["on"]]$tot$cum[d, rp, ] <- cumsum(bt)
    A[["on"]]$took[d, rp] <- took

    ## WITHOUT OBV ("off") = realised + prevented index (direct counterfactual)
    bh <- bin_counts(res$without_hcw,   bin_width, n_bins)
    bt <- bin_counts(res$without_total, bin_width, n_bins)
    A[["off"]]$hcw$inc[d, rp, ] <- bh; A[["off"]]$hcw$cum[d, rp, ] <- cumsum(bh)
    A[["off"]]$tot$inc[d, rp, ] <- bt; A[["off"]]$tot$cum[d, rp, ] <- cumsum(bt)
    A[["off"]]$took[d, rp] <- took
  }

  ## ---- magnitude check (mean over taken-off reps, per arm) ----
  for (ak in names(ARM_LABELS)) {
    tk <- A[[ak]]$took
    ft <- A[[ak]]$tot$cum[, , n_bins]; fh <- A[[ak]]$hcw$cum[, , n_bins]
    cat(sprintf("[%s] %-12s: take-off %.1f%% of reps; mean over taken-off reps -> total = %.0f, HCW = %.1f.\n",
                scn$label, ARM_LABELS[[ak]], 100 * mean(tk),
                if (any(tk)) mean(ft[tk]) else NA_real_,
                if (any(tk)) mean(fh[tk]) else NA_real_))
  }

  ## ---- per-draw central lines -> long form ----
  draws_list <- list(); li <- 1L
  for (ak in names(ARM_LABELS)) for (mk in names(METRIC_KEYS)) for (meas in names(MEASURES_K)) {
    mt <- per_draw_central(A[[ak]][[mk]][[meas]], A[[ak]]$took)   # [N_POST x n_bins]
    draws_list[[li]] <- data.frame(
      scenario = scn$label,
      draw     = rep(seq_len(N_POST), times = n_bins),
      week     = rep(week, each = N_POST),
      value    = as.vector(mt),
      arm      = ARM_LABELS[[ak]], metric = METRIC_KEYS[[mk]], measure = MEASURES_K[[meas]],
      stringsAsFactors = FALSE
    )
    li <- li + 1L
  }
  do.call(rbind, draws_list)
}


# -----------------------------------------------------------------------------
# 5. RUN BOTH SCENARIOS
# -----------------------------------------------------------------------------
N_CLUSTER <- if (parallel::detectCores() > 120) {
  min(120L, parallel::detectCores() - 10L)
} else {
  min(12L, parallel::detectCores() - 4L)
}
N_WORKERS <- N_WORKERS_OVERRIDE %||% max(1L, min(N_CLUSTER, future::availableCores()))

cat(sprintf("OBV: efficacy %.0f%%, coverage %.1f, adherence %.1f, target %s/%s\n",
            100 * OBV$efficacy, OBV$coverage, OBV$adherence, OBV$target_class, OBV$target_locations))
cat(sprintf("Per scenario: %d draws x %d reps = %d single OBV runs, on %d workers.\n",
            N_POST, N_REPS, N_POST * N_REPS, N_WORKERS))

plan(multisession, workers = N_WORKERS)
draws_long <- do.call(rbind, lapply(SCENARIOS, compute_scenario_draws))
plan(sequential)

draws_long$scenario <- factor(draws_long$scenario, levels = SCEN_LEVELS)


# -----------------------------------------------------------------------------
# 6. CROSS-DRAW BANDS  (central line + IQR + 95% CrI, per scenario/metric/measure/week)
# -----------------------------------------------------------------------------
# `draws_long`: one row per scenario x draw x week (value = that draw's central
# across the N_REPS reps). Summarise ACROSS the draws into a band per cell.
grp <- with(draws_long, interaction(scenario, arm, metric, measure, week,
                                    drop = TRUE, sep = "\r"))
band95 <- do.call(rbind, lapply(split(draws_long, grp), function(d) data.frame(
  scenario = as.character(d$scenario[1]),
  arm      = as.character(d$arm[1]),
  metric   = as.character(d$metric[1]),
  measure  = as.character(d$measure[1]),
  week     = d$week[1],
  med      = central_fun(d$value),
  q25      = qf(d$value, 0.25),  q75  = qf(d$value, 0.75),    # IQR
  lo95     = qf(d$value, 0.025), hi95 = qf(d$value, 0.975),   # 95% CrI
  stringsAsFactors = FALSE
)))
band95$scenario <- factor(band95$scenario, levels = SCEN_LEVELS)
band95$arm      <- factor(band95$arm,      levels = ARM_LEVELS)
band95$measure  <- factor(band95$measure,  levels = MEASURE_LEVELS)

## Headline: cross-draw central cumulative deaths at the final week, per scenario.
cat(sprintf("\n==== OBV direct effect (cross-draw %s cumulative deaths at final week) ====\n", BAND_CENTRAL))
for (sc in SCEN_LEVELS) for (mlab in unname(METRIC_KEYS)) {
  d <- band95[band95$scenario == sc & band95$metric == mlab & band95$measure == "Cumulative", ]
  if (!nrow(d)) next
  fw <- max(d$week)
  no <- d$med[d$arm == "Without OBV" & d$week == fw]; no <- if (length(no)) no[1] else NA_real_
  ob <- d$med[d$arm == "With OBV"    & d$week == fw]; ob <- if (length(ob)) ob[1] else NA_real_
  cat(sprintf("  %-12s | %-12s: without OBV = %7.1f, with OBV = %7.1f (averted %.1f%%)\n",
              sc, mlab, no, ob, if (is.finite(no) && no > 0) 100 * (no - ob) / no else NA_real_))
}


# -----------------------------------------------------------------------------
# 7. PER-SCENARIO BAND FIGURES (raw ggplot; weekly incidence on top, cumulative below)
# -----------------------------------------------------------------------------
make_band_fig <- function(metric_label, scen_label) {
  d <- band95[band95$metric == metric_label & band95$scenario == scen_label, , drop = FALSE]
  ggplot(d, aes(week, med, colour = arm, fill = arm)) +
    geom_ribbon(aes(ymin = lo95, ymax = hi95), alpha = CRI_ALPHA, colour = NA) +   # 95% CrI (lighter)
    geom_ribbon(aes(ymin = q25,  ymax = q75),  alpha = IQR_ALPHA, colour = NA) +   # IQR (darker)
    geom_line(linewidth = BAND_LWD, na.rm = TRUE) +
    facet_wrap(~ measure, ncol = 1, scales = "free_y") +
    scale_colour_manual(values = ARM_COLS) +
    scale_fill_manual(values = ARM_COLS) +
    coord_cartesian(xlim = c(0, XLIM)) +
    labs(x = "Days since outbreak seeding", y = metric_label,
         colour = NULL, fill = NULL,
         title = sprintf("%s decoupled fit: %s, obeldesivir effect",
                         scen_label, tolower(metric_label)),
         subtitle = sprintf(paste0("Single OBV run per draw x rep; 'Without OBV' = realised + prevented_completed ",
                                   "(direct within-run counterfactual).  Line: cross-draw %s; dark band: IQR; light band: 95%% CrI."),
                            BAND_CENTRAL)) +
    theme_bw(base_size = 12) +
    theme(plot.subtitle = element_text(size = 7), legend.position = "top")
}

for (sc in SCEN_LEVELS) {
  print(make_band_fig("HCW deaths",   sc))
  print(make_band_fig("Total deaths", sc))
}


# -----------------------------------------------------------------------------
# 8. FINAL 2x2 PANEL: HCW deaths -- weekly incidence (top) / cumulative (bottom),
#    one column per scenario, rel_heights = c(1, 2).
# -----------------------------------------------------------------------------
# One un-faceted cell. Top-row cells carry the scenario title and hide the x-axis
# (shared with the cumulative cell beneath); only the left column keeps a y title.
make_cell <- function(scen_label, measure_label, y_lab = NULL, top = FALSE) {
  d <- band95[band95$metric == "HCW deaths" &
              band95$scenario == scen_label &
              band95$measure  == measure_label, , drop = FALSE]
  g <- ggplot(d, aes(week, med, colour = arm, fill = arm)) +
    geom_ribbon(aes(ymin = lo95, ymax = hi95), alpha = CRI_ALPHA, colour = NA) +
    geom_ribbon(aes(ymin = q25,  ymax = q75),  alpha = IQR_ALPHA, colour = NA) +
    geom_line(linewidth = BAND_LWD, na.rm = TRUE) +
    scale_colour_manual(values = ARM_COLS) +
    scale_fill_manual(values = ARM_COLS) +
    coord_cartesian(xlim = c(0, XLIM)) +
    labs(x = "Days since outbreak seeding", y = y_lab,
         colour = NULL, fill = NULL,
         title = if (top) scen_label else NULL) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold"))
  if (top) g <- g + theme(axis.title.x = element_blank(),
                          axis.text.x  = element_blank(),
                          axis.ticks.x = element_blank())
  g
}

# Cells listed row-major: top row (incidence) then bottom row (cumulative).
cells <- list(
  make_cell("DRC",         "Weekly incidence", y_lab = "Weekly HCW deaths",     top = TRUE),
  make_cell("West Africa", "Weekly incidence", y_lab = NULL,                    top = TRUE),
  make_cell("DRC",         "Cumulative",       y_lab = "Cumulative HCW deaths", top = FALSE),
  make_cell("West Africa", "Cumulative",       y_lab = NULL,                    top = FALSE)
)

# The 2x2 grid exactly as requested (incidence row over cumulative row).
grid_2x2 <- plot_grid(plotlist = cells, nrow = 2, ncol = 2,
                      rel_heights = c(1, 2), align = "hv", axis = "tblr")

# Shared arm legend on top (kept as a separate strip so the 2x2 itself keeps
# rel_heights = c(1, 2)). If your cowplot/ggplot2 versions return an empty
# legend here, swap get_legend(...) for:
#   get_plot_component(make_band_fig("HCW deaths","DRC"), "guide-box-top")
shared_legend <- get_legend(make_band_fig("HCW deaths", "DRC") +
                              theme(legend.position = "top"))
panel_hcw <- plot_grid(shared_legend, grid_2x2, ncol = 1, rel_heights = c(0.08, 1))

print(panel_hcw)