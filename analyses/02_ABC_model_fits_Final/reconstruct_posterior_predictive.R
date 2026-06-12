# reconstruct_posterior_predictive.R
# =============================================================================
# Read the two FINAL decoupled ABC-SMC fits, reconstruct posterior trajectories,
# and draw posterior-predictive checks (summary-stat histograms + fit-ratio), a
# la sections 9-11 of the *_run_abc_calibration_decoupled.R scripts -- but driven
# straight off the saved .RDS (each file IS the `result`, with $param/$stats/
# $weights and a self-describing $run_metadata).
#
#   STAGE A (always; needs only base R): posterior-predictive summary statistics
#     - per-stat histograms with the 2.5/50/97.5% posterior-predictive lines and
#       the observed target overlaid;
#     - a single-panel "simulated / observed" fit-ratio plot.
#   STAGE B (DO_TRAJECTORIES; needs the `fiber` package + the functions/): re-
#     simulate N_TRAJ posterior draws and draw spaghetti + median/95% CrI bands.
#
# Observed targets + stat names are pulled from each file's $run_metadata, so the
# only per-scenario inputs are the .RDS, the scenario CSV, and check_final_size.
# =============================================================================

suppressPackageStartupMessages({ library(here); library(ggplot2) })
`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- the two fits (as supplied) --------------------------------------------
SCENARIOS <- list(
  WestAfrica = list(
    id           = "Worst_WestAfrica",
    rds          = here("outputs", "02_ABC_model_fits_Final",
                        "fiber_ABC_SMC_Worst_WestAfrica_Decoupled_20260608_162044_check_NP5_NS4_NBREPS_30_NBSIMUL_472.RDS"),
    scenario_csv = here("data-processed", "final_six_scenario_values_original_approach.csv"),
    param_names  = c("R0", "prop_funeral", "etu_efficacy", "ppe_efficacy", "hcw_risk_scalar"),
    check_final_size = 40000
  ),
  DRC = list(
    id           = "Middle_DRC_ConflictSmoothed_PlusPlus",
    rds          = here("outputs", "02_ABC_model_fits_Final",
                        "fiber_ABC_SMC_Middle_DRC_ConflictSmoothed_PlusPlus_Decoupled_20260607_215621_check_NP5_NS4_NBREPS_30_NBSIMUL_590.RDS"),
    scenario_csv = here("data-processed", "final_six_scenario_values_original_approach.csv"),
    param_names  = c("R0", "prop_funeral", "etu_efficacy", "ppe_efficacy", "hcw_risk_scalar"),
    check_final_size = 10000
  )
)

# ---- knobs ------------------------------------------------------------------
DO_TRAJECTORIES <- TRUE       # STAGE B: re-simulate (needs `fiber`); the slow part
SAVE_PDF        <- FALSE      # FALSE = plot to the active device; TRUE = one PDF per scenario
N_POST          <- 10000L     # posterior-predictive resample size (summary stats)
N_TRAJ          <- 200L       # trajectories re-simulated per scenario (lower = faster)
BIN_WIDTH       <- 7L         # weekly bins (matches PEAK_BIN_WIDTH in the fits)
SEEDING_CASES   <- 25L        # matches SEEDING_CASES in the fits (not stored in metadata)
SETUP_R0_N      <- 100000L    # R0-invariant calibration draws (matches the fits)
SETUP_R0_SEED   <- 42L
PP_SEED         <- 1L
TRAJ_SEED       <- 100L
TRAJ_WORKERS    <- max(1L, min(parallel::detectCores() - 1L, 8L))

FUNCTIONS_DIR <- here("functions")
OUT_DIR       <- here("outputs", "02_ABC_model_fits_Final")

# =============================================================================
# STAGE A -- posterior predictive on the fitted summary statistics
# =============================================================================
load_fit <- function(sc) {
  stopifnot(file.exists(sc$rds))
  res <- readRDS(sc$rds)
  md  <- res$run_metadata %||% list()
  stat_names <- md$summary_stats %||% colnames(res$stats)
  observed   <- md$observed_summaries
  if (is.null(stat_names) || is.null(observed))
    stop("No summary_stats/observed_summaries in $run_metadata for ", sc$id,
         " -- supply them in SCENARIOS to proceed.")
  stats <- as.data.frame(res$stats); colnames(stats) <- stat_names
  set.seed(PP_SEED)
  idx  <- sample(seq_len(nrow(stats)), size = N_POST, replace = TRUE, prob = res$weights)
  list(res = res, md = md, stat_names = stat_names, observed = observed,
       post = stats[idx, , drop = FALSE])
}

# per-stat histograms: PP distribution + 2.5/50/97.5% lines (blue) + observed (red)
plot_pp_hist <- function(name, fit) {
  obs <- fit$observed; ns <- length(obs)
  op <- par(mfrow = c(ceiling(ns / 2), 2), mar = c(4, 4, 3, 1)); on.exit(par(op))
  for (s in names(obs)) {
    x <- fit$post[[s]]
    hist(x, breaks = 12, main = paste0(name, "  |  PP: ", s), xlab = s,
         col = adjustcolor("steelblue", 0.6), border = "white")
    abline(v = quantile(x, c(0.025, 0.5, 0.975)), col = "darkblue",
           lty = c(2, 1, 2), lwd = c(1, 2, 1))
    abline(v = obs[s], col = "red", lwd = 2.5)
  }
}

# single panel: simulated / observed for every fitted summary (fitted scale, so
# for the log_* summaries this is a ratio of logs -- ~1 when well fit).
plot_fit_ratio <- function(name, fit) {
  obs <- fit$observed; ns <- length(obs)
  qs <- sapply(names(obs), function(s) quantile(fit$post[[s]] / obs[s], c(0.025, 0.5, 0.975)))
  op <- par(mfrow = c(1, 1), mar = c(8, 4, 3, 1)); on.exit(par(op))
  plot(NA, xlim = c(0.5, ns + 0.5), ylim = c(min(0, qs), max(1.5, qs) * 1.02),
       xaxt = "n", xlab = "", ylab = "Simulated / Observed",
       main = paste0(name, ": posterior-predictive fit ratio"))
  axis(1, at = seq_len(ns), labels = names(obs), las = 2, cex.axis = 0.8)
  abline(h = 1, lty = 2, col = "red")
  for (i in seq_len(ns)) {
    segments(i, qs[1, i], i, qs[3, i], lwd = 2, col = "darkblue")
    points(i, qs[2, i], pch = 16, cex = 1.5, col = "darkblue")
  }
}

# =============================================================================
# STAGE B -- reconstruct posterior trajectories by re-simulating the model
# =============================================================================
reconstruct_trajectories <- function(sc, fit) {
  suppressPackageStartupMessages({
    library(fiber); library(future); library(future.apply)
  })
  source(file.path(FUNCTIONS_DIR, "setup_model_parameters.R"))
  source(file.path(FUNCTIONS_DIR, "abc_calibration_functions_common.R"))
  source(file.path(FUNCTIONS_DIR, "abc_calibration_functions_decoupled.R"))
  source(file.path(FUNCTIONS_DIR, "calculate_model_approx_r0.R"))
  source(file.path(FUNCTIONS_DIR, "simulation_helpers.R"))   # bin_counts()

  res <- fit$res; md <- fit$md
  fit_params   <- md$fit_params %||% sc$param_names
  fixed_values <- md$fixed_values
  ghqe <- md$fixed_efficacies$general_hospital_quarantine_efficacy
  sfe  <- md$fixed_efficacies$safe_funeral_efficacy
  hcw_base_prob <- md$hcw_base_prob %||% 0.25

  # rebuild the scenario args exactly as the calibration script's section 4
  scenario_matrix <- read_scenario_matrix(sc$scenario_csv)
  mp <- make_model_parameters(scenario_id = sc$id, scenario_matrix = scenario_matrix,
                              overrides = list(check_final_size = sc$check_final_size))
  R0_invariants <- compute_R0_invariants(args = mp$args, n = SETUP_R0_N, seed = SETUP_R0_SEED)

  posterior <- as.data.frame(res$param); colnames(posterior) <- fit_params
  set.seed(TRAJ_SEED)
  traj_idx <- sample(seq_len(nrow(posterior)), size = N_TRAJ, replace = TRUE, prob = res$weights)

  plan(multisession, workers = TRAJ_WORKERS); on.exit(plan(sequential), add = TRUE)
  traj_runs <- future_lapply(seq_along(traj_idx), function(k) {
    i    <- traj_idx[k]
    full <- assemble_decoupled_theta(as.numeric(posterior[i, fit_params]), fit_params, fixed_values)
    a <- build_abc_model_args_decoupled(
      R0 = full$R0, prop_funeral = full$prop_funeral, etu_efficacy = full$etu_efficacy,
      ppe_efficacy = full$ppe_efficacy, hcw_risk_scalar = full$hcw_risk_scalar,
      base = mp$base_args, tv = mp$tv_args, invariants = R0_invariants,
      general_hospital_quarantine_efficacy = ghqe, safe_funeral_efficacy = sfe,
      hcw_base_prob = hcw_base_prob, seeding_cases = SEEDING_CASES)
    a$seed <- TRAJ_SEED + k
    tdf  <- do.call(branching_process_main, a)$tdf
    tdf  <- tdf[!is.na(tdf$time_infection_absolute), , drop = FALSE]
    died <- !is.na(tdf$outcome) & tdf$outcome
    list(case_days  = as.integer(floor(tdf$time_infection_absolute)),
         death_days = as.integer(floor(tdf$time_outcome_absolute[died])))
  }, future.packages = "fiber", future.seed = TRUE)

  max_day <- max(unlist(lapply(traj_runs, function(r) c(r$case_days, r$death_days))), 0)
  n_bins  <- max(1L, ceiling((max_day + 1L) / BIN_WIDTH))
  week    <- ((seq_len(n_bins) - 1L) + 0.5) * BIN_WIDTH
  do.call(rbind, lapply(seq_along(traj_runs), function(k) {
    r <- traj_runs[[k]]
    rbind(data.frame(scenario = sc$id, draw = k, week = week, metric = "Cases",
                     incidence = bin_counts(r$case_days,  BIN_WIDTH, n_bins)),
          data.frame(scenario = sc$id, draw = k, week = week, metric = "Deaths",
                     incidence = bin_counts(r$death_days, BIN_WIDTH, n_bins)))
  }))
}

plot_trajectories <- function(name, traj_long) {
  print(ggplot(traj_long, aes(week, incidence, group = draw)) +
          geom_line(alpha = 0.15, colour = "steelblue", linewidth = 0.3) +
          facet_wrap(~ metric, scales = "free_y") +
          labs(x = "Time since outbreak start (days)", y = "Count per week",
               title = sprintf("Posterior trajectories: %s", name)) + theme_bw())
  band <- do.call(rbind, lapply(split(traj_long, list(traj_long$metric, traj_long$week), drop = TRUE),
    function(d) data.frame(metric = d$metric[1], week = d$week[1],
                           lo  = quantile(d$incidence, 0.025, names = FALSE),
                           med = quantile(d$incidence, 0.500, names = FALSE),
                           hi  = quantile(d$incidence, 0.975, names = FALSE))))
  print(ggplot(band, aes(week, med)) +
          geom_ribbon(aes(ymin = lo, ymax = hi), fill = "steelblue", alpha = 0.3) +
          geom_line(colour = "steelblue", linewidth = 0.8) +
          facet_wrap(~ metric, scales = "free_y") +
          labs(x = "Time since outbreak start (days)", y = "Count per week",
               title = sprintf("Posterior median + 95%% CrI: %s", name)) + theme_bw())
}

# =============================================================================
# DRIVER -- loop the scenarios
# =============================================================================
for (name in names(SCENARIOS)) {
  sc  <- SCENARIOS[[name]]
  message("\n==== ", name, "  (", sc$id, ") ====")
  fit <- load_fit(sc)

  cat(sprintf("  particles: %d   summaries: %s\n",
              nrow(fit$post), paste(fit$stat_names, collapse = ", ")))

  traj_long <- if (DO_TRAJECTORIES) {
    message("  re-simulating ", N_TRAJ, " posterior trajectories ...")
    reconstruct_trajectories(sc, fit)
  } else NULL

  draw_all <- function() {
    plot_pp_hist(name, fit)
    plot_fit_ratio(name, fit)
    if (!is.null(traj_long)) plot_trajectories(name, traj_long)
  }

  if (SAVE_PDF) {
    pdf_path <- file.path(OUT_DIR, sprintf("posterior_predictive_%s.pdf", name))
    pdf(pdf_path, width = 9, height = 7); draw_all(); dev.off()
    message("  wrote ", pdf_path)
  } else {
    draw_all()
  }
}
