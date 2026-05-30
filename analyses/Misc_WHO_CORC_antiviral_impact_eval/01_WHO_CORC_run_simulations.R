# 01_run_simulations.R  (analyses/04_obeldesivir_impact)
# =============================================================================
# Obeldesivir (OBV PEP) impact analysis — simulation step.
#
# What this script does, end to end:
#   1. Reads the posterior of the 3 fitted parameters (R0, prop_funeral,
#      hcw_risk_scalar) from the DRC ABC-SMC run in analyses/02_ABC_model_fits_HCWrisk,
#      using the final completed step (step 7 by default).
#   2. Converts the 3 fitted parameters into the 4 fiber model parameters
#      (mn_offspring_genPop, mn_offspring_funeral, and the two
#      prob_hcw_cond_*_hospital probabilities), using exactly the calibration's
#      mapping (build_abc_model_args()).
#   3. Downsamples the weighted posterior to N_SETS = 100 parameter sets and
#      runs N_REPS = 5 stochastic replicates for each set.
#   4. Runs the fiber branching-process model twice for every (set, rep):
#         - WITHOUT obeldesivir
#         - WITH    obeldesivir (80% efficacy, 100% coverage, 100% adherence)
#      Work is parallelised with the `future` framework (multisession backend,
#      Windows-compatible) across N_WORKERS cores.
#   5. Computes (a) the total number of deaths and the total number of HCW
#      deaths for each arm, and (b) the deaths averted by obeldesivir (total and
#      HCW), then saves:
#         - per-replicate results,
#         - per-set mean epidemic-curve trajectories (deaths/week and HCW
#           deaths/week), for the plotting script,
#         - headline summary tables (CSV) in outputs/.
#
# Run it with your working directory set to either the repo root or this folder.
# Heavy compute lives here; 02_plot.R only reads the saved intermediate.
# =============================================================================


# -----------------------------------------------------------------------------
# 0. Locate the repo root (works however you launch the script)
# -----------------------------------------------------------------------------
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

# Path to the repo root (the folder containing obv_hcw_paper.Rproj). Known
# machines are baked in (mirroring analyses/02_ABC_model_fits_HCWrisk/DRC_run_abc_calibration.R);
# on any other machine we fall back to auto-detection. If neither works, just set
# REPO_ROOT directly below, e.g.
#   REPO_ROOT <- "C:/Users/cwhittaker/Documents/Research Projects/obv_hcw_paper"
REPO_ROOT <- switch(
  Sys.info()[["user"]],
  "cwhittaker" = "C:/Users/cwhittaker/Documents/Research Projects/obv_hcw_paper",
  "PETAL_WS_2" = "C:/Users/PETAL_WS_2/Documents/obv_hcw_paper",
  "PETAL_WS_1" = "C:/Users/PETAL_WS_1/Documents/obv_hcw_paper",
  NA_character_
)

# Fallback: walk up from the script location (Rscript / RStudio "Source") or,
# failing that, the current working directory, until the .Rproj is found.
if (is.na(REPO_ROOT)) {
  get_script_dir <- function() {
    args     <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)))
    for (i in rev(seq_len(sys.nframe()))) {
      of <- sys.frame(i)$ofile
      if (!is.null(of)) return(dirname(normalizePath(of, mustWork = FALSE)))
    }
    NULL
  }
  find_repo_root_local <- function(start, marker = "obv_hcw_paper.Rproj") {
    d <- normalizePath(start, winslash = "/", mustWork = FALSE)
    repeat {
      if (file.exists(file.path(d, marker))) return(d)
      parent <- dirname(d)
      if (identical(parent, d)) return(NULL)
      d <- parent
    }
  }
  REPO_ROOT <- find_repo_root_local(get_script_dir() %||% getwd()) %||%
    find_repo_root_local(getwd())
}

if (is.null(REPO_ROOT) || is.na(REPO_ROOT) ||
    !file.exists(file.path(REPO_ROOT, "obv_hcw_paper.Rproj"))) {
  stop("Could not locate the repo root. Set REPO_ROOT at the top of this script, e.g.\n",
       '  REPO_ROOT <- "C:/Users/cwhittaker/Documents/Research Projects/obv_hcw_paper"',
       call. = FALSE)
}

ANALYSIS_DIR <- file.path(REPO_ROOT, "analyses", "Misc_WHO_CORC_obeldesivir_impact_eval")


# -----------------------------------------------------------------------------
# 1. CONFIGURATION
# -----------------------------------------------------------------------------
SCENARIO_ID <- "Middle_DRC_ConflictSmoothed"

# Posterior downsampling and replication.
N_SETS <- 100L   # parameter sets drawn from the (weighted) posterior
N_REPS <- 5L     # stochastic replicates per parameter set per arm
ABC_STEP <- 7L   # which ABC step to read the posterior from (the final one); NULL = latest

# Obeldesivir PEP settings for the "with obeldesivir" arm.
#   80% efficacy, 100% coverage, 100% adherence.
# By default obeldesivir is modelled (as in fiber) as post-exposure prophylaxis
# for HCWs exposed in the hospital setting (obv_pep_target_class = "HCW",
# obv_pep_target_locations = "hospital").
OBV_CFG <- list(
  efficacy         = 0.80,
  coverage         = 1.00,
  adherence        = 1.00,
  dpc              = 1,            # days post-exposure to first dose (unused when efficacy is a scalar)
  target_class     = "HCW",
  target_locations = "hospital"
)

# Parallelisation. The user runs this on Windows with 10 cores available;
# `multisession` is the Windows-compatible future backend. Capped at the number
# of cores actually available so the script is portable.
N_WORKERS <- 10L

# Model / calibration settings. These reproduce the settings the DRC ABC used so
# the posterior -> model-parameter conversion is identical to the calibration.
HCW_BASE_PROB    <- 0.25       # symmetric base for prob_hcw_cond_*_hospital
SEEDING_CASES    <- 25L        # seeding cases used during the ABC
CHECK_FINAL_SIZE <- 10000L     # MODEL_OVERRIDES$check_final_size in the DRC run
SETUP_FUNERAL_SHARE <- 0.5     # used only to pin the D / F units (as in the ABC)
SETUP_R0_N    <- 100000L
SETUP_R0_SEED <- 42L

# Seeds.
RESAMPLE_SEED <- 1L            # reproducible posterior downsampling
SEED_BASE     <- 20260528L     # base for per-(set, rep) simulation seeds

# Epidemic-curve binning (days per bin). 7 = weekly, which is the standard,
# readable unit for a multi-month outbreak; set to 1 for daily incidence.
BIN_WIDTH_DAYS <- 7L

# Output locations.
OUTPUT_DIR  <- file.path(REPO_ROOT, "outputs/misc")
RESULTS_RDS <- file.path(ANALYSIS_DIR, "WHO_CORC_prelim_obeldesivir_simulation_results.rds")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)


# -----------------------------------------------------------------------------
# 2. LIBRARIES + SOURCES
# -----------------------------------------------------------------------------
library(fiber)
library(future)
library(future.apply)
library(progressr)

# Calibration helpers (parameter setup, build_abc_model_args, R0 solver).
CALIB_HELPERS <- file.path(REPO_ROOT, "analyses", "02_ABC_model_fits_HCWrisk", "helper_functions")
source(file.path(CALIB_HELPERS, "setup_model_parameters.R"))
source(file.path(CALIB_HELPERS, "abc_calibration_functions.R"))
source(file.path(CALIB_HELPERS, "calculate_model_approx_r0.R"))

# This analysis's helpers.
source(file.path(ANALYSIS_DIR, "helper_functions.R"))

handlers("progress")


# -----------------------------------------------------------------------------
# 3. BUILD SCENARIO MODEL ARGS + SOLVE FOR D / F MULTIPLIERS
# -----------------------------------------------------------------------------
# Reproduce section 3 of DRC_run_abc_calibration.R so the conversion of the
# fitted (R0, prop_funeral) into (mn_offspring_genPop, mn_offspring_funeral) is
# identical to the calibration.
SCENARIO_CSV <- file.path(REPO_ROOT, "analyses", "02_ABC_model_fits_HCWrisk", "final_four_scenario_values.csv")
scenario_matrix <- read_scenario_matrix(SCENARIO_CSV)

mp <- make_model_parameters(
  scenario_id     = SCENARIO_ID,
  scenario_matrix = scenario_matrix,
  overrides       = list(check_final_size = CHECK_FINAL_SIZE)
)
base_args <- mp$base_args
tv_args   <- mp$tv_args

setup_solve <- solve_offspring_means_for_R0(
  R0   = 1.0,
  args = mp$args,
  proportion_transmission_from_funerals = SETUP_FUNERAL_SHARE,
  n    = SETUP_R0_N,
  seed = SETUP_R0_SEED
)
D_direct_multiplier  <- setup_solve$D_direct_multiplier
F_funeral_multiplier <- setup_solve$F_funeral_multiplier
message(sprintf("Scenario %s: D = %.5f, F = %.5f",
                SCENARIO_ID, D_direct_multiplier, F_funeral_multiplier))


# -----------------------------------------------------------------------------
# 4. READ POSTERIOR (STEP 7) + DOWNSAMPLE + DERIVE MODEL PARAMETERS
# -----------------------------------------------------------------------------
ABC_OUTPUTS_DIR <- file.path(REPO_ROOT, "analyses", "02_ABC_model_fits_HCWrisk", "abc_outputs")
ABC_RUN_DIR     <- find_latest_abc_run_dir(ABC_OUTPUTS_DIR, SCENARIO_ID)
message("Reading ABC posterior from: ", ABC_RUN_DIR)

posterior <- read_abc_posterior_step(ABC_RUN_DIR, step = ABC_STEP)
message(sprintf("Posterior step %d: %d particles.", attr(posterior, "step"), nrow(posterior)))

# Weighted resample to N_SETS parameter sets.
theta <- downsample_posterior(posterior, n_sets = N_SETS, seed = RESAMPLE_SEED)

# Transparent record of the 4 derived model parameters for each set.
derived_params <- derive_model_parameters(
  theta, D = D_direct_multiplier, F_fun = F_funeral_multiplier, hcw_base_prob = HCW_BASE_PROB
)
write.csv(derived_params,
          file.path(OUTPUT_DIR, "obeldesivir_downsampled_posterior_parameters.csv"),
          row.names = FALSE)

# Pre-build one full fiber argument list per parameter set, in the main session,
# using the calibration's build_abc_model_args(). Doing the conversion here (not
# on the workers) keeps it provably identical to the calibration and means the
# workers only need fiber. seed is set per-replicate inside simulate_one().
model_args_list <- lapply(seq_len(N_SETS), function(i) {
  build_abc_model_args(
    R0              = theta$R0[i],
    prop_funeral    = theta$prop_funeral[i],
    hcw_risk_scalar = theta$hcw_risk_scalar[i],
    base            = base_args,
    tv              = tv_args,
    D               = D_direct_multiplier,
    F_fun           = F_funeral_multiplier,
    seeding_cases   = SEEDING_CASES,
    hcw_base_prob   = HCW_BASE_PROB
  )
})


# -----------------------------------------------------------------------------
# 5. BUILD THE JOB LIST (set x rep x arm) AND RUN IN PARALLEL WITH future
# -----------------------------------------------------------------------------
# One job = one fiber simulation. 100 sets x 5 reps x 2 arms = 1000 jobs.
# The same seed is reused across arms for a (set, rep) pair (paired comparison).
arms <- c("no_obv", "obv")
jobs <- list()
k <- 0L
for (s in seq_len(N_SETS)) {
  for (r in seq_len(N_REPS)) {
    seed_sr <- SEED_BASE + (s - 1L) * N_REPS + (r - 1L)
    for (arm in arms) {
      k <- k + 1L
      jobs[[k]] <- list(set_id = s, rep_id = r, arm = arm, seed = seed_sr)
    }
  }
}
message(sprintf("Total simulations to run: %d (%d sets x %d reps x %d arms).",
                length(jobs), N_SETS, N_REPS, length(arms)))

n_workers <- min(N_WORKERS, future::availableCores())
message(sprintf("Parallelising over %d worker(s) via future::multisession.", n_workers))
plan(multisession, workers = n_workers)
on.exit(plan(sequential), add = TRUE)

start_time <- Sys.time()
with_progress({
  p <- progressor(along = jobs)
  results <- future_lapply(
    jobs,
    function(job, args_list, obv_cfg) {
      out <- simulate_one(job, args_list = args_list, obv_cfg = obv_cfg)
      p()
      out
    },
    args_list       = model_args_list,
    obv_cfg         = OBV_CFG,
    future.packages = "fiber",
    future.seed     = TRUE
  )
})
elapsed <- Sys.time() - start_time
message(sprintf("Simulations complete in %.1f %s.", as.numeric(elapsed), attr(elapsed, "units")))

plan(sequential)


# -----------------------------------------------------------------------------
# 6. PER-REPLICATE SUMMARY TABLE
# -----------------------------------------------------------------------------
per_rep <- data.frame(
  set_id       = vapply(results, `[[`, integer(1), "set_id"),
  rep_id       = vapply(results, `[[`, integer(1), "rep_id"),
  arm          = vapply(results, `[[`, character(1), "arm"),
  seed         = vapply(results, `[[`, numeric(1), "seed"),
  n_cases      = vapply(results, `[[`, integer(1), "n_cases"),
  n_deaths     = vapply(results, `[[`, integer(1), "n_deaths"),
  n_hcw_deaths = vapply(results, `[[`, integer(1), "n_hcw_deaths"),
  stringsAsFactors = FALSE
)
# Attach the fitted + derived parameters for each set.
per_rep <- merge(per_rep, derived_params, by = "set_id", all.x = TRUE)
per_rep <- per_rep[order(per_rep$set_id, per_rep$rep_id, per_rep$arm), ]


# -----------------------------------------------------------------------------
# 7. EPIDEMIC-CURVE TRAJECTORIES
# -----------------------------------------------------------------------------
# Bin each replicate's deaths into BIN_WIDTH_DAYS-day bins, then summarise each
# parameter set's N_REPS replicates by their MEAN trajectory (per arm). The
# plotting script turns these per-set mean trajectories into (i) individual
# lines and (ii) a median + 25%/75% interval band.
all_days <- unlist(lapply(results, function(r) c(r$death_days, r$hcw_death_days)), use.names = FALSE)
max_day  <- if (length(all_days)) max(all_days) else 0L
n_bins   <- max(1L, ceiling((max_day + 1L) / BIN_WIDTH_DAYS))
bin_mid  <- ((seq_len(n_bins) - 1L) + 0.5) * BIN_WIDTH_DAYS   # bin midpoint (days)

# Accumulate summed weekly incidence per (set, arm) and the replicate count.
sum_all <- array(0, dim = c(n_bins, N_SETS, length(arms)))
sum_hcw <- array(0, dim = c(n_bins, N_SETS, length(arms)))
rep_n   <- matrix(0L, N_SETS, length(arms))
for (r in results) {
  ai <- match(r$arm, arms)
  sum_all[, r$set_id, ai] <- sum_all[, r$set_id, ai] + bin_counts(r$death_days,     BIN_WIDTH_DAYS, n_bins)
  sum_hcw[, r$set_id, ai] <- sum_hcw[, r$set_id, ai] + bin_counts(r$hcw_death_days, BIN_WIDTH_DAYS, n_bins)
  rep_n[r$set_id, ai]     <- rep_n[r$set_id, ai] + 1L
}

# Long data frame of per-set MEAN trajectories.
traj_list <- vector("list", length(arms) * N_SETS)
ix <- 0L
for (ai in seq_along(arms)) {
  for (s in seq_len(N_SETS)) {
    ix <- ix + 1L
    denom <- max(rep_n[s, ai], 1L)
    traj_list[[ix]] <- data.frame(
      set_id          = s,
      arm             = arms[ai],
      time_days       = bin_mid,
      deaths_per_bin  = sum_all[, s, ai] / denom,
      hcw_deaths_per_bin = sum_hcw[, s, ai] / denom,
      stringsAsFactors = FALSE
    )
  }
}
trajectories <- do.call(rbind, traj_list)


# -----------------------------------------------------------------------------
# 8. HEADLINE OUTPUTS: TOTAL DEATHS / HCW DEATHS, AND DEATHS AVERTED
# -----------------------------------------------------------------------------
# (a) Distribution of the total number of deaths and HCW deaths per arm
#     (across all 100 x 5 = 500 replicates per arm).
arm_totals <- do.call(rbind, lapply(arms, function(arm) {
  sub <- per_rep[per_rep$arm == arm, ]
  rbind(
    data.frame(arm = arm, outcome = "total_deaths", t(q_summary(sub$n_deaths))),
    data.frame(arm = arm, outcome = "hcw_deaths",   t(q_summary(sub$n_hcw_deaths)))
  )
}))
rownames(arm_totals) <- NULL

# (b) Deaths averted by obeldesivir. Pair at the parameter-set level: average the
#     N_REPS replicates within each arm, then avert = (no obeldesivir) - (with).
#     This yields one averted value (and % averted) per parameter set; we then
#     summarise across the 100 sets.
agg <- aggregate(cbind(n_deaths, n_hcw_deaths) ~ set_id + arm, data = per_rep, FUN = mean)
wide <- reshape(agg, idvar = "set_id", timevar = "arm", direction = "wide")
# columns: n_deaths.no_obv, n_deaths.obv, n_hcw_deaths.no_obv, n_hcw_deaths.obv
per_set_averted <- data.frame(
  set_id                 = wide$set_id,
  deaths_no_obv          = wide$n_deaths.no_obv,
  deaths_obv             = wide$n_deaths.obv,
  hcw_deaths_no_obv      = wide$n_hcw_deaths.no_obv,
  hcw_deaths_obv         = wide$n_hcw_deaths.obv,
  deaths_averted         = wide$n_deaths.no_obv      - wide$n_deaths.obv,
  hcw_deaths_averted     = wide$n_hcw_deaths.no_obv  - wide$n_hcw_deaths.obv,
  stringsAsFactors       = FALSE
)
per_set_averted$pct_deaths_averted     <- with(per_set_averted,
                                               ifelse(deaths_no_obv     > 0, 100 * deaths_averted     / deaths_no_obv,     NA_real_))
per_set_averted$pct_hcw_deaths_averted <- with(per_set_averted,
                                               ifelse(hcw_deaths_no_obv > 0, 100 * hcw_deaths_averted / hcw_deaths_no_obv, NA_real_))

averted_summary <- rbind(
  data.frame(quantity = "deaths_averted",         t(q_summary(per_set_averted$deaths_averted))),
  data.frame(quantity = "hcw_deaths_averted",     t(q_summary(per_set_averted$hcw_deaths_averted))),
  data.frame(quantity = "pct_deaths_averted",     t(q_summary(per_set_averted$pct_deaths_averted))),
  data.frame(quantity = "pct_hcw_deaths_averted", t(q_summary(per_set_averted$pct_hcw_deaths_averted)))
)
rownames(averted_summary) <- NULL

# Diagnostic: fraction of replicates that "took off" (>= 100 deaths), the same
# threshold the calibration used; lets the reader see how many reps fizzled.
takeoff_rate <- aggregate(n_deaths ~ arm, data = per_rep,
                          FUN = function(x) mean(x >= 100))
names(takeoff_rate)[2] <- "prop_takeoff"

# Save CSV tables.
write.csv(arm_totals,      file.path(OUTPUT_DIR, "obeldesivir_total_deaths_by_arm.csv"),      row.names = FALSE)
write.csv(averted_summary, file.path(OUTPUT_DIR, "obeldesivir_deaths_averted_summary.csv"),   row.names = FALSE)
write.csv(per_set_averted, file.path(OUTPUT_DIR, "obeldesivir_per_set_deaths_averted.csv"),   row.names = FALSE)


# -----------------------------------------------------------------------------
# 9. SAVE EVERYTHING THE PLOTTING SCRIPT NEEDS
# -----------------------------------------------------------------------------
saveRDS(
  list(
    config = list(
      scenario_id = SCENARIO_ID, n_sets = N_SETS, n_reps = N_REPS,
      abc_step = attr(posterior, "step"), abc_run_dir = ABC_RUN_DIR,
      obv = OBV_CFG, bin_width_days = BIN_WIDTH_DAYS, arms = arms,
      D = D_direct_multiplier, F = F_funeral_multiplier,
      seeding_cases = SEEDING_CASES, check_final_size = CHECK_FINAL_SIZE,
      hcw_base_prob = HCW_BASE_PROB
    ),
    theta            = theta,
    derived_params   = derived_params,
    per_rep          = per_rep,
    trajectories     = trajectories,
    arm_totals       = arm_totals,
    per_set_averted  = per_set_averted,
    averted_summary  = averted_summary,
    takeoff_rate     = takeoff_rate
  ),
  RESULTS_RDS
)
message("Saved simulation results to: ", RESULTS_RDS)


# -----------------------------------------------------------------------------
# 10. CONSOLE SUMMARY
# -----------------------------------------------------------------------------
fmt <- function(v) sprintf("%s [%s, %s]",
                           formatC(v["q50"], format = "f", digits = 0, big.mark = ","),
                           formatC(v["q25"], format = "f", digits = 0, big.mark = ","),
                           formatC(v["q75"], format = "f", digits = 0, big.mark = ","))
cat("\n================ Obeldesivir impact (median [IQR] across the posterior) ================\n")
for (arm in arms) {
  td <- arm_totals[arm_totals$arm == arm & arm_totals$outcome == "total_deaths", ]
  hd <- arm_totals[arm_totals$arm == arm & arm_totals$outcome == "hcw_deaths", ]
  lbl <- if (arm == "no_obv") "WITHOUT obeldesivir" else "WITH    obeldesivir"
  cat(sprintf("  %s : total deaths = %s | HCW deaths = %s\n", lbl,
              fmt(unlist(td[paste0("q", c(25, 50, 75))])),
              fmt(unlist(hd[paste0("q", c(25, 50, 75))]))))
}
pa  <- averted_summary[averted_summary$quantity == "pct_deaths_averted", ]
pha <- averted_summary[averted_summary$quantity == "pct_hcw_deaths_averted", ]
cat(sprintf("  Deaths averted        : %.1f%% [%.1f, %.1f]\n", pa[["q50"]],  pa[["q25"]],  pa[["q75"]]))
cat(sprintf("  HCW deaths averted    : %.1f%% [%.1f, %.1f]\n", pha[["q50"]], pha[["q25"]], pha[["q75"]]))
cat(sprintf("  Take-off rate (>=100 deaths): no_obv = %.1f%%, obv = %.1f%%\n",
            100 * takeoff_rate$prop_takeoff[takeoff_rate$arm == "no_obv"],
            100 * takeoff_rate$prop_takeoff[takeoff_rate$arm == "obv"]))
cat("========================================================================================\n\n")
cat("Next: run 02_plot.R to produce the figures.\n")
