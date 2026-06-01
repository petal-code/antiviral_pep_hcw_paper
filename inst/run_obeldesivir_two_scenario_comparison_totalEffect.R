# =============================================================================
# run_obeldesivir_two_scenario_comparison_totalEffect.R
#
# Estimates the % of healthcare-worker (HCW) deaths that obeldesivir (OBV)
# post-exposure prophylaxis averts, for BOTH fitted scenarios (DRC and West
# Africa) and at TWO OBV efficacies (40% and 80%), then draws one stratified plot.
#
# TOTAL vs DIRECT effect
# ----------------------
# The companion script ..._directEffect.R uses fiber's within-run counterfactual
# (prevented_deaths), which counts ONLY the would-be death of each prevented
# index HCW infection -- it excludes the onward transmission those infections
# would have seeded. That is a conservative DIRECT estimate.
#
# This script measures the FULL (TOTAL) effect, including averted onward chains,
# by actually simulating both worlds and differencing their MEANS:
#   * For each posterior draw, run many reps WITHOUT OBV and many reps WITH OBV.
#   * Average HCW deaths over reps in each arm  -> mu_no(draw), mu_obv(draw).
#   * Averted(draw)   = mu_no - mu_obv ;  %Averted = 100 * (mu_no - mu_obv)/mu_no.
#
# WHY MEANS (and why this avoids the earlier negative-% artefact): enabling OBV
# consumes extra RNG, so a WITH-OBV run and a same-seed WITHOUT-OBV run desync --
# you cannot pair them rep-to-rep. But the MEAN over reps estimates the arm's
# EXPECTED deaths, and expectations do not depend on coupling. So mu_no - mu_obv
# is an UNBIASED estimate of the true averted burden regardless of desync; the
# desync only adds variance, never a biased sign. Because OBV can only remove
# deaths, the true per-draw effect is >= 0 -- so with ENOUGH reps each estimate
# sits at a genuine non-negative number. The bimodal fizzle/takeoff split means
# "enough" is large: prefer N_REPS >= 50, not 5.
#
# Efficiency: the no-OBV arm does not depend on efficacy, so it is run ONCE per
# (scenario, draw, rep) and reused as the shared baseline for every efficacy.
#
# Defines NO new functions -- composes existing functions/ helpers with base R +
# ggplot. Run from the repo root:
#   Rscript inst/run_obeldesivir_two_scenario_comparison_totalEffect.R
# =============================================================================


# ---------------------------------------------------------------------------
# 0. CONFIGURATION (everything tweakable lives here)
# ---------------------------------------------------------------------------
# NOTE: differencing arm means needs many reps to settle (bimodal outbreak
# sizes). 50 is a sensible floor; raise for a smoother estimate.
N_REPS  <- 50L    # stochastic replicates per particle per arm
N_CORES <- 12L    # parallel workers (capped at availableCores() below)

# OBV efficacies to test. Add/remove freely.
OBV_EFFICACIES <- c(obv_40 = 0.40,
                    obv_80 = 0.80)
# Fixed OBV delivery settings shared by every efficacy: full coverage/adherence,
# modelled as PEP for HCWs exposed in hospital.
OBV_BASE <- list(coverage = 1.0, adherence = 1.0, dpc = 1,
                 target_class = "HCW", target_locations = "hospital")

# The two scenarios: a friendly label -> (posterior RDS in inst/, scenario id).
SCENARIOS <- list(
  DRC        = list(rds = "fiber_ABC_SMC_Middle_DRC_ConflictSmoothed_2026-05-26.rds",
                    id  = "Middle_DRC_ConflictSmoothed"),
  WestAfrica = list(rds = "fiber_ABC_SMC_Worst_WestAfrica_2026-05-26.rds",
                    id  = "Worst_WestAfrica")
)

# Calibration settings reproduced so the posterior -> model-args conversion is
# identical to how each ABC was run (see 01_WHO_CORC_run_simulations.R).
HCW_BASE_PROB <- 0.25 ## relates a fitted parameter scaling healthcare transmission;
                      ## will be replaced by a different fitted parameter in the final fits
SEEDING_CASES <- 25L
CHECK_FINAL_SIZE <- 10000L
SETUP_FUNERAL_SHARE <- 0.5
SETUP_R0_N <- 100000
SETUP_R0_SEED <- 42L
RESAMPLE_SEED <- 1L; SEED_BASE <- 20260601L


# ---------------------------------------------------------------------------
# 1. LIBRARIES + SHARED HELPERS
# ---------------------------------------------------------------------------
library(fiber); library(future); library(future.apply); library(progressr); library(ggplot2)
FN <- here::here("functions")
source(file.path(FN, "setup_model_parameters.R"))            # scenario -> model args
source(file.path(FN, "calculate_model_approx_r0.R"))         # D/F solver
source(file.path(FN, "abc_calibration_functions_common.R"))  # generic ABC helpers
source(file.path(FN, "abc_calibration_functions_hcwRisk.R")) # build_abc_model_args()
source(file.path(FN, "abc_posterior.R"))                     # downsample_posterior()
source(file.path(FN, "simulation_helpers.R"))                # simulate_one()
handlers("progress")
check_model_function_version()   # fail fast if fiber predates the NPI interface


# ---------------------------------------------------------------------------
# 2. READ BOTH POSTERIORS, THEN PICK THE COMMON (SMALLER) SAMPLE SIZE
# ---------------------------------------------------------------------------
# The inst/ files are EasyABC result objects: $param (3 columns: R0,
# prop_funeral, hcw_risk_scalar) and $weights. Rebuild the weighted-posterior
# data frame downsample_posterior() expects.
posteriors <- lapply(SCENARIOS, function(sc) {
  res <- readRDS(here::here("inst", sc$rds))
  data.frame(weight = res$weights,
             R0 = res$param[, 1],               ## fitted: baseline R0
             prop_funeral = res$param[, 2],     ## fitted: funeral transmission share
             hcw_risk_scalar = res$param[, 3])  ## fitted: HCW in-hospital transmission scaler
})
N_SETS <- min(50, ## lower this for a quicker run
              min(vapply(posteriors, nrow, integer(1))))   # the scenario with the fewest draws
message(sprintf("Particles available: %s; using N_SETS = %d.",
                paste(sprintf("%s=%d", names(posteriors), vapply(posteriors, nrow, integer(1))), collapse = ", "), N_SETS))


# ---------------------------------------------------------------------------
# 3. PER SCENARIO: DOWNSAMPLE, SOLVE D/F, PRE-BUILD ONE ARG LIST PER PARTICLE
# ---------------------------------------------------------------------------
# Build the fiber argument lists in the main session (not on the workers) so the
# conversion stays provably identical to the calibration; workers only run fiber.
scenario_matrix  <- read_scenario_matrix(here::here("data-processed", "final_four_scenario_values.csv"))
args_by_scenario <- list()
for (sc in names(SCENARIOS)) {
  mp <- make_model_parameters(SCENARIOS[[sc]]$id, scenario_matrix,
                              overrides = list(check_final_size = CHECK_FINAL_SIZE))
  # D and F are the t=0 direct/funeral R0 multipliers for THIS scenario.
  sol <- solve_offspring_means_for_R0(R0 = 1.0, args = mp$args,
                                      proportion_transmission_from_funerals = SETUP_FUNERAL_SHARE,
                                      n = SETUP_R0_N, seed = SETUP_R0_SEED)
  theta <- downsample_posterior(posteriors[[sc]], n_sets = N_SETS, seed = RESAMPLE_SEED)
  # One ready-to-run fiber arg list per particle (OBV is OFF here; the worker
  # switches it on per arm below).
  args_by_scenario[[sc]] <- lapply(seq_len(N_SETS), function(i)
    build_abc_model_args(R0 = theta$R0[i], prop_funeral = theta$prop_funeral[i],
                         hcw_risk_scalar = theta$hcw_risk_scalar[i],
                         base = mp$base_args, tv = mp$tv_args,
                         D = sol$D_direct_multiplier, F_fun = sol$F_funeral_multiplier,
                         seeding_cases = SEEDING_CASES, hcw_base_prob = HCW_BASE_PROB))
}


# ---------------------------------------------------------------------------
# 4. BUILD THE JOB LIST (scenario x particle x rep x arm) AND RUN IN PARALLEL
# ---------------------------------------------------------------------------
# Arms = "no_obv" (the shared baseline) + one per efficacy. simulate_one() turns
# OBV on only when arm == "obv", so we run each arm via simulate_one() with the
# matching config and tag the descriptive arm name onto the result.
arms_cfg <- c(list(no_obv = NULL),
              lapply(OBV_EFFICACIES, function(e) c(OBV_BASE, list(efficacy = e))))

jobs <- list(); k <- 0L
for (sc in names(SCENARIOS)) {
  sc_i <- match(sc, names(SCENARIOS))
  for (s in seq_len(N_SETS)) for (r in seq_len(N_REPS)) {
    # Seed depends on (scenario, set, rep) only -- arms within a rep share it.
    seed_sr <- SEED_BASE + ((sc_i - 1L) * N_SETS + (s - 1L)) * N_REPS + (r - 1L)
    for (arm in names(arms_cfg)) {
      k <- k + 1L
      jobs[[k]] <- list(scenario = sc, set_id = s, rep_id = r, arm = arm, seed = seed_sr)
    }
  }
}
message(sprintf("Total simulations: %d (%d scenarios x %d sets x %d reps x %d arms).",
                length(jobs), length(SCENARIOS), N_SETS, N_REPS, length(arms_cfg)))

plan(multisession, workers = min(N_CORES, future::availableCores()))
on.exit(plan(sequential), add = TRUE)

with_progress({
  p <- progressor(along = jobs)
  results <- future_lapply(jobs, function(job, args_by_scenario, arms_cfg) {
    # Collapse the descriptive arm (no_obv / obv_40 / obv_80) to simulate_one()'s
    # "no_obv"/"obv" switch + matching config, then restore the descriptive label.
    cfg <- arms_cfg[[job$arm]]
    j   <- list(set_id = job$set_id, rep_id = job$rep_id, seed = job$seed,
                arm = if (is.null(cfg)) "no_obv" else "obv")
    out <- simulate_one(j, args_list = args_by_scenario[[job$scenario]],
                        obv_cfg = if (is.null(cfg)) list() else cfg)
    out$scenario <- job$scenario
    out$arm      <- job$arm
    p(); out
  }, args_by_scenario = args_by_scenario, arms_cfg = arms_cfg,
     future.packages = "fiber", future.seed = TRUE)
})
plan(sequential)


# ---------------------------------------------------------------------------
# 5. % HCW DEATHS AVERTED -- DIFFERENCE OF ARM MEANS, PER POSTERIOR DRAW
# ---------------------------------------------------------------------------
per_rep <- data.frame(
  scenario     = vapply(results, `[[`, character(1), "scenario"),
  arm          = vapply(results, `[[`, character(1), "arm"),
  set_id       = vapply(results, `[[`, integer(1),   "set_id"),
  n_hcw_deaths = vapply(results, `[[`, integer(1),   "n_hcw_deaths"),
  stringsAsFactors = FALSE
)

# Mean HCW deaths over reps, per (scenario, set, arm) -- the per-arm EXPECTED
# burden for that posterior draw (unbiased despite OBV's RNG desync; see header).
mu <- aggregate(n_hcw_deaths ~ scenario + set_id + arm, data = per_rep, FUN = mean)

# Pull out the shared no-OBV baseline and merge onto each OBV arm.
base_mu <- mu[mu$arm == "no_obv", c("scenario", "set_id", "n_hcw_deaths")]
names(base_mu)[3] <- "mu_no"
averted <- merge(mu[mu$arm != "no_obv", ], base_mu, by = c("scenario", "set_id"))
names(averted)[names(averted) == "n_hcw_deaths"] <- "mu_obv"

# Per-draw TOTAL effect: averted burden and % averted (>= 0 in expectation; tiny
# negatives can still occur from finite-rep noise where the baseline is small).
averted$hcw_averted     <- averted$mu_no - averted$mu_obv
averted$pct_hcw_averted <- ifelse(averted$mu_no > 0,
                                  100 * averted$hcw_averted / averted$mu_no, NA_real_)

# Friendly labels.
averted$scenario_lbl <- c(DRC = "DRC", WestAfrica = "West Africa")[averted$scenario]
averted$efficacy_lbl <- sprintf("%d%% efficacy", round(100 * OBV_EFFICACIES[averted$arm]))

# Headline: per-draw median [IQR] % averted, plus the burden-weighted pooled %
# (sum of averted over sum of baseline -- weights draws by outbreak size).
summ <- do.call(rbind, lapply(split(averted, list(averted$scenario_lbl, averted$efficacy_lbl), drop = TRUE),
  function(d) {
    q <- quantile(d$pct_hcw_averted, c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE)
    pooled <- 100 * sum(d$hcw_averted) / sum(d$mu_no)
    data.frame(scenario = d$scenario_lbl[1], efficacy = d$efficacy_lbl[1],
               median = round(q[2], 1), iqr_lo = round(q[1], 1), iqr_hi = round(q[3], 1),
               pooled = round(pooled, 1), stringsAsFactors = FALSE)
  }))
print(summ, row.names = FALSE)


# ---------------------------------------------------------------------------
# 6. PLOT: distribution of % HCW deaths averted, by scenario and efficacy
# ---------------------------------------------------------------------------
ggplot(averted[!is.na(averted$pct_hcw_averted), ],
       aes(scenario_lbl, pct_hcw_averted, fill = efficacy_lbl)) +
  geom_boxplot(outlier.size = 0.5, position = position_dodge(0.8), width = 0.7) +
  labs(x = NULL, y = "HCW deaths averted (%)", fill = "Obeldesivir",
       title = "Obeldesivir impact on HCW deaths, by scenario and efficacy",
       subtitle = sprintf("TOTAL effect (with vs without OBV, difference of arm means); %d posterior draws x %d reps per arm",
                          N_SETS, N_REPS)) +
  theme_bw()
