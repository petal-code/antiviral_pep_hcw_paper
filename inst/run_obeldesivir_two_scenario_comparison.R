# =============================================================================
# run_obeldesivir_two_scenario_comparison.R
#
# Estimates the % of healthcare-worker (HCW) deaths that obeldesivir (OBV)
# post-exposure prophylaxis would avert, for BOTH fitted scenarios (DRC and
# West Africa) and at TWO OBV efficacies (40% and 80%), then draws a single
# stratified comparison plot.
#
# Workflow (reusing the project's existing helpers; NO new functions defined):
#   1. Read the two fitted ABC posteriors shipped in inst/.
#   2. Downsample BOTH to the same size -- the smaller of the two particle
#      counts -- so the scenarios are compared on an equal footing.
#   3. For each scenario, convert every posterior draw (R0, prop_funeral,
#      hcw_risk_scalar) into a fiber argument list, using that scenario's own
#      time-varying inputs and the calibration's exact mapping.
#   4. Run every (scenario x particle x replicate) under 3 arms -- no OBV,
#      OBV@40%, OBV@80% -- in parallel, pairing arms on a shared seed.
#   5. Compare ARM-LEVEL DISTRIBUTIONS (not paired reps -- OBV decouples the
#      per-arm RNG stream, so per-rep differences are meaningless): % HCW deaths
#      averted = 100*(median_no_obv - median_obv)/median_no_obv per scenario,
#      with a bootstrap 95% CI; plot by scenario x efficacy.
#
# Run from the repo root with:  Rscript inst/run_obeldesivir_two_scenario_comparison.R
# =============================================================================


# ---------------------------------------------------------------------------
# 0. CONFIGURATION (everything tweakable lives here)
# ---------------------------------------------------------------------------
N_REPS  <- 5L     # stochastic replicates per particle per arm (1 = one draw per particle)
N_CORES <- 14L    # parallel workers (capped at availableCores() below)

<<<<<<< Updated upstream
# A parameter set enters the % averted plot only if its no-OBV baseline outbreak
# TOOK OFF, i.e. produced at least this many TOTAL deaths. The % is
# 100*(baseline - obv)/baseline, so a small/fizzled outbreak makes it explode
# (and, because OBV decouples the per-arm RNG stream, occasionally go sharply
# negative). Gating on total deaths is a clean "did the outbreak happen" filter;
# the ABSOLUTE averted (computed alongside) is robust and never blows up.
MIN_BASELINE_DEATHS <- 50

=======
>>>>>>> Stashed changes
# OBV efficacies to test (the "with obeldesivir" arms). Add/remove freely.
OBV_EFFICACIES <- c(obv_40 = 0.40, 
                    obv_80 = 0.80)
# Fixed OBV delivery settings shared by every OBV arm: full coverage/adherence,
# modelled as PEP for HCWs exposed in hospital (matches the WHO CORC analysis).
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
HCW_BASE_PROB <- 0.25 ## note this relates a fitted parameter, which scales this to influence amount of healthcare transmission
                      ## this will be replaced with a different fitted parameter in the final version of the fits
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

# Per-arm OBV config: no_obv = NULL (gate stays off); each OBV arm = base + efficacy.
## Creating a list of model parameters for each of the two different obeldesivir parameterisations
arms_cfg <- c(list(no_obv = NULL),
              lapply(OBV_EFFICACIES, function(e) c(OBV_BASE, list(efficacy = e))))


# ---------------------------------------------------------------------------
# 2. READ BOTH POSTERIORS, THEN PICK THE COMMON (SMALLER) SAMPLE SIZE
# ---------------------------------------------------------------------------
# The inst/ files are EasyABC result objects: $param (3 columns: R0,
# prop_funeral, hcw_risk_scalar) and $weights. Rebuild the weighted-posterior
# data frame downsample_posterior() expects.
posteriors <- lapply(SCENARIOS, function(sc) {
  res <- readRDS(here::here("inst", sc$rds))  ## reading in the fitted results for each scenario
  data.frame(weight = res$weights, 
             R0 = res$param[, 1],               ## first fitted parameter is R0 - this will stay the same in the new fitting version currently underway
             prop_funeral = res$param[, 2],     ## second fitted parameter is % of transmission attributable to funerals - this will stay the same in the new fitting version currently underway
             hcw_risk_scalar = res$param[, 3])  ## third fitted parameter is a parameter that scales in-hospital transmission towards HCW - this will change in the new fitting version
})
N_SETS <- min(50, ## change this number if you want to do fewer that the number of posterior draws available
              min(vapply(posteriors, nrow, integer(1))))   # calculate which output has the least number of posterior draws
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
  # Scenario-specific parameters + time-varying curves (prop_etu, ppe_coverage etc).
  mp <- make_model_parameters(SCENARIOS[[sc]]$id, scenario_matrix,
                              overrides = list(check_final_size = CHECK_FINAL_SIZE))
  # D and F are the t=0 direct/funeral R0 multipliers for THIS scenario (they
  # depend on its prop_etu(0) and efficacies), so they are solved per scenario.
  sol <- solve_offspring_means_for_R0(R0 = 1.0, args = mp$args,
                                      proportion_transmission_from_funerals = SETUP_FUNERAL_SHARE,
                                      n = SETUP_R0_N, seed = SETUP_R0_SEED)
  theta <- downsample_posterior(posteriors[[sc]], n_sets = N_SETS, seed = RESAMPLE_SEED)
  # One ready-to-run fiber arg list per particle (default scalars + scenario tv
  # + this particle's R0 / prop_funeral / hcw_risk_scalar via the ABC mapping).
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
# One job = one fiber simulation. The seed depends on (scenario, set, rep) but
# NOT the arm, so the 3 arms share a stochastic history -> a paired comparison.
jobs <- list()
k <- 0L
for (sc in names(SCENARIOS)) {
  sc_i <- match(sc, names(SCENARIOS))
  for (s in seq_len(N_SETS)) for (r in seq_len(N_REPS)) {
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
    # simulate_one() switches OBV on only when arm == "obv", reading efficacy
    # from obv_cfg. We have several OBV arms, so collapse the descriptive label
    # to "obv"/"no_obv" for the call, pass the matching config, then restore the
    # descriptive label on the result.
    cfg <- arms_cfg[[job$arm]] ## get the specific model parameters for this scenario and posterior draw 
    j   <- list(set_id = job$set_id, rep_id = job$rep_id, seed = job$seed, arm = if (is.null(cfg)) "no_obv" else "obv")
    out <- simulate_one(j, args_list = args_by_scenario[[job$scenario]],
                        obv_cfg = if (is.null(cfg)) list() else cfg)
    out$scenario <- job$scenario
    out$arm      <- job$arm          # overwrite "obv" with "obv_40" / "obv_80"
    p(); out
  }, args_by_scenario = args_by_scenario, arms_cfg = arms_cfg,
     future.packages = "fiber", future.seed = TRUE)
})
plan(sequential)


# ---------------------------------------------------------------------------
# 5. % HCW DEATHS AVERTED -- DISTRIBUTION-LEVEL (unpaired) ESTIMAND
# ---------------------------------------------------------------------------
# WHY NOT PAIR PER (set, rep): enabling OBV consumes extra rbinom() draws inside
# branching_process_main(), so an OBV arm and its same-seed no-OBV arm DESYNC as
# soon as the first OBV draw fires (documented fiber caveat). Near R0 = 1 with
# heavy overdispersion, that desync routinely flips a rep between fizzle and full
# take-off -- so a per-rep difference compares two *independent* outbreaks and its
# sign is meaningless (the diagnostic showed total deaths rose under OBV in ~48%
# of reps, which is impossible in a truly paired run). fiber's own guidance:
# "use many replicates and compare distributions." So we do exactly that.
per_rep <- data.frame(
  scenario     = vapply(results, `[[`, character(1), "scenario"),
  arm          = vapply(results, `[[`, character(1), "arm"),
  n_deaths     = vapply(results, `[[`, integer(1),   "n_deaths"),
  n_hcw_deaths = vapply(results, `[[`, integer(1),   "n_hcw_deaths"),
  stringsAsFactors = FALSE
)

# Compare ARM-LEVEL DISTRIBUTIONS, not paired reps. Within each (scenario, arm)
# we summarise all N_SETS x N_REPS draws by their MEDIAN HCW deaths, then
#     % averted = 100 * (median_no_obv - median_arm) / median_no_obv.
# A CI comes from a non-parametric bootstrap: resample the reps within each arm,
# recompute the two medians and the % -- no pairing, so it cannot blow up.
N_BOOT <- 2000L

pct_averted_ci <- do.call(rbind, lapply(names(SCENARIOS), function(sc) {
  d   <- per_rep[per_rep$scenario == sc, ]
  obs <- tapply(d$n_hcw_deaths, d$arm, median)   # median HCW deaths per arm
  do.call(rbind, lapply(names(OBV_EFFICACIES), function(e) {
    # Point estimate from the observed medians.
    point <- 100 * (obs[["no_obv"]] - obs[[e]]) / obs[["no_obv"]]
    # Bootstrap: resample reps within the no_obv and this-efficacy arm separately.
    base_v <- d$n_hcw_deaths[d$arm == "no_obv"]
    obv_v  <- d$n_hcw_deaths[d$arm == e]
    boot <- replicate(N_BOOT, {
      mb <- median(sample(base_v, length(base_v), replace = TRUE))
      mo <- median(sample(obv_v,  length(obv_v),  replace = TRUE))
      100 * (mb - mo) / mb
    })
    data.frame(scenario = sc, arm = e,
               pct_averted = point,
               lo = quantile(boot, 0.025, names = FALSE),
               hi = quantile(boot, 0.975, names = FALSE),
               stringsAsFactors = FALSE)
  }))
}))

# Friendly labels for the plot.
pct_averted_ci$scenario_lbl <- c(DRC = "DRC", WestAfrica = "West Africa")[pct_averted_ci$scenario]
pct_averted_ci$efficacy_lbl <- sprintf("%d%% efficacy", round(100 * OBV_EFFICACIES[pct_averted_ci$arm]))
print(pct_averted_ci[, c("scenario_lbl", "efficacy_lbl", "pct_averted", "lo", "hi")], row.names = FALSE)

# ---------------------------------------------------------------------------
# 6. PLOT: % HCW deaths averted (distribution-level), by scenario and efficacy
# ---------------------------------------------------------------------------
# A point estimate per scenario x efficacy with its bootstrap 95% interval.
ggplot(pct_averted_ci, aes(scenario_lbl, pct_averted, colour = efficacy_lbl)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60") +
  geom_pointrange(aes(ymin = lo, ymax = hi),
                  position = position_dodge(0.4), linewidth = 0.8) +
  labs(x = NULL, y = "HCW deaths averted (%)", colour = "Obeldesivir",
       title = "Obeldesivir impact on HCW deaths, by scenario and efficacy",
       subtitle = sprintf("Distribution-level (unpaired) estimate across %d posterior draws x %d reps per arm; bars = bootstrap 95%% CI",
                          N_SETS, N_REPS)) +
  theme_bw()
