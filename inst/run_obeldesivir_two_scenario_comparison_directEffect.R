# =============================================================================
# run_obeldesivir_two_scenario_comparison_directEffect.R
#
# Estimates the % of healthcare-worker (HCW) deaths that obeldesivir (OBV)
# post-exposure prophylaxis averts, for BOTH fitted scenarios (DRC and West
# Africa) and at TWO OBV efficacies (40% and 80%), then draws one stratified plot.
#
# METHOD (requires the updated fiber, petal-code/fiber PR #70):
#   Each OBV run carries its OWN internal counterfactual. After the outbreak is
#   simulated, fiber replays the infections the OBV gate prevented through the
#   same outcome model and reports their would-be deaths as
#       out$sim_info$obv_pep_num_treated$prevented_deaths
#   This is computed AFTER the trajectory (drawing no RNG inside it), so it does
#   not perturb the simulation. Therefore a SINGLE OBV run gives both:
#       realised HCW deaths (in tdf, with OBV in place), and
#       HCW deaths OBV prevented (prevented_deaths),
#   from which, per run:
#       counterfactual HCW deaths (no OBV) = realised + prevented
#       HCW deaths averted                 = prevented
#       % HCW deaths averted = 100 * prevented / (realised + prevented)   in [0,100]
#   No separate no-OBV arm is needed, the two arms can't desync, and the % can
#   never be negative -- which fixes the artefacts the paired approach produced.
#
#   CAVEAT: prevented_deaths is the DIRECT would-be death of each prevented index
#   infection ONLY; it excludes the onward transmission chains those infections
#   would have seeded. So this is a CONSERVATIVE (direct) estimate of % averted.
#
# Defines NO new functions -- composes existing functions/ helpers with base R +
# ggplot. Run from the repo root:
#   Rscript inst/run_obeldesivir_two_scenario_comparison.R
# =============================================================================


# ---------------------------------------------------------------------------
# 0. CONFIGURATION (everything tweakable lives here)
# ---------------------------------------------------------------------------
N_REPS  <- 5L     # stochastic replicates per particle per efficacy
N_CORES <- 14L    # parallel workers (capped at availableCores() below)

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
  # switches it on per efficacy below).
  args_by_scenario[[sc]] <- lapply(seq_len(N_SETS), function(i)
    build_abc_model_args(R0 = theta$R0[i], prop_funeral = theta$prop_funeral[i],
                         hcw_risk_scalar = theta$hcw_risk_scalar[i],
                         base = mp$base_args, tv = mp$tv_args,
                         D = sol$D_direct_multiplier, F_fun = sol$F_funeral_multiplier,
                         seeding_cases = SEEDING_CASES, hcw_base_prob = HCW_BASE_PROB))
}


# ---------------------------------------------------------------------------
# 4. BUILD THE JOB LIST (scenario x particle x rep x efficacy) AND RUN IN PARALLEL
# ---------------------------------------------------------------------------
# One job = one OBV simulation. No no-OBV arm: each run is its own counterfactual.
jobs <- list(); k <- 0L
for (sc in names(SCENARIOS)) {
  sc_i <- match(sc, names(SCENARIOS))
  for (s in seq_len(N_SETS)) for (r in seq_len(N_REPS)) {
    seed_sr <- SEED_BASE + ((sc_i - 1L) * N_SETS + (s - 1L)) * N_REPS + (r - 1L)
    for (e in names(OBV_EFFICACIES)) {
      k <- k + 1L
      jobs[[k]] <- list(scenario = sc, set_id = s, efficacy = e,
                        eff_value = OBV_EFFICACIES[[e]], seed = seed_sr)
    }
  }
}
message(sprintf("Total simulations: %d (%d scenarios x %d sets x %d reps x %d efficacies).",
                length(jobs), length(SCENARIOS), N_SETS, N_REPS, length(OBV_EFFICACIES)))

plan(multisession, workers = min(N_CORES, future::availableCores()))
on.exit(plan(sequential), add = TRUE)

with_progress({
  p <- progressor(along = jobs)
  results <- future_lapply(jobs, function(job, args_by_scenario, OBV_BASE) {
    # Switch OBV on for this efficacy, then run the model directly so we can read
    # the deferred prevented-deaths counterfactual (simulate_one() doesn't expose it).
    a <- args_by_scenario[[job$scenario]][[job$set_id]]
    a$seed                     <- job$seed
    a$obv_pep_enabled          <- TRUE
    a$obv_pep_coverage         <- OBV_BASE$coverage
    a$obv_pep_adherence        <- OBV_BASE$adherence
    a$obv_pep_efficacy         <- job$eff_value
    a$obv_pep_dpc              <- OBV_BASE$dpc
    a$obv_pep_target_class     <- OBV_BASE$target_class
    a$obv_pep_target_locations <- OBV_BASE$target_locations

    out <- do.call(fiber::branching_process_main, a)

    # Realised HCW deaths in this (with-OBV) outbreak.
    tdf <- out$tdf
    tdf <- tdf[!is.na(tdf$time_infection_absolute), , drop = FALSE]
    realised_hcw <- sum(!is.na(tdf$outcome) & tdf$outcome &
                        !is.na(tdf$class) & tdf$class == "HCW")
    # Direct HCW deaths OBV prevented (deferred counterfactual; >= 0 by construction).
    prevented_hcw <- out$sim_info$obv_pep_num_treated$prevented_deaths
    if (is.null(prevented_hcw) || length(prevented_hcw) == 0L || is.na(prevented_hcw)) prevented_hcw <- 0

    p()
    list(scenario = job$scenario, efficacy = job$efficacy, set_id = job$set_id,
         realised_hcw = realised_hcw, prevented_hcw = prevented_hcw)
  }, args_by_scenario = args_by_scenario, OBV_BASE = OBV_BASE,
     future.packages = "fiber", future.seed = TRUE)
})
plan(sequential)


# ---------------------------------------------------------------------------
# 5. % HCW DEATHS AVERTED (per-run counterfactual; pooled per particle)
# ---------------------------------------------------------------------------
per_run <- data.frame(
  scenario      = vapply(results, `[[`, character(1), "scenario"),
  efficacy      = vapply(results, `[[`, character(1), "efficacy"),
  set_id        = vapply(results, `[[`, integer(1),   "set_id"),
  realised_hcw  = vapply(results, `[[`, numeric(1),   "realised_hcw"),
  prevented_hcw = vapply(results, `[[`, numeric(1),   "prevented_hcw"),
  stringsAsFactors = FALSE
)
# Counterfactual (no-OBV) HCW deaths = realised + prevented, per run.
per_run$counterfactual_hcw <- per_run$realised_hcw + per_run$prevented_hcw

if (sum(per_run$prevented_hcw) == 0)
  warning("No prevented HCW deaths recorded in any run -- check fiber is the updated ",
          "version (PR #70) and OBV coverage > 0.", call. = FALSE)

# Per particle: SUM prevented and counterfactual across its reps (burden-weighted,
# so the denominator is stable), then % averted. The distribution of these
# per-particle values across the posterior is the parameter-variation story --
# and, with this method, every value is well-defined and in [0, 100].
agg <- aggregate(cbind(prevented_hcw, counterfactual_hcw) ~ scenario + efficacy + set_id,
                 data = per_run, FUN = sum)
agg$pct_averted <- ifelse(agg$counterfactual_hcw > 0,
                          100 * agg$prevented_hcw / agg$counterfactual_hcw, NA_real_)

# Friendly labels.
agg$scenario_lbl <- c(DRC = "DRC", WestAfrica = "West Africa")[agg$scenario]
agg$efficacy_lbl <- sprintf("%d%% efficacy", round(100 * OBV_EFFICACIES[agg$efficacy]))

# Headline summaries: per-particle median [IQR], and the pooled burden % averted.
summ <- do.call(rbind, lapply(split(agg, list(agg$scenario_lbl, agg$efficacy_lbl), drop = TRUE),
  function(d) {
    q <- quantile(d$pct_averted, c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE)
    pooled <- 100 * sum(d$prevented_hcw) / sum(d$counterfactual_hcw)  # burden-weighted point
    data.frame(scenario = d$scenario_lbl[1], efficacy = d$efficacy_lbl[1],
               median = round(q[2], 1), iqr_lo = round(q[1], 1), iqr_hi = round(q[3], 1),
               pooled = round(pooled, 1), stringsAsFactors = FALSE)
  }))
print(summ, row.names = FALSE)


# ---------------------------------------------------------------------------
# 6. PLOT: distribution of % HCW deaths averted, by scenario and efficacy
# ---------------------------------------------------------------------------
ggplot(agg[!is.na(agg$pct_averted), ],
       aes(scenario_lbl, pct_averted, fill = efficacy_lbl)) +
  geom_boxplot(outlier.size = 0.5, position = position_dodge(0.8), width = 0.7) +
  labs(x = NULL, y = "HCW deaths averted (%)", fill = "Obeldesivir",
       title = "Obeldesivir impact on HCW deaths, by scenario and efficacy",
       subtitle = sprintf("Per-run counterfactual (direct prevented deaths); %d posterior draws x %d reps per efficacy",
                          N_SETS, N_REPS)) +
  theme_bw()
