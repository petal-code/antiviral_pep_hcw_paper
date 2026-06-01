# =============================================================================
# run_obeldesivir_two_scenario_comparison_directVsTotal.R
#
# Plots the DIRECT vs TOTAL obeldesivir (OBV) effect side by side, for BOTH
# scenarios (DRC, West Africa) and TWO efficacies (40%, 80%). The gap between
# them is exactly the benefit OBV gets from averting ONWARD transmission chains
# (not just the index HCW infections it blocks).
#
# Two estimands, BOTH extracted from ONE set of simulations:
#   DIRECT (within-run counterfactual; conservative): of the HCW infections OBV
#     blocked, how many would themselves have died -- excludes their onward
#     chains. Per OBV run, fiber reports these as
#       out$sim_info$obv_pep_num_treated$prevented_deaths
#     %direct = 100 * prevented / (realised + prevented).
#   TOTAL (with vs without OBV, difference of arm MEANS; full effect): includes
#     averted onward chains. Needs a no-OBV baseline arm.
#     %total = 100 * (mean_no_obv - mean_obv) / mean_no_obv.
#   (Means are unbiased despite OBV's RNG desync -- expectations don't depend on
#   coupling -- so the difference is valid; see the totalEffect script header.)
#
# Each OBV run yields BOTH its realised HCW deaths (-> total) and its prevented
# deaths (-> direct), so no extra simulation is needed for the comparison; we
# only add the shared no-OBV arm the total effect requires.
#
# Defines NO new functions -- composes existing functions/ helpers with base R +
# ggplot. Run from the repo root:
#   Rscript inst/run_obeldesivir_two_scenario_comparison_directVsTotal.R
# =============================================================================


# ---------------------------------------------------------------------------
# 0. CONFIGURATION (everything tweakable lives here)
# ---------------------------------------------------------------------------
# The TOTAL effect differences arm means of a bimodal outbreak-size distribution,
# so it needs many reps to settle. 50 is a sensible floor.
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
# Arms = "no_obv" (shared baseline for TOTAL) + one per efficacy. Every run is
# executed directly (not via simulate_one) so we can read both realised HCW
# deaths AND the deferred prevented-deaths counterfactual from the same call.
arms <- c("no_obv", names(OBV_EFFICACIES))

jobs <- list(); k <- 0L
for (sc in names(SCENARIOS)) {
  sc_i <- match(sc, names(SCENARIOS))
  for (s in seq_len(N_SETS)) for (r in seq_len(N_REPS)) {
    # Seed depends on (scenario, set, rep) only -- arms within a rep share it.
    seed_sr <- SEED_BASE + ((sc_i - 1L) * N_SETS + (s - 1L)) * N_REPS + (r - 1L)
    for (arm in arms) {
      k <- k + 1L
      jobs[[k]] <- list(scenario = sc, set_id = s, rep_id = r, arm = arm,
                        eff_value = if (arm == "no_obv") NA_real_ else OBV_EFFICACIES[[arm]],
                        seed = seed_sr)
    }
  }
}
message(sprintf("Total simulations: %d (%d scenarios x %d sets x %d reps x %d arms).",
                length(jobs), length(SCENARIOS), N_SETS, N_REPS, length(arms)))

plan(multisession, workers = min(N_CORES, future::availableCores()))
on.exit(plan(sequential), add = TRUE)

with_progress({
  p <- progressor(along = jobs)
  results <- future_lapply(jobs, function(job, args_by_scenario, OBV_BASE) {
    a <- args_by_scenario[[job$scenario]][[job$set_id]]
    a$seed <- job$seed
    if (job$arm != "no_obv") {                 # OBV arms: turn the gate on
      a$obv_pep_enabled          <- TRUE
      a$obv_pep_coverage         <- OBV_BASE$coverage
      a$obv_pep_adherence        <- OBV_BASE$adherence
      a$obv_pep_efficacy         <- job$eff_value
      a$obv_pep_dpc              <- OBV_BASE$dpc
      a$obv_pep_target_class     <- OBV_BASE$target_class
      a$obv_pep_target_locations <- OBV_BASE$target_locations
    }
    out <- do.call(fiber::branching_process_main, a)

    tdf <- out$tdf
    tdf <- tdf[!is.na(tdf$time_infection_absolute), , drop = FALSE]
    realised_hcw <- sum(!is.na(tdf$outcome) & tdf$outcome &
                        !is.na(tdf$class) & tdf$class == "HCW")
    prevented_hcw <- out$sim_info$obv_pep_num_treated$prevented_deaths   # 0 for no_obv
    if (is.null(prevented_hcw) || length(prevented_hcw) == 0L || is.na(prevented_hcw)) prevented_hcw <- 0

    p()
    list(scenario = job$scenario, arm = job$arm, set_id = job$set_id,
         realised_hcw = realised_hcw, prevented_hcw = prevented_hcw)
  }, args_by_scenario = args_by_scenario, OBV_BASE = OBV_BASE,
     future.packages = "fiber", future.seed = TRUE)
})
plan(sequential)


# ---------------------------------------------------------------------------
# 5. COMPUTE BOTH ESTIMANDS PER PARTICLE
# ---------------------------------------------------------------------------
per_run <- data.frame(
  scenario      = vapply(results, `[[`, character(1), "scenario"),
  arm           = vapply(results, `[[`, character(1), "arm"),
  set_id        = vapply(results, `[[`, integer(1),   "set_id"),
  realised_hcw  = vapply(results, `[[`, numeric(1),   "realised_hcw"),
  prevented_hcw = vapply(results, `[[`, numeric(1),   "prevented_hcw"),
  stringsAsFactors = FALSE
)
if (sum(per_run$prevented_hcw) == 0)
  warning("No prevented HCW deaths recorded -- check fiber is the updated version (PR #70).",
          call. = FALSE)

# Mean realised HCW deaths over reps, per (scenario, set, arm).
mu <- aggregate(realised_hcw ~ scenario + set_id + arm, data = per_run, FUN = mean)
base_mu <- mu[mu$arm == "no_obv", c("scenario", "set_id", "realised_hcw")]
names(base_mu)[3] <- "mu_no"

# --- TOTAL: (mean_no_obv - mean_obv) / mean_no_obv, per particle x efficacy ----
total <- merge(mu[mu$arm != "no_obv", ], base_mu, by = c("scenario", "set_id"))
names(total)[names(total) == "realised_hcw"] <- "mu_obv"
total$pct <- ifelse(total$mu_no > 0, 100 * (total$mu_no - total$mu_obv) / total$mu_no, NA_real_)
total$estimand <- "Total (incl. averted chains)"

# --- DIRECT: prevented / (realised + prevented), per particle x efficacy -------
# Sum across reps within the OBV arm (burden-weighted, stable denominator).
obv_runs <- per_run[per_run$arm != "no_obv", ]
obv_runs$counterfactual_hcw <- obv_runs$realised_hcw + obv_runs$prevented_hcw
direct_agg <- aggregate(cbind(prevented_hcw, counterfactual_hcw) ~ scenario + set_id + arm,
                        data = obv_runs, FUN = sum)
direct_agg$pct <- ifelse(direct_agg$counterfactual_hcw > 0,
                         100 * direct_agg$prevented_hcw / direct_agg$counterfactual_hcw, NA_real_)
direct_agg$estimand <- "Direct (index deaths only)"

# Stack the two estimands into one long frame for a faceted comparison.
combined <- rbind(
  total[,      c("scenario", "set_id", "arm", "pct", "estimand")],
  direct_agg[, c("scenario", "set_id", "arm", "pct", "estimand")]
)
combined$scenario_lbl <- c(DRC = "DRC", WestAfrica = "West Africa")[combined$scenario]
combined$efficacy_lbl <- sprintf("%d%% efficacy", round(100 * OBV_EFFICACIES[combined$arm]))
# Order so Direct is drawn left of Total.
combined$estimand <- factor(combined$estimand,
                            levels = c("Direct (index deaths only)", "Total (incl. averted chains)"))

# Headline table: median % averted per scenario x efficacy x estimand.
summ <- aggregate(pct ~ scenario_lbl + efficacy_lbl + estimand, data = combined,
                  FUN = function(x) round(median(x, na.rm = TRUE), 1))
print(summ[order(summ$scenario_lbl, summ$efficacy_lbl, summ$estimand), ], row.names = FALSE)


# ---------------------------------------------------------------------------
# 6. PLOT: Direct vs Total side by side, by scenario and efficacy
# ---------------------------------------------------------------------------
# Two panels (Direct | Total); within each, scenario on x, efficacy by fill. The
# vertical step from the Direct panel to the Total panel is the onward-chain gain.
ggplot(combined[!is.na(combined$pct), ],
       aes(scenario_lbl, pct, fill = efficacy_lbl)) +
  geom_boxplot(outlier.size = 0.5, position = position_dodge(0.8), width = 0.7) +
  facet_wrap(~ estimand) +
  labs(x = NULL, y = "HCW deaths averted (%)", fill = "Obeldesivir",
       title = "Obeldesivir impact on HCW deaths: direct vs total effect",
       subtitle = sprintf("Direct = index deaths only; Total = incl. averted onward chains. %d posterior draws x %d reps per arm",
                          N_SETS, N_REPS)) +
  theme_bw()
