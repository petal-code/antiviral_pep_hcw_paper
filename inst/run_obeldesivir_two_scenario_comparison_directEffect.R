# =============================================================================
# run_obeldesivir_two_scenario_comparison_directEffect.R
#
# Estimates the % of healthcare-worker (HCW) deaths that obeldesivir (OBV)
# post-exposure prophylaxis averts, for each fitted scenario and at TWO OBV
# efficacies (40% and 80%), then draws one stratified plot.
#
# >>> NPI-EFFICACY FITS <<<  This script now consumes the NEW fitting scheme that
# fits (R0, prop_funeral, npi_scaler) instead of (R0, prop_funeral,
# hcw_risk_scalar). The single npi_scaler in [-1, 1] positions the two fitted
# conditional efficacies (ppe_efficacy, etu_efficacy) on their [min, max]
# intervals (the "npi_spec"); see functions/abc_calibration_functions_npi.R and
# npi_efficacy_from_scaler(). Consequences for this script:
#   * the posterior's 3rd column is npi_scaler, not hcw_risk_scalar;
#   * because etu_efficacy is FITTED, the direct multiplier D depends on the
#     particle, so we cache the efficacy-independent R0 invariants ONCE per
#     scenario (compute_R0_invariants()) and let the NPI build_abc_model_args()
#     recompute D/F per particle (it also sets ppe/etu efficacy on the args);
#   * mapping npi_scaler -> efficacies REQUIRES the npi_spec the fit used. We read
#     it from the fit itself when present (newer fits embed it / ship a sidecar),
#     else fall back to a per-scenario constant below (see SCENARIOS).
# Currently only the DRC (PlusPlus) NPI fit exists, so the "comparison" runs on
# that one scenario; West Africa is a ready-to-uncomment placeholder for when its
# NPI fit lands. The machinery handles any number of scenarios.
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
#   Rscript inst/run_obeldesivir_two_scenario_comparison_directEffect.R
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

# The scenarios: a friendly label -> (NPI-efficacy posterior RDS [repo-relative],
# scenario id, scenario CSV, plot label, and a npi_spec fallback used ONLY if the
# fit does not carry its own npi_spec). The fallback below is the npi_spec the
# 154658 DRC PlusPlus fit was run with (confirmed; see the fit's _metadata.txt).
SCENARIOS <- list(
  DRC_PlusPlus = list(
    rds   = file.path("outputs", "02_ABC_model_fits_NPI_Eff",
                      "fiber_ABC_SMC_Middle_DRC_ConflictSmoothed_PlusPlus_NPIeff_20260601_154658.rds"),
    id    = "Middle_DRC_ConflictSmoothed_PlusPlus",
    csv   = "final_six_scenario_values_original_approach.csv",
    label = "DRC (PlusPlus)",
    npi_spec_fallback = list(ppe_efficacy = list(min = 0.30, max = 0.90),
                             etu_efficacy = list(min = 0.60, max = 0.95))
  )
  # ,
  # To add West Africa once its NPI fit exists, point `rds` at it and set the
  # matching `csv`/`npi_spec_fallback`:
  # WestAfrica = list(
  #   rds   = file.path("outputs", "02_ABC_model_fits_NPI_Eff", "<WA NPI fit>.rds"),
  #   id    = "Worst_WestAfrica",
  #   csv   = "final_four_scenario_values.csv",
  #   label = "West Africa",
  #   npi_spec_fallback = list(ppe_efficacy = list(min = 0.20, max = 0.90),
  #                            etu_efficacy = list(min = 0.50, max = 0.95)))
)

# Calibration settings reproduced so the posterior -> model-args conversion is
# identical to how each NPI ABC was run (see analyses/02_ABC_model_fits_NPI_Eff/).
# The fixed (not fitted) efficacies general_hospital_quarantine_efficacy and
# safe_funeral_efficacy come from DEFAULT_SCALAR_INPUTS, exactly as the fit did.
SEEDING_CASES <- 25L
CHECK_FINAL_SIZE <- 10000L
SETUP_R0_N <- 100000     # MC draws for compute_R0_invariants() (match calibration)
SETUP_R0_SEED <- 42L
RESAMPLE_SEED <- 1L; SEED_BASE <- 20260601L

# Where to save results (plot + tables + the npi_spec actually used).
OUTPUT_DIR <- file.path("outputs", "03_obeldesivir_impact_NPI")


# ---------------------------------------------------------------------------
# 1. LIBRARIES + SHARED HELPERS
# ---------------------------------------------------------------------------
library(fiber); library(future); library(future.apply); library(progressr); library(ggplot2)
FN <- here::here("functions")
source(file.path(FN, "setup_model_parameters.R"))            # scenario -> model args; DEFAULT_SCALAR_INPUTS
source(file.path(FN, "calculate_model_approx_r0.R"))         # compute_R0_invariants(); `%||%`
source(file.path(FN, "abc_calibration_functions_common.R"))  # generic ABC helpers
source(file.path(FN, "abc_calibration_functions_npi.R"))     # NPI build_abc_model_args(); npi_efficacy_from_scaler()
source(file.path(FN, "abc_posterior.R"))                     # downsample_posterior()
handlers("progress")
check_model_function_version()   # fail fast if fiber predates the NPI interface


# ---------------------------------------------------------------------------
# 2. READ POSTERIORS + RESOLVE EACH FIT'S npi_spec, THEN PICK THE COMMON SIZE
# ---------------------------------------------------------------------------
# The fit files are EasyABC result objects: $param (3 columns: R0, prop_funeral,
# npi_scaler) and $weights. Rebuild the weighted-posterior data frame
# downsample_posterior() expects, and resolve the npi_spec for each scenario in
# priority order: embedded in the .rds -> source-able "<stem>_metadata.txt"
# sidecar -> the per-scenario script fallback. npi_scaler is uninterpretable
# without it, so we record which source was used.
posteriors      <- list()
npi_specs       <- list()
npi_spec_source <- list()
for (sc in names(SCENARIOS)) {
  rds_path <- here::here(SCENARIOS[[sc]]$rds)
  if (!file.exists(rds_path)) stop("Fit RDS not found: ", rds_path, call. = FALSE)
  res <- readRDS(rds_path)

  posteriors[[sc]] <- data.frame(
    weight       = res$weights,
    R0           = res$param[, 1],   ## fitted: baseline R0
    prop_funeral = res$param[, 2],   ## fitted: funeral transmission share
    npi_scaler   = res$param[, 3]    ## fitted: shared NPI-efficacy position in [-1, 1]
  )

  spec <- res$npi_spec                       # (1) embedded by newer calibration runs
  src  <- "embedded in .rds"
  if (is.null(spec)) {                       # (2) source-able sidecar next to the .rds
    side <- sub("\\.rds$", "_metadata.txt", rds_path, ignore.case = TRUE)
    if (file.exists(side)) {
      e  <- new.env()
      ok <- tryCatch({ sys.source(side, envir = e); TRUE }, error = function(...) FALSE)
      if (ok && !is.null(e$npi_run_metadata$npi_spec)) {
        spec <- e$npi_run_metadata$npi_spec
        src  <- "sidecar _metadata.txt"
      }
    }
  }
  if (is.null(spec)) {                       # (3) per-scenario script fallback
    spec <- SCENARIOS[[sc]]$npi_spec_fallback
    src  <- "script fallback (user-confirmed)"
  }
  if (is.null(spec)) stop("No npi_spec available for scenario '", sc, "'.", call. = FALSE)

  npi_specs[[sc]]       <- spec
  npi_spec_source[[sc]] <- src
  message(sprintf("[%s] npi_spec from %s: ppe[%.2f, %.2f], etu[%.2f, %.2f]",
                  sc, src,
                  spec$ppe_efficacy$min, spec$ppe_efficacy$max,
                  spec$etu_efficacy$min, spec$etu_efficacy$max))
}
N_SETS <- min(50, ## lower this for a quicker run
              min(vapply(posteriors, nrow, integer(1))))   # the scenario with the fewest draws
message(sprintf("Particles available: %s; using N_SETS = %d.",
                paste(sprintf("%s=%d", names(posteriors), vapply(posteriors, nrow, integer(1))), collapse = ", "), N_SETS))


# ---------------------------------------------------------------------------
# 3. PER SCENARIO: INVARIANTS, DOWNSAMPLE, PRE-BUILD ONE ARG LIST PER PARTICLE
# ---------------------------------------------------------------------------
# Build the fiber argument lists in the main session (not on the workers) so the
# conversion stays provably identical to the calibration; workers only run fiber.
# Each scenario may use a different scenario CSV, so read it inside the loop.
args_by_scenario <- list()
for (sc in names(SCENARIOS)) {
  scenario_matrix <- read_scenario_matrix(here::here("data-processed", SCENARIOS[[sc]]$csv))
  mp <- make_model_parameters(SCENARIOS[[sc]]$id, scenario_matrix,
                              overrides = list(check_final_size = CHECK_FINAL_SIZE))
  # Efficacy-INDEPENDENT R0 invariants for THIS scenario (computed once). The NPI
  # build_abc_model_args() recomputes the per-particle D/F from these + the
  # particle's efficacies (which it derives from npi_scaler via npi_spec).
  inv <- compute_R0_invariants(args = mp$args, n = SETUP_R0_N, seed = SETUP_R0_SEED)
  theta <- downsample_posterior(posteriors[[sc]], n_sets = N_SETS, seed = RESAMPLE_SEED,
                                param_names = c("R0", "prop_funeral", "npi_scaler"))
  # One ready-to-run fiber arg list per particle (OBV is OFF here; the worker
  # switches it on per efficacy below). The fixed efficacies match the fit.
  args_by_scenario[[sc]] <- lapply(seq_len(N_SETS), function(i)
    build_abc_model_args(
      R0 = theta$R0[i], prop_funeral = theta$prop_funeral[i],
      npi_scaler = theta$npi_scaler[i],
      base = mp$base_args, tv = mp$tv_args, invariants = inv,
      npi_spec = npi_specs[[sc]],
      general_hospital_quarantine_efficacy = DEFAULT_SCALAR_INPUTS$general_hospital_quarantine_efficacy,
      safe_funeral_efficacy = DEFAULT_SCALAR_INPUTS$safe_funeral_efficacy,
      seeding_cases = SEEDING_CASES))
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

# Friendly labels (scenario label from the SCENARIOS config; efficacy from value).
# as.character() guards against aggregate() returning these grouping cols as factors.
agg$scenario_lbl <- vapply(as.character(agg$scenario),
                           function(sc) SCENARIOS[[sc]]$label, character(1))
agg$efficacy_lbl <- sprintf("%d%% efficacy",
                            round(100 * OBV_EFFICACIES[as.character(agg$efficacy)]))

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
p_plot <- ggplot(agg[!is.na(agg$pct_averted), ],
       aes(scenario_lbl, pct_averted, fill = efficacy_lbl)) +
  geom_boxplot(outlier.size = 0.5, position = position_dodge(0.8), width = 0.7) +
  labs(x = NULL, y = "HCW deaths averted (%)", fill = "Obeldesivir",
       title = "Obeldesivir impact on HCW deaths, by scenario and efficacy",
       subtitle = sprintf("NPI-efficacy fits; per-run counterfactual (direct prevented deaths); %d posterior draws x %d reps per efficacy",
                          N_SETS, N_REPS)) +
  theme_bw()
print(p_plot)


# ---------------------------------------------------------------------------
# 7. SAVE OUTPUTS (plot + tables + the npi_spec actually used, for provenance)
# ---------------------------------------------------------------------------
out_dir <- here::here(OUTPUT_DIR)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(out_dir, "obv_pct_hcw_deaths_averted.png"),
       plot = p_plot, width = 8, height = 5, dpi = 150)
write.csv(agg,  file.path(out_dir, "obv_per_particle_pct_averted.csv"), row.names = FALSE)
write.csv(summ, file.path(out_dir, "obv_summary_pct_averted.csv"),      row.names = FALSE)

# Record exactly how these numbers were produced -- crucially, the npi_spec each
# scenario was interpreted with (and where it came from). source()-able.
obv_run_metadata <- list(
  generated_at     = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
  scenarios        = setNames(lapply(names(SCENARIOS), function(sc) list(
                        id              = SCENARIOS[[sc]]$id,
                        rds             = SCENARIOS[[sc]]$rds,
                        npi_spec        = npi_specs[[sc]],
                        npi_spec_source = npi_spec_source[[sc]])), names(SCENARIOS)),
  obv_efficacies   = OBV_EFFICACIES,
  obv_base         = OBV_BASE,
  n_sets           = N_SETS,
  n_reps           = N_REPS,
  seeding_cases    = SEEDING_CASES,
  check_final_size = CHECK_FINAL_SIZE,
  fixed_efficacies = list(
    general_hospital_quarantine_efficacy = DEFAULT_SCALAR_INPUTS$general_hospital_quarantine_efficacy,
    safe_funeral_efficacy                = DEFAULT_SCALAR_INPUTS$safe_funeral_efficacy)
)
con <- file(file.path(out_dir, "obv_run_metadata.txt"), open = "wt")
writeLines("# Provenance for the OBV direct-effect results in this directory. source()-able.", con)
cat("obv_run_metadata <- ", file = con); dput(obv_run_metadata, file = con)
close(con)

message("Saved plot + tables + provenance to ", out_dir)
