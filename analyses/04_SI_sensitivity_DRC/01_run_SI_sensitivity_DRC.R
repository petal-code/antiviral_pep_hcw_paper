# =============================================================================
# 01_run_SI_sensitivity_DRC.R
#
# Supplementary-information (SI) sensitivity analyses for the DRC-like archetype.
#
#   (1) "Vaccine-free stress test" -- increase transmissibility.
#       Vaccination in the observed DRC outbreak likely suppressed transmission,
#       so a posterior calibrated to the realised dynamics may be conservative
#       for a vaccine-UNAVAILABLE setting. We rerun the DRC-like posterior
#       simulations with baseline transmissibility scaled UP by +10/+20/+30%
#       (via R0, equivalently the offspring means -- see note), keeping the PEP
#       scenarios fixed.
#       Report: HCW deaths in the no-PEP baseline, HCW deaths averted, % reduction.
#
#   (2) "HCW-exposure upshift" -- increase the HCW exposure scalar.
#       Vaccination of frontline responders in North Kivu/Ituri may mean the
#       realised HCW burden partly reflects protection of response personnel. We
#       rerun with the healthcare-worker exposure scalar (hcw_risk_scalar) scaled
#       UP by +10/+25/+50%.
#       Report: baseline HCW deaths, HCW deaths averted, HCW-days lost averted,
#       % reduction.
#
# This script reuses the EXACT machinery, posterior, and scenario inputs that
# generated the main figures (cf. analyses/03_figure_template/01_analysis_figure2.R):
#   * the same DRC posterior .rds  (Middle_DRC_ConflictSmoothed_PlusPlus, the
#     20260607_215621 NS4 fit wired into figures 1-5),
#   * the same original-approach scenario CSV + scenario id,
#   * the same build_abc_model_args_decoupled() argument builder,
#   * the same obv_pep_* PEP settings, takeoff-retry loop, seeds, N_PARTICLES,
#     N_REPS, and HCW_BASE_PROB.
# The ONLY change is a per-particle multiplicative scaling of R0 (analysis 1) or
# hcw_risk_scalar (analysis 2), injected at the build step. Everything else --
# priors, fixed efficacies, the PEP settings -- is identical to the main
# analysis, so the results are directly comparable to the figures.
#
# WHY scaling R0 == scaling the offspring means.
#   solve_offspring_means(R0, prop_funeral, D, F) returns
#       mn_offspring_genPop  = (1 - prop_funeral) * R0 / D
#       mn_offspring_funeral =      prop_funeral  * R0 / F
#   and the multipliers D, F do NOT depend on R0 (calculate_model_approx_r0.R).
#   So multiplying a particle's R0 by `f` multiplies BOTH offspring means by
#   exactly `f` -- a clean `f`x increase in baseline transmissibility that holds
#   the funeral share (prop_funeral) and the natural-history structure fixed.
#
# HOW the no-PEP baseline is obtained (explicit paired runs).
#   For each (scaling, particle, rep) we run a no-PEP simulation
#   (obv_pep_enabled = FALSE) and the with-PEP arms at the SAME takeoff seed (the
#   variance-reduction pairing in functions/simulation_helpers.R), then take HCW
#   deaths averted = baseline - with-PEP as a matched difference. We do NOT
#   reconstruct the baseline from out$prevented_completed (as the figures do):
#   that field is empty under some fiber builds and silently collapses the
#   baseline onto the with-PEP tdf (averted = 0). prevented_completed is still
#   recorded as a cross-check.
#
# Outputs:
#   outputs/04_SI_sensitivity_DRC/SI_sensitivity_DRC_run_summary.csv  (run-level;
#       one row per scaling x arm x particle x rep -- heavy, gitignored)
#   output_figgen/SI_sensitivity_DRC_particle_df.csv  (compact, per-particle means;
#       the analysable unit, tracked in git)
#   output_figgen/SI_sensitivity_DRC_hcw_saturation.csv  (diagnostic: fraction of
#       particles whose scaled prob_hcw saturates the pmin(., 1) cap)
# Summary tables + SI figures are produced by 02_summarise_SI_sensitivity_DRC.R.
#
# Requires the `fiber` package: devtools::install_github("petal-code/fiber").
# =============================================================================

library(here)
library(future)
library(future.apply)
library(dplyr)
library(fiber)

source(here("functions", "setup_model_parameters.R"))
source(here("functions", "calculate_model_approx_r0.R"))
source(here("functions", "abc_calibration_functions_common.R"))
source(here("functions", "abc_calibration_functions_decoupled.R"))
source(here("functions", "abc_posterior.R"))

# =============================================================================
# Configuration  (defaults match the main figure analysis exactly)
# =============================================================================
N_PARTICLES   <- 200L     # posterior particles (resampled); main analysis = 200
N_REPS        <- 10L      # stochastic replicates per particle/arm; main = 10
N_WORKERS     <- min(120L, future::availableCores())
SEEDING_CASES <- 25L
RESAMPLE_SEED <- 42L      # posterior downsample seed (matches figures)
SEED_BASE     <- 20260801L
HCW_BASE_PROB <- 0.25     # hcw_risk_scalar multiplies this (matches figures)

TAKEOFF_DEATH_THRESHOLD <- 100L
MAX_RETRIES             <- 50L

# HCW-days lost counts every infected HCW as absent from symptom onset to the end
# of that run's epidemic (the figure-2/3 obv_return = FALSE convention), computed
# per run from its own tdf.

# --- Quick smoke-test override (uncomment to sanity-check the plumbing fast) ---
# N_PARTICLES <- 20L; N_REPS <- 2L

# =============================================================================
# Sensitivity scalings.
# Each entry scales one calibrated quantity per particle, leaving everything else
# at its posterior value. The `reference` cell (both factors = 1) is the
# as-fitted DRC archetype; it is shared by BOTH analyses in the summary step, so
# it is only simulated once here.
# =============================================================================
SENS_SCENARIOS <- list(
  list(key = "reference",          analysis = "reference",       label = "DRC as-fitted (x1.00)",
       r0_factor = 1.00, hcw_factor = 1.00),

  list(key = "transmissibility_110", analysis = "transmissibility", label = "Transmissibility +10%",
       r0_factor = 1.10, hcw_factor = 1.00),
  list(key = "transmissibility_120", analysis = "transmissibility", label = "Transmissibility +20%",
       r0_factor = 1.20, hcw_factor = 1.00),
  list(key = "transmissibility_130", analysis = "transmissibility", label = "Transmissibility +30%",
       r0_factor = 1.30, hcw_factor = 1.00),

  list(key = "hcw_exposure_125",   analysis = "hcw_exposure",     label = "HCW exposure +25%",
       r0_factor = 1.00, hcw_factor = 1.25),
  list(key = "hcw_exposure_150",   analysis = "hcw_exposure",     label = "HCW exposure +50%",
       r0_factor = 1.00, hcw_factor = 1.50),
  list(key = "hcw_exposure_200",   analysis = "hcw_exposure",     label = "HCW exposure +100%",
       r0_factor = 1.00, hcw_factor = 2.00)
)

# =============================================================================
# PEP scenarios (arms).
# "Keeping the same PEP scenarios": the full-coverage antiviral-efficacy sweep
# used in Figure 2 (the headline HCW-targeted PEP results: 100% coverage,
# 50-90% efficacy). full_obv80 is the headline arm. To extend to the Figure-3
# coverage scenarios, add ramp_high / ramp_low entries here and supply the
# matching coverage curves in COVERAGE_FNS (see helper_functions_figure_1to4.R
# COVERAGE_SPECS for the canonical splines).
# =============================================================================
ARM_EFFICACIES <- 0.80
ARMS <- data.frame(
  arm_name = sprintf("full_obv%02d", round(ARM_EFFICACIES * 100)),
  coverage = "full",
  efficacy = ARM_EFFICACIES,
  stringsAsFactors = FALSE
)

COVERAGE_FNS <- list(
  full = function(t) rep(1.0, length(t))
)

# Resolve a path, falling back to a case-insensitive match in its directory.
# The figures reference the posterior as ".RDS" but it is saved as ".rds"; this
# keeps the script working on case-sensitive (Linux) and case-insensitive (macOS)
# filesystems alike.
resolve_existing <- function(path) {
  if (file.exists(path)) return(path)
  hits <- list.files(dirname(path), full.names = TRUE)
  m <- hits[tolower(basename(hits)) == tolower(basename(path))]
  if (length(m) >= 1L) return(m[[1L]])
  stop("Cannot find file: ", path, call. = FALSE)
}

# =============================================================================
# DRC scenario definition (identical to the figure scripts).
# =============================================================================
DRC <- list(
  id           = "Middle_DRC_ConflictSmoothed_PlusPlus",
  rds          = resolve_existing(here("outputs", "02_ABC_model_fits_Final",
                      "fiber_ABC_SMC_Middle_DRC_ConflictSmoothed_PlusPlus_Decoupled_20260607_215621_check_NP5_NS4_NBREPS_30_NBSIMUL_590.rds")),
  scenario_csv = here("data-processed", "final_six_scenario_values_original_approach.csv"),
  param_names  = c("R0", "prop_funeral", "etu_efficacy", "ppe_efficacy", "hcw_risk_scalar"),
  check_final_size = 15000
)

OUT_DIR_HEAVY <- here("outputs", "04_SI_sensitivity_DRC")
OUT_DIR_FIG   <- here("output_figgen")
dir.create(OUT_DIR_HEAVY, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR_FIG,   recursive = TRUE, showWarnings = FALSE)

check_model_function_version()

# Absolute paths to the function files, so PSOCK workers can re-source them
# (here::here() may not resolve inside a fresh worker session).
FN_PATHS <- list(
  setup     = here("functions", "setup_model_parameters.R"),
  r0        = here("functions", "calculate_model_approx_r0.R"),
  common    = here("functions", "abc_calibration_functions_common.R"),
  decoupled = here("functions", "abc_calibration_functions_decoupled.R")
)

# =============================================================================
# Per-scenario setup: load posterior, downsample, build base args + invariants.
# (Mirrors analyses/03_figure_template/01_analysis_figure2.R lines 109-130.)
# =============================================================================
message("Loading DRC posterior + building model parameters ...")
res       <- readRDS(DRC$rds)
posterior <- as.data.frame(res$param)
colnames(posterior) <- DRC$param_names
posterior$weight    <- res$weights

theta <- downsample_posterior(posterior, n_sets = N_PARTICLES,
                              seed = RESAMPLE_SEED, param_names = DRC$param_names)

scenario_matrix <- read_scenario_matrix(DRC$scenario_csv)
mp <- make_model_parameters(
  scenario_id     = DRC$id,
  scenario_matrix = scenario_matrix,
  overrides       = list(seeding_cases    = SEEDING_CASES,
                         check_final_size = DRC$check_final_size)
)
# Efficacy-independent R0 invariants -- do NOT depend on R0 or hcw_risk_scalar,
# so they are computed once and reused across all scalings.
inv <- compute_R0_invariants(args = mp$args, n = 50000, seed = 42L)

base_args <- mp$base_args
tv_args   <- mp$tv_args

# =============================================================================
# hcw_risk_scalar saturation diagnostic.
# prob_hcw = pmin(HCW_BASE_PROB * hcw_risk_scalar * factor, 1). For the +50% arm,
# particles with hcw_risk_scalar above ~ 1 / (HCW_BASE_PROB * 1.5) = 2.67 hit the
# cap. We record how many do, so the interpretation can flag any saturation.
# =============================================================================
hcw_factors <- sort(unique(vapply(SENS_SCENARIOS, function(s) s$hcw_factor, numeric(1))))
saturation_df <- do.call(rbind, lapply(hcw_factors, function(hf) {
  raw    <- HCW_BASE_PROB * theta$hcw_risk_scalar * hf
  capped <- pmin(raw, 1.0)
  data.frame(
    hcw_factor               = hf,
    frac_particles_saturated = mean(raw >= 1.0),
    median_prob_hcw          = median(capped),
    max_uncapped_prob_hcw    = max(raw),
    stringsAsFactors = FALSE
  )
}))
write.csv(saturation_df,
          file.path(OUT_DIR_FIG, "SI_sensitivity_DRC_hcw_saturation.csv"),
          row.names = FALSE)
message("hcw_risk_scalar saturation by factor:")
print(saturation_df, row.names = FALSE)

# =============================================================================
# Per-run HCW outcomes from a single fiber output (tdf only -- no reliance on
# out$prevented_completed). HCW-days lost uses the figure-2/3 convention
# (obv_return = FALSE): every infected HCW is absent from symptom onset to the
# end of that run's epidemic.
# =============================================================================
hcw_metrics <- function(out) {
  tdf <- out$tdf
  tdf <- tdf[!is.na(tdf$time_infection_absolute), , drop = FALSE]
  duration <- max(tdf$time_outcome_absolute, na.rm = TRUE)

  is_hcw <- !is.na(tdf$class)   & tdf$class == "HCW"
  died   <- !is.na(tdf$outcome) & tdf$outcome

  hcw <- tdf[is_hcw, , drop = FALSE]
  days_lost <- if (nrow(hcw) == 0) 0 else
    sum(pmax(duration - hcw$time_symptom_onset_absolute, 0), na.rm = TRUE)

  list(
    n_hcw_deaths = sum(is_hcw & died),
    days_lost    = days_lost,
    total_deaths = sum(died)
  )
}

# Diagnostic: how many cases OBV recorded as prevented (line-list + counter).
# Reported alongside the paired estimates so the prevented_completed channel can
# be cross-checked, but it is NOT used to compute the headline outcomes.
obv_prevented_diag <- function(out) {
  pc <- out$prevented_completed
  n_rows <- if (!is.null(pc) && nrow(pc) > 0) nrow(pc) else 0L
  n_treated <- tryCatch(out$sim_info$obv_pep_num_treated$prevented, error = function(e) NA_real_)
  list(n_rows = n_rows, n_treated = if (is.null(n_treated)) NA_real_ else n_treated)
}

# =============================================================================
# Job list: one job per worker, each handling a contiguous block of particles.
# =============================================================================
particle_chunks <- split(seq_len(N_PARTICLES),
                         cut(seq_len(N_PARTICLES), N_WORKERS, labels = FALSE))
particle_chunks <- Filter(function(x) length(x) > 0, particle_chunks)

n_arms <- nrow(ARMS)
n_sens <- length(SENS_SCENARIOS)
# Per (scaling, particle, rep): 1 no-PEP baseline + n_arms with-PEP runs.
total_runs <- n_sens * N_PARTICLES * N_REPS * (n_arms + 1L)
message(sprintf(
  "%d scalings x %d particles x %d reps x (1 baseline + %d arms) = %d simulations on %d workers.",
  n_sens, N_PARTICLES, N_REPS, n_arms, total_runs, length(particle_chunks)))

# =============================================================================
# Worker: for a block of particles, run -- per (scaling, rep) -- a no-PEP
# baseline and the with-PEP arms PAIRED at the same takeoff seed, extract HCW
# outcomes inline, and return a tidy data.frame.
#
# Why paired explicit runs (not tdf + prevented_completed). The no-PEP
# counterfactual is obtained from its OWN run (obv_pep_enabled = FALSE) rather
# than reconstructed from out$prevented_completed, which is empty under some
# fiber builds and silently collapses the baseline onto the with-PEP tdf. The
# baseline and every arm share the same seed (the variance-reduction pairing in
# simulation_helpers.R), so HCW deaths averted = baseline - with-PEP is a matched
# difference on the same stochastic realisation.
# =============================================================================
run_particle_block <- function(particle_ids) {
  # Re-source the model code once per fresh worker session.
  if (!isTRUE(get0(".si_worker_ready", envir = globalenv()))) {
    source(FN_PATHS$setup)
    source(FN_PATHS$r0)
    source(FN_PATHS$common)
    source(FN_PATHS$decoupled)
    if (!"fiber" %in% loadedNamespaces()) library(fiber)
    assign(".si_worker_ready", TRUE, envir = globalenv())
  }

  rows <- vector("list", 0L)

  for (p in particle_ids) {
    for (s_idx in seq_along(SENS_SCENARIOS)) {
      s <- SENS_SCENARIOS[[s_idx]]

      # Build the particle args with the scaled calibrated quantity. Identical to
      # the figure scripts, except R0 and hcw_risk_scalar carry the scaling.
      args_base <- build_abc_model_args_decoupled(
        R0              = theta$R0[p]              * s$r0_factor,
        prop_funeral    = theta$prop_funeral[p],
        etu_efficacy    = theta$etu_efficacy[p],
        ppe_efficacy    = theta$ppe_efficacy[p],
        hcw_risk_scalar = theta$hcw_risk_scalar[p] * s$hcw_factor,
        base            = base_args,
        tv              = tv_args,
        invariants      = inv,
        general_hospital_quarantine_efficacy = DEFAULT_SCALAR_INPUTS$general_hospital_quarantine_efficacy,
        safe_funeral_efficacy                = DEFAULT_SCALAR_INPUTS$safe_funeral_efficacy,
        hcw_base_prob   = HCW_BASE_PROB,
        seeding_cases   = SEEDING_CASES
      )
      args_base$check_final_size <- DRC$check_final_size

      for (r in seq_len(N_REPS)) {
        # One seed per (scaling, particle, rep); shared by baseline and all arms.
        seed0 <- SEED_BASE +
          (s_idx - 1L) * (N_PARTICLES * N_REPS) +
          (p     - 1L) * N_REPS +
          (r     - 1L)

        # --- no-PEP baseline: retry until takeoff, then fix the seed ---
        a0 <- args_base
        a0$obv_pep_enabled <- FALSE
        retry <- 0L; takeoff_seed <- seed0; out0 <- NULL
        repeat {
          a0$seed <- takeoff_seed
          out0    <- do.call(fiber::branching_process_main, a0)
          tdf0    <- out0$tdf[!is.na(out0$tdf$time_infection_absolute), ]
          base_deaths <- sum(!is.na(tdf0$outcome) & tdf0$outcome)
          if (base_deaths >= TAKEOFF_DEATH_THRESHOLD) break
          retry        <- retry + 1L
          takeoff_seed <- takeoff_seed + 1L
          if (retry >= MAX_RETRIES) break
        }
        base_m <- hcw_metrics(out0)
        took_off <- base_m$total_deaths >= TAKEOFF_DEATH_THRESHOLD

        # --- with-PEP arms at the SAME takeoff seed (paired) ---
        for (arm_idx in seq_len(nrow(ARMS))) {
          arm <- ARMS[arm_idx, ]
          a1  <- args_base
          a1$obv_pep_enabled          <- TRUE
          a1$obv_pep_coverage         <- COVERAGE_FNS[[arm$coverage]]
          a1$obv_pep_adherence        <- 1.0
          a1$obv_pep_dpc              <- 0
          a1$obv_pep_efficacy         <- arm$efficacy
          a1$obv_pep_target_class     <- "HCW"
          a1$obv_pep_target_locations <- c("hospital", "community", "funeral")
          a1$seed                     <- takeoff_seed

          out1  <- do.call(fiber::branching_process_main, a1)
          arm_m <- hcw_metrics(out1)
          diag  <- obv_prevented_diag(out1)

          rows[[length(rows) + 1L]] <- data.frame(
            analysis           = s$analysis,
            scenario_key       = s$key,
            scenario_label     = s$label,
            r0_factor          = s$r0_factor,
            hcw_factor         = s$hcw_factor,
            arm                = arm$arm_name,
            efficacy           = arm$efficacy,
            particle_id        = p,
            rep                = r,
            took_off           = took_off,
            n_hcw_deaths       = arm_m$n_hcw_deaths,                 # WITH PEP
            counterfactual_hcw = base_m$n_hcw_deaths,                # no-PEP baseline
            prevented_hcw      = base_m$n_hcw_deaths - arm_m$n_hcw_deaths,
            hcw_days_lost      = arm_m$days_lost,                    # WITH PEP
            baseline_days_lost = base_m$days_lost,                   # no-PEP baseline
            total_deaths_obv   = arm_m$total_deaths,
            total_deaths_base  = base_m$total_deaths,
            obv_prevented_rows = diag$n_rows,                        # diagnostic only
            obv_num_treated    = diag$n_treated,                     # diagnostic only
            stringsAsFactors   = FALSE
          )
        }
      }
    }
  }
  do.call(rbind, rows)
}

# =============================================================================
# Run in parallel over particle blocks.
# =============================================================================
plan(multisession, workers = length(particle_chunks))
t_start <- proc.time()

results <- future_lapply(
  particle_chunks, run_particle_block,
  future.globals = list(
    theta = theta, base_args = base_args, tv_args = tv_args, inv = inv,
    ARMS = ARMS, COVERAGE_FNS = COVERAGE_FNS, SENS_SCENARIOS = SENS_SCENARIOS,
    DRC = DRC, FN_PATHS = FN_PATHS,
    hcw_metrics = hcw_metrics, obv_prevented_diag = obv_prevented_diag,
    run_particle_block = run_particle_block,
    N_REPS = N_REPS, N_PARTICLES = N_PARTICLES, SEED_BASE = SEED_BASE,
    SEEDING_CASES = SEEDING_CASES, HCW_BASE_PROB = HCW_BASE_PROB,
    TAKEOFF_DEATH_THRESHOLD = TAKEOFF_DEATH_THRESHOLD,
    MAX_RETRIES = MAX_RETRIES
  ),
  future.packages = c("fiber"),
  future.seed     = TRUE
)

plan(sequential)
run_df <- do.call(rbind, results)
elapsed <- (proc.time() - t_start)["elapsed"]
message(sprintf("Done %d simulations in %.1f minutes.", nrow(run_df), elapsed / 60))

# =============================================================================
# Sanity check: OBV impact must be non-zero, and report whether the
# prevented_completed channel was populated (cross-check only -- the headline
# numbers come from the paired no-PEP baseline, not from prevented_completed).
# =============================================================================
sanity <- run_df %>%
  group_by(analysis) %>%
  summarise(
    mean_baseline_hcw_deaths = mean(counterfactual_hcw),
    mean_pep_hcw_deaths      = mean(n_hcw_deaths),
    mean_hcw_deaths_averted  = mean(prevented_hcw),
    mean_prevented_rows      = mean(obv_prevented_rows),
    .groups = "drop"
  )
message("Per-analysis sanity (paired no-PEP vs with-PEP):")
print(as.data.frame(sanity), row.names = FALSE)
if (max(sanity$mean_hcw_deaths_averted, na.rm = TRUE) <= 0)
  warning("HCW deaths averted is <= 0 everywhere -- check the OBV/PEP settings.")
if (all(sanity$mean_prevented_rows == 0))
  message("NB: out$prevented_completed was empty in every run (expected with this ",
          "fiber build); the paired-baseline method does not rely on it.")

# =============================================================================
# Save run-level results (heavy; gitignored under outputs/).
# =============================================================================
run_csv <- file.path(OUT_DIR_HEAVY, "SI_sensitivity_DRC_run_summary.csv")
write.csv(run_df, run_csv, row.names = FALSE)
message("Saved run-level results: ", run_csv)

# =============================================================================
# Aggregate over replicates -> one row per (scaling x arm x particle).
# Means over reps put every quantity on the per-epidemic scale; the percentage
# reductions are ratios of those means and so match the figures exactly
# (make_particle_df()'s ratio-of-sums == ratio-of-means).
# =============================================================================
particle_df <- run_df %>%
  group_by(analysis, scenario_key, scenario_label, r0_factor, hcw_factor,
           arm, efficacy, particle_id) %>%
  summarise(
    n_reps                 = dplyr::n(),
    baseline_hcw_deaths    = mean(counterfactual_hcw),
    pep_hcw_deaths         = mean(n_hcw_deaths),
    hcw_deaths_averted     = mean(prevented_hcw),
    baseline_hcw_days_lost = mean(baseline_days_lost),
    pep_hcw_days_lost      = mean(hcw_days_lost),
    hcw_days_lost_averted  = mean(baseline_days_lost) - mean(hcw_days_lost),
    .groups = "drop"
  ) %>%
  mutate(
    pct_hcw_deaths_averted = ifelse(baseline_hcw_deaths > 0,
                                    100 * hcw_deaths_averted / baseline_hcw_deaths,
                                    NA_real_),
    pct_days_lost_averted  = ifelse(baseline_hcw_days_lost > 0,
                                    100 * hcw_days_lost_averted / baseline_hcw_days_lost,
                                    NA_real_)
  )

particle_csv <- file.path(OUT_DIR_FIG, "SI_sensitivity_DRC_particle_df.csv")
write.csv(particle_df, particle_csv, row.names = FALSE)
message("Saved per-particle results: ", particle_csv)
message("Run 02_summarise_SI_sensitivity_DRC.R to build SI tables + figures.")
