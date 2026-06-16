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
# HOW the no-PEP baseline is obtained (same convention as the figures).
#   Every simulation is run with OBV/PEP enabled; the matched no-PEP
#   counterfactual is reconstructed as cases_base = tdf + prevented_completed (the
#   cases OBV prevented, added back) -- exactly as helper_functions_figure_1to4.R
#   / extract_run_summary() do. prevented_completed MUST be populated; the run
#   records its row count and fiber's obv_pep_num_treated$prevented and errors if
#   it is empty everywhere (rather than silently reporting averted = 0).
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
# Pre-build the REFERENCE per-particle args in the MAIN process, exactly as the
# figure scripts do (build_abc_model_args_decoupled pulled from a dedicated
# environment via sys.source). These carry no scaling (r0_factor = hcw_factor = 1)
# and no OBV settings; the worker derives the scaled args from them by arithmetic
# and switches OBV on per arm.
# =============================================================================
env_decoupled <- new.env(parent = globalenv())
sys.source(here("functions", "abc_calibration_functions_decoupled.R"), envir = env_decoupled)
build_decoupled <- env_decoupled$build_abc_model_args_decoupled

base_particle_args <- lapply(seq_len(N_PARTICLES), function(p) {
  args <- build_decoupled(
    R0              = theta$R0[p],
    prop_funeral    = theta$prop_funeral[p],
    etu_efficacy    = theta$etu_efficacy[p],
    ppe_efficacy    = theta$ppe_efficacy[p],
    hcw_risk_scalar = theta$hcw_risk_scalar[p],
    base            = base_args,
    tv              = tv_args,
    invariants      = inv,
    general_hospital_quarantine_efficacy = DEFAULT_SCALAR_INPUTS$general_hospital_quarantine_efficacy,
    safe_funeral_efficacy                = DEFAULT_SCALAR_INPUTS$safe_funeral_efficacy,
    hcw_base_prob   = HCW_BASE_PROB,
    seeding_cases   = SEEDING_CASES
  )
  args$check_final_size <- DRC$check_final_size
  args
})

# Scale a reference particle's args by (r0_factor, hcw_factor). Both are exact:
#   mn_offspring_{genPop,funeral} are linear in R0, so x r0_factor;
#   prob_hcw_cond_*_hospital = pmin(hcw_base_prob * hcw_risk_scalar, 1), so the
#   already-built (reference, unclamped) value x hcw_factor, re-clamped at 1.
scale_particle_args <- function(args, r0_factor, hcw_factor) {
  args$mn_offspring_genPop  <- args$mn_offspring_genPop  * r0_factor
  args$mn_offspring_funeral <- args$mn_offspring_funeral * r0_factor
  args$prob_hcw_cond_genPop_hospital <- pmin(args$prob_hcw_cond_genPop_hospital * hcw_factor, 1.0)
  args$prob_hcw_cond_hcw_hospital    <- pmin(args$prob_hcw_cond_hcw_hospital    * hcw_factor, 1.0)
  args
}

# =============================================================================
# Per-run HCW outcomes (in-memory analogue of extract_run_summary()).
# The no-PEP baseline is reconstructed as cases_base = tdf + prevented_completed,
# EXACTLY as the figures do (helper_functions_figure_1to4.R). We also return the
# prevented_completed row count and fiber's own obv_pep_num_treated$prevented so
# an empty channel is flagged rather than silently collapsing the baseline.
# =============================================================================
summarise_one_run <- function(out) {
  tdf       <- out$tdf
  prevented <- out$prevented_completed
  n_prevented_rows <- if (!is.null(prevented) && nrow(prevented) > 0) nrow(prevented) else 0L
  n_treated <- tryCatch({
    v <- out$sim_info$obv_pep_num_treated$prevented
    if (is.null(v)) NA_real_ else as.numeric(v)
  }, error = function(e) NA_real_)

  duration <- max(tdf$time_outcome_absolute, na.rm = TRUE)
  if (n_prevented_rows > 0) {
    missing_cols <- setdiff(names(tdf), names(prevented))
    for (col in missing_cols) prevented[[col]] <- NA
    prevented  <- prevented[, names(tdf), drop = FALSE]
    cases_base <- rbind(tdf, prevented)              # no-PEP counterfactual
  } else {
    cases_base <- tdf
  }

  .hcw_deaths <- function(cases) {
    is_hcw <- !is.na(cases$class)   & cases$class == "HCW"
    died   <- !is.na(cases$outcome) & cases$outcome
    sum(died & is_hcw)
  }
  .days_lost <- function(cases) {
    hcw <- cases[!is.na(cases$class) & cases$class == "HCW", ]
    if (nrow(hcw) == 0) return(0)
    sum(pmax(duration - hcw$time_symptom_onset_absolute, 0), na.rm = TRUE)
  }

  n_base <- .hcw_deaths(cases_base)
  n_obv  <- .hcw_deaths(tdf)
  list(
    n_hcw_deaths       = n_obv,                  # HCW deaths WITH PEP
    hcw_days_lost      = .days_lost(tdf),        # HCW-days lost WITH PEP
    counterfactual_hcw = n_base,                 # HCW deaths in no-PEP baseline
    prevented_hcw      = n_base - n_obv,         # HCW deaths averted
    baseline_days_lost = .days_lost(cases_base), # HCW-days lost in no-PEP baseline
    total_deaths       = sum(!is.na(tdf$outcome) & tdf$outcome),
    n_prevented_rows   = n_prevented_rows,       # diagnostic: prevented_completed rows
    obv_num_treated    = n_treated               # diagnostic: fiber's prevented count
  )
}

# =============================================================================
# Job list: one job per worker, each handling a contiguous block of particles.
# =============================================================================
particle_chunks <- split(seq_len(N_PARTICLES),
                         cut(seq_len(N_PARTICLES), N_WORKERS, labels = FALSE))
particle_chunks <- Filter(function(x) length(x) > 0, particle_chunks)

n_arms <- nrow(ARMS)
n_sens <- length(SENS_SCENARIOS)
total_runs <- n_sens * n_arms * N_PARTICLES * N_REPS
message(sprintf(
  "%d scalings x %d arms x %d particles x %d reps = %d simulations on %d workers.",
  n_sens, n_arms, N_PARTICLES, N_REPS, total_runs, length(particle_chunks)))

# =============================================================================
# Worker: for a block of particles, run every (scaling x arm x rep) with OBV
# enabled, reconstruct the no-PEP baseline from tdf + prevented_completed (the
# figure convention), and return a tidy data.frame. Mirrors the figure worker:
# it consumes the MAIN-process-built reference args and only scales them, sets
# OBV, and runs -- no model code is sourced inside the worker.
# =============================================================================
run_particle_block <- function(particle_ids) {
  rows <- vector("list", 0L)

  for (p in particle_ids) {
    ref <- base_particle_args[[p]]

    for (s_idx in seq_along(SENS_SCENARIOS)) {
      s    <- SENS_SCENARIOS[[s_idx]]
      args <- scale_particle_args(ref, s$r0_factor, s$hcw_factor)

      for (arm_idx in seq_len(nrow(ARMS))) {
        arm <- ARMS[arm_idx, ]

        # OBV / PEP settings (identical to 01_analysis_figure2.R).
        args$obv_pep_enabled          <- TRUE
        args$obv_pep_coverage         <- COVERAGE_FNS[[arm$coverage]]
        args$obv_pep_adherence        <- 1.0
        args$obv_pep_dpc              <- 0
        args$obv_pep_efficacy         <- arm$efficacy
        args$obv_pep_target_class     <- "HCW"
        args$obv_pep_target_locations <- c("hospital", "community", "funeral")

        for (r in seq_len(N_REPS)) {
          # Deterministic, collision-free seed across (scaling, arm, particle, rep).
          seed <- SEED_BASE +
            (s_idx   - 1L) * (nrow(ARMS) * N_PARTICLES * N_REPS) +
            (arm_idx - 1L) * (N_PARTICLES * N_REPS) +
            (p       - 1L) * N_REPS +
            (r       - 1L)

          # Retry until the epidemic takes off (matches the figure pipeline).
          out          <- NULL
          retry        <- 0L
          current_seed <- seed
          repeat {
            args$seed <- current_seed
            out       <- do.call(fiber::branching_process_main, args)
            tdf_all   <- out$tdf[!is.na(out$tdf$time_infection_absolute), ]
            n_deaths  <- sum(!is.na(tdf_all$outcome) & tdf_all$outcome)
            if (n_deaths >= TAKEOFF_DEATH_THRESHOLD) break
            retry        <- retry + 1L
            current_seed <- current_seed + 1L
            if (retry >= MAX_RETRIES) break
          }

          m <- summarise_one_run(out)
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
            took_off           = n_deaths >= TAKEOFF_DEATH_THRESHOLD,
            n_hcw_deaths       = m$n_hcw_deaths,                     # WITH PEP
            counterfactual_hcw = m$counterfactual_hcw,              # no-PEP baseline
            prevented_hcw      = m$prevented_hcw,
            hcw_days_lost      = m$hcw_days_lost,
            baseline_days_lost = m$baseline_days_lost,
            total_deaths       = m$total_deaths,
            obv_prevented_rows = m$n_prevented_rows,                # diagnostic
            obv_num_treated    = m$obv_num_treated,                 # diagnostic
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
    base_particle_args = base_particle_args, scale_particle_args = scale_particle_args,
    summarise_one_run = summarise_one_run, run_particle_block = run_particle_block,
    ARMS = ARMS, COVERAGE_FNS = COVERAGE_FNS, SENS_SCENARIOS = SENS_SCENARIOS,
    N_REPS = N_REPS, N_PARTICLES = N_PARTICLES, SEED_BASE = SEED_BASE,
    TAKEOFF_DEATH_THRESHOLD = TAKEOFF_DEATH_THRESHOLD, MAX_RETRIES = MAX_RETRIES
  ),
  future.packages = c("fiber"),
  future.seed     = TRUE
)

plan(sequential)
run_df <- do.call(rbind, results)
elapsed <- (proc.time() - t_start)["elapsed"]
message(sprintf("Done %d simulations in %.1f minutes.", nrow(run_df), elapsed / 60))

# =============================================================================
# Save run-level results first (heavy; gitignored), so the per-run diagnostics
# are on disk even if the prevented_completed health check below halts.
# =============================================================================
run_csv <- file.path(OUT_DIR_HEAVY, "SI_sensitivity_DRC_run_summary.csv")
write.csv(run_df, run_csv, row.names = FALSE)
message("Saved run-level results: ", run_csv)

# =============================================================================
# Sanity + prevented_completed health check.
# The no-PEP baseline is reconstructed from out$prevented_completed, so that
# channel MUST be populated. Report, per analysis, the mean prevented_completed
# row count and fiber's own obv_pep_num_treated$prevented, and STOP loudly if it
# is empty everywhere (which would mean OBV prevention was not recorded and the
# baseline has collapsed onto the with-PEP tdf, giving averted = 0).
# =============================================================================
sanity <- run_df %>%
  group_by(analysis) %>%
  summarise(
    mean_baseline_hcw_deaths = mean(counterfactual_hcw),
    mean_pep_hcw_deaths      = mean(n_hcw_deaths),
    mean_hcw_deaths_averted  = mean(prevented_hcw),
    mean_prevented_rows      = mean(obv_prevented_rows),
    mean_obv_num_treated     = mean(obv_num_treated, na.rm = TRUE),
    frac_runs_empty_prev     = mean(obv_prevented_rows == 0),
    .groups = "drop"
  )
message("Per-analysis sanity (reconstruction from prevented_completed):")
print(as.data.frame(sanity), row.names = FALSE)

frac_empty_overall <- mean(run_df$obv_prevented_rows == 0)
message(sprintf("prevented_completed empty in %.1f%% of runs (mean rows = %.1f, mean obv_num_treated = %.1f).",
                100 * frac_empty_overall, mean(run_df$obv_prevented_rows),
                mean(run_df$obv_num_treated, na.rm = TRUE)))
if (all(run_df$obv_prevented_rows == 0)) {
  stop("out$prevented_completed is EMPTY in every run: OBV prevention was not ",
       "recorded, so the no-PEP baseline collapses onto the with-PEP tdf and ",
       "HCW deaths averted = 0. Run-level diagnostics were saved to\n  ", run_csv,
       "\nInvestigate the OBV/PEP arguments or the fiber build BEFORE trusting ",
       "these results (cf. fiber branching_process_main.R prevented_completed ",
       "population conditions).", call. = FALSE)
}
if (max(sanity$mean_hcw_deaths_averted, na.rm = TRUE) <= 0)
  warning("HCW deaths averted is <= 0 everywhere -- check the OBV/PEP settings.")

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
