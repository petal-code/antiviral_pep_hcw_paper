# ============================================================================
# time_final_size_runs.R
# ----------------------------------------------------------------------------
# Quick benchmark: how long does ONE fiber run with check_final_size = 75000
# take, exactly as configured in 05_projection_6_qcurves.R?
#
# It reuses 05's setup (config -> NPI matrices -> R0 invariants -> args_grid),
# then times one *established* (taken-off) run per R0 value and reports the
# per-run elapsed time plus the average. Bump REPS_PER_R0 to average more.
#
#   Rscript analyses/dose_estimation_subanalysis/time_final_size_runs.R
#
# Nothing is written to disk; results are printed to the console.
# ============================================================================

suppressPackageStartupMessages({
  library(here)
  library(fiber)
})

source(here("analyses", "dose_estimation_subanalysis", "helpers.R"))
source(here("functions", "setup_model_parameters.R"))
source(here("functions", "calculate_model_approx_r0.R"))

set.seed(123)

# ============================================================================
# Timing knobs
# ============================================================================
# Which NPI scenario(s) to time. Default is "conflict" — the lowest-control
# scenario, so its epidemics are the most likely to hit the 75000 cap and give
# the *heaviest* (slowest) single-run time. Set to scen_names (loaded below) to
# time all six scenarios, or to e.g. "logistic" for the current-trajectory case.
TIMING_SCENARIOS <- "conflict"

# Runs to time per (R0, scenario). The user asked for ~1 per R0; raise this for
# a smoother average.
REPS_PER_R0 <- 1L

# ============================================================================
# Config (copied verbatim from 05_projection_6_qcurves.R, section 1)
# ============================================================================
R0_GRID <- c(1.45, 1.50, 1.55, 1.60, 1.65)

SCALAR_OVERRIDES <- list(
  etu_efficacy                         = 0.84,
  general_hospital_quarantine_efficacy = 0.3,
  ppe_efficacy                         = 0.84,
  safe_funeral_efficacy                = 0.88
)

NPI_SPECS <- list(
  prob_hosp                = list(q0 = 0.00, q1 = 0.80),
  delay_hosp               = list(q0 = 6.00, q1 = 1.50),
  prop_etu                 = list(q0 = 0.00, q1 = 0.90),
  safe_funeral_prop        = list(q0 = 0.00, q1 = 0.90),
  unsafe_funeral_prop_hosp = list(q0 = 0.05, q1 = 0.01),
  ppe_coverage             = list(q0 = 0.00, q1 = 0.90)
)
UNSAFE_FUNERAL_ETU <- 0.0

FUNERAL_FRAC        <- 0.25
SEEDING_CASES       <- 5L
CHECK_FINAL_SIZE    <- 75000L          # <-- the value being benchmarked

# Takeoff criterion (so the timed run is an established epidemic, not an early
# stochastic extinction) — identical to 05.
TAKEOFF_N             <- 250L
TAKEOFF_DEADLINE_DATE <- as.Date("2026-06-15")
MAX_RETRIES           <- 50L

EPIDEMIC_START_DATE <- as.Date("2026-02-17")
MATRIX_HORIZON      <- 730L
SEED_BASE           <- 20260619L

epi_start <- as.Date(EPIDEMIC_START_DATE)

# ============================================================================
# Load Q-curve scenarios (05 section 3)
# ============================================================================
scen_path <- file.path(DIR_OUT, "dose_q_curve_extrapolation_scenarios.rds")
if (!file.exists(scen_path))
  stop("Run 01_fit_dose_q_curve.R first to generate: ", scen_path, call. = FALSE)

all_scen      <- readRDS(scen_path)
all_scen$date <- as.Date(all_scen$date)
scen_names    <- if (is.factor(all_scen$scenario)) levels(all_scen$scenario) else
                   as.character(unique(all_scen$scenario))

q_first_date <- min(all_scen$date)
offset_days  <- as.integer(q_first_date - epi_start)

# ============================================================================
# Build NPI scenario matrices (05 section 4, in-memory)
# ============================================================================
lin         <- function(spec, q) spec$q0 + (spec$q1 - spec$q0) * q
matrix_days <- 0:MATRIX_HORIZON

build_matrix_df <- function(sn) {
  sel    <- all_scen[as.character(all_scen$scenario) == sn, , drop = FALSE]
  q_sim  <- as.integer(sel$date - epi_start)
  q_vals <- approx(q_sim, sel$q, xout = matrix_days, rule = 2)$y
  q_vals[matrix_days < offset_days] <- 0
  q_vals <- clip01(q_vals)

  safe_fp   <- clip01(lin(NPI_SPECS$safe_funeral_prop,        q_vals))
  unsafe_fh <- clip01(lin(NPI_SPECS$unsafe_funeral_prop_hosp, q_vals))

  data.frame(
    scenario                 = sn,
    scenario_label           = sn,
    relative_day             = matrix_days,
    prob_hosp                = clip01(lin(NPI_SPECS$prob_hosp,    q_vals)),
    delay_hosp               = pmax(lin(NPI_SPECS$delay_hosp,     q_vals), 0.01),
    prob_unsafe_funeral_comm = clip01(1 - safe_fp),
    prob_unsafe_funeral_hosp = unsafe_fh,
    prob_unsafe_funeral_etu  = clip01(rep(UNSAFE_FUNERAL_ETU, length(matrix_days))),
    prop_etu                 = clip01(lin(NPI_SPECS$prop_etu,     q_vals)),
    ipc_helper               = clip01(lin(NPI_SPECS$ppe_coverage, q_vals)),
    q_value                  = q_vals,
    stringsAsFactors = FALSE
  )
}
scen_mats <- setNames(lapply(scen_names, build_matrix_df), scen_names)

# ============================================================================
# Assemble FIBER args per (R0, scenario) (05 section 6)
# ============================================================================
overrides_base <- modifyList(
  list(seeding_cases = SEEDING_CASES, check_final_size = CHECK_FINAL_SIZE),
  SCALAR_OVERRIDES
)

# R0 invariants (expensive Monte Carlo step, done once).
mp_ref <- make_model_parameters(scen_names[1], scen_mats[[scen_names[1]]],
                                overrides = overrides_base)
inv    <- compute_R0_invariants(mp_ref$args, n = 50000L, seed = 42L)
D      <- D_from_invariants(inv, mp_ref$args$etu_efficacy,
                            mp_ref$args$general_hospital_quarantine_efficacy)
F_fun  <- F_from_invariants(inv, mp_ref$args$safe_funeral_efficacy)

args_grid <- list()
for (r0 in R0_GRID) {
  means <- solve_offspring_means(r0, FUNERAL_FRAC, D, F_fun)
  for (sn in scen_names) {
    key <- paste(sprintf("R0_%.2f", r0), sn, sep = "__")
    mp  <- make_model_parameters(sn, scen_mats[[sn]], overrides = overrides_base)
    a   <- mp$args
    a$mn_offspring_genPop  <- means$mn_genPop
    a$mn_offspring_funeral <- means$mn_funeral
    a$seed                 <- NULL
    args_grid[[key]] <- a
  }
}

# ============================================================================
# Timing
# ============================================================================
takeoff_day <- as.integer(TAKEOFF_DEADLINE_DATE - epi_start)

# Time ONE established run: advance the seed until the epidemic takes off (as 05
# does), but only record the elapsed time of that final, taken-off run. Early
# extinctions are reported via `retries` and are NOT counted in `elapsed_s`.
time_one_run <- function(args, base_seed) {
  seed <- base_seed; retry <- 0L
  repeat {
    args$seed <- seed
    tt  <- system.time(out <- do.call(fiber::branching_process_main, args))
    it  <- out$tdf$time_infection_absolute
    it  <- it[!is.na(it)]
    took_off <- sum(it <= takeoff_day) >= TAKEOFF_N
    if (took_off || retry >= MAX_RETRIES)
      return(list(elapsed = unname(tt["elapsed"]), n_cases = length(it),
                  hit_cap = length(it) >= CHECK_FINAL_SIZE,
                  took_off = took_off, retries = retry))
    retry <- retry + 1L; seed <- seed + 1L
  }
}

if (identical(TIMING_SCENARIOS, "all")) TIMING_SCENARIOS <- scen_names
stopifnot(all(TIMING_SCENARIOS %in% scen_names))

message(sprintf("Timing %d run(s) per R0 x scenario | check_final_size = %d",
                REPS_PER_R0, CHECK_FINAL_SIZE))
message(sprintf("Scenarios: %s | R0: %s\n",
                paste(TIMING_SCENARIOS, collapse = ", "),
                paste(sprintf("%.2f", R0_GRID), collapse = ", ")))

rows <- list()
for (sn in TIMING_SCENARIOS) {
  for (r0 in R0_GRID) {
    args <- args_grid[[paste(sprintf("R0_%.2f", r0), sn, sep = "__")]]
    for (rep in seq_len(REPS_PER_R0)) {
      res <- time_one_run(args, SEED_BASE + rep)
      rows[[length(rows) + 1L]] <- data.frame(
        scenario = sn, r0 = r0, rep = rep, elapsed_s = res$elapsed,
        n_cases = res$n_cases, hit_cap = res$hit_cap, retries = res$retries)
      message(sprintf(
        "  %-9s R0=%.2f rep %d: %6.2f s | final size %6d%s | %d retr%s to take off",
        sn, r0, rep, res$elapsed, res$n_cases,
        if (res$hit_cap) " (hit cap)" else "",
        res$retries, if (res$retries == 1L) "y" else "ies"))
    }
  }
}
timing <- do.call(rbind, rows)

cat("\n--- Per-run timing ---\n")
print(timing, row.names = FALSE)

cat("\n--- Mean elapsed per run, by R0 ---\n")
print(aggregate(elapsed_s ~ r0, timing, mean), row.names = FALSE)

cat(sprintf(
  "\nOverall mean: %.2f s/run over %d run(s)  [scenario(s): %s, check_final_size = %d]\n",
  mean(timing$elapsed_s), nrow(timing),
  paste(TIMING_SCENARIOS, collapse = ", "), CHECK_FINAL_SIZE))
