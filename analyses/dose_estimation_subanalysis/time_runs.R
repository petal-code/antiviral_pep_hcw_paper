# ============================================================================
# time_runs.R -- quick wall-clock timing of single FIBER replicates at the 05
# configuration (CHECK_FINAL_SIZE = 75000), one run per R0 value.
#
# Mirrors 05_projection_6_qcurves.R's model setup exactly (NPI matrix -> R0
# invariants -> offspring means -> run_one_takeoff) so the per-run time is
# representative of the real run. Times run_one_takeoff() (the unit of work in
# 05, including any takeoff retries) for each R0 and prints the average.
#
# Run from the repo root:  Rscript analyses/dose_estimation_subanalysis/time_runs.R
# ============================================================================

suppressPackageStartupMessages({ library(here); library(fiber) })
source(here("analyses", "dose_estimation_subanalysis", "helpers.R"))        # DIR_OUT, clip01
source(here("functions", "setup_model_parameters.R"))                       # read_scenario_matrix, make_model_parameters
source(here("functions", "calculate_model_approx_r0.R"))                    # *_invariants, solve_offspring_means

# ---- config (must match 05) ------------------------------------------------
R0_GRID          <- c(1.45, 1.50, 1.55, 1.60, 1.65)
REPS_PER_R0      <- 1L            # runs to time per R0 (bump to 2-3 for a less noisy mean)
# Which NPI scenario(s) to time. Cost varies A LOT by scenario: well-controlled
# curves stay small (fast); poorly-controlled ones grow to the cap (slow). The
# slow, cap-hitting ones dominate the real run, so "conflict" is a conservative
# default. Set to a full vector e.g. c("logistic","flat","conflict") -- or to
# all six -- for a truer average (more runs = longer test).
TEST_SCENARIOS   <- c("conflict")

CHECK_FINAL_SIZE <- 75000L
SEEDING_CASES    <- 5L
FUNERAL_FRAC     <- 0.25
TAKEOFF_N             <- 250L
TAKEOFF_DEADLINE_DATE <- as.Date("2026-06-15")
MAX_RETRIES           <- 50L
EPIDEMIC_START_DATE   <- as.Date("2026-02-17")
SEED_BASE        <- 20260619L
MATRIX_HORIZON   <- 730L
N_WORKERS_EST    <- 106L          # workers you expect for the full run (for the projection)

SCALAR_OVERRIDES <- list(
  etu_efficacy = 0.84, general_hospital_quarantine_efficacy = 0.3,
  ppe_efficacy = 0.84, safe_funeral_efficacy = 0.88
)
NPI_SPECS <- list(
  prob_hosp = list(q0 = 0.00, q1 = 0.80), delay_hosp = list(q0 = 6.00, q1 = 1.50),
  prop_etu = list(q0 = 0.00, q1 = 0.90), safe_funeral_prop = list(q0 = 0.00, q1 = 0.90),
  unsafe_funeral_prop_hosp = list(q0 = 0.05, q1 = 0.01), ppe_coverage = list(q0 = 0.00, q1 = 0.90)
)
UNSAFE_FUNERAL_ETU <- 0.0

epi_start   <- EPIDEMIC_START_DATE
takeoff_day <- as.integer(TAKEOFF_DEADLINE_DATE - epi_start)
lin         <- function(spec, q) spec$q0 + (spec$q1 - spec$q0) * q

# ---- build a scenario's NPI matrix (verbatim from 05 section 4) -------------
all_scen      <- readRDS(file.path(DIR_OUT, "dose_q_curve_extrapolation_scenarios.rds"))
all_scen$date <- as.Date(all_scen$date)
offset_days   <- as.integer(min(all_scen$date) - epi_start)
matrix_days   <- 0:MATRIX_HORIZON

build_matrix_csv <- function(sn) {
  sel    <- all_scen[as.character(all_scen$scenario) == sn, , drop = FALSE]
  q_vals <- approx(as.integer(sel$date - epi_start), sel$q, xout = matrix_days, rule = 2)$y
  q_vals[matrix_days < offset_days] <- 0
  q_vals <- clip01(q_vals)
  safe_fp <- clip01(lin(NPI_SPECS$safe_funeral_prop, q_vals))
  df <- data.frame(
    scenario = sn, scenario_label = sn, relative_day = matrix_days,
    prob_hosp                = clip01(lin(NPI_SPECS$prob_hosp,  q_vals)),
    delay_hosp               = pmax(lin(NPI_SPECS$delay_hosp,    q_vals), 0.01),
    prob_unsafe_funeral_comm = clip01(1 - safe_fp),
    prob_unsafe_funeral_hosp = clip01(lin(NPI_SPECS$unsafe_funeral_prop_hosp, q_vals)),
    prob_unsafe_funeral_etu  = clip01(rep(UNSAFE_FUNERAL_ETU, length(matrix_days))),
    prop_etu                 = clip01(lin(NPI_SPECS$prop_etu,    q_vals)),
    ipc_helper               = clip01(lin(NPI_SPECS$ppe_coverage, q_vals)),
    q_value                  = q_vals, stringsAsFactors = FALSE)
  pth <- tempfile(fileext = ".csv"); write.csv(df, pth, row.names = FALSE); pth
}

overrides_base <- modifyList(
  list(seeding_cases = SEEDING_CASES, check_final_size = CHECK_FINAL_SIZE), SCALAR_OVERRIDES)

# R0 invariants (efficacy-only; computed once, outside the timed loop) --------
ref_sm <- read_scenario_matrix(build_matrix_csv(TEST_SCENARIOS[1]))
mp_ref <- make_model_parameters(TEST_SCENARIOS[1], ref_sm, overrides = overrides_base)
inv    <- compute_R0_invariants(mp_ref$args, n = 50000L, seed = 42L)
D      <- D_from_invariants(inv, mp_ref$args$etu_efficacy,
                            mp_ref$args$general_hospital_quarantine_efficacy)
F_fun  <- F_from_invariants(inv, mp_ref$args$safe_funeral_efficacy)

# ---- runner (verbatim from 05 section 8) -----------------------------------
run_one_takeoff <- function(args, base_seed, takeoff_n, max_retries, takeoff_day) {
  seed <- base_seed; retry <- 0L; inf_times <- numeric(0); n_cases <- 0L; n_by_deadline <- 0L
  repeat {
    args$seed <- seed
    out <- do.call(fiber::branching_process_main, args)
    it  <- out$tdf$time_infection_absolute; it <- it[!is.na(it)]
    n_cases <- length(it); n_by_deadline <- sum(it <= takeoff_day)
    if (n_by_deadline >= takeoff_n || retry >= max_retries) { inf_times <- it; break }
    retry <- retry + 1L; seed <- seed + 1L
  }
  list(n_cases = n_cases, took_off = n_by_deadline >= takeoff_n, retries = retry)
}

# ---- time one run per (scenario, R0) ---------------------------------------
res <- data.frame()
for (sn in TEST_SCENARIOS) {
  sm <- read_scenario_matrix(build_matrix_csv(sn))
  mp <- make_model_parameters(sn, sm, overrides = overrides_base)
  for (i in seq_along(R0_GRID)) {
    r0    <- R0_GRID[i]
    means <- solve_offspring_means(r0, FUNERAL_FRAC, D, F_fun)
    a <- mp$args
    a$mn_offspring_genPop  <- means$mn_genPop
    a$mn_offspring_funeral <- means$mn_funeral
    a$seed <- NULL
    for (rep in seq_len(REPS_PER_R0)) {
      seed0 <- SEED_BASE + i * 1000L + rep
      t <- system.time(run <- run_one_takeoff(a, seed0, TAKEOFF_N, MAX_RETRIES, takeoff_day))[["elapsed"]]
      res <- rbind(res, data.frame(scenario = sn, r0 = r0, rep = rep,
                                   secs = round(t, 2), n_cases = run$n_cases,
                                   retries = run$retries, hit_cap = run$n_cases >= CHECK_FINAL_SIZE))
      cat(sprintf("%-12s R0 %.2f rep %d: %6.1f s | n_cases = %6d | retries = %d | hit_cap = %s\n",
                  sn, r0, rep, t, run$n_cases, run$retries, run$hit_cap))
    }
  }
}

# ---- summary + projection to the full run ----------------------------------
cat("\n--- per run ---\n"); print(res, row.names = FALSE)
m <- mean(res$secs)
cat(sprintf("\nMean per-run wall time: %.1f s (%.2f min) over %d runs | cap %d | scenario(s): %s\n",
            m, m / 60, nrow(res), CHECK_FINAL_SIZE, paste(TEST_SCENARIOS, collapse = ", ")))

total_sims <- length(R0_GRID) * 6L * 330L         # 5 R0 x 6 scenarios x 330 reps = 9,900
core_h     <- total_sims * m / 3600
cat(sprintf("Full run = %d sims -> %.1f core-hours -> ~%.1f h wall on %d workers.\n",
            total_sims, core_h, core_h / N_WORKERS_EST, N_WORKERS_EST))
cat("NB: per-run cost varies strongly by scenario (cap-hitting = slow). If you\n",
    "   timed only slow/fast scenarios this over/under-states the true average;\n",
    "   set TEST_SCENARIOS <- c(\"linear_to_90\",\"logistic\",\"flat\",\"conflict\", ...) for a mix.\n", sep = "")
