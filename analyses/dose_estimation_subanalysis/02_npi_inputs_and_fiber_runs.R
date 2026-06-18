# ============================================================================
# 02_npi_inputs_and_fiber_runs.R
# ----------------------------------------------------------------------------
# WHAT THIS SCRIPT DOES
#   1. Reads the dose Q curve fitted in 01_fit_dose_q_curve.R.
#   2. Turns it into TIME-VARYING NPI inputs: the user defines a min and a max
#      value for each NPI parameter the Q curve drives (prob_hosp, delay to
#      hosp, ETU proportion, safe-funeral proportion, PPE coverage, ...), and
#      each parameter is linearly mapped along the Q curve between those two
#      values. The resulting scenario matrix is saved to outputs/.
#   3. Uses those NPI curves as FIXED inputs to fiber, run in parallel for
#      N_STOCH stochastic replicates across a grid of baseline R0 values
#      (1.3 -> 1.6), seeded with 5 infections, re-running any replicate that
#      fails to reach TAKEOFF_N (= 100) infections.
#   4. For each R0, summarises FROM THE INDIVIDUAL RAW RUNS:
#        * the median cumulative number of cases at a set of timepoints
#          (day 10, 20, 30, ...), and
#        * the median time to reach a set of cumulative case amounts.
#
# HOW Q MAPS TO EACH NPI VALUE
#   For an NPI parameter with worst-response value q0 (at Q = 0) and best-
#   response value q1 (at Q = 1):   value(t) = q0 + (q1 - q0) * Q(t).
#   Direction is encoded purely by whether q0 < q1 (improves upward, e.g. PPE
#   coverage) or q0 > q1 (improves downward, e.g. delay to hospitalisation).
#   "Safe-funeral proportion" is mapped this way and then converted to the
#   model's unsafe-funeral probability as (1 - safe-funeral proportion).
#
# Inputs : outputs/dose_q_curve.rds          (from 01)
# Output : outputs/dose_npi_scenario_matrix.csv / .rds   (the NPI inputs)
#          outputs/dose_npi_timevarying_long.csv         (tidy, for plotting)
#          outputs/dose_r0_grid_per_run.rds              (per-replicate metrics)
#          outputs/dose_r0_grid_cumulative_cases.csv     (median cum cases)
#          outputs/dose_r0_grid_time_to_amounts.csv      (median time-to-amount)
#          outputs/dose_r0_grid_results.rds              (everything bundled)
# ============================================================================

suppressPackageStartupMessages({
  library(here)
  library(future)
  library(future.apply)
  library(fiber)
  library(dplyr)
  library(ggplot2)
})

source(here("analyses", "dose_estimation_subanalysis", "helpers.R"))
source(here("functions", "setup_model_parameters.R"))
source(here("functions", "calculate_model_approx_r0.R"))

set.seed(123)

# ----------------------------------------------------------------------------
# 1. Configuration  <-- the knobs to change
# ----------------------------------------------------------------------------

# --- NPI parameters the Q curve drives. q0 = value at Q = 0 (worst response),
#     q1 = value at Q = 1 (best response). Edit these min/max values freely.
NPI_SPECS <- list(
  prob_hosp         = list(q0 = 0.00, q1 = 0.80),   # P(hospitalised | symptomatic)
  delay_hosp        = list(q0 = 6.00, q1 = 1.50),   # onset->hosp delay factor (days); improves DOWN
  prop_etu          = list(q0 = 0.00, q1 = 0.90),   # proportion of hospitalised cases in an ETU
  safe_funeral_prop = list(q0 = 0.00, q1 = 0.90),   # proportion of funerals that are safe
  ppe_coverage      = list(q0 = 0.00, q1 = 0.90)    # PPE coverage lever (-> ppe_coverage_hcw)
)
# Unsafe-funeral probability for ETU deaths (kept separate; ETU deaths are
# managed safely, so 0 by default).
UNSAFE_FUNERAL_ETU <- 0.0

# --- Fixed efficacy scalars (NOT time-varying; override DEFAULT_SCALAR_INPUTS).
SCALAR_OVERRIDES <- list(
  etu_efficacy                         = 0.90,
  general_hospital_quarantine_efficacy = 0.30,
  ppe_efficacy                         = 0.70,
  safe_funeral_efficacy                = 0.80
)

# --- Simulation grid + controls.
R0_GRID          <- seq(1.30, 1.60, by = 0.05)   # baseline (t=0) R0 grid
FUNERAL_FRAC     <- 0.25                          # share of t=0 transmission via funerals
SEEDING_CASES    <- 5L                            # initial seeding infections
N_STOCH          <- 100L                          # stochastic replicates per R0
TAKEOFF_N        <- 100L                          # re-run any outbreak < this many infections
MAX_RETRIES      <- 50L                           # cap on re-runs per replicate
CHECK_FINAL_SIZE <- 50000L                        # stop a run once this many cases exist

# --- Summary grids (computed per run, then median across runs).
TIMEPOINTS <- c(seq(10L, 360L, by = 10L), 365L)   # days at which to read cumulative cases
AMOUNTS    <- c(10L, 25L, 50L, 100L, 250L, 500L,  # cumulative case amounts to time
                1000L, 2500L, 5000L, 10000L)

# --- Scenario identity + horizon over which the NPI matrix is defined.
SCENARIO_ID     <- "dose_npi"
SCENARIO_LABEL  <- "Dose-driven NPIs"
MATRIX_HORIZON  <- max(730L, max(TIMEPOINTS))     # days; Q held flat past its grid

# --- Epidemic start date. The epidemic (relative day 0) is seeded on this
#     calendar date, which MUST be on or before the Q curve's first date (from
#     01). If it is earlier, the gap before the dose roll-out is filled with
#     Q = 0, so the NPIs sit at their worst-response (q0) values and the epidemic
#     runs unmitigated until the dose curve begins. Set to NA to start exactly at
#     the Q curve's first date (no prepended zeros).
EPIDEMIC_START_DATE <- as.Date("2026-05-01")

# --- Parallel + RNG.
N_WORKERS <- min(future::availableCores(), 50L)
SEED_BASE <- 20260617L

# Fail fast if the installed fiber predates the time-varying NPI interface.
check_model_function_version()

if (CHECK_FINAL_SIZE < max(AMOUNTS)) {
  warning("CHECK_FINAL_SIZE (", CHECK_FINAL_SIZE, ") < max(AMOUNTS) (",
          max(AMOUNTS), "): time-to for the largest amounts will be NA.",
          call. = FALSE)
}

# ----------------------------------------------------------------------------
# 2. Read the Q curve and map it onto the NPI-matrix day grid
# ----------------------------------------------------------------------------
# The Q curve carries calendar dates (the `date` column from 01). The epidemic
# is seeded at EPIDEMIC_START_DATE, which may be EARLIER than the Q curve's first
# date: the days in between (before the dose roll-out) get Q = 0 prepended, so
# the NPIs sit at their worst-response values until the curve begins.
q_curve <- readRDS(file.path(DIR_OUT, "dose_q_curve.rds"))
if (is.list(q_curve) && !is.data.frame(q_curve) && !is.null(q_curve$q_curve)) {
  q_curve <- q_curve$q_curve   # tolerate being handed the full fit list
}
if (is.null(q_curve$date)) {
  stop("The Q curve has no `date` column; re-run 01_fit_dose_q_curve.R so it ",
       "stores calendar dates.", call. = FALSE)
}
q_curve$date <- as.Date(q_curve$date)
q_first_date <- min(q_curve$date)

# Resolve + validate the epidemic start date (day 0 of the simulation).
epi_start <- if (is.null(EPIDEMIC_START_DATE) || all(is.na(EPIDEMIC_START_DATE)))
  q_first_date else as.Date(EPIDEMIC_START_DATE)
if (epi_start > q_first_date) {
  stop("EPIDEMIC_START_DATE (", epi_start, ") must be on or before the Q curve's ",
       "first date (", q_first_date, ").", call. = FALSE)
}
offset_days <- as.integer(q_first_date - epi_start)   # days of Q = 0 prepended
message(sprintf("Epidemic starts %s; Q curve starts %s -> %d day(s) of Q = 0 prepended.",
                as.character(epi_start), as.character(q_first_date), offset_days))

# Map the Q curve onto the simulation's relative-day axis (day 0 = epi_start):
# each Q point's sim day = its calendar date - epi_start. rule = 2 holds Q flat
# beyond the curve's last day; days before the roll-out are set to exactly 0.
matrix_days <- 0:MATRIX_HORIZON
q_sim_day   <- as.integer(q_curve$date - epi_start)
q_on_grid   <- approx(q_sim_day, q_curve$q_mean, xout = matrix_days, rule = 2)$y
q_on_grid[matrix_days < offset_days] <- 0
q_on_grid   <- clip01(q_on_grid)

# ----------------------------------------------------------------------------
# 3. Build + save the time-varying NPI inputs (the scenario matrix)
# ----------------------------------------------------------------------------
lin <- function(spec, q) spec$q0 + (spec$q1 - spec$q0) * q

safe_funeral_prop <- clip01(lin(NPI_SPECS$safe_funeral_prop, q_on_grid))

scenario_matrix_df <- data.frame(
  scenario                 = SCENARIO_ID,
  scenario_label           = SCENARIO_LABEL,
  relative_day             = matrix_days,
  prob_hosp                = clip01(lin(NPI_SPECS$prob_hosp, q_on_grid)),
  delay_hosp               = pmax(lin(NPI_SPECS$delay_hosp, q_on_grid), 0.01),
  prob_unsafe_funeral_comm = clip01(1 - safe_funeral_prop),
  prob_unsafe_funeral_hosp = clip01(1 - safe_funeral_prop),
  prob_unsafe_funeral_etu  = clip01(rep(UNSAFE_FUNERAL_ETU, length(matrix_days))),
  prop_etu                 = clip01(lin(NPI_SPECS$prop_etu, q_on_grid)),
  ipc_helper               = clip01(lin(NPI_SPECS$ppe_coverage, q_on_grid)),
  q_value                  = q_on_grid,
  stringsAsFactors = FALSE
)

matrix_csv <- file.path(DIR_OUT, "dose_npi_scenario_matrix.csv")
write.csv(scenario_matrix_df, matrix_csv, row.names = FALSE)
saveRDS(scenario_matrix_df, file.path(DIR_OUT, "dose_npi_scenario_matrix.rds"))

# Tidy long version (handy for plotting the inputs).
tv_long <- scenario_matrix_df %>%
  select(relative_day, q_value, prob_hosp, delay_hosp,
         prob_unsafe_funeral_comm, prob_unsafe_funeral_hosp,
         prop_etu, ipc_helper) %>%
  tidyr::pivot_longer(-relative_day, names_to = "input", values_to = "value")
write.csv(tv_long, file.path(DIR_OUT, "dose_npi_timevarying_long.csv"), row.names = FALSE)

message("Wrote NPI scenario matrix (", nrow(scenario_matrix_df), " daily rows) to outputs/.")

# ----------------------------------------------------------------------------
# 4. Assemble the fiber model args per R0
# ----------------------------------------------------------------------------
# Read the matrix back through the validated reader, build the base + time-
# varying args once, then invert each grid R0 into offspring means at t = 0.
scenario_matrix <- read_scenario_matrix(matrix_csv)

mp <- make_model_parameters(
  scenario_id     = SCENARIO_ID,
  scenario_matrix = scenario_matrix,
  overrides       = modifyList(
    list(seeding_cases = SEEDING_CASES, check_final_size = CHECK_FINAL_SIZE),
    SCALAR_OVERRIDES
  )
)

# Efficacy-independent R0 invariants at t = 0 (Q(0) = worst-response state).
inv   <- compute_R0_invariants(mp$args, n = 50000L, seed = 42L)
D     <- D_from_invariants(inv, mp$args$etu_efficacy,
                           mp$args$general_hospital_quarantine_efficacy)
F_fun <- F_from_invariants(inv, mp$args$safe_funeral_efficacy)

# One ready-to-run args list per R0 (offspring means solved from the target R0).
args_by_r0 <- lapply(R0_GRID, function(R0) {
  means <- solve_offspring_means(R0, FUNERAL_FRAC, D, F_fun)
  a <- mp$args
  a$mn_offspring_genPop  <- means$mn_genPop
  a$mn_offspring_funeral <- means$mn_funeral
  a$seeding_cases        <- SEEDING_CASES
  a$check_final_size     <- CHECK_FINAL_SIZE
  a$seed                 <- NULL
  a
})
names(args_by_r0) <- sprintf("R0_%.2f", R0_GRID)

# ----------------------------------------------------------------------------
# 5. Per-replicate runner + summariser (used on the parallel workers)
# ----------------------------------------------------------------------------
# Run ONE replicate, re-running (advancing the seed) until the outbreak reaches
# `takeoff_n` infections or `max_retries` is hit. Returns the realised infection
# times plus bookkeeping.
run_one_takeoff <- function(args, base_seed, takeoff_n, max_retries) {
  seed <- base_seed; retry <- 0L; inf_times <- numeric(0); n_cases <- 0L
  repeat {
    args$seed <- seed
    out <- do.call(fiber::branching_process_main, args)
    it  <- out$tdf$time_infection_absolute
    it  <- it[!is.na(it)]
    n_cases <- length(it)
    if (n_cases >= takeoff_n || retry >= max_retries) { inf_times <- it; break }
    retry <- retry + 1L; seed <- seed + 1L
  }
  list(inf_times = inf_times, n_cases = n_cases,
       took_off = n_cases >= takeoff_n, final_seed = seed, retries = retry)
}

# Per-run summaries computed from the raw infection times:
#   cum_at[t]  = cumulative cases by day t
#   time_to[a] = day on which the a-th infection occurred (time to `a` cases),
#                or NA if the run never reached `a` cases.
summarise_run <- function(inf_times, timepoints, amounts) {
  st <- sort(inf_times); n <- length(st)
  list(
    cum_at  = vapply(timepoints, function(t) sum(st <= t), numeric(1)),
    time_to = vapply(amounts,    function(a) if (n >= a) st[a] else NA_real_, numeric(1)),
    n_cases = n
  )
}

# ----------------------------------------------------------------------------
# 6. Run the grid in parallel, one R0 at a time, with progress + ETA
# ----------------------------------------------------------------------------
# Reps for a given R0 run in parallel across the workers; we loop over the R0
# grid in the MAIN process so we can report how far through we are (and a rough
# ETA) after each R0 finishes. A "running ..." line is printed before each batch
# so progress is visible even while the first (cold) batch is in flight.
n_r0  <- length(R0_GRID)
total <- n_r0 * N_STOCH
message(sprintf("Running %d R0 x %d reps = %d replicates on %d workers...",
                n_r0, N_STOCH, total, N_WORKERS))

plan(multisession, workers = N_WORKERS)
t_start <- proc.time()
results <- vector("list", n_r0)

for (i in seq_along(R0_GRID)) {
  args_i <- args_by_r0[[i]]
  r0_i   <- R0_GRID[i]
  message(sprintf("  [%d/%d] R0 = %.2f: running %d reps ...", i, n_r0, r0_i, N_STOCH))

  results[[i]] <- future_lapply(seq_len(N_STOCH), function(rep_id) {
    base_seed <- SEED_BASE +
      (i      - 1L) * N_STOCH * (MAX_RETRIES + 1L) +
      (rep_id - 1L) *           (MAX_RETRIES + 1L)
    run <- run_one_takeoff(args_i, base_seed, TAKEOFF_N, MAX_RETRIES)
    s   <- summarise_run(run$inf_times, TIMEPOINTS, AMOUNTS)
    list(r0 = r0_i, r0_idx = i, rep_id = rep_id,
         n_cases = run$n_cases, took_off = run$took_off, retries = run$retries,
         cum_at = s$cum_at, time_to = s$time_to)
  },
  future.globals = list(
    args_i = args_i, i = i, r0_i = r0_i, N_STOCH = N_STOCH,
    MAX_RETRIES = MAX_RETRIES, SEED_BASE = SEED_BASE, TAKEOFF_N = TAKEOFF_N,
    TIMEPOINTS = TIMEPOINTS, AMOUNTS = AMOUNTS,
    run_one_takeoff = run_one_takeoff, summarise_run = summarise_run
  ),
  future.packages = "fiber", future.seed = TRUE)

  # Progress + rough ETA. Later R0s tend to cost a little more (bigger
  # outbreaks), so the early ETA is somewhat optimistic.
  n_takeoff <- sum(vapply(results[[i]], function(r) isTRUE(r$took_off), logical(1)))
  elapsed_s <- (proc.time() - t_start)[["elapsed"]]
  eta_s     <- elapsed_s / i * (n_r0 - i)
  message(sprintf(
    "  [%d/%d] R0 = %.2f done | %d/%d reps reached takeoff | elapsed %.1f min | rough ETA %.1f min",
    i, n_r0, r0_i, n_takeoff, N_STOCH, elapsed_s / 60, eta_s / 60))
}

plan(sequential)

# Flatten the per-R0 lists into one flat list of per-replicate results.
results <- do.call(c, results)

elapsed <- proc.time() - t_start
message(sprintf("Done: %d replicates in %.1f min.", length(results), elapsed["elapsed"] / 60))

saveRDS(results, file.path(DIR_OUT, "dose_r0_grid_per_run.rds"))

# ----------------------------------------------------------------------------
# 7. Summarise across replicates (medians from the individual raw runs)
# ----------------------------------------------------------------------------
safe_median <- function(x) if (all(is.na(x))) NA_real_ else stats::median(x, na.rm = TRUE)
safe_q <- function(x, p) if (all(is.na(x))) NA_real_ else
  as.numeric(stats::quantile(x, p, na.rm = TRUE, names = FALSE))

took   <- Filter(function(r) isTRUE(r$took_off), results)
n_fail <- length(results) - length(took)
if (n_fail > 0L)
  message(sprintf("Note: %d/%d replicates failed to reach %d infections after %d retries (excluded from summaries).",
                  n_fail, length(results), TAKEOFF_N, MAX_RETRIES))
if (length(took) == 0L)
  stop("No replicate reached takeoff (", TAKEOFF_N, " infections). ",
       "Check R0_GRID, the NPI settings, MAX_RETRIES or CHECK_FINAL_SIZE.",
       call. = FALSE)
by_r0 <- split(took, vapply(took, function(r) r$r0_idx, integer(1)))

# (a) Median cumulative cases at each timepoint.
cumulative_cases <- do.call(rbind, lapply(seq_along(R0_GRID), function(i) {
  rs <- by_r0[[as.character(i)]]
  if (is.null(rs)) return(NULL)
  M <- do.call(rbind, lapply(rs, function(r) r$cum_at))   # n_runs x n_timepoints
  data.frame(
    r0               = R0_GRID[i],
    timepoint_day    = TIMEPOINTS,
    n_runs           = nrow(M),
    median_cum_cases = apply(M, 2, stats::median),
    q25_cum_cases    = apply(M, 2, safe_q, 0.25),
    q75_cum_cases    = apply(M, 2, safe_q, 0.75),
    stringsAsFactors = FALSE
  )
}))

# (b) Median time to reach each cumulative case amount.
time_to_amounts <- do.call(rbind, lapply(seq_along(R0_GRID), function(i) {
  rs <- by_r0[[as.character(i)]]
  if (is.null(rs)) return(NULL)
  M <- do.call(rbind, lapply(rs, function(r) r$time_to))  # n_runs x n_amounts
  data.frame(
    r0              = R0_GRID[i],
    case_amount     = AMOUNTS,
    n_runs_reaching = apply(M, 2, function(x) sum(!is.na(x))),
    median_time     = apply(M, 2, safe_median),
    q25_time        = apply(M, 2, safe_q, 0.25),
    q75_time        = apply(M, 2, safe_q, 0.75),
    stringsAsFactors = FALSE
  )
}))

write.csv(cumulative_cases, file.path(DIR_OUT, "dose_r0_grid_cumulative_cases.csv"), row.names = FALSE)
write.csv(time_to_amounts,  file.path(DIR_OUT, "dose_r0_grid_time_to_amounts.csv"),  row.names = FALSE)

saveRDS(list(
  cumulative_cases = cumulative_cases,
  time_to_amounts  = time_to_amounts,
  per_run          = results,
  config = list(
    npi_specs = NPI_SPECS, unsafe_funeral_etu = UNSAFE_FUNERAL_ETU,
    scalar_overrides = SCALAR_OVERRIDES, r0_grid = R0_GRID,
    funeral_frac = FUNERAL_FRAC, seeding_cases = SEEDING_CASES,
    n_stoch = N_STOCH, takeoff_n = TAKEOFF_N, max_retries = MAX_RETRIES,
    check_final_size = CHECK_FINAL_SIZE, timepoints = TIMEPOINTS, amounts = AMOUNTS,
    epidemic_start_date = epi_start, q_first_date = q_first_date,
    offset_days = offset_days
  )
), file.path(DIR_OUT, "dose_r0_grid_results.rds"))

message("\n02_npi_inputs_and_fiber_runs.R complete. Results in outputs/.")

# ----------------------------------------------------------------------------
# 8. Quick diagnostic plots (display + saved)
# ----------------------------------------------------------------------------
# (i) The time-varying NPI inputs.
p_inputs <- ggplot(tv_long, aes(relative_day, value)) +
  geom_line(colour = "#1f77b4", linewidth = 0.8) +
  facet_wrap(~ input, scales = "free_y") +
  labs(title = "Time-varying NPI inputs driven by the dose Q curve",
       x = "Relative day", y = "Input value") +
  theme_bw(base_size = 10)
ggsave(file.path(DIR_OUT, "dose_npi_timevarying.png"), p_inputs, width = 9, height = 6, dpi = 150)
print(p_inputs)

# (ii) Median cumulative cases over time, by R0.
p_cum <- ggplot(cumulative_cases,
                aes(timepoint_day, median_cum_cases, colour = factor(r0))) +
  geom_ribbon(aes(ymin = q25_cum_cases, ymax = q75_cum_cases, fill = factor(r0)),
              colour = NA, alpha = 0.12) +
  geom_line(linewidth = 0.8) +
  labs(title = "Median cumulative cases by baseline R0",
       subtitle = "Lines = median across replicates; bands = 25-75%",
       x = "Day", y = "Cumulative cases", colour = "R0", fill = "R0") +
  theme_bw(base_size = 11)
ggsave(file.path(DIR_OUT, "dose_r0_grid_cumulative_cases.png"), p_cum, width = 9, height = 5.5, dpi = 150)
print(p_cum)
