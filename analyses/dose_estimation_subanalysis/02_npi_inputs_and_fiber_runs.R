# ============================================================================
# 02_npi_inputs_and_fiber_runs.R
# ----------------------------------------------------------------------------
# WHAT THIS SCRIPT DOES
#   1. Reads ONE forward-extrapolation scenario of the dose Q curve (chosen via
#      EXTRAP_SCENARIO) from 01_fit_dose_q_curve.R.
#   2. Turns it into TIME-VARYING NPI inputs: the user defines a min and a max
#      value for each NPI parameter the Q curve drives (prob_hosp, delay to
#      hosp, ETU proportion, safe-funeral proportion, PPE coverage, ...), and
#      each parameter is linearly mapped along the Q curve between those two
#      values. The resulting scenario matrix is saved to outputs/.
#   3. BEFORE the long runs, computes + overlays an analytic single-type Rt
#      profile per R0 (Rt_curve_single_type), as a quick sanity check that Rt
#      starts near each target R0 and falls below 1 as the response ramps up.
#   4. Uses those NPI curves as FIXED inputs to fiber, run in parallel for
#      N_STOCH stochastic replicates across a grid of baseline R0 values
#      (1.3 -> 1.6), seeded with 5 infections, re-running any replicate that
#      fails to reach TAKEOFF_N cumulative infections BY the takeoff-deadline
#      date (TAKEOFF_DEADLINE_DATE, relative to EPIDEMIC_START_DATE).
#   5. For each R0, summarises FROM THE INDIVIDUAL RAW RUNS:
#        * the mean (and median) cumulative number of cases at a set of timepoints
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
# Inputs : outputs/dose_q_curve_extrapolation_scenarios.rds  (from 01; pick one
#          via EXTRAP_SCENARIO). Falls back to dose_q_curve.rds (the logistic
#          projection) only if the scenarios file is absent and EXTRAP_SCENARIO
#          is "logistic".
#          data-processed/onsets_daily_incidence.csv  (optional; from
#          analyses/onset_incidence/compute_daily_incidence_from_onsets.R --
#          its observed cumulative-onset curve is overlaid on the trajectory plots).
#          data-processed/insp_sitrep__national_cumulative_confirmed_cases__daily.csv
#          (optional; observed national cumulative confirmed cases, also overlaid).
# Output : outputs/dose_npi_scenario_matrix.csv / .rds   (the NPI inputs)
#          outputs/dose_npi_timevarying_long.csv         (tidy, for plotting)
#          outputs/dose_r0_grid_rt_profiles.csv / .png   (analytic Rt per R0, pre-run)
#          outputs/dose_r0_grid_per_run.rds              (per-replicate metrics)
#          outputs/dose_r0_grid_cumulative_cases.csv     (mean + median cum cases)
#          outputs/dose_r0_grid_cumulative_trajectories.png (per-replicate curves by R0)
#          outputs/dose_r0_grid_daily_incidence.png      (implied daily incidence by R0)
#          outputs/dose_r0_grid_snapshot_cumulative.png  (cum infections at 14 May / 14 Jun)
#          outputs/dose_r0_grid_time_to_amounts.csv      (median time-to-amount)
#          outputs/dose_growth_rates_doubling_times.csv / .png (growth rate + doubling
#                                          time per period: fiber R0s vs data sources)
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
source(here("functions", "calculate_model_approx_rt.R"))

set.seed(123)

# ----------------------------------------------------------------------------
# 1. Configuration  <-- the knobs to change
# ----------------------------------------------------------------------------

# --- Which forward-extrapolation scenario of the dose Q curve to run. One of the
#     names produced by 01_fit_dose_q_curve.R:
#       "linear_to_90" | "logistic" | "flat" | "conflict".
EXTRAP_SCENARIO <- "logistic"

# --- NPI parameters the Q curve drives. q0 = value at Q = 0 (worst response),
#     q1 = value at Q = 1 (best response). Edit these min/max values freely.
NPI_SPECS <- list(
  prob_hosp         = list(q0 = 0.00, q1 = 0.80),   # P(hospitalised | symptomatic)
  delay_hosp        = list(q0 = 6.00, q1 = 1.50),   # onset->hosp delay factor (days); improves DOWN
  prop_etu          = list(q0 = 0.00, q1 = 0.90),   # proportion of hospitalised cases in an ETU
  safe_funeral_prop = list(q0 = 0.00, q1 = 0.90),   # COMMUNITY: proportion of funerals that are SAFE (-> 1-this = unsafe comm)
  unsafe_funeral_prop_hosp = list(q0 = 0.90, q1 = 1.00),  # HOSPITAL (non-ETU): UNSAFE-funeral probability, given directly
  ppe_coverage      = list(q0 = 0.00, q1 = 0.90)    # PPE coverage lever (-> ppe_coverage_hcw)
)
# Unsafe-funeral probability for ETU deaths (kept separate; ETU deaths are
# managed safely, so 0 by default).
UNSAFE_FUNERAL_ETU <- 0.0

# --- Fixed efficacy scalars (NOT time-varying; override DEFAULT_SCALAR_INPUTS).
SCALAR_OVERRIDES <- list(
  etu_efficacy                         = 0.60,
  general_hospital_quarantine_efficacy = 0.20,
  ppe_efficacy                         = 0.60,
  safe_funeral_efficacy                = 0.70
)

# --- Simulation grid + controls.
R0_GRID          <- seq(1.45, 1.75, by = 0.05)   # baseline (t=0) R0 grid
FUNERAL_FRAC     <- 0.25                          # share of t=0 transmission via funerals
SEEDING_CASES    <- 5L                            # initial seeding infections
N_STOCH          <- 25L                          # stochastic replicates per R0
# Takeoff condition: an outbreak counts as "taken off" only if it has reached at
# least TAKEOFF_N cumulative infections BY the deadline date (relative to
# EPIDEMIC_START_DATE); otherwise it is re-run (seed advanced).
TAKEOFF_N             <- 250L                     # cumulative infections required ...
TAKEOFF_DEADLINE_DATE <- as.Date("2026-06-15")     # ... by this calendar date
MAX_RETRIES      <- 50L                           # cap on re-runs per replicate
CHECK_FINAL_SIZE <- 20000L                        # stop a run once this many cases exist

# --- Summary grids (computed per run, then median across runs).
# Daily grid: cumulative cases (and the incidence we derive from them) are read
# at every day, so even very short growth-rate windows resolve. cum_at is just
# count(infections <= t), so a daily grid is cheap. Push the upper bound out if
# you set GROWTH_LATE_END beyond day 365.
TIMEPOINTS <- 0:365                               # days at which to read cumulative cases
AMOUNTS    <- c(10L, 25L, 50L, 100L, 250L, 500L,  # cumulative case amounts to time
                1000L, 2000L, 2500L, 4000L, 5000L)

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
EPIDEMIC_START_DATE <- as.Date("2026-02-17")

# --- Parallel + RNG.
N_WORKERS <- min(future::availableCores() - 3, 50L)
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
# Pick the chosen forward-extrapolation scenario (saved by 01). Each scenario
# shares the fitted curve up to the last data point and differs only afterwards;
# we use its point Q series `q` as the Q curve here.
scen_path <- file.path(DIR_OUT, "dose_q_curve_extrapolation_scenarios.rds")
if (file.exists(scen_path)) {
  all_scen <- readRDS(scen_path)
  avail    <- unique(as.character(all_scen$scenario))
  if (!EXTRAP_SCENARIO %in% avail) {
    stop("EXTRAP_SCENARIO '", EXTRAP_SCENARIO, "' not found in ",
         basename(scen_path), ". Available: ", paste(avail, collapse = ", "),
         ".", call. = FALSE)
  }
  sel     <- all_scen[as.character(all_scen$scenario) == EXTRAP_SCENARIO, , drop = FALSE]
  q_curve <- data.frame(relative_day = sel$relative_day,
                        date = as.Date(sel$date), q_mean = sel$q)
  message("Using dose Q-curve extrapolation scenario: '", EXTRAP_SCENARIO, "'.")
} else {
  # Fallback: scenarios file absent (older 01 run). Only the logistic projection
  # is available, as dose_q_curve.rds.
  if (EXTRAP_SCENARIO != "logistic") {
    stop("Extrapolation-scenarios file not found (", scen_path, "). Re-run ",
         "01_fit_dose_q_curve.R to use EXTRAP_SCENARIO = '", EXTRAP_SCENARIO,
         "'.", call. = FALSE)
  }
  q_curve <- readRDS(file.path(DIR_OUT, "dose_q_curve.rds"))
  if (is.list(q_curve) && !is.data.frame(q_curve) && !is.null(q_curve$q_curve)) {
    q_curve <- q_curve$q_curve   # tolerate being handed the full fit list
  }
  message("Scenarios file not found; using dose_q_curve.rds (logistic projection).")
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

# All time-axis plots below are drawn in CALENDAR time: a relative day `d` maps
# to the date epi_start + d. The dose-data window (first and last observation)
# is marked on every such plot with dashed verticals.
day_to_date     <- function(d) epi_start + d
DATA_FIRST_DATE <- min(DOSE_OBS$date)   # 18 May 2026 (first dose observation)
DATA_LAST_DATE  <- max(DOSE_OBS$date)   # 14 Jun 2026 (last  dose observation)
data_window_vlines <- list(
  geom_vline(xintercept = DATA_FIRST_DATE, linetype = "dashed", colour = "grey40"),
  geom_vline(xintercept = DATA_LAST_DATE,  linetype = "dashed", colour = "grey40")
)

# Takeoff deadline as a relative day (since epi_start): an outbreak must reach
# TAKEOFF_N cumulative infections on or before this day to count as taken off.
takeoff_day <- as.integer(TAKEOFF_DEADLINE_DATE - epi_start)
if (takeoff_day <= 0L) {
  stop("TAKEOFF_DEADLINE_DATE (", TAKEOFF_DEADLINE_DATE, ") must be after the ",
       "epidemic start (", epi_start, ").", call. = FALSE)
}
message(sprintf("Takeoff: >= %d infections by %s (relative day %d).",
                TAKEOFF_N, as.character(TAKEOFF_DEADLINE_DATE), takeoff_day))

# Snapshot dates for the end-of-script cumulative-infection comparison. Add their
# relative days to TIMEPOINTS so each replicate's cumulative is read off EXACTLY
# at these dates (rather than interpolated off the 10-day grid).
# Snapshot dates for the per-R0 cumulative-infection plot, paired with which
# observed reference line(s) to draw at each: 14 May (confirmed-series start ->
# both onsets and confirmed), 7 Jun (onsets finish, but confirmed is also still
# available -> both), 16 Jun (confirmed finish; onsets ended -> confirmed only).
# SNAPSHOT_REFS is matched to SNAPSHOT_DATES element-by-element; edit together.
SNAPSHOT_DATES <- as.Date(c("2026-05-14", "2026-06-07", "2026-06-16"))
SNAPSHOT_REFS  <- list(c("onsets", "confirmed cases"),  # 14 May: both available
                       c("onsets", "confirmed cases"),  # 07 Jun: onsets finish + confirmed
                       "confirmed cases")                # 16 Jun: confirmed finish (onsets ended)
snapshot_days  <- as.integer(SNAPSHOT_DATES - epi_start)
TIMEPOINTS     <- sort(unique(c(TIMEPOINTS, snapshot_days[snapshot_days > 0])))

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

safe_funeral_prop   <- clip01(lin(NPI_SPECS$safe_funeral_prop, q_on_grid))       # community
unsafe_funeral_hosp <- clip01(lin(NPI_SPECS$unsafe_funeral_prop_hosp, q_on_grid)) # hospital (non-ETU), direct

scenario_matrix_df <- data.frame(
  scenario                 = SCENARIO_ID,
  scenario_label           = SCENARIO_LABEL,
  relative_day             = matrix_days,
  prob_hosp                = clip01(lin(NPI_SPECS$prob_hosp, q_on_grid)),
  delay_hosp               = pmax(lin(NPI_SPECS$delay_hosp, q_on_grid), 0.01),
  prob_unsafe_funeral_comm = clip01(1 - safe_funeral_prop),
  prob_unsafe_funeral_hosp = unsafe_funeral_hosp,
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
# 5. Rt profile per R0 -- computed and plotted BEFORE the long fiber runs
# ----------------------------------------------------------------------------
# Quick analytic sanity check (NO simulation, so it's fast): for each starting
# R0, compute the single-type reproduction-number profile over time from that
# R0's own args -- the same NPI curves the simulations use -- with
# Rt_curve_single_type() (functions/calculate_model_approx_rt.R), and overlay all
# R0s on one graph. We expect Rt to start near each target R0 and fall below 1 as
# the dose-driven response ramps up. The function returns both the instantaneous
# (frozen-conditions) and case (cohort) reproduction numbers; both are plotted.
RT_TIMES <- 0:max(TIMEPOINTS)   # day grid (since epidemic start) for the profiles
RT_MC_N  <- 50000L              # Monte-Carlo parents (so R_inst(0) ~ target R0)
RT_SEED  <- 1L

message("Computing analytic Rt profiles for each R0 (before the runs)...")
rt_profiles <- do.call(rbind, lapply(seq_along(R0_GRID), function(i) {
  rt <- Rt_curve_single_type(args_by_r0[[i]], times = RT_TIMES,
                             n = RT_MC_N, seed = RT_SEED)
  rt$r0 <- R0_GRID[i]
  rt
}))
write.csv(rt_profiles, file.path(DIR_OUT, "dose_r0_grid_rt_profiles.csv"), row.names = FALSE)

# Long form: one row per (R0, day, mode) for the two reproduction-number curves.
rt_long <- rbind(
  data.frame(r0 = rt_profiles$r0, day = rt_profiles$time,
             mode = "instantaneous", Rt = rt_profiles$R_inst),
  data.frame(r0 = rt_profiles$r0, day = rt_profiles$time,
             mode = "case",          Rt = rt_profiles$R_case)
)
rt_long$mode <- factor(rt_long$mode, levels = c("instantaneous", "case"))
rt_long$date <- day_to_date(rt_long$day)

p_rt <- ggplot(rt_long, aes(date, Rt, colour = factor(r0), linetype = mode,
                            group = interaction(r0, mode))) +
  geom_hline(yintercept = 1, colour = "grey55", linewidth = 0.4) +
  data_window_vlines +
  geom_line(linewidth = 0.8) +
  scale_linetype_manual(values = c(instantaneous = "solid", case = "22"),
                        name = "Rt type") +
  scale_colour_viridis_d(option = "C", end = 0.9, name = expression(R[0] ~ "(t=0)")) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
  labs(title = "Analytic Rt profile by starting R0 (before the fiber runs)",
       subtitle = sprintf("Single-type approximation; solid = instantaneous, dashed = case; vertical dashes = dose-data window (%s, %s)",
                          format(DATA_FIRST_DATE, "%d %b"), format(DATA_LAST_DATE, "%d %b")),
       x = "Date", y = expression(R[t])) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(DIR_OUT, "dose_r0_grid_rt_profiles.png"), p_rt, width = 9, height = 5.5, dpi = 150)
print(p_rt)

# ----------------------------------------------------------------------------
# 6. Per-replicate runner + summariser (used on the parallel workers)
# ----------------------------------------------------------------------------
# Run ONE replicate, re-running (advancing the seed) until the outbreak has at
# least `takeoff_n` cumulative infections BY day `takeoff_day` (or `max_retries`
# is hit). Returns the realised infection times plus bookkeeping.
run_one_takeoff <- function(args, base_seed, takeoff_n, max_retries, takeoff_day) {
  seed <- base_seed; retry <- 0L; inf_times <- numeric(0)
  n_cases <- 0L; n_by_deadline <- 0L
  repeat {
    args$seed <- seed
    out <- do.call(fiber::branching_process_main, args)
    it  <- out$tdf$time_infection_absolute
    it  <- it[!is.na(it)]
    n_cases       <- length(it)
    n_by_deadline <- sum(it <= takeoff_day)        # infections on/before the deadline
    if (n_by_deadline >= takeoff_n || retry >= max_retries) { inf_times <- it; break }
    retry <- retry + 1L; seed <- seed + 1L
  }
  list(inf_times = inf_times, n_cases = n_cases, n_by_deadline = n_by_deadline,
       took_off = n_by_deadline >= takeoff_n, final_seed = seed, retries = retry)
}

# Per-run summaries computed from the raw infection times:
#   cum_at[t]  = cumulative cases by day t
#   time_to[a] = day on which the a-th infection occurred (time to `a` cases),
#                or NA if the run never reached `a` cases.
summarise_run <- function(inf_times, timepoints, amounts) {
  st <- sort(inf_times); n <- length(st)
  list(
    # cum_at[k] = count of infections on/before timepoints[k]; findInterval on the
    # sorted times is O(n + m) (matters now that timepoints is a daily grid).
    cum_at  = findInterval(timepoints, st),
    time_to = vapply(amounts,    function(a) if (n >= a) st[a] else NA_real_, numeric(1)),
    n_cases = n
  )
}

# ----------------------------------------------------------------------------
# 7. Run the grid in parallel, one R0 at a time, with progress + ETA
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
    run <- run_one_takeoff(args_i, base_seed, TAKEOFF_N, MAX_RETRIES, takeoff_day)
    s   <- summarise_run(run$inf_times, TIMEPOINTS, AMOUNTS)
    list(r0 = r0_i, r0_idx = i, rep_id = rep_id,
         n_cases = run$n_cases, n_by_deadline = run$n_by_deadline,
         took_off = run$took_off, retries = run$retries,
         cum_at = s$cum_at, time_to = s$time_to)
  },
  future.globals = list(
    args_i = args_i, i = i, r0_i = r0_i, N_STOCH = N_STOCH,
    MAX_RETRIES = MAX_RETRIES, SEED_BASE = SEED_BASE, TAKEOFF_N = TAKEOFF_N,
    takeoff_day = takeoff_day,
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
# 8. Summarise across replicates (medians from the individual raw runs)
# ----------------------------------------------------------------------------
safe_median <- function(x) if (all(is.na(x))) NA_real_ else stats::median(x, na.rm = TRUE)
safe_q <- function(x, p) if (all(is.na(x))) NA_real_ else
  as.numeric(stats::quantile(x, p, na.rm = TRUE, names = FALSE))

took   <- Filter(function(r) isTRUE(r$took_off), results)
n_fail <- length(results) - length(took)
if (n_fail > 0L)
  message(sprintf("Note: %d/%d replicates failed to reach %d infections by %s after %d retries (excluded from summaries).",
                  n_fail, length(results), TAKEOFF_N,
                  as.character(TAKEOFF_DEADLINE_DATE), MAX_RETRIES))
if (length(took) == 0L)
  stop("No replicate reached takeoff (", TAKEOFF_N, " infections by ",
       as.character(TAKEOFF_DEADLINE_DATE), "). This is a demanding bar for low ",
       "R0 / few seeds -- consider a higher R0_GRID/SEEDING_CASES, a later ",
       "TAKEOFF_DEADLINE_DATE, or a lower TAKEOFF_N.", call. = FALSE)
by_r0 <- split(took, vapply(took, function(r) r$r0_idx, integer(1)))

# (a) Mean + median cumulative cases at each timepoint.
cumulative_cases <- do.call(rbind, lapply(seq_along(R0_GRID), function(i) {
  rs <- by_r0[[as.character(i)]]
  if (is.null(rs)) return(NULL)
  M <- do.call(rbind, lapply(rs, function(r) r$cum_at))   # n_runs x n_timepoints
  data.frame(
    r0               = R0_GRID[i],
    timepoint_day    = TIMEPOINTS,
    n_runs           = nrow(M),
    mean             = apply(M, 2, mean),
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
    takeoff_deadline_date = TAKEOFF_DEADLINE_DATE, takeoff_day = takeoff_day,
    check_final_size = CHECK_FINAL_SIZE, timepoints = TIMEPOINTS, amounts = AMOUNTS,
    extrap_scenario = EXTRAP_SCENARIO,
    epidemic_start_date = epi_start, q_first_date = q_first_date,
    offset_days = offset_days
  )
), file.path(DIR_OUT, "dose_r0_grid_results.rds"))

message("\n02_npi_inputs_and_fiber_runs.R complete. Results in outputs/.")

# ----------------------------------------------------------------------------
# 9. Quick diagnostic plots (display + saved)
# ----------------------------------------------------------------------------
# Observed cumulative-onset curve (from analyses/onset_incidence/
# compute_daily_incidence_from_onsets.R) to overlay on the fiber TRAJECTORY
# plots for comparison. Those plots are CUMULATIVE, so we compare against the
# observed CUMULATIVE onsets on the same calendar-date axis. `onset_overlay` is
# NULL (a no-op when added to a ggplot) if the file has not been generated yet.
onset_csv     <- here("data-processed", "onsets_daily_incidence.csv")
onset_overlay <- NULL
if (file.exists(onset_csv)) {
  onsets_obs       <- read.csv(onset_csv, stringsAsFactors = FALSE)
  onsets_obs$date  <- as.Date(onsets_obs$date)
  onset_overlay <- geom_line(data = onsets_obs, aes(date, cumulative_onsets),
                             inherit.aes = FALSE, colour = "#d62728", linewidth = 1.1)
  message(sprintf("Overlaying observed cumulative onsets (%s to %s, max %.0f).",
                  as.character(min(onsets_obs$date)), as.character(max(onsets_obs$date)),
                  max(onsets_obs$cumulative_onsets, na.rm = TRUE)))
} else {
  message("Note: ", basename(onset_csv), " not found -- skipping observed-onset ",
          "overlay. Run analyses/onset_incidence/compute_daily_incidence_from_onsets.R.")
}

# Observed national cumulative CONFIRMED cases (INSP sitrep), overlaid on the
# same cumulative trajectory plots in a distinct colour. NULL no-op if absent.
confirmed_csv     <- here("data-processed", "insp_sitrep__national_cumulative_confirmed_cases__daily.csv")
confirmed_overlay <- NULL
if (file.exists(confirmed_csv)) {
  confirmed_obs      <- read.csv(confirmed_csv, stringsAsFactors = FALSE)
  confirmed_obs$date <- as.Date(confirmed_obs$date, format = "%d/%m/%Y")
  confirmed_overlay <- geom_line(
    data = confirmed_obs, aes(date, national_cumulative_confirmed_cases),
    inherit.aes = FALSE, colour = "#1a9850", linewidth = 1.1)
  message(sprintf("Overlaying observed cumulative confirmed cases (%s to %s, max %.0f).",
                  as.character(min(confirmed_obs$date)), as.character(max(confirmed_obs$date)),
                  max(confirmed_obs$national_cumulative_confirmed_cases, na.rm = TRUE)))
} else {
  message("Note: ", basename(confirmed_csv), " not found -- skipping confirmed-cases overlay.")
}

# McCabe external estimate: total cumulative cases by a given date, given as a
# low/high pair. Drawn on the cumulative-total comparison plots as a vertical
# range (a line) at that date, capped with a marker at each value, so it can be
# read off against where the model / data curves sit on that date. Edit freely.
MCCABE_DATE   <- as.Date("2026-05-27")
MCCABE_VALUES <- c(600, 1000)
mccabe_df   <- data.frame(date = MCCABE_DATE,
                          lo = min(MCCABE_VALUES), hi = max(MCCABE_VALUES))
mccabe_pts  <- data.frame(date = MCCABE_DATE, value = MCCABE_VALUES)
mccabe_layer <- list(
  geom_linerange(data = mccabe_df, aes(x = date, ymin = lo, ymax = hi),
                 inherit.aes = FALSE, colour = "#984ea3", linewidth = 1.0),
  geom_point(data = mccabe_pts, aes(x = date, y = value),
             inherit.aes = FALSE, colour = "#984ea3", size = 2.4, shape = 18)
)
mccabe_desc <- sprintf("purple = McCabe est. total (%s by %s)",
                       paste(MCCABE_VALUES, collapse = " & "), format(MCCABE_DATE, "%d %b"))
# Same estimate as horizontal reference lines, for the R0-vs-cumulative snapshot
# plot (x-axis is R0, not date): one dashed purple line at each value.
mccabe_hlines <- geom_hline(yintercept = MCCABE_VALUES, colour = "#984ea3",
                            linetype = "dashed", linewidth = 0.7)

# (i) The time-varying NPI inputs.
tv_long$date <- day_to_date(tv_long$relative_day)
p_inputs <- ggplot(tv_long, aes(date, value)) +
  data_window_vlines +
  geom_line(colour = "#1f77b4", linewidth = 0.8) +
  facet_wrap(~ input, scales = "free_y") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
  labs(title = "Time-varying NPI inputs driven by the dose Q curve",
       subtitle = sprintf("Vertical dashes = dose-data window (%s, %s)",
                          format(DATA_FIRST_DATE, "%d %b"), format(DATA_LAST_DATE, "%d %b")),
       x = "Date", y = "Input value") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))
ggsave(file.path(DIR_OUT, "dose_npi_timevarying.png"), p_inputs, width = 9, height = 6, dpi = 150)
print(p_inputs)

# (ii) Mean cumulative cases over time, by R0.
cumulative_cases$date <- day_to_date(cumulative_cases$timepoint_day)
p_cum <- ggplot(cumulative_cases,
                aes(date, mean , colour = factor(r0))) +
  data_window_vlines +
  geom_ribbon(aes(ymin = q25_cum_cases, ymax = q75_cum_cases, fill = factor(r0)),
              colour = NA, alpha = 0.12) +
  geom_line(linewidth = 0.8) +
  onset_overlay +
  confirmed_overlay +
  mccabe_layer +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
  labs(title = "Mean cumulative cases by baseline R0",
       subtitle = paste0("Mean (lines) + 25-75% (bands); red = observed cumulative onsets; green = confirmed cases; ", mccabe_desc),
       x = "Date", y = "Cumulative cases", colour = "R0", fill = "R0") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  +
  coord_cartesian(xlim = c(epi_start, as.Date("2026-06-20")),
                  ylim = c(0, NA))
ggsave(file.path(DIR_OUT, "dose_r0_grid_cumulative_cases.png"), p_cum, width = 9, height = 5.5, dpi = 150)
print(p_cum)

# (iii) Per-replicate cumulative-incidence trajectories, one panel per R0.
# Each thin line is ONE stochastic replicate's cumulative cases over time (read
# at TIMEPOINTS); the bold line is the mean. One facet per starting R0.
traj_long <- do.call(rbind, lapply(took, function(r)
  data.frame(r0 = r$r0, rep_id = r$rep_id, day = TIMEPOINTS, cum = r$cum_at)))
traj_long$date <- day_to_date(traj_long$day)

# Fade individual lines more when there are many replicates, less when few.
traj_alpha <- max(0.06, min(0.6, 25 / N_STOCH))

p_traj <- ggplot(traj_long, aes(date, cum, group = interaction(r0, rep_id))) +
  data_window_vlines +
  geom_line(colour = "#1f77b4", alpha = traj_alpha, linewidth = 0.35) +
  geom_line(data = cumulative_cases, aes(date, mean),
            inherit.aes = FALSE, colour = "black", linewidth = 0.9) +
  onset_overlay +
  confirmed_overlay +
  mccabe_layer +
  facet_wrap(~ r0, scales = "free_y", labeller = label_both) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
  labs(title = "Cumulative-incidence trajectories by replicate, per R0",
       subtitle = sprintf("Thin = replicates; black = mean; red = onsets; green = confirmed cases; %s; scenario '%s', %d reps",
                          mccabe_desc, EXTRAP_SCENARIO, N_STOCH),
       x = "Date", y = "Cumulative cases") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))  +
  coord_cartesian(xlim = c(epi_start, as.Date("2026-06-20")),
                  ylim = c(0, 3000))
ggsave(file.path(DIR_OUT, "dose_r0_grid_cumulative_trajectories.png"), p_traj,
       width = 10, height = 7, dpi = 150)
print(p_traj)

p_traj_log10 <- ggplot(traj_long, aes(date, cum, group = interaction(r0, rep_id))) +
  data_window_vlines +
  geom_line(colour = "#1f77b4", alpha = traj_alpha, linewidth = 0.35) +
  geom_line(data = cumulative_cases, aes(date, mean),
            inherit.aes = FALSE, colour = "black", linewidth = 0.9) +
  onset_overlay +
  confirmed_overlay +
  mccabe_layer +
  facet_wrap(~ r0, scales = "free_y", labeller = label_both) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
  scale_y_log10(labels = scales::label_number()) +
  labs(title = "Cumulative-incidence trajectories by replicate, per R0",
       subtitle = sprintf("Thin = replicates; black = mean; red = onsets; green = confirmed cases; %s; scenario '%s', %d reps",
                          mccabe_desc, EXTRAP_SCENARIO, N_STOCH),
       x = "Date", y = "Cumulative cases (log10)") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))  +
  coord_cartesian(xlim = c(epi_start, as.Date("2026-06-20")),
                  ylim = c(1, NA))
print(p_traj_log10)

# (iv) Daily incidence implied by each replicate's cumulative curve, per R0.
# Incidence over each interval = diff(cumulative) / diff(day), plotted at the
# interval midpoint. Thin = replicates; black = mean across replicates; the
# observed daily incidence is overlaid (red = onsets, green = confirmed).
inc_mids <- (TIMEPOINTS[-1] + TIMEPOINTS[-length(TIMEPOINTS)]) / 2
inc_long <- do.call(rbind, lapply(took, function(r)
  data.frame(r0 = r$r0, rep_id = r$rep_id, day = inc_mids,
             incidence = diff(r$cum_at) / diff(TIMEPOINTS))))
inc_long$date <- day_to_date(inc_long$day)
inc_mean <- aggregate(incidence ~ r0 + day, data = inc_long, FUN = mean)
inc_mean$date <- day_to_date(inc_mean$day)

onset_inc_overlay <- if (exists("onsets_obs"))
  geom_line(data = onsets_obs, aes(date, daily_incidence), inherit.aes = FALSE,
            colour = "#d62728", linewidth = 1.0) else NULL
confirmed_inc_overlay <- NULL
if (exists("confirmed_obs")) {
  co <- confirmed_obs[order(confirmed_obs$date), ]
  ci <- data.frame(date = co$date[-1],
                   incidence = diff(co$national_cumulative_confirmed_cases) /
                               as.numeric(diff(co$date)))
  confirmed_inc_overlay <- geom_line(data = ci, aes(date, incidence),
                                     inherit.aes = FALSE, colour = "#1a9850", linewidth = 1.0)
}

p_traj_inc <- ggplot(inc_long, aes(date, incidence, group = interaction(r0, rep_id))) +
  data_window_vlines +
  geom_line(colour = "#1f77b4", alpha = traj_alpha, linewidth = 0.35) +
  geom_line(data = inc_mean, aes(date, incidence), inherit.aes = FALSE,
            colour = "black", linewidth = 0.9) +
  onset_inc_overlay +
  confirmed_inc_overlay +
  facet_wrap(~ r0, scales = "free_y", labeller = label_both) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
  labs(title = "Daily incidence implied by the cumulative curves, per R0",
       subtitle = sprintf("diff(cumulative)/diff(day); thin = replicates, black = mean; red = onset, green = confirmed daily incidence; scenario '%s'",
                          EXTRAP_SCENARIO),
       x = "Date", y = "Incidence (per day)") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) +
  coord_cartesian(xlim = c(epi_start, as.Date("2026-06-20")))
ggsave(file.path(DIR_OUT, "dose_r0_grid_daily_incidence.png"), p_traj_inc,
       width = 10, height = 7, dpi = 150)
print(p_traj_inc)

# ----------------------------------------------------------------------------
# 10. Growth rate + doubling time by period (fiber curves vs data sources)
# ----------------------------------------------------------------------------
# For each curve we fit log(incidence) ~ day over a date window; the slope is
# the exponential growth rate r (per day) and the doubling time is log(2)/r.
# Working from DAILY INCIDENCE (not the cumulative curve) is the standard
# estimator: it reflects the CURRENT growth rate, so it is not contaminated by
# the seeding backlog early on and it can turn negative once incidence peaks and
# falls (a cumulative curve keeps climbing and would overstate r). Incidence:
#   * fiber     -- diff(mean cumulative) / diff(day) at the timepoint midpoints
#   * onsets    -- the daily_incidence column
#   * confirmed -- diff(cumulative confirmed) / diff(day)
# Windows (all boundaries AND labels come from the GROWTH_* config dates):
#   * "pre <MID_START>"      : start .. day before MID_START
#   * "<MID_START>-<MID_END>": the main window
#   * "<MID_END>-<LATE_END>" : the window after MID_END (fiber only)
#   * "<X>-<MID_END>"        : adjustable alternative start to MID_END
#   * "full"                 : a data source's full observed range (confirmed)
# Doubling time is only defined for r > 0 (a flat/declining curve has r <= 0);
# we always report r and set doubling time to NA when r <= 0.

# --- Growth-rate analysis (section 10) period boundaries. Edit these freely;
#     ALL period labels in the table/plot are generated automatically from them.
GROWTH_MID_START <- as.Date("2026-05-15")  # start of the main window / "pre" cutoff
GROWTH_MID_END   <- as.Date("2026-06-14")  # shared END date of the windows
GROWTH_LATE_END  <- as.Date("2026-07-14")  # end of the "after MID_END" window
GROWTH_X_DATE    <- as.Date("2026-05-26")  # adjustable alternative window start (to MID_END)

# Window boundaries + auto-generated labels from the configured GROWTH_* dates.
pre_end  <- GROWTH_MID_START - 1L
lab_pre  <- sprintf("pre %s", format(GROWTH_MID_START, "%d %b"))
lab_mid  <- sprintf("%s-%s",  format(GROWTH_MID_START, "%d %b"), format(GROWTH_MID_END, "%d %b"))
lab_late <- sprintf("%s-%s",  format(GROWTH_MID_END,   "%d %b"), format(GROWTH_LATE_END, "%d %b"))
lab_x    <- sprintf("%s-%s",  format(GROWTH_X_DATE,    "%d %b"), format(GROWTH_MID_END, "%d %b"))

# Fit r (per day) + doubling time by regressing log(incidence) on day over the
# window [lo, hi]. Non-positive / non-finite points are dropped (log undefined);
# NA is returned when fewer than two distinct days remain in the window.
fit_growth <- function(dates, values, lo, hi) {
  dates <- as.Date(dates); values <- suppressWarnings(as.numeric(values))
  keep <- is.finite(values) & values > 0 & !is.na(dates) & dates >= lo & dates <= hi
  d <- as.numeric(dates[keep]); v <- values[keep]
  if (length(unique(d)) < 2L)
    return(data.frame(growth_rate = NA_real_, doubling_time = NA_real_,
                      n_points = sum(keep), date_from = as.Date(NA), date_to = as.Date(NA)))
  r <- unname(coef(stats::lm(log(v) ~ d))[2])
  data.frame(growth_rate   = r,
             doubling_time = if (is.finite(r) && r > 0) log(2) / r else NA_real_,
             n_points      = length(v),
             date_from     = min(dates[keep]), date_to = max(dates[keep]))
}

# --- Fiber: per-R0 daily incidence = diff(mean cumulative)/diff(day) at the
#     timepoint interval midpoints, fit over each window. ------------------------
fiber_periods <- list()
fiber_periods[[lab_pre]]  <- c(epi_start,        pre_end)
fiber_periods[[lab_x]]    <- c(GROWTH_X_DATE,    GROWTH_MID_END)
fiber_periods[[lab_mid]]  <- c(GROWTH_MID_START, GROWTH_MID_END)
fiber_periods[[lab_late]] <- c(GROWTH_MID_END,   GROWTH_LATE_END)
growth_fiber <- do.call(rbind, lapply(R0_GRID, function(r0v) {
  cc <- cumulative_cases[cumulative_cases$r0 == r0v, ]
  cc <- cc[order(cc$timepoint_day), ]
  inc_day  <- (cc$timepoint_day[-1] + cc$timepoint_day[-nrow(cc)]) / 2  # interval midpoints
  inc_val  <- diff(cc$mean) / diff(cc$timepoint_day)                    # cases per day
  inc_date <- day_to_date(inc_day)
  do.call(rbind, lapply(names(fiber_periods), function(pn) {
    rng <- fiber_periods[[pn]]
    data.frame(series = sprintf("R0=%.2f", r0v), source_type = "fiber", r0 = r0v,
               period = pn, fit_growth(inc_date, inc_val, rng[1], rng[2]),
               row.names = NULL)
  }))
}))

# --- Data sources (reuse the curves read for the overlays, if present). ------
add_data_growth <- function(acc, label, dates, values, periods) {
  rows <- do.call(rbind, lapply(names(periods), function(pn) {
    rng <- periods[[pn]]
    data.frame(series = label, source_type = "data", r0 = NA_real_,
               period = pn, fit_growth(dates, values, rng[1], rng[2]), row.names = NULL)
  }))
  rbind(acc, rows)
}
growth_data <- NULL
if (exists("confirmed_obs")) {
  conf_periods <- list("full" = c(min(confirmed_obs$date), max(confirmed_obs$date)))
  conf_periods[[lab_mid]] <- c(GROWTH_MID_START, GROWTH_MID_END)
  conf_periods[[lab_x]]   <- c(GROWTH_X_DATE,    GROWTH_MID_END)
  co       <- confirmed_obs[order(confirmed_obs$date), ]
  conf_inc <- diff(co$national_cumulative_confirmed_cases) / as.numeric(diff(co$date))
  growth_data <- add_data_growth(growth_data, "confirmed cases",
    co$date[-1], conf_inc, conf_periods)
}
if (exists("onsets_obs")) {
  onset_periods <- list()
  onset_periods[[lab_pre]] <- c(min(onsets_obs$date), pre_end)
  onset_periods[[lab_mid]] <- c(GROWTH_MID_START, GROWTH_MID_END)
  onset_periods[[lab_x]]   <- c(GROWTH_X_DATE,    GROWTH_MID_END)
  growth_data <- add_data_growth(growth_data, "onsets",
    onsets_obs$date, onsets_obs$daily_incidence, onset_periods)
}

growth_tbl <- rbind(growth_fiber, growth_data)
growth_tbl$period <- factor(growth_tbl$period,
  levels = unique(c(lab_pre, lab_x, lab_mid, lab_late, "full")))

write.csv(growth_tbl, file.path(DIR_OUT, "dose_growth_rates_doubling_times.csv"), row.names = FALSE)
disp <- growth_tbl
disp$growth_rate   <- round(disp$growth_rate, 4)
disp$doubling_time <- round(disp$doubling_time, 1)
cat("\nGrowth rate (per day) + doubling time (days) by period:\n")
print(disp[order(disp$period, disp$series),
           c("series", "period", "growth_rate", "doubling_time", "n_points",
             "date_from", "date_to")], row.names = FALSE)

# --- Graphical comparison: bars per series, faceted by metric x period. ------
series_levels <- c(sprintf("R0=%.2f", R0_GRID), "confirmed cases", "onsets")
series_cols   <- c(setNames(viridisLite::viridis(length(R0_GRID), option = "C", end = 0.9),
                            sprintf("R0=%.2f", R0_GRID)),
                   "confirmed cases" = "#1a9850", "onsets" = "#d62728")
growth_tbl$series <- factor(growth_tbl$series, levels = series_levels)

growth_long <- tidyr::pivot_longer(growth_tbl, c(growth_rate, doubling_time),
                                   names_to = "metric", values_to = "value")
# Hide very long doubling times (near-zero growth) so that facet stays readable.
growth_long$value[growth_long$metric == "doubling_time" &
                  is.finite(growth_long$value) & growth_long$value > 120] <- NA
growth_long$metric <- factor(growth_long$metric, levels = c("growth_rate", "doubling_time"),
                             labels = c("Growth rate (per day)", "Doubling time (days)"))

p_growth <- ggplot(growth_long, aes(series, value, fill = series)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_col() +
  facet_grid(metric ~ period, scales = "free_y") +
  scale_fill_manual(values = series_cols, guide = "none") +
  labs(title = "Growth rate and doubling time by period (from daily incidence)",
       subtitle = sprintf("log-linear fit to incidence; fiber = per-R0 mean trajectory vs data sources; scenario '%s' (doubling times > 120 d hidden)",
                          EXTRAP_SCENARIO),
       x = NULL, y = NULL) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
ggsave(file.path(DIR_OUT, "dose_growth_rates_doubling_times.png"), p_growth,
       width = 13, height = 6.5, dpi = 150)
print(p_growth)

# ----------------------------------------------------------------------------
# 11. Cumulative infections at two snapshot dates (per R0) vs data sources
# ----------------------------------------------------------------------------
# For each replicate, its cumulative infections at each SNAPSHOT_DATES entry
# (read EXACTLY off the daily TIMEPOINTS): one dot per replicate at its R0, with
# the mean drawn as a large dot. Observed reference lines are overlaid, but only
# the source(s) listed in SNAPSHOT_REFS for that date (so 7 Jun shows onsets only
# and 16 Jun shows confirmed only).
snap_idx  <- match(snapshot_days, TIMEPOINTS)
snap_labs <- format(SNAPSHOT_DATES, "%d %b %Y")
snap_long <- do.call(rbind, lapply(took, function(r)
  data.frame(r0 = r$r0, rep_id = r$rep_id, snapshot = snap_labs,
             cum = r$cum_at[snap_idx])))
snap_long$snapshot <- factor(snap_long$snapshot, levels = snap_labs)
snap_mean <- aggregate(cum ~ r0 + snapshot, data = snap_long, FUN = mean)
snap_median <- aggregate(cum ~ r0 + snapshot, data = snap_long, FUN = stats::median)

# Observed cumulative value at a date (interpolated; held flat past the data end).
obs_at <- function(d, dates, values) {
  dates <- as.Date(dates); values <- as.numeric(values)
  if (d <= max(dates)) approx(as.numeric(dates), values, as.numeric(d))$y
  else values[which.max(dates)]
}
# Requestable observed sources (only those whose file was read are available).
obs_sources <- list()
if (exists("onsets_obs"))
  obs_sources[["onsets"]] <- list(d = onsets_obs$date, v = onsets_obs$cumulative_onsets)
if (exists("confirmed_obs"))
  obs_sources[["confirmed cases"]] <- list(
    d = confirmed_obs$date, v = confirmed_obs$national_cumulative_confirmed_cases)

# One reference row per (snapshot, requested-and-available source).
ref <- do.call(rbind, lapply(seq_along(SNAPSHOT_DATES), function(k) {
  do.call(rbind, lapply(SNAPSHOT_REFS[[k]], function(sr) {
    s <- obs_sources[[sr]]
    if (is.null(s)) return(NULL)
    data.frame(snapshot = factor(snap_labs[k], levels = snap_labs), source = sr,
               value = obs_at(SNAPSHOT_DATES[k], s$d, s$v))
  }))
}))
if (is.null(ref)) ref <- data.frame()

ref_layer <- if (nrow(ref) > 0)
  geom_hline(data = ref, aes(yintercept = value, colour = source), linewidth = 0.9) else NULL

# Describe the per-date reference scheme in the subtitle (generated from config).
ref_desc <- paste(mapply(function(lab, srcs) sprintf("%s = %s", lab, paste(srcs, collapse = " + ")),
                         snap_labs, SNAPSHOT_REFS), collapse = "; ")

p_snap <- ggplot(snap_long, aes(factor(r0), cum)) +
  geom_jitter(width = 0.12, height = 0, colour = "#1f77b4", alpha = 0.5, size = 1.5) +
  geom_point(data = snap_median, aes(factor(r0), cum), colour = "black", size = 4) +
  ref_layer +
  mccabe_hlines +
  facet_wrap(~ snapshot, scales = "free_y") +
  scale_colour_manual(values = c("onsets" = "#d62728", "confirmed cases" = "#1a9850"),
                      name = "Observed") +
  labs(title = "Cumulative infections by R0 at snapshot dates",
       subtitle = sprintf("Blue = replicates; large black dot = mean; reference lines: %s; %s (dashed); scenario '%s'",
                          ref_desc, mccabe_desc, EXTRAP_SCENARIO),
       x = "Baseline R0", y = "Cumulative infections") +
  theme_bw(base_size = 11)
ggsave(file.path(DIR_OUT, "dose_r0_grid_snapshot_cumulative.png"), p_snap,
       width = 9, height = 5, dpi = 150)
print(p_snap)

# ----------------------------------------------------------------------------
# 12. Cumulative cases RE-ZEROED to the confirmed-series start date
# ----------------------------------------------------------------------------
# Put the fiber trajectories, the observed onsets and the observed confirmed
# cases ON THE SAME BASELINE so their accumulation rates can be compared head to
# head. We take t0 = the first date of the confirmed-case series (mid-May), call
# that "day 0", and for every series plot the cumulative cases ACCRUED SINCE t0,
# i.e. value(t) - value(t0), against the number of days since t0. Every curve
# therefore departs from the origin (0, 0): the figure answers "starting from
# mid-May, how fast does each series pile up cases?" regardless of how many had
# already accumulated before t0.
#
# t0 defaults to the confirmed-series start; if that file is absent we fall back
# to the onset-series start, then to the epidemic start.
REBASE_END_DATE <- as.Date("2026-06-20")   # right-hand x limit (set NA for auto)

rebase_t0 <- if (exists("confirmed_obs") && any(!is.na(as.Date(confirmed_obs$date))))
               min(as.Date(confirmed_obs$date), na.rm = TRUE) else
             if (exists("onsets_obs") && any(!is.na(as.Date(onsets_obs$date))))
               min(as.Date(onsets_obs$date), na.rm = TRUE) else epi_start
message(sprintf("Re-zeroing cumulative curves at t0 = %s (day 0).",
                as.character(rebase_t0)))

# Re-zero one cumulative series at t0: keep dates >= t0, subtract the (linearly
# interpolated) value at t0, and guarantee an explicit (0, 0) anchor so every
# curve starts at the origin. rule = 2 holds the series flat outside its range,
# so t0 need not be one of the series' own observation dates.
rebase_cum <- function(dates, values, t0) {
  # Coerce defensively: `dates` may arrive as character (e.g. an un-converted
  # CSV column), in which case as.numeric() would silently give all-NA and break
  # approx(); as.Date() turns ISO strings into real dates (and is a no-op on a
  # Date). `values` is coerced to numeric for the same reason.
  dates  <- as.Date(dates)
  values <- suppressWarnings(as.numeric(values))
  t0     <- as.Date(t0)
  ok <- !is.na(dates) & is.finite(values); dates <- dates[ok]; values <- values[ok]
  o  <- order(dates); dates <- dates[o]; values <- values[o]
  if (length(unique(dates)) < 2L || length(t0) != 1L || is.na(t0)) return(NULL)
  v0  <- stats::approx(as.numeric(dates), values, xout = as.numeric(t0), rule = 2)$y
  sel <- dates >= t0
  out <- data.frame(days_since = as.numeric(dates[sel] - t0),
                    cum_since  = values[sel] - v0)
  if (!any(out$days_since == 0))
    out <- rbind(data.frame(days_since = 0, cum_since = 0), out)
  out[order(out$days_since), , drop = FALSE]
}

# Fiber: one re-zeroed mean trajectory per R0.
rebased_fiber <- do.call(rbind, lapply(R0_GRID, function(r0v) {
  cc <- cumulative_cases[cumulative_cases$r0 == r0v, ]
  rb <- rebase_cum(cc$date, cc$mean, rebase_t0)
  if (is.null(rb)) return(NULL)
  data.frame(series = sprintf("R0=%.2f", r0v), source_type = "fiber", rb, row.names = NULL)
}))

# Data sources: onsets and confirmed cases, re-zeroed the same way.
rebased_data <- data.frame()
if (exists("onsets_obs")) {
  rb <- rebase_cum(onsets_obs$date, onsets_obs$cumulative_onsets, rebase_t0)
  if (!is.null(rb)) rebased_data <- rbind(rebased_data,
    data.frame(series = "onsets", source_type = "data", rb, row.names = NULL))
}
if (exists("confirmed_obs")) {
  rb <- rebase_cum(confirmed_obs$date, confirmed_obs$national_cumulative_confirmed_cases, rebase_t0)
  if (!is.null(rb)) rebased_data <- rbind(rebased_data,
    data.frame(series = "confirmed cases", source_type = "data", rb, row.names = NULL))
}

rebased_all <- rbind(rebased_fiber, rebased_data)
rebased_all$series <- factor(rebased_all$series,
  levels = c(sprintf("R0=%.2f", R0_GRID), "onsets", "confirmed cases"))
rebased_cols <- c(setNames(viridisLite::viridis(length(R0_GRID), option = "C", end = 0.9),
                           sprintf("R0=%.2f", R0_GRID)),
                  "onsets" = "#d62728", "confirmed cases" = "#1a9850")

write.csv(rebased_all, file.path(DIR_OUT, "dose_r0_grid_cumulative_rebased.csv"),
          row.names = FALSE)

x_hi <- if (length(REBASE_END_DATE) && is.na(REBASE_END_DATE)) max(rebased_all$days_since) else as.numeric(REBASE_END_DATE - rebase_t0)

p_rebased <- ggplot(rebased_all, aes(days_since, cum_since, colour = series)) +
  geom_line(aes(linewidth = source_type)) +
  scale_colour_manual(values = rebased_cols, name = NULL) +
  scale_linewidth_manual(values = c(fiber = 0.7, data = 1.3), guide = "none") +
  scale_x_continuous(breaks = seq(0, ceiling(x_hi / 7) * 7, by = 7)) +
  coord_cartesian(xlim = c(0, x_hi), ylim = c(0, NA)) +
  labs(title = "Cumulative cases accrued since the confirmed-series start date",
       subtitle = sprintf("Day 0 = %s (confirmed-series start); fiber = per-R0 mean; data = thick lines; scenario '%s'",
                          format(rebase_t0, "%d %b %Y"), EXTRAP_SCENARIO),
       x = sprintf("Days since %s", format(rebase_t0, "%d %b %Y")),
       y = "Cumulative cases since day 0") +
  theme_bw(base_size = 11)
ggsave(file.path(DIR_OUT, "dose_r0_grid_cumulative_rebased.png"), p_rebased,
       width = 9, height = 5.5, dpi = 150)
print(p_rebased)
