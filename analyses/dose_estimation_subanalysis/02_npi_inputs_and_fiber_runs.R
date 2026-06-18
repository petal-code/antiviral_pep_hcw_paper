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
# Inputs : outputs/dose_q_curve_extrapolation_scenarios.rds  (from 01; pick one
#          via EXTRAP_SCENARIO). Falls back to dose_q_curve.rds (the logistic
#          projection) only if the scenarios file is absent and EXTRAP_SCENARIO
#          is "logistic".
#          data-processed/onsets_daily_incidence.csv  (optional; from
#          analyses/onset_incidence/compute_daily_incidence_from_onsets.R --
#          its observed cumulative-onset curve is overlaid on the trajectory plots).
# Output : outputs/dose_npi_scenario_matrix.csv / .rds   (the NPI inputs)
#          outputs/dose_npi_timevarying_long.csv         (tidy, for plotting)
#          outputs/dose_r0_grid_rt_profiles.csv / .png   (analytic Rt per R0, pre-run)
#          outputs/dose_r0_grid_per_run.rds              (per-replicate metrics)
#          outputs/dose_r0_grid_cumulative_cases.csv     (median cum cases)
#          outputs/dose_r0_grid_cumulative_trajectories.png (per-replicate curves by R0)
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
  safe_funeral_prop = list(q0 = 0.00, q1 = 0.90),   # proportion of funerals that are safe
  ppe_coverage      = list(q0 = 0.00, q1 = 0.90)    # PPE coverage lever (-> ppe_coverage_hcw)
)
# Unsafe-funeral probability for ETU deaths (kept separate; ETU deaths are
# managed safely, so 0 by default).
UNSAFE_FUNERAL_ETU <- 0.0

# --- Fixed efficacy scalars (NOT time-varying; override DEFAULT_SCALAR_INPUTS).
SCALAR_OVERRIDES <- list(
  etu_efficacy                         = 0.70,
  general_hospital_quarantine_efficacy = 0.30,
  ppe_efficacy                         = 0.60,
  safe_funeral_efficacy                = 0.80
)

# --- Simulation grid + controls.
R0_GRID          <- seq(1.35, 1.55, by = 0.05)   # baseline (t=0) R0 grid
FUNERAL_FRAC     <- 0.25                          # share of t=0 transmission via funerals
SEEDING_CASES    <- 5L                            # initial seeding infections
N_STOCH          <- 10L                          # stochastic replicates per R0
# Takeoff condition: an outbreak counts as "taken off" only if it has reached at
# least TAKEOFF_N cumulative infections BY the deadline date (relative to
# EPIDEMIC_START_DATE); otherwise it is re-run (seed advanced).
TAKEOFF_N             <- 1000L                     # cumulative infections required ...
TAKEOFF_DEADLINE_DATE <- as.Date("2026-06-15")     # ... by this calendar date
MAX_RETRIES      <- 50L                           # cap on re-runs per replicate
CHECK_FINAL_SIZE <- 10000L                        # stop a run once this many cases exist

# --- Summary grids (computed per run, then median across runs).
TIMEPOINTS <- c(seq(10L, 360L, by = 10L), 365L)   # days at which to read cumulative cases
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
EPIDEMIC_START_DATE <- as.Date("2026-03-01")

# --- Parallel + RNG.
N_WORKERS <- min(future::availableCores() - 4, 50L)
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
  scale_x_date(date_labels = "%b %Y") +
  labs(title = "Analytic Rt profile by starting R0 (before the fiber runs)",
       subtitle = sprintf("Single-type approximation; solid = instantaneous, dashed = case; vertical dashes = dose-data window (%s, %s)",
                          format(DATA_FIRST_DATE, "%d %b"), format(DATA_LAST_DATE, "%d %b")),
       x = "Date", y = expression(R[t])) +
  theme_bw(base_size = 11)
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
    cum_at  = vapply(timepoints, function(t) sum(st <= t), numeric(1)),
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

# (i) The time-varying NPI inputs.
tv_long$date <- day_to_date(tv_long$relative_day)
p_inputs <- ggplot(tv_long, aes(date, value)) +
  data_window_vlines +
  geom_line(colour = "#1f77b4", linewidth = 0.8) +
  facet_wrap(~ input, scales = "free_y") +
  scale_x_date(date_labels = "%b %Y") +
  labs(title = "Time-varying NPI inputs driven by the dose Q curve",
       subtitle = sprintf("Vertical dashes = dose-data window (%s, %s)",
                          format(DATA_FIRST_DATE, "%d %b"), format(DATA_LAST_DATE, "%d %b")),
       x = "Date", y = "Input value") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(DIR_OUT, "dose_npi_timevarying.png"), p_inputs, width = 9, height = 6, dpi = 150)
print(p_inputs)

# (ii) Median cumulative cases over time, by R0.
cumulative_cases$date <- day_to_date(cumulative_cases$timepoint_day)
p_cum <- ggplot(cumulative_cases,
                aes(date, median_cum_cases, colour = factor(r0))) +
  data_window_vlines +
  geom_ribbon(aes(ymin = q25_cum_cases, ymax = q75_cum_cases, fill = factor(r0)),
              colour = NA, alpha = 0.12) +
  geom_line(linewidth = 0.8) +
  onset_overlay +
  scale_x_date(date_labels = "%b %Y") +
  labs(title = "Median cumulative cases by baseline R0",
       subtitle = "Lines = median across replicates; bands = 25-75%; dashes = dose-data window; red = observed cumulative onsets",
       x = "Date", y = "Cumulative cases", colour = "R0", fill = "R0") +
  theme_bw(base_size = 11)
ggsave(file.path(DIR_OUT, "dose_r0_grid_cumulative_cases.png"), p_cum, width = 9, height = 5.5, dpi = 150)
print(p_cum)

# (iii) Per-replicate cumulative-incidence trajectories, one panel per R0.
# Each thin line is ONE stochastic replicate's cumulative cases over time (read
# at TIMEPOINTS); the bold line is the median. One facet per starting R0.
traj_long <- do.call(rbind, lapply(took, function(r)
  data.frame(r0 = r$r0, rep_id = r$rep_id, day = TIMEPOINTS, cum = r$cum_at)))
traj_long$date <- day_to_date(traj_long$day)

# Fade individual lines more when there are many replicates, less when few.
traj_alpha <- max(0.06, min(0.6, 25 / N_STOCH))

p_traj <- ggplot(traj_long, aes(date, cum, group = interaction(r0, rep_id))) +
  data_window_vlines +
  geom_line(colour = "#1f77b4", alpha = traj_alpha, linewidth = 0.35) +
  geom_line(data = cumulative_cases, aes(date, median_cum_cases),
            inherit.aes = FALSE, colour = "black", linewidth = 0.9) +
  onset_overlay +
  facet_wrap(~ r0, scales = "free_y", labeller = label_both) +
  scale_x_date(date_labels = "%b %Y") +
  labs(title = "Cumulative-incidence trajectories by replicate, per R0",
       subtitle = sprintf("Thin = replicates; black = median; red = observed cumulative onsets; scenario '%s', %d reps",
                          EXTRAP_SCENARIO, N_STOCH),
       x = "Date", y = "Cumulative cases") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(DIR_OUT, "dose_r0_grid_cumulative_trajectories.png"), p_traj,
       width = 10, height = 7, dpi = 150)
print(p_traj)
