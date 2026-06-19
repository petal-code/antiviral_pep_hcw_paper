# ============================================================================
# 05_projection_5_qcurves.R
# ============================================================================
#
# INPUTS (must exist before running)
# -----------------------------------
# outputs/dose_q_curve_extrapolation_scenarios.rds
#     Six forward-extrapolation Q curves from 01_fit_dose_q_curve.R:
#       1. linear_to_90      -- ramps to 90% dose coverage by day 100
#       2. linear_to_95      -- ramps to 95% dose coverage by day 100
#       3. linear_to_95_dec  -- ramps slowly to 95% by end of Dec 2026
#       4. logistic          -- fitted logistic projection (current trajectory)
#       5. flat              -- held flat at the last observed value (~63%)
#       6. conflict          -- drops to ~20% during a conflict episode
#
# data-processed/onsets_daily_incidence.csv            (optional, for overlays)
# data-processed/insp_sitrep__national_cumulative_     (optional, for overlays)
#     confirmed_cases__daily.csv
#
# OUTPUTS (all written to outputs/05_outputs/)
# --------------------------------------------
# 05_projections_raw.rds          Raw replicate results — saved immediately
#                                 after the parallel run, before any plotting,
#                                 so data is safe even if a plot step errors.
# 05_projections.rds              Full bundle: raw results + summaries + config
# 05_npi_inputs.png               Time-varying NPI parameter curves by scenario
# 05_rt_profiles.png              Analytic Rt curves (all R0 x scenario combos)
# 05_cumulative_trajectories.png  Individual replicate traces, faceted by scen.
# 05_cumulative_bands.png         Median + 50% / 95% CI ribbons, faceted
# 05_daily_incidence.png          Implied daily incidence
# 05_cumulative_summary.csv       Aggregated cumulative incidence (median + CIs)
# 05_summary_table.csv            Cumulative cases at Sep 1 / Dec 31 + final size
# 05_rt_profiles.csv              Analytic Rt values (numeric)
#
#
# SCRIPT STRUCTURE
# ----------------
#  1. Configuration        <- knobs: R0 values, replicates, efficacy params
#  2. Setup                <- output folder, date helpers, colour palettes
#  3. Load Q-curve scenarios
#  4. Build NPI matrices   <- translate Q -> time-varying NPI inputs
#  5. Plot NPI inputs
#  6. Assemble FIBER args  <- one args list per (R0 x scenario)
#  7. Analytic Rt profiles <- fast check Rt drops below 1 for each scenario
#  8. Runner helpers       <- run_one_takeoff(), summarise_run()
#  9. Parallel FIBER runs  <- flat future_lapply over all 4,500 combinations
# 10. Aggregate results    <- median, 25-75%, 2.5-97.5% across replicates
# 11. Observed data overlay (optional)
# 12. Build trajectory table
# 13. Plot: individual trajectories
# 14. Plot: median + CI bands
# 15. Plot: daily incidence
# 16. Summary table        <- cumulative cases at Sep 1, Dec 31, final size
# 17. Bundle and save
# ============================================================================

suppressPackageStartupMessages({
  library(here)
  library(future)
  library(future.apply)
  library(fiber)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

source(here("analyses", "dose_estimation_subanalysis", "helpers.R"))
source(here("functions", "setup_model_parameters.R"))
source(here("functions", "calculate_model_approx_r0.R"))
source(here("functions", "calculate_model_approx_rt.R"))

set.seed(123)

# ============================================================================
# 1. Configuration  <-- the knobs to change
# ============================================================================

# --- R0 grid -----------------------------------------------------------------
# Basic reproductive number: average secondary infections per case in a fully
# susceptible population with NO interventions. We run three scenarios:
# 1.50 (moderate), 1.55 (medium), 1.60 (higher). The effect of interventions
# is captured by the time-varying NPI parameters below, not by changing R0.
R0_GRID <- c(1.40, 1.45, 1.50)

# --- Replicates --------------------------------------------------------------
# Number of stochastic replicates per (R0, scenario) combination.
N_STOCH <- 3 # 10 # 250L

# --- Efficacy parameters (PI-confirmed values, scalar = 0.6) -----------------
# These are FIXED efficacies for each intervention type. They represent how
# effective each intervention is *when fully implemented* (i.e. at Q = 1).
# Values were confirmed by the PI and correspond to a scalar of 0.6 applied to
# the baseline values from 02_npi_inputs_and_fiber_runs.R:
#   ETU efficacy      : 0.60 + (1-0.60)*0.6 = 0.84
#   Gen hosp efficacy : 0.20 + (1-0.20)*0.6 = 0.68
#   PPE efficacy      : 0.60 + (1-0.60)*0.6 = 0.84
#   Funeral efficacy  : 0.70 + (1-0.70)*0.6 = 0.88
SCALAR_OVERRIDES <- list(
  etu_efficacy                         = 0.84,  # ETU isolation efficacy
  general_hospital_quarantine_efficacy = 0.3,  # general hospital quarantine
  ppe_efficacy                         = 0.84,  # HCW PPE efficacy
  safe_funeral_efficacy                = 0.88   # safe burial efficacy
)

# --- NPI specs: Q -> time-varying NPI parameters ----------------------------
# For each NPI parameter, q0 is the value when the operational response is at
# its worst (Q = 0) and q1 is the value at full response (Q = 1). The actual
# value on any given day is linearly interpolated based on the Q curve value.
#
# NOTE: unsafe_funeral_prop_hosp is the probability a hospital death leads to
# an unsafe community funeral. It is set very low (5% at worst, 1% at best)
# because hospital deaths are managed safely by the response team.
NPI_SPECS <- list(
  prob_hosp                = list(q0 = 0.00, q1 = 0.80),  # P(hospitalised | symptomatic)
  delay_hosp               = list(q0 = 6.00, q1 = 1.50),  # onset->hosp delay (days); LOWER is better
  prop_etu                 = list(q0 = 0.00, q1 = 0.90),  # fraction of hosp cases going to ETU vs general ward
  safe_funeral_prop        = list(q0 = 0.00, q1 = 0.90),  # fraction of community funerals that are safe
  unsafe_funeral_prop_hosp = list(q0 = 0.05, q1 = 0.01),  # fraction of hospital deaths with unsafe funeral
  ppe_coverage             = list(q0 = 0.00, q1 = 0.90)   # HCW PPE coverage fraction
)
UNSAFE_FUNERAL_ETU <- 0.0  # ETU deaths are always managed safely (0% unsafe)

# --- Other simulation parameters --------------------------------------------
# Fraction of baseline transmission (at t=0) attributable to funerals.
# The remainder is community/household transmission.
FUNERAL_FRAC <- 0.25

# Number of index cases to seed the epidemic with at t=0.
SEEDING_CASES <- 5L

# Hard cap on epidemic size per replicate (L = integer type, required by FIBER).
# Set to 5000L for a quick test run; 60000L for the full production run.
# Replicates that hit this cap are flagged in the summary table.
CHECK_FINAL_SIZE <- 60000L

# Takeoff criterion: a replicate must accumulate at least TAKEOFF_N infections
# by TAKEOFF_DEADLINE_DATE to be counted as a "successful" (established)
# epidemic. Replicates that fail this threshold after MAX_RETRIES attempts
# are excluded from all outputs (they represent stochastic extinctions).
TAKEOFF_N             <- 250L
TAKEOFF_DEADLINE_DATE <- as.Date("2026-06-15")
MAX_RETRIES           <- 50L  

# --- Date / time parameters -------------------------------------------------
# Day 0 of the epidemic (first known case cluster).
EPIDEMIC_START_DATE <- as.Date("2026-02-17")

# Calendar dates at which we read off cumulative case counts for the summary
# table: "how big is the epidemic by Sep 1?" and "by Dec 31?"
SUMMARY_DATES <- as.Date(c("2026-09-01", "2026-12-31"))

# PHEIC-related dates: first PHEIC alert and formal PHEIC declaration.
# These are drawn as dashed vertical lines on every calendar-axis plot.
PHEIC_DATES <- as.Date(c("2026-05-14", "2026-05-18"))

# McCabe external estimate: TOTAL cumulative cases by MCCABE_DATE, given as a
# low/high pair. Drawn on the cumulative-total plots as a vertical range at that
# date so it can be read off against the model trajectories and the data.
MCCABE_DATE   <- as.Date("2026-05-27")
MCCABE_VALUES <- c(600, 1000)

# How far ahead (in days from epidemic day 0) the NPI matrices extend.
MATRIX_HORIZON <- 730L  # ~2 years

# Days (relative to epidemic day 0) at which cumulative cases are recorded per
# replicate. A DAILY grid (matching 02) so the implied daily-incidence plot and
# the comparison to daily data are crisp; cum_at is count(infections <= t) via
# findInterval, so a daily grid is cheap.
TIMEPOINTS <- 0:730L

# Case-count thresholds for the time_to metric: "on what day did cumulative
# infections first reach 100 cases? 500 cases? ..." This captures epidemic
# speed, not just final size.
AMOUNTS <- c(100L, 500L, 1000L, 5000L, 10000L, 25000L, 50000L, 60000L)

# --- Parallelisation --------------------------------------------------------
# Use at most (total cores - 1) workers so the machine stays responsive.
# On Linux/macOS the future package uses fork-based multicore (lower overhead);
# on Windows it falls back to socket-based multisession automatically.
N_WORKERS <- min(future::availableCores() - 4L, 50L)

# Base random seed. Each replicate gets a unique derived seed so results are
# exactly reproducible regardless of how many cores are used.
SEED_BASE <- 20260619L

# ============================================================================
# 2. Setup: dates, output folder, colour palettes
# ============================================================================

out_dir <- file.path(DIR_OUT, "05_outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

epi_start   <- as.Date(EPIDEMIC_START_DATE)
# Helper: convert relative day (integer from epidemic day 0) to calendar date.
day_to_date <- function(d) epi_start + d

# Add dashed vertical lines at both PHEIC dates to every calendar-axis plot.
pheic_vlines <- lapply(PHEIC_DATES, function(d)
  geom_vline(xintercept = d, linetype = "dashed", colour = "grey40"))

# Add summary dates as timepoints so they appear exactly in the output table.
snapshot_days <- as.integer(SUMMARY_DATES - epi_start)
TIMEPOINTS    <- sort(unique(c(TIMEPOINTS, snapshot_days[snapshot_days > 0L])))
takeoff_day   <- as.integer(TAKEOFF_DEADLINE_DATE - epi_start)

# Colour palette for all 6 scenarios (Dark2-based, matches 004 script).
# Each scenario has a fixed colour so figures are visually consistent across
# scripts and the PI can compare plots directly.
scen_palette <- c(
  linear_to_90     = "#1b9e77",  # teal
  linear_to_95     = "#66a61e",  # green
  linear_to_95_dec = "#e6ab02",  # yellow/gold
  logistic         = "#7570b3",  # purple
  flat             = "#d95f02",  # orange
  conflict         = "#e7298a"   # pink
)

# Colour palette for R0 values (dark-to-light viridis so they are
# distinguishable in both colour and greyscale).
r0_palette <- setNames(
  viridisLite::viridis(length(R0_GRID), option = "C", end = 0.85),
  sprintf("%.2f", R0_GRID)
)

# ============================================================================
# 3. Load Q-curve scenarios
# ============================================================================

scen_path <- file.path(DIR_OUT, "dose_q_curve_extrapolation_scenarios.rds")
if (!file.exists(scen_path))
  stop("Run 01_fit_dose_q_curve.R first to generate: ", scen_path, call. = FALSE)

all_scen      <- readRDS(scen_path)
all_scen$date <- as.Date(all_scen$date)

# Preserve the canonical scenario order set by script 01 (factor levels).
scen_names <- if (is.factor(all_scen$scenario)) levels(all_scen$scenario) else
                as.character(unique(all_scen$scenario))

# The Q curves start on the date of the first dose observation (~May 14 2026).
# The epidemic started Feb 27 2026, so there is an offset period during which
# Q = 0 (no operational response yet).
q_first_date <- min(all_scen$date)
if (epi_start > q_first_date)
  stop("EPIDEMIC_START_DATE must be <= Q curve start (", q_first_date, ").", call. = FALSE)
offset_days <- as.integer(q_first_date - epi_start)

message(sprintf("Loaded %d scenarios: %s", length(scen_names),
                paste(scen_names, collapse = ", ")))
message(sprintf("Epidemic starts %s; Q curve starts %s -> %d days of Q=0 prepended.",
                epi_start, q_first_date, offset_days))

# Descriptive labels matching the scenario descriptions in 01_fit_dose_q_curve.R.
# Used in plot legends and facet labels throughout.
scen_label_lookup <- c(
  linear_to_90     = "1. Linear to 90% by day 100",
  linear_to_95     = "2. Linear to 95% by day 100",
  linear_to_95_dec = "3. Linear to 95% by Dec 2026",
  logistic         = "4. Logistic projection (current)",
  flat             = "5. Flat at last value",
  conflict         = "6. Conflict at day 100"
)
scen_labels <- setNames(
  ifelse(scen_names %in% names(scen_label_lookup),
         scen_label_lookup[scen_names],
         # Fallback for any unexpected scenario name not in the lookup.
         sprintf("%d. %s", seq_along(scen_names),
                 tools::toTitleCase(gsub("_", " ", scen_names)))),
  scen_names
)

# Subset palette to the scenarios actually present (named lookup, not positional).
scen_colours <- scen_palette[scen_names]

# ============================================================================
# 4. Build NPI scenario matrices
# ============================================================================
# For each scenario, translate the daily Q value into all six time-varying NPI
# parameters using the linear interpolation defined in NPI_SPECS above.
# Each matrix is one row per day (relative day from epi_start) and is saved as
# a CSV that read_scenario_matrix() can parse.

lin         <- function(spec, q) spec$q0 + (spec$q1 - spec$q0) * q
matrix_days <- 0:MATRIX_HORIZON

build_matrix_df <- function(sn) {
  # Interpolate Q values from the scenario data onto the daily grid.
  sel    <- all_scen[as.character(all_scen$scenario) == sn, , drop = FALSE]
  q_sim  <- as.integer(sel$date - epi_start)
  q_vals <- approx(q_sim, sel$q, xout = matrix_days, rule = 2)$y
  # Pre-Q-curve period (epidemic start -> first dose data): set Q = 0 (no response).
  q_vals[matrix_days < offset_days] <- 0
  q_vals <- clip01(q_vals)

  # Derive funeral-related parameters from Q.
  safe_fp   <- clip01(lin(NPI_SPECS$safe_funeral_prop,        q_vals))
  unsafe_fh <- clip01(lin(NPI_SPECS$unsafe_funeral_prop_hosp, q_vals))

  data.frame(
    scenario                 = sn,
    scenario_label           = sn,
    relative_day             = matrix_days,
    # Probability a symptomatic case is hospitalised (improves with Q).
    prob_hosp                = clip01(lin(NPI_SPECS$prob_hosp,    q_vals)),
    # Days from symptom onset to hospitalisation (decreases with Q = improves).
    delay_hosp               = pmax(lin(NPI_SPECS$delay_hosp,     q_vals), 0.01),
    # Probability a community death has an unsafe funeral (decreases with Q).
    prob_unsafe_funeral_comm = clip01(1 - safe_fp),
    # Probability a hospital death has an unsafe funeral (low, improves slowly).
    prob_unsafe_funeral_hosp = unsafe_fh,
    # ETU deaths: always managed safely, so always 0% unsafe.
    prob_unsafe_funeral_etu  = clip01(rep(UNSAFE_FUNERAL_ETU, length(matrix_days))),
    # Fraction of hospitalised cases admitted to an ETU (vs general ward).
    prop_etu                 = clip01(lin(NPI_SPECS$prop_etu,     q_vals)),
    # HCW PPE coverage (feeds into ppe_coverage_hcw in FIBER).
    ipc_helper               = clip01(lin(NPI_SPECS$ppe_coverage, q_vals)),
    q_value                  = q_vals,
    stringsAsFactors = FALSE
  )
}

mat_dir <- file.path(out_dir, "scenario_matrices")
dir.create(mat_dir, recursive = TRUE, showWarnings = FALSE)

matrix_csvs <- setNames(vapply(scen_names, function(sn) {
  df  <- build_matrix_df(sn)
  pth <- file.path(mat_dir, sprintf("matrix_%s.csv", sn))
  write.csv(df, pth, row.names = FALSE)
  pth
}, character(1)), scen_names)

message("Built NPI matrices for ", length(scen_names), " scenarios.")

# ============================================================================
# 5. Plot NPI input curves
# ============================================================================
# Visualise how each NPI parameter varies over time under each scenario.
# This is a useful sanity check: the curves should diverge only AFTER the last
# observed Q data point (around Jun 14 2026).

npi_long <- do.call(rbind, lapply(scen_names, function(sn) {
  read.csv(matrix_csvs[[sn]], stringsAsFactors = FALSE) %>%
    select(relative_day, q_value, prob_hosp, delay_hosp,
           prob_unsafe_funeral_comm, prob_unsafe_funeral_hosp,
           prop_etu, ipc_helper) %>%
    pivot_longer(-relative_day, names_to = "input", values_to = "value") %>%
    mutate(scenario = sn, date = day_to_date(relative_day))
}))
npi_long$scenario <- factor(npi_long$scenario, levels = scen_names)

p_inputs <- ggplot(npi_long, aes(date, value, colour = scenario)) +
  pheic_vlines +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ input, scales = "free_y", ncol = 3) +
  scale_colour_manual(values = scen_colours, labels = scen_labels, name = "Scenario") +
  scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") +
  labs(title = "Time-varying NPI inputs by scenario",
       subtitle = "Scenarios are identical up to the last Q observation (~Jun 14 2026), then diverge",
       x = "Date", y = "Value") +
  theme_bw(base_size = 10) +
  theme(axis.text.x   = element_text(angle = 45, hjust = 1, size = 6),
        legend.position = "bottom")
print(p_inputs)
message("Saved 05_npi_inputs.png")

# ============================================================================
# 6. Assemble FIBER args: one per (R0, scenario) combination
# ============================================================================
# FIBER's branching_process_main() takes a large list of named parameters.
# make_model_parameters() builds that list from a scenario matrix + overrides.
#
# R0 INVARIANTS: The relationship between R0 and the offspring mean parameters
# (mn_offspring_genPop, mn_offspring_funeral) is mediated by D (the
# generation-time distribution) and F_fun (the funeral-transmission function).
# Both depend only on the fixed efficacy scalars (SCALAR_OVERRIDES), NOT on
# the time-varying Q scenario, so we compute them once from any scenario and
# reuse them across all (R0, scenario) combinations.

overrides_base <- modifyList(
  list(seeding_cases = SEEDING_CASES, check_final_size = CHECK_FINAL_SIZE),
  SCALAR_OVERRIDES
)

# Use the first scenario as the reference for computing R0 invariants.
ref_sm <- read_scenario_matrix(matrix_csvs[[scen_names[1]]])
mp_ref <- make_model_parameters(scen_names[1], ref_sm, overrides = overrides_base)

# Compute the R0 invariants (expensive Monte Carlo step, done once).
inv   <- compute_R0_invariants(mp_ref$args, n = 50000L, seed = 42L)
D     <- D_from_invariants(inv, mp_ref$args$etu_efficacy,
                           mp_ref$args$general_hospital_quarantine_efficacy)
F_fun <- F_from_invariants(inv, mp_ref$args$safe_funeral_efficacy)

# For each R0 value, solve for the community and funeral offspring means that
# are consistent with that R0 given D and F_fun. Then build a full args list
# for every (R0, scenario) combination.
args_grid <- setNames(
  do.call(c, lapply(seq_along(R0_GRID), function(ri) {
    r0    <- R0_GRID[ri]
    # solve_offspring_means finds mn_genPop and mn_funeral such that the
    # implied R0 (computed from the branching process) equals r0.
    means <- solve_offspring_means(r0, FUNERAL_FRAC, D, F_fun)
    message(sprintf("R0 = %.2f: mn_genPop = %.4f, mn_funeral = %.4f",
                    r0, means$mn_genPop, means$mn_funeral))
    lapply(scen_names, function(sn) {
      sm <- read_scenario_matrix(matrix_csvs[[sn]])
      mp <- make_model_parameters(sn, sm, overrides = overrides_base)
      a  <- mp$args
      a$mn_offspring_genPop  <- means$mn_genPop
      a$mn_offspring_funeral <- means$mn_funeral
      a$seed                 <- NULL  # seed is set per-replicate in Section 9
      a
    })
  })),
  # Key format: "R0_1.50__linear_to_90", "R0_1.50__logistic", etc.
  as.vector(outer(sprintf("R0_%.2f", R0_GRID), scen_names,
                  function(r, s) paste(r, s, sep = "__")))
)
message(sprintf("Assembled %d (R0 x scenario) fiber arg sets.", length(args_grid)))

# ============================================================================
# 7. Analytic Rt profiles (fast sanity check before running FIBER)
# ============================================================================
# Rt_curve_single_type() computes the ANALYTIC (not simulated) effective
# reproductive number Rt at each point in time, given the time-varying NPI
# parameters. This is much faster than simulation and serves as a crucial
# sanity check: if Rt never drops below 1 for any scenario, the epidemic will
# always grow to the CHECK_FINAL_SIZE cap regardless of R0. The plot below
# lets the PI verify that at least some scenarios achieve Rt < 1.
#
# Rt is shown for all (R0, scenario) combinations. Facets = R0 values;
# colours = scenarios.

RT_TIMES <- 0:365L   # compute Rt for the first year of the epidemic
RT_MC_N  <- 50000L   # Monte Carlo draws for the analytic calculation
RT_SEED  <- 1L

message("Computing analytic Rt profiles (sanity check before FIBER)...")
rt_profiles <- do.call(rbind, lapply(R0_GRID, function(r0) {
  do.call(rbind, lapply(scen_names, function(sn) {
    key <- paste(sprintf("R0_%.2f", r0), sn, sep = "__")
    rt  <- Rt_curve_single_type(args_grid[[key]], times = RT_TIMES,
                                n = RT_MC_N, seed = RT_SEED)
    data.frame(r0 = r0, scenario = sn, day = rt$time,
               R_inst = rt$R_inst, R_case = rt$R_case,
               stringsAsFactors = FALSE)
  }))
}))
rt_profiles$date     <- day_to_date(rt_profiles$day)
rt_profiles$scenario <- factor(rt_profiles$scenario, levels = scen_names)
rt_profiles$r0_label <- sprintf("R0 = %.2f", rt_profiles$r0)

write.csv(rt_profiles, file.path(out_dir, "05_rt_profiles.csv"), row.names = FALSE)

p_rt <- ggplot(rt_profiles, aes(date, R_inst, colour = scenario, group = scenario)) +
  # Dashed black line at Rt = 1: above = epidemic growing, below = controlled.
  geom_hline(yintercept = 1, linetype = "dashed", colour = "black", linewidth = 0.5) +
  pheic_vlines +
  geom_line(linewidth = 0.75) +
  facet_wrap(~ r0_label, ncol = 3) +
  scale_colour_manual(values = scen_colours, labels = scen_labels, name = "Scenario") +
  scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") +
  labs(
    title    = "Analytic instantaneous Rt by R0 and NPI scenario",
    subtitle = paste0(
      "ETU=0.84, GenHosp=0.68, PPE=0.84, SafeFuneral=0.88\n",
      "Black dashed = Rt 1 (above = growing, below = controlled). ",
      "Grey dashed = PHEIC dates."
    ),
    x = "Date", y = expression(R[t] ~ "(instantaneous)")
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x   = element_text(angle = 45, hjust = 1, size = 7),
        legend.position = "bottom")
ggsave(file.path(out_dir, "05_rt_profiles.png"), p_rt,
       width = 13, height = 5, dpi = 150)
print(p_rt)
message("Saved 05_rt_profiles.png")

# ============================================================================
# 8. Runner helpers
# ============================================================================

# run_one_takeoff()
# -----------------
# Runs a single FIBER replicate, retrying with a new seed if the epidemic
# fails to "take off" (i.e. does not reach TAKEOFF_N cumulative infections by
# TAKEOFF_DEADLINE_DATE). Stochastic epidemics can go extinct early by chance,
# so retrying is standard practice to obtain a sample of "established" outbreaks.
#
# Arguments:
#   args          -- FIBER args list (from args_grid above)
#   base_seed     -- starting random seed for this replicate
#   takeoff_n     -- minimum infections required by the deadline
#   max_retries   -- maximum number of seed increments before giving up
#   takeoff_day   -- relative day corresponding to TAKEOFF_DEADLINE_DATE
#
# Returns a list with:
#   inf_times     -- sorted infection times (relative days) for all cases
#   n_cases       -- total epidemic size
#   took_off      -- TRUE if the takeoff criterion was met
#   retries       -- how many seed increments were needed
run_one_takeoff <- function(args, base_seed, takeoff_n, max_retries, takeoff_day) {
  seed <- base_seed; retry <- 0L
  inf_times <- numeric(0); n_cases <- 0L; n_by_deadline <- 0L
  repeat {
    args$seed <- seed
    out <- do.call(fiber::branching_process_main, args)
    it  <- out$tdf$time_infection_absolute
    it  <- it[!is.na(it)]
    n_cases       <- length(it)
    n_by_deadline <- sum(it <= takeoff_day)
    if (n_by_deadline >= takeoff_n || retry >= max_retries) { inf_times <- it; break }
    retry <- retry + 1L; seed <- seed + 1L
  }
  list(inf_times = inf_times, n_cases = n_cases, n_by_deadline = n_by_deadline,
       took_off = n_by_deadline >= takeoff_n, final_seed = seed, retries = retry)
}

# summarise_run()
# ---------------
# Given a vector of individual infection times (relative days), computes:
#   cum_at  -- cumulative case count at each timepoint in TIMEPOINTS
#   time_to -- relative day when cumulative infections first reached each
#              threshold in AMOUNTS (NA if that threshold was never reached)
#   n_cases -- total epidemic size
summarise_run <- function(inf_times, timepoints, amounts) {
  st <- sort(inf_times); n <- length(st)
  list(
    # count(infections <= t) via findInterval on the sorted times: O(n+m),
    # which matters now that timepoints is a daily grid.
    cum_at  = findInterval(timepoints, st),
    time_to = vapply(amounts,    function(a) if (n >= a) st[a] else NA_real_, numeric(1)),
    n_cases = n
  )
}

# ============================================================================
# 9. Parallel FIBER runs: all (R0 x scenario x replicate) in one flat batch
# ============================================================================
# All 4,500 runs are submitted to the parallel scheduler in a single
# future_lapply() call. This is more efficient than running one scenario at a
# time because the scheduler can pack every available core continuously with no
# idle gaps between scenario batches.
#
# SEED SCHEME: each (R0, scenario, replicate) combination has a unique seed
# block of (MAX_RETRIES + 1) consecutive seeds reserved for it, so seeds never
# collide even if the maximum number of retries is reached.
#
# DATA SAFETY: results_flat (raw output) is saved to disk IMMEDIATELY after
# the parallel run completes, before any plotting code runs. If a downstream
# plot step fails, the data is already safe on disk.

n_r0   <- length(R0_GRID)
n_scen <- length(scen_names)
total  <- n_r0 * n_scen * N_STOCH

# Build the flat run grid: one row per (R0, scenario, replicate).
run_grid <- do.call(rbind, lapply(seq_along(R0_GRID), function(ri)
  do.call(rbind, lapply(seq_along(scen_names), function(si)
    data.frame(ri = ri, r0 = R0_GRID[ri], si = si, sn = scen_names[si],
               rep_id = seq_len(N_STOCH), stringsAsFactors = FALSE)))))

if (future::supportsMulticore()) {
  plan(multicore, workers = N_WORKERS)
  message(sprintf("Using multicore (fork-based) on %d workers.", N_WORKERS))
} else {
  plan(multisession, workers = N_WORKERS)
  message(sprintf("Using multisession (socket-based) on %d workers.", N_WORKERS))
}

message(sprintf(
  "Starting FIBER: %d R0 values x %d scenarios x %d reps = %d total runs.",
  n_r0, n_scen, N_STOCH, total))
t_start <- proc.time()

results_flat <- future_lapply(seq_len(nrow(run_grid)), function(idx) {
  ri     <- run_grid$ri[idx]
  r0     <- run_grid$r0[idx]
  si     <- run_grid$si[idx]
  sn     <- run_grid$sn[idx]
  rep_id <- run_grid$rep_id[idx]
  key    <- paste(sprintf("R0_%.2f", r0), sn, sep = "__")

  # Unique seed block for this (R0, scenario, replicate).
  base_seed <- SEED_BASE +
    ((ri - 1L) * n_scen + (si - 1L)) * N_STOCH * (MAX_RETRIES + 1L) +
    (rep_id - 1L) * (MAX_RETRIES + 1L)

  run <- run_one_takeoff(args_grid[[key]], base_seed, TAKEOFF_N, MAX_RETRIES, takeoff_day)
  s   <- summarise_run(run$inf_times, TIMEPOINTS, AMOUNTS)

  list(r0 = r0, r0_idx = ri, scenario = sn, scen_idx = si, rep_id = rep_id,
       n_cases = run$n_cases, n_by_deadline = run$n_by_deadline,
       took_off = run$took_off, retries = run$retries,
       cum_at = s$cum_at, time_to = s$time_to)
},
future.globals = list(
  run_grid = run_grid, args_grid = args_grid,
  n_scen = n_scen, N_STOCH = N_STOCH,
  MAX_RETRIES = MAX_RETRIES, SEED_BASE = SEED_BASE,
  TAKEOFF_N = TAKEOFF_N, takeoff_day = takeoff_day,
  TIMEPOINTS = TIMEPOINTS, AMOUNTS = AMOUNTS,
  run_one_takeoff = run_one_takeoff, summarise_run = summarise_run
),
future.packages = "fiber", future.seed = TRUE)

plan(sequential)
elapsed <- proc.time() - t_start
message(sprintf("FIBER complete: %d replicates in %.1f min.",
                length(results_flat), elapsed["elapsed"] / 60))

# --- CHECKPOINT SAVE (raw results, before any plotting) --------------------
# Saved immediately so the raw data is safe even if a downstream plot fails.
# To reload: results_flat <- readRDS(file.path(out_dir, "05_projections_raw.rds"))
saveRDS(results_flat, file.path(out_dir, "05_projections_raw.rds"))
message("Checkpoint: raw results saved to 05_projections_raw.rds")

# ============================================================================
# 10. Aggregate replicate results
# ============================================================================
# Compute median and uncertainty intervals across replicates for each
# (R0, scenario, timepoint) combination.

safe_q <- function(x, p) if (all(is.na(x))) NA_real_ else
  as.numeric(stats::quantile(x, p, na.rm = TRUE, names = FALSE))

# Exclude replicates that failed the takeoff criterion (stochastic extinctions).
took   <- Filter(function(r) isTRUE(r$took_off), results_flat)
n_fail <- length(results_flat) - length(took)
if (n_fail > 0L)
  message(sprintf(
    "Note: %d/%d replicates failed the takeoff criterion (excluded from all plots/tables).\n",
    n_fail, length(results_flat)))

# Split the successful replicates by (R0, scenario) key.
combo_key <- function(r) paste(sprintf("%.2f", r$r0), r$scenario, sep = "__")
by_combo  <- split(took, vapply(took, combo_key, character(1)))

cumulative_summary <- do.call(rbind, lapply(R0_GRID, function(r0) {
  do.call(rbind, lapply(scen_names, function(sn) {
    key <- paste(sprintf("%.2f", r0), sn, sep = "__")
    rs  <- by_combo[[key]]
    if (is.null(rs) || length(rs) == 0L) return(NULL)
    # Matrix: rows = replicates, columns = timepoints.
    M <- do.call(rbind, lapply(rs, function(r) r$cum_at))
    data.frame(r0 = r0, scenario = sn, timepoint_day = TIMEPOINTS,
               n_runs = nrow(M),
               mean   = apply(M, 2, mean),
               median = apply(M, 2, stats::median),
               q025   = apply(M, 2, safe_q, 0.025),  # 2.5th percentile
               q25    = apply(M, 2, safe_q, 0.25),   # 25th percentile
               q75    = apply(M, 2, safe_q, 0.75),   # 75th percentile
               q975   = apply(M, 2, safe_q, 0.975),  # 97.5th percentile
               stringsAsFactors = FALSE)
  }))
}))
cumulative_summary$date     <- day_to_date(cumulative_summary$timepoint_day)
cumulative_summary$scenario <- factor(cumulative_summary$scenario, levels = scen_names)
cumulative_summary$r0_label <- sprintf("R0 = %.2f", cumulative_summary$r0)

write.csv(cumulative_summary, file.path(out_dir, "05_cumulative_summary.csv"),
          row.names = FALSE)
message("Saved 05_cumulative_summary.csv")

# ============================================================================
# 11. Optional observed data overlays
# ============================================================================
# If the processed data CSVs exist, overlay observed cumulative onsets and
# confirmed cases on the projection plots. These are not required to run the
# script; the overlays are simply skipped if the files are absent.

onset_overlay <- confirmed_overlay <- onset_inc_overlay <- NULL

onset_csv <- here("data-processed", "onsets_daily_incidence.csv")
if (file.exists(onset_csv)) {
  onsets_obs      <- read.csv(onset_csv, stringsAsFactors = FALSE)
  onsets_obs$date <- as.Date(onsets_obs$date)
  onset_overlay   <- geom_line(data = onsets_obs,
                               aes(date, cumulative_onsets),
                               inherit.aes = FALSE,
                               colour = "#d62728", linewidth = 1.1, linetype = "solid")
  if ("daily_incidence" %in% names(onsets_obs))
    onset_inc_overlay <- geom_line(data = onsets_obs,
                                   aes(date, daily_incidence),
                                   inherit.aes = FALSE,
                                   colour = "#d62728", linewidth = 1.0)
  message(sprintf("Overlaying observed onsets (max %.0f).",
                  max(onsets_obs$cumulative_onsets, na.rm = TRUE)))
}

confirmed_csv <- here("data-processed",
                      "insp_sitrep__national_cumulative_confirmed_cases__daily.csv")
if (file.exists(confirmed_csv)) {
  confirmed_obs      <- read.csv(confirmed_csv, stringsAsFactors = FALSE)
  confirmed_obs$date <- as.Date(confirmed_obs$date, format = "%d/%m/%Y")
  confirmed_overlay  <- geom_line(
    data = confirmed_obs,
    aes(date, national_cumulative_confirmed_cases),
    inherit.aes = FALSE, colour = "#1a9850", linewidth = 1.1)
  message(sprintf("Overlaying confirmed cases (max %.0f).",
                  max(confirmed_obs$national_cumulative_confirmed_cases, na.rm = TRUE)))
}

# Confirmed DAILY incidence (diff of cumulative) for the daily-incidence plot.
confirmed_inc_overlay <- NULL
if (exists("confirmed_obs")) {
  cc_ord <- confirmed_obs[order(confirmed_obs$date), ]
  confirmed_inc_df <- data.frame(
    date = cc_ord$date[-1],
    inc  = diff(cc_ord$national_cumulative_confirmed_cases) / as.numeric(diff(cc_ord$date)))
  confirmed_inc_overlay <- geom_line(data = confirmed_inc_df, aes(date, inc),
                                     inherit.aes = FALSE, colour = "#1a9850", linewidth = 1.0)
}

# McCabe external estimate as a vertical range (+ end points) at MCCABE_DATE,
# for the CUMULATIVE-total plots (a single absolute total, so not the daily plot).
mccabe_df    <- data.frame(date = MCCABE_DATE, lo = min(MCCABE_VALUES), hi = max(MCCABE_VALUES))
mccabe_pts   <- data.frame(date = MCCABE_DATE, value = MCCABE_VALUES)
mccabe_layer <- list(
  geom_linerange(data = mccabe_df, aes(x = date, ymin = lo, ymax = hi),
                 inherit.aes = FALSE, colour = "#984ea3", linewidth = 1.0),
  geom_point(data = mccabe_pts, aes(x = date, y = value),
             inherit.aes = FALSE, colour = "#984ea3", size = 2.4, shape = 18)
)
mccabe_desc  <- sprintf("Purple = McCabe est. total (%s by %s)",
                        paste(MCCABE_VALUES, collapse = " & "), format(MCCABE_DATE, "%d %b"))

# ============================================================================
# 12. Build replicate trajectory long table
# ============================================================================
traj_long <- do.call(rbind, lapply(took, function(r)
  data.frame(r0 = r$r0, scenario = r$scenario, rep_id = r$rep_id,
             day = TIMEPOINTS, cum = r$cum_at, stringsAsFactors = FALSE)))
traj_long$date     <- day_to_date(traj_long$day)
traj_long$scenario <- factor(traj_long$scenario, levels = scen_names)
traj_long$r0_label <- sprintf("R0 = %.2f", traj_long$r0)

# Reduce alpha (transparency) for individual lines as N_STOCH increases to
# prevent overplotting from hiding the density distribution.
traj_alpha <- max(0.03, min(0.25, 15 / N_STOCH))

# ============================================================================
# 13. Plot: cumulative trajectories (individual replicates), faceted by scenario
# ============================================================================
# Each thin line = one replicate. Bold line = median across all replicates.
# Facets = scenario; colour = R0 value.
# Red line = observed onsets; green line = confirmed cases (if data loaded).

p_traj <- ggplot(traj_long,
                 aes(date, cum, colour = r0_label,
                     group = interaction(r0, rep_id))) +
  pheic_vlines +
  geom_line(alpha = traj_alpha, linewidth = 0.25) +
  geom_line(data = cumulative_summary,
            aes(date, median, colour = r0_label, group = r0_label),
            inherit.aes = FALSE, linewidth = 1.0) +
  onset_overlay + confirmed_overlay + mccabe_layer +
  facet_wrap(~ scenario, ncol = 2, scales = "free_y",
             labeller = labeller(scenario = scen_labels)) +
  scale_colour_manual(values = r0_palette, name = expression(R[0])) +
  scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") +
  labs(
    title    = sprintf("Cumulative incidence: %d replicates per R0 x scenario", N_STOCH),
    subtitle = paste0(
      "Thin lines = replicates  |  Bold lines = median  |  ",
      "Red = observed onsets  |  Green = confirmed cases  |  ", mccabe_desc
    ),
    x = "Date", y = "Cumulative cases"
  ) +
  theme_bw(base_size = 10) +
  theme(axis.text.x   = element_text(angle = 45, hjust = 1, size = 6),
        legend.position = "bottom")
ggsave(file.path(out_dir, "05_cumulative_trajectories.png"), p_traj,
       width = 12, height = 10, dpi = 150)
print(p_traj)
message("Saved 05_cumulative_trajectories.png")

# ============================================================================
# 14. Plot: median + 95% CI bands, faceted by scenario
# ============================================================================
# Inner ribbon = interquartile range (25-75%).
# Outer ribbon = 95% CI (2.5-97.5%).
# Bold line    = median.

p_bands <- ggplot(cumulative_summary,
                  aes(date, colour = r0_label, fill = r0_label)) +
  pheic_vlines +
  geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.10, colour = NA) +
  geom_ribbon(aes(ymin = q25,  ymax = q75),  alpha = 0.20, colour = NA) +
  geom_line(aes(y = median), linewidth = 1.0) +
  onset_overlay + confirmed_overlay + mccabe_layer +
  facet_wrap(~ scenario, ncol = 2, scales = "free_y",
             labeller = labeller(scenario = scen_labels)) +
  scale_colour_manual(values = r0_palette, name = expression(R[0])) +
  scale_fill_manual(values   = r0_palette, name = expression(R[0])) +
  scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") +
  labs(
    title    = "Cumulative incidence: median + 50% and 95% CI",
    subtitle = paste0(
      "Inner ribbon = 25-75%  |  Outer ribbon = 2.5-97.5%  |  ",
      "Red = onsets  |  Green = confirmed cases  |  ", mccabe_desc
    ),
    x = "Date", y = "Cumulative cases"
  ) +
  theme_bw(base_size = 10) +
  theme(axis.text.x   = element_text(angle = 45, hjust = 1, size = 6),
        legend.position = "bottom")
ggsave(file.path(out_dir, "05_cumulative_bands.png"), p_bands,
       width = 12, height = 10, dpi = 150)
print(p_bands)
message("Saved 05_cumulative_bands.png")

# ============================================================================
# 15. Plot: implied daily incidence
# ============================================================================
# Derived from the cumulative trajectories by differencing consecutive
# timepoints and dividing by the number of days between them.

inc_long <- do.call(rbind, lapply(took, function(r) {
  cum <- r$cum_at
  data.frame(
    r0       = r$r0,
    scenario = r$scenario,
    rep_id   = r$rep_id,
    # Use mid-point of each interval as the "date" for the incidence value.
    day      = (TIMEPOINTS[-1] + TIMEPOINTS[-length(TIMEPOINTS)]) / 2,
    inc      = diff(cum) / diff(TIMEPOINTS),
    stringsAsFactors = FALSE
  )
}))
inc_long$date     <- day_to_date(inc_long$day)
inc_long$scenario <- factor(inc_long$scenario, levels = scen_names)
inc_long$r0_label <- sprintf("R0 = %.2f", inc_long$r0)

inc_mean <- aggregate(inc ~ r0 + r0_label + scenario + day, data = inc_long, FUN = mean)
inc_mean$date     <- day_to_date(inc_mean$day)
inc_mean$scenario <- factor(inc_mean$scenario, levels = scen_names)

p_inc <- ggplot(inc_long,
                aes(date, inc, colour = r0_label,
                    group = interaction(r0, rep_id))) +
  pheic_vlines +
  geom_line(alpha = traj_alpha, linewidth = 0.25) +
  geom_line(data = inc_mean,
            aes(date, inc, colour = r0_label, group = r0_label),
            inherit.aes = FALSE, linewidth = 1.0) +
  onset_inc_overlay + confirmed_inc_overlay +
  facet_wrap(~ scenario, ncol = 2, scales = "free_y",
             labeller = labeller(scenario = scen_labels)) +
  scale_colour_manual(values = r0_palette, name = expression(R[0])) +
  scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") +
  labs(
    title    = "Daily incidence by R0 and NPI scenario",
    subtitle = "Thin lines = replicates  |  Bold lines = mean  |  Red = observed onsets  |  Green = confirmed daily incidence",
    x = "Date", y = "Daily incidence (cases/day)"
  ) +
  theme_bw(base_size = 10) +
  theme(axis.text.x   = element_text(angle = 45, hjust = 1, size = 6),
        legend.position = "bottom")
ggsave(file.path(out_dir, "05_daily_incidence.png"), p_inc,
       width = 12, height = 10, dpi = 150)
print(p_inc)
message("Saved 05_daily_incidence.png")

# ============================================================================
# 16. Summary table: cumulative cases at key dates + total epidemic size
# ============================================================================
# For each (R0, scenario) combination, reports:
#   - Median and IQR/95% CI of cumulative cases at each SUMMARY_DATE.
#   - Final epidemic size (total infections; flagged if many hit the cap).

summary_rows <- do.call(rbind, lapply(R0_GRID, function(r0) {
  do.call(rbind, lapply(scen_names, function(sn) {
    key <- paste(sprintf("%.2f", r0), sn, sep = "__")
    rs  <- by_combo[[key]]
    if (is.null(rs) || length(rs) == 0L) return(NULL)
    do.call(rbind, lapply(seq_along(SUMMARY_DATES), function(di) {
      day <- as.integer(SUMMARY_DATES[di] - epi_start)
      # Find the closest timepoint to the summary date.
      idx <- match(day, TIMEPOINTS)
      if (is.na(idx)) idx <- which.min(abs(TIMEPOINTS - day))
      vals <- vapply(rs, function(r) r$cum_at[idx], numeric(1))
      data.frame(r0 = r0, scenario = sn,
                 date   = as.character(SUMMARY_DATES[di]),
                 n_reps = length(vals),
                 median = round(stats::median(vals, na.rm = TRUE)),
                 q25    = round(safe_q(vals, 0.25)),
                 q75    = round(safe_q(vals, 0.75)),
                 q025   = round(safe_q(vals, 0.025)),
                 q975   = round(safe_q(vals, 0.975)),
                 stringsAsFactors = FALSE)
    }))
  }))
}))

final_rows <- do.call(rbind, lapply(R0_GRID, function(r0) {
  do.call(rbind, lapply(scen_names, function(sn) {
    key  <- paste(sprintf("%.2f", r0), sn, sep = "__")
    rs   <- by_combo[[key]]
    if (is.null(rs) || length(rs) == 0L) return(NULL)
    vals    <- vapply(rs, function(r) as.numeric(r$n_cases), numeric(1))
    # Count replicates that hit the cap — a high fraction suggests the true
    # final size is larger than CHECK_FINAL_SIZE.
    n_trunc <- sum(vals >= CHECK_FINAL_SIZE)
    data.frame(r0 = r0, scenario = sn,
               date   = sprintf("final_size [%d/%d hit %d-case cap]",
                                n_trunc, length(rs), CHECK_FINAL_SIZE),
               n_reps = length(rs),
               median = round(stats::median(vals, na.rm = TRUE)),
               q25    = round(safe_q(vals, 0.25)),
               q75    = round(safe_q(vals, 0.75)),
               q025   = round(safe_q(vals, 0.025)),
               q975   = round(safe_q(vals, 0.975)),
               stringsAsFactors = FALSE)
  }))
}))

summary_tbl          <- rbind(summary_rows, final_rows)
summary_tbl          <- summary_tbl[order(summary_tbl$r0, summary_tbl$scenario), ]
summary_tbl$scenario <- factor(summary_tbl$scenario, levels = scen_names)
write.csv(summary_tbl, file.path(out_dir, "05_summary_table.csv"), row.names = FALSE)
cat("\n--- Summary table ---\n")
print(summary_tbl, row.names = FALSE)

# ============================================================================
# 17. Bundle and save all results
# ============================================================================
# Saves a single .rds with all results + configuration for downstream use
# (e.g. by a paper-figure script). The config block records every parameter
# so the run is fully reproducible.

saveRDS(list(
  results_flat       = results_flat,        # raw per-replicate output
  cumulative_summary = cumulative_summary,   # aggregated cumulative incidence
  summary_tbl        = summary_tbl,          # snapshot table (Sep 1, Dec 31)
  rt_profiles        = rt_profiles,          # analytic Rt curves
  config = list(
    r0_grid             = R0_GRID,
    scenario_names      = scen_names,
    scalar_overrides    = SCALAR_OVERRIDES,
    npi_specs           = NPI_SPECS,
    funeral_frac        = FUNERAL_FRAC,
    seeding_cases       = SEEDING_CASES,
    n_stoch             = N_STOCH,
    takeoff_n           = TAKEOFF_N,
    max_retries         = MAX_RETRIES,
    check_final_size    = CHECK_FINAL_SIZE,
    epidemic_start_date = epi_start,
    pheic_dates         = PHEIC_DATES,
    summary_dates       = SUMMARY_DATES,
    seed_base           = SEED_BASE,
    n_workers           = N_WORKERS,
    run_date            = Sys.Date()   # record when the run was executed
  )
), file.path(out_dir, "05_projections.rds"))

message(sprintf(
  "\n05_projection_5_qcurves.R complete.\n  %d R0 x %d scenarios x %d reps = %d runs.\n  All outputs saved to: %s",
  n_r0, n_scen, N_STOCH, total, out_dir))
