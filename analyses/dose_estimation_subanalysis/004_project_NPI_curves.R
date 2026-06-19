# ============================================================================
# 004_project_NPI_curves.R
# ----------------------------------------------------------------------------
# WHAT THIS SCRIPT DOES
#   Runs FIBER stochastic epidemic simulations for ALL SIX forward-extrapolation
#   NPI scenarios produced by 01_fit_dose_q_curve.R, at a single fixed R0 value
#   (default: 1.65), with 100 stochastic replicates per scenario (600 total).
#
#   The six NPI scenarios are:
#     1. linear_to_90      -- dose coverage ramps linearly to 90% by day 100
#     2. linear_to_95      -- dose coverage ramps linearly to 95% by day 100
#     3. linear_to_95_dec  -- slow ramp to 95% reached only by end of Dec 2026
#     4. logistic          -- the fitted logistic projection (business as usual)
#     5. flat              -- held flat at the last observed Q value (~63%)
#     6. conflict          -- flat, then a conflict episode drops Q to ~20%, reverts
#
#   NPI parameters are carried over from 02_npi_inputs_and_fiber_runs.R unchanged.
#
#   Q curves only start from ~May 18 2026. The epidemic begins Feb 27 2026, so
#   the pre-intervention gap (Feb 27 -> May 18) is padded with Q = 0.
#
# Inputs  : outputs/dose_q_curve_extrapolation_scenarios.rds  (from 01)
#           data-processed/onsets_daily_incidence.csv          (optional overlay)
#           data-processed/insp_sitrep__national_cumulative_confirmed_cases__daily.csv
#                                                              (optional overlay)
# Outputs : outputs/004_scenario_matrices/          (per-scenario NPI CSVs)
#           outputs/004_npi_scenario_inputs.png
#           outputs/004_npi_scenario_cumulative_all.png
#           outputs/004_npi_scenario_cumulative_bands.png
#           outputs/004_npi_scenario_summary_table.csv
#           outputs/004_npi_scenario_projections.rds
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

# ----------------------------------------------------------------------------
# 1. Configuration  <-- the knobs to change
# ----------------------------------------------------------------------------

# Fixed R0 for all scenarios. Change this one number to explore a different R0.
R0_TARGET <- 1.65

# NPI parameters the Q curve drives (q0 = value at Q=0 [worst], q1 = at Q=1 [best]).
# Identical to 02 EXCEPT unsafe_funeral_prop_hosp: hospital deaths are managed
# safely, so only ~3% result in community funerals at Q=0, declining to ~1% at Q=1.
NPI_SPECS <- list(
  prob_hosp                = list(q0 = 0.00, q1 = 0.80),  # P(hospitalised | symptomatic)
  delay_hosp               = list(q0 = 6.00, q1 = 1.50),  # onset->hosp delay (days); improves DOWN
  prop_etu                 = list(q0 = 0.00, q1 = 0.90),  # fraction of hosp cases in ETU
  safe_funeral_prop        = list(q0 = 0.00, q1 = 0.90),  # community: fraction of funerals that are safe
  unsafe_funeral_prop_hosp = list(q0 = 0.95, q1 = 1.00),  # hospital (non-ETU): unsafe funeral prob (direct; matches 02)
  ppe_coverage             = list(q0 = 0.00, q1 = 0.90)   # PPE coverage lever -> ppe_coverage_hcw
)
UNSAFE_FUNERAL_ETU <- 0.0  # ETU deaths managed safely

# Fixed efficacy scalars (identical to 02).
SCALAR_OVERRIDES <- list(
  etu_efficacy                         = 0.60,
  general_hospital_quarantine_efficacy = 0.20,
  ppe_efficacy                         = 0.60,
  safe_funeral_efficacy                = 0.70
)

# Simulation controls.
FUNERAL_FRAC     <- 0.25    # share of t=0 transmission via funerals
SEEDING_CASES    <- 3L
N_STOCH          <- 1L      # replicates per scenario (set to 100 for full run)

TAKEOFF_N             <- 250L
TAKEOFF_DEADLINE_DATE <- as.Date("2026-06-15")
MAX_RETRIES           <- 50L
CHECK_FINAL_SIZE      <- 40000L  # cap each FIBER run at this many infections

# Epidemic start date: Feb 27 is day 0. Q curves start May 18; the gap is Q=0.
EPIDEMIC_START_DATE <- as.Date("2026-02-27")

# Summary dates for the output table.
SUMMARY_DATES <- as.Date(c("2026-09-01", "2026-12-31"))

# Timepoints at which to read cumulative cases per replicate (10-day grid + extras).
TIMEPOINTS <- sort(unique(c(seq(10L, 730L, by = 10L), 730L)))

# PHEIC-related dates for vertical lines on every calendar plot.
PHEIC_DATES <- as.Date(c("2026-05-14", "2026-05-18"))

# R0 grid for the analytic Rt profiles (section 6). Independent of R0_TARGET
# (which drives the actual FIBER runs). Change the range/step freely here.
R0_RT_GRID <- seq(1.45, 1.65, by = 0.05)

# Scenario matrix horizon (days since epidemic start).
MATRIX_HORIZON <- 730L

# Parallel workers.
N_WORKERS <- min(future::availableCores() - 1L, 50L)
SEED_BASE  <- 20260618L

if (CHECK_FINAL_SIZE < max(TIMEPOINTS))
  warning("CHECK_FINAL_SIZE < max(TIMEPOINTS): some cumulative readings will be truncated.",
          call. = FALSE)

# ----------------------------------------------------------------------------
# 2. Load the four Q-curve scenarios from 01
# ----------------------------------------------------------------------------
scen_path <- file.path(DIR_OUT, "dose_q_curve_extrapolation_scenarios.rds")
if (!file.exists(scen_path))
  stop("Scenarios file not found:\n  ", scen_path,
       "\nRun 01_fit_dose_q_curve.R first.", call. = FALSE)

all_scen      <- readRDS(scen_path)
all_scen$date <- as.Date(all_scen$date)

# Preserve the canonical scenario order set by script 01.
scen_names <- if (is.factor(all_scen$scenario)) levels(all_scen$scenario) else as.character(unique(all_scen$scenario))
expected <- c("linear_to_90", "linear_to_95", "linear_to_95_dec",
              "logistic", "flat", "conflict")
missing  <- setdiff(expected, scen_names)
if (length(missing) > 0L)
  stop("Expected scenario(s) not found in the file: ", paste(missing, collapse = ", "),
       "\nAvailable: ", paste(scen_names, collapse = ", "), call. = FALSE)
scen_names <- expected   # enforce canonical order

message("Loaded scenarios: ", paste(scen_names, collapse = ", "))

q_first_date <- min(all_scen$date)
epi_start    <- as.Date(EPIDEMIC_START_DATE)
if (epi_start > q_first_date)
  stop("EPIDEMIC_START_DATE (", epi_start, ") must be <= Q curve start (", q_first_date, ").",
       call. = FALSE)
offset_days <- as.integer(q_first_date - epi_start)
message(sprintf("Epidemic starts %s; Q curve starts %s -> %d day(s) of Q=0 prepended.",
                as.character(epi_start), as.character(q_first_date), offset_days))

day_to_date <- function(d) epi_start + d

pheic_vlines <- lapply(PHEIC_DATES, function(d)
  geom_vline(xintercept = d, linetype = "dashed", colour = "grey40"))

takeoff_day <- as.integer(TAKEOFF_DEADLINE_DATE - epi_start)
if (takeoff_day <= 0L)
  stop("TAKEOFF_DEADLINE_DATE must be after EPIDEMIC_START_DATE.", call. = FALSE)
message(sprintf("Takeoff: >= %d infections by %s (relative day %d).",
                TAKEOFF_N, as.character(TAKEOFF_DEADLINE_DATE), takeoff_day))

snapshot_days <- as.integer(SUMMARY_DATES - epi_start)
TIMEPOINTS    <- sort(unique(c(TIMEPOINTS, snapshot_days[snapshot_days > 0L])))

# ----------------------------------------------------------------------------
# 3. Build NPI scenario matrices (one per scenario)
# ----------------------------------------------------------------------------
lin         <- function(spec, q) spec$q0 + (spec$q1 - spec$q0) * q
matrix_days <- 0:MATRIX_HORIZON

build_scenario_matrix_df <- function(sn) {
  sel    <- all_scen[as.character(all_scen$scenario) == sn, , drop = FALSE]
  q_sim  <- as.integer(sel$date - epi_start)
  q_vals <- approx(q_sim, sel$q, xout = matrix_days, rule = 2)$y
  q_vals[matrix_days < offset_days] <- 0
  q_vals <- clip01(q_vals)

  safe_fp   <- clip01(lin(NPI_SPECS$safe_funeral_prop, q_vals))
  unsafe_fh <- clip01(lin(NPI_SPECS$unsafe_funeral_prop_hosp, q_vals))

  data.frame(
    scenario                 = sn,
    scenario_label           = sn,
    relative_day             = matrix_days,
    prob_hosp                = clip01(lin(NPI_SPECS$prob_hosp, q_vals)),
    delay_hosp               = pmax(lin(NPI_SPECS$delay_hosp, q_vals), 0.01),
    prob_unsafe_funeral_comm = clip01(1 - safe_fp),
    prob_unsafe_funeral_hosp = unsafe_fh,
    prob_unsafe_funeral_etu  = clip01(rep(UNSAFE_FUNERAL_ETU, length(matrix_days))),
    prop_etu                 = clip01(lin(NPI_SPECS$prop_etu, q_vals)),
    ipc_helper               = clip01(lin(NPI_SPECS$ppe_coverage, q_vals)),
    q_value                  = q_vals,
    stringsAsFactors = FALSE
  )
}

mat_dir <- file.path(DIR_OUT, "004_scenario_matrices")
dir.create(mat_dir, recursive = TRUE, showWarnings = FALSE)

scenario_matrix_dfs  <- setNames(lapply(scen_names, build_scenario_matrix_df), scen_names)
scenario_matrix_csvs <- setNames(
  vapply(scen_names, function(sn) {
    p <- file.path(mat_dir, sprintf("dose_npi_matrix_%s.csv", sn))
    write.csv(scenario_matrix_dfs[[sn]], p, row.names = FALSE)
    p
  }, character(1)),
  scen_names
)
message("Wrote NPI scenario matrices to ", mat_dir)

# ----------------------------------------------------------------------------
# 4. Plot NPI input curves for all scenarios
# ----------------------------------------------------------------------------
scenario_colours <- c(
  linear_to_90     = "#1b9e77",
  linear_to_95     = "#66a61e",
  linear_to_95_dec = "#e6ab02",
  logistic         = "#7570b3",
  flat             = "#d95f02",
  conflict         = "#e7298a"
)
scenario_labels <- c(
  linear_to_90     = "1. Linear to 90% by day 100",
  linear_to_95     = "2. Linear to 95% by day 100",
  linear_to_95_dec = "3. Linear to 95% by Dec 2026",
  logistic         = "4. Logistic projection (current)",
  flat             = "5. Flat at last value",
  conflict         = "6. Conflict at day 100"
)

npi_long <- do.call(rbind, lapply(scen_names, function(sn) {
  scenario_matrix_dfs[[sn]] %>%
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
  scale_colour_manual(values = scenario_colours, labels = scenario_labels,
                      name = "NPI scenario") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
  labs(
    title    = "Time-varying NPI inputs by scenario",
    subtitle = sprintf("Dashed verticals = PHEIC dates (%s, %s); R0 = %.2f",
                       format(PHEIC_DATES[1], "%d %b"), format(PHEIC_DATES[2], "%d %b"),
                       R0_TARGET),
    x = "Date", y = "Input value"
  ) +
  theme_bw(base_size = 10) +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 6),
        legend.position = "bottom")
ggsave(file.path(DIR_OUT, "004_npi_scenario_inputs.png"), p_inputs,
       width = 11, height = 8, dpi = 150)
print(p_inputs)

# ----------------------------------------------------------------------------
# 5. Assemble fiber args per scenario at R0_TARGET
# ----------------------------------------------------------------------------
# R0 invariants depend only on the fixed scalar overrides, not the time-varying
# trajectory, so we compute them once from an arbitrary scenario's matrix.
ref_sm <- read_scenario_matrix(scenario_matrix_csvs[[scen_names[1]]])
mp_ref <- make_model_parameters(
  scenario_id     = scen_names[1],
  scenario_matrix = ref_sm,
  overrides       = modifyList(
    list(seeding_cases = SEEDING_CASES, check_final_size = CHECK_FINAL_SIZE),
    SCALAR_OVERRIDES
  )
)

inv   <- compute_R0_invariants(mp_ref$args, n = 50000L, seed = 42L)
D     <- D_from_invariants(inv, mp_ref$args$etu_efficacy,
                           mp_ref$args$general_hospital_quarantine_efficacy)
F_fun <- F_from_invariants(inv, mp_ref$args$safe_funeral_efficacy)
means <- solve_offspring_means(R0_TARGET, FUNERAL_FRAC, D, F_fun)
message(sprintf("R0 = %.2f: mn_genPop = %.4f, mn_funeral = %.4f",
                R0_TARGET, means$mn_genPop, means$mn_funeral))

args_by_scenario <- setNames(
  lapply(scen_names, function(sn) {
    sm <- read_scenario_matrix(scenario_matrix_csvs[[sn]])
    mp <- make_model_parameters(
      scenario_id     = sn,
      scenario_matrix = sm,
      overrides       = modifyList(
        list(seeding_cases = SEEDING_CASES, check_final_size = CHECK_FINAL_SIZE),
        SCALAR_OVERRIDES
      )
    )
    a <- mp$args
    a$mn_offspring_genPop  <- means$mn_genPop
    a$mn_offspring_funeral <- means$mn_funeral
    a$seeding_cases        <- SEEDING_CASES
    a$check_final_size     <- CHECK_FINAL_SIZE
    a$seed                 <- NULL
    a
  }),
  scen_names
)

# ----------------------------------------------------------------------------
# 6. Analytic Rt profiles: R0_RT_GRID x all 4 NPI scenarios
# ----------------------------------------------------------------------------
# For every combination of R0 (1.45 -> 1.65) and NPI scenario, compute the
# single-type instantaneous and case Rt over time. D and F_fun are reused from
# section 5 (they depend only on the fixed scalar overrides, not on R0 or the
# time-varying trajectory), so we only need to re-solve the offspring means per R0.
# Adapted from 02_npi_inputs_and_fiber_runs.R section 5.
RT_TIMES <- 0:730L
RT_MC_N  <- 50000L
RT_SEED  <- 1L

message(sprintf("Computing analytic Rt profiles: %d R0 values x %d scenarios...",
                length(R0_RT_GRID), length(scen_names)))

# Build a named list of args for every (R0, scenario) pair.
rt_args_grid <- do.call(c, lapply(R0_RT_GRID, function(r0) {
  m <- solve_offspring_means(r0, FUNERAL_FRAC, D, F_fun)
  setNames(lapply(scen_names, function(sn) {
    a <- args_by_scenario[[sn]]
    a$mn_offspring_genPop  <- m$mn_genPop
    a$mn_offspring_funeral <- m$mn_funeral
    a
  }), paste(sprintf("R0_%.2f", r0), scen_names, sep = "__"))
}))

rt_profiles <- do.call(rbind, lapply(R0_RT_GRID, function(r0) {
  do.call(rbind, lapply(scen_names, function(sn) {
    key <- paste(sprintf("R0_%.2f", r0), sn, sep = "__")
    rt  <- Rt_curve_single_type(rt_args_grid[[key]], times = RT_TIMES,
                                n = RT_MC_N, seed = RT_SEED)
    rt$r0       <- r0
    rt$scenario <- sn
    rt
  }))
}))
write.csv(rt_profiles, file.path(DIR_OUT, "004_npi_scenario_rt_profiles.csv"),
          row.names = FALSE)

rt_long <- rbind(
  data.frame(r0 = rt_profiles$r0, scenario = rt_profiles$scenario,
             day = rt_profiles$time, mode = "instantaneous", Rt = rt_profiles$R_inst),
  data.frame(r0 = rt_profiles$r0, scenario = rt_profiles$scenario,
             day = rt_profiles$time, mode = "case",          Rt = rt_profiles$R_case)
)
rt_long$mode     <- factor(rt_long$mode, levels = c("instantaneous", "case"))
rt_long$scenario <- factor(rt_long$scenario, levels = scen_names)
rt_long$date     <- day_to_date(rt_long$day)

p_rt <- ggplot(rt_long, aes(date, Rt, colour = scenario, linetype = mode,
                            group = interaction(r0, scenario, mode))) +
  geom_hline(yintercept = 1, colour = "grey40", linewidth = 0.5) +
  pheic_vlines +
  geom_line(linewidth = 0.75) +
  facet_wrap(~ r0, ncol = 3,
             labeller = label_bquote(R[0] == .(r0))) +
  scale_colour_manual(values = scenario_colours, labels = scenario_labels,
                      name = "Scenario") +
  scale_linetype_manual(values = c(instantaneous = "solid", case = "22"),
                        name = "Rt type") +
  scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") +
  labs(
    title    = "Analytic Rt by R0 and NPI scenario",
    subtitle = sprintf("Solid = instantaneous Rt; dashed = case Rt; Rt < 1 = controlled\nDashed verticals = PHEIC dates (%s, %s)",
                       format(PHEIC_DATES[1], "%d %b"), format(PHEIC_DATES[2], "%d %b")),
    x = "Date", y = expression(R[t])
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        legend.position = "bottom")
ggsave(file.path(DIR_OUT, "004_npi_scenario_rt_profiles.png"), p_rt,
       width = 11, height = 8, dpi = 150)
print(p_rt)

# ----------------------------------------------------------------------------
# 7. Per-replicate runner + summariser (identical pattern to 02)
# ----------------------------------------------------------------------------
run_one_takeoff <- function(args, base_seed, takeoff_n, max_retries, takeoff_day) {
  seed <- base_seed; retry <- 0L; inf_times <- numeric(0)
  n_cases <- 0L; n_by_deadline <- 0L
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

AMOUNTS <- c(100L, 500L, 1000L, 5000L, 10000L, 20000L, 40000L)

summarise_run <- function(inf_times, timepoints, amounts) {
  st <- sort(inf_times); n <- length(st)
  list(
    cum_at  = vapply(timepoints, function(t) sum(st <= t), numeric(1)),
    time_to = vapply(amounts,    function(a) if (n >= a) st[a] else NA_real_, numeric(1)),
    n_cases = n
  )
}

# ----------------------------------------------------------------------------
# 8. Run FIBER: all 400 runs in one parallel batch
# ----------------------------------------------------------------------------
# All (scenario, replicate) combinations are submitted at once so the scheduler
# can pack every available core without the sequential scenario-loop bottleneck.
# On Linux/macOS multicore (fork-based, lower overhead) is used automatically;
# on Windows we fall back to multisession (socket-based).
n_scen <- length(scen_names)
total  <- n_scen * N_STOCH

run_grid <- do.call(rbind, lapply(seq_along(scen_names), function(si)
  data.frame(si = si, sn = scen_names[si], rep_id = seq_len(N_STOCH),
             stringsAsFactors = FALSE)))

if (future::supportsMulticore()) {
  plan(multicore, workers = N_WORKERS)
  message(sprintf("Using multicore (fork) on %d workers.", N_WORKERS))
} else {
  plan(multisession, workers = N_WORKERS)
  message(sprintf("Using multisession on %d workers.", N_WORKERS))
}

message(sprintf("Running %d scenarios x %d reps = %d total, R0 = %.2f...",
                n_scen, N_STOCH, total, R0_TARGET))
t_start <- proc.time()

results_flat <- future_lapply(seq_len(nrow(run_grid)), function(idx) {
  si     <- run_grid$si[idx]
  sn     <- run_grid$sn[idx]
  rep_id <- run_grid$rep_id[idx]
  base_seed <- SEED_BASE +
    (si     - 1L) * N_STOCH * (MAX_RETRIES + 1L) +
    (rep_id - 1L) *           (MAX_RETRIES + 1L)
  run <- run_one_takeoff(args_by_scenario[[si]], base_seed, TAKEOFF_N, MAX_RETRIES, takeoff_day)
  s   <- summarise_run(run$inf_times, TIMEPOINTS, AMOUNTS)
  list(scenario = sn, scen_idx = si, rep_id = rep_id,
       n_cases = run$n_cases, n_by_deadline = run$n_by_deadline,
       took_off = run$took_off, retries = run$retries,
       cum_at = s$cum_at, time_to = s$time_to)
},
future.globals = list(
  run_grid = run_grid, scen_names = scen_names,
  args_by_scenario = args_by_scenario,
  N_STOCH = N_STOCH, MAX_RETRIES = MAX_RETRIES, SEED_BASE = SEED_BASE,
  TAKEOFF_N = TAKEOFF_N, takeoff_day = takeoff_day,
  TIMEPOINTS = TIMEPOINTS, AMOUNTS = AMOUNTS,
  run_one_takeoff = run_one_takeoff, summarise_run = summarise_run
),
future.packages = "fiber", future.seed = TRUE)

plan(sequential)
elapsed <- proc.time() - t_start
message(sprintf("Done: %d replicates in %.1f min.", length(results_flat), elapsed["elapsed"] / 60))
saveRDS(results_flat, file.path(DIR_OUT, "004_npi_scenario_projections.rds"))

# ----------------------------------------------------------------------------
# 8. Summarise across replicates per scenario
# ----------------------------------------------------------------------------
safe_median <- function(x) if (all(is.na(x))) NA_real_ else stats::median(x, na.rm = TRUE)
safe_q <- function(x, p)   if (all(is.na(x))) NA_real_ else
  as.numeric(stats::quantile(x, p, na.rm = TRUE, names = FALSE))

took   <- Filter(function(r) isTRUE(r$took_off), results_flat)
n_fail <- length(results_flat) - length(took)
if (n_fail > 0L)
  message(sprintf("Note: %d/%d replicates failed to reach takeoff (excluded from summaries).",
                  n_fail, length(results_flat)))
if (length(took) == 0L)
  stop("No replicate reached takeoff. Check TAKEOFF_N, SEEDING_CASES, TAKEOFF_DEADLINE_DATE.",
       call. = FALSE)

by_scen <- split(took, vapply(took, function(r) r$scen_idx, integer(1)))

cumulative_by_scenario <- do.call(rbind, lapply(seq_along(scen_names), function(si) {
  rs <- by_scen[[as.character(si)]]
  if (is.null(rs)) return(NULL)
  M <- do.call(rbind, lapply(rs, function(r) r$cum_at))
  data.frame(
    scenario      = scen_names[si],
    timepoint_day = TIMEPOINTS,
    n_runs        = nrow(M),
    mean          = apply(M, 2, mean),
    median        = apply(M, 2, stats::median),
    q10           = apply(M, 2, safe_q, 0.10),
    q25           = apply(M, 2, safe_q, 0.25),
    q75           = apply(M, 2, safe_q, 0.75),
    q90           = apply(M, 2, safe_q, 0.90),
    stringsAsFactors = FALSE
  )
}))
cumulative_by_scenario$date <- day_to_date(cumulative_by_scenario$timepoint_day)
cumulative_by_scenario$scenario <- factor(cumulative_by_scenario$scenario, levels = scen_names)

write.csv(cumulative_by_scenario,
          file.path(DIR_OUT, "004_npi_scenario_cumulative.csv"), row.names = FALSE)

# ----------------------------------------------------------------------------
# 9. Optional observed-data overlays (reused from 02)
# ----------------------------------------------------------------------------
onset_overlay <- confirmed_overlay <- NULL
onset_csv <- here("data-processed", "onsets_daily_incidence.csv")
if (file.exists(onset_csv)) {
  onsets_obs      <- read.csv(onset_csv, stringsAsFactors = FALSE)
  onsets_obs$date <- as.Date(onsets_obs$date)
  onset_overlay   <- geom_line(data = onsets_obs, aes(date, cumulative_onsets),
                               inherit.aes = FALSE, colour = "#d62728", linewidth = 1.1)
  message(sprintf("Overlaying observed cumulative onsets (%s to %s, max %.0f).",
                  as.character(min(onsets_obs$date)), as.character(max(onsets_obs$date)),
                  max(onsets_obs$cumulative_onsets, na.rm = TRUE)))
}
confirmed_csv <- here("data-processed",
                      "insp_sitrep__national_cumulative_confirmed_cases__daily.csv")
if (file.exists(confirmed_csv)) {
  confirmed_obs      <- read.csv(confirmed_csv, stringsAsFactors = FALSE)
  confirmed_obs$date <- as.Date(confirmed_obs$date, format = "%d/%m/%Y")
  confirmed_overlay  <- geom_line(data = confirmed_obs,
                                  aes(date, national_cumulative_confirmed_cases),
                                  inherit.aes = FALSE, colour = "#1a9850", linewidth = 1.1)
  message(sprintf("Overlaying observed cumulative confirmed cases (%s to %s, max %.0f).",
                  as.character(min(confirmed_obs$date)), as.character(max(confirmed_obs$date)),
                  max(confirmed_obs$national_cumulative_confirmed_cases, na.rm = TRUE)))
}

# ----------------------------------------------------------------------------
# 10. Plot: individual stochastic trajectories for all 4 scenarios
# ----------------------------------------------------------------------------
traj_long <- do.call(rbind, lapply(took, function(r)
  data.frame(scenario = r$scenario, rep_id = r$rep_id,
             day = TIMEPOINTS, cum = r$cum_at, stringsAsFactors = FALSE)))
traj_long$date     <- day_to_date(traj_long$day)
traj_long$scenario <- factor(traj_long$scenario, levels = scen_names)

traj_alpha <- max(0.04, min(0.35, 20 / N_STOCH))

p_all_traj <- ggplot(traj_long,
                     aes(date, cum, colour = scenario,
                         group = interaction(scenario, rep_id))) +
  pheic_vlines +
  geom_line(alpha = traj_alpha, linewidth = 0.3) +
  geom_line(data = cumulative_by_scenario,
            aes(date, median, colour = scenario, group = scenario),
            inherit.aes = FALSE, linewidth = 1.1) +
  onset_overlay +
  confirmed_overlay +
  scale_colour_manual(values = scenario_colours, labels = scenario_labels,
                      name = "Scenario") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
  labs(
    title    = sprintf("Cumulative incidence: all %d NPI scenarios (R0 = %.2f, n = %d reps each)",
                       length(scen_names), R0_TARGET, N_STOCH),
    subtitle = sprintf("Thin lines = stochastic replicates; bold = median; red = onsets; green = confirmed\nDashed verticals = PHEIC dates (%s, %s)",
                       format(PHEIC_DATES[1], "%d %b"), format(PHEIC_DATES[2], "%d %b")),
    x = "Date", y = "Cumulative cases"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")
ggsave(file.path(DIR_OUT, "004_npi_scenario_cumulative_all.png"), p_all_traj,
       width = 11, height = 6.5, dpi = 150)
print(p_all_traj)

# ----------------------------------------------------------------------------
# 11. Plot: median + 10-90% and 25-75% uncertainty bands, all scenarios
# ----------------------------------------------------------------------------
p_bands <- ggplot(cumulative_by_scenario,
                  aes(date, colour = scenario, fill = scenario)) +
  pheic_vlines +
  geom_ribbon(aes(ymin = q10, ymax = q90), alpha = 0.10, colour = NA) +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.18, colour = NA) +
  geom_line(aes(y = median), linewidth = 1.0) +
  onset_overlay +
  confirmed_overlay +
  scale_colour_manual(values = scenario_colours, labels = scenario_labels,
                      name = "Scenario") +
  scale_fill_manual(values   = scenario_colours, labels = scenario_labels,
                    name = "Scenario") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
  labs(
    title    = sprintf("Cumulative incidence by NPI scenario: median + uncertainty (R0 = %.2f)",
                       R0_TARGET),
    subtitle = sprintf("Solid = median; inner band = 25-75%%; outer band = 10-90%%\nDashed verticals = PHEIC dates (%s, %s)",
                       format(PHEIC_DATES[1], "%d %b"), format(PHEIC_DATES[2], "%d %b")),
    x = "Date", y = "Cumulative cases"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")
ggsave(file.path(DIR_OUT, "004_npi_scenario_cumulative_bands.png"), p_bands,
       width = 11, height = 6.5, dpi = 150)
print(p_bands)

# ----------------------------------------------------------------------------
# 12. Summary table: cumulative cases at key dates + total epidemic size
# ----------------------------------------------------------------------------
summary_date_rows <- do.call(rbind, lapply(seq_along(scen_names), function(si) {
  sn <- scen_names[si]
  rs <- by_scen[[as.character(si)]]
  if (is.null(rs)) return(NULL)
  do.call(rbind, lapply(seq_along(SUMMARY_DATES), function(di) {
    sd  <- SUMMARY_DATES[di]
    day <- as.integer(sd - epi_start)
    idx <- match(day, TIMEPOINTS)
    if (is.na(idx)) idx <- which.min(abs(TIMEPOINTS - day))
    vals <- vapply(rs, function(r) r$cum_at[idx], numeric(1))
    data.frame(
      scenario = sn,
      date     = as.character(sd),
      n_reps   = length(vals),
      median   = round(safe_median(vals)),
      q10      = round(safe_q(vals, 0.10)),
      q25      = round(safe_q(vals, 0.25)),
      q75      = round(safe_q(vals, 0.75)),
      q90      = round(safe_q(vals, 0.90)),
      stringsAsFactors = FALSE
    )
  }))
}))

final_size_rows <- do.call(rbind, lapply(seq_along(scen_names), function(si) {
  sn  <- scen_names[si]
  rs  <- by_scen[[as.character(si)]]
  if (is.null(rs)) return(NULL)
  vals        <- vapply(rs, function(r) as.numeric(r$n_cases), numeric(1))
  n_truncated <- sum(vals >= CHECK_FINAL_SIZE)
  data.frame(
    scenario = sn,
    date     = sprintf("total_final_size [%d/%d reps hit %d-case cap]",
                       n_truncated, length(rs), CHECK_FINAL_SIZE),
    n_reps   = length(vals),
    median   = round(safe_median(vals)),
    q10      = round(safe_q(vals, 0.10)),
    q25      = round(safe_q(vals, 0.25)),
    q75      = round(safe_q(vals, 0.75)),
    q90      = round(safe_q(vals, 0.90)),
    stringsAsFactors = FALSE
  )
}))

summary_tbl          <- rbind(summary_date_rows, final_size_rows)
summary_tbl$scenario <- factor(summary_tbl$scenario, levels = scen_names)
summary_tbl          <- summary_tbl[order(summary_tbl$scenario, summary_tbl$date), ]

write.csv(summary_tbl, file.path(DIR_OUT, "004_npi_scenario_summary_table.csv"),
          row.names = FALSE)
cat("\n--- Summary: cumulative cases by scenario and date ---\n")
print(summary_tbl, row.names = FALSE)

# ----------------------------------------------------------------------------
# 13. Bundle and save everything
# ----------------------------------------------------------------------------
saveRDS(list(
  results_flat       = results_flat,
  cumulative_by_scen = cumulative_by_scenario,
  summary_tbl        = summary_tbl,
  config = list(
    r0_target             = R0_TARGET,
    scenario_names        = scen_names,
    npi_specs             = NPI_SPECS,
    unsafe_funeral_etu    = UNSAFE_FUNERAL_ETU,
    scalar_overrides      = SCALAR_OVERRIDES,
    funeral_frac          = FUNERAL_FRAC,
    seeding_cases         = SEEDING_CASES,
    n_stoch               = N_STOCH,
    takeoff_n             = TAKEOFF_N,
    max_retries           = MAX_RETRIES,
    takeoff_deadline_date = TAKEOFF_DEADLINE_DATE,
    takeoff_day           = takeoff_day,
    check_final_size      = CHECK_FINAL_SIZE,
    timepoints            = TIMEPOINTS,
    amounts               = AMOUNTS,
    epidemic_start_date   = epi_start,
    q_first_date          = q_first_date,
    offset_days           = offset_days,
    pheic_dates           = PHEIC_DATES,
    summary_dates         = SUMMARY_DATES
  )
), file.path(DIR_OUT, "004_npi_scenario_projections.rds"))

message("\n004_project_NPI_curves.R complete. Results in outputs/.")
