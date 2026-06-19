# ============================================================================
# 05b_projection_logistic_good.R
# ----------------------------------------------------------------------------
# Companion to 05_projection_6_qcurves.R that runs the SAME pipeline but for the
# SINGLE new scenario `logistic_projection_good` only. It writes a separate
# results bundle in the IDENTICAL format to 05's, so the new scenario can be
# appended to an existing 6-scenario 05 run WITHOUT re-running those 6 scenarios.
#
# Everything that affects the numbers (R0 grid, efficacies, NPI specs, takeoff
# rules, CHECK_FINAL_SIZE, TIMEPOINTS, AMOUNTS, R0-invariant calibration) is
# copied verbatim from 05 so the results are directly comparable / mergeable.
# Plotting is omitted — regenerate the combined figures from the merged results.
#
# PREREQUISITE
#   Re-run 01_fit_dose_q_curve.R first so the scenarios .rds contains the
#   `logistic_projection_good` curve.
#
# OUTPUTS (to outputs/05_outputs/, alongside the 05 outputs)
#   05b_logistic_good_projections_raw.rds   raw results_flat list (clean append)
#   05b_logistic_good_projections.rds        bundle: results_flat + summaries + config
#
# ----------------------------------------------------------------------------
# HOW TO MERGE WITH THE MAIN 6-SCENARIO 05 RUN
# ----------------------------------------------------------------------------
#   main <- readRDS(file.path(out_dir, "05_projections.rds"))
#   good <- readRDS(file.path(out_dir, "05b_logistic_good_projections.rds"))
#
#   # 1) RAW per-replicate results -- scenario is plain character here, so this
#   #    is a clean concatenation (this is the statistically correct merge: pool
#   #    the raw replicates, then re-run 05's aggregation/plots on the result):
#   results_flat <- c(main$results_flat, good$results_flat)
#
#   # 2) Pre-computed summary tables -- coerce the scenario column on BOTH sides
#   #    to character before rbind (the main run saved it as a 6-level factor),
#   #    then optionally re-factor with the full 7-scenario ordering:
#   fix <- function(df) { df$scenario <- as.character(df$scenario); df }
#   cumulative_summary <- rbind(fix(main$cumulative_summary), fix(good$cumulative_summary))
#   summary_tbl        <- rbind(fix(main$summary_tbl),        fix(good$summary_tbl))
#   rt_profiles        <- rbind(fix(main$rt_profiles),        fix(good$rt_profiles))
#   lev <- c("linear_to_90","linear_to_95","linear_to_95_dec",
#            "logistic","flat","conflict","logistic_projection_good")
#   for (nm in c("scenario"))
#     cumulative_summary[[nm]] <- factor(cumulative_summary[[nm]], levels = lev)
# ============================================================================

suppressPackageStartupMessages({
  library(here)
  library(future)
  library(future.apply)
  library(fiber)
})

source(here("analyses", "dose_estimation_subanalysis", "helpers.R"))
source(here("functions", "setup_model_parameters.R"))
source(here("functions", "calculate_model_approx_r0.R"))
source(here("functions", "calculate_model_approx_rt.R"))

set.seed(123)

# ============================================================================
# 0. Which scenario to run
# ============================================================================
SCENARIO_TO_RUN <- "logistic_projection_good"

# ============================================================================
# 1. Configuration  <-- MUST match 05_projection_6_qcurves.R for mergeability
# ============================================================================
R0_GRID <- c(1.50, 1.55, 1.60)

N_STOCH <- 330  # stochastic replicates per (R0, scenario)

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

FUNERAL_FRAC     <- 0.25
SEEDING_CASES    <- 5L
CHECK_FINAL_SIZE <- 75000L

TAKEOFF_N             <- 250L
TAKEOFF_DEADLINE_DATE <- as.Date("2026-06-15")
MAX_RETRIES           <- 50L

EPIDEMIC_START_DATE <- as.Date("2026-02-17")
SUMMARY_DATES       <- as.Date(c("2026-09-01", "2026-12-31"))
MATRIX_HORIZON      <- 730L
TIMEPOINTS          <- 0:730L
AMOUNTS <- c(100L, 500L, 1000L, 5000L, 10000L, 25000L, 50000L, 60000L)

N_WORKERS <- min(future::availableCores() - 4L, 120)

# Distinct seed base (offset well clear of 05's SEED_BASE = 20260619) so this
# scenario's replicates are independent of the original 6-scenario run.
SEED_BASE_NEW <- 20260619L + 100000000L

# Analytic Rt settings (identical to 05).
RT_TIMES <- 0:365L
RT_MC_N  <- 5000L
RT_SEED  <- 1L

# ============================================================================
# 2. Setup: dates, output folder
# ============================================================================
out_dir <- file.path(DIR_OUT, "05_outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

epi_start   <- as.Date(EPIDEMIC_START_DATE)
day_to_date <- function(d) epi_start + d

# Replicate 05's TIMEPOINTS augmentation EXACTLY so cum_at vectors align.
snapshot_days <- as.integer(SUMMARY_DATES - epi_start)
TIMEPOINTS    <- sort(unique(c(TIMEPOINTS, snapshot_days[snapshot_days > 0L])))
takeoff_day   <- as.integer(TAKEOFF_DEADLINE_DATE - epi_start)

# ============================================================================
# 3. Load Q-curve scenarios; select the new one
# ============================================================================
scen_path <- file.path(DIR_OUT, "dose_q_curve_extrapolation_scenarios.rds")
if (!file.exists(scen_path))
  stop("Run 01_fit_dose_q_curve.R first to generate: ", scen_path, call. = FALSE)

all_scen      <- readRDS(scen_path)
all_scen$date <- as.Date(all_scen$date)
scen_full     <- if (is.factor(all_scen$scenario)) levels(all_scen$scenario) else
                   as.character(unique(all_scen$scenario))

if (!SCENARIO_TO_RUN %in% scen_full)
  stop("Scenario '", SCENARIO_TO_RUN, "' not found in ", scen_path,
       ".\n  Re-run the updated 01_fit_dose_q_curve.R to add it. Available: ",
       paste(scen_full, collapse = ", "), call. = FALSE)

# Index of this scenario in the full (combined) scenario ordering, recorded on
# each result for parity with 05's scen_idx (not used by the aggregation, which
# keys on the scenario NAME).
scen_idx_full <- match(SCENARIO_TO_RUN, scen_full)

q_first_date <- min(all_scen$date)
offset_days  <- as.integer(q_first_date - epi_start)

message(sprintf("Running scenario '%s' (index %d of %d) over R0 = %s, %d reps each.",
                SCENARIO_TO_RUN, scen_idx_full, length(scen_full),
                paste(sprintf("%.2f", R0_GRID), collapse = ", "), N_STOCH))

# ============================================================================
# 4. Build the NPI scenario matrix (05 section 4, in-memory, single scenario)
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
scen_mat <- build_matrix_df(SCENARIO_TO_RUN)

# ============================================================================
# 5. Assemble FIBER args per R0 (05 section 6)
# ============================================================================
# D and F_fun depend only on the fixed efficacies and the t=0 (Q=0) scenario
# inputs, which are identical across all scenarios, so the offspring means solved
# here match 05's for each R0 exactly.
overrides_base <- modifyList(
  list(seeding_cases = SEEDING_CASES, check_final_size = CHECK_FINAL_SIZE),
  SCALAR_OVERRIDES
)

mp_ref <- make_model_parameters(SCENARIO_TO_RUN, scen_mat, overrides = overrides_base)
inv    <- compute_R0_invariants(mp_ref$args, n = 50000L, seed = 42L)
D      <- D_from_invariants(inv, mp_ref$args$etu_efficacy,
                            mp_ref$args$general_hospital_quarantine_efficacy)
F_fun  <- F_from_invariants(inv, mp_ref$args$safe_funeral_efficacy)

args_grid <- list()
for (r0 in R0_GRID) {
  means <- solve_offspring_means(r0, FUNERAL_FRAC, D, F_fun)
  message(sprintf("R0 = %.2f: mn_genPop = %.4f, mn_funeral = %.4f",
                  r0, means$mn_genPop, means$mn_funeral))
  key <- paste(sprintf("R0_%.2f", r0), SCENARIO_TO_RUN, sep = "__")
  mp  <- make_model_parameters(SCENARIO_TO_RUN, scen_mat, overrides = overrides_base)
  a   <- mp$args
  a$mn_offspring_genPop  <- means$mn_genPop
  a$mn_offspring_funeral <- means$mn_funeral
  a$seed                 <- NULL
  args_grid[[key]] <- a
}

# ============================================================================
# 6. Analytic Rt profiles for this scenario (05 section 7)
# ============================================================================
rt_profiles <- do.call(rbind, lapply(R0_GRID, function(r0) {
  key <- paste(sprintf("R0_%.2f", r0), SCENARIO_TO_RUN, sep = "__")
  rt  <- Rt_curve_single_type(args_grid[[key]], times = RT_TIMES,
                              n = RT_MC_N, seed = RT_SEED)
  data.frame(r0 = r0, scenario = SCENARIO_TO_RUN, day = rt$time,
             R_inst = rt$R_inst, R_case = rt$R_case, stringsAsFactors = FALSE)
}))
rt_profiles$date     <- day_to_date(rt_profiles$day)
rt_profiles$r0_label <- sprintf("R0 = %.2f", rt_profiles$r0)

p_rt <- ggplot(rt_profiles, aes(date, R_inst, colour = scenario,
                                group = interaction(r0_label, scenario))) +
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
      "ETU=0.84, GenHosp=0.30, PPE=0.84, SafeFuneral=0.88\n",
      "Black dashed = Rt 1 (above = growing, below = controlled). ",
      "Grey dashed = PHEIC dates."
    ),
    x = "Date", y = expression(R[t] ~ "(instantaneous)")
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x   = element_text(angle = 45, hjust = 1, size = 7),
        legend.position = "bottom")

# ============================================================================
# 7. Runner helpers (copied verbatim from 05 section 8)
# ============================================================================
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

summarise_run <- function(inf_times, timepoints, amounts) {
  st <- sort(inf_times); n <- length(st)
  list(
    cum_at  = findInterval(timepoints, st),
    time_to = vapply(amounts, function(a) if (n >= a) st[a] else NA_real_, numeric(1)),
    n_cases = n
  )
}

# ============================================================================
# 8. Parallel FIBER runs: (R0 x replicate) for this one scenario (05 section 9)
# ============================================================================
n_r0  <- length(R0_GRID)
total <- n_r0 * N_STOCH

run_grid <- do.call(rbind, lapply(seq_along(R0_GRID), function(ri)
  data.frame(ri = ri, r0 = R0_GRID[ri], rep_id = seq_len(N_STOCH),
             stringsAsFactors = FALSE)))

if (future::supportsMulticore()) {
  plan(multicore, workers = N_WORKERS)
  message(sprintf("Using multicore (fork-based) on %d workers.", N_WORKERS))
} else {
  plan(multisession, workers = N_WORKERS)
  message(sprintf("Using multisession (socket-based) on %d workers.", N_WORKERS))
}

message(sprintf("Starting FIBER: %d R0 values x 1 scenario x %d reps = %d total runs.",
                n_r0, N_STOCH, total))
t_start <- proc.time()

results_flat <- future_lapply(seq_len(nrow(run_grid)), function(idx) {
  ri     <- run_grid$ri[idx]
  r0     <- run_grid$r0[idx]
  rep_id <- run_grid$rep_id[idx]
  key    <- paste(sprintf("R0_%.2f", r0), SCENARIO_TO_RUN, sep = "__")

  # Unique seed block for this (R0, replicate), disjoint from the 05 run.
  base_seed <- SEED_BASE_NEW +
    ((ri - 1L) * N_STOCH + (rep_id - 1L)) * (MAX_RETRIES + 1L)

  run <- run_one_takeoff(args_grid[[key]], base_seed, TAKEOFF_N, MAX_RETRIES, takeoff_day)
  s   <- summarise_run(run$inf_times, TIMEPOINTS, AMOUNTS)

  # Same field layout as 05's results_flat so the lists concatenate cleanly.
  list(r0 = r0, r0_idx = ri, scenario = SCENARIO_TO_RUN, scen_idx = scen_idx_full,
       rep_id = rep_id, n_cases = run$n_cases, n_by_deadline = run$n_by_deadline,
       took_off = run$took_off, retries = run$retries,
       cum_at = s$cum_at, time_to = s$time_to)
},
future.globals = list(
  run_grid = run_grid, args_grid = args_grid, SCENARIO_TO_RUN = SCENARIO_TO_RUN,
  scen_idx_full = scen_idx_full, N_STOCH = N_STOCH,
  MAX_RETRIES = MAX_RETRIES, SEED_BASE_NEW = SEED_BASE_NEW,
  TAKEOFF_N = TAKEOFF_N, takeoff_day = takeoff_day,
  TIMEPOINTS = TIMEPOINTS, AMOUNTS = AMOUNTS,
  run_one_takeoff = run_one_takeoff, summarise_run = summarise_run
),
future.packages = "fiber", future.seed = TRUE)

plan(sequential)
elapsed <- proc.time() - t_start
message(sprintf("FIBER complete: %d replicates in %.1f min.",
                length(results_flat), elapsed["elapsed"] / 60))

# --- CHECKPOINT SAVE (raw results, clean to c() onto 05's results_flat) ------
saveRDS(results_flat, file.path(out_dir, "05b_logistic_good_projections_raw.rds"))
message("Checkpoint: raw results saved to 05b_logistic_good_projections_raw.rds")

# ============================================================================
# 9. Aggregate replicate results (05 section 10) -- scenario kept as character
# ============================================================================
safe_q <- function(x, p) if (all(is.na(x))) NA_real_ else
  as.numeric(stats::quantile(x, p, na.rm = TRUE, names = FALSE))

took   <- Filter(function(r) isTRUE(r$took_off), results_flat)
n_fail <- length(results_flat) - length(took)
if (n_fail > 0L)
  message(sprintf("Note: %d/%d replicates failed the takeoff criterion (excluded).",
                  n_fail, length(results_flat)))

combo_key <- function(r) paste(sprintf("%.2f", r$r0), r$scenario, sep = "__")
by_combo  <- split(took, vapply(took, combo_key, character(1)))

cumulative_summary <- do.call(rbind, lapply(R0_GRID, function(r0) {
  key <- paste(sprintf("%.2f", r0), SCENARIO_TO_RUN, sep = "__")
  rs  <- by_combo[[key]]
  if (is.null(rs) || length(rs) == 0L) return(NULL)
  M <- do.call(rbind, lapply(rs, function(r) r$cum_at))
  data.frame(r0 = r0, scenario = SCENARIO_TO_RUN, timepoint_day = TIMEPOINTS,
             n_runs = nrow(M),
             mean   = apply(M, 2, mean),
             median = apply(M, 2, stats::median),
             q025   = apply(M, 2, safe_q, 0.025),
             q25    = apply(M, 2, safe_q, 0.25),
             q75    = apply(M, 2, safe_q, 0.75),
             q975   = apply(M, 2, safe_q, 0.975),
             stringsAsFactors = FALSE)
}))
cumulative_summary$date     <- day_to_date(cumulative_summary$timepoint_day)
cumulative_summary$r0_label <- sprintf("R0 = %.2f", cumulative_summary$r0)

# ============================================================================
# 10. Summary table (05 section 16) -- scenario kept as character
# ============================================================================
summary_rows <- do.call(rbind, lapply(R0_GRID, function(r0) {
  key <- paste(sprintf("%.2f", r0), SCENARIO_TO_RUN, sep = "__")
  rs  <- by_combo[[key]]
  if (is.null(rs) || length(rs) == 0L) return(NULL)
  do.call(rbind, lapply(seq_along(SUMMARY_DATES), function(di) {
    day <- as.integer(SUMMARY_DATES[di] - epi_start)
    idx <- match(day, TIMEPOINTS)
    if (is.na(idx)) idx <- which.min(abs(TIMEPOINTS - day))
    vals <- vapply(rs, function(r) r$cum_at[idx], numeric(1))
    data.frame(r0 = r0, scenario = SCENARIO_TO_RUN,
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

final_rows <- do.call(rbind, lapply(R0_GRID, function(r0) {
  key  <- paste(sprintf("%.2f", r0), SCENARIO_TO_RUN, sep = "__")
  rs   <- by_combo[[key]]
  if (is.null(rs) || length(rs) == 0L) return(NULL)
  vals    <- vapply(rs, function(r) as.numeric(r$n_cases), numeric(1))
  n_trunc <- sum(vals >= CHECK_FINAL_SIZE)
  data.frame(r0 = r0, scenario = SCENARIO_TO_RUN,
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

summary_tbl <- rbind(summary_rows, final_rows)
summary_tbl <- summary_tbl[order(summary_tbl$r0), ]
cat("\n--- Summary table (", SCENARIO_TO_RUN, ") ---\n", sep = "")
print(summary_tbl, row.names = FALSE)

# ============================================================================
# 11. Bundle and save (mirrors 05's 05_projections.rds, minus plots)
# ============================================================================
saveRDS(list(
  results_flat       = results_flat,
  cumulative_summary = cumulative_summary,
  summary_tbl        = summary_tbl,
  rt_profiles        = rt_profiles,
  config = list(
    scenario            = SCENARIO_TO_RUN,
    scen_idx_full       = scen_idx_full,
    r0_grid             = R0_GRID,
    scalar_overrides    = SCALAR_OVERRIDES,
    npi_specs           = NPI_SPECS,
    funeral_frac        = FUNERAL_FRAC,
    seeding_cases       = SEEDING_CASES,
    n_stoch             = N_STOCH,
    takeoff_n           = TAKEOFF_N,
    max_retries         = MAX_RETRIES,
    check_final_size    = CHECK_FINAL_SIZE,
    epidemic_start_date = epi_start,
    summary_dates       = SUMMARY_DATES,
    seed_base           = SEED_BASE_NEW,
    n_workers           = N_WORKERS,
    run_date            = Sys.Date()
  )
), file.path(out_dir, "05b_logistic_good_projections.rds"))

message(sprintf(
  "\n05b_projection_logistic_good.R complete.\n  %d R0 x 1 scenario x %d reps = %d runs.\n  Saved to: %s",
  n_r0, N_STOCH, total, out_dir))
