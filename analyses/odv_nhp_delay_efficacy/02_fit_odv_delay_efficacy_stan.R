# ============================================================================
# ODV NHP delay-to-initiation efficacy analysis -- FULL Bayesian (Stan) version
# ============================================================================
#
# Purpose
# -------
# This is the Stan counterpart to 01_fit_odv_delay_efficacy.R. It performs the
# SAME data processing (parse the hand-extracted NHP survival data, expand to
# animal-level records, compute empirical hazard-scale efficacy points, and
# build the split survival data), but then fits the COMPLETE survival model in
# Stan instead of profiling out the baseline hazards.
#
# What "complete" means here
# ---------------------------
# Script 01 profiles the K piecewise-constant baseline interval hazards out of
# the likelihood analytically and approximates uncertainty in the two curve
# parameters (E0, d50) with a 2-parameter Laplace / multivariate-normal step.
# This script instead estimates ALL parameters jointly -- the K baseline hazards
# AND the curve parameters -- by HMC/NUTS, so baseline-hazard uncertainty,
# parameter correlation and posterior skew are propagated into the fitted curve.
#
# Two user choices (set in the block below)
# -----------------------------------------
#   fit_k            : FALSE -> fix the logistic steepness k at k_fixed (mirrors
#                      script 01, k = 1).  TRUE -> estimate k with a weakly-
#                      informative lognormal prior centred on k_fixed.
#   use_hazard_prior : 0 -> flat (improper) prior on each baseline hazard
#                      lambda_k >= 0 (the most uninformative option; posterior is
#                      still proper here). 1 -> a very diffuse but proper
#                      exponential prior, as a fallback.
#
# Intended use
# ------------
# As in script 01, the fitted curve is an empirical guide for sensitivity
# analyses of delayed antiviral PEP. It is NOT a precise estimate of human
# efficacy.
#
# Expected working directory
# ---------------------------
# Paths are resolved with here::here(), so the script can be run from anywhere
# inside the repository, for example:
#   Rscript analyses/odv_nhp_delay_efficacy/02_fit_odv_delay_efficacy_stan.R
#
# Output
# ------
# Writes ONE processed output:
#   data-processed/odv_nhp_delay/odv_ebov_rhesus_delay_efficacy_fit_stan.rds
# (a sibling of script 01's output; this script does not overwrite it). The .rds
# contains raw data, processed survival data, empirical points, the posterior
# fitted curve, parameter summaries, sampler diagnostics, a compact draws table,
# the Stan data list and metadata.
#
# ============================================================================


# ------------------------------------------------------------
# 0) Load packages
# ------------------------------------------------------------
# survival : Kaplan-Meier estimates and survival-data splitting (as in 01).
# cmdstanr : the repo's Stan interface (see analyses/01_latent_response_*).
# Package startup messages are suppressed so automated runs remain readable.

suppressPackageStartupMessages({
  library(survival)
  library(cmdstanr)
})


# ------------------------------------------------------------
# 1) USER CHOICES
# ------------------------------------------------------------
# ---- k (logistic steepness): fix it, or fit it with a prior ----------------
# fit_k = FALSE -> fix k at k_fixed (identical treatment to script 01).
# fit_k = TRUE  -> estimate k with a lognormal prior whose median is k_fixed.
fit_k           <- FALSE
k_fixed         <- 1
k_prior_logmean <- log(k_fixed)  # lognormal median = k_fixed
k_prior_logsd   <- 0.5           # ~95% of prior mass within roughly [0.37, 2.7] * k_fixed

# ---- Baseline interval hazards: prior --------------------------------------
# Per request, the default is the most uninformative option: a flat (improper)
# prior on each lambda_k >= 0. Stan adds the lower-bound Jacobian automatically,
# so omitting a sampling statement yields a uniform prior on the constrained
# (>= 0) scale. The posterior remains proper because every retained interval has
# positive total exposure (and most have >= 1 death).
#
# Set use_hazard_prior = 1L only if you want a guaranteed-proper prior or see
# sampling pathologies far out in the tail; the exponential mean is then huge
# (1 / hazard_prior_rate), so it stays extremely diffuse.
use_hazard_prior  <- 0L
hazard_prior_rate <- 1e-3

# ---- d50 (decline location): weakly-informative prior ----------------------
# The likelihood only sees efficacy at the observed initiation days (1-4 dpc),
# and for d50 beyond ~5-6 the curve is flat across them -- so the decline
# LOCATION is not identified by these data (a flat prior lets the posterior
# wander out toward d_zero, giving a flat-then-late-drop median curve unlike the
# profiled fit). We therefore centre d50 in the observed window and cap it at
# d_zero (in the Stan model). Set d50_prior_sd very large to recover a flat prior.
d50_prior_mean <- 4
d50_prior_sd   <- 2


# ------------------------------------------------------------
# 2) Define paths and analysis settings
# ------------------------------------------------------------
# Paths are resolved with here::here() so the script runs from anywhere in repo.

raw_path  <- here::here("data-raw", "odv_nhp_delay", "odv_ebov_rhesus_delay_survival.csv")
out_dir   <- here::here("data-processed", "odv_nhp_delay")
out_path  <- file.path(out_dir, "odv_ebov_rhesus_delay_efficacy_fit_stan.rds")
stan_file <- here::here("analyses", "odv_nhp_delay_efficacy",
                        "stan-models", "odv_delay_efficacy.stan")

# Stop early with an informative message if inputs are missing.
if (!file.exists(raw_path)) {
  stop(
    "Could not find raw data file: ", raw_path, "\n",
    "Check that data-raw is present in the cloned repository."
  )
}
if (!file.exists(stan_file)) {
  stop("Could not find Stan model file: ", stan_file)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Settings are stored together so they are saved into the output .rds and can be
# inspected later without reading the script. (Shared values mirror script 01.)
settings <- list(
  # data / curve settings (shared with script 01)
  t_censor_plot = 28,   # follow-up/censoring day read from the slide
  t_fit_end     = 15,   # final event time included in the likelihood fit
  dpc_zero      = 15,   # initiation day by which efficacy is constrained to zero
  eps_hr        = 1e-6, # lower numerical bound for hazard ratios
  B_emp         = 500,  # bootstrap replicates for empirical point uncertainty
  seed          = 1,    # fixed seed for reproducibility

  # k handling
  fit_k           = fit_k,
  k_fixed         = k_fixed,
  k_prior_logmean = k_prior_logmean,
  k_prior_logsd   = k_prior_logsd,

  # baseline-hazard prior
  use_hazard_prior  = use_hazard_prior,
  hazard_prior_rate = hazard_prior_rate,

  # d50 prior
  d50_prior_mean = d50_prior_mean,
  d50_prior_sd   = d50_prior_sd,

  # Stan sampler settings
  chains        = 4,
  iter_warmup   = 1000,
  iter_sampling = 1000,
  adapt_delta   = 0.95, # a little above default given bounded params / small data
  max_treedepth = 12
)

set.seed(settings$seed)


# ------------------------------------------------------------
# 3) Helper functions for raw-data parsing and simple transforms
# ------------------------------------------------------------
# Identical to script 01: the raw CSV stores death times as a semicolon-
# separated string because the original data are arm-level hand extractions.

parse_death_times <- function(x) {
  # Convert entries such as "7;8;8;9;11" into numeric death times.
  x <- trimws(as.character(x))
  if (is.na(x) || x == "") return(numeric(0))
  as.numeric(strsplit(x, ";", fixed = TRUE)[[1]])
}

clamp01 <- function(x) {
  # Keep efficacy values within [0, 1] after numerical approximation.
  pmin(1, pmax(0, x))
}

km_surv_at <- function(dat, t_eval) {
  # Kaplan-Meier survival at a fixed evaluation time; extend = TRUE carries the
  # final estimate forward to t_eval.
  sf <- survfit(Surv(time, event) ~ 1, data = dat)
  as.numeric(summary(sf, times = t_eval, extend = TRUE)$surv[1])
}


# ------------------------------------------------------------
# 4) Read raw arm-level data and expand to animal-level data
# ------------------------------------------------------------
# This block reproduces steps 3-4 of script 01 so the two fits consume exactly
# the same processed data.

raw_dat <- read.csv(raw_path, stringsAsFactors = FALSE)

required_cols <- c(
  "arm", "dpc", "n_total", "n_survivors_reported",
  "death_times", "censor_day"
)
missing_cols <- setdiff(required_cols, names(raw_dat))
if (length(missing_cols) > 0) {
  stop("Raw data file is missing required columns: ", paste(missing_cols, collapse = ", "))
}

# Derive deaths/survivors from the death-time strings and check internal
# consistency with the reported survival counts.
raw_dat$n_deaths_extracted <- vapply(
  raw_dat$death_times,
  function(x) length(parse_death_times(x)),
  integer(1)
)
raw_dat$n_survivors_extracted <- raw_dat$n_total - raw_dat$n_deaths_extracted

if (any(raw_dat$n_survivors_extracted != raw_dat$n_survivors_reported)) {
  stop("Extracted death times do not match reported survivor counts in raw CSV.")
}

# Expand one raw-data row into animal-level records for that arm.
make_arm_df <- function(row) {
  death_times <- parse_death_times(row$death_times)
  n_total <- as.integer(row$n_total)
  n_deaths <- length(death_times)
  n_cens <- n_total - n_deaths
  if (n_cens < 0) stop("More death times than animals in arm: ", row$arm)

  data.frame(
    subject_id = sprintf("%s_%02d", row$arm, seq_len(n_total)),
    arm        = row$arm,
    dpc        = as.integer(row$dpc),
    time       = c(death_times, rep(as.numeric(row$censor_day), n_cens)),
    event      = c(rep(1L, n_deaths), rep(0L, n_cens)),
    stringsAsFactors = FALSE
  )
}

individual_dat <- do.call(
  rbind,
  lapply(seq_len(nrow(raw_dat)), function(i) make_arm_df(raw_dat[i, ]))
)
individual_dat$time  <- as.numeric(individual_dat$time)
individual_dat$event <- as.integer(individual_dat$event)

vehicle_dat <- individual_dat[individual_dat$arm == "vehicle", ]
odv_dat     <- individual_dat[individual_dat$arm != "vehicle", ]


# ------------------------------------------------------------
# 5) Empirical hazard-scale efficacy at observed dpc values
# ------------------------------------------------------------
# Reproduces step 5 of script 01 (descriptive only) so the saved .rds carries
# the same empirical points to plot against the Bayesian posterior curve.
#   efficacy = 1 - H_treated / H_vehicle,  H = -log(S_KM) at t_censor_plot.

S_vehicle <- km_surv_at(vehicle_dat, settings$t_censor_plot)
H_vehicle <- -log(S_vehicle)

observed_dpc <- sort(unique(odv_dat$dpc))

empirical_points <- do.call(rbind, lapply(observed_dpc, function(d) {
  dat_d <- odv_dat[odv_dat$dpc == d, ]
  S_d <- km_surv_at(dat_d, settings$t_censor_plot)
  H_d <- -log(S_d)
  data.frame(
    dpc = d,
    n_total = nrow(dat_d),
    n_survivors = sum(dat_d$event == 0L),
    survival = S_d,
    efficacy_hazard_scale = 1 - H_d / H_vehicle,
    stringsAsFactors = FALSE
  )
}))

# Non-parametric bootstrap for empirical-point uncertainty (plotting only).
boot_eff <- matrix(NA_real_, nrow = settings$B_emp, ncol = length(observed_dpc))
colnames(boot_eff) <- paste0("dpc", observed_dpc)

for (b in seq_len(settings$B_emp)) {
  vehicle_b <- vehicle_dat[sample(seq_len(nrow(vehicle_dat)), replace = TRUE), ]
  H_vehicle_b <- -log(km_surv_at(vehicle_b, settings$t_censor_plot))
  if (!is.finite(H_vehicle_b) || H_vehicle_b <= 0) next

  for (j in seq_along(observed_dpc)) {
    d <- observed_dpc[j]
    dat_d <- odv_dat[odv_dat$dpc == d, ]
    dat_b <- dat_d[sample(seq_len(nrow(dat_d)), replace = TRUE), ]
    H_d_b <- -log(km_surv_at(dat_b, settings$t_censor_plot))
    if (is.finite(H_d_b)) boot_eff[b, j] <- 1 - H_d_b / H_vehicle_b
  }
}

empirical_points$efficacy_lo <- NA_real_
empirical_points$efficacy_hi <- NA_real_
for (j in seq_along(observed_dpc)) {
  d <- observed_dpc[j]
  vals <- boot_eff[, j]
  vals <- vals[is.finite(vals)]
  if (length(vals) >= 50) {
    empirical_points$efficacy_lo[empirical_points$dpc == d] <- quantile(vals, 0.025, na.rm = TRUE)
    empirical_points$efficacy_hi[empirical_points$dpc == d] <- quantile(vals, 0.975, na.rm = TRUE)
  }
}
empirical_points$efficacy_hazard_scale <- clamp01(empirical_points$efficacy_hazard_scale)
empirical_points$efficacy_lo <- clamp01(empirical_points$efficacy_lo)
empirical_points$efficacy_hi <- clamp01(empirical_points$efficacy_hi)


# ------------------------------------------------------------
# 6) Prepare split survival data for the likelihood
# ------------------------------------------------------------
# Reproduces step 6 of script 01: piecewise-constant baseline hazards over
# intervals defined by death times and treatment-initiation times, truncated at
# t_fit_end. The resulting split_dat is exactly the design fed to Stan.

all_deaths <- sort(unique(individual_dat$time[individual_dat$event == 1L]))
all_deaths <- all_deaths[all_deaths <= settings$t_fit_end]

breaks_full <- sort(unique(c(0, all_deaths, observed_dpc, settings$t_fit_end)))
cuts_for_split <- breaks_full[-c(1, length(breaks_full))]

fit_dat <- individual_dat
fit_dat$time_fit <- pmin(fit_dat$time, settings$t_fit_end)
fit_dat$event_fit <- as.integer(
  fit_dat$event == 1L & fit_dat$time <= settings$t_fit_end
)

split_dat <- survSplit(
  Surv(time_fit, event_fit) ~ .,
  data = fit_dat,
  cut = cuts_for_split,
  episode = "interval"
)

split_dat$tstart   <- as.numeric(split_dat$tstart)
split_dat$tstop    <- as.numeric(split_dat$time_fit)
split_dat$exposure <- split_dat$tstop - split_dat$tstart
split_dat <- split_dat[split_dat$exposure > 0, ]

split_dat$interval <- factor(split_dat$interval)
split_dat$interval_idx <- as.integer(split_dat$interval)
K <- length(levels(split_dat$interval))

# Treatment only affects intervals after the relevant ODV initiation time.
split_dat$post_treatment <- as.integer(
  split_dat$dpc > 0 & split_dat$tstart >= split_dat$dpc
)


# ------------------------------------------------------------
# 7) Assemble the Stan data list
# ------------------------------------------------------------
# The curve grid mirrors the one saved by script 01 (0 .. dpc_zero by 0.02).

D <- length(observed_dpc)

# Map each split row to its arm's efficacy entry. Vehicle rows have dpc = 0 and
# no match; set those to 1 (a valid, but unused, index -- guarded in the model
# by post_treatment == 0).
dpc_idx <- match(split_dat$dpc, observed_dpc)
dpc_idx[is.na(dpc_idx)] <- 1L

curve_grid <- seq(0, settings$dpc_zero, by = 0.02)

stan_data <- list(
  N = nrow(split_dat),
  K = K,
  D = D,
  interval_idx   = as.integer(split_dat$interval_idx),
  event          = as.integer(split_dat$event_fit),
  exposure       = as.numeric(split_dat$exposure),
  post_treatment = as.integer(split_dat$post_treatment),
  dpc_idx        = as.integer(dpc_idx),
  dpc_obs        = as.numeric(observed_dpc),
  d_zero         = settings$dpc_zero,
  eps_hr         = settings$eps_hr,
  fit_k             = as.integer(settings$fit_k),
  k_fixed           = settings$k_fixed,
  k_prior_logmean   = settings$k_prior_logmean,
  k_prior_logsd     = settings$k_prior_logsd,
  use_hazard_prior  = as.integer(settings$use_hazard_prior),
  hazard_prior_rate = settings$hazard_prior_rate,
  d50_prior_mean    = settings$d50_prior_mean,
  d50_prior_sd      = settings$d50_prior_sd,
  n_grid   = length(curve_grid),
  grid_dpc = as.numeric(curve_grid)
)


# ------------------------------------------------------------
# 8) Compile and fit the Stan model
# ------------------------------------------------------------

mod <- cmdstan_model(stan_file)

fit <- mod$sample(
  data            = stan_data,
  seed            = settings$seed,
  chains          = settings$chains,
  parallel_chains = settings$chains,
  iter_warmup     = settings$iter_warmup,
  iter_sampling   = settings$iter_sampling,
  adapt_delta     = settings$adapt_delta,
  max_treedepth   = settings$max_treedepth,
  refresh         = 200
)

cat("\nSampler diagnostics:\n")
diag_summ <- fit$diagnostic_summary()
print(diag_summ)


# ------------------------------------------------------------
# 8b) MAP check: does the mode reproduce the profiled fit (script 01)?
# ------------------------------------------------------------
# Penalised MLE (jacobian = FALSE => mode on the constrained scale, matching a
# classical optimum rather than the unconstrained-scale mode). With flat priors
# this should land on script 01's optim result; with the weakly-informative d50
# prior it is the MAP. A quick confirmation that the likelihood/model match and a
# central estimate that is directly comparable to script 01.
map_estimate <- tryCatch({
  opt <- mod$optimize(data = stan_data, jacobian = FALSE, seed = settings$seed)
  os <- opt$summary(c("E0", "d50", "k"))
  cat("\nMAP (penalised MLE) estimate:\n"); print(os)
  as.data.frame(os)
}, error = function(e) {
  message("optimize() failed (", conditionMessage(e), "); skipping MAP check.")
  NULL
})


# ------------------------------------------------------------
# 9) Posterior summaries and fitted curve
# ------------------------------------------------------------
# The fitted curve uses pointwise posterior quantiles of efficacy_grid, replacing
# the Laplace / MVN-draw uncertainty of script 01.

# Curve parameters (k is reported whether fixed or fitted -- it is a transformed
# parameter in the model, so it is always available here).
param_summ <- fit$summary(c("E0", "d50", "k"))

# Pointwise efficacy curve with a 95% credible ribbon. Columns of the draws
# matrix are efficacy_grid[1..n_grid], in the same order as curve_grid.
grid_draws <- fit$draws(variables = "efficacy_grid", format = "draws_matrix")
grid_q <- apply(grid_draws, 2, quantile, probs = c(0.025, 0.5, 0.975), names = FALSE)
curve_dat <- data.frame(
  dpc         = curve_grid,
  efficacy    = grid_q[2, ],  # posterior median
  efficacy_lo = grid_q[1, ],  # 2.5%
  efficacy_hi = grid_q[3, ]   # 97.5%
)

# Posterior efficacy at the observed initiation days, to compare with empirical.
eff_summ <- fit$summary("eff")
pred_at_observed <- eff_summ$median
rmse_emp <- sqrt(mean((pred_at_observed - empirical_points$efficacy_hazard_scale)^2))

# Helper to pull a single parameter's posterior summary row.
get_par <- function(tab, nm, col) tab[[col]][tab$variable == nm]

fit_summary <- data.frame(
  E0_median  = get_par(param_summ, "E0", "median"),
  E0_lo      = get_par(param_summ, "E0", "q5"),
  E0_hi      = get_par(param_summ, "E0", "q95"),
  d50_median = get_par(param_summ, "d50", "median"),
  d50_lo     = get_par(param_summ, "d50", "q5"),
  d50_hi     = get_par(param_summ, "d50", "q95"),
  k_median   = get_par(param_summ, "k", "median"),
  k_estimated = as.logical(settings$fit_k),
  dpc_zero   = settings$dpc_zero,
  rmse_empirical_points = rmse_emp,
  n_animals  = nrow(individual_dat),
  n_events   = sum(individual_dat$event),
  n_divergent = sum(diag_summ$num_divergent),
  max_rhat    = max(param_summ$rhat, na.rm = TRUE),
  min_ess_bulk = min(param_summ$ess_bulk, na.rm = TRUE),
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 10) Save one reusable processed output
# ------------------------------------------------------------
# Structure mirrors script 01's output, with the Stan posterior replacing the
# optim object. A compact draws table for the curve parameters is stored (the
# full cmdstanr object is not, as it references temporary CSV files).

draws_df <- as.data.frame(fit$draws(c("E0", "d50", "k"), format = "df"))

output <- list(
  metadata = list(
    analysis = "ODV NHP delay-to-initiation efficacy curve (full Bayesian, Stan)",
    source = if ("source_note" %in% names(raw_dat)) unique(raw_dat$source_note) else NA,
    method = paste(
      "Full Bayesian piecewise-exponential survival model fit in Stan (cmdstanr);",
      "baseline interval hazards estimated jointly with the delay-efficacy curve",
      "parameters rather than profiled out."
    ),
    interpretation = paste(
      "Parsimonious empirical guide for sensitivity analyses of delayed PEP;",
      "not a precise estimate of human efficacy."
    ),
    raw_path = raw_path,
    output_path = out_path,
    stan_model = stan_file,
    created_by = "analyses/odv_nhp_delay_efficacy/02_fit_odv_delay_efficacy_stan.R",
    settings = settings,
    k_treatment = if (settings$fit_k) "estimated" else "fixed"
  ),
  raw_data = raw_dat,
  individual_survival_data = individual_dat,
  split_survival_data = split_dat,
  empirical_points = empirical_points,
  stan_data = stan_data,
  fitted_curve = curve_dat,
  param_summary = as.data.frame(param_summ),
  fit_summary = fit_summary,
  diagnostics = diag_summ,
  map_estimate = map_estimate,
  draws = draws_df
)

saveRDS(output, out_path)

message("Saved ODV NHP delay efficacy (Stan) fit to: ", out_path)
print(fit_summary)
