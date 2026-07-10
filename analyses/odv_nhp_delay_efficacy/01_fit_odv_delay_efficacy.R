# ============================================================
# ODV NHP delay-to-initiation efficacy analysis
# ============================================================
#
# Purpose
# -------
# This script reads hand-extracted NHP survival data from data-raw,
# converts the arm-level entries into animal-level survival data, estimates
# empirical hazard-scale efficacy at each treatment-initiation day, fits a
# simple delay-efficacy curve, and saves one reusable .rds object.
#
# Intended use
# ------------
# The fitted curve is intended as an empirical guide for sensitivity analyses
# of delayed antiviral post-exposure prophylaxis (PEP). It should not be
# interpreted as a precise estimate of human efficacy.
#
# Expected working directory
# --------------------------
# Run this script from the repository root, for example:
#   Rscript analyses/odv_nhp_delay_efficacy/01_fit_odv_delay_efficacy.R
#
# Output
# ------
# This script writes ONE processed output only:
#   data-processed/odv_nhp_delay/odv_ebov_rhesus_delay_efficacy_fit.rds
#
# The .rds contains raw data, processed survival data, empirical points,
# fitted curve, optimisation results, fit summary, and metadata.
#
# ============================================================


# ------------------------------------------------------------
# 0) Load packages
# ------------------------------------------------------------
# survival: Kaplan-Meier estimates and survival-data splitting.
# MASS: multivariate normal draws for approximate parameter uncertainty.
# Package startup messages are suppressed so automated runs remain readable.

suppressPackageStartupMessages({
  library(survival)
  library(MASS)
})


# ------------------------------------------------------------
# 1) Define paths and analysis settings
# ------------------------------------------------------------
# Paths are deliberately relative to the repository root so that the script
# can be run reproducibly on another machine after cloning the repo.

raw_path <- "data-raw/odv_nhp_delay/odv_ebov_rhesus_delay_survival.csv"
out_dir  <- "data-processed/odv_nhp_delay"
out_path <- file.path(out_dir, "odv_ebov_rhesus_delay_efficacy_fit.rds")

# Stop early with an informative message if the script is not being run from
# the repository root or if the raw-data file has not been copied into place.
if (!file.exists(raw_path)) {
  stop(
    "Could not find raw data file: ", raw_path, "\n",
    "Run this script from the repository root and check that data-raw is present."
  )
}

# Create the processed-data output directory if it does not already exist.
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Settings are stored together so they are saved into the output .rds and can
# be inspected later without reading the script.
settings <- list(
  t_censor_plot = 28,  # follow-up/censoring day read from the slide
  t_fit_end     = 15,  # final event time included in the likelihood fit
  dpc_zero      = 15,  # initiation day by which efficacy is constrained to zero
  k_fixed       = 1,   # fixed logistic steepness to avoid overfitting
  eps_hr        = 1e-6,# lower numerical bound for hazard ratios
  B_emp         = 500, # bootstrap replicates for empirical point uncertainty
  nsim_param    = 3000,# approximate parameter draws for fitted-curve uncertainty
  seed          = 1    # fixed seed for reproducibility
)

set.seed(settings$seed)


# ------------------------------------------------------------
# 2) Helper functions for raw-data parsing and simple transforms
# ------------------------------------------------------------
# The raw CSV stores death times as a semicolon-separated string because the
# original data are arm-level hand extractions from Kaplan-Meier curves.

parse_death_times <- function(x) {
  # Convert entries such as "7;8;8;9;11" into numeric death times.
  # Blank or missing strings are treated as arms with no observed deaths.
  x <- trimws(as.character(x))
  if (is.na(x) || x == "") return(numeric(0))
  as.numeric(strsplit(x, ";", fixed = TRUE)[[1]])
}

clamp01 <- function(x) {
  # Keep efficacy values within [0, 1] after numerical approximation.
  pmin(1, pmax(0, x))
}

km_surv_at <- function(dat, t_eval) {
  # Estimate Kaplan-Meier survival at a fixed evaluation time.
  # extend = TRUE carries the final survival estimate forward to t_eval.
  sf <- survfit(Surv(time, event) ~ 1, data = dat)
  as.numeric(summary(sf, times = t_eval, extend = TRUE)$surv[1])
}


# ------------------------------------------------------------
# 3) Read the raw arm-level data
# ------------------------------------------------------------
# Each row in the raw CSV represents one treatment arm: vehicle or ODV started
# at 1, 2, 3, or 4 days post-challenge (dpc).

raw_dat <- read.csv(raw_path, stringsAsFactors = FALSE)

# Check that the expected columns are present before doing any processing.
required_cols <- c(
  "arm", "dpc", "n_total", "n_survivors_reported",
  "death_times", "censor_day"
)
missing_cols <- setdiff(required_cols, names(raw_dat))
if (length(missing_cols) > 0) {
  stop("Raw data file is missing required columns: ", paste(missing_cols, collapse = ", "))
}

# Derive deaths and survivors from the death-time strings so the script checks
# that the hand extraction is internally consistent with the reported survival.
raw_dat$n_deaths_extracted <- vapply(
  raw_dat$death_times,
  function(x) length(parse_death_times(x)),
  integer(1)
)
raw_dat$n_survivors_extracted <- raw_dat$n_total - raw_dat$n_deaths_extracted

# If this fails, either the death times or the reported survivor counts in the
# raw CSV need to be corrected before fitting the curve.
if (any(raw_dat$n_survivors_extracted != raw_dat$n_survivors_reported)) {
  stop("Extracted death times do not match reported survivor counts in raw CSV.")
}


# ------------------------------------------------------------
# 4) Expand arm-level data to animal-level survival data
# ------------------------------------------------------------
# The survival likelihood is easiest to write using one row per animal.
# Deaths receive their extracted death time and event = 1.
# Survivors are censored at the common follow-up day and receive event = 0.

make_arm_df <- function(row) {
  # Expand one raw-data row into animal-level records for that arm.
  death_times <- parse_death_times(row$death_times)
  n_total <- as.integer(row$n_total)
  n_deaths <- length(death_times)
  n_cens <- n_total - n_deaths

  # Guard against impossible hand-entered data.
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

# Ensure numeric/integer types after expansion.
individual_dat$time  <- as.numeric(individual_dat$time)
individual_dat$event <- as.integer(individual_dat$event)

# Split the animal-level data into comparator and treated datasets for the
# empirical efficacy calculations below.
vehicle_dat <- individual_dat[individual_dat$arm == "vehicle", ]
odv_dat     <- individual_dat[individual_dat$arm != "vehicle", ]


# ------------------------------------------------------------
# 5) Estimate empirical hazard-scale efficacy at observed dpc values
# ------------------------------------------------------------
# Empirical efficacy is summarised as:
#   1 - H_treated / H_vehicle
# where H = -log(S_KM) at 28 days.
#
# This puts the empirical points on the same hazard-ratio scale used by the
# fitted model. It is descriptive and should not be over-interpreted given the
# small NHP sample size.

S_vehicle <- km_surv_at(vehicle_dat, settings$t_censor_plot)
H_vehicle <- -log(S_vehicle)

# Estimate one empirical point for each observed ODV initiation day.
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

# Use a simple non-parametric bootstrap to describe empirical-point uncertainty.
# This uncertainty is saved for plotting/inspection, but the model fit itself
# uses the full split survival likelihood below.
boot_eff <- matrix(NA_real_, nrow = settings$B_emp, ncol = length(observed_dpc))
colnames(boot_eff) <- paste0("dpc", observed_dpc)

for (b in seq_len(settings$B_emp)) {
  # Resample vehicle animals with replacement and recompute vehicle hazard.
  vehicle_b <- vehicle_dat[sample(seq_len(nrow(vehicle_dat)), replace = TRUE), ]
  H_vehicle_b <- -log(km_surv_at(vehicle_b, settings$t_censor_plot))

  # Skip bootstrap replicates with degenerate vehicle hazards.
  if (!is.finite(H_vehicle_b) || H_vehicle_b <= 0) next

  # Resample each treated arm independently and recompute hazard-scale efficacy.
  for (j in seq_along(observed_dpc)) {
    d <- observed_dpc[j]
    dat_d <- odv_dat[odv_dat$dpc == d, ]
    dat_b <- dat_d[sample(seq_len(nrow(dat_d)), replace = TRUE), ]
    H_d_b <- -log(km_surv_at(dat_b, settings$t_censor_plot))

    if (is.finite(H_d_b)) boot_eff[b, j] <- 1 - H_d_b / H_vehicle_b
  }
}

# Add percentile bootstrap intervals to the empirical-point table.
empirical_points$efficacy_lo <- NA_real_
empirical_points$efficacy_hi <- NA_real_

for (j in seq_along(observed_dpc)) {
  d <- observed_dpc[j]
  vals <- boot_eff[, j]
  vals <- vals[is.finite(vals)]

  # Require enough non-degenerate bootstrap replicates to report an interval.
  if (length(vals) >= 50) {
    empirical_points$efficacy_lo[empirical_points$dpc == d] <-
      quantile(vals, 0.025, na.rm = TRUE)
    empirical_points$efficacy_hi[empirical_points$dpc == d] <-
      quantile(vals, 0.975, na.rm = TRUE)
  }
}

# Clamp empirical values to the interpretable efficacy range.
empirical_points$efficacy_hazard_scale <- clamp01(empirical_points$efficacy_hazard_scale)
empirical_points$efficacy_lo <- clamp01(empirical_points$efficacy_lo)
empirical_points$efficacy_hi <- clamp01(empirical_points$efficacy_hi)


# ------------------------------------------------------------
# 6) Prepare split survival data for the likelihood fit
# ------------------------------------------------------------
# The likelihood uses piecewise-constant baseline hazards over intervals defined
# by observed death times and treatment initiation times. Treatment changes the
# post-initiation hazard in each ODV arm via a delay-specific hazard multiplier.

all_deaths <- sort(unique(individual_dat$time[individual_dat$event == 1L]))
all_deaths <- all_deaths[all_deaths <= settings$t_fit_end]

# Include death times, observed ODV initiation times, and the fitting horizon as
# interval breakpoints.
breaks_full <- sort(unique(c(0, all_deaths, observed_dpc, settings$t_fit_end)))
cuts_for_split <- breaks_full[-c(1, length(breaks_full))]

# Truncate follow-up at t_fit_end for the fitted likelihood.
fit_dat <- individual_dat
fit_dat$time_fit <- pmin(fit_dat$time, settings$t_fit_end)
fit_dat$event_fit <- as.integer(
  fit_dat$event == 1L & fit_dat$time <= settings$t_fit_end
)

# survSplit converts each animal into one row per risk interval.
split_dat <- survSplit(
  Surv(time_fit, event_fit) ~ .,
  data = fit_dat,
  cut = cuts_for_split,
  episode = "interval"
)

# Compute interval start, stop, and exposure time for each split row.
split_dat$tstart   <- as.numeric(split_dat$tstart)
split_dat$tstop    <- as.numeric(split_dat$time_fit)
split_dat$exposure <- split_dat$tstop - split_dat$tstart

# Remove any zero-length rows introduced by the splitting procedure.
split_dat <- split_dat[split_dat$exposure > 0, ]

# Store integer interval indices for fast likelihood calculations.
split_dat$interval <- factor(split_dat$interval)
split_dat$interval_idx <- as.integer(split_dat$interval)
K <- length(levels(split_dat$interval))

# Treatment only affects intervals after the relevant ODV initiation time.
# Vehicle rows have dpc = 0 and therefore remain untreated throughout.
split_dat$post_treatment <- as.integer(
  split_dat$dpc > 0 & split_dat$tstart >= split_dat$dpc
)


# ------------------------------------------------------------
# 7) Define the delay-efficacy curve
# ------------------------------------------------------------
# The fitted curve is a scaled logistic decline:
#   efficacy(d) = E0 * {g(d) - g(d_zero)} / {g(0) - g(d_zero)}
#   g(d) = 1 / [1 + exp(k * (d - d50))]
#
# E0 is the fitted efficacy at dpc = 0 on the hazard scale.
# d50 controls where the decline occurs.
# k is fixed to avoid overfitting the small NHP dataset.
# d_zero fixes efficacy to zero by a late initiation day.

logistic_g <- function(d, k, d50) {
  # Basic logistic component used to define the decline in efficacy with delay.
  1 / (1 + exp(k * (d - d50)))
}

efficacy_curve <- function(d, E0, d50, k = settings$k_fixed,
                           d_zero = settings$dpc_zero) {
  # Evaluate the scaled delay-efficacy curve at initiation day d.
  g0 <- logistic_g(0, k, d50)
  gz <- logistic_g(d_zero, k, d50)
  gd <- logistic_g(d, k, d50)

  # Scaling ensures efficacy is E0 at d = 0 and 0 at d = d_zero.
  eff <- E0 * (gd - gz) / (g0 - gz)

  # Replace numerical failures with NA and force zero at/after d_zero.
  eff[!is.finite(eff)] <- NA_real_
  eff[d >= d_zero] <- 0
  eff
}


# ------------------------------------------------------------
# 8) Define the profiled negative log-likelihood
# ------------------------------------------------------------
# For a given delay-efficacy curve, the baseline interval hazards are profiled
# out analytically. This keeps optimisation to two shape parameters: E0 and d50.

row_hazard_multiplier <- function(par, dat) {
  # Convert delay-specific efficacies into hazard multipliers for split rows.
  E0 <- par[1]
  d50 <- par[2]

  # Evaluate efficacy only at observed ODV initiation days.
  eff_d <- efficacy_curve(observed_dpc, E0, d50)
  names(eff_d) <- as.character(observed_dpc)

  # Hazard ratio = 1 - efficacy on the hazard scale.
  hr_d <- 1 - eff_d

  # Reject invalid or near-zero hazard ratios.
  if (any(!is.finite(hr_d)) || any(hr_d <= settings$eps_hr)) return(NULL)

  # Vehicle and pre-treatment rows have multiplier 1.
  mult <- rep(1, nrow(dat))

  # Post-treatment rows get the multiplier corresponding to their dpc arm.
  idx <- which(dat$post_treatment == 1L)
  if (length(idx) > 0) {
    mult[idx] <- hr_d[as.character(dat$dpc[idx])]
  }

  mult
}

profile_nll <- function(par, dat, K) {
  # Negative log-likelihood for the piecewise exponential survival model.
  E0 <- par[1]
  d50 <- par[2]

  # Constrain fitted parameters to interpretable ranges.
  if (!is.finite(E0) || E0 < 0 || E0 > 1 - settings$eps_hr) return(1e30)
  if (!is.finite(d50) || d50 < 0 || d50 > 30) return(1e30)

  # Reject curves with invalid observed-day efficacies.
  eff_d <- efficacy_curve(observed_dpc, E0, d50)
  if (any(!is.finite(eff_d)) || any(eff_d < 0) || any(eff_d > 1 - settings$eps_hr)) {
    return(1e30)
  }

  # Convert the curve into row-level hazard multipliers.
  mult <- row_hazard_multiplier(par, dat)
  if (is.null(mult)) return(1e30)

  # Profile the baseline hazard separately in each interval.
  idx_fac <- factor(dat$interval_idx, levels = seq_len(K))
  exposure_k <- tapply(dat$exposure * mult, idx_fac, sum)
  deaths_k   <- tapply(dat$event_fit, idx_fac, sum)

  # Replace empty intervals with zeros for safe indexing.
  exposure_k[is.na(exposure_k)] <- 0
  deaths_k[is.na(deaths_k)] <- 0

  # Any interval with deaths but no exposure is impossible.
  if (any(exposure_k == 0 & deaths_k > 0)) return(1e30)

  # Maximum-likelihood baseline hazard for each interval.
  lambda_k <- rep(0, K)
  ok <- which(exposure_k > 0)
  lambda_k[ok] <- deaths_k[ok] / exposure_k[ok]

  # Expected events for each split row.
  mu <- dat$exposure * mult * lambda_k[dat$interval_idx]
  y <- dat$event_fit

  # Guard against invalid Poisson means.
  if (any(y == 1L & mu <= 0) || any(mu < 0) || any(!is.finite(mu))) return(1e30)

  # Piecewise-exponential log-likelihood is equivalent to a Poisson likelihood
  # for deaths with exposure offsets; constants are omitted.
  -sum(ifelse(y == 1L, log(mu), 0) - mu)
}


# ------------------------------------------------------------
# 9) Fit the curve
# ------------------------------------------------------------
# Optimise E0 and d50 by minimising the profiled negative log-likelihood.
# The starting values are deliberately simple because there are few data points.

fit <- optim(
  par = c(E0 = 0.95, d50 = 6),
  fn = profile_nll,
  method = "L-BFGS-B",
  lower = c(0, 0),
  upper = c(1 - settings$eps_hr, 30),
  dat = split_dat,
  K = K,
  control = list(maxit = 5000)
)

# Extract fitted parameters.
E0_hat <- fit$par[["E0"]]
d50_hat <- fit$par[["d50"]]

# Compare fitted values with the empirical points at the observed initiation days.
pred_at_observed <- clamp01(efficacy_curve(observed_dpc, E0_hat, d50_hat))

# Summarise goodness-of-fit using quantities that are easy to inspect later.
n_shape_par <- 2
n_par_total <- K + n_shape_par
nll <- fit$value
aic <- 2 * n_par_total + 2 * nll
rmse_emp <- sqrt(mean((pred_at_observed - empirical_points$efficacy_hazard_scale)^2))


# ------------------------------------------------------------
# 10) Approximate uncertainty in the fitted curve
# ------------------------------------------------------------
# Approximate uncertainty is obtained by numerically estimating the Hessian of
# the profiled negative log-likelihood and drawing parameters from the local
# asymptotic normal approximation. This is descriptive, not definitive.

hessian2 <- function(f, par, step = c(1e-3, 1e-2), ...) {
  # Numerical Hessian for the two fitted parameters, E0 and d50.
  p <- as.numeric(par)
  h <- as.numeric(step)
  H <- matrix(0, 2, 2)
  f0 <- f(p, ...)

  # Diagonal second derivatives.
  for (i in 1:2) {
    pp <- p
    pm <- p
    pp[i] <- p[i] + h[i]
    pm[i] <- p[i] - h[i]
    H[i, i] <- (f(pp, ...) - 2 * f0 + f(pm, ...)) / h[i]^2
  }

  # Off-diagonal mixed derivative.
  ppp <- p
  ppm <- p
  pmp <- p
  pmm <- p
  ppp[1] <- p[1] + h[1]; ppp[2] <- p[2] + h[2]
  ppm[1] <- p[1] + h[1]; ppm[2] <- p[2] - h[2]
  pmp[1] <- p[1] - h[1]; pmp[2] <- p[2] + h[2]
  pmm[1] <- p[1] - h[1]; pmm[2] <- p[2] - h[2]

  H[1, 2] <- (f(ppp, ...) - f(ppm, ...) - f(pmp, ...) + f(pmm, ...)) / (4 * h[1] * h[2])
  H[2, 1] <- H[1, 2]

  # Symmetrise to reduce numerical noise.
  (H + t(H)) / 2
}

# Curve grid to save for downstream use.
curve_grid <- seq(0, settings$dpc_zero, by = 0.02)
curve_dat <- data.frame(
  dpc = curve_grid,
  efficacy = clamp01(efficacy_curve(curve_grid, E0_hat, d50_hat)),
  efficacy_lo = NA_real_,
  efficacy_hi = NA_real_
)

# Estimate the covariance matrix for E0 and d50 if the Hessian is invertible.
H <- hessian2(profile_nll, c(E0_hat, d50_hat), dat = split_dat, K = K)
V <- tryCatch(solve(H), error = function(e) NULL)

if (!is.null(V) && all(is.finite(V))) {
  # Force the covariance matrix to be positive semi-definite if numerical noise
  # creates tiny negative eigenvalues.
  eig <- eigen(V, symmetric = TRUE)
  eig$values[eig$values < 1e-8] <- 1e-8
  V <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)

  # Draw plausible parameter values from the local normal approximation.
  draws <- MASS::mvrnorm(settings$nsim_param, mu = c(E0_hat, d50_hat), Sigma = V)

  # Keep only draws within the same parameter bounds used in optimisation.
  keep <- draws[, 1] >= 0 & draws[, 1] <= 1 - settings$eps_hr &
    draws[, 2] >= 0 & draws[, 2] <= 30
  draws <- draws[keep, , drop = FALSE]

  # Summarise pointwise curve uncertainty if enough valid draws remain.
  if (nrow(draws) >= 200) {
    eff_mat <- vapply(
      curve_grid,
      function(d) efficacy_curve(d, draws[, 1], draws[, 2]),
      FUN.VALUE = numeric(nrow(draws))
    )

    curve_dat$efficacy_lo <- clamp01(apply(eff_mat, 2, quantile, 0.025, na.rm = TRUE))
    curve_dat$efficacy_hi <- clamp01(apply(eff_mat, 2, quantile, 0.975, na.rm = TRUE))
  }
}


# ------------------------------------------------------------
# 11) Build a compact fit summary
# ------------------------------------------------------------
# This table is saved inside the .rds so that downstream scripts can inspect
# the fitted parameters without digging through the optimisation object.

fit_summary <- data.frame(
  E0_hat = E0_hat,
  d50_hat = d50_hat,
  k_fixed = settings$k_fixed,
  dpc_zero = settings$dpc_zero,
  nll = nll,
  aic = aic,
  rmse_empirical_points = rmse_emp,
  n_animals = nrow(individual_dat),
  n_events = sum(individual_dat$event),
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 12) Save one reusable processed output
# ------------------------------------------------------------
# Charlie requested an .rds object for downstream analyses. To keep the analysis
# footprint minimal, this script deliberately does not write separate CSV files.
# Everything needed later is stored in the single output list below.

output <- list(
  metadata = list(
    analysis = "ODV NHP delay-to-initiation efficacy curve",
    source = unique(raw_dat$source_note),
    interpretation = paste(
      "Parsimonious empirical guide for sensitivity analyses of delayed PEP;",
      "not a precise estimate of human efficacy."
    ),
    raw_path = raw_path,
    output_path = out_path,
    created_by = "analyses/odv_nhp_delay_efficacy/01_fit_odv_delay_efficacy.R",
    settings = settings
  ),
  raw_data = raw_dat,
  individual_survival_data = individual_dat,
  split_survival_data = split_dat,
  empirical_points = empirical_points,
  fitted_curve = curve_dat,
  fit_summary = fit_summary,
  optim_fit = fit
)

saveRDS(output, out_path)

message("Saved ODV NHP delay efficacy fit to: ", out_path)
print(fit_summary)
