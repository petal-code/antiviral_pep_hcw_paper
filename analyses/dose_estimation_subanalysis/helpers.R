# ============================================================================
# helpers.R  (dose_estimation_subanalysis)
# ----------------------------------------------------------------------------
# Small shared utilities and path constants for the two scripts in this
# subanalysis:
#   01_fit_dose_q_curve.R              -- fits + extrapolates the dose Q curve
#   02_npi_inputs_and_fiber_runs.R     -- turns the Q curve into time-varying
#                                         NPI inputs and runs fiber across R0
#
# Paths are resolved from the repository root with here::here() so the scripts
# run the same regardless of the working directory they are launched from.
# ============================================================================

# ---- Project paths ---------------------------------------------------------
ANALYSIS_DIR <- here::here("analyses", "dose_estimation_subanalysis")
DIR_STAN     <- file.path(ANALYSIS_DIR, "stan-models")
DIR_OUT      <- file.path(ANALYSIS_DIR, "outputs")

# Ensure the outputs folder exists (it is created on first run; the generated
# .rds / .csv artefacts are regenerable and are not tracked in git).
dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)

# ---- The raw dose-coverage observations ------------------------------------
# Five (calendar date, percentage) pairs. Percentages are converted to
# proportions ( /100 ) downstream. The first point is 0% on 18 May 2026.
DOSE_OBS <- data.frame(
  date       = as.Date(c("2026-05-18", "2026-05-24", "2026-05-31",
                         "2026-06-07", "2026-06-14")),
  percentage = c(0.0, 19.3, 30.2, 64.4, 63.1),
  stringsAsFactors = FALSE
)

# ---- Tiny numeric helper ----------------------------------------------------
clip01 <- function(x) pmin(1, pmax(0, x))

# ---- Build the relative-day observation table (with an optional start anchor) -
# Converts the calendar-date observations to a RELATIVE-DAY axis measured from
# `start_date` (day 0), and turns percentages into proportions.
#
# Start anchor: if `start_date` is earlier than the first observation, a SINGLE
# zero-valued point is added at the start date (relative day 0) to anchor the
# curve at the 0% floor from that day. (Earlier days are NOT each padded with a
# zero.) With start_date == the first observation date (the default) nothing is
# added and the relative days are simply {0, 6, 13, 20, 27}.
#
# Returns a data.frame with columns: date, relative_day, proportion, padded.
build_dose_obs <- function(dose_obs = DOSE_OBS,
                           start_date = min(dose_obs$date)) {
  start_date <- as.Date(start_date)
  first_obs  <- min(dose_obs$date)
  if (start_date > first_obs) {
    stop("`start_date` (", start_date, ") must be on or before the first ",
         "observation (", first_obs, ").", call. = FALSE)
  }

  obs <- data.frame(
    date         = dose_obs$date,
    relative_day = as.integer(dose_obs$date - start_date),
    proportion   = dose_obs$percentage / 100,
    padded       = FALSE,
    stringsAsFactors = FALSE
  )

  # If the start date is earlier than the first observation, add a SINGLE zero
  # at the start date (relative day 0) to anchor the curve at the 0% floor.
  if (start_date < first_obs) {
    pad <- data.frame(
      date         = start_date,
      relative_day = 0L,
      proportion   = 0,
      padded       = TRUE,
      stringsAsFactors = FALSE
    )
    obs <- rbind(pad, obs)
  }

  obs[order(obs$relative_day), , drop = FALSE]
}

# ---- Fit + extrapolate the logistic dose Q curve ---------------------------
# Shared fitting routine used by BOTH 01_fit_dose_q_curve.R (a single fit) and
# 03_compare_start_dates.R (one fit per candidate start date), so they use an
# IDENTICAL model + priors and cannot silently drift apart. It fits
# stan-models/logistic_qcurve_single.stan to the relative-day observations for a
# given `start_date` (front padding handled by build_dose_obs()), then returns
# the fitted/extrapolated Q curve on a daily grid out to `forward_days` past the
# last observation. The shape priors (t50, log k) are derived from the data
# window exactly as before; the endpoint pins (min 0, max 1) and the noise scale
# are the tunable knobs.
#
# Requires cmdstanr to be loaded by the caller. Pass a pre-compiled model in
# `mod` to avoid recompiling when fitting many start dates.
#
# Returns a list:
#   q_curve       data.frame(relative_day, date, q_mean, q_median, q_lower,
#                            q_upper, segment)
#   observations  the relative-day obs table that was fit
#   param_summary posterior summary of L, U, t50, k, sigma
#   priors_used   the prior settings actually used (incl. the data-derived t50)
#   fit           the CmdStanMCMC object
#   last_obs_day, horizon_day, start_date
fit_dose_q_curve <- function(start_date,
                             dose_obs       = DOSE_OBS,
                             forward_days   = 365L,
                             mod            = NULL,
                             lower_prior    = c(mean = 0.0, sd = 0.01),
                             upper_prior    = c(mean = 1.0, sd = 0.02),
                             logk_prior     = c(mean = log(0.2), sd = 0.7),
                             sigma_prior_sd = 0.15,
                             chains         = 4L,
                             iter_warmup    = 1500L,
                             iter_sampling  = 1500L,
                             adapt_delta    = 0.95,
                             max_treedepth  = 12L,
                             seed           = 123L,
                             refresh        = 200L) {
  start_date <- as.Date(start_date)
  obs <- build_dose_obs(dose_obs, start_date = start_date)

  last_obs_day <- max(obs$relative_day)
  horizon_day  <- last_obs_day + forward_days
  day_pred     <- 0:horizon_day

  # Weakly-informative shape priors derived from the observation window so the
  # data dominate (the data window shifts with the start date / padding).
  data_span      <- diff(range(obs$relative_day))
  t50_prior_mean <- mean(range(obs$relative_day))
  t50_prior_sd   <- max(data_span, 10)

  stan_data <- list(
    N = nrow(obs), day = as.numeric(obs$relative_day), y = as.numeric(obs$proportion),
    lower_prior_mean = lower_prior[["mean"]], lower_prior_sd = lower_prior[["sd"]],
    upper_prior_mean = upper_prior[["mean"]], upper_prior_sd = upper_prior[["sd"]],
    t50_prior_mean = t50_prior_mean, t50_prior_sd = t50_prior_sd,
    logk_prior_mean = logk_prior[["mean"]], logk_prior_sd = logk_prior[["sd"]],
    sigma_prior_sd = sigma_prior_sd,
    N_pred = length(day_pred), day_pred = as.numeric(day_pred)
  )

  if (is.null(mod)) {
    mod <- cmdstanr::cmdstan_model(file.path(DIR_STAN, "logistic_qcurve_single.stan"))
  }

  # Explicit, in-range initial values per chain. These keep early warmup away
  # from the awkward boundaries (sigma -> 0, U pinned at 1, t50/k far off) that
  # otherwise trigger sporadic "student_t scale = 0" proposal-rejection messages.
  init_fun <- function() list(
    L     = stats::runif(1, 0.00, 0.05),
    U     = stats::runif(1, 0.85, 0.99),
    t50   = stats::rnorm(1, t50_prior_mean, 5),
    log_k = stats::rnorm(1, logk_prior[["mean"]], 0.3),
    sigma = stats::runif(1, 0.05, 0.20)
  )

  fit <- mod$sample(
    data = stan_data, seed = seed,
    chains = chains, parallel_chains = chains,
    iter_warmup = iter_warmup, iter_sampling = iter_sampling,
    adapt_delta = adapt_delta, max_treedepth = max_treedepth, refresh = refresh,
    init = init_fun
  )

  # Parse q_pred[m] -> m (base R; no stringr dependency in helpers).
  qp <- fit$summary(variables = "q_pred")
  qp$grid_id <- as.integer(sub("^q_pred\\[([0-9]+)\\]$", "\\1", qp$variable))
  qp <- qp[order(qp$grid_id), , drop = FALSE]

  q_curve <- data.frame(
    relative_day = day_pred[qp$grid_id],
    date         = start_date + day_pred[qp$grid_id],
    q_mean       = qp$mean,
    q_median     = qp$median,
    q_lower      = qp$q5,
    q_upper      = qp$q95,
    segment      = ifelse(day_pred[qp$grid_id] <= last_obs_day,
                          "observed_window", "extrapolated"),
    stringsAsFactors = FALSE
  )

  list(
    q_curve       = q_curve,
    observations  = obs,
    param_summary = fit$summary(variables = c("L", "U", "t50", "k", "sigma")),
    priors_used   = list(
      lower = lower_prior, upper = upper_prior, logk = logk_prior,
      t50   = c(mean = t50_prior_mean, sd = t50_prior_sd),
      sigma = c(sd = sigma_prior_sd)
    ),
    fit          = fit,
    last_obs_day = last_obs_day,
    horizon_day  = horizon_day,
    start_date   = start_date
  )
}

# ---- Smooth a daily-incidence series + its instantaneous growth rate -------
# Used to estimate a robust growth rate from noisy daily incidence (in 02 for the
# fiber + confirmed r(t) comparison; 04 carries its own equivalent copy). Prefer
# a quasi-Poisson GAM (log link) spline of time -- it handles count noise, zeros
# and overdispersion -- and fall back to a loess on log(incidence) if mgcv is
# unavailable. `day` is numeric time; `value` the (possibly fractional) incidence.
fit_incidence_smoother <- function(day, value) {
  ok <- is.finite(day) & is.finite(value)
  df <- data.frame(day = day[ok], value = pmax(value[ok], 0))
  if (requireNamespace("mgcv", quietly = TRUE) && length(unique(df$day)) >= 6) {
    k <- max(4L, min(10L, length(unique(df$day)) %/% 3L))
    list(type = "gam",
         model = mgcv::gam(value ~ s(day, k = k), family = quasipoisson(), data = df))
  } else {
    list(type = "loess", model = stats::loess(log(value + 0.5) ~ day, data = df, span = 0.6))
  }
}
# Fitted incidence (response scale) at arbitrary days.
smooth_incidence_at <- function(fit, days) {
  if (fit$type == "gam")
    as.numeric(predict(fit$model, data.frame(day = days), type = "response"))
  else exp(as.numeric(predict(fit$model, data.frame(day = days))))
}
# Instantaneous growth rate r(t) = d/dt log(fitted incidence), by central
# difference (for the log-link GAM this is the derivative of the fitted smooth).
smooth_r_at <- function(fit, days, h = 0.5) {
  (log(smooth_incidence_at(fit, days + h)) -
   log(smooth_incidence_at(fit, days - h))) / (2 * h)
}
