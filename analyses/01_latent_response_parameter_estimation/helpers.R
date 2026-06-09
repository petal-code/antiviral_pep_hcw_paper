# ============================================================================
# helpers.R
# ----------------------------------------------------------------------------
# Small shared utilities used across the latent-response-parameter pipeline.
# Keeping them here (rather than copy-pasting into every script) means there is
# exactly ONE definition of each, so they cannot silently drift apart.
#
#   source("helpers.R")
#
# The heavy, scenario-specific data wrangling (the SDB reconstruction) lives in
# 00_DataPreparation_and_Cleaning.R, not here - this file is only generic glue.
# ============================================================================

# ---- The shared model "vocabulary" -----------------------------------------

# The six latent response parameters, in a fixed canonical order. Every script
# uses this ordering so that parameter id 1..6 means the same thing everywhere.
PARAM_LEVELS <- c(
  "delay_hosp",            # mean delay (days) from symptom onset to hospitalisation
  "p_hosp",                # probability an infected person is hospitalised
  "p_ETU",                 # proportion of hospitalised cases managed in an ETU/ETC
  "latent_IPC",            # latent infection-prevention-and-control / PPE index
  "p_unsafe_funeral_comm", # probability of an unsafe funeral after a community death
  "p_unsafe_funeral_hosp"  # probability of an unsafe funeral after a hospital death
)

# Human-readable panel titles (used when plotting), keyed by parameter.
PANEL_LOOKUP <- c(
  delay_hosp            = "Delay to hospitalisation",
  p_hosp                = "Probability of hospitalisation",
  p_ETU                 = "Proportion in ETU / ETC",
  latent_IPC            = "Latent IPC index",
  p_unsafe_funeral_comm = "Unsafe funeral after community death",
  p_unsafe_funeral_hosp = "Unsafe funeral after hospital death"
)

# The final scenario matrices are all reported on a common 0..730 day horizon
# (731 daily rows), with tau = relative_day / HORIZON_DAYS.
HORIZON_DAYS <- 730L

# ---- Tiny numeric helpers ---------------------------------------------------

# Clamp a value (or vector) into the closed unit interval [0, 1].
clip01 <- function(x) pmin(1, pmax(0, x))

# Min-max rescale to [0, 1]. Used only for the diagnostic q_value column.
rescale_01 <- function(x) {
  r <- range(x, na.rm = TRUE)
  if (!all(is.finite(r))) stop("rescale_01(): non-finite range.")
  if (diff(r) <= 0)       stop("rescale_01(): zero-width range.")
  (x - r[1]) / diff(r)
}

# ---- Locating raw input files ----------------------------------------------

# Resolve a raw input workbook by name inside data-raw/. Accepts several
# candidate filenames (the SDB workbook has been distributed under a couple of
# slightly different names) and returns the first that exists.
resolve_input_file <- function(candidates, description = "input file",
                               data_raw_dir = "data-raw") {
  for (nm in candidates) {
    p <- file.path(data_raw_dir, nm)
    if (file.exists(p)) return(p)
  }
  stop(
    "Could not find ", description, " in '", data_raw_dir, "'. Looked for: ",
    paste(candidates, collapse = ", ")
  )
}

# ---- Interpolation ----------------------------------------------------------

# Build a linear interpolation function from (x, y) support points. Extrapolates
# by holding the nearest endpoint flat (rule = 2), which is what we want for the
# response curves (before the first / after the last support point the value
# simply holds). Duplicated x values are averaged.
make_interp <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  stats::approxfun(x[ok], y[ok], method = "linear", rule = 2, ties = mean)
}

# Centred rolling mean over an ordered vector (used to smooth the weekly SDB
# reconstruction). For an even window of width k (default 4) the window is
# necessarily slightly asymmetric: floor((k-1)/2) points before, the current
# point, and the remainder after. Positions with no finite values stay NA and
# are filled by the caller.
rolling_mean_centered <- function(y, k = 4) {
  n <- length(y)
  if (n == 0) return(numeric(0))
  if (k < 1) stop("rolling window k must be >= 1.")
  out <- rep(NA_real_, n)
  left_n  <- floor((k - 1) / 2)
  right_n <- (k - 1) - left_n
  for (i in seq_len(n)) {
    idx <- seq.int(max(1, i - left_n), min(n, i + right_n))
    yy  <- y[idx]
    ok  <- is.finite(yy)
    if (any(ok)) out[i] <- mean(yy[ok])
  }
  out
}

# ---- Shared construction of the per-parameter endpoint supports ------------

# Both Model A and Model B estimate per-parameter lower/upper magnitude
# endpoints inside an admissible window, with normal priors centred on the
# workbook bounds. This helper builds that support table (the lower_floor /
# lower_cap / upper_cap transform limits and the prior SDs) identically for both
# models, so the two fits use exactly the same endpoint geometry.
#
# Inputs:
#   param_meta  : a data.frame/tibble with columns
#                   parameter, lb_prior_mean, ub_prior_mean, abs_min, abs_max
#   support_half_width_mult, bound_sd_frac_of_initial_span, bound_sd_frac_of_domain
#                 : the three tuning constants (defaults match the originals)
# Returns param_meta with span0, domain_w, lower_floor, lower_cap, upper_cap,
# lb_prior_sd, ub_prior_sd added.
build_param_support <- function(param_meta,
                                support_half_width_mult       = 0.75,
                                bound_sd_frac_of_initial_span = 0.15,
                                bound_sd_frac_of_domain       = 0.03) {

  pm <- within(param_meta, {
    span0    <- ub_prior_mean - lb_prior_mean   # prior width of each parameter
    domain_w <- abs_max - abs_min               # full admissible width
  })

  # Lower endpoint may sit anywhere in [lower_floor, lower_cap]; upper endpoint
  # anywhere in [lower_est, upper_cap]. The caps keep lower below upper and keep
  # both inside the hard admissible domain.
  pm$lower_floor <- pmax(pm$abs_min, pm$lb_prior_mean - support_half_width_mult * pm$span0)

  lower_cap_raw <- pmin(
    pm$lb_prior_mean + support_half_width_mult * pm$span0,
    pm$ub_prior_mean - 0.05 * pm$span0,
    pm$abs_max       - 0.05 * pm$domain_w
  )
  pm$lower_cap <- pmax(pm$lower_floor + 1e-3, lower_cap_raw)

  upper_cap_raw <- pmin(pm$abs_max, pm$ub_prior_mean + support_half_width_mult * pm$span0)
  pm$upper_cap  <- pmax(pm$lower_cap + 1e-3, upper_cap_raw)

  # Prior SDs on the endpoints: a fraction of the prior span or of the domain,
  # whichever is larger, with a tiny floor to avoid degenerate priors.
  pm$lb_prior_sd <- pmax(bound_sd_frac_of_initial_span * pm$span0,
                         bound_sd_frac_of_domain * pm$domain_w, 1e-3)
  pm$ub_prior_sd <- pmax(bound_sd_frac_of_initial_span * pm$span0,
                         bound_sd_frac_of_domain * pm$domain_w, 1e-3)

  # Sanity checks: the admissible windows must be non-empty and ordered.
  if (any(pm$lower_cap <= pm$lower_floor)) stop("build_param_support(): lower_cap <= lower_floor.")
  if (any(pm$upper_cap <= pm$lower_cap))   stop("build_param_support(): upper_cap <= lower_cap.")
  pm
}

# ---- Mapping a fitted tau-curve onto the reported 0..730 day grid ----------

# A fit produces a parameter's mean trajectory on a normalised tau in [0, 1].
# The final scenario matrices are reported on integer days 0..730 with
# tau = day / 730. This helper linearly interpolates a (tau_in, value_in) curve
# onto that daily grid and returns a tibble with day, relative_day, tau, value.
curve_to_day_grid <- function(tau_in, value_in, horizon_days = HORIZON_DAYS) {
  f <- make_interp(tau_in, value_in)
  days <- 0:horizon_days
  tau_out <- days / horizon_days
  data.frame(
    time_index   = days + 1L,
    relative_day = days,
    tau          = tau_out,
    value        = as.numeric(f(tau_out))
  )
}
