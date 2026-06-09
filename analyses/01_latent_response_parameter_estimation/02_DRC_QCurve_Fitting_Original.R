# ============================================================================
# 02_DRC_QCurve_Fitting_Original.R
# ----------------------------------------------------------------------------
# Fit the DRC "conflict" and "conflict++" response curves using MODEL B (the
# bounds-only model: the response-quality curve Q is supplied as fixed data from
# the empirical SDB reconstruction, and only the per-parameter magnitude
# endpoints are estimated).
#
# Inputs : data-processed/drc_anchors.csv
#          data-processed/drc_conflict_qseries.csv
#          data-processed/drc_conflict_plusplus_qseries.csv
# Stan   : stan-models/modelB_fixedQ_boundsOnly.stan
# Output : data-processed/drc_conflict_fit.rds
#          data-processed/drc_conflict_plusplus_fit.rds
#
# NOTE on the no-conflict scenario: by design it is NOT fit here. Its data are
# still prepared in 00, and it is explored separately in DRC_no_conflict_checking.R
# (Model B vs Model A). Its horizon is used by 03 to time-stretch the West Africa
# conflict scenarios.
# ============================================================================

suppressPackageStartupMessages({
  library(cmdstanr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(tibble)
})

source("helpers.R")

set.seed(123)

# ----------------------------------------------------------------------------
# Per-parameter metadata (DRC admissible domains)
# ----------------------------------------------------------------------------
# Note these differ slightly from West Africa: the delay ceiling is 12 days and
# hospital unsafe funerals are capped much tighter at 0.010.
domain_meta <- tribble(
  ~parameter,                ~abs_min, ~abs_max,
  "delay_hosp",              0.5,      12.0,
  "p_hosp",                  0.0,      1.0,
  "p_ETU",                   0.0,      1.0,
  "latent_IPC",              0.0,      1.0,
  "p_unsafe_funeral_comm",   0.0,      1.0,
  "p_unsafe_funeral_hosp",   0.0,      0.010
)

drc_anchors <- read_csv("data-processed/drc_anchors.csv", show_col_types = FALSE)

param_meta <- drc_anchors %>%
  group_by(parameter) %>%
  summarise(
    lb_prior_mean = first(lower_bound),
    ub_prior_mean = first(upper_bound),
    direction     = first(direction),
    .groups = "drop"
  ) %>%
  left_join(domain_meta, by = "parameter") %>%
  mutate(
    param_id  = match(parameter, PARAM_LEVELS),
    increases = if_else(direction == "up", 1L, 0L)
  ) %>%
  arrange(param_id) %>%
  build_param_support()

J <- nrow(param_meta)

# Shared hyperprior for the observation scale.
SIGMA_FRAC_PRIOR_MEANLOG <- log(0.12)
SIGMA_FRAC_PRIOR_SDLOG   <- 0.60

mod <- cmdstan_model("stan-models/modelB_fixedQ_boundsOnly.stan")

# ----------------------------------------------------------------------------
# Fit one DRC scenario
# ----------------------------------------------------------------------------
# Given a Q series (the shared empirical curve) and the DRC anchors, this:
#   1. maps each anchor onto the shared Q (q_obs = Q at the anchor's day),
#   2. builds a dense prediction grid of the shared Q,
#   3. fits Model B (estimating only the endpoints), and
#   4. returns the per-parameter native-unit curves, with the community
#      unsafe-funeral curve replaced by its ABSOLUTE Warsame proxy (1 - success),
#      because that parameter is on the absolute scale, not the relative Q scale.
fit_drc_scenario <- function(qseries, label) {

  message("\n==== Fitting DRC scenario: ", label, " ====")

  # The shared horizon spans the anchors and the Q series.
  max_day <- max(c(drc_anchors$relative_day, qseries$relative_day), na.rm = TRUE)

  # The largest empirical success value sets the absolute floor of community
  # unsafe funerals (1 - max success) and is used in the tweak prior below.
  success_max <- max(qseries$success_smoothed, na.rm = TRUE)

  # ---- 1. Map each anchor onto the shared Q ----
  q_at_day <- make_interp(qseries$relative_day, qseries$q_value)
  fit_df <- drc_anchors %>%
    left_join(select(param_meta, parameter, param_id, increases), by = "parameter") %>%
    mutate(
      q_obs       = clip01(q_at_day(relative_day)),
      y_obs       = value_used,
      weight      = if_else(is.na(weight), 1, weight),
      obs_sd_mult = 1 / pmax(weight, 0.25)
    ) %>%
    arrange(param_id)

  # ---- 2. Dense prediction grid of the shared Q (regular grid + exact support) ----
  tau_grid <- sort(unique(c(seq(0, 1, length.out = 250), qseries$tau_q)))
  q_on_tau <- make_interp(qseries$tau_q, qseries$q_value)
  q_pred   <- clip01(q_on_tau(tau_grid))

  # ---- 3. Endpoint tweak priors (per parameter, Model B exposes upper/lower) ----
  # Community unsafe funerals: pin the high end (at Q=0) to 1, and the low end
  # (at Q=1) to the absolute Warsame floor 1 - max(success).
  # Hospital unsafe funerals: strongly regularise to a near-zero band.
  zero_flags <- rep(0L, J)
  tw <- list(
    use_upper = zero_flags, upper_mean = param_meta$ub_prior_mean, upper_sd = rep(1.0, J),
    use_lower = zero_flags, lower_mean = param_meta$lb_prior_mean, lower_sd = rep(1.0, J)
  )
  set_tweak <- function(tw, p, upper = NULL, lower = NULL) {
    j <- param_meta$param_id[param_meta$parameter == p]
    if (!is.null(upper)) { tw$use_upper[j] <- 1L; tw$upper_mean[j] <- upper$mean; tw$upper_sd[j] <- upper$sd }
    if (!is.null(lower)) { tw$use_lower[j] <- 1L; tw$lower_mean[j] <- lower$mean; tw$lower_sd[j] <- lower$sd }
    tw
  }
  tw <- set_tweak(tw, "p_unsafe_funeral_comm",
                  upper = list(mean = 1.00,               sd = 0.02),
                  lower = list(mean = clip01(1 - success_max), sd = 0.02))
  tw <- set_tweak(tw, "p_unsafe_funeral_hosp",
                  upper = list(mean = 0.010,  sd = 0.003),
                  lower = list(mean = 0.0005, sd = 0.001))

  stan_data <- list(
    N = nrow(fit_df), J = J, N_pred = length(tau_grid),
    param_id = fit_df$param_id, q_obs = fit_df$q_obs,
    y_obs = fit_df$y_obs, obs_sd_mult = fit_df$obs_sd_mult,
    abs_min = param_meta$abs_min, abs_max = param_meta$abs_max,
    lower_floor = param_meta$lower_floor, lower_cap = param_meta$lower_cap,
    upper_cap = param_meta$upper_cap,
    lb_prior_mean = param_meta$lb_prior_mean, lb_prior_sd = param_meta$lb_prior_sd,
    ub_prior_mean = param_meta$ub_prior_mean, ub_prior_sd = param_meta$ub_prior_sd,
    increases = param_meta$increases,
    sigma_frac_prior_meanlog = SIGMA_FRAC_PRIOR_MEANLOG,
    sigma_frac_prior_sdlog   = SIGMA_FRAC_PRIOR_SDLOG,
    use_upper_tweak = tw$use_upper, upper_tweak_mean = tw$upper_mean, upper_tweak_sd = tw$upper_sd,
    use_lower_tweak = tw$use_lower, lower_tweak_mean = tw$lower_mean, lower_tweak_sd = tw$lower_sd,
    q_pred = q_pred
  )

  fit <- mod$sample(
    data = stan_data, seed = 123,
    chains = 4, parallel_chains = 4,
    iter_warmup = 1500, iter_sampling = 1500,
    adapt_delta = 0.98, max_treedepth = 13, refresh = 200
  )
  cat("\nSampler diagnostics (", label, "):\n", sep = ""); print(fit$diagnostic_summary())

  # ---- 4. Per-parameter native-unit curves ----
  curve_summ <- fit$summary(variables = "theta_pred") %>%
    mutate(
      param_id = as.integer(str_match(variable, "\\[([0-9]+),([0-9]+)\\]")[, 2]),
      grid_id  = as.integer(str_match(variable, "\\[([0-9]+),([0-9]+)\\]")[, 3]),
      tau          = tau_grid[grid_id],
      relative_day = tau * max_day
    ) %>%
    left_join(select(param_meta, param_id, parameter), by = "param_id") %>%
    select(parameter, param_id, tau, relative_day, mean, q5, q95)

  # Deterministic override: community unsafe funerals on the ABSOLUTE scale.
  ufc_on_tau <- make_interp(qseries$tau_q, qseries$unsafe_funeral_comm_proxy)
  curve_summ <- curve_summ %>%
    mutate(mean = if_else(parameter == "p_unsafe_funeral_comm",
                          clip01(ufc_on_tau(tau)), mean))

  # The shared Q on the prediction grid (used as the scenario's q_value in 03).
  q_grid <- tibble(tau = tau_grid, relative_day = tau_grid * max_day, q_value = q_pred)

  list(curve_summ = curve_summ, q_grid = q_grid, param_meta = param_meta, max_day = max_day)
}

# ----------------------------------------------------------------------------
# Run both conflict scenarios and save
# ----------------------------------------------------------------------------
drc_conflict_qseries          <- read_csv("data-processed/drc_conflict_qseries.csv", show_col_types = FALSE)
drc_conflict_plusplus_qseries <- read_csv("data-processed/drc_conflict_plusplus_qseries.csv", show_col_types = FALSE)

drc_conflict_fit          <- fit_drc_scenario(drc_conflict_qseries,          "drc_conflict")
drc_conflict_plusplus_fit <- fit_drc_scenario(drc_conflict_plusplus_qseries, "drc_conflict_plusplus")

saveRDS(drc_conflict_fit,          "data-processed/drc_conflict_fit.rds")
saveRDS(drc_conflict_plusplus_fit, "data-processed/drc_conflict_plusplus_fit.rds")

message("\n02_DRC_QCurve_Fitting_Original.R complete. Saved DRC conflict + conflict++ fits.")
