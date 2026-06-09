# ============================================================================
# 01_WestAfrica_QCurve_Fitting_Original.R
# ----------------------------------------------------------------------------
# Fit the West Africa worst-case response curves using MODEL A (the
# partial-pooled model that ESTIMATES both the curve shape and the magnitude
# endpoints), WITH the targeted "tweak" priors switched on.
#
# Input  : data-processed/wa_anchors.csv   (from 00)
# Stan   : stan-models/modelA_partialpool_estimateQ_withTweaks.stan
# Output : data-processed/wa_fit.rds        (curves + endpoints, consumed by 03)
#          data-processed/wa_fit_curve_summaries.csv  (human-readable)
#
# The tweaks are isolated in one clearly-labelled block below so they are easy to
# audit. west_africa_checking.R fits this same data with and without the tweaks
# and overlays the curves so their effect can be seen directly.
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
# 1. Read the cleaned West Africa anchors
# ----------------------------------------------------------------------------
wa_anchors <- read_csv("data-processed/wa_anchors.csv", show_col_types = FALSE)

# ----------------------------------------------------------------------------
# 2. Per-parameter metadata: priors, direction, and admissible support
# ----------------------------------------------------------------------------
# Hard admissible domains for the West Africa parameters. These are the physical
# limits within which each parameter's estimated endpoints must lie. (Note the
# delay ceiling is 10 days and hospital unsafe funerals are capped at 0.10 here.)
domain_meta <- tribble(
  ~parameter,                ~abs_min, ~abs_max,
  "delay_hosp",              0.5,      10.0,
  "p_hosp",                  0.0,      1.0,
  "p_ETU",                   0.0,      1.0,
  "latent_IPC",              0.0,      1.0,
  "p_unsafe_funeral_comm",   0.0,      1.0,
  "p_unsafe_funeral_hosp",   0.0,      0.10
)

# One row per parameter: take the workbook lower/upper bounds as the prior
# centres for the endpoints, and the direction (up = improves upward).
param_meta <- wa_anchors %>%
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
  arrange(param_id)

if (!all(PARAM_LEVELS %in% param_meta$parameter)) {
  stop("West Africa anchors are missing parameters: ",
       paste(setdiff(PARAM_LEVELS, param_meta$parameter), collapse = ", "))
}

# Add the endpoint transform limits and prior SDs (shared logic with Model B).
param_meta <- build_param_support(param_meta)
J <- nrow(param_meta)

# ----------------------------------------------------------------------------
# 3. Assemble the per-observation fit table
# ----------------------------------------------------------------------------
# Normalise outbreak time to tau in [0,1] over the West Africa anchor window, and
# turn each anchor's weight into a noise multiplier (higher weight = tighter fit).
max_day <- max(wa_anchors$relative_day, na.rm = TRUE)

fit_df <- wa_anchors %>%
  left_join(select(param_meta, parameter, param_id, increases), by = "parameter") %>%
  mutate(
    tau         = relative_day / max_day,
    y_obs       = value_used,
    weight      = if_else(is.na(weight), 1, weight),
    obs_sd_mult = 1 / pmax(weight, 0.25)
  ) %>%
  arrange(param_id, tau)

# ----------------------------------------------------------------------------
# 4. TARGETED TWEAK PRIORS  <-- the part to scrutinise
# ----------------------------------------------------------------------------
# Every tweak below targets ONE parameter: community unsafe funerals
# (p_unsafe_funeral_comm). They encode an external belief that, in West Africa,
# community unsafe funerals stayed near-universal early in the response and then
# fell sharply and relatively late. Each is an extra informative prior on top of
# the hierarchical model, and each should be justifiable from evidence.
#
#   - t50  -> 0.72 : the decline's midpoint sits ~72% through the response window
#   - logk -> log(11) : a steep (rather than gradual) decline
#   - upper -> 0.95 : the early level is pinned near the ceiling
#
# In addition, one anchor (WA_UFC_01) that is discordant with this shape is
# DOWN-WEIGHTED (its observation SD is inflated) so it is fit more loosely.
#
# To reproduce the clean no-tweak baseline, set TWEAKS_ON <- FALSE (or use
# west_africa_checking.R, which fits both).
TWEAKS_ON <- TRUE

tweak_param <- "p_unsafe_funeral_comm"

# Tweak values (prior mean, prior sd) for the targeted parameter.
t50_tweak   <- list(mean = 0.72,     sd = 0.07)
logk_tweak  <- list(mean = log(11),  sd = 0.35)
upper_tweak <- list(mean = 0.95,     sd = 0.04)

# Discordant-anchor down-weighting.
anchor_to_downweight      <- "WA_UFC_01"
anchor_downweight_sd_mult <- 1.75   # > 1 loosens the fit to that single anchor

if (TWEAKS_ON) {
  fit_df <- fit_df %>%
    mutate(obs_sd_mult = if_else(
      anchor_id == anchor_to_downweight & parameter == tweak_param,
      obs_sd_mult * anchor_downweight_sd_mult,
      obs_sd_mult
    ))
}

# Expand the tweak settings into the per-parameter vectors the Stan model expects.
# Default everything OFF; switch on only the targeted parameter's tweaks.
j_tweak <- param_meta$param_id[param_meta$parameter == tweak_param]

zero_flags <- rep(0L, J)
tweak <- list(
  use_t50   = zero_flags, t50_mean   = rep(0.5, J),                  t50_sd   = rep(1.0, J),
  use_logk  = zero_flags, logk_mean  = rep(log(8), J),               logk_sd  = rep(1.0, J),
  use_upper = zero_flags, upper_mean = param_meta$ub_prior_mean,     upper_sd = rep(1.0, J),
  use_lower = zero_flags, lower_mean = param_meta$lb_prior_mean,     lower_sd = rep(1.0, J)
)

if (TWEAKS_ON) {
  tweak$use_t50[j_tweak]    <- 1L; tweak$t50_mean[j_tweak]   <- t50_tweak$mean;   tweak$t50_sd[j_tweak]   <- t50_tweak$sd
  tweak$use_logk[j_tweak]   <- 1L; tweak$logk_mean[j_tweak]  <- logk_tweak$mean;  tweak$logk_sd[j_tweak]  <- logk_tweak$sd
  tweak$use_upper[j_tweak]  <- 1L; tweak$upper_mean[j_tweak] <- upper_tweak$mean; tweak$upper_sd[j_tweak] <- upper_tweak$sd
  # The lower-endpoint tweak hook exists in the model but is unused for West Africa.
}

# ----------------------------------------------------------------------------
# 5. Hyperpriors and the prediction grid
# ----------------------------------------------------------------------------
# Hyperpriors for the partially-pooled curve shape and the observation scale.
hyper <- list(
  mu_t50_raw_prior_mean = 0.0,  mu_t50_raw_prior_sd = 1.25, sigma_t50_prior_sd  = 0.50,
  mu_log_k_prior_mean   = log(8), mu_log_k_prior_sd = 0.80, sigma_log_k_prior_sd = 0.40,
  sigma_frac_prior_meanlog = log(0.12), sigma_frac_prior_sdlog = 0.60
)

# Predict the curves on a regular 100-point grid over normalised time.
tau_pred <- seq(0, 1, length.out = 100)

# ----------------------------------------------------------------------------
# 6. Assemble the Stan data list
# ----------------------------------------------------------------------------
stan_data <- c(
  list(
    N = nrow(fit_df), J = J,
    param_id = fit_df$param_id, tau = fit_df$tau,
    y_obs = fit_df$y_obs, obs_sd_mult = fit_df$obs_sd_mult,
    abs_min = param_meta$abs_min, abs_max = param_meta$abs_max,
    lower_floor = param_meta$lower_floor, lower_cap = param_meta$lower_cap,
    upper_cap = param_meta$upper_cap,
    lb_prior_mean = param_meta$lb_prior_mean, lb_prior_sd = param_meta$lb_prior_sd,
    ub_prior_mean = param_meta$ub_prior_mean, ub_prior_sd = param_meta$ub_prior_sd,
    increases = param_meta$increases,
    use_t50_tweak = tweak$use_t50,   t50_tweak_mean = tweak$t50_mean,   t50_tweak_sd = tweak$t50_sd,
    use_logk_tweak = tweak$use_logk, logk_tweak_mean = tweak$logk_mean, logk_tweak_sd = tweak$logk_sd,
    use_upper_tweak = tweak$use_upper, upper_tweak_mean = tweak$upper_mean, upper_tweak_sd = tweak$upper_sd,
    use_lower_tweak = tweak$use_lower, lower_tweak_mean = tweak$lower_mean, lower_tweak_sd = tweak$lower_sd,
    N_pred = length(tau_pred), tau_pred = tau_pred
  ),
  hyper
)

# ----------------------------------------------------------------------------
# 7. Compile and sample
# ----------------------------------------------------------------------------
mod <- cmdstan_model("stan-models/modelA_partialpool_estimateQ_withTweaks.stan")

fit <- mod$sample(
  data = stan_data, seed = 123,
  chains = 4, parallel_chains = 4,
  iter_warmup = 1500, iter_sampling = 1500,
  adapt_delta = 0.97, max_treedepth = 13, refresh = 200
)

cat("\nSampler diagnostics:\n"); print(fit$diagnostic_summary())

# ----------------------------------------------------------------------------
# 8. Extract tidy posterior-mean curves
# ----------------------------------------------------------------------------
# Pull a [J, N_pred] matrix variable (Q_pred or theta_pred) out of the fit and
# return a tidy table with one row per (parameter, grid point).
tidy_matrix <- function(fit, varname) {
  fit$summary(variables = varname) %>%
    mutate(
      param_id = as.integer(str_match(variable, "\\[([0-9]+),([0-9]+)\\]")[, 2]),
      grid_id  = as.integer(str_match(variable, "\\[([0-9]+),([0-9]+)\\]")[, 3]),
      tau          = tau_pred[grid_id],
      relative_day = tau * max_day
    ) %>%
    left_join(select(param_meta, param_id, parameter), by = "param_id") %>%
    select(parameter, param_id, grid_id, tau, relative_day, mean, median, sd, q5, q95)
}

q_summ     <- tidy_matrix(fit, "Q_pred")      # per-parameter response-quality curve Q_j(tau)
curve_summ <- tidy_matrix(fit, "theta_pred")  # per-parameter value curve theta_j(tau)

# ----------------------------------------------------------------------------
# 9. Save
# ----------------------------------------------------------------------------
# wa_fit.rds carries everything script 03 needs to (a) report the standalone
# West Africa scenario and (b) build the West-Africa-with-conflict scenarios.
wa_fit <- list(
  q_summ      = q_summ,        # Q_j(tau): per-parameter shape, used to modulate by DRC conflict Q
  curve_summ  = curve_summ,    # theta_j(tau): per-parameter native-unit trajectory
  param_meta  = param_meta,
  tau_pred    = tau_pred,
  max_day     = max_day,
  tweaks_on   = TWEAKS_ON
)
saveRDS(wa_fit, "data-processed/wa_fit.rds")
write_csv(curve_summ, "data-processed/wa_fit_curve_summaries.csv")

message("\n01_WestAfrica_QCurve_Fitting_Original.R complete. Saved data-processed/wa_fit.rds")
