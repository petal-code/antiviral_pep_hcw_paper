# ============================================================================
# 01_WestAfrica_QCurve_Fitting_Original.R
# ----------------------------------------------------------------------------
# WHAT THIS SCRIPT DOES
#   Fits the six West Africa "worst-case" response-parameter curves using
#   MODEL A, then saves the fitted curves for the combine step (03).
#
# THE MODEL IN ONE PARAGRAPH
#   Outbreak progress is measured on a normalised clock tau in [0, 1]. Each
#   response parameter j (probability of hospitalisation, delay to
#   hospitalisation, ...) is assumed to move along a smooth S-shaped
#   "response-quality" curve Q_j(tau) that runs from 0 (worst response) to 1
#   (best response), scaled between that parameter's own low and high magnitude
#   endpoints. MODEL A ESTIMATES BOTH (a) the shape of each Q_j (its midpoint t50
#   and steepness k, partially pooled across parameters) AND (b) the two
#   magnitude endpoints per parameter. The only data are a handful of literature
#   "anchor" points per parameter. (Contrast Model B, used for DRC, where the Q
#   curve is supplied as fixed data instead of being estimated.)
#
# THE "TWEAKS"
#   This version switches ON a set of targeted, informative priors ("tweaks")
#   that nudge ONE parameter (community unsafe funerals) toward an externally
#   believed shape. They are isolated in section 4 so they are easy to audit.
#   west_africa_checking.R fits the same data with and without them and overlays
#   the curves, so their effect is visible.
#
# Input  : data-processed/wa_prep.rds   (from 00; uses its $anchors element)
# Stan   : stan-models/modelA_partialpool_estimateQ_withTweaks.stan
# Output : data-processed/wa_fit.rds              (curves + endpoints; read by 03)
#          data-processed/wa_fit_curve_summaries.csv   (human-readable curves)
# ============================================================================

suppressPackageStartupMessages({
  library(cmdstanr)   # interface to CmdStan (compiles + samples the Stan model)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(tibble)
})

# helpers.R defines the path constants (DIR_PROCESSED, DIR_STAN), the canonical
# parameter order (PARAM_LEVELS), small utilities (clip01, make_interp, ...) and
# build_param_support() used below.
source(here::here("analyses", "01_latent_response_parameter_estimation", "helpers.R"))

set.seed(123)  # reproducible R-side randomness (the sampler seed is set separately)

# ----------------------------------------------------------------------------
# 1. Read the cleaned West Africa anchors
# ----------------------------------------------------------------------------
# wa_prep.rds was produced by 00. Its $anchors element is the cleaned anchor
# table (one row per literature data point; see read_anchor_sheet() in 00 for the
# column meanings, e.g. parameter, relative_day, value_used, weight, direction,
# lower_bound, upper_bound).
wa_anchors <- readRDS(file.path(DIR_PROCESSED, "wa_prep.rds"))$anchors

# ----------------------------------------------------------------------------
# 2. Per-parameter metadata: priors, direction, and admissible support
# ----------------------------------------------------------------------------
# domain_meta gives the HARD admissible domain [abs_min, abs_max] of each
# parameter: the physically sensible limits its estimated endpoints may never
# leave. (Note these are West-Africa-specific: the delay-to-hospitalisation
# ceiling is 10 days, and hospital unsafe funerals are capped at 0.10.)
domain_meta <- tribble(
  ~parameter,                ~abs_min, ~abs_max,
  "delay_hosp",              0.5,      10.0,
  "p_hosp",                  0.0,      1.0,
  "p_ETU",                   0.0,      1.0,
  "latent_IPC",              0.0,      1.0,
  "p_unsafe_funeral_comm",   0.0,      1.0,
  "p_unsafe_funeral_hosp",   0.0,      0.10
)

# Collapse the anchor table to ONE row per parameter. The workbook lower/upper
# bounds become the PRIOR CENTRES for the two magnitude endpoints, and the
# direction tells us whether the parameter improves upward or downward.
param_meta <- wa_anchors %>%
  group_by(parameter) %>%
  summarise(
    lb_prior_mean = first(lower_bound),   # prior centre for the LOWER endpoint
    ub_prior_mean = first(upper_bound),   # prior centre for the UPPER endpoint
    direction     = first(direction),     # "up" or "down"
    .groups = "drop"
  ) %>%
  left_join(domain_meta, by = "parameter") %>%
  mutate(
    param_id  = match(parameter, PARAM_LEVELS),     # integer id 1..6, canonical order
    increases = if_else(direction == "up", 1L, 0L)  # 1 = improves upward (Stan flag)
  ) %>%
  arrange(param_id)

# Guard: every modelled parameter must be present, or the fit is mis-specified.
if (!all(PARAM_LEVELS %in% param_meta$parameter)) {
  stop("West Africa anchors are missing parameters: ",
       paste(setdiff(PARAM_LEVELS, param_meta$parameter), collapse = ", "))
}

# build_param_support() (in helpers.R) adds the endpoint transform limits and the
# prior SDs, identically to the DRC fit. After this call, param_meta has one row
# per parameter with the following columns ------------------------------------
#   parameter      parameter name
#   param_id       integer id 1..6 in canonical PARAM_LEVELS order
#   direction      "up" / "down"   (how the parameter improves over the response)
#   increases      1 if direction == "up", else 0   (passed to the Stan model)
#   abs_min/abs_max  hard admissible domain (physical limits)
#   lb_prior_mean  prior CENTRE for the lower magnitude endpoint (workbook low end)
#   ub_prior_mean  prior CENTRE for the upper magnitude endpoint (workbook high end)
#   span0          ub_prior_mean - lb_prior_mean   (prior width of the parameter)
#   domain_w       abs_max - abs_min               (full admissible width)
#   lower_floor    lowest value the fitted LOWER endpoint may take
#   lower_cap      highest value the fitted LOWER endpoint may take
#   upper_cap      highest value the fitted UPPER endpoint may take
#                  (the fitted upper endpoint lives in [lower endpoint, upper_cap])
#   lb_prior_sd    prior SD on the lower endpoint
#   ub_prior_sd    prior SD on the upper endpoint
param_meta <- build_param_support(param_meta)
J <- nrow(param_meta)   # number of parameters (6)

# ----------------------------------------------------------------------------
# 3. Assemble the per-observation fit table
# ----------------------------------------------------------------------------
# Stan needs one row per anchor OBSERVATION. We (a) put every anchor on the
# normalised clock tau = relative_day / max_day, and (b) turn each anchor's
# "weight" into an observation-noise multiplier: higher weight -> smaller noise
# -> the curve is pulled more tightly through that point.
max_day <- max(wa_anchors$relative_day, na.rm = TRUE)   # the West Africa anchor horizon

# fit_df: one row per anchor observation fed to Stan, with columns -------------
#   param_id      which parameter the observation belongs to (1..6)
#   tau           normalised outbreak time of the observation, in [0, 1]
#   y_obs         the observed parameter value at that time (the anchor value)
#   obs_sd_mult   per-observation noise multiplier (= 1 / weight, floored);
#                 larger -> looser fit to that point
#   (anchor_id, parameter, etc. are carried through for bookkeeping/plots)
fit_df <- wa_anchors %>%
  left_join(select(param_meta, parameter, param_id, increases), by = "parameter") %>%
  mutate(
    tau         = relative_day / max_day,
    y_obs       = value_used,
    weight      = if_else(is.na(weight), 1, weight),
    obs_sd_mult = 1 / pmax(weight, 0.25)   # floor weight at 0.25 so noise can't explode
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
# west_africa_checking.R, which fits both and overlays them).
TWEAKS_ON <- TRUE

tweak_param <- "p_unsafe_funeral_comm"   # the single parameter the tweaks target

# Tweak values, as (prior mean, prior sd) on the targeted parameter's shape/level.
t50_tweak   <- list(mean = 0.72,     sd = 0.07)   # later midpoint
logk_tweak  <- list(mean = log(11),  sd = 0.35)   # steeper decline
upper_tweak <- list(mean = 0.95,     sd = 0.04)   # high early level (near ceiling)

# Discordant-anchor down-weighting (inflate this one anchor's observation SD).
anchor_to_downweight      <- "WA_UFC_01"
anchor_downweight_sd_mult <- 1.75   # > 1 loosens the fit to that single anchor

# Apply the down-weight to the matching row of fit_df (only when tweaks are on).
if (TWEAKS_ON) {
  fit_df <- fit_df %>%
    mutate(obs_sd_mult = if_else(
      anchor_id == anchor_to_downweight & parameter == tweak_param,
      obs_sd_mult * anchor_downweight_sd_mult,
      obs_sd_mult
    ))
}

# The Stan model expects the tweak settings as length-J vectors (one entry per
# parameter), with a 0/1 switch plus a (mean, sd) for each of t50, logk, upper,
# lower. Start with EVERYTHING OFF (switches = 0, harmless placeholder means/sds),
# then turn on only the targeted parameter's tweaks below.
j_tweak <- param_meta$param_id[param_meta$parameter == tweak_param]   # row to switch on

zero_flags <- rep(0L, J)
tweak <- list(
  use_t50   = zero_flags, t50_mean   = rep(0.5, J),                  t50_sd   = rep(1.0, J),
  use_logk  = zero_flags, logk_mean  = rep(log(8), J),               logk_sd  = rep(1.0, J),
  use_upper = zero_flags, upper_mean = param_meta$ub_prior_mean,     upper_sd = rep(1.0, J),
  use_lower = zero_flags, lower_mean = param_meta$lb_prior_mean,     lower_sd = rep(1.0, J)
)

# Switch on the three tweaks for the targeted parameter only.
if (TWEAKS_ON) {
  tweak$use_t50[j_tweak]    <- 1L; tweak$t50_mean[j_tweak]   <- t50_tweak$mean;   tweak$t50_sd[j_tweak]   <- t50_tweak$sd
  tweak$use_logk[j_tweak]   <- 1L; tweak$logk_mean[j_tweak]  <- logk_tweak$mean;  tweak$logk_sd[j_tweak]  <- logk_tweak$sd
  tweak$use_upper[j_tweak]  <- 1L; tweak$upper_mean[j_tweak] <- upper_tweak$mean; tweak$upper_sd[j_tweak] <- upper_tweak$sd
  # The lower-endpoint tweak hook exists in the model but is unused for West Africa.
}

# ----------------------------------------------------------------------------
# 5. Hyperpriors and the prediction grid
# ----------------------------------------------------------------------------
# These hyperpriors govern the PARTIAL POOLING of the curve shape across the six
# parameters, and the observation-noise scale:
#   mu_t50_raw_*  : group-mean midpoint (on the logit scale) and its prior spread
#   sigma_t50_*   : how much the parameters' midpoints may differ (pooling strength)
#   mu_log_k_*    : group-mean log-steepness and its prior spread
#   sigma_log_k_* : how much steepness may differ across parameters
#   sigma_frac_*  : lognormal prior on the per-parameter observation noise fraction
hyper <- list(
  mu_t50_raw_prior_mean = 0.0,  mu_t50_raw_prior_sd = 1.25, sigma_t50_prior_sd  = 0.50,
  mu_log_k_prior_mean   = log(8), mu_log_k_prior_sd = 0.80, sigma_log_k_prior_sd = 0.40,
  sigma_frac_prior_meanlog = log(0.12), sigma_frac_prior_sdlog = 0.60
)

# Grid of normalised times at which to output the fitted curves (100 points,
# evenly spaced over the response window). 03 later stretches this onto 0..730 days.
tau_pred <- seq(0, 1, length.out = 100)

# ----------------------------------------------------------------------------
# 6. Assemble the Stan data list
# ----------------------------------------------------------------------------
# This list must contain EXACTLY the variables declared in the model's `data`
# block (sizes, the observations, the per-parameter support/priors, the tweak
# vectors, and the prediction grid). It is concatenated with `hyper` above.
stan_data <- c(
  list(
    N = nrow(fit_df), J = J,                                  # sizes
    param_id = fit_df$param_id, tau = fit_df$tau,             # observation indexing + time
    y_obs = fit_df$y_obs, obs_sd_mult = fit_df$obs_sd_mult,   # observed values + noise mult
    abs_min = param_meta$abs_min, abs_max = param_meta$abs_max,                 # hard domain
    lower_floor = param_meta$lower_floor, lower_cap = param_meta$lower_cap,     # endpoint limits
    upper_cap = param_meta$upper_cap,
    lb_prior_mean = param_meta$lb_prior_mean, lb_prior_sd = param_meta$lb_prior_sd,  # endpoint priors
    ub_prior_mean = param_meta$ub_prior_mean, ub_prior_sd = param_meta$ub_prior_sd,
    increases = param_meta$increases,                         # direction flags
    # The four optional tweak-prior blocks (all OFF except the targeted parameter):
    use_t50_tweak = tweak$use_t50,   t50_tweak_mean = tweak$t50_mean,   t50_tweak_sd = tweak$t50_sd,
    use_logk_tweak = tweak$use_logk, logk_tweak_mean = tweak$logk_mean, logk_tweak_sd = tweak$logk_sd,
    use_upper_tweak = tweak$use_upper, upper_tweak_mean = tweak$upper_mean, upper_tweak_sd = tweak$upper_sd,
    use_lower_tweak = tweak$use_lower, lower_tweak_mean = tweak$lower_mean, lower_tweak_sd = tweak$lower_sd,
    N_pred = length(tau_pred), tau_pred = tau_pred           # prediction grid
  ),
  hyper
)

# ----------------------------------------------------------------------------
# 7. Compile and sample
# ----------------------------------------------------------------------------
# Compile the Stan model (cached after the first run) and draw posterior samples.
# adapt_delta/max_treedepth are raised above their defaults because the
# hierarchical geometry can be awkward; 4 chains x 1500 post-warmup draws.
mod <- cmdstan_model(file.path(DIR_STAN, "modelA_partialpool_estimateQ_withTweaks.stan"))

fit <- mod$sample(
  data = stan_data, seed = 123,
  chains = 4, parallel_chains = 4,
  iter_warmup = 1500, iter_sampling = 1500,
  adapt_delta = 0.97, max_treedepth = 13, refresh = 200
)

# Always glance at the diagnostics (divergences, tree-depth hits, E-BFMI).
cat("\nSampler diagnostics:\n"); print(fit$diagnostic_summary())

# ----------------------------------------------------------------------------
# 8. Extract tidy posterior-mean curves
# ----------------------------------------------------------------------------
# The model returns two [J, N_pred] matrices: Q_pred (the response-quality curve
# per parameter) and theta_pred (the parameter's value per grid point). cmdstanr
# reports them with names like "theta_pred[3,57]"; this helper parses the [j,m]
# indices back out and attaches the parameter name, time, and day.
tidy_matrix <- function(fit, varname) {
  fit$summary(variables = varname) %>%
    mutate(
      param_id = as.integer(str_match(variable, "\\[([0-9]+),([0-9]+)\\]")[, 2]),  # j
      grid_id  = as.integer(str_match(variable, "\\[([0-9]+),([0-9]+)\\]")[, 3]),  # m
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
# West Africa scenario and (b) build the West-Africa-with-conflict scenarios (it
# uses q_summ to modulate by the DRC conflict curve, and curve_summ for the
# native-unit start/end magnitudes).
wa_fit <- list(
  q_summ      = q_summ,        # Q_j(tau): per-parameter shape (cols: parameter, tau, relative_day, mean, q5, q95, ...)
  curve_summ  = curve_summ,    # theta_j(tau): per-parameter native-unit trajectory (same columns)
  param_meta  = param_meta,    # the per-parameter metadata table documented in section 2
  tau_pred    = tau_pred,      # the 100-point normalised-time grid the curves are on
  max_day     = max_day,       # the West Africa anchor horizon (days) behind tau
  tweaks_on   = TWEAKS_ON      # whether the targeted tweaks were applied
)
saveRDS(wa_fit, file.path(DIR_PROCESSED, "wa_fit.rds"))

# Also write the fitted parameter curves as a plain CSV for quick human inspection.
write_csv(curve_summ, file.path(DIR_PROCESSED, "wa_fit_curve_summaries.csv"))

message("\n01_WestAfrica_QCurve_Fitting_Original.R complete. Saved data-processed/wa_fit.rds")
