# ============================================================================
# 01_WestAfrica_QCurve_Fitting(revisedMethod).R
# ----------------------------------------------------------------------------
# WHAT THIS SCRIPT DOES
#   Fits the six West Africa "worst-case" response-parameter curves using
#   MODEL C (the REVISED-methodology, endpoint-constrained model), then saves the
#   fitted curves for the combine step (03).
#
# HOW THIS DIFFERS FROM THE ORIGINAL-METHODOLOGY 01 (Model A)
#   Model A estimated BOTH the curve shape (t50, k) AND each parameter's two
#   magnitude endpoints, and exposed optional "tweak" priors. Model C instead:
#     * LOCKS each parameter's endpoints to early/late literature-window extrema
#       (computed by lock_endpoints() in helpers, the same rule used for DRC),
#       so the curve begins and ends at intended literature-supported values;
#     * estimates ONLY the curve SHAPE (t50, k), still partially pooled across
#       the six parameters; and
#     * has NO tweak priors at all - there is nothing to tweak, so there is also
#       no "with vs without tweaks" comparison in the revised methodology.
#   Everything else (the normalised logistic Q, the partial pooling, the
#   student-t likelihood, the 100-point prediction grid) is the same as Model A.
#
# Input  : data-processed/WestAfrica_QCurve_revisedMethod/WestAfrica_QCurve_PreppedData(revisedMethod).rds
#          (from 00; uses its $anchors element)
# Stan   : stan-models/modelC_endpointConstrained_estimateShape(revisedMethod).stan
# Output : data-processed/WestAfrica_QCurve_revisedMethod/WestAfrica_QCurve_Fit(revisedMethod).rds
#          data-processed/WestAfrica_QCurve_revisedMethod/WestAfrica_QCurve_Summaries(revisedMethod).csv
# ============================================================================

suppressPackageStartupMessages({
  library(cmdstanr)   # interface to CmdStan (compiles + samples the Stan model)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(tibble)
  library(ggplot2)   # display-only plots of the fitted curves at the end
})

source(here::here("analyses", "01_latent_response_parameter_estimation_revisedMethodology",
                  "helpers(revisedMethod).R"))

set.seed(123)  # reproducible R-side randomness (the sampler seed is set separately)

# ----------------------------------------------------------------------------
# 1. Read the cleaned West Africa anchors
# ----------------------------------------------------------------------------
# The $anchors element is the cleaned anchor table (one row per literature data
# point; see read_anchor_sheet() in 00 for the column meanings, e.g. parameter,
# relative_day, value_used, weight, direction, lower_bound, upper_bound).
wa_anchors <- readRDS(file.path(DIR_PROCESSED,
  "WestAfrica_QCurve_revisedMethod/WestAfrica_QCurve_PreppedData(revisedMethod).rds"))$anchors

# ----------------------------------------------------------------------------
# 2. Per-parameter metadata: direction and canonical ordering
# ----------------------------------------------------------------------------
# Model C does not estimate endpoints, so (unlike Model A) there is no admissible
# domain or endpoint prior to set up here - all we need per parameter is its
# direction and its canonical id. The magnitude endpoints come next, from the
# literature-window locking rule.
param_meta <- wa_anchors %>%
  group_by(parameter) %>%
  summarise(
    direction   = first(direction),                 # "up" or "down"
    lower_bound = first(lower_bound),               # workbook range, used only as
    upper_bound = first(upper_bound),               #   the endpoint-locking fallback
    .groups = "drop"
  ) %>%
  mutate(param_id = match(parameter, PARAM_LEVELS)) %>%
  arrange(param_id)

# Guard: every modelled parameter must be present, or the fit is mis-specified.
if (!all(PARAM_LEVELS %in% param_meta$parameter)) {
  stop("West Africa anchors are missing parameters: ",
       paste(setdiff(PARAM_LEVELS, param_meta$parameter), collapse = ", "))
}

J       <- nrow(param_meta)                          # number of parameters (6)
max_day <- max(wa_anchors$relative_day, na.rm = TRUE)  # the West Africa anchor horizon (~357)

# ----------------------------------------------------------------------------
# 3. LOCK the magnitude endpoints (this is the heart of the revised methodology)
# ----------------------------------------------------------------------------
# lock_endpoints() (in helpers) picks each parameter's start/end magnitude from
# early/late literature-window extrema, forces unsafe-funeral parameters with an
# observed zero to end at zero, and falls back to the workbook range when a
# window is empty. See its long comment for the full rule. The result is one row
# per parameter with:
#   theta_start          value at Q = 0 (worst response; direction already baked in)
#   theta_end            value at Q = 1 (best response)
#   endpoint_day_for_tau the day at which Q reaches 1 (scenario end, or earlier
#                        for a terminal-zero parameter)
#   start_source/end_source  audit labels saying where each endpoint came from
endpoint_table <- lock_endpoints(
  wa_anchors,
  scenario_duration_days = max_day,
  early_window_day       = 50,
  late_start_day         = 325
)

param_meta <- param_meta %>%
  left_join(endpoint_table, by = c("parameter", "direction")) %>%
  arrange(param_id)

message("West Africa locked endpoints:")
print(param_meta %>% select(parameter, direction, theta_start, theta_end,
                            endpoint_day_for_tau, start_source, end_source))

# ----------------------------------------------------------------------------
# 4. Assemble the per-observation fit table
# ----------------------------------------------------------------------------
# Stan needs one row per anchor OBSERVATION. We (a) put every anchor on the
# normalised clock tau = relative_day / endpoint_day_for_tau (so a parameter that
# reaches its endpoint early - a terminal-zero funeral - still has tau = 1 at its
# endpoint day), and (b) build a per-observation noise scale obs_sd.
#
# obs_sd combines two things, both in the parameter's native units:
#   * a baseline fraction (OBS_SD_FRAC_OF_SPAN) of the parameter's endpoint SPAN
#     |theta_end - theta_start| - i.e. wider-ranging parameters get more slack; and
#   * the anchor's trust weight, as the multiplier 1/weight (higher weight ->
#     tighter), exactly the obs_sd_mult convention used by Models A and B.
OBS_SD_FRAC_OF_SPAN <- 0.08      # baseline obs sd as a fraction of the endpoint span
MIN_OBS_SD          <- 1e-4      # floor so obs_sd can never be zero (Stan needs > 0)

# fit_df: one row per anchor observation fed to Stan, with columns -------------
#   param_id      which parameter the observation belongs to (1..6)
#   tau           normalised outbreak time of the observation, in [0, 1]
#   y_obs         the observed parameter value at that time (the anchor value)
#   obs_sd        per-observation observation scale (native units, > 0)
fit_df <- wa_anchors %>%
  left_join(select(param_meta, parameter, param_id, theta_start, theta_end,
                   endpoint_day_for_tau), by = "parameter") %>%
  mutate(
    tau         = clip01(relative_day / endpoint_day_for_tau),
    y_obs       = value_used,
    weight      = if_else(is.na(weight), 1, weight),
    obs_sd_mult = 1 / pmax(weight, 0.25),                       # looser when smaller weight
    span        = abs(theta_end - theta_start),
    obs_sd      = pmax(OBS_SD_FRAC_OF_SPAN * span * obs_sd_mult, MIN_OBS_SD)
  ) %>%
  arrange(param_id, tau)

# ----------------------------------------------------------------------------
# 5. Hyperpriors and the prediction grid
# ----------------------------------------------------------------------------
# These hyperpriors govern the PARTIAL POOLING of the curve shape across the six
# parameters (identical to Model A's shape hyperpriors):
#   mu_t50_raw_*  : group-mean midpoint (on the logit scale) and its prior spread
#   sigma_t50_*   : how much the parameters' midpoints may differ (pooling strength)
#   mu_log_k_*    : group-mean log-steepness and its prior spread
#   sigma_log_k_* : how much steepness may differ across parameters
hyper <- list(
  mu_t50_raw_prior_mean = 0.0,    mu_t50_raw_prior_sd = 1.25, sigma_t50_prior_sd   = 0.50,
  mu_log_k_prior_mean   = log(8), mu_log_k_prior_sd   = 0.80, sigma_log_k_prior_sd = 0.40
)

# Grid of normalised times at which to output the fitted curves (100 points,
# evenly spaced over the response window). 03 later places this on the 0..730 day
# axis (on each parameter's own endpoint_day_for_tau) and holds it flat to day 730.
tau_pred <- seq(0, 1, length.out = 100)

# ----------------------------------------------------------------------------
# 6. Assemble the Stan data list
# ----------------------------------------------------------------------------
# This list must contain EXACTLY the variables declared in Model C's `data`
# block: sizes, the observations, the FIXED endpoints, the shape hyperpriors and
# the prediction grid. Note theta_start/theta_end are DATA here - only the shape
# (t50, k) is estimated.
stan_data <- c(
  list(
    N = nrow(fit_df), J = J,
    param_id = fit_df$param_id, tau = fit_df$tau,
    y_obs = fit_df$y_obs, obs_sd = fit_df$obs_sd,
    theta_start = param_meta$theta_start, theta_end = param_meta$theta_end,
    N_pred = length(tau_pred), tau_pred = tau_pred
  ),
  hyper
)

# ----------------------------------------------------------------------------
# 7. Compile and sample
# ----------------------------------------------------------------------------
# Compile Model C (cached after the first run) and draw posterior samples.
# adapt_delta/max_treedepth are raised above their defaults because the
# hierarchical geometry can be awkward; 4 chains x 1500 post-warmup draws.
mod <- cmdstan_model(file.path(DIR_STAN,
  "modelC_endpointConstrained_estimateShape(revisedMethod).stan"))

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
# indices back out and attaches the parameter name, time, and day. Each
# parameter's day axis uses ITS OWN endpoint_day_for_tau (so a terminal-zero
# parameter's curve spans only up to its zero day, then is plateaued in step 9).
endpoint_day <- setNames(param_meta$endpoint_day_for_tau, param_meta$parameter)

tidy_matrix <- function(fit, varname) {
  fit$summary(variables = varname) %>%
    mutate(
      param_id = as.integer(str_match(variable, "\\[([0-9]+),([0-9]+)\\]")[, 2]),  # j
      grid_id  = as.integer(str_match(variable, "\\[([0-9]+),([0-9]+)\\]")[, 3])   # m
    ) %>%
    left_join(select(param_meta, param_id, parameter), by = "param_id") %>%
    mutate(
      tau          = tau_pred[grid_id],
      relative_day = tau * endpoint_day[parameter]
    ) %>%
    select(parameter, param_id, grid_id, tau, relative_day, mean, median, sd, q5, q95)
}

q_summ     <- tidy_matrix(fit, "Q_pred")      # per-parameter response-quality curve Q_j(tau)
curve_summ <- tidy_matrix(fit, "theta_pred")  # per-parameter value curve theta_j(tau)

# ----------------------------------------------------------------------------
# 9. Plateau any terminal-zero parameter out to the scenario end
# ----------------------------------------------------------------------------
# A terminal-zero parameter (e.g. an unsafe-funeral curve that reaches 0 before
# day 357) has curve points only up to its endpoint_day_for_tau. Append one
# plateau row at max_day holding the endpoint value, so the curve (and the
# combine step) carry the endpoint forward rather than stopping early. (03 also
# holds the last value flat, so this is belt-and-braces, but it makes the saved
# curve explicit and the display plot correct.)
plateau_rows <- curve_summ %>%
  group_by(parameter) %>%
  filter(max(relative_day, na.rm = TRUE) < max_day - 1e-8) %>%
  arrange(relative_day, .by_group = TRUE) %>%
  slice_tail(n = 1) %>%
  ungroup() %>%
  mutate(tau = 1, relative_day = max_day)

if (nrow(plateau_rows) > 0) {
  curve_summ <- bind_rows(curve_summ, plateau_rows) %>% arrange(parameter, relative_day)
}

# ----------------------------------------------------------------------------
# 10. Save
# ----------------------------------------------------------------------------
# WestAfrica_QCurve_Fit(revisedMethod).rds carries everything 03 needs to report
# the West Africa scenario. Its structure mirrors the original-methodology fit
# object (q_summ + curve_summ + param_meta + tau_pred + max_day), with the locked
# endpoint_table added and methodology tagged "revised".
wa_fit <- list(
  q_summ         = q_summ,         # Q_j(tau): per-parameter shape
  curve_summ     = curve_summ,     # theta_j(tau): per-parameter native-unit trajectory
  param_meta     = param_meta,     # per-parameter metadata + locked endpoints
  endpoint_table = endpoint_table, # the endpoint-locking audit table
  tau_pred       = tau_pred,       # the 100-point normalised-time grid
  max_day        = max_day,        # the West Africa anchor horizon (days)
  methodology    = "revised"
)
saveRDS(wa_fit, file.path(DIR_PROCESSED,
        "WestAfrica_QCurve_revisedMethod/WestAfrica_QCurve_Fit(revisedMethod).rds"))

write_csv(curve_summ, file.path(DIR_PROCESSED,
          "WestAfrica_QCurve_revisedMethod/WestAfrica_QCurve_Summaries(revisedMethod).csv"))

# ----------------------------------------------------------------------------
# 11. Plot the fitted curves with the anchor data on top  (display only)
# ----------------------------------------------------------------------------
# One facet per parameter: the posterior-MEAN fitted curve (blue) with its 90%
# credible-interval ribbon, the literature anchors the model was fit to (orange
# points), and a dashed grey marker at each locked endpoint value, so the fit and
# the endpoint locking can both be eyeballed. The plot is printed to the active
# graphics device and is deliberately NOT saved.
wa_curve_plot_df <- curve_summ %>%
  mutate(panel = factor(PANEL_LOOKUP[parameter], levels = unname(PANEL_LOOKUP)))
wa_anchor_plot_df <- fit_df %>%
  mutate(panel = factor(PANEL_LOOKUP[parameter], levels = unname(PANEL_LOOKUP)))
wa_endpoint_plot_df <- param_meta %>%
  transmute(parameter, theta_start, theta_end) %>%
  pivot_longer(c(theta_start, theta_end), names_to = "which", values_to = "value") %>%
  mutate(panel = factor(PANEL_LOOKUP[parameter], levels = unname(PANEL_LOOKUP)))

ggplot(wa_curve_plot_df, aes(relative_day, mean)) +
  geom_ribbon(aes(ymin = q5, ymax = q95), fill = "#1f77b4", alpha = 0.20) +   # 90% interval
  geom_line(colour = "#1f77b4", linewidth = 0.9) +                            # posterior mean
  geom_hline(data = wa_endpoint_plot_df, aes(yintercept = value),             # locked endpoints
             linetype = "dashed", colour = "grey55", linewidth = 0.3) +
  geom_point(data = wa_anchor_plot_df, aes(relative_day, y_obs),              # the anchor data
             inherit.aes = FALSE, colour = "#ff7f0e", size = 2) +
  facet_wrap(~ panel, scales = "free_y", ncol = 2) +
  labs(title = "West Africa Model C (revised): endpoint-constrained fitted curves",
       subtitle = "Blue = posterior mean + 90% interval; orange = literature anchors; dashed grey = locked endpoints",
       x = "Relative outbreak day", y = NULL) +
  theme_bw(base_size = 11) +
  theme(strip.text = element_text(face = "bold"))
