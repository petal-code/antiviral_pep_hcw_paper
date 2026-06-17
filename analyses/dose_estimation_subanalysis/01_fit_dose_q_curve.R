# ============================================================================
# 01_fit_dose_q_curve.R
# ----------------------------------------------------------------------------
# WHAT THIS SCRIPT DOES
#   Fits a single logistic "Q curve" to five (date, dose-coverage %) pairs, of
#   the SAME family as the previous Q-curve analyses (Model A: a logistic scaled
#   between a lower and an upper endpoint, robust Student-t likelihood), and
#   projects it forward up to a year past the data. The fitted/extrapolated Q
#   curve is saved to outputs/ as the headline artefact.
#
# THE DATA
#   {18 May 2026, 0%}, {24 May 2026, 19.3%}, {31 May 2026, 30.2%},
#   {7 Jun 2026, 64.4%}, {14 Jun 2026, 63.1%}
#   Percentages -> proportions; calendar dates -> RELATIVE days from a start
#   date (day 0). The front of the series can be padded with zeros by choosing a
#   start date earlier than 18 May (see START_DATE below).
#
# THE MODEL  (stan-models/logistic_qcurve_single.stan)
#   q(day) = L + (U - L) * inv_logit( k * (day - t50) )
#   with VERY informative priors pinning L ~ 0 (min 0%) and U ~ 1 (max 100%);
#   the data inform only the shape (midpoint t50, growth rate k). The logistic is
#   in raw day units (NOT renormalised to the window), so projecting it forward
#   simply means evaluating it on later days at the estimated growth rate.
#
# Stan   : stan-models/logistic_qcurve_single.stan
# Output : outputs/dose_q_curve.rds       <- the Q curve itself (daily data.frame)
#          outputs/dose_q_curve_fit.rds   <- full fit (params, draws summary, data)
#          outputs/dose_q_curve.csv       <- the Q curve as plain CSV
#          outputs/dose_q_curve_fit.png   <- diagnostic fit + extrapolation plot
# ============================================================================

suppressPackageStartupMessages({
  library(cmdstanr)
  library(dplyr)
  library(ggplot2)
})

source(here::here("analyses", "dose_estimation_subanalysis", "helpers.R"))

set.seed(123)

# ----------------------------------------------------------------------------
# 1. Configuration  <-- the knobs to change
# ----------------------------------------------------------------------------
# START_DATE: day 0 of the relative-day axis. Set EARLIER than 18 May 2026 to
#   pad the front with zeros (one zero per day from the start date up to 18 May).
#   Default = 18 May 2026 (no padding; relative days {0, 6, 13, 20, 27}).
START_DATE   <- as.Date("2026-05-18")

# How far PAST the last observation to project the curve forward (days). The
# curve is reported on a daily grid from day 0 out to (last obs day + this).
FORWARD_DAYS <- 365L

# Endpoint priors -- "very informative" pins on the min (0) and max (1) levels.
# Tighten/loosen these SDs to harden/soften how strictly the curve is held to
# the 0% floor and 100% ceiling.
LOWER_PRIOR_MEAN <- 0.0; LOWER_PRIOR_SD <- 0.01   # min level: 0%
UPPER_PRIOR_MEAN <- 1.0; UPPER_PRIOR_SD <- 0.02   # max level: 100%

# Observation-noise prior scale (half-normal); proportions live in [0, 1].
SIGMA_PRIOR_SD <- 0.15

# Sampler settings.
N_CHAINS <- 4L; ITER_WARMUP <- 1500L; ITER_SAMPLING <- 1500L

# ----------------------------------------------------------------------------
# 2. Assemble the relative-day observations
# ----------------------------------------------------------------------------
obs <- build_dose_obs(DOSE_OBS, start_date = START_DATE)

message("Dose observations (relative-day axis, day 0 = ", START_DATE, "):")
print(obs, row.names = FALSE)

last_obs_day <- max(obs$relative_day)
horizon_day  <- last_obs_day + FORWARD_DAYS         # last day of the reported curve
day_pred     <- 0:horizon_day                       # daily prediction grid

# Weakly-informative shape priors, derived from the observation window so the
# data dominate: midpoint centred on the middle of the window with an SD as wide
# as the window; growth rate around ~0.2 / day (rises over a few weeks).
data_span    <- diff(range(obs$relative_day))
T50_PRIOR_MEAN <- mean(range(obs$relative_day))
T50_PRIOR_SD   <- max(data_span, 10)
LOGK_PRIOR_MEAN <- log(0.2)
LOGK_PRIOR_SD   <- 0.7

# ----------------------------------------------------------------------------
# 3. Fit the logistic Q curve in Stan
# ----------------------------------------------------------------------------
stan_data <- list(
  N = nrow(obs), day = as.numeric(obs$relative_day), y = as.numeric(obs$proportion),
  lower_prior_mean = LOWER_PRIOR_MEAN, lower_prior_sd = LOWER_PRIOR_SD,
  upper_prior_mean = UPPER_PRIOR_MEAN, upper_prior_sd = UPPER_PRIOR_SD,
  t50_prior_mean = T50_PRIOR_MEAN, t50_prior_sd = T50_PRIOR_SD,
  logk_prior_mean = LOGK_PRIOR_MEAN, logk_prior_sd = LOGK_PRIOR_SD,
  sigma_prior_sd = SIGMA_PRIOR_SD,
  N_pred = length(day_pred), day_pred = as.numeric(day_pred)
)

mod <- cmdstan_model(file.path(DIR_STAN, "logistic_qcurve_single.stan"))

fit <- mod$sample(
  data = stan_data, seed = 123,
  chains = N_CHAINS, parallel_chains = N_CHAINS,
  iter_warmup = ITER_WARMUP, iter_sampling = ITER_SAMPLING,
  adapt_delta = 0.95, max_treedepth = 12, refresh = 200
)

cat("\nSampler diagnostics:\n"); print(fit$diagnostic_summary())

# ----------------------------------------------------------------------------
# 4. Extract the fitted / extrapolated Q curve
# ----------------------------------------------------------------------------
# q_pred[m] is the posterior Q at day_pred[m]. Attach the calendar date and a
# flag separating the observed-window segment from the forward extrapolation.
q_curve <- fit$summary(variables = "q_pred") %>%
  mutate(grid_id = as.integer(stringr::str_match(variable, "\\[([0-9]+)\\]")[, 2])) %>%
  arrange(grid_id) %>%
  transmute(
    relative_day = day_pred[grid_id],
    date         = START_DATE + day_pred[grid_id],
    q_mean       = mean,
    q_median     = median,
    q_lower      = q5,
    q_upper      = q95,
    segment      = ifelse(day_pred[grid_id] <= last_obs_day,
                          "observed_window", "extrapolated")
  )

# Posterior summaries of the logistic parameters.
param_summary <- fit$summary(variables = c("L", "U", "t50", "k", "sigma"))
cat("\nFitted logistic parameters:\n"); print(param_summary)

# ----------------------------------------------------------------------------
# 5. Save  (the Q curve is the headline output)
# ----------------------------------------------------------------------------
# dose_q_curve.rds is LITERALLY the Q curve (one row per day, mean + 90% band,
# observed-vs-extrapolated flag). dose_q_curve_fit.rds carries the full fit for
# provenance / reuse.
saveRDS(q_curve, file.path(DIR_OUT, "dose_q_curve.rds"))
write.csv(q_curve, file.path(DIR_OUT, "dose_q_curve.csv"), row.names = FALSE)

q_curve_fit <- list(
  q_curve       = q_curve,
  param_summary = param_summary,
  observations  = obs,
  config = list(
    start_date = START_DATE, forward_days = FORWARD_DAYS,
    last_obs_day = last_obs_day, horizon_day = horizon_day,
    priors = list(
      lower = c(mean = LOWER_PRIOR_MEAN, sd = LOWER_PRIOR_SD),
      upper = c(mean = UPPER_PRIOR_MEAN, sd = UPPER_PRIOR_SD),
      t50   = c(mean = T50_PRIOR_MEAN,   sd = T50_PRIOR_SD),
      logk  = c(mean = LOGK_PRIOR_MEAN,  sd = LOGK_PRIOR_SD),
      sigma = c(sd = SIGMA_PRIOR_SD)
    )
  )
)
saveRDS(q_curve_fit, file.path(DIR_OUT, "dose_q_curve_fit.rds"))

message("\n01_fit_dose_q_curve.R complete. Wrote dose_q_curve.rds (",
        nrow(q_curve), " daily rows, ", last_obs_day, " observed + ",
        FORWARD_DAYS, " projected days) to outputs/.")

# ----------------------------------------------------------------------------
# 6. Diagnostic plot: fit + extrapolation vs data
# ----------------------------------------------------------------------------
p <- ggplot(q_curve, aes(relative_day, q_mean)) +
  geom_ribbon(aes(ymin = q_lower, ymax = q_upper), fill = "#1f77b4", alpha = 0.18) +
  geom_line(aes(linetype = segment), colour = "#1f77b4", linewidth = 0.9) +
  geom_vline(xintercept = last_obs_day, linetype = "dotted", colour = "grey50") +
  geom_point(data = obs, aes(relative_day, proportion,
                             shape = ifelse(padded, "padded zero", "observation")),
             inherit.aes = FALSE, colour = "#ff7f0e", size = 2.4) +
  scale_linetype_manual(values = c(observed_window = "solid", extrapolated = "22"),
                        name = NULL) +
  scale_shape_manual(values = c("observation" = 16, "padded zero" = 1), name = NULL) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(title = "Dose-coverage Q curve: logistic fit + forward projection",
       subtitle = sprintf("Day 0 = %s; dotted line = last observation (day %d); curve pinned to 0%% / 100%% endpoints",
                          as.character(START_DATE), last_obs_day),
       x = "Relative day", y = "Dose coverage (Q)") +
  theme_bw(base_size = 11)

ggsave(file.path(DIR_OUT, "dose_q_curve_fit.png"), p, width = 9, height = 5, dpi = 150)
print(p)   # also display
