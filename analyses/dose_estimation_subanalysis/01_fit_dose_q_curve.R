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
#          outputs/dose_q_curve_extrapolation_scenarios.{rds,csv,png}
#                                          <- 4 forward-extrapolation scenarios + plot
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
#   add a single anchoring zero at the start date (relative day 0), so the curve
#   starts from 0% there. Default = 18 May 2026 (no added zero; relative days
#   {0, 6, 13, 20, 27}).
START_DATE   <- as.Date("2026-05-14")

# How far PAST the last observation to project the curve forward (days). The
# curve is reported on a daily grid from day 0 out to (last obs day + this).
FORWARD_DAYS <- 365L

# Endpoint priors -- "very informative" pins on the min (0) and max (1) levels.
# Tighten/loosen these SDs to harden/soften how strictly the curve is held to
# the 0% floor and 100% ceiling.
LOWER_PRIOR_MEAN <- 0.0; LOWER_PRIOR_SD <- 0.01   # min level: 0%
UPPER_PRIOR_MEAN <- 1.0; UPPER_PRIOR_SD <- 0.3   # max level: 100%

# Observation-noise prior scale (half-normal); proportions live in [0, 1].
SIGMA_PRIOR_SD <- 0.15

# Sampler settings.
N_CHAINS <- 4L; ITER_WARMUP <- 1500L; ITER_SAMPLING <- 1500L

# Forward-extrapolation scenario knobs (section 5). Each scenario continues the
# curve AFTER the last data point, starting from the fitted value there (q_end).
LINEAR_TARGET_Q      <- 0.90   # scenario 1: straight line up to this Q ...
LINEAR_TARGET_DAY    <- 100L   #             ... reached by this day, flat after
CONFLICT_START_DAY   <- 100L   # scenario 4: conflict begins here
CONFLICT_DROP_DAYS   <- 25L    #   drop from q_end to CONFLICT_LOW_Q over this many days
CONFLICT_LOW_DAYS    <- 50L    #   then hold at CONFLICT_LOW_Q this many days
CONFLICT_REVERT_DAYS <- 25L    #   then revert to q_end over this many days
CONFLICT_LOW_Q       <- 0.20   #   the depressed Q during conflict

# ----------------------------------------------------------------------------
# 2. Fit + extrapolate the logistic Q curve  (fit_dose_q_curve() in helpers.R)
# ----------------------------------------------------------------------------
# The fitting itself -- padding, data-derived shape priors, the Stan model and
# the extraction of the daily Q curve out to FORWARD_DAYS past the last
# observation -- lives in helpers.R, so 03_compare_start_dates.R fits in exactly
# the same way (no drift).
res <- fit_dose_q_curve(
  start_date     = START_DATE,
  forward_days   = FORWARD_DAYS,
  lower_prior    = c(mean = LOWER_PRIOR_MEAN, sd = LOWER_PRIOR_SD),
  upper_prior    = c(mean = UPPER_PRIOR_MEAN, sd = UPPER_PRIOR_SD),
  sigma_prior_sd = SIGMA_PRIOR_SD,
  chains = N_CHAINS, iter_warmup = ITER_WARMUP, iter_sampling = ITER_SAMPLING
)

obs           <- res$observations
q_curve       <- res$q_curve
param_summary <- res$param_summary
last_obs_day  <- res$last_obs_day
horizon_day   <- res$horizon_day

message("Dose observations (relative-day axis, day 0 = ", START_DATE, "):")
print(obs, row.names = FALSE)
cat("\nSampler diagnostics:\n");        print(res$fit$diagnostic_summary())
cat("\nFitted logistic parameters:\n"); print(param_summary)

# ----------------------------------------------------------------------------
# 3. Save  (the Q curve is the headline output)
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
    priors = res$priors_used
  )
)
saveRDS(q_curve_fit, file.path(DIR_OUT, "dose_q_curve_fit.rds"))

message("\n01_fit_dose_q_curve.R complete. Wrote dose_q_curve.rds (",
        nrow(q_curve), " daily rows, ", last_obs_day, " observed + ",
        FORWARD_DAYS, " projected days) to outputs/.")

# ----------------------------------------------------------------------------
# 4. Diagnostic plot: fit + extrapolation vs data
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

# ----------------------------------------------------------------------------
# 5. Forward-extrapolation scenarios (applied AFTER the last data point)
# ----------------------------------------------------------------------------
# All four scenarios share the fitted curve up to last_obs_day and start their
# continuation from the fitted value there (q_end). They differ only afterwards:
#   1. linear_to_90  -- straight line from (last_obs_day, q_end) to
#                        (LINEAR_TARGET_DAY, LINEAR_TARGET_Q); flat thereafter.
#   2. logistic      -- the fitted logistic projected forward (current behaviour).
#   3. flat          -- held flat at q_end (no further improvement).
#   4. conflict      -- flat at q_end until CONFLICT_START_DAY, then a conflict
#                        drops Q to CONFLICT_LOW_Q over CONFLICT_DROP_DAYS, holds
#                        it for CONFLICT_LOW_DAYS, reverts to q_end over
#                        CONFLICT_REVERT_DAYS; flat at q_end afterwards.
# dose_q_curve.rds (read by 02) stays the logistic projection; these scenarios
# are saved separately so they can be selected downstream if desired.

# Fitted value at the last data point -- the common start of every extrapolation.
q_end <- q_curve$q_mean[match(last_obs_day, q_curve$relative_day)]

# Common daily grid, long enough to show the whole conflict episode.
conflict_end <- CONFLICT_START_DAY + CONFLICT_DROP_DAYS + CONFLICT_LOW_DAYS +
  CONFLICT_REVERT_DAYS
scen_horizon <- max(horizon_day, conflict_end + 60L)
scen_days    <- 0:scen_horizon
post         <- scen_days > last_obs_day                 # the extrapolation region

# The fitted logistic on the grid (flat-held beyond its own horizon, where it has
# saturated). This IS scenario 2 and the shared [0, last_obs_day] segment.
fit_on <- clip01(approx(q_curve$relative_day, q_curve$q_mean,
                        xout = scen_days, rule = 2)$y)

# 1. Linear ramp to LINEAR_TARGET_Q by LINEAR_TARGET_DAY, flat thereafter.
s_linear <- fit_on
if (LINEAR_TARGET_DAY > last_obs_day) {
  lin_fun <- approxfun(c(last_obs_day, LINEAR_TARGET_DAY),
                       c(q_end, LINEAR_TARGET_Q), rule = 2)
  s_linear[post] <- lin_fun(scen_days[post])
} else {
  s_linear[post] <- LINEAR_TARGET_Q
}

# 2. Logistic forward projection (the fit itself).
s_logistic <- fit_on

# 3. Flat at q_end.
s_flat <- fit_on
s_flat[post] <- q_end

# 4. Flat, then a conflict episode, then revert to q_end. A single piecewise-
#    linear function over the four knots does it: rule = 2 holds q_end before the
#    drop and after the reversion.
conflict_fun <- approxfun(
  x = c(CONFLICT_START_DAY,
        CONFLICT_START_DAY + CONFLICT_DROP_DAYS,
        CONFLICT_START_DAY + CONFLICT_DROP_DAYS + CONFLICT_LOW_DAYS,
        conflict_end),
  y = c(q_end, CONFLICT_LOW_Q, CONFLICT_LOW_Q, q_end),
  rule = 2
)
s_conflict <- fit_on
s_conflict[post] <- conflict_fun(scen_days[post])

scenarios <- rbind(
  data.frame(relative_day = scen_days, q = clip01(s_linear),   scenario = "linear_to_90"),
  data.frame(relative_day = scen_days, q = clip01(s_linear),   scenario = "linear_to_95"),
  data.frame(relative_day = scen_days, q = clip01(s_logistic), scenario = "logistic"),
  data.frame(relative_day = scen_days, q = clip01(s_flat),     scenario = "flat"),
  data.frame(relative_day = scen_days, q = clip01(s_conflict), scenario = "conflict")
)
scenarios$date <- START_DATE + scenarios$relative_day

scen_levels <- c("linear_to_90", "logistic", "flat", "conflict")
scen_labels <- c(linear_to_90 = sprintf("1. Linear to %d%% by day %d",
                                        round(100 * LINEAR_TARGET_Q), LINEAR_TARGET_DAY),
                 logistic     = "2. Logistic projection (current)",
                 flat         = "3. Flat at last value",
                 conflict     = sprintf("4. Conflict at day %d", CONFLICT_START_DAY))
scenarios$scenario <- factor(scenarios$scenario, levels = scen_levels)

saveRDS(scenarios, file.path(DIR_OUT, "dose_q_curve_extrapolation_scenarios.rds"))
write.csv(scenarios, file.path(DIR_OUT, "dose_q_curve_extrapolation_scenarios.csv"),
          row.names = FALSE)

# Plot: the fit (90% band + data) plus the four extrapolations on one graph, in
# calendar time with the dose-data window marked (matching 02).
data_first_date <- min(DOSE_OBS$date)   # 18 May 2026 (first dose observation)
data_last_date  <- max(DOSE_OBS$date)   # 14 Jun 2026 (last  dose observation)

p_scen <- ggplot(scenarios, aes(date, q, colour = scenario)) +
  geom_ribbon(data = q_curve[q_curve$relative_day <= last_obs_day, ],
              aes(x = date, ymin = q_lower, ymax = q_upper),
              inherit.aes = FALSE, fill = "grey70", alpha = 0.35) +
  geom_vline(xintercept = data_first_date, linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = data_last_date,  linetype = "dashed", colour = "grey40") +
  geom_line(linewidth = 0.9) +
  geom_point(data = obs, aes(date, proportion),
             inherit.aes = FALSE, colour = "black", size = 2.2) +
  scale_colour_brewer(palette = "Dark2", labels = scen_labels, name = "Extrapolation") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(title = "Dose Q curve: fit + forward-extrapolation scenarios",
       subtitle = sprintf("Grey = 90%% CI, points = data; dashed verticals = dose-data window (%s, %s); scenarios diverge after the last data point",
                          format(data_first_date, "%d %b"), format(data_last_date, "%d %b")),
       x = "Date", y = "Dose coverage (Q)") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(DIR_OUT, "dose_q_curve_extrapolation_scenarios.png"), p_scen,
       width = 9.5, height = 5.5, dpi = 150)
print(p_scen)   # also display
