# ============================================================================
# 03_compare_start_dates.R
# ----------------------------------------------------------------------------
# WHAT THIS SCRIPT DOES
#   Refits the dose-coverage logistic Q curve for several candidate START DATES
#   (i.e. several amounts of front zero-padding before 18 May), then overlays all
#   of the fitted/extrapolated curves on ONE plot so the effect of the start-date
#   choice is visible at a glance.
#
#   Each curve is fit with EXACTLY the same model/priors as 01_fit_dose_q_curve.R
#   (both call fit_dose_q_curve() in helpers.R). The only thing that changes
#   between curves is the start date, which sets how many leading zeros are added
#   between the start date and 18 May.
#
#   The curves are drawn against the CALENDAR DATE (not relative day) so they are
#   directly comparable: the real observations sit at fixed calendar dates, and
#   each curve simply begins at its own start date and shares the same forward
#   horizon.
#
# Stan   : stan-models/logistic_qcurve_single.stan  (compiled once, reused)
# Output : outputs/dose_q_curve_startdate_comparison.rds / .csv  (all curves)
#          outputs/dose_q_curve_startdate_params.csv             (fitted params)
#          outputs/dose_q_curve_startdate_comparison.png         (overlay plot)
#          outputs/dose_q_curve_startdate_comparison_zoom.png    (zoomed overlay)
# ============================================================================

suppressPackageStartupMessages({
  library(cmdstanr)
  library(ggplot2)
})

source(here::here("analyses", "dose_estimation_subanalysis", "helpers.R"))

set.seed(123)

# ----------------------------------------------------------------------------
# 1. Configuration
# ----------------------------------------------------------------------------
# Candidate start dates: 1, 7, 13 and 18 May 2026. An earlier start date adds a
# single anchoring zero at that date (relative day 0), shifting where the curve
# leaves the 0% floor.
START_DATES  <- as.Date(c("2026-05-01", "2026-05-07", "2026-05-13", "2026-05-18"))
FORWARD_DAYS <- 365L   # forward projection past the last observation (matches 01)

# Priors -- the SAME knobs as 01_fit_dose_q_curve.R, applied to every start-date
# fit. Edit here to tune the comparison (keep in sync with 01 if you want the
# 18 May curve to match the headline fit exactly).
LOWER_PRIOR_MEAN <- 0.0; LOWER_PRIOR_SD <- 0.01   # min level: 0%
UPPER_PRIOR_MEAN <- 1.0; UPPER_PRIOR_SD <- 0.02   # max level: 100%
SIGMA_PRIOR_SD   <- 0.15

# ----------------------------------------------------------------------------
# 2. Fit one Q curve per start date (compile the model once, reuse it)
# ----------------------------------------------------------------------------
mod <- cmdstan_model(file.path(DIR_STAN, "logistic_qcurve_single.stan"))

fits <- lapply(START_DATES, function(sd) {
  message("Fitting start date ", sd, " ...")
  fit_dose_q_curve(
    start_date     = sd,
    forward_days   = FORWARD_DAYS,
    lower_prior    = c(mean = LOWER_PRIOR_MEAN, sd = LOWER_PRIOR_SD),
    upper_prior    = c(mean = UPPER_PRIOR_MEAN, sd = UPPER_PRIOR_SD),
    sigma_prior_sd = SIGMA_PRIOR_SD,
    mod = mod, refresh = 0
  )
})
names(fits) <- as.character(START_DATES)

# ----------------------------------------------------------------------------
# 3. Combine the curves and the fitted parameters
# ----------------------------------------------------------------------------
curves <- do.call(rbind, lapply(fits, function(f)
  cbind(start_date = as.character(f$start_date), f$q_curve)))
curves$start_date <- factor(curves$start_date, levels = as.character(START_DATES))

params <- do.call(rbind, lapply(fits, function(f) {
  d <- as.data.frame(f$param_summary)[, c("variable", "mean", "q5", "q95")]
  cbind(start_date = as.character(f$start_date), d)
}))
cat("\nFitted logistic parameters by start date:\n"); print(params, row.names = FALSE)

# ----------------------------------------------------------------------------
# 4. Save  (compact per-start-date details; not the heavy CmdStan fit objects,
#           which hold references to temp CSVs and don't reload cleanly)
# ----------------------------------------------------------------------------
details <- lapply(fits, function(f) list(
  observations  = f$observations,
  param_summary = as.data.frame(f$param_summary),
  priors_used   = f$priors_used,
  last_obs_day  = f$last_obs_day,
  horizon_day   = f$horizon_day,
  start_date    = f$start_date
))
saveRDS(list(curves = curves, params = params, details = details,
             start_dates = START_DATES, forward_days = FORWARD_DAYS),
        file.path(DIR_OUT, "dose_q_curve_startdate_comparison.rds"))
write.csv(curves[, c("start_date", "relative_day", "date",
                     "q_mean", "q_lower", "q_upper", "segment")],
          file.path(DIR_OUT, "dose_q_curve_startdate_comparison.csv"), row.names = FALSE)
write.csv(params, file.path(DIR_OUT, "dose_q_curve_startdate_params.csv"), row.names = FALSE)

# ----------------------------------------------------------------------------
# 5. Overlay plot (all start dates on one calendar axis)
# ----------------------------------------------------------------------------
# The observations are fixed in calendar time and shared by every fit.
obs_cal <- data.frame(date = DOSE_OBS$date, proportion = DOSE_OBS$percentage / 100)
last_obs_date <- max(DOSE_OBS$date)   # 14 Jun 2026

base_plot <- function(dat) {
  ggplot(dat, aes(date, q_mean, colour = start_date, fill = start_date,
                  group = interaction(start_date, segment))) +
    geom_ribbon(aes(ymin = q_lower, ymax = q_upper), colour = NA, alpha = 0.10) +
    geom_line(aes(linetype = segment), linewidth = 0.9) +
    geom_vline(xintercept = last_obs_date, linetype = "dotted", colour = "grey50") +
    geom_point(data = obs_cal, aes(date, proportion), inherit.aes = FALSE,
               colour = "black", size = 2.2) +
    scale_linetype_manual(values = c(observed_window = "solid", extrapolated = "22"),
                          guide = "none") +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    labs(colour = "Start date", fill = "Start date",
         x = NULL, y = "Dose coverage (Q)") +
    theme_bw(base_size = 11)
}

# (i) Full horizon.
p_full <- base_plot(curves) +
  labs(title = "Dose-coverage Q curve: sensitivity to start date (single anchoring zero)",
       subtitle = "Solid = fit over the observed window; dashed = forward projection; black points = observations; dotted = last observation")
ggsave(file.path(DIR_OUT, "dose_q_curve_startdate_comparison.png"),
       p_full, width = 10, height = 5.5, dpi = 150)
print(p_full)

# (ii) Zoom on the informative region (start of the earliest curve -> ~6 weeks
# past the last observation), where the curves actually differ; far out they all
# saturate near 100%.
zoom_to <- last_obs_date + 45L
p_zoom <- base_plot(curves[curves$date <= zoom_to, , drop = FALSE]) +
  labs(title = "Dose-coverage Q curve by start date (zoom on the data window)",
       subtitle = "Same fits, zoomed: an earlier start date adds one anchoring zero, shifting where the curve leaves 0%")
ggsave(file.path(DIR_OUT, "dose_q_curve_startdate_comparison_zoom.png"),
       p_zoom, width = 10, height = 5.5, dpi = 150)
print(p_zoom)

message("\n03_compare_start_dates.R complete. Wrote ", length(START_DATES),
        " start-date curves + overlay plots to outputs/.")
