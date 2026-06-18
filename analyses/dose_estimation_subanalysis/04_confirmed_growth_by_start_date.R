# ============================================================================
# 04_confirmed_growth_by_start_date.R  (dose_estimation_subanalysis)
# ----------------------------------------------------------------------------
# PURPOSE
#   Sanity-check the growth-rate estimate for the observed national cumulative
#   CONFIRMED cases by scanning EVERY possible window START date against a FIXED
#   window END date (default 16 Jun 2026, the last confirmed observation). For
#   each start date we fit the exponential growth rate r (per day) and doubling
#   time (log(2)/r) over [start, END] THREE ways and compare them:
#
#     * daily (raw)      -- log-linear fit of log(daily incidence) ~ day. Targets
#                           the right quantity (incidence growth) but is NOISY
#                           (small counts + reporting artefacts: spikes, gaps).
#     * daily (smoothed) -- the same fit, but on a SMOOTHED daily-incidence curve
#                           (a quasi-Poisson GAM spline of time; loess fallback).
#                           Keeps the incidence target while removing day-to-day
#                           noise -- the recommended estimate.
#     * cumulative       -- log-linear fit of log(cumulative) ~ day. SMOOTH but
#                           BIASED: once incidence decelerates the cumulative
#                           log-slope lags and stays ABOVE the true (instantaneous)
#                           rate, and being monotone it can never go negative.
#
#   The companion context figure shows the raw incidence, the smoother, and the
#   instantaneous growth rate r(t) = d/dt log(smoothed incidence), which makes
#   the late-May spike-then-dip (the source of negative raw r) explicit.
#
# OUTPUT (in outputs/)
#   dose_confirmed_growth_by_start.csv          -- the full table (3 methods)
#   dose_confirmed_growth_by_start.png          -- r + doubling time vs start
#   dose_confirmed_incidence_context.png        -- incidence + smoother + r(t)
#
# RUN (from the repo root)
#   Rscript analyses/dose_estimation_subanalysis/04_confirmed_growth_by_start_date.R
# ============================================================================

suppressPackageStartupMessages({
  library(here)
  library(ggplot2)
})
source(here::here("analyses", "dose_estimation_subanalysis", "helpers.R"))  # DIR_OUT

# ---- config ----------------------------------------------------------------
CONFIRMED_CSV  <- here::here("data-processed",
                             "insp_sitrep__national_cumulative_confirmed_cases__daily.csv")
END_DATE       <- as.Date("2026-06-16")   # fixed window END (last confirmed obs)
MIN_WIN_POINTS <- 2L                       # min INCIDENCE points (>=2) for a fit

# ---- read + derive daily incidence -----------------------------------------
co <- utils::read.csv(CONFIRMED_CSV, stringsAsFactors = FALSE)
co$date <- as.Date(co$date, format = "%d/%m/%Y")
co <- co[!is.na(co$date) & co$date <= END_DATE, ]
co <- co[order(co$date), ]
cum_val <- co$national_cumulative_confirmed_cases

# Daily incidence = diff(cumulative)/diff(day), attributed to the later date.
inc <- data.frame(date      = co$date[-1],
                  incidence = diff(cum_val) / as.numeric(diff(co$date)))

# Day axis (numeric) measured from the first confirmed date, shared by all fits.
d0 <- min(co$date)
to_day <- function(d) as.numeric(d - d0)

# ---- smoother for the daily incidence --------------------------------------
# Prefer a quasi-Poisson GAM (log link) spline of time: it handles count noise,
# zeros and overdispersion and gives the underlying incidence trend without the
# cumulative bias. Falls back to a loess on log(incidence) if mgcv is absent.
fit_incidence_smoother <- function(day, value) {
  ok <- is.finite(day) & is.finite(value)
  df <- data.frame(day = day[ok], value = pmax(value[ok], 0))
  if (requireNamespace("mgcv", quietly = TRUE) && length(unique(df$day)) >= 6) {
    k <- max(4L, min(10L, length(unique(df$day)) %/% 3L))
    list(type = "gam",
         model = mgcv::gam(value ~ s(day, k = k), family = quasipoisson(), data = df))
  } else {
    message("Note: mgcv unavailable (or too few points) -- using a loess fallback.")
    list(type = "loess", model = stats::loess(log(value + 0.5) ~ day, data = df, span = 0.6))
  }
}
# Fitted incidence (response scale) at arbitrary days.
smooth_incidence_at <- function(fit, days) {
  if (fit$type == "gam")
    as.numeric(predict(fit$model, data.frame(day = days), type = "response"))
  else exp(as.numeric(predict(fit$model, data.frame(day = days))))
}
# Instantaneous growth rate r(t) = d/dt log(fitted incidence) by central
# difference (for the log-link GAM this is the derivative of the fitted smooth).
smooth_r_at <- function(fit, days, h = 0.5) {
  (log(smooth_incidence_at(fit, days + h)) -
   log(smooth_incidence_at(fit, days - h))) / (2 * h)
}

inc_fit       <- fit_incidence_smoother(to_day(inc$date), inc$incidence)
inc$smoothed  <- smooth_incidence_at(inc_fit, to_day(inc$date))

# ---- log-linear growth-rate helper -----------------------------------------
# r = slope of log(value) ~ day over the supplied points; NA if < 2 positive pts.
loglin_r <- function(day, value) {
  ok  <- is.finite(day) & is.finite(value) & value > 0
  day <- day[ok]; value <- value[ok]
  if (length(unique(day)) < 2L) return(c(r = NA_real_, n = length(value)))
  c(r = unname(coef(stats::lm(log(value) ~ day))[2]), n = length(value))
}
dbl <- function(r) if (is.finite(r) && r > 0) log(2) / r else NA_real_

# ---- scan every start date --------------------------------------------------
# Candidate starts: every confirmed observation date that still leaves at least
# MIN_WIN_POINTS incidence points up to END_DATE.
starts <- co$date[co$date < END_DATE]
rows <- list()
for (s in starts) {
  s <- as.Date(s, origin = "1970-01-01")

  inc_w <- inc[inc$date >= s & inc$date <= END_DATE, ]
  cum_w <- co[co$date  >= s & co$date  <= END_DATE, ]
  if (nrow(inc_w) < MIN_WIN_POINTS) next

  fi <- loglin_r(to_day(inc_w$date), inc_w$incidence)   # raw daily incidence
  fs <- loglin_r(to_day(inc_w$date), inc_w$smoothed)    # smoothed daily incidence
  fc <- loglin_r(to_day(cum_w$date), cum_w$national_cumulative_confirmed_cases)  # cumulative

  rows[[length(rows) + 1L]] <- data.frame(
    start_date = s, end_date = END_DATE,
    method        = c("daily (raw)", "daily (smoothed)", "cumulative"),
    growth_rate   = c(fi[["r"]], fs[["r"]], fc[["r"]]),
    doubling_time = c(dbl(fi[["r"]]), dbl(fs[["r"]]), dbl(fc[["r"]])),
    n_points      = c(fi[["n"]], fs[["n"]], fc[["n"]]),
    row.names     = NULL)
}
tbl <- do.call(rbind, rows)
method_levels <- c("daily (raw)", "daily (smoothed)", "cumulative")
tbl$method <- factor(tbl$method, levels = method_levels)

write.csv(tbl, file.path(DIR_OUT, "dose_confirmed_growth_by_start.csv"), row.names = FALSE)

cat(sprintf("\nConfirmed-case growth rate by window start (END = %s):\n",
            format(END_DATE, "%d %b %Y")))
disp <- tbl
disp$growth_rate   <- round(disp$growth_rate, 4)
disp$doubling_time <- round(disp$doubling_time, 1)
print(disp[order(disp$method, disp$start_date), ], row.names = FALSE)

for (m in method_levels) {
  nn <- sum(tbl$method == m & is.finite(tbl$growth_rate) & tbl$growth_rate < 0)
  cat(sprintf("  %-18s negative-growth windows: %d of %d\n",
              m, nn, sum(tbl$method == m)))
}

# ---- figure 1: growth rate + doubling time vs start date --------------------
cols <- c("daily (raw)"      = "#969696",   # noisy baseline
          "daily (smoothed)" = "#1a9850",   # recommended
          "cumulative"       = "#762a83")   # biased
long <- tidyr::pivot_longer(tbl, c(growth_rate, doubling_time),
                            names_to = "metric", values_to = "value")
long$metric <- factor(long$metric, levels = c("growth_rate", "doubling_time"),
                       labels = c("Growth rate r (per day)", "Doubling time (days)"))
# Hide very long doubling times (near-zero growth) so the panel stays readable.
long$value[long$metric == "Doubling time (days)" &
           is.finite(long$value) & long$value > 120] <- NA

p1 <- ggplot(long, aes(start_date, value, colour = method)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_line(linewidth = 0.6) +
  geom_point(aes(size = n_points), alpha = 0.85) +
  facet_wrap(~ metric, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = cols, name = "Fit on") +
  scale_size_continuous(range = c(0.8, 3.2), name = "n obs in window") +
  scale_x_date(date_breaks = "1 week", date_labels = "%d %b") +
  labs(title = "Confirmed-case growth rate vs window start date",
       subtitle = sprintf("Window = [start, %s]. Raw daily = noisy; cumulative = biased high (& never < 0); smoothed = recommended. Doubling times > 120 d hidden.",
                          format(END_DATE, "%d %b %Y")),
       x = "Window start date", y = NULL) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(DIR_OUT, "dose_confirmed_growth_by_start.png"), p1,
       width = 9, height = 7, dpi = 150)
print(p1)

# ---- figure 2: incidence + smoother + instantaneous r(t) --------------------
# Fine grid for the smoother / its derivative.
grid_day <- seq(min(to_day(inc$date)), max(to_day(inc$date)), by = 0.5)
grid_dat <- d0 + grid_day
inc_panel <- data.frame(date = inc$date,  value = inc$incidence,
                        panel = "Daily incidence")
fit_panel <- data.frame(date = grid_dat,  value = smooth_incidence_at(inc_fit, grid_day),
                        panel = "Daily incidence")
r_panel   <- data.frame(date = grid_dat,  value = smooth_r_at(inc_fit, grid_day),
                        panel = "Instantaneous growth r(t) (per day)")
panel_lvls <- c("Daily incidence", "Instantaneous growth r(t) (per day)")
inc_panel$panel <- factor(inc_panel$panel, panel_lvls)
fit_panel$panel <- factor(fit_panel$panel, panel_lvls)
r_panel$panel   <- factor(r_panel$panel,   panel_lvls)

p2 <- ggplot(mapping = aes(date, value)) +
  geom_col(data = inc_panel, fill = "#969696", alpha = 0.55) +
  geom_line(data = fit_panel, colour = "#1a9850", linewidth = 1.0) +
  geom_line(data = r_panel,   colour = "#08519c", linewidth = 1.0) +
  geom_hline(data = data.frame(panel = factor("Instantaneous growth r(t) (per day)", panel_lvls),
                               y0 = 0),
             aes(yintercept = y0), inherit.aes = FALSE, colour = "grey60", linewidth = 0.3) +
  facet_wrap(~ panel, ncol = 1, scales = "free_y") +
  scale_x_date(date_breaks = "1 week", date_labels = "%d %b") +
  labs(title = "Confirmed daily incidence, smoother, and instantaneous growth rate",
       subtitle = sprintf("Grey bars = raw daily incidence; green = %s smoother; blue = d/dt log(smoothed). r(t) < 0 where incidence is falling.",
                          inc_fit$type),
       x = "Date", y = NULL) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(DIR_OUT, "dose_confirmed_incidence_context.png"), p2,
       width = 9, height = 6.5, dpi = 150)
print(p2)

message("\n04_confirmed_growth_by_start_date.R complete. Outputs in outputs/.")
