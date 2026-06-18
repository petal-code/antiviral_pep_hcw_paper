# ============================================================================
# 04_confirmed_growth_by_start_date.R  (dose_estimation_subanalysis)
# ----------------------------------------------------------------------------
# PURPOSE
#   Sanity-check the growth-rate estimate for the observed national cumulative
#   CONFIRMED cases by scanning EVERY possible window START date against a FIXED
#   window END date (default 16 Jun 2026, the last confirmed observation). For
#   each start date we fit the exponential growth rate r (per day) and doubling
#   time (log(2)/r) over [start, END] two ways:
#
#     * incidence-based   -- log-linear fit of log(daily incidence) ~ day
#                            (this is what 02_npi_inputs_and_fiber_runs.R uses)
#     * cumulative-based  -- log-linear fit of log(cumulative) ~ day
#
#   WHY: 02 reports some NEGATIVE growth rates for confirmed cases, which feels
#   wrong. This shows it is expected: the confirmed daily incidence is noisy
#   (a late-May spike, then dips), so windows that start after the spike capture
#   a flat/declining incidence trend and give r <= 0. The CUMULATIVE curve is
#   monotone, so its log-linear slope is essentially always >= 0 -- hence the two
#   estimators diverge, and the incidence one (correctly) can be negative. Short
#   late-start windows are dominated by day-to-day noise and swing wildly; point
#   size in the figure encodes how many observations each fit used.
#
# OUTPUT (in outputs/)
#   dose_confirmed_growth_by_start.csv          -- the full table
#   dose_confirmed_growth_by_start.png          -- r + doubling time vs start
#   dose_confirmed_incidence_context.png        -- the incidence/cumulative curves
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

# Day axis (numeric) measured from the first confirmed date, shared by both fits.
d0 <- min(co$date)
to_day <- function(d) as.numeric(d - d0)

# ---- log-linear growth-rate helper -----------------------------------------
# r = slope of log(value) ~ day over the supplied points; NA if < 2 positive pts.
loglin_r <- function(day, value) {
  ok  <- is.finite(day) & is.finite(value) & value > 0
  day <- day[ok]; value <- value[ok]
  if (length(unique(day)) < 2L) return(c(r = NA_real_, n = length(value)))
  c(r = unname(coef(stats::lm(log(value) ~ day))[2]), n = length(value))
}

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

  fi <- loglin_r(to_day(inc_w$date), inc_w$incidence)
  fc <- loglin_r(to_day(cum_w$date), cum_w$national_cumulative_confirmed_cases)

  rows[[length(rows) + 1L]] <- data.frame(
    start_date = s, end_date = END_DATE,
    method     = c("incidence", "cumulative"),
    growth_rate   = c(fi[["r"]], fc[["r"]]),
    doubling_time = c(if (is.finite(fi[["r"]]) && fi[["r"]] > 0) log(2) / fi[["r"]] else NA_real_,
                      if (is.finite(fc[["r"]]) && fc[["r"]] > 0) log(2) / fc[["r"]] else NA_real_),
    n_points   = c(fi[["n"]], fc[["n"]]),
    row.names  = NULL)
}
tbl <- do.call(rbind, rows)
tbl$method <- factor(tbl$method, levels = c("incidence", "cumulative"))

write.csv(tbl, file.path(DIR_OUT, "dose_confirmed_growth_by_start.csv"), row.names = FALSE)

cat(sprintf("\nConfirmed-case growth rate by window start (END = %s):\n",
            format(END_DATE, "%d %b %Y")))
disp <- tbl
disp$growth_rate   <- round(disp$growth_rate, 4)
disp$doubling_time <- round(disp$doubling_time, 1)
print(disp[order(disp$method, disp$start_date), ], row.names = FALSE)

n_neg <- sum(tbl$method == "incidence" & is.finite(tbl$growth_rate) & tbl$growth_rate < 0)
cat(sprintf("\nIncidence-based fits with NEGATIVE growth: %d of %d start dates.\n",
            n_neg, sum(tbl$method == "incidence")))

# ---- figure 1: growth rate + doubling time vs start date --------------------
cols <- c(incidence = "#1a9850", cumulative = "#762a83")
long <- tidyr::pivot_longer(tbl, c(growth_rate, doubling_time),
                            names_to = "metric", values_to = "value")
long$metric <- factor(long$metric, levels = c("growth_rate", "doubling_time"),
                       labels = c("Growth rate r (per day)", "Doubling time (days)"))

p1 <- ggplot(long, aes(start_date, value, colour = method)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_line(linewidth = 0.5) +
  geom_point(aes(size = n_points), alpha = 0.8) +
  facet_wrap(~ metric, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = cols, name = "Fit on") +
  scale_size_continuous(range = c(0.8, 3.2), name = "n obs in window") +
  scale_x_date(date_breaks = "1 week", date_labels = "%d %b") +
  labs(title = "Confirmed-case growth rate vs window start date",
       subtitle = sprintf("Window = [start, %s]. Incidence-based r can be negative (declining/noisy tail); cumulative-based r stays >= 0 (monotone curve). Larger points use more data.",
                          format(END_DATE, "%d %b %Y")),
       x = "Window start date", y = NULL) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(DIR_OUT, "dose_confirmed_growth_by_start.png"), p1,
       width = 9, height = 7, dpi = 150)
print(p1)

# ---- figure 2: the confirmed incidence + cumulative (explains the negatives) -
p2 <- ggplot(inc, aes(date, incidence)) +
  geom_col(fill = "#1a9850", alpha = 0.55) +
  geom_smooth(se = FALSE, colour = "#0b3d0b", linewidth = 0.7, method = "loess",
              formula = y ~ x, span = 0.5) +
  scale_x_date(date_breaks = "1 week", date_labels = "%d %b") +
  labs(title = "Observed confirmed daily incidence",
       subtitle = "diff(cumulative)/diff(day). The late-May spike then dip is what drives negative incidence-based growth over windows that start afterwards.",
       x = "Date", y = "Confirmed cases per day") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(DIR_OUT, "dose_confirmed_incidence_context.png"), p2,
       width = 9, height = 4.5, dpi = 150)
print(p2)

message("\n04_confirmed_growth_by_start_date.R complete. Outputs in outputs/.")
