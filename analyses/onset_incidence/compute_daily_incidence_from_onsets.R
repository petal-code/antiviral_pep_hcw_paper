# ============================================================================
# compute_daily_incidence_from_onsets.R
# ----------------------------------------------------------------------------
# PURPOSE
#   Derive a daily incidence curve from the cumulative-onset quantile bands in
#   data-raw/onsets_over_time.csv. Three steps, in the order requested:
#
#     1. cumulative_onsets = midpoint of the cumulative_onsets_lower_30 and
#        cumulative_onsets_lower_60 columns, i.e.
#        (cumulative_onsets_lower_30 + cumulative_onsets_lower_60) / 2.
#
#     2. Shift the date backwards in time by SHIFT_DAYS (= 7) days, so each
#        cumulative value is re-attributed to a date 7 days earlier.
#
#     3. daily_incidence = the day-to-day difference of cumulative_onsets
#        (cumulative[t] - cumulative[t-1]), which is the incidence implied by
#        the cumulative curve.
#
# OUTPUT
#   data-processed/onsets_daily_incidence.csv with columns:
#     date              -- original date minus SHIFT_DAYS
#     cumulative_onsets -- midpoint of the two lower bounds (step 1)
#     daily_incidence   -- day-to-day difference of cumulative_onsets (step 3)
#
# NOTE on the first day
#   The day-to-day difference is undefined for the very first date (no previous
#   day). The cumulative series starts from ~0 (the first cumulative value is
#   ~3e-6), so we treat the pre-series cumulative as 0 and set the first day's
#   incidence equal to its own cumulative value. This makes the daily incidences
#   sum back to the final cumulative value. To instead leave the first day
#   blank, change `default = 0` to `default = NA_real_` in the lag() call below.
#
# RUN (from the repo root)
#   Rscript analyses/onset_incidence/compute_daily_incidence_from_onsets.R
# ============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

# ---- config -----------------------------------------------------------------
SHIFT_DAYS  <- 7L
INPUT_CSV   <- here::here("data-raw",       "onsets_over_time.csv")
OUTPUT_CSV  <- here::here("data-processed", "onsets_daily_incidence.csv")

# ---- read -------------------------------------------------------------------
onsets <- read_csv(
  INPUT_CSV,
  col_types = cols(date = col_date(format = "%Y-%m-%d"), .default = col_double())
)

# ---- transform --------------------------------------------------------------
daily <- onsets %>%
  arrange(date) %>%                                  # ensure chronological order
  mutate(
    # Step 1: midpoint of the two lower bounds.
    cumulative_onsets = (cumulative_onsets_lower_30 + cumulative_onsets_lower_60) / 2,
    # Step 2: throw the series back in time by SHIFT_DAYS days.
    date = date - SHIFT_DAYS
  ) %>%
  arrange(date) %>%                                  # order is unchanged, but be safe
  mutate(
    # Step 3: daily incidence = difference between subsequent days. The first
    # day takes the cumulative value itself (pre-series cumulative assumed 0).
    daily_incidence = cumulative_onsets - lag(cumulative_onsets, default = 0)
  ) %>%
  select(date, cumulative_onsets, daily_incidence)

# ---- write ------------------------------------------------------------------
dir.create(dirname(OUTPUT_CSV), showWarnings = FALSE, recursive = TRUE)
write_csv(daily, OUTPUT_CSV)

# ---- report -----------------------------------------------------------------
message("Wrote ", nrow(daily), " rows to ", OUTPUT_CSV)
message("Date range (after ", SHIFT_DAYS, "-day shift): ",
        min(daily$date), " to ", max(daily$date))
print(head(daily))
print(tail(daily))
