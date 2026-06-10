# ============================================================================
# 00b_DataPreparation_CommunityDeaths.R
# ----------------------------------------------------------------------------
# Prep the DRC community-death-proportion series for Model C (reach estimation).
#
# Source = data-processed/ebola_drc_community_deaths.csv (one row per reporting
# period; community_deaths_pct = community deaths / all deaths, with numerator /
# denominator_n where reported).
#
# Rules (per CW):
#   * where numerator + denominator_n are reported, use them directly;
#   * where only the % is reported, impute the denominator as the MEAN of the
#     reported denominators, and back out numerator = round(pct * denom).
#
# Output: a tidy (relative_day, n_comm, N_deaths) table for the binomial reach
# likelihood, saved to data-processed/DRC_QCurve/DRC_CommunityDeaths_Prepped.rds.
#
# >>> CONFIRM <<< DAY0 (the calendar date that the scenario's relative_day = 0
#     maps to). It MUST match the day-0 anchor used to build the scenario curves,
#     or the reach curve will be time-shifted relative to everything else.
# ============================================================================

suppressPackageStartupMessages({ library(dplyr); library(readr); library(stringr) })
source(here::here("analyses", "01_latent_response_parameter_estimation", "helpers.R"))

# ---- day-0 anchor (EDIT/CONFIRM) -------------------------------------------
DAY0 <- as.Date("2018-08-01")   # DRC outbreak declaration; confirm vs scenario day 0

csv_path <- file.path(DIR_PROCESSED, "ebola_drc_community_deaths.csv")
stopifnot(file.exists(csv_path))
raw <- readr::read_csv(csv_path, show_col_types = FALSE)

cd <- raw %>%
  mutate(
    pct  = as.numeric(str_remove(community_deaths_pct, "%")) / 100,
    date = as.Date(date_midpoint, format = "%m/%d/%y")
  ) %>%
  filter(!is.na(pct), !is.na(date))

# Impute missing denominators with the mean of the reported ones, then back out
# the numerator from the percentage where it is missing.
denom_mean <- round(mean(cd$denominator_n, na.rm = TRUE))
message("Mean reported denominator (used for imputation): ", denom_mean)

cd <- cd %>%
  mutate(
    N_deaths = if_else(is.na(denominator_n), denom_mean, as.numeric(denominator_n)),
    n_comm   = if_else(is.na(numerator), round(pct * N_deaths), as.numeric(numerator)),
    # keep counts internally consistent (n_comm <= N_deaths)
    n_comm   = pmin(n_comm, N_deaths),
    relative_day = as.integer(date - DAY0)
  ) %>%
  filter(relative_day >= 0) %>%
  transmute(relative_day, date, n_comm = as.integer(n_comm), N_deaths = as.integer(N_deaths),
            pct_reported = pct) %>%
  arrange(relative_day)

# Empirical worst/best community-death proportions -> reach-clock normalisation
# constants c_hi (worst) and c_lo (best). Used by Model C (4a, normalised reach).
c_hi <- max(cd$n_comm / cd$N_deaths)
c_lo <- min(cd$n_comm / cd$N_deaths)
message(sprintf("Observed community-death proportion range: c_lo=%.3f .. c_hi=%.3f (n=%d periods)",
                c_lo, c_hi, nrow(cd)))

out <- list(obs = cd, c_hi = c_hi, c_lo = c_lo, day0 = DAY0)
dir.create(file.path(DIR_PROCESSED, "DRC_QCurve"), recursive = TRUE, showWarnings = FALSE)
saveRDS(out, file.path(DIR_PROCESSED, "DRC_QCurve", "DRC_CommunityDeaths_Prepped.rds"))
message("Wrote DRC_CommunityDeaths_Prepped.rds  (", nrow(cd), " periods)")

print(cd, n = 50)
