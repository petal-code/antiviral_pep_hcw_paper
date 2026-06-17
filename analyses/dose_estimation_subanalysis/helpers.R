# ============================================================================
# helpers.R  (dose_estimation_subanalysis)
# ----------------------------------------------------------------------------
# Small shared utilities and path constants for the two scripts in this
# subanalysis:
#   01_fit_dose_q_curve.R              -- fits + extrapolates the dose Q curve
#   02_npi_inputs_and_fiber_runs.R     -- turns the Q curve into time-varying
#                                         NPI inputs and runs fiber across R0
#
# Paths are resolved from the repository root with here::here() so the scripts
# run the same regardless of the working directory they are launched from.
# ============================================================================

# ---- Project paths ---------------------------------------------------------
ANALYSIS_DIR <- here::here("analyses", "dose_estimation_subanalysis")
DIR_STAN     <- file.path(ANALYSIS_DIR, "stan-models")
DIR_OUT      <- file.path(ANALYSIS_DIR, "outputs")

# Ensure the outputs folder exists (it is created on first run; the generated
# .rds / .csv artefacts are regenerable and are not tracked in git).
dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)

# ---- The raw dose-coverage observations ------------------------------------
# Five (calendar date, percentage) pairs. Percentages are converted to
# proportions ( /100 ) downstream. The first point is 0% on 18 May 2026.
DOSE_OBS <- data.frame(
  date       = as.Date(c("2026-05-18", "2026-05-24", "2026-05-31",
                         "2026-06-07", "2026-06-14")),
  percentage = c(0.0, 19.3, 30.2, 64.4, 63.1),
  stringsAsFactors = FALSE
)

# ---- Tiny numeric helper ----------------------------------------------------
clip01 <- function(x) pmin(1, pmax(0, x))

# ---- Build the relative-day observation table (with optional front padding) -
# Converts the calendar-date observations to a RELATIVE-DAY axis measured from
# `start_date` (day 0), and turns percentages into proportions.
#
# Front padding ("set a start date and make all days between the start date and
# 18 May set to 0"): if `start_date` is earlier than the first observation, one
# zero-valued point is added for every day from day 0 up to (but not including)
# the first real observation. These extra zeros anchor the early curve at the 0%
# floor. With start_date == the first observation date (the default) there is no
# padding and the relative days are simply {0, 6, 13, 20, 27}.
#
# Returns a data.frame with columns: date, relative_day, proportion, padded.
build_dose_obs <- function(dose_obs = DOSE_OBS,
                           start_date = min(dose_obs$date),
                           pad_step = 1L) {
  start_date <- as.Date(start_date)
  first_obs  <- min(dose_obs$date)
  if (start_date > first_obs) {
    stop("`start_date` (", start_date, ") must be on or before the first ",
         "observation (", first_obs, ").", call. = FALSE)
  }

  obs <- data.frame(
    date         = dose_obs$date,
    relative_day = as.integer(dose_obs$date - start_date),
    proportion   = dose_obs$percentage / 100,
    padded       = FALSE,
    stringsAsFactors = FALSE
  )

  offset <- as.integer(first_obs - start_date)   # days of padding requested
  if (offset >= 1L) {
    pad_days <- seq.int(0L, offset - 1L, by = pad_step)
    pad <- data.frame(
      date         = start_date + pad_days,
      relative_day = pad_days,
      proportion   = 0,
      padded       = TRUE,
      stringsAsFactors = FALSE
    )
    obs <- rbind(pad, obs)
  }

  obs[order(obs$relative_day), , drop = FALSE]
}
