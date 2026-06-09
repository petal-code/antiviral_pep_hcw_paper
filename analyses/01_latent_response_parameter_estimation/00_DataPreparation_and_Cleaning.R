# ============================================================================
# 00_DataPreparation_and_Cleaning.R
# ----------------------------------------------------------------------------
# PURPOSE
#   Turn the raw inputs in data-raw/ into clean, analysis-ready objects in
#   data-processed/. Nothing here fits a model; this script only reads, cleans,
#   reconstructs and saves. The two fitting scripts (01, 02), the combine script
#   (03) and the two checking scripts all read ONLY from data-processed/.
#
# WHAT IT PRODUCES
#   Two bundled R objects (named lists), so data-processed/ stays tidy instead of
#   filling up with many loose CSVs. Every downstream script reads one of these:
#
#   WestAfrica_QCurve_PreppedData.rds   list(
#                   anchors = <cleaned West Africa literature anchors>
#                 )
#
#   DRC_QCurve_PreppedData.rds  list(
#                   anchors                   = <cleaned DRC literature anchors>,
#                   conflict_qseries          = <empirical "conflict" Q(t) curve>,
#                   conflict_plusplus_qseries = <conflict Q with a forced collapse
#                                                (success -> 0) over days 200-300>,
#                   no_conflict_qseries       = <counterfactual no-conflict Q(t):
#                                                monthly, plateaued, tail removed>,
#                   durations                 = <scenario horizons (max day)>,
#                   province_weekly_qc        = <per-province weekly reconstruction;
#                                                diagnostic only, not used in fits>
#                 )
#
#   The anchor tables and the Q series are documented column-by-column where they
#   are built, below.
#
# KEY IDEA - the "response-quality curve" Q(t)
#   For the DRC scenarios, Q(t) is NOT estimated. It is the empirical fraction of
#   safe-and-dignified burials (SDB) that were successful over time (the Warsame
#   et al. line-list), reconstructed here and then supplied as fixed data to
#   Model B. Q is a relative 0-1 index: Q = success(t) / max(success). The
#   community-unsafe-funeral parameter, by contrast, is kept on the ABSOLUTE
#   Warsame scale (1 - success(t)), because its floor is 1 - max(success), not 0.
# ============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
})

source(here::here("analyses", "01_latent_response_parameter_estimation", "helpers.R"))

dir.create(DIR_PROCESSED, showWarnings = FALSE, recursive = TRUE)

# ----------------------------------------------------------------------------
# Settings that define the SDB reconstruction
# These are the choices actually used to produce the published curves. They are
# written out explicitly (rather than buried as optional toggles) so the reader
# can see exactly how the empirical curve is built.
# ----------------------------------------------------------------------------

# Provinces whose province-specific SDB success curves are averaged to form the
# shared DRC response curve.
PROVINCES_TO_AVERAGE <- c("North Kivu", "Ituri")

# "Successful response" follows Warsame Fig. 3: a burial counts as a successful
# response if it was a literal success OR an SDB turned out not to be needed.
# The performance denominator adds failures. False alerts / unclear are excluded.
SUCCESS_OUTCOMES <- c("success", "sdb not needed")
FAILURE_OUTCOMES <- c("failure")

# Early, tiny-sample success spikes are artefacts. A weekly bin in the first
# INITIAL_SPIKE_MAX_DAY days, with fewer than INITIAL_SPIKE_MIN_ELIGIBLE eligible
# burials and an implausibly high success proportion, is set to zero success.
INITIAL_SPIKE_MAX_DAY        <- 75
INITIAL_SPIKE_MIN_ELIGIBLE   <- 10
INITIAL_SPIKE_SUCCESS_THRESH <- 0.50

# A bin must have at least this many eligible burials to contribute to Q.
MIN_ELIGIBLE_FOR_Q <- 1

# Rolling-average window (weeks) used to smooth the weekly conflict reconstruction.
ROLLING_WINDOW_WEEKS <- 4

# The DRC "++" scenario forces a complete response collapse (SDB success -> 0,
# and therefore community unsafe funerals -> 1) over this day window. Because the
# collapse is applied to the shared Q itself, every response parameter collapses
# with it.
PLUSPLUS_WINDOW_DAY        <- c(200, 300)
PLUSPLUS_SUCCESS_VALUE     <- 0
PLUSPLUS_UNSAFE_FUNERAL    <- 1

# The no-conflict counterfactual is truncated at the first data-rich maximum of
# the bridged success line. A bin counts as "data-rich" if it has at least this
# many eligible burials.
NO_CONFLICT_END_MIN_ELIGIBLE <- 25

# Anchors superseded by the line-list reconstruction (a community-unsafe-funeral
# summary that would otherwise double-count the SDB evidence, and a start anchor
# inconsistent with "begin at 1"). Held out of all DRC fits.
DRC_SUPERSEDED_ANCHORS <- c("DRC_UFC_00", "DRC_UFC_01")

# ----------------------------------------------------------------------------
# Locate the raw workbooks
# ----------------------------------------------------------------------------
curve_workbook <- resolve_input_file(
  "filovirus_three_scenario_curve_inputs_bestcase_recreate.xlsx",
  "curve-anchor workbook"
)
sdb_workbook <- resolve_input_file(
  c("evd_drc_sdb_performance_datasets_pub.xlsx",
    "evd_drc_sdb_performance_datasets_pub(1).xlsx",
    "evd_drc_sdb_performance_datasets_pub (1).xlsx"),
  "Warsame SDB workbook"
)

# ============================================================================
# PART 1 - Clean the literature anchor tables
# ============================================================================
# Both the West Africa and DRC sheets share the same layout (two header rows to
# skip, then one row per anchor). We keep only the rows flagged for fitting and
# coerce the columns to sensible types.

REQUIRED_ANCHOR_COLS <- c(
  "anchor_id", "parameter", "relative_day", "value_used",
  "fit_role", "include_in_fit", "weight", "direction",
  "lower_bound", "upper_bound"
)

# read_anchor_sheet() returns a cleaned anchor table: one row per literature data
# point kept for fitting, with these columns ---------------------------------
#   anchor_id      unique id of the anchor (e.g. "WA_UFC_01")
#   parameter      which response parameter it informs (one of PARAM_LEVELS)
#   relative_day   outbreak day the anchor refers to (days since response start)
#   value_used     the observed / derived value of that parameter at that day
#   fit_role       the anchor's role label from the workbook (e.g. "main_fit")
#   include_in_fit "YES" for rows used in the fit (others are filtered out here)
#   weight         relative trust in the anchor; higher = fit it more tightly
#                  (becomes obs_sd_mult = 1 / weight in the fitting scripts)
#   direction      "up" if the parameter improves upward over the response,
#                  "down" if it improves downward
#   lower_bound    workbook low end of the parameter -> prior CENTRE for that
#                  parameter's LOWER magnitude endpoint
#   upper_bound    workbook high end of the parameter -> prior CENTRE for the
#                  UPPER magnitude endpoint
read_anchor_sheet <- function(path, sheet) {
  raw <- read_excel(path, sheet = sheet, skip = 2)

  missing <- setdiff(REQUIRED_ANCHOR_COLS, names(raw))
  if (length(missing) > 0) {
    stop("Sheet '", sheet, "' is missing columns: ", paste(missing, collapse = ", "))
  }

  raw %>%
    transmute(
      anchor_id      = as.character(anchor_id),
      parameter      = as.character(parameter),
      relative_day   = as.numeric(relative_day),
      value_used     = as.numeric(value_used),
      fit_role       = trimws(as.character(fit_role)),
      include_in_fit = toupper(trimws(as.character(include_in_fit))),
      weight         = as.numeric(weight),
      direction      = tolower(trimws(as.character(direction))),
      lower_bound    = as.numeric(lower_bound),
      upper_bound    = as.numeric(upper_bound)
    ) %>%
    # Keep only the parameters we model, the rows marked for fitting, and rows
    # with a usable value and time.
    filter(parameter %in% PARAM_LEVELS,
           include_in_fit == "YES",
           !is.na(value_used),
           !is.na(relative_day))
}

wa_anchors  <- read_anchor_sheet(curve_workbook, "Worst_WestAfrica")

drc_anchors <- read_anchor_sheet(curve_workbook, "Middle_DRC_2018_2019") %>%
  # Drop the SDB-summary anchors that the line-list reconstruction supersedes.
  filter(!(anchor_id %in% DRC_SUPERSEDED_ANCHORS))

# (wa_anchors and drc_anchors are bundled into the .rds objects at the end.)
message("West Africa anchors per parameter:")
print(table(wa_anchors$parameter))
message("DRC anchors per parameter:")
print(table(drc_anchors$parameter))

# ============================================================================
# PART 2 - Reconstruct the empirical SDB success time series
# ============================================================================

# ---- 2a. Read and filter the SDB line-list ---------------------------------
admin_units <- read_excel(sdb_workbook, sheet = "admin_units") %>%
  mutate(hz = as.character(hz), province = as.character(province)) %>%
  distinct(hz, province)

sdb_linelist <- read_excel(sdb_workbook, sheet = "sdb_dataset") %>%
  mutate(
    hz            = as.character(hz),
    outcome_lshtm = tolower(trimws(as.character(outcome_lshtm))),
    origin_cat    = tolower(trimws(as.character(origin_cat))),
    # readxl usually returns date-formatted cells as Date already; guard against
    # the case where the column arrives as a raw Excel serial number.
    date          = if (is.numeric(date)) as.Date(date, origin = "1899-12-30") else as.Date(date)
  ) %>%
  left_join(admin_units, by = "hz") %>%
  filter(province %in% PROVINCES_TO_AVERAGE) %>%
  # Some line-list rows carry epi_year/epi_week but no actual date. If kept, they
  # form NA-date bins that can spuriously become the empirical maximum used to
  # scale Q. Drop them so the scaling maximum matches the dated, plotted curve.
  filter(!is.na(date))

if (nrow(sdb_linelist) == 0) stop("No SDB rows after filtering to the requested provinces.")

# The outbreak clock for the SDB series starts at the first dated burial.
sdb_start_date <- min(sdb_linelist$date, na.rm = TRUE)

# ---- 2b. Bin into time periods and compute per-province success ------------
# aggregation_unit: "epi_week" (Monday-start weeks; recreates the step-like
# Warsame curve) for the conflict scenarios, or "monthly" for the smoother
# no-conflict counterfactual.
bin_provinces <- function(linelist, aggregation_unit) {
  # Assign each burial to a time bin and timestamp it with the bin's mid-date.
  if (aggregation_unit == "epi_week") {
    # Monday-start week; use the mid-week (Thursday) as the bin date.
    monday <- linelist$date - (as.integer(format(linelist$date, "%u")) - 1L)
    linelist$bin_mid_date <- monday + 3L
  } else if (aggregation_unit == "monthly") {
    linelist$bin_mid_date <- as.Date(format(linelist$date, "%Y-%m-01")) + 14L
  } else {
    stop("aggregation_unit must be 'epi_week' or 'monthly'.")
  }

  binned <- linelist %>%
    group_by(province, bin_mid_date) %>%
    summarise(
      n_successful_response = sum(outcome_lshtm %in% SUCCESS_OUTCOMES, na.rm = TRUE),
      n_failure             = sum(outcome_lshtm %in% FAILURE_OUTCOMES, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      relative_day = as.numeric(bin_mid_date - sdb_start_date),
      n_eligible   = n_successful_response + n_failure,
      prop_success = if_else(n_eligible > 0, n_successful_response / n_eligible, NA_real_)
    ) %>%
    # Suppress early tiny-sample success spikes (set them to zero success).
    mutate(
      is_initial_spike = relative_day <= INITIAL_SPIKE_MAX_DAY &
                         n_eligible   <  INITIAL_SPIKE_MIN_ELIGIBLE &
                         prop_success >= INITIAL_SPIKE_SUCCESS_THRESH,
      prop_success = if_else(is_initial_spike, 0, prop_success)
    ) %>%
    filter(!is.na(prop_success)) %>%
    arrange(province, relative_day)

  binned
}

# Per-province "no-conflict bridge": each province rises to its data-rich maximum
# success and then plateaus (holds flat). Removing the later decline is what
# turns the observed curve into a no-conflict counterfactual.
apply_province_plateau <- function(binned, min_eligible) {
  binned %>%
    group_by(province) %>%
    arrange(relative_day, .by_group = TRUE) %>%
    group_modify(function(g, ...) {
      rich <- g$n_eligible >= min_eligible & is.finite(g$prop_success)
      if (!any(rich)) { g$prop_success_plateau <- g$prop_success; return(g) }
      max_success <- max(g$prop_success[rich], na.rm = TRUE)
      first_max   <- min(which(g$prop_success >= max_success - 1e-8 & rich))
      pv <- g$prop_success
      pv[first_max:length(pv)] <- max_success
      g$prop_success_plateau <- clip01(pv)
      g
    }) %>%
    ungroup()
}

# Average the province curves (unweighted mean) into the single shared line.
average_provinces <- function(binned, value_col) {
  binned %>%
    group_by(bin_mid_date) %>%
    summarise(
      relative_day  = first(relative_day),
      success_avg   = mean(.data[[value_col]], na.rm = TRUE),
      n_eligible_sum = sum(n_eligible, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(relative_day) %>%
    filter(n_eligible_sum >= MIN_ELIGIBLE_FOR_Q, is.finite(success_avg))
}

# ---- 2c. Turn an averaged success line into a finalised Q series -----------
# Steps: (optionally smooth), prepend a forced Q(0)=0 start, max-scale to get the
# relative 0-1 Q, and compute the absolute community-unsafe-funeral proxy.
finalise_q_series <- function(avg, smooth, drc_anchor_max_day,
                              collapse_window = NULL, end_day = NULL) {

  s <- avg %>% arrange(relative_day)

  # Smooth (centred 4-week rolling mean) or pass through unchanged.
  s$success_smoothed <- if (smooth) {
    sm <- clip01(rolling_mean_centered(s$success_avg, k = ROLLING_WINDOW_WEEKS))
    sm[!is.finite(sm)] <- s$success_avg[!is.finite(sm)]
    sm
  } else {
    s$success_avg
  }

  # "++" collapse: force success to zero across the collapse window. Applied
  # before normalisation so the shared Q (and therefore every parameter) drops.
  if (!is.null(collapse_window)) {
    in_win <- s$relative_day >= collapse_window[1] & s$relative_day <= collapse_window[2]
    s$success_smoothed[in_win] <- PLUSPLUS_SUCCESS_VALUE
  }

  # No-conflict truncation: cut the conflict-deterioration tail at the supplied
  # end day, then renormalise the horizon to that endpoint.
  if (!is.null(end_day)) {
    s <- s %>% filter(relative_day <= end_day)
  }

  # The shared horizon spans the workbook anchors and the SDB series.
  max_day <- max(c(drc_anchor_max_day, s$relative_day), na.rm = TRUE)

  # Prepend the forced outbreak-start support point: success = 0 at day 0.
  s <- bind_rows(
    tibble(relative_day = 0, success_smoothed = 0, n_eligible_sum = 0),
    s %>% select(relative_day, success_smoothed, n_eligible_sum)
  ) %>%
    arrange(relative_day) %>%
    distinct(relative_day, .keep_all = TRUE)

  scale_max <- max(s$success_smoothed, na.rm = TRUE)
  if (!is.finite(scale_max) || scale_max <= 0) stop("Q scaling maximum is not positive.")

  # The returned Q series has one row per support point (the forced day-0 start
  # first, then each retained SDB bin), with these columns -------------------
  #   relative_day               outbreak day of the support point
  #   tau_q                      relative_day / max_day  (normalised time, 0-1)
  #   q_value                    the shared response-quality index Q in [0,1]
  #                              (= success / max success); supplied to Model B
  #                              as DATA and used as the scenario's q_value
  #   unsafe_funeral_comm_proxy  ABSOLUTE community unsafe-funeral level (1 - success);
  #                              floor is 1 - max(success), not 0
  #   success_smoothed           the smoothed SDB success proportion behind Q
  #   n_eligible_sum             eligible burials behind the bin (0 at the start row)
  s %>%
    mutate(
      tau_q = relative_day / max_day,
      # Q is RELATIVE: divide by the maximum only (no min subtraction). The best
      # achieved response is Q = 1; later dips stay partial unless success hits 0.
      q_value = clip01(success_smoothed / scale_max),
      # The community unsafe-funeral proxy is ABSOLUTE: 1 - success (its floor is
      # 1 - max(success), not 0).
      unsafe_funeral_comm_proxy = clip01(1 - success_smoothed)
    ) %>%
    select(relative_day, tau_q, q_value, unsafe_funeral_comm_proxy,
           success_smoothed, n_eligible_sum)
}

# ---- 2d. Build the three DRC Q series --------------------------------------
drc_anchor_max_day <- max(drc_anchors$relative_day, na.rm = TRUE)

# CONFLICT and CONFLICT++ : weekly reconstruction, no plateau.
weekly_binned <- bin_provinces(sdb_linelist, "epi_week")
weekly_avg    <- average_provinces(weekly_binned, "prop_success")

drc_conflict_qseries <- finalise_q_series(
  weekly_avg, smooth = TRUE, drc_anchor_max_day = drc_anchor_max_day
)

drc_conflict_plusplus_qseries <- finalise_q_series(
  weekly_avg, smooth = TRUE, drc_anchor_max_day = drc_anchor_max_day,
  collapse_window = PLUSPLUS_WINDOW_DAY
)

# NO-CONFLICT : monthly reconstruction, per-province plateau, conflict tail
# removed at the first data-rich maximum of the averaged plateaued line.
monthly_binned <- bin_provinces(sdb_linelist, "monthly") %>%
  apply_province_plateau(min_eligible = NO_CONFLICT_END_MIN_ELIGIBLE)
monthly_avg    <- average_provinces(monthly_binned, "prop_success_plateau")

no_conflict_end_day <- monthly_avg %>%
  filter(n_eligible_sum >= NO_CONFLICT_END_MIN_ELIGIBLE) %>%
  { if (nrow(.) == 0) stop("No data-rich monthly bins to define the no-conflict endpoint.")
    .$relative_day[which.max(.$success_avg)] }

drc_no_conflict_qseries <- finalise_q_series(
  monthly_avg, smooth = FALSE, drc_anchor_max_day = drc_anchor_max_day,
  end_day = no_conflict_end_day
)

# ---- 2e. Bundle everything into two tidy .rds objects ----------------------
# drc_durations: one row per scenario giving the largest day on its curve. The
# ratio duration(conflict) / duration(no-conflict) is used in 03 to stretch the
# West-Africa-with-conflict timeline (conflict drags the response out over a
# longer window). The wa_anchor_max_day row is carried for context only.
drc_durations <- tibble(
  scenario = c("drc_conflict", "drc_conflict_plusplus", "drc_no_conflict", "wa_anchor_max_day"),
  max_day  = c(
    max(drc_conflict_qseries$relative_day),
    max(drc_conflict_plusplus_qseries$relative_day),
    max(drc_no_conflict_qseries$relative_day),
    max(wa_anchors$relative_day, na.rm = TRUE)
  )
)

# Save exactly two bundled lists (see the header for their full contents). Each
# element is a data frame documented where it was built above; the trailing
# comments say which downstream script consumes it.
wa_prep <- list(
  anchors = wa_anchors                          # cleaned WA anchors   (01, west_africa_checking)
)
saveRDS(wa_prep, file.path(DIR_PROCESSED, "WestAfrica_QCurve_PreppedData.rds"))

drc_prep <- list(
  anchors                   = drc_anchors,                   # cleaned DRC anchors          (02, no-conflict check)
  conflict_qseries          = drc_conflict_qseries,          # shared Q, "conflict"         (02)
  conflict_plusplus_qseries = drc_conflict_plusplus_qseries, # shared Q, "conflict++"       (02)
  no_conflict_qseries       = drc_no_conflict_qseries,       # shared Q, no-conflict        (no-conflict check)
  durations                 = drc_durations,                 # scenario horizons            (03 time-stretch)
  province_weekly_qc        = weekly_binned                  # per-province reconstruction  (diagnostic only)
)
saveRDS(drc_prep, file.path(DIR_PROCESSED, "DRC_QCurve_PreppedData.rds"))

message("\nDRC scenario horizons (max day):")
print(drc_durations)
message("\n00_DataPreparation_and_Cleaning.R complete. ",
        "Wrote data-processed/WestAfrica_QCurve_PreppedData.rds and data-processed/DRC_QCurve_PreppedData.rds")
