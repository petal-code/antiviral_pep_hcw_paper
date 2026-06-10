# ============================================================================
# 00_DataPreparation_and_Cleaning(revisedMethod).R
# ----------------------------------------------------------------------------
# PURPOSE
#   Turn the raw inputs in data-raw/ into clean, analysis-ready objects in
#   data-processed/. Nothing here fits a model; this script only reads, cleans,
#   reconstructs and saves. The two fitting scripts (01, 02) and the combine
#   script (03) all read ONLY from data-processed/.
#
# HOW THIS DIFFERS FROM THE ORIGINAL-METHODOLOGY 00
#   The SDB reconstruction (PART 2) is IDENTICAL - the Warsame line-list is the
#   same for both methodologies. What differs is PART 1, the literature anchors:
#   the REVISED methodology reads its anchors and per-parameter literature ranges
#   from the parameter-table workbook
#       filovirus_model_parameter_table_4_scenarios_with_etu_baseline.xlsx
#   (sheets "Worst West Africa" and "DRC conflict-smoothed"), exactly as the
#   github_upload revised-methodology scripts do. That workbook is the authoritative
#   source for the revised endpoints (e.g. the DRC IPC/PPE q-scaling range 0.071-
#   0.746 is its DRC "ipc_helper" summary range). The original methodology instead
#   reads filovirus_three_scenario_curve_inputs_bestcase_recreate.xlsx; the two
#   workbooks carry different anchor sets, so the locked endpoints differ.
#
#   The parameter-table workbook has a very different LAYOUT from the flat anchor
#   sheets, so PART 1 below has its own parser (read_parameter_table_sheet):
#     * a "Time-varying response" section gives each parameter's summary
#       literature low-high range (-> lower_bound / upper_bound, the endpoint
#       fallbacks and, for DRC IPC, the q-scaling range); and
#     * a "Time-varying curve anchors" section lists one row per orange-point
#       anchor, with the day:value encoded in the "Value / range" cell and the
#       fit role / parameter encoded in the description.
#   Parameters are identified by the clean Symbol column (d_hosp(t), p_hosp(t), ...).
#
# WHAT IT PRODUCES  (same bundles + paths as before, revised-method subfolders)
#   WestAfrica_QCurve_revisedMethod/WestAfrica_QCurve_PreppedData(revisedMethod).rds
#                 list(anchors = <cleaned WA anchors>)
#   DRC_QCurve_revisedMethod/DRC_QCurve_PreppedData(revisedMethod).rds  list(
#                 anchors, conflict_qseries, conflict_plusplus_qseries,
#                 no_conflict_qseries, durations, province_weekly_qc)
#
# KEY IDEA - the "response-quality curve" Q(t)  (unchanged from the original 00)
#   For the DRC scenarios Q(t) is NOT estimated; it is the empirical fraction of
#   safe-and-dignified burials that were successful over time (the Warsame
#   line-list), reconstructed here and (in the revised methodology) mapped onto
#   each parameter's LOCKED endpoints. Community unsafe funerals stay on the
#   ABSOLUTE Warsame scale (1 - success), floor 1 - max(success).
# ============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(tibble)
})

source(here::here("analyses", "01_latent_response_parameter_estimation_revisedMethodology",
                  "helpers(revisedMethod).R"))

dir.create(DIR_PROCESSED, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(DIR_PROCESSED, "WestAfrica_QCurve_revisedMethod"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(DIR_PROCESSED, "DRC_QCurve_revisedMethod"),        showWarnings = FALSE, recursive = TRUE)

# ----------------------------------------------------------------------------
# Settings that define the SDB reconstruction  (identical to the original 00)
# ----------------------------------------------------------------------------
PROVINCES_TO_AVERAGE <- c("North Kivu", "Ituri")
SUCCESS_OUTCOMES <- c("success", "sdb not needed")
FAILURE_OUTCOMES <- c("failure")
INITIAL_SPIKE_MAX_DAY        <- 75
INITIAL_SPIKE_MIN_ELIGIBLE   <- 10
INITIAL_SPIKE_SUCCESS_THRESH <- 0.50
MIN_ELIGIBLE_FOR_Q <- 1
ROLLING_WINDOW_WEEKS <- 4
PLUSPLUS_WINDOW_DAY        <- c(200, 300)
PLUSPLUS_SUCCESS_VALUE     <- 0
PLUSPLUS_UNSAFE_FUNERAL    <- 1
NO_CONFLICT_END_MIN_ELIGIBLE <- 25

# ----------------------------------------------------------------------------
# Locate the raw workbooks
# ----------------------------------------------------------------------------
# REVISED methodology: anchors come from the parameter-table workbook.
parameter_workbook <- resolve_input_file(
  "filovirus_model_parameter_table_4_scenarios_with_etu_baseline.xlsx",
  "parameter-table workbook"
)
sdb_workbook <- resolve_input_file(
  c("evd_drc_sdb_performance_datasets_pub.xlsx",
    "evd_drc_sdb_performance_datasets_pub(1).xlsx",
    "evd_drc_sdb_performance_datasets_pub (1).xlsx"),
  "Warsame SDB workbook"
)

# ============================================================================
# PART 1 - Read the literature anchors from the parameter-table workbook
# ============================================================================

# Each of the six modelled parameters is identified by its Symbol-column tag in
# the workbook (the cleanest, most stable key). p_UF,ETU(t) (ETU funerals) is
# intentionally omitted - it is fixed at 0 and not modelled.
SYMBOL_MAP <- c(
  "d_hosp(t)"    = "delay_hosp",
  "p_hosp(t)"    = "p_hosp",
  "p_ETU(t)"     = "p_ETU",
  "I_IPC(t)"     = "latent_IPC",
  "p_UF,comm(t)" = "p_unsafe_funeral_comm",
  "p_UF,hosp(t)" = "p_unsafe_funeral_hosp"
)

# Direction each parameter improves over the response (revised-methodology
# convention: "up" = increasing toward a better response, "down" = decreasing).
DIRECTION_MAP <- c(
  delay_hosp            = "down",
  p_hosp                = "up",
  p_ETU                 = "up",
  latent_IPC            = "up",
  p_unsafe_funeral_comm = "down",
  p_unsafe_funeral_hosp = "down"
)

# Parse a numeric (low, high) pair out of a "low-high [unit]" summary range cell,
# e.g. "0.324-0.905" or "1.843-6.316 days" or "0-0". The range separator is an
# en-dash; split on it first, then read the first number from each side (so a
# leading minus in subtraction-looking text is never mistaken for a negative).
extract_range <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "[−—]", "–")     # minus / em-dash -> en-dash
  first_num <- function(s) {
    n <- str_extract_all(s, "[0-9]*\\.?[0-9]+")[[1]]
    if (length(n) == 0) NA_real_ else as.numeric(n[1])
  }
  if (str_detect(x, "–")) {
    parts <- str_split_fixed(x, "–", 2)
    return(c(first_num(parts[1]), first_num(parts[2])))
  }
  if (str_detect(x, "[0-9]\\s*-\\s*[0-9]")) {
    parts <- str_split_fixed(x, "\\s*-\\s*", 2)
    return(c(first_num(parts[1]), first_num(parts[2])))
  }
  v <- first_num(x)
  c(v, v)
}

# Parse "relative day NNN: VVV [unit]" -> (relative_day, value_used).
parse_anchor_value <- function(x) {
  x   <- as.character(x)
  day <- as.numeric(str_match(x, "relative\\s+day\\s+([0-9.]+)\\s*:")[, 2])
  val <- as.numeric(str_match(x, ":\\s*([0-9.]+)")[, 2])
  list(relative_day = day, value_used = val)
}

# Parse a "key: value;" field out of an anchor's description (e.g. fit role).
parse_field <- function(x, key) {
  trimws(str_match(as.character(x), paste0(key, ":\\s*([^;]+)"))[, 2])
}

# read_parameter_table_sheet() returns a cleaned anchor table with EXACTLY the
# columns the rest of the revised pipeline expects (the same schema the original
# methodology's read_anchor_sheet produced), so 01/02/lock_endpoints are unchanged:
#   anchor_id, parameter (internal name), relative_day, value_used, fit_role,
#   weight, direction, lower_bound, upper_bound
# The parameter-table format carries no per-anchor weight, so weight is uniform 1
# (this reproduces the revised reference, which scales observation noise by the
# endpoint span only, not by anchor weights).
read_parameter_table_sheet <- function(path, sheet) {
  raw <- read_excel(path, sheet = sheet, col_names = FALSE)
  nm  <- c("Section", "Parameter", "Symbol", "Description", "Value_range", "Reference", "URL")
  names(raw)[seq_along(nm)] <- nm
  raw <- raw %>% mutate(across(everything(), as.character))

  # Per-parameter summary literature range (the "Time-varying response" rows).
  summary_ranges <- raw %>%
    filter(Section == "Time-varying response", Symbol %in% names(SYMBOL_MAP)) %>%
    mutate(
      parameter   = unname(SYMBOL_MAP[Symbol]),
      rng         = lapply(Value_range, extract_range),
      lower_bound = vapply(rng, function(z) min(z, na.rm = TRUE), numeric(1)),
      upper_bound = vapply(rng, function(z) max(z, na.rm = TRUE), numeric(1))
    ) %>%
    select(parameter, lower_bound, upper_bound) %>%
    distinct(parameter, .keep_all = TRUE)

  # One row per orange-point anchor (the "Time-varying curve anchors" rows).
  anchors <- raw %>%
    filter(Section == "Time-varying curve anchors",
           Symbol %in% names(SYMBOL_MAP),
           !is.na(Parameter), Parameter != "Orange-point fitted anchors") %>%
    mutate(
      anchor_id    = Parameter,
      parameter    = unname(SYMBOL_MAP[Symbol]),
      fit_role     = parse_field(Description, "fit role"),
      av           = lapply(Value_range, parse_anchor_value),
      relative_day = vapply(av, function(z) z$relative_day, numeric(1)),
      value_used   = vapply(av, function(z) z$value_used,   numeric(1))
    ) %>%
    filter(is.finite(relative_day), is.finite(value_used)) %>%
    left_join(summary_ranges, by = "parameter") %>%
    mutate(direction = unname(DIRECTION_MAP[parameter]), weight = 1) %>%
    select(anchor_id, parameter, relative_day, value_used, fit_role,
           weight, direction, lower_bound, upper_bound)

  if (nrow(anchors) == 0) stop("No usable anchors parsed from sheet '", sheet, "'.")
  anchors
}

wa_anchors  <- read_parameter_table_sheet(parameter_workbook, "Worst West Africa")
drc_anchors <- read_parameter_table_sheet(parameter_workbook, "DRC conflict-smoothed")

message("West Africa anchors per parameter:")
print(table(wa_anchors$parameter))
message("DRC anchors per parameter:")
print(table(drc_anchors$parameter))

# ============================================================================
# PART 2 - Reconstruct the empirical SDB success time series
# (IDENTICAL to the original-methodology 00: same Warsame line-list, same rules)
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
    date          = if (is.numeric(date)) as.Date(date, origin = "1899-12-30") else as.Date(date)
  ) %>%
  left_join(admin_units, by = "hz") %>%
  filter(province %in% PROVINCES_TO_AVERAGE) %>%
  filter(!is.na(date))

if (nrow(sdb_linelist) == 0) stop("No SDB rows after filtering to the requested provinces.")

sdb_start_date <- min(sdb_linelist$date, na.rm = TRUE)

# ---- 2b. Bin into time periods and compute per-province success ------------
bin_provinces <- function(linelist, aggregation_unit) {
  if (aggregation_unit == "epi_week") {
    monday <- linelist$date - (as.integer(format(linelist$date, "%u")) - 1L)
    linelist$bin_mid_date <- monday + 3L
  } else if (aggregation_unit == "monthly") {
    linelist$bin_mid_date <- as.Date(format(linelist$date, "%Y-%m-01")) + 14L
  } else {
    stop("aggregation_unit must be 'epi_week' or 'monthly'.")
  }

  linelist %>%
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
    mutate(
      is_initial_spike = relative_day <= INITIAL_SPIKE_MAX_DAY &
                         n_eligible   <  INITIAL_SPIKE_MIN_ELIGIBLE &
                         prop_success >= INITIAL_SPIKE_SUCCESS_THRESH,
      prop_success = if_else(is_initial_spike, 0, prop_success)
    ) %>%
    filter(!is.na(prop_success)) %>%
    arrange(province, relative_day)
}

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

average_provinces <- function(binned, value_col) {
  binned %>%
    group_by(bin_mid_date) %>%
    summarise(
      relative_day   = first(relative_day),
      success_avg    = mean(.data[[value_col]], na.rm = TRUE),
      n_eligible_sum = sum(n_eligible, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(relative_day) %>%
    filter(n_eligible_sum >= MIN_ELIGIBLE_FOR_Q, is.finite(success_avg))
}

# ---- 2c. Turn an averaged success line into a finalised Q series -----------
finalise_q_series <- function(avg, smooth, drc_anchor_max_day,
                              collapse_window = NULL, end_day = NULL) {

  s <- avg %>% arrange(relative_day)

  s$success_smoothed <- if (smooth) {
    sm <- clip01(rolling_mean_centered(s$success_avg, k = ROLLING_WINDOW_WEEKS))
    sm[!is.finite(sm)] <- s$success_avg[!is.finite(sm)]
    sm
  } else {
    s$success_avg
  }

  if (!is.null(collapse_window)) {
    in_win <- s$relative_day >= collapse_window[1] & s$relative_day <= collapse_window[2]
    s$success_smoothed[in_win] <- PLUSPLUS_SUCCESS_VALUE
  }

  if (!is.null(end_day)) {
    s <- s %>% filter(relative_day <= end_day)
  }

  max_day <- max(c(drc_anchor_max_day, s$relative_day), na.rm = TRUE)

  s <- bind_rows(
    tibble(relative_day = 0, success_smoothed = 0, n_eligible_sum = 0),
    s %>% select(relative_day, success_smoothed, n_eligible_sum)
  ) %>%
    arrange(relative_day) %>%
    distinct(relative_day, .keep_all = TRUE)

  scale_max <- max(s$success_smoothed, na.rm = TRUE)
  if (!is.finite(scale_max) || scale_max <= 0) stop("Q scaling maximum is not positive.")

  # Columns: relative_day, tau_q, q_value (= success / max success, the shared Q),
  # unsafe_funeral_comm_proxy (= 1 - success, absolute), success_smoothed, n_eligible_sum.
  s %>%
    mutate(
      tau_q = relative_day / max_day,
      q_value = clip01(success_smoothed / scale_max),
      unsafe_funeral_comm_proxy = clip01(1 - success_smoothed)
    ) %>%
    select(relative_day, tau_q, q_value, unsafe_funeral_comm_proxy,
           success_smoothed, n_eligible_sum)
}

# ---- 2d. Build the three DRC Q series --------------------------------------
drc_anchor_max_day <- max(drc_anchors$relative_day, na.rm = TRUE)

weekly_binned <- bin_provinces(sdb_linelist, "epi_week")
weekly_avg    <- average_provinces(weekly_binned, "prop_success")

drc_conflict_qseries <- finalise_q_series(
  weekly_avg, smooth = TRUE, drc_anchor_max_day = drc_anchor_max_day
)
drc_conflict_plusplus_qseries <- finalise_q_series(
  weekly_avg, smooth = TRUE, drc_anchor_max_day = drc_anchor_max_day,
  collapse_window = PLUSPLUS_WINDOW_DAY
)

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
drc_durations <- tibble(
  scenario = c("drc_conflict", "drc_conflict_plusplus", "drc_no_conflict", "wa_anchor_max_day"),
  max_day  = c(
    max(drc_conflict_qseries$relative_day),
    max(drc_conflict_plusplus_qseries$relative_day),
    max(drc_no_conflict_qseries$relative_day),
    max(wa_anchors$relative_day, na.rm = TRUE)
  )
)

wa_prep <- list(anchors = wa_anchors)
saveRDS(wa_prep, file.path(DIR_PROCESSED,
        "WestAfrica_QCurve_revisedMethod/WestAfrica_QCurve_PreppedData(revisedMethod).rds"))

drc_prep <- list(
  anchors                   = drc_anchors,
  conflict_qseries          = drc_conflict_qseries,
  conflict_plusplus_qseries = drc_conflict_plusplus_qseries,
  no_conflict_qseries       = drc_no_conflict_qseries,
  durations                 = drc_durations,
  province_weekly_qc        = weekly_binned
)
saveRDS(drc_prep, file.path(DIR_PROCESSED,
        "DRC_QCurve_revisedMethod/DRC_QCurve_PreppedData(revisedMethod).rds"))

message("\nDRC scenario horizons (max day):")
print(drc_durations)
message("\n00_DataPreparation_and_Cleaning(revisedMethod).R complete. ",
        "Anchors read from the parameter-table workbook; SDB Q reconstructed as in the original 00.")
