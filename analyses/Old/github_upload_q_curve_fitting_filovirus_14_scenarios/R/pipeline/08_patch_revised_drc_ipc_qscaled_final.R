# ============================================================
# Patch revised DRC conflict / DRC conflict++ IPC/PPE so that
# it follows the revised-methodology q-scaled rule throughout:
#   ipc_helper(t) = ipc_low + (ipc_high - ipc_low) * q_value(t)
# with the DRC++ collapse window forced to the poor-response endpoint.
#
# This script does NOT rerun Stan, does NOT rebuild scenarios, and does
# NOT overwrite earlier final CSVs. It writes new *_ipcQscaled.csv files.
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(readxl)
  library(dplyr)
  library(stringr)
})

# ------------------------------------------------------------
# Locate folders
# ------------------------------------------------------------
base_dir <- getwd()
out_dir <- file.path(base_dir, "v6_run_outputs", "final_730day_outputs")

if (!dir.exists(out_dir)) {
  stop(
    "Cannot find final output folder. Run this from the filovirus_14scenario_v8 folder, e.g.\n",
    "setwd('C:/Users/jnstapley/Documents/Efficacy_curves/filovirus_14scenario_v8')"
  )
}

# Prefer the latest checked files as inputs; fall back sensibly.
choose_existing <- function(paths, label) {
  paths <- paths[file.exists(paths)]
  if (length(paths) == 0) stop("Could not find input file for: ", label)
  paths[[1]]
}

orig_in <- choose_existing(c(
  file.path(out_dir, "final_original_methodology_7_scenarios_730day_matrix_ipcLitMax.csv"),
  file.path(out_dir, "final_original_methodology_7_scenarios_730day_matrix_terminalIPC.csv"),
  file.path(out_dir, "final_original_methodology_7_scenarios_730day_matrix.csv")
), "original matrix")

rev_in <- choose_existing(c(
  file.path(out_dir, "final_revised_methodology_7_scenarios_730day_matrix_ipcLitMax.csv"),
  file.path(out_dir, "final_revised_methodology_7_scenarios_730day_matrix_terminalIPC.csv"),
  file.path(out_dir, "final_revised_methodology_7_scenarios_730day_matrix.csv")
), "revised matrix")

orig_out <- file.path(out_dir, "final_original_methodology_7_scenarios_730day_matrix_ipcQscaled.csv")
rev_out  <- file.path(out_dir, "final_revised_methodology_7_scenarios_730day_matrix_ipcQscaled.csv")

audit_path <- file.path(out_dir, "scenario_manifest_output_resolution_audit.csv")

message("Original input: ", orig_in)
message("Revised input:  ", rev_in)

orig <- readr::read_csv(orig_in, show_col_types = FALSE)
rev  <- readr::read_csv(rev_in,  show_col_types = FALSE)

# ------------------------------------------------------------
# Repair any missing q_value values from q summary files if needed
# ------------------------------------------------------------
repair_missing_q <- function(df, audit_path) {
  if (!"q_value" %in% names(df) || !file.exists(audit_path)) return(df)
  audit <- readr::read_csv(audit_path, show_col_types = FALSE)
  if (!all(c("methodology", "scenario_key", "q_path") %in% names(audit))) return(df)

  for (i in seq_len(nrow(audit))) {
    row <- audit[i, ]
    if (!all(c("methodology", "scenario_key") %in% names(df))) next

    idx <- df$methodology == row$methodology & df$scenario_key == row$scenario_key
    if (!any(idx)) next
    if (sum(is.na(df$q_value[idx])) == 0) next
    if (is.na(row$q_path) || !file.exists(row$q_path)) next

    qdat <- readr::read_csv(row$q_path, show_col_types = FALSE)
    value_col <- dplyr::case_when(
      "mean" %in% names(qdat) ~ "mean",
      "q_value" %in% names(qdat) ~ "q_value",
      "median" %in% names(qdat) ~ "median",
      TRUE ~ NA_character_
    )
    if (is.na(value_col) || !"relative_day" %in% names(qdat)) next

    qday <- qdat %>%
      dplyr::mutate(
        relative_day = as.numeric(relative_day),
        q_value = as.numeric(.data[[value_col]])
      ) %>%
      dplyr::filter(is.finite(relative_day), is.finite(q_value)) %>%
      dplyr::group_by(relative_day) %>%
      dplyr::summarise(q_value = mean(q_value, na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(relative_day)

    if (nrow(qday) < 1) next

    df$q_value[idx] <- stats::approx(
      x = qday$relative_day,
      y = qday$q_value,
      xout = df$relative_day[idx],
      rule = 2,
      ties = mean
    )$y
    df$q_value[idx] <- pmin(1, pmax(0, df$q_value[idx]))
    message("Repaired q_value for: ", row$methodology, " / ", row$scenario_key)
  }
  df
}

orig <- repair_missing_q(orig, audit_path)
rev  <- repair_missing_q(rev,  audit_path)

# ------------------------------------------------------------
# Read IPC endpoint range from workbook, with explicit fallback
# ------------------------------------------------------------
extract_range <- function(x) {
  x <- as.character(x)
  x <- stringr::str_replace_all(x, "[−—]", "–")
  nums <- stringr::str_extract_all(x, "[-+]?\\d*\\.?\\d+(?:[eE][-+]?\\d+)?")[[1]]
  if (length(nums) == 0) return(c(NA_real_, NA_real_))
  if (length(nums) == 1) return(rep(as.numeric(nums[1]), 2))
  as.numeric(nums[1:2])
}

find_parameter_workbook <- function() {
  candidates <- c(
    file.path(base_dir, "v6_run_outputs", "filovirus_model_parameter_table_4_scenarios_with_etu_baseline.xlsx"),
    file.path(base_dir, "filovirus_model_parameter_table_4_scenarios_with_etu_baseline.xlsx"),
    file.path(dirname(base_dir), "filovirus_model_parameter_table_4_scenarios_with_etu_baseline.xlsx"),
    file.path(getwd(), "filovirus_model_parameter_table_4_scenarios_with_etu_baseline.xlsx")
  )
  candidates[file.exists(candidates)][1]
}

parameter_workbook <- find_parameter_workbook()

ipc_low <- 0.071
ipc_high <- 0.746
source_range <- "fallback hard-coded 0.071-0.746"
source_workbook <- NA_character_

if (!is.na(parameter_workbook) && file.exists(parameter_workbook)) {
  raw <- tryCatch(
    readxl::read_excel(parameter_workbook, sheet = "DRC conflict-smoothed", col_names = FALSE),
    error = function(e) NULL
  )
  if (!is.null(raw) && ncol(raw) >= 5) {
    names(raw)[1:7] <- c("Section", "Parameter", "Symbol", "Description", "Value_range", "Reference", "URL")
    ipc_row <- raw %>%
      dplyr::filter(.data$Section == "Time-varying response") %>%
      dplyr::filter(stringr::str_detect(as.character(.data$Parameter), "ipc_helper|ppe_efficacy|IPC|PPE")) %>%
      dplyr::slice(1)
    if (nrow(ipc_row) == 1) {
      rng <- extract_range(ipc_row$Value_range[[1]])
      if (all(is.finite(rng))) {
        ipc_low <- min(rng)
        ipc_high <- max(rng)
        source_range <- as.character(ipc_row$Value_range[[1]])
        source_workbook <- parameter_workbook
      }
    }
  }
}

message("IPC endpoint range used: ", ipc_low, " to ", ipc_high)
message("Source range: ", source_range)

clip01 <- function(x) pmin(1, pmax(0, x))
ipc_from_q <- function(q) clip01(ipc_low + (ipc_high - ipc_low) * q)

# ------------------------------------------------------------
# Apply q-scaled IPC rule only to revised DRC conflict scenarios
# ------------------------------------------------------------
required_cols <- c("scenario_key", "relative_day", "q_value", "ipc_helper")
missing_cols <- setdiff(required_cols, names(rev))
if (length(missing_cols) > 0) stop("Revised matrix is missing columns: ", paste(missing_cols, collapse = ", "))

rev_qscaled <- rev %>%
  dplyr::mutate(
    ipc_helper = dplyr::case_when(
      scenario_key == "middle_drc_conflict" ~
        ipc_from_q(q_value),

      scenario_key == "middle_drc_conflict_plusplus" & relative_day >= 200 & relative_day <= 300 ~
        ipc_low,

      scenario_key == "middle_drc_conflict_plusplus" ~
        ipc_from_q(q_value),

      TRUE ~ ipc_helper
    )
  )

# Original methodology is deliberately unchanged except for any q_value repair.
readr::write_csv(orig, orig_out)
readr::write_csv(rev_qscaled, rev_out)

message("Wrote unchanged original paired matrix: ", orig_out)
message("Wrote revised q-scaled IPC matrix: ", rev_out)

# ------------------------------------------------------------
# Checks
# ------------------------------------------------------------
check_tbl <- dplyr::bind_rows(orig, rev_qscaled) %>%
  dplyr::group_by(methodology, scenario_key) %>%
  dplyr::summarise(
    min_day = min(relative_day, na.rm = TRUE),
    max_day = max(relative_day, na.rm = TRUE),
    min_q = min(q_value, na.rm = TRUE),
    max_q = max(q_value, na.rm = TRUE),
    n_missing_q = sum(is.na(q_value)),
    .groups = "drop"
  )

print(check_tbl, n = Inf)

ipc_check <- rev_qscaled %>%
  dplyr::filter(scenario_key %in% c("middle_drc_conflict", "middle_drc_conflict_plusplus")) %>%
  dplyr::group_by(scenario_key) %>%
  dplyr::summarise(
    min_ipc = min(ipc_helper, na.rm = TRUE),
    max_ipc = max(ipc_helper, na.rm = TRUE),
    final_ipc = ipc_helper[which.max(relative_day)],
    .groups = "drop"
  )

print(ipc_check)

message("Done. Use *_ipcQscaled.csv as the final files for the consistent revised methodology.")
