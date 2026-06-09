# ============================================================
# 12_patch_WA_conflict_plusplus_direction_final.R
#
# Fix West Africa conflict++ directionality in the final 730-day matrices.
#
# Problem fixed:
#   The West Africa conflict++ collapse window must represent response
#   deterioration, not improvement.
#
# Logic:
#   For each methodology, construct worst_west_africa_conflict_plusplus from
#   the existing worst_west_africa_conflict trajectory.
#
#   Outside the ++ collapse window: keep West Africa with conflict.
#   Inside the ++ collapse window: force parameters to poor-response endpoints.
#
#   Increasing/good-response parameters:
#     prob_hosp, prop_etu, ipc_helper -> poor endpoint = minimum value
#
#   Adverse parameters:
#     delay_hosp, prob_unsafe_funeral_comm, prob_unsafe_funeral_hosp
#     -> poor endpoint = maximum value
#
#   q_value is set to 0 during the collapse window.
#
# This script does not rerun Stan, rebuild all scenarios, or overwrite the
# current final files. It writes *_ipcQscaled_WAplusplusFixed.csv.
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

out_dir <- "C:/Users/jnstapley/Documents/Efficacy_curves/filovirus_14scenario_v8/v6_run_outputs/final_730day_outputs"

orig_in <- file.path(out_dir, "final_original_methodology_7_scenarios_730day_matrix_ipcQscaled.csv")
rev_in  <- file.path(out_dir, "final_revised_methodology_7_scenarios_730day_matrix_ipcQscaled.csv")

orig_out <- file.path(out_dir, "final_original_methodology_7_scenarios_730day_matrix_ipcQscaled_WAplusplusFixed.csv")
rev_out  <- file.path(out_dir, "final_revised_methodology_7_scenarios_730day_matrix_ipcQscaled_WAplusplusFixed.csv")

if (!file.exists(orig_in)) stop("Cannot find original input: ", orig_in)
if (!file.exists(rev_in))  stop("Cannot find revised input: ", rev_in)

plusplus_window <- c(200, 300)

good_response_params <- c(
  "prob_hosp",
  "prop_etu",
  "ipc_helper"
)

adverse_params <- c(
  "delay_hosp",
  "prob_unsafe_funeral_comm",
  "prob_unsafe_funeral_hosp"
)

all_patch_params <- c(good_response_params, adverse_params)

poor_endpoint_values <- function(base_df) {
  # Derive poor-response endpoints from the baseline WA-with-conflict
  # trajectory, using extrema over the whole trajectory.
  #
  # This avoids relying on possibly incorrect existing ++ values.
  out <- list()

  for (p in good_response_params) {
    if (p %in% names(base_df)) {
      out[[p]] <- min(base_df[[p]], na.rm = TRUE)
    }
  }

  for (p in adverse_params) {
    if (p %in% names(base_df)) {
      out[[p]] <- max(base_df[[p]], na.rm = TRUE)
    }
  }

  out
}

patch_one_methodology <- function(df, methodology_label) {
  required_scenarios <- c(
    "worst_west_africa_conflict",
    "worst_west_africa_conflict_plusplus"
  )

  missing_scenarios <- setdiff(required_scenarios, unique(df$scenario_key))
  if (length(missing_scenarios) > 0) {
    stop(
      "Missing required scenario(s) in ", methodology_label, ": ",
      paste(missing_scenarios, collapse = ", ")
    )
  }

  base <- df %>%
    filter(scenario_key == "worst_west_africa_conflict") %>%
    arrange(relative_day)

  target_template <- df %>%
    filter(scenario_key == "worst_west_africa_conflict_plusplus") %>%
    arrange(relative_day)

  if (nrow(base) != nrow(target_template)) {
    stop("Base and target row counts differ for ", methodology_label)
  }

  if (!all(base$relative_day == target_template$relative_day)) {
    stop("Base and target relative_day grids differ for ", methodology_label)
  }

  poor <- poor_endpoint_values(base)

  patched_target <- base

  # Preserve target scenario identifiers/names from existing ++ scenario where possible.
  patched_target$scenario_key <- target_template$scenario_key

  if ("scenario" %in% names(patched_target) && "scenario" %in% names(target_template)) {
    patched_target$scenario <- target_template$scenario
  }

  if ("methodology" %in% names(patched_target)) {
    patched_target$methodology <- methodology_label
  }

  collapse_idx <- patched_target$relative_day >= plusplus_window[1] &
    patched_target$relative_day <= plusplus_window[2]

  for (p in names(poor)) {
    patched_target[[p]][collapse_idx] <- poor[[p]]
  }

  if ("q_value" %in% names(patched_target)) {
    patched_target$q_value[collapse_idx] <- 0
  }

  # Replace target scenario in full data frame.
  out <- df %>%
    filter(scenario_key != "worst_west_africa_conflict_plusplus") %>%
    bind_rows(patched_target) %>%
    arrange(methodology, scenario_key, relative_day)

  # Diagnostics for this scenario.
  diag <- patched_target %>%
    mutate(in_collapse = collapse_idx) %>%
    group_by(in_collapse) %>%
    summarise(
      min_day = min(relative_day, na.rm = TRUE),
      max_day = max(relative_day, na.rm = TRUE),
      min_prob_hosp = min(prob_hosp, na.rm = TRUE),
      max_prob_hosp = max(prob_hosp, na.rm = TRUE),
      min_delay_hosp = min(delay_hosp, na.rm = TRUE),
      max_delay_hosp = max(delay_hosp, na.rm = TRUE),
      min_ufc = min(prob_unsafe_funeral_comm, na.rm = TRUE),
      max_ufc = max(prob_unsafe_funeral_comm, na.rm = TRUE),
      min_ufh = min(prob_unsafe_funeral_hosp, na.rm = TRUE),
      max_ufh = max(prob_unsafe_funeral_hosp, na.rm = TRUE),
      min_etu = min(prop_etu, na.rm = TRUE),
      max_etu = max(prop_etu, na.rm = TRUE),
      min_ipc = min(ipc_helper, na.rm = TRUE),
      max_ipc = max(ipc_helper, na.rm = TRUE),
      min_q = min(q_value, na.rm = TRUE),
      max_q = max(q_value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(methodology = methodology_label, .before = 1)

  attr(out, "wa_plusplus_diagnostic") <- diag
  attr(out, "wa_plusplus_poor_endpoints") <- tibble(
    methodology = methodology_label,
    parameter = names(poor),
    poor_endpoint = unlist(poor)
  )

  out
}

orig <- readr::read_csv(orig_in, show_col_types = FALSE)
rev  <- readr::read_csv(rev_in, show_col_types = FALSE)

orig_fixed <- patch_one_methodology(orig, "original")
rev_fixed  <- patch_one_methodology(rev, "revised")

readr::write_csv(orig_fixed, orig_out)
readr::write_csv(rev_fixed, rev_out)

message("Wrote WA conflict++ direction-fixed original matrix:")
message(orig_out)
message("Wrote WA conflict++ direction-fixed revised matrix:")
message(rev_out)

diagnostics <- bind_rows(
  attr(orig_fixed, "wa_plusplus_diagnostic"),
  attr(rev_fixed, "wa_plusplus_diagnostic")
)

poor_endpoints <- bind_rows(
  attr(orig_fixed, "wa_plusplus_poor_endpoints"),
  attr(rev_fixed, "wa_plusplus_poor_endpoints")
)

message("\nPoor-response endpoints used during WA conflict++ collapse:")
print(poor_endpoints)

message("\nWA conflict++ collapse diagnostics:")
print(diagnostics)

# Whole-file sanity checks.
combined_check <- bind_rows(orig_fixed, rev_fixed) %>%
  group_by(methodology, scenario_key) %>%
  summarise(
    n_rows = n(),
    n_days = n_distinct(relative_day),
    min_day = min(relative_day, na.rm = TRUE),
    max_day = max(relative_day, na.rm = TRUE),
    min_q = min(q_value, na.rm = TRUE),
    max_q = max(q_value, na.rm = TRUE),
    n_missing_q = sum(is.na(q_value)),
    .groups = "drop"
  )

message("\nFull final-file sanity check:")
print(combined_check)

bad <- combined_check %>%
  filter(n_rows != 731 | n_days != 731 | min_day != 0 | max_day != 730 |
           n_missing_q != 0 | min_q < -1e-8 | max_q > 1 + 1e-8)

if (nrow(bad) > 0) {
  print(bad)
  stop("Sanity check failed.")
}

message("\nDone. Use *_ipcQscaled_WAplusplusFixed.csv as the final matrices if the plots now look correct.")
