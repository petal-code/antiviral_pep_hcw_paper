
# ============================================================
# Run all 14 scripts and print local plots without saving by default
# ============================================================
V6_BUNDLE_ROOT <- normalizePath(dirname(sys.frame(1)$ofile), winslash = "/", mustWork = TRUE)
source(file.path(V6_BUNDLE_ROOT, "_00_v6_helpers.R"))

options(scenario.save_plots = getOption("scenario.save_plots", FALSE))
v6_disable_saving_plots()

run_dir <- file.path(V6_BUNDLE_ROOT, "v6_run_outputs")
v6_prepare_run_dir(run_dir)
oldwd <- getwd()
setwd(run_dir)
on.exit(setwd(oldwd), add = TRUE)

message("Bundle root: ", V6_BUNDLE_ROOT)
message("Run/output directory: ", run_dir)
message("Plot saving: ", ifelse(isTRUE(getOption("scenario.save_plots")), "ON", "OFF: plots will print locally only"))

# Scenario run plan.  Each script is followed by a capture step so that, even
# if the native script writes to an unexpected directory or leaves only the
# in-memory object, the expected CSV is materialised in v6_run_outputs.
orig_dir <- file.path(V6_BUNDLE_ROOT, "scripts", "original_methodology")
rev_dir  <- file.path(V6_BUNDLE_ROOT, "scripts", "revised_methodology")

run_plan <- tibble::tribble(
  ~methodology, ~script_dir, ~script_file, ~matrix_file, ~curve_file, ~q_file, ~anchor_file, ~sdb_file,
  "original", orig_dir, "West_Africa_USE_partial_pooling_independent_bounds_UFaCD_tweaked.R", "west_africa_partial_pool_normalisedQ_estimated_bounds_tweaked_ufc_bp_input_matrix.csv", "west_africa_partial_pool_normalisedQ_estimated_bounds_tweaked_ufc_curve_summaries.csv", "west_africa_partial_pool_normalisedQ_estimated_bounds_tweaked_ufc_Q_summaries.csv", "west_africa_partial_pool_normalisedQ_estimated_bounds_tweaked_ufc_anchor_rows_used.csv", NA_character_,
  "original", orig_dir, "best_composite_partial_pool_normalisedQ_estimated_bounds_tuned.R", "best_composite_partial_pool_normalisedQ_estimated_bounds_tuned_bp_input_matrix.csv", "best_composite_partial_pool_normalisedQ_estimated_bounds_tuned_curve_summaries.csv", "best_composite_partial_pool_normalisedQ_estimated_bounds_tuned_Q_summaries.csv", "best_composite_partial_pool_normalisedQ_estimated_bounds_tuned_anchor_rows_used.csv", NA_character_,
  "original", orig_dir, "drc_conflict_warsame_weekly_Q_relative_absolute_UFC_rolling4_sparse_greydots.R", "drc_conflict_warsame_weekly_Q_relative_absolute_UFC_rolling4_sparse_greydots_bp_input_matrix.csv", "drc_conflict_warsame_weekly_Q_relative_absolute_UFC_rolling4_sparse_greydots_curve_summaries.csv", "drc_conflict_warsame_weekly_Q_relative_absolute_UFC_rolling4_sparse_greydots_Q_summaries.csv", "drc_conflict_warsame_weekly_Q_relative_absolute_UFC_rolling4_sparse_greydots_anchor_rows_used.csv", "drc_conflict_warsame_weekly_Q_relative_absolute_UFC_rolling4_sparse_greydots_averaged_conflict_binned_sdb_success.csv",
  "original", orig_dir, "drc_non_conflict_warsame_only_Q_minus_conflict_time_original_plateau_selectfix.R", "drc_non_conflict_warsame_only_Q_minus_conflict_time_original_plateau_bp_input_matrix.csv", "drc_non_conflict_warsame_only_Q_minus_conflict_time_original_plateau_curve_summaries.csv", "drc_non_conflict_warsame_only_Q_minus_conflict_time_original_plateau_Q_summaries.csv", "drc_non_conflict_warsame_only_Q_minus_conflict_time_original_plateau_anchor_rows_used.csv", "drc_non_conflict_warsame_only_Q_minus_conflict_time_original_plateau_averaged_no_conflict_binned_sdb_success.csv",
  "original", orig_dir, "drc_conflict_plusplus_original_method_directCollapse.R", "drc_conflict_plusplus_original_method_directCollapse_bp_input_matrix.csv", "drc_conflict_plusplus_original_method_directCollapse_curve_summaries.csv", "drc_conflict_plusplus_original_method_directCollapse_Q_summaries.csv", "drc_conflict_plusplus_original_method_directCollapse_anchor_rows_used.csv", "drc_conflict_plusplus_original_method_directCollapse_averaged_conflict_binned_sdb_success.csv",
  "original", orig_dir, "west_africa_with_conflict_extended_duration_update_matrix_FIX3.R", "west_africa_with_conflict_extended_DRC_Q_modulated_bp_input_matrix.csv", "west_africa_with_conflict_extended_DRC_Q_modulated_curve_summaries.csv", "west_africa_with_conflict_extended_DRC_Q_modulated_Q_summaries.csv", NA_character_, NA_character_,
  "original", orig_dir, "west_africa_conflict_plusplus_original_method_from_conflict.R", "west_africa_conflict_plusplus_original_method_directCollapse_bp_input_matrix.csv", "west_africa_conflict_plusplus_original_method_directCollapse_curve_summaries.csv", "west_africa_conflict_plusplus_original_method_directCollapse_Q_summaries.csv", NA_character_, NA_character_,
  "revised", rev_dir, "01_refit_worst_west_africa_endpoint_constrained_WINDOW_EXTREMA_ZERO_PLATEAU_LABELS.R", "worst_west_africa_endpoint_constrained_zero_plateau_matrix.csv", "worst_west_africa_endpoint_constrained_zero_plateau_curve_summaries.csv", "worst_west_africa_endpoint_constrained_zero_plateau_Q_summaries.csv", "worst_west_africa_endpoint_constrained_zero_plateau_anchor_table_used.csv", NA_character_,
  "revised", rev_dir, "04_refit_best_east_africa_endpoint_constrained_WINDOW_EXTREMA_ZERO_PLATEAU_LABELS.R", "best_east_africa_endpoint_constrained_zero_plateau_matrix.csv", "best_east_africa_endpoint_constrained_zero_plateau_curve_summaries.csv", "best_east_africa_endpoint_constrained_zero_plateau_Q_summaries.csv", "best_east_africa_endpoint_constrained_zero_plateau_anchor_table_used.csv", NA_character_,
  "revised", rev_dir, "02_refit_drc_no_conflict_endpoint_constrained_SHORT_HORIZON_FILTERED.R", "drc_no_conflict_endpoint_constrained_zero_plateau_matrix.csv", "drc_no_conflict_endpoint_constrained_zero_plateau_curve_summaries.csv", "drc_no_conflict_endpoint_constrained_zero_plateau_Q_summaries.csv", "drc_no_conflict_endpoint_constrained_zero_plateau_anchor_table_used.csv", NA_character_,
  "revised", rev_dir, "03_rebuild_drc_conflict_endpoint_constrained_PRESERVE_Q_DIRECT_UFC_GREY_POINTS_FIRST_GREY_START.R", "drc_conflict_endpoint_constrained_preserveQ_directUFC_firstGreyStart_matrix.csv", "drc_conflict_endpoint_constrained_preserveQ_directUFC_firstGreyStart_matrix_with_uncertainty_long.csv", NA_character_, "drc_conflict_endpoint_constrained_preserveQ_directUFC_firstGreyStart_anchor_table_used.csv", "drc_conflict_endpoint_constrained_preserveQ_directUFC_firstGreyStart_warsame_sdb_points_used.csv",
  "revised", rev_dir, "03_DRC_conflict_PLUSPLUS_direct_Warsame_UFC_ALLPARAMS_COLLAPSE_PATCHED_FIRST_GREY_START_IPC275_PROBAXIS01.R", "drc_conflict_plusplus_allparams_plateau_directWarsameUFC_firstGreyStart_IPC275_probAxis01_matrix.csv", "drc_conflict_plusplus_allparams_plateau_directWarsameUFC_firstGreyStart_IPC275_probAxis01_matrix_with_uncertainty_long.csv", NA_character_, "drc_conflict_plusplus_allparams_plateau_directWarsameUFC_firstGreyStart_IPC275_probAxis01_anchor_table_used.csv", "drc_conflict_plusplus_allparams_plateau_directWarsameUFC_firstGreyStart_IPC275_probAxis01_warsame_sdb_points_used.csv",
  "revised", rev_dir, "05_construct_west_africa_conflict_ORIGINAL_HYBRID_Q_latestWA_completeCI.R", "west_africa_with_conflict_originalHybridQ_latestWA_completeCI_bp_input_matrix.csv", "west_africa_with_conflict_originalHybridQ_latestWA_completeCI_curve_summaries.csv", "west_africa_with_conflict_originalHybridQ_latestWA_completeCI_Q_summaries.csv", NA_character_, NA_character_,
  "revised", rev_dir, "06_construct_west_africa_conflict_PLUSPLUS_originalHybridQ_latestWA_completeCI.R", "west_africa_with_conflict_plusplus_originalHybridQ_latestWA_completeCI_bp_input_matrix.csv", "west_africa_with_conflict_plusplus_originalHybridQ_latestWA_completeCI_curve_summaries.csv", "west_africa_with_conflict_plusplus_originalHybridQ_latestWA_completeCI_Q_summaries.csv", NA_character_, NA_character_
)

for (i in seq_len(nrow(run_plan))) {
  row <- run_plan[i,]
  script_path <- file.path(row$script_dir, row$script_file)
  v6_run_script(script_path, paste0(row$methodology, " / ", row$script_file))
  v6_capture_expected_outputs(
    matrix_file = row$matrix_file,
    curve_file = row$curve_file,
    q_file = row$q_file,
    anchor_file = row$anchor_file,
    sdb_file = row$sdb_file,
    label = paste0(row$methodology, " / ", row$script_file)
  )
}

source(file.path(V6_BUNDLE_ROOT, "01_build_final_730day_csvs.R"), local = .GlobalEnv)

message("\nRun complete. Final 730-day CSVs are in: ", file.path(run_dir, "final_730day_outputs"))
message("To reprint all 14 extended audit plots without rerunning fits, source 02_print_all_14_plots_local_only.R")
