# =============================================================================
# 02_extract_figure_supp_A.R
# Extract and save run summaries for Figure Supp A.
# Run this once after simulations are complete.
# Output: output_figgen/figure_supp_A_run_summary.csv
# =============================================================================
source(here::here("analyses", "03_figure_template", "helper_functions_figure_1to4.R"))
FIGSUPPA_EFFICACY_LEVELS <- c("obv_50", "obv_60", "obv_70", "obv_80", "obv_90")
message("Extracting run summaries for figure supp A...")
run_df <- do.call(rbind, lapply(FIGSUPPA_EFFICACY_LEVELS, function(eff_name) {
  arm_dir <- sprintf("full_obv%02d", round(OBV_EFFICACY_VALUES[[eff_name]] * 100))
  extract_run_summary(arm_dir, arm_label = eff_name, n_workers = 10L)
}))
save_figure_data(run_df, "figure_supp_A_run_summary.csv")
message("Figure Supp A data extraction complete.")