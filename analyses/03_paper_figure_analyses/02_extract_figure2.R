# =============================================================================
# 02_extract_figure2.R
# Extract and save run summaries for Figure 2.
# Run this once after simulations are complete.
# Output: output_figgen/figure_2_run_summary.csv
# =============================================================================
source(here::here("analyses", "03_paper_figure_analyses", "helper_functions_figure_1to4.R"))

FIG2_EFFICACY_LEVELS <- c("obv_50", "obv_60", "obv_70", "obv_80", "obv_90")

message("Extracting run summaries for figure 2...")

run_df <- do.call(rbind, lapply(FIG2_EFFICACY_LEVELS, function(eff_name) {
  arm_dir <- sprintf("full_obv%02d", round(OBV_EFFICACY_VALUES[[eff_name]] * 100))
  extract_run_summary(arm_dir, arm_label = eff_name, n_workers = 10L)
}))

save_figure_data(run_df, "figure_2_run_summary.csv")
message("Figure 2 data extraction complete.")