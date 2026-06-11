# =============================================================================
# 02_extract_figure3.R
# Extract and save run summaries for Figure 3.
# Output: output_figgen/figure_3_run_summary.csv
# =============================================================================
source(here::here("analyses", "03_figure_template", "helper_functions_figure_1to4.R"))

FIG3_EFFICACY_LEVELS <- c("obv_50", "obv_60", "obv_70", "obv_80", "obv_90")

message("Extracting run summaries for figure 3...")

run_df <- do.call(rbind, lapply(COVERAGE_LEVELS, function(cov) {
  do.call(rbind, lapply(FIG3_EFFICACY_LEVELS, function(eff_name) {
    arm_dir   <- sprintf("%s_obv%02d", cov, round(OBV_EFFICACY_VALUES[[eff_name]] * 100))
    arm_label <- sprintf("%s__%s", cov, eff_name)
    extract_run_summary(arm_dir, arm_label = arm_label, n_workers = 10L)
  }))
}))

save_figure_data(run_df, "figure_3_run_summary.csv")
message("Figure 3 data extraction complete.")