# =============================================================================
# 02_extract_figure4.R
# Extract and save run summaries for Figure 4.
# Output: output_figgen/figure_4_run_summary.csv
# =============================================================================
source(here::here("analyses", "03_paper_figure_analyses", "helper_functions_figure_1to4.R"))

COVERAGE_GRID <- c(0.10, 0.30, 0.50, 0.70, 0.90)
EFFICACY_GRID <- c(0.50, 0.60, 0.70, 0.80, 0.90)

message("Extracting run summaries for figure 4...")

run_df <- do.call(rbind, lapply(COVERAGE_GRID, function(cov) {
  do.call(rbind, lapply(EFFICACY_GRID, function(eff) {
    arm_dir   <- sprintf("const%02d_obv%02d", round(cov * 100), round(eff * 100))
    arm_label <- sprintf("cov%02d_obv%02d", round(cov * 100), round(eff * 100))
    df <- extract_run_summary(arm_dir, arm_label = arm_label, n_workers = 10L, obv_return = FALSE)
    df$obv_efficacy <- eff
    df$obv_coverage <- cov
    df
  }))
}))

save_figure_data(run_df, "figure_4_run_summary.csv")
message("Figure 4 data extraction complete.")