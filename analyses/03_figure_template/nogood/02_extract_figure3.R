# =============================================================================
# 02_extract_figure3.R
# Extract and save run summaries for Figure 3.
# Output: output_figgen/figure_3_run_summary.csv
#
# Also extracts weekly averted HCW deaths incidence (efficacy 80%) for the
# cumulative averted-deaths overlay panel.
# Output: output_figgen/figure_3_averted_weekly_ts.csv
# =============================================================================
source(here::here("analyses", "03_figure_template", "helper_functions_figure_1to4.R"))
FIG3_EFFICACY_LEVELS <- c("obv_50", "obv_60", "obv_70", "obv_80", "obv_90")

# -----------------------------------------------------------------------
# Run summary (deaths averted / days lost averted boxplots)
# -----------------------------------------------------------------------
message("Extracting run summaries for figure 3...")
run_df <- do.call(rbind, lapply(COVERAGE_LEVELS, function(cov) {
  do.call(rbind, lapply(FIG3_EFFICACY_LEVELS, function(eff_name) {
    arm_dir   <- sprintf("%s_obv%02d", cov, round(OBV_EFFICACY_VALUES[[eff_name]] * 100))
    arm_label <- sprintf("%s__%s", cov, eff_name)
    extract_run_summary(arm_dir, arm_label = arm_label, n_workers = 10L, obv_return = FALSE)
  }))
}))
save_figure_data(run_df, "figure_3_run_summary.csv")

# -----------------------------------------------------------------------
# Weekly averted HCW deaths incidence (for cumulative averted-deaths panel)
# at efficacy 80% only
# -----------------------------------------------------------------------
EFF <- "obv_80"

message("Extracting averted weekly ts for figure 3 (efficacy 80%)...")
ts_list <- lapply(COVERAGE_LEVELS, function(cov) {
  arm_dir <- sprintf("%s_obv%02d", cov, round(OBV_EFFICACY_VALUES[[EFF]] * 100))
  df <- extract_weekly_ts(arm_dir, n_workers = 10L)
  df <- df[df$metric == "averted_hcw_deaths_incidence", ]
  df$coverage_name <- cov
  df
})
ts_all <- do.call(rbind, ts_list)
save_figure_data(ts_all, "figure_3_averted_weekly_ts.csv")

message("Figure 3 data extraction complete.")