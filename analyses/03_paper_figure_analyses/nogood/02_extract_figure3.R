# =============================================================================
# 02_extract_figure3.R
# Extract and save run summaries for Figure 3.
# Output: output_figgen/figure_3_run_summary.csv
#
<<<<<<< HEAD
# Also extracts weekly averted HCW-deaths and HCW-days-lost incidence across
# all coverage x efficacy combinations, for cumulative averted-burden panels.
=======
# Also extracts weekly averted HCW deaths incidence (efficacy 80%) for the
# cumulative averted-deaths overlay panel.
>>>>>>> 873ecc12b709e76c5085cd6ebf2f57c289f1da8c
# Output: output_figgen/figure_3_averted_weekly_ts.csv
# =============================================================================
source(here::here("analyses", "03_paper_figure_analyses", "helper_functions_figure_1to4.R"))
FIG3_EFFICACY_LEVELS <- c("obv_50", "obv_60", "obv_70", "obv_80", "obv_90")

# -----------------------------------------------------------------------
# Run summary (deaths averted / days lost averted boxplots)
# -----------------------------------------------------------------------
message("Extracting run summaries for figure 3...")
run_df <- do.call(rbind, lapply(COVERAGE_LEVELS, function(cov) {
  do.call(rbind, lapply(FIG3_EFFICACY_LEVELS, function(eff_name) {
    arm_dir   <- sprintf("%s_obv%02d", cov, round(OBV_EFFICACY_VALUES[[eff_name]] * 100))
    arm_label <- sprintf("%s__%s", cov, eff_name)
<<<<<<< HEAD
    extract_run_summary(arm_dir, arm_label = arm_label, n_workers = 14L, obv_return = FALSE)
=======
    extract_run_summary(arm_dir, arm_label = arm_label, n_workers = 10L, obv_return = FALSE)
>>>>>>> 873ecc12b709e76c5085cd6ebf2f57c289f1da8c
  }))
}))
save_figure_data(run_df, "figure_3_run_summary.csv")

# -----------------------------------------------------------------------
<<<<<<< HEAD
# Weekly averted HCW deaths / days lost incidence
# across all coverage x efficacy combinations
# -----------------------------------------------------------------------
message("Extracting averted weekly ts for figure 3 (all coverage x efficacy)...")
ts_list <- do.call(rbind, lapply(COVERAGE_LEVELS, function(cov) {
  do.call(rbind, lapply(FIG3_EFFICACY_LEVELS, function(eff_name) {
    arm_dir <- sprintf("%s_obv%02d", cov, round(OBV_EFFICACY_VALUES[[eff_name]] * 100))
    df <- extract_weekly_ts(arm_dir, n_workers = 14L)
    df <- df[df$metric %in% c("averted_hcw_deaths_incidence",
                              "averted_hcw_days_lost_incidence"), ]
    df$coverage_name <- cov
    df$eff_name       <- eff_name
    df
  }))
}))
save_figure_data(ts_list, "figure_3_averted_weekly_ts.csv")

# -----------------------------------------------------------------------
# Weekly HCW deaths incidence at 80% efficacy, for each coverage scenario
# (baseline + full/ramp_high/ramp_low), for the epi-curve overlay panel
# -----------------------------------------------------------------------
EFF80 <- "obv_80"

message("Extracting weekly HCW deaths incidence at 80% efficacy...")
weekly_80_list <- lapply(COVERAGE_LEVELS, function(cov) {
  arm_dir <- sprintf("%s_obv%02d", cov, round(OBV_EFFICACY_VALUES[[EFF80]] * 100))
  df <- extract_weekly_ts(arm_dir, n_workers = 14L)
  df <- df[df$metric == "hcw_deaths_incidence", ]
  df$coverage_name <- cov
  df
})
weekly_80_all <- do.call(rbind, weekly_80_list)
save_figure_data(weekly_80_all, "figure_3_weekly_hcw_deaths_80.csv")
=======
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
>>>>>>> 873ecc12b709e76c5085cd6ebf2f57c289f1da8c

message("Figure 3 data extraction complete.")