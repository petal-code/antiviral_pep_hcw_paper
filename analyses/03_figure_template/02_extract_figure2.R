# =============================================================================
# 02_extract_figure2.R
# Extract and save run summaries for Figure 2.
# Output: output_figgen/figure_2_run_summary.csv
#
# Also extracts weekly incident HCW deaths at 80% antiviral efficacy, for
# baseline (no antiviral) and each coverage scenario (full/ramp_high/ramp_low).
#
# Weekly incidence is now summarised using the two-step particle-quantile
# method (matching Figure 1 and Figure 3new):
#   Step 1: mean over reps within each particle
#   Step 2: quantiles across particles (q025/q25/q50/q75/q975)
# The old mean +/- 1.96*SE approach pooled particle x rep directly,
# producing near-invisible ribbons due to overly small SE.
#
# Output: output_figgen/figure_2_weekly_hcw_deaths_80.csv
# =============================================================================
source(here::here("analyses", "03_figure_template", "helper_functions_figure_1to4.R"))
FIG2_EFFICACY_LEVELS <- c("obv_50", "obv_60", "obv_70", "obv_80", "obv_90")
# -----------------------------------------------------------------------
# Run summary (deaths averted / days lost averted boxplots)
# -----------------------------------------------------------------------
message("Extracting run summaries for figure 2...")
run_df <- do.call(rbind, lapply(COVERAGE_LEVELS, function(cov) {
  do.call(rbind, lapply(FIG2_EFFICACY_LEVELS, function(eff_name) {
    arm_dir   <- sprintf("%s_obv%02d", cov, round(OBV_EFFICACY_VALUES[[eff_name]] * 100))
    arm_label <- sprintf("%s__%s", cov, eff_name)
    extract_run_summary(arm_dir, arm_label = arm_label, n_workers = 14L, obv_return = FALSE)
  }))
}))
save_figure_data(run_df, "figure_2_run_summary.csv")
# -----------------------------------------------------------------------
# Weekly incident HCW deaths at 80% efficacy
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
weekly_80_raw <- do.call(rbind, weekly_80_list)
# Baseline (no antiviral) is identical regardless of which arm it came from,
# since OBV doesn't affect transmission dynamics. Keep baseline only from the
# "full" arm and label it separately from the OBV arms.
weekly_80_clean <- weekly_80_raw %>%
  filter(
    (arm == "obv") |
      (arm == "baseline" & coverage_name == "ramp_low")
  ) %>%
  mutate(line_group = ifelse(arm == "baseline", "baseline", coverage_name))
# Two-step aggregation matching Figure 1 / Figure 3new convention:
#   Step 1: mean over reps within each particle (reduces 200x10 to 200)
#   Step 2: quantiles across the 200 particles
# Incident values do NOT receive cumsum (unlike hcw_deaths cumulative metric).
weekly_80_q <- weekly_80_clean %>%
  mutate(week = week / 7) %>%
  # Step 1: rep average within each particle
  group_by(scenario, line_group, particle_id, week) %>%
  summarise(value = mean(value), .groups = "drop") %>%
  # Step 2: quantiles across particles
  group_by(scenario, line_group, week) %>%
  summarise(
    q025 = quantile(value, 0.025, na.rm = TRUE),
    q25  = quantile(value, 0.25,  na.rm = TRUE),
    q50  = quantile(value, 0.50,  na.rm = TRUE),
    q75  = quantile(value, 0.75,  na.rm = TRUE),
    q975 = quantile(value, 0.975, na.rm = TRUE),
    .groups = "drop"
  )
save_figure_data(weekly_80_q, "figure_2_weekly_hcw_deaths_80.csv")
message("Figure 2 data extraction complete.")