# =============================================================================
# 02_extract_figure3.R
# Extract and save run summaries for Figure 3.
# Output: output_figgen/figure_3_run_summary.csv
#
# Also extracts weekly incident HCW deaths at 80% antiviral efficacy, for
# baseline (no antiviral) and each coverage scenario (full/ramp_high/ramp_low),
# pre-aggregated to quantiles across particle x rep (same convention as
# figure 1's weekly ts CSV).
# Output: output_figgen/figure_3_weekly_hcw_deaths_80.csv
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
    extract_run_summary(arm_dir, arm_label = arm_label, n_workers = 14L, obv_return = FALSE)
  }))
}))
save_figure_data(run_df, "figure_3_run_summary.csv")

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
      (arm == "baseline" & coverage_name == "full")
  ) %>%
  mutate(line_group = ifelse(arm == "baseline", "baseline", coverage_name))

# Aggregate to quantiles across particle x rep
# weekly_80_q <- weekly_80_clean %>%
#   mutate(week = week / 7) %>%
#   group_by(scenario, line_group, week) %>%
#   summarise(
#     q025 = quantile(value, 0.025),
#     q25  = quantile(value, 0.25),
#     q50  = quantile(value, 0.50),
#     q75  = quantile(value, 0.75),
#     q975 = quantile(value, 0.975),
#     .groups = "drop"
#   )
weekly_80_q <- weekly_80_clean %>%
  mutate(week = week / 7) %>%
  group_by(scenario, line_group, week) %>%
  summarise(
    mean_val = mean(value),
    sd_val   = sd(value),
    n        = n(),
    .groups = "drop"
  ) %>%
  mutate(
    se_val = sd_val / sqrt(n),
    q025   = mean_val - 1.96 * se_val,
    q975   = mean_val + 1.96 * se_val,
    q25    = mean_val - se_val,
    q75    = mean_val + se_val,
    q50    = mean_val
  )

save_figure_data(weekly_80_q, "figure_3_weekly_hcw_deaths_80.csv")

message("Figure 3 data extraction complete.")