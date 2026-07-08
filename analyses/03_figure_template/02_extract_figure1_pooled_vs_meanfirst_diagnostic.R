# =============================================================================
# 01_extract_figure1_pooled_test.R
#
# PURPOSE
#   Reviewer comment: averaging the stochastic replicates before taking
#   quantiles across posterior draws compresses variance, because it collapses
#   replicate-level (stochastic) uncertainty before the quantile step. Only
#   parameter (posterior-draw) uncertainty survives.
#
#   This script computes BOTH versions on the same raw extraction and compares
#   interval widths directly, so we can see how much the current approach
#   (mean-first) shrinks the credible interval relative to full pooling.
#
#   Version A "mean_first"  (current 01_extract_figure1.R approach):
#     1. average the N_REPS replicates within each particle_id  -> 200 values
#     2. quantiles across those 200 particle-level means
#
#   Version B "pooled"  (reviewer-suggested approach):
#     1. do NOT average reps away
#     2. quantiles directly across all (particle_id, rep) pairs
#        (200 particles x 10 reps = 2000 values per week/metric)
#
# Both use the exact same raw_ts extraction so the only difference is the
# aggregation step. Nothing here is written back into the production
# figure_1_weekly_ts.csv -- this is a diagnostic comparison only.
# =============================================================================
source(here::here("analyses", "03_figure_template", "helper_functions_figure_1to4.R"))

raw_ts <- extract_weekly_ts("cov80_obv80", bin_width = 7, n_workers = 10L)

# -----------------------------------------------------------------------------
# Version A: mean-first (current approach)
# -----------------------------------------------------------------------------
ts_mean_first <- raw_ts %>%
  # Step 1: mean over reps per particle
  group_by(scenario, arm, particle_id, week, metric) %>%
  summarise(value = mean(value), .groups = "drop") %>%
  # Step 2: cumsum for hcw_deaths per particle (after rep averaging)
  arrange(scenario, arm, particle_id, metric, week) %>%
  group_by(scenario, arm, particle_id, metric) %>%
  mutate(value = if (unique(metric) == "hcw_deaths") cumsum(value) else value) %>%
  ungroup() %>%
  # Step 3: quantiles over particles (n = 200 per week/metric)
  group_by(scenario, arm, week, metric) %>%
  summarise(
    n    = n(),
    q025 = quantile(value, 0.025),
    q25  = quantile(value, 0.25),
    q50  = quantile(value, 0.50),
    q75  = quantile(value, 0.75),
    q975 = quantile(value, 0.975),
    .groups = "drop"
  ) %>%
  mutate(week = week / 7, method = "mean_first")

# -----------------------------------------------------------------------------
# Version B: pooled (no rep-averaging before quantiles)
# -----------------------------------------------------------------------------
ts_pooled <- raw_ts %>%
  # Step 1: cumsum for hcw_deaths per (particle, rep) -- reps kept separate
  arrange(scenario, arm, particle_id, rep, metric, week) %>%
  group_by(scenario, arm, particle_id, rep, metric) %>%
  mutate(value = if (unique(metric) == "hcw_deaths") cumsum(value) else value) %>%
  ungroup() %>%
  # Step 2: quantiles over ALL (particle, rep) pairs (n = 200 x N_REPS per week/metric)
  group_by(scenario, arm, week, metric) %>%
  summarise(
    n    = n(),
    q025 = quantile(value, 0.025),
    q25  = quantile(value, 0.25),
    q50  = quantile(value, 0.50),
    q75  = quantile(value, 0.75),
    q975 = quantile(value, 0.975),
    .groups = "drop"
  ) %>%
  mutate(week = week / 7, method = "pooled")

# -----------------------------------------------------------------------------
# Comparison
# -----------------------------------------------------------------------------
comparison <- bind_rows(ts_mean_first, ts_pooled) %>%
  mutate(iqr_width = q75 - q25, ci95_width = q975 - q025)

save_figure_data(comparison, "figure_1_weekly_ts_pooled_vs_meanfirst.csv")

# Summary: how much narrower is mean_first vs pooled, on average, by metric?
width_summary <- comparison %>%
  select(scenario, arm, week, metric, method, ci95_width, iqr_width) %>%
  pivot_wider(names_from = method, values_from = c(ci95_width, iqr_width)) %>%
  mutate(
    ci95_ratio = ci95_width_mean_first / ci95_width_pooled,   # < 1 => mean_first is narrower
    iqr_ratio  = iqr_width_mean_first  / iqr_width_pooled
  ) %>%
  group_by(scenario, metric) %>%
  summarise(
    mean_ci95_ratio = mean(ci95_ratio, na.rm = TRUE),
    mean_iqr_ratio  = mean(iqr_ratio,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(scenario, metric)

message("Average width ratio (mean_first / pooled) across weeks, by scenario x metric:")
message("(values well below 1 = mean-first compresses the interval a lot)")
print(as.data.frame(width_summary))

save_figure_data(width_summary, "figure_1_pooled_vs_meanfirst_width_summary.csv")

message("Done. See figure_1_weekly_ts_pooled_vs_meanfirst.csv and figure_1_pooled_vs_meanfirst_width_summary.csv")
