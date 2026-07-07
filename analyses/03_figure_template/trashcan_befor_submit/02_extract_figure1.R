# =============================================================================
# 01_extract_figure1.R
# Extract and save weekly time series data for Figure 1.
# Run this once after simulations are complete.
# Output saved to output_figgen/figure_1_weekly_ts.csv
# =============================================================================
source(here::here("analyses", "03_figure_template", "helper_functions_figure_1to4.R"))

raw_ts <- extract_weekly_ts("cov80_obv80", bin_width = 7, n_workers = 10L)

INCIDENCE_METRICS <- c("deaths", "infections", "hcw_deaths_incidence")
# hcw_deaths is cumulative; all others are incidence

ts_quantiles <- raw_ts %>%
  # Step 1: mean over reps per particle
  group_by(scenario, arm, particle_id, week, metric) %>%
  summarise(value = mean(value), .groups = "drop") %>%
  # Step 2: cumsum for hcw_deaths per particle (after rep averaging)
  arrange(scenario, arm, particle_id, metric, week) %>%
  group_by(scenario, arm, particle_id, metric) %>%
  mutate(value = if (unique(metric) == "hcw_deaths") cumsum(value) else value) %>%
  ungroup() %>%
  # Step 3: quantiles over particles
  group_by(scenario, arm, week, metric) %>%
  summarise(
    q025 = quantile(value, 0.025),
    q25  = quantile(value, 0.25),
    q50  = quantile(value, 0.50),
    q75  = quantile(value, 0.75),
    q975 = quantile(value, 0.975),
    .groups = "drop"
  ) %>%
  mutate(week = week / 7)

save_figure_data(ts_quantiles, "figure_1_weekly_ts.csv")

# Particle-level final cumulative HCW deaths (baseline vs obv), for computing
# averted (baseline - obv) and % reduction with proper CIs at the particle level.
particle_cum_hcw_deaths <- raw_ts %>%
  filter(metric == "hcw_deaths") %>%
  group_by(scenario, arm, particle_id, week) %>%
  summarise(value = mean(value), .groups = "drop") %>%
  arrange(scenario, arm, particle_id, week) %>%
  group_by(scenario, arm, particle_id) %>%
  mutate(value = cumsum(value)) %>%
  filter(week == max(week)) %>%
  ungroup() %>%
  select(scenario, arm, particle_id, cum_hcw_deaths = value)

save_figure_data(particle_cum_hcw_deaths, "figure_1_particle_cum_hcw_deaths.csv")

message("Figure 1 data extraction complete.")