# =============================================================================
# 02_extract_figure1.R
# Extract and save weekly time series data for Figure 1.
# Run this once after simulations are complete.
# Output: output_figgen/figure_1_weekly_ts.csv
#         output_figgen/figure_1_particle_cum_hcw_deaths.csv
#
# Pooled aggregation (variance-compression fix):
#   Quantiles are taken directly across all (particle_id, rep) pairs -- every
#   stochastic replicate is its own observation -- instead of first averaging
#   the N_REPS replicates within each particle and then taking quantiles
#   across only the resulting particle-level means. The old mean-first
#   approach compressed replicate-level (stochastic) uncertainty out of the
#   reported interval and left only parameter (posterior-draw) uncertainty.
#   See 02_extract_figure1_pooled_vs_meanfirst_diagnostic.R for the earlier
#   side-by-side comparison of the two methods.
# =============================================================================
source(here::here("analyses", "03_figure_template", "helper_functions_figure_1to4.R"))

raw_ts <- extract_weekly_ts("cov80_obv80", bin_width = 7, n_workers = 10L)
INCIDENCE_METRICS <- c("deaths", "infections", "hcw_deaths_incidence")

# hcw_deaths is cumulative; all others are incidence.
# cumsum is computed per (particle_id, rep) so each replicate keeps its own
# trajectory, then quantiles are taken across ALL (particle_id, rep) pairs.
ts_quantiles <- raw_ts %>%
  arrange(scenario, arm, particle_id, rep, metric, week) %>%
  group_by(scenario, arm, particle_id, rep, metric) %>%
  mutate(value = if (unique(metric) == "hcw_deaths") cumsum(value) else value) %>%
  ungroup() %>%
  group_by(scenario, arm, week, metric) %>%
  summarise(
    q025 = quantile(value, 0.025),
    q5 = quantile(value, 0.05),
    q10 = quantile(value, 0.10),
    q20 = quantile(value, 0.20),
    q25  = quantile(value, 0.25),
    q50  = quantile(value, 0.50),
    mean  = mean(value),
    q75  = quantile(value, 0.75),
    q80  = quantile(value, 0.80),
    q90  = quantile(value, 0.90),
    q95  = quantile(value, 0.95),
    q975 = quantile(value, 0.975),
    .groups = "drop"
  ) %>%
  mutate(week = week / 7)
save_figure_data(ts_quantiles, "figure_1_weekly_ts.csv")

# Replicate-level final cumulative HCW deaths (baseline vs obv), pooled across
# (particle_id, rep), for computing averted (baseline - obv) and % reduction
# with credible intervals that reflect both parameter and stochastic
# uncertainty.
particle_cum_hcw_deaths <- raw_ts %>%
  filter(metric == "hcw_deaths") %>%
  arrange(scenario, arm, particle_id, rep, week) %>%
  group_by(scenario, arm, particle_id, rep) %>%
  mutate(value = cumsum(value)) %>%
  filter(week == max(week)) %>%
  ungroup() %>%
  select(scenario, arm, particle_id, rep, cum_hcw_deaths = value)
save_figure_data(particle_cum_hcw_deaths, "figure_1_particle_cum_hcw_deaths.csv")

message("Figure 1 data extraction complete (pooled particle x rep quantiles).")
