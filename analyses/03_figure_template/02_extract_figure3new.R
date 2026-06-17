# =============================================================================
# 02_extract_figure3new.R
# Extract weekly HCW death time series and particle-level summaries for the
# four coverage/DPC scenarios x three efficacy arms (figure3_new).
#
# Outputs:
#   figure_3new_weekly_ts.csv       : weekly incident + cumulative HCW deaths
#   figure_3new_particle_summary.csv : particle-level deaths averted % per arm
# =============================================================================
source(here::here("analyses", "03_figure_template", "helper_functions_figure_1to4.R"))

ARM_NAMES_3NEW <- c(
  "with_conflict_mid", "with_conflict_lo", "with_conflict_hi",
  "without_conflict_mid", "without_conflict_lo", "without_conflict_hi",
  "optimistic_mid", "optimistic_lo", "optimistic_hi",
  "no_pep"
)

# ---- 1. Weekly time series (incident + cumulative HCW deaths) ----
message("Extracting weekly time series for figure3_new arms...")

# raw_ts_3new <- do.call(rbind, lapply(ARM_NAMES_3NEW, function(arm_name) {
#   extract_weekly_ts(
#     arm_dir   = file.path("conflict_dpc", arm_name),
#     bin_width = 7,
#     n_workers = 10L
#   )
# }))

raw_ts_3new <- do.call(rbind, lapply(ARM_NAMES_3NEW, function(arm_name) {
  df <- extract_weekly_ts(
    arm_dir   = file.path("conflict_dpc", arm_name),
    bin_width = 7,
    n_workers = 10L
  )
  # extract_weekly_ts always labels rows "baseline"/"obv" internally;
  # overwrite with the actual arm name so multiple arms can be combined.
  # Keep only the "obv" rows, since each arm_name already corresponds to
  # a specific coverage/DPC/efficacy setting (no need for its own baseline).
  df <- df[df$arm == "obv", ]
  df$arm <- arm_name
  df
}))

ts_quantiles_3new <- raw_ts_3new %>%
  group_by(scenario, arm, particle_id, week, metric) %>%
  summarise(value = mean(value), .groups = "drop") %>%
  arrange(scenario, arm, particle_id, metric, week) %>%
  group_by(scenario, arm, particle_id, metric) %>%
  mutate(value = if (unique(metric) == "hcw_deaths") cumsum(value) else value) %>%
  ungroup() %>%
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

save_figure_data(ts_quantiles_3new, "figure_3new_weekly_ts.csv")

# ---- 2. Particle-level run summary (for deaths averted % per arm) ----
message("Extracting particle-level run summaries for figure3_new arms...")

run_df_3new <- do.call(rbind, lapply(ARM_NAMES_3NEW, function(arm_name) {
  extract_run_summary(
    arm_dir   = file.path("conflict_dpc", arm_name),
    arm_label = arm_name,
    n_workers = 10L
  )
}))

particle_df_3new <- make_particle_df(run_df_3new)

save_figure_data(particle_df_3new, "figure_3new_particle_summary.csv")

message("Figure 3new extraction complete.")