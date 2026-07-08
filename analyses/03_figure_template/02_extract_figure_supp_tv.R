# =============================================================================
# 02_extract_figure_supp_tv.R
# Extract particle-level run summaries for the time-varying-parameter
# sensitivity sweep (01_analysis_figure3new_tv_param_sensitivity.R), plus
# the baseline (unperturbed) condition for comparison.
#
#   baseline : outputs/simulation/conflict_dpc_max5/   (existing baseline run)
#   good     : outputs/simulation/tv_sens_good/
#   goodgood : outputs/simulation/tv_sens_goodgood/
#   bad      : outputs/simulation/tv_sens_bad/
#   badbad   : outputs/simulation/tv_sens_badbad/
#
# For each condition, reads the three with_conflict arms (mid/lo/hi),
# computes:
#   - particle-level pct_hcw_deaths_averted (via make_particle_df; no
#     separate no_pep arm needed -- averted outcomes come from each run's
#     own prevented_completed data)
#   - particle-level raw HCW death counts (mean over reps within particle)
# and tags the result with the condition label.
#
# Output: output_figgen/figure_supp_tv_particle_summary.csv
# =============================================================================
source(here::here("analyses", "03_figure_template", "helper_functions_figure_1to4.R"))

CONDITIONS <- list(
  baseline = "conflict_dpc_max5",
  good     = "tv_sens_good",
  goodgood = "tv_sens_goodgood",
  bad      = "tv_sens_bad",
  badbad   = "tv_sens_badbad"
)

ARM_NAMES <- c("with_conflict_mid", "with_conflict_lo", "with_conflict_hi")

message("Extracting run summaries for TV-parameter sensitivity sweep...")

particle_df_tv <- do.call(rbind, lapply(names(CONDITIONS), function(cond_name) {
  message(sprintf("  condition = %s...", cond_name))
  cond_dir <- CONDITIONS[[cond_name]]
  
  run_df <- do.call(rbind, lapply(ARM_NAMES, function(arm_name) {
    df <- extract_run_summary(
      arm_dir    = file.path(cond_dir, arm_name),
      arm_label  = arm_name,
      n_workers  = 10L,
      obv_return = FALSE
    )
    df$scenario <- "DRC"  # normalise: DRC_conflict -> DRC
    df
  }))
  
  # Averted % (per-run prevented/counterfactual, no no_pep arm needed)
  particle_df <- make_particle_df(run_df) %>%
    select(scenario, particle_id, arm, pct_hcw_deaths_averted, pct_days_lost_averted)
  
  # Raw HCW death counts, mean over reps within each particle x arm
  deaths_df <- run_df %>%
    group_by(scenario, particle_id, arm) %>%
    summarise(n_hcw_deaths = mean(n_hcw_deaths, na.rm = TRUE), .groups = "drop")
  
  out <- left_join(particle_df, deaths_df, by = c("scenario", "particle_id", "arm"))
  out$condition <- cond_name
  out
}))

save_figure_data(particle_df_tv, "figure_supp_tv_particle_summary.csv")
message("TV-parameter sensitivity extraction complete.")