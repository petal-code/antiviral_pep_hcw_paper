# =============================================================================
# 02_extract_figure_supp_conflict_intensity.R
# Extract particle-level run summaries for the conflict-INTENSITY sensitivity
# sweep (01_analysis_figure3new_conflict_dpc_sensitivity_intensity.R), plus
# the baseline condition for comparison.
#
#   baseline : outputs/simulation/conflict_dpc_max5/      (coverage 80, dpc 4)
#   weak     : outputs/simulation/conflict_dpc_sens_weak/   (coverage 100, dpc 1)
#   strong   : outputs/simulation/conflict_dpc_sens_strong/ (coverage 60, dpc 7)
#
# For each condition, reads all 13 arms, computes particle-level
# pct_hcw_deaths_averted (via make_particle_df, using "no_pep" as the
# counterfactual denominator within that condition), and tags the result
# with the condition label.
#
# Output: output_figgen/figure_supp_conflict_intensity_particle_summary.csv
# =============================================================================
source(here::here("analyses", "03_figure_template", "helper_functions_figure_1to4.R"))

CONDITIONS <- list(
  baseline = "conflict_dpc_max5",
  weak     = "conflict_dpc_sens_weak",
  strong   = "conflict_dpc_sens_strong"
)

ARM_NAMES <- c(
  "no_pep",
  "with_conflict_mid", "with_conflict_lo", "with_conflict_hi",
  "cov_conflict_mid",  "cov_conflict_lo",  "cov_conflict_hi",
  "dpc_conflict_mid",  "dpc_conflict_lo",  "dpc_conflict_hi",
  "optimistic_mid",    "optimistic_lo",    "optimistic_hi"
)

message("Extracting run summaries for conflict-intensity sensitivity sweep...")

particle_df_intensity <- do.call(rbind, lapply(names(CONDITIONS), function(cond_name) {
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
  
  particle_df <- make_particle_df(run_df)
  particle_df$condition <- cond_name
  particle_df
}))

save_figure_data(particle_df_intensity, "figure_supp_conflict_intensity_particle_summary.csv")
message("Conflict-intensity sensitivity extraction complete.")