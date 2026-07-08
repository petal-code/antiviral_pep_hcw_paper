# =============================================================================
# 02_extract_figure_supp_dpcdelay.R
# Extract particle-level run summaries for the DPC-delay sensitivity sweep
# (01_analysis_figure3new_conflict_dpc_sensitivity_delay.R), plus shift=0
# (our default/main setting, no delay) for comparison.
#
# For each shift level (0 [baseline], 2, 4, 6, 8, 10 days), reads all 13 arms
# from outputs/simulation/conflict_dpc_sens_delay{shift}/ (shift=0 instead
# reuses outputs/simulation/conflict_dpc_max5/, the existing baseline run),
# computes particle-level pct_hcw_deaths_averted (via make_particle_df, using
# "no_pep" as the counterfactual denominator within that shift), and tags
# the result with the shift level.
#
# Output: output_figgen/figure_supp_delay_particle_summary.csv
# =============================================================================
source(here::here("analyses", "03_figure_template", "helper_functions_figure_1to4.R"))

SHIFTS <- c(0, 2, 4, 6, 8, 10)

# shift=0 is our default/main setting -- not a separately simulated
# "conflict_dpc_sens_delay00" folder, but the existing baseline run.
SHIFT_DIRS <- setNames(
  c("conflict_dpc_max5", sprintf("conflict_dpc_sens_delay%02d", setdiff(SHIFTS, 0))),
  as.character(SHIFTS)
)

ARM_NAMES <- c(
  "no_pep",
  "with_conflict_mid", "with_conflict_lo", "with_conflict_hi",
  "cov_conflict_mid",  "cov_conflict_lo",  "cov_conflict_hi",
  "dpc_conflict_mid",  "dpc_conflict_lo",  "dpc_conflict_hi",
  "optimistic_mid",    "optimistic_lo",    "optimistic_hi"
)

message("Extracting run summaries for DPC-delay sensitivity sweep...")

particle_df_delay <- do.call(rbind, lapply(SHIFTS, function(shift) {
  message(sprintf("  shift = +%d days...", shift))
  shift_dir <- SHIFT_DIRS[[as.character(shift)]]
  
  run_df <- do.call(rbind, lapply(ARM_NAMES, function(arm_name) {
    df <- extract_run_summary(
      arm_dir    = file.path(shift_dir, arm_name),
      arm_label  = arm_name,
      n_workers  = 10L,
      obv_return = FALSE
    )
    df$scenario <- "DRC"  # normalise: DRC_conflict -> DRC
    df
  }))
  
  particle_df <- make_particle_df(run_df)
  particle_df$shift <- shift
  particle_df
}))

save_figure_data(particle_df_delay, "figure_supp_delay_particle_summary.csv")
message("DPC-delay sensitivity extraction complete.")