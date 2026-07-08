# =============================================================================
# 02_extract_figure_supp_C.R
# Extract run summaries for Figure Supp C (Figure 3 sensitivity analysis).
# Reads RDS files from outputs/simulation/fig3sens/
# Output: output_figgen/figure_supp_C_run_summary.csv
# =============================================================================
source(here::here("analyses", "03_figure_template", "helper_functions_figure_1to4.R"))
SENS_BASE <- here("outputs", "simulation", "fig3sens")
ONSET_DAYS <- seq(0L, 100L, by = 20L)
MAX_COVS   <- seq(0.20, 1.00, by = 0.20)
FIXED_EFFICACY <- 0.80
# Build arm table matching 01_analysis_figure3_sensitivity.R
ARMS_SENS <- do.call(rbind, lapply(ONSET_DAYS, function(onset) {
  do.call(rbind, lapply(MAX_COVS, function(max_cov) {
    data.frame(
      arm_name  = sprintf("sens_t%03d_cov%03d_obv80", onset, round(max_cov * 100)),
      onset_day = onset,
      max_cov   = max_cov,
      efficacy  = FIXED_EFFICACY,
      stringsAsFactors = FALSE
    )
  }))
}))
message(sprintf("Extracting %d sensitivity arms...", nrow(ARMS_SENS)))
sens_df <- do.call(rbind, lapply(seq_len(nrow(ARMS_SENS)), function(i) {
  arm      <- ARMS_SENS[i, ]
  arm_dir  <- file.path(SENS_BASE, arm$arm_name)
  
  # Use existing extract_run_summary from helper, pointing to fig3sens subdir
  df <- extract_run_summary(
    arm_dir    = file.path("fig3sens", arm$arm_name),
    arm_label  = arm$arm_name,
    n_workers  = 10L,
    obv_return = FALSE
  )
  df$onset_day <- arm$onset_day
  df$max_cov   <- arm$max_cov
  df$efficacy  <- arm$efficacy
  df
}))
save_figure_data(sens_df, "figure_supp_C_run_summary.csv")
message("Figure Supp C extraction complete.")