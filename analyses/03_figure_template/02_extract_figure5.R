# =============================================================================
# 02_extract_figure5.R
# Extract and save dose summaries for Figure 5.
# Output: output_figgen/figure_5_dose_summary.csv
#         output_figgen/figure_5_ppe_by_particle.csv
# =============================================================================
source(here::here("analyses", "03_figure_template", "helper_functions_figure_1to4.R"))
source(here("functions", "abc_posterior.R"))

FIG5_EFFICACY_LEVELS <- c("obv_20", "obv_30", "obv_40",
                          "obv_50", "obv_60", "obv_70", "obv_80")

# Load PPE efficacy per particle from posterior RDS
SCENARIOS_RDS <- list(
  WestAfrica = here("outputs", "02_ABC_model_fits_Final",
                    "fiber_ABC_SMC_Worst_WestAfrica_Decoupled_20260608_162044_check_NP5_NS4_NBREPS_30_NBSIMUL_472.RDS"),
  DRC        = here("outputs", "02_ABC_model_fits_Final",
                    "fiber_ABC_SMC_Middle_DRC_ConflictSmoothed_PlusPlus_Decoupled_20260607_215621_check_NP5_NS4_NBREPS_30_NBSIMUL_590.RDS")
)
PARAM_NAMES   <- c("R0", "prop_funeral", "etu_efficacy", "ppe_efficacy", "hcw_risk_scalar")
RESAMPLE_SEED <- 42L
N_PARTICLES   <- 200L

ppe_by_particle <- do.call(rbind, lapply(names(SCENARIOS_RDS), function(sc_name) {
  res       <- readRDS(SCENARIOS_RDS[[sc_name]])
  posterior <- as.data.frame(res$param)
  colnames(posterior) <- PARAM_NAMES
  posterior$weight    <- res$weights
  theta <- downsample_posterior(posterior, n_sets = N_PARTICLES,
                                seed = RESAMPLE_SEED, param_names = PARAM_NAMES)
  data.frame(scenario = sc_name, particle_id = seq_len(N_PARTICLES),
             ppe_efficacy = theta$ppe_efficacy, stringsAsFactors = FALSE)
}))

save_figure_data(ppe_by_particle, "figure_5_ppe_by_particle.csv")

message("Extracting dose summaries for figure 5...")

dose_df <- do.call(rbind, lapply(FIG5_EFFICACY_LEVELS, function(eff_name) {
  arm_dir <- sprintf("full_obv%02d", round(OBV_EFFICACY_VALUES[[eff_name]] * 100))
  extract_dose_summary(arm_dir, eff_name = eff_name, n_workers = 10L)
}))

save_figure_data(dose_df, "figure_5_dose_summary.csv")
message("Figure 5 data extraction complete.")