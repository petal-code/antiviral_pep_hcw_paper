# =============================================================================
# 02_extract_figure4.R
#
# WestAfrica and DRC processed completely independently.
# Stockpile grid is demand-adaptive: covers up to 4x median demand per
# scenario so that Policy B always reaches plateau within the grid.
# =============================================================================

library(here)
library(dplyr)

source(here("analyses", "03_figure_template", "helper_functions_figure_1to4.R"))
source(here("analyses", "03_figure_template", "helper_functions_figure4.R"))

CURVE_DAT    <- readRDS(here("data-processed", "DPC_fixed_efficacy_varied_d50.rds"))
N_WORKERS    <- 10L
DPC_PANEL_AB <- c(0, 5)
DOSES_PER_COURSE <- 20L

for (sc_name in c("WestAfrica", "DRC")) {
  message(sprintf("\n===== Processing %s =====", sc_name))
  
  # ---- Step 1: Extract baseline data ----
  message("Step 1: Extracting per-run baseline data...")
  run_df <- extract_figure4_posthoc(sc_name = sc_name, n_workers = N_WORKERS)
  message(sprintf("  Extracted %d runs.", nrow(run_df)))
  
  # Stockpile grid: fine enough to give smooth curves
  # Units: doses (1 course = 20 doses)
  stockpile_seq <- seq(20, 200000, by = 20)
  message(sprintf("  stockpile grid: %d to %d doses (%d steps)",
                  min(stockpile_seq), max(stockpile_seq), length(stockpile_seq)))
  
  # ---- Step 2: Stockpile sweep (panels a and b) ----
  message("Step 2: Stockpile sweep...")
  stockpile_raw <- do.call(rbind, lapply(DPC_PANEL_AB, function(dpc) {
    message(sprintf("  DPC = %d ...", dpc))
    apply_stockpile_posthoc(
      run_df        = run_df,
      stockpile_seq = stockpile_seq,
      dpc           = dpc,
      curve_dat     = CURVE_DAT
    )
  }))
  
  stockpile_panels <- summarise_stockpile_panel(stockpile_raw, run_df)
  
  save_figure_data(stockpile_panels$panel_a,
                   sprintf("figure4_%s_panel_a_summary.csv", sc_name))
  save_figure_data(stockpile_panels$panel_b,
                   sprintf("figure4_%s_panel_b_summary.csv", sc_name))
  message("  Stockpile summaries saved.")
  
  # ---- Step 3: Doses per death (panel c) ----
  message("Step 3: Doses per death averted...")
  doses_summary <- compute_doses_per_death(
    run_df          = run_df,
    efficacy_scales = seq(0.2, 0.9, by = 0.1),
    dpc_vals        = c(0, 5),
    curve_dat       = CURVE_DAT
  )
  
  save_figure_data(doses_summary,
                   sprintf("figure4_%s_doses_per_death.csv", sc_name))
  message("  Doses per death summary saved.")
}

message("\nFigure 4 extraction complete.")