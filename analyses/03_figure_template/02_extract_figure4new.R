# =============================================================================
# 02_extract_figure4.R
#
# Post-hoc extraction for Figure 4 panels a, b, c.
# Reads baseline (no antiviral) simulation RDS files and applies
# Policy A / B logic entirely in R.
#
# Outputs (saved to output_figgen/):
#   figure4_stockpile_summary.csv   : panels a and b
#   figure4_doses_per_death.csv     : panel c
# =============================================================================

library(here)
library(dplyr)

source(here("analyses", "03_figure_template", "helper_functions_figure_1to4.R"))
source(here("analyses", "03_figure_template", "helper_functions_figure4.R"))

SCENARIO   <- "WestAfrica"
CURVE_DAT  <- readRDS(here("data-processed", "DPC_fixed_efficacy_varied_d50.rds"))
N_WORKERS  <- 10L

# =============================================================================
# 1. Extract per-run baseline data
# =============================================================================
message("Step 1: Extracting per-run baseline data...")

run_df <- extract_figure4_posthoc(
  sc_name   = SCENARIO,
  n_workers = N_WORKERS
)

message(sprintf("  Extracted %d runs (%d particles x %d reps).",
                nrow(run_df),
                length(unique(run_df$particle_id)),
                length(unique(run_df$rep))))

# =============================================================================
# 2. Panel a / b: stockpile sweep
#    DPC 0 and DPC 5, efficacy from curve mid at each DPC
# =============================================================================
message("Step 2: Stockpile sweep for panels a and b...")

STOCKPILE_SEQ <- seq(1000, 100000, by = 1000)  # doses (1 course = 20 doses, up to 5000 courses)
DPC_PANEL_AB  <- c(0, 5)

stockpile_list <- lapply(DPC_PANEL_AB, function(dpc) {
  message(sprintf("  DPC = %d ...", dpc))
  apply_stockpile_posthoc(
    run_df        = run_df,
    stockpile_seq = STOCKPILE_SEQ,
    dpc           = dpc,
    curve_dat     = CURVE_DAT,
    seed          = 42L
  )
})

stockpile_raw <- do.call(rbind, stockpile_list)

stockpile_summary <- summarise_stockpile_panel(stockpile_raw, run_df)

save_figure_data(stockpile_summary, "figure4_stockpile_summary.csv")
message("  Stockpile summary saved.")

# =============================================================================
# 3. Panel c: doses per death averted
#    Efficacy scales 0.2 - 0.9, DPC 0 and 5, Policy A and B
# =============================================================================
message("Step 3: Doses per death averted for panel c...")

doses_summary <- compute_doses_per_death(
  run_df          = run_df,
  efficacy_scales = seq(0.2, 0.9, by = 0.1),
  dpc_vals        = c(0, 5),
  curve_dat       = CURVE_DAT,
  seed            = 42L
)

save_figure_data(doses_summary, "figure4_doses_per_death.csv")
message("  Doses per death summary saved.")

message("Figure 4 extraction complete.")