# =============================================================================
# 02_extract_figure_supp_sensitivity_numbers.R
#
# Pulls every number needed for the Section 3.9 sensitivity-analysis text
# (Figures S13-S15). No text-generation step -- just the extracted values,
# printed with enough precision to verify by hand and paste into the
# manuscript in place of the [xx] placeholders.
# =============================================================================
library(here)
library(dplyr)

options(pillar.sigfig = 6)

summarise_pct <- function(df, group_vars) {
  df %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(
      med = median(pct_hcw_deaths_averted, na.rm = TRUE),
      lo  = quantile(pct_hcw_deaths_averted, 0.025, na.rm = TRUE),
      hi  = quantile(pct_hcw_deaths_averted, 0.975, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(across(c(med, lo, hi), ~round(.x, 1)))
}

# =============================================================================
# Figure S13: conflict-intensity sensitivity
# =============================================================================
cat("\n============================================================\n")
cat("Figure S13: conflict-intensity sensitivity\n")
cat("============================================================\n")

intensity_df <- read.csv(here("output_figgen", "figure_supp_conflict_intensity_particle_summary.csv"),
                         stringsAsFactors = FALSE)

# Central efficacy, both scenarios, all three conditions
cat("\n-- with_conflict_mid ('Delayed coverage + dosing', central efficacy) by condition --\n")
intensity_df %>%
  filter(arm == "with_conflict_mid") %>%
  summarise_pct(c("condition")) %>%
  print()

cat("\n-- dpc_conflict_mid ('Delayed dosing', central efficacy) by condition --\n")
intensity_df %>%
  filter(arm == "dpc_conflict_mid") %>%
  summarise_pct(c("condition")) %>%
  print()

# Full breakdown (all efficacy levels x both scenarios x all conditions),
# for the paragraph's supporting detail / SI table if needed
cat("\n-- Full breakdown: all arms x conditions --\n")
intensity_df %>%
  summarise_pct(c("arm", "condition")) %>%
  arrange(arm, condition) %>%
  print(n = Inf)

# =============================================================================
# Figure S14: DPC-delay (efficacy shift) sensitivity
# =============================================================================
cat("\n============================================================\n")
cat("Figure S14: DPC-delay sensitivity\n")
cat("============================================================\n")

delay_df <- read.csv(here("output_figgen", "figure_supp_delay_particle_summary.csv"),
                     stringsAsFactors = FALSE)

cat("\n-- with_conflict_mid ('Delayed coverage + dosing', central efficacy) by shift --\n")
delay_df %>%
  filter(arm == "with_conflict_mid") %>%
  summarise_pct(c("shift")) %>%
  arrange(shift) %>%
  print()

cat("\n-- dpc_conflict_mid ('Delayed dosing', central efficacy) by shift --\n")
delay_df %>%
  filter(arm == "dpc_conflict_mid") %>%
  summarise_pct(c("shift")) %>%
  arrange(shift) %>%
  print()

# Difference between shift=0 and shift=8, central efficacy, both scenarios
cat("\n-- Change from shift=0 to shift=8 (median pct averted, central efficacy) --\n")
delay_df %>%
  filter(arm %in% c("with_conflict_mid", "dpc_conflict_mid"), shift %in% c(0, 8)) %>%
  group_by(arm, shift) %>%
  summarise(med = median(pct_hcw_deaths_averted, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = shift, values_from = med, names_prefix = "shift") %>%
  mutate(delta_pp = round(shift8 - shift0, 1)) %>%
  print()

cat("\n-- Full breakdown: all arms x shifts --\n")
delay_df %>%
  summarise_pct(c("arm", "shift")) %>%
  arrange(arm, shift) %>%
  print(n = Inf)

# =============================================================================
# Figure S15: TV-parameter sensitivity
# =============================================================================
cat("\n============================================================\n")
cat("Figure S15: TV-parameter sensitivity\n")
cat("============================================================\n")

tv_df <- read.csv(here("output_figgen", "figure_supp_tv_particle_summary.csv"),
                  stringsAsFactors = FALSE)

# Raw HCW deaths (absolute counts) by condition, central efficacy
cat("\n-- Raw HCW deaths, with_conflict_mid, by condition --\n")
tv_df %>%
  filter(arm == "with_conflict_mid") %>%
  group_by(condition) %>%
  summarise(
    med = median(n_hcw_deaths, na.rm = TRUE),
    lo  = quantile(n_hcw_deaths, 0.025, na.rm = TRUE),
    hi  = quantile(n_hcw_deaths, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(c(med, lo, hi), ~round(.x, 0))) %>%
  arrange(match(condition, c("goodgood", "good", "baseline", "bad", "badbad"))) %>%
  print()

# Percentage averted by condition, central efficacy
cat("\n-- Pct HCW deaths averted, with_conflict_mid, by condition --\n")
tv_df %>%
  filter(arm == "with_conflict_mid") %>%
  summarise_pct(c("condition")) %>%
  arrange(match(condition, c("goodgood", "good", "baseline", "bad", "badbad"))) %>%
  print()

cat("\n-- Full breakdown: all arms x conditions (deaths + pct averted) --\n")
tv_df %>%
  group_by(arm, condition) %>%
  summarise(
    med_deaths  = round(median(n_hcw_deaths, na.rm = TRUE), 0),
    lo_deaths   = round(quantile(n_hcw_deaths, 0.025, na.rm = TRUE), 0),
    hi_deaths   = round(quantile(n_hcw_deaths, 0.975, na.rm = TRUE), 0),
    med_averted = round(median(pct_hcw_deaths_averted, na.rm = TRUE), 1),
    lo_averted  = round(quantile(pct_hcw_deaths_averted, 0.025, na.rm = TRUE), 1),
    hi_averted  = round(quantile(pct_hcw_deaths_averted, 0.975, na.rm = TRUE), 1),
    .groups = "drop"
  ) %>%
  arrange(arm, match(condition, c("goodgood", "good", "baseline", "bad", "badbad"))) %>%
  print(n = Inf)

message("\nDone. Paste the printed values into the [xx] placeholders in the Section 3.9 text.")