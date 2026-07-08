# =============================================================================
# 02_extract_figure4.R
#
# Pulls every number needed for the Figure 4 results paragraph.
# No text-generation step -- just the extracted values, printed with enough
# precision to verify by hand. Re-run any time the underlying CSVs change.
#
# Panel letters below match figure_4_alt3's actual tags:
#   a = West Africa panel a (deaths averted vs stockpile)
#   b = DRC panel a (deaths averted vs stockpile)
#   c = West Africa panel b (Policy A vs B, % HCW deaths averted)
#   d = West Africa doses/death (dpc 0 only)
# =============================================================================

library(here)
library(dplyr)
library(tidyr)

options(pillar.sigfig = 6)   # show enough decimals to verify differences by hand

# ---- Load data (same CSVs as 02_plot_figure4.R) -----------------------------
read_sc <- function(sc, file_suffix) {
  read.csv(here("output_figgen",
                sprintf("figure4_%s_%s", sc, file_suffix)),
           stringsAsFactors = FALSE) %>%
    mutate(scenario = sc,
           dpc_chr  = as.character(dpc))
}

SC_ORDER <- c("WestAfrica", "DRC")

panel_a_raw <- bind_rows(
  read_sc("WestAfrica", "panel_a_summary.csv"),
  read_sc("DRC",        "panel_a_summary.csv")
) %>% mutate(scenario = factor(scenario, levels = SC_ORDER))

panel_b_raw <- bind_rows(
  read_sc("WestAfrica", "panel_b_summary.csv"),
  read_sc("DRC",        "panel_b_summary.csv")
) %>% mutate(scenario = factor(scenario, levels = SC_ORDER))

doses_df <- bind_rows(
  read_sc("WestAfrica", "doses_per_death.csv"),
  read_sc("DRC",        "doses_per_death.csv")
) %>% mutate(scenario = factor(scenario, levels = SC_ORDER))

# Helper: pull the value at the closest available x to a target, warning if
# it isn't an exact grid match.
nearest_value <- function(df, xcol, target, ycol) {
  i <- which.min(abs(df[[xcol]] - target))
  if (df[[xcol]][i] != target) {
    warning(sprintf("No exact match for %s = %s; using nearest value %s",
                    xcol, target, df[[xcol]][i]))
  }
  df[[ycol]][i]
}

# =============================================================================
# 1) Deaths averted at 30,000 doses, dpc0 vs dpc5, Policy B
#    (Panels a/b) -- absolute counts AND the unrounded difference, both
#    kept, so you can see exactly how the rounded difference was derived.
# =============================================================================
stockpile_target <- 30000

deaths_at_30k <- panel_a_raw %>%
  filter(policy == "B") %>%
  group_by(scenario, dpc) %>%
  slice_min(abs(stockpile_doses - stockpile_target), n = 1) %>%
  ungroup() %>%
  select(scenario, dpc, stockpile_doses, deaths_averted_med) %>%
  pivot_wider(id_cols = scenario, names_from = dpc,
              values_from = deaths_averted_med, names_prefix = "dpc") %>%
  mutate(
    dpc0_rounded         = round(dpc0),
    dpc5_rounded         = round(dpc5),
    diff_unrounded        = dpc0 - dpc5,             # exact, pre-rounding
    diff_of_rounded       = dpc0_rounded - dpc5_rounded,  # what you'd get rounding first
    deaths_lost_to_delay  = round(diff_unrounded)    # the value to actually use
  )

# =============================================================================
# 2) Fold-more-doses for Policy A vs B at matched % HCW deaths averted
#    (Panel c, West Africa only, restricted to the plotted x-range)
# =============================================================================
SUPPLY_XLIM <- 2.5   # keep in sync with the panel's x-axis limit

panel_c_data <- panel_b_raw %>%
  filter(scenario == "WestAfrica", dpc == 0, supply_ratio <= SUPPLY_XLIM)

target_pct <- nearest_value(panel_c_data %>% filter(policy == "B"),
                            "supply_ratio", 1, "pct_averted_med")

policyA_match <- panel_c_data %>%
  filter(policy == "A") %>%
  slice_min(abs(pct_averted_med - target_pct), n = 1)

fold_more_doses     <- round(policyA_match$supply_ratio / 1, 1)
policyA_pct_reached  <- round(policyA_match$pct_averted_med, 1)
pct_gap              <- round(target_pct - policyA_pct_reached, 1)

if (abs(pct_gap) > 2) {
  warning(sprintf(
    paste("Policy A reaches only %.1f%% (vs target %.1f%%) within",
          "supply_ratio <= %s -- fold_more_doses is a lower bound, not an",
          "exact match."),
    policyA_pct_reached, target_pct, SUPPLY_XLIM))
}




# =============================================================================
# 3) Doses per HCW death averted at 80% intrinsic efficacy
#    (Panel d, West Africa)
# =============================================================================
eff_target <- 0.8

doses_at_80pct <- doses_df %>%
  filter(scenario == "WestAfrica") %>%
  group_by(dpc, policy) %>%
  slice_min(abs(intrinsic_efficacy - eff_target), n = 1) %>%
  ungroup() %>%
  select(dpc, policy, intrinsic_efficacy,
         doses_per_death_med_unrounded = doses_per_death_med) %>%
  mutate(doses_per_death_med = round(doses_per_death_med_unrounded))



# =============================================================================
# Print everything
# =============================================================================
cat("=== 1) Deaths averted at 30,000 doses (Policy B) ===\n")
print(deaths_at_30k)

cat("\n=== 2) Fold-more-doses, Policy A vs B (West Africa) ===\n")
print(tibble(target_pct, policyA_pct_reached, pct_gap, fold_more_doses))

cat("\n=== 3) Doses per HCW death averted at 80% efficacy (West Africa) ===\n")
print(doses_at_80pct)