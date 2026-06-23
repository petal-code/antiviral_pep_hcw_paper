# =============================================================================
# 03_plot_figure4.R
# =============================================================================

library(here)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(grid)

source(here("analyses", "03_figure_template", "helper_functions_figure_1to4.R"))

OUT_DIR <- here("figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

DPC_COLORS <- c("0" = "#1D9E75", "5" = "#C0392B")
DPC_LABELS <- c("0" = "prompt delivery (dpc 0)", "5" = "delayed delivery (dpc 5)")
POLICY_LINETYPES <- c(B = "solid", A = "dashed")

XLIM_A <- c(WestAfrica = 30000, DRC = 10000)

SC_ORDER  <- c("WestAfrica", "DRC")
SC_LABELS <- c(WestAfrica = "West Africa (Worst)", DRC = "DRC (Middle, PlusPlus)")

save_fig <- function(filename_base, plot, width, height) {
  ggsave(file.path(OUT_DIR, paste0(filename_base, ".png")),
         plot, width = width, height = height, dpi = 400, units = "in")
  ggsave(file.path(OUT_DIR, paste0(filename_base, ".pdf")),
         plot, width = width, height = height, units = "in")
}

# =============================================================================
# Load data
# =============================================================================
read_sc <- function(sc, file_suffix) {
  read.csv(here("output_figgen",
                sprintf("figure4_%s_%s", sc, file_suffix)),
           stringsAsFactors = FALSE) %>%
    mutate(scenario       = sc,
           dpc_chr        = as.character(dpc),
           scenario_label = SC_LABELS[sc])
}

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

# =============================================================================
# Panel builders
# =============================================================================

make_panel_a <- function(sc) {
  df <- panel_a_raw %>%
    filter(scenario == sc, policy == "B",
           stockpile_doses <= XLIM_A[sc]) %>%
    mutate(dpc_chr = as.character(dpc))
  
  wide <- df %>%
    select(stockpile_doses, dpc_chr, deaths_averted_med) %>%
    pivot_wider(id_cols = stockpile_doses,
                names_from = dpc_chr, values_from = deaths_averted_med,
                names_prefix = "dpc") %>%
    filter(!is.na(dpc0) & !is.na(dpc5))
  
  ggplot() +
    geom_ribbon(data = wide,
                aes(x = stockpile_doses, ymin = dpc5, ymax = dpc0),
                fill = "#C0392B", alpha = 0.12, color = NA) +
    geom_line(data = df,
              aes(x = stockpile_doses, y = deaths_averted_med, color = dpc_chr),
              linewidth = 1.0) +
    scale_color_manual(values = DPC_COLORS, labels = DPC_LABELS, name = NULL) +
    scale_x_continuous(labels = scales::comma,
                       expand = expansion(mult = c(0, 0.02)),
                       name   = "stockpile size (doses)") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                       limits = c(0, NA),
                       name   = "HCW deaths averted") +
    theme_fig() +
    theme(legend.position   = "bottom",
          legend.background = element_blank(),
          legend.key        = element_blank(),
          legend.text       = element_text(size = 8))
}

make_panel_b <- function(sc, xlim_max = 3.3, x_breaks = 0:3,
                         x_labels = c("0", "1×", "2×", "3×")) {
  df <- panel_b_raw %>%
    filter(scenario == sc, dpc == 0, supply_ratio <= xlim_max)
  
  wide <- df %>%
    select(supply_ratio, policy, pct_averted_med) %>%
    pivot_wider(id_cols = supply_ratio,
                names_from = policy, values_from = pct_averted_med) %>%
    filter(!is.na(A) & !is.na(B))
  
  ggplot(df, aes(x = supply_ratio, y = pct_averted_med, linetype = policy)) +
    geom_ribbon(data = wide,
                aes(x = supply_ratio, ymin = A, ymax = B),
                inherit.aes = FALSE,
                fill = DPC_COLORS[["0"]], alpha = 0.12) +
    geom_line(color = DPC_COLORS[["0"]], linewidth = 1.0) +
    geom_vline(xintercept = 1, linetype = "dashed",
               color = "grey60", linewidth = 0.4) +
    scale_linetype_manual(values = POLICY_LINETYPES,
                          labels = c(B = "Policy B (targeted)",
                                     A = "Policy A (broad)"),
                          name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.02)),
                       breaks = x_breaks,
                       labels = x_labels,
                       name   = "supply (× targeted demand)") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                       limits = c(0, NA),
                       name   = "% HCW deaths averted") +
    theme_fig() +
    theme(legend.position   = "bottom",
          legend.background = element_blank(),
          legend.key        = element_blank(),
          legend.key.width  = unit(1.5, "cm"),
          legend.text       = element_text(size = 8)) +
    guides(linetype = guide_legend(
      override.aes = list(linewidth = 1.2, color = DPC_COLORS[["0"]])
    ))
}

make_panel_c <- function(sc) {
  df <- doses_df %>% filter(scenario == sc)
  
  ribbon <- df %>%
    select(dpc_chr, efficacy_scale, intrinsic_efficacy,
           policy, doses_per_death_med) %>%
    pivot_wider(id_cols = c(dpc_chr, efficacy_scale, intrinsic_efficacy),
                names_from = policy, values_from = doses_per_death_med) %>%
    filter(!is.na(A) & !is.na(B))
  
  ggplot() +
    geom_ribbon(data = ribbon,
                aes(x = intrinsic_efficacy * 100, ymin = B, ymax = A,
                    fill = dpc_chr),
                alpha = 0.18) +
    geom_line(data = df %>% filter(policy == "B"),
              aes(x = intrinsic_efficacy * 100, y = doses_per_death_med,
                  color = dpc_chr),
              linetype = "solid", linewidth = 1.0) +
    geom_line(data = df %>% filter(policy == "A"),
              aes(x = intrinsic_efficacy * 100, y = doses_per_death_med,
                  color = dpc_chr),
              linetype = "dashed", linewidth = 1.0) +
    scale_color_manual(values = DPC_COLORS, guide = "none") +
    scale_fill_manual(values  = DPC_COLORS, guide = "none") +
    scale_x_continuous(breaks = seq(20, 90, by = 10),
                       labels = function(x) paste0(x, "%"),
                       name   = "intrinsic efficacy (%)") +
    scale_y_log10(name   = "doses / death averted (log)",
                  breaks = c(20, 50, 100, 300, 1000, 3000),
                  labels = scales::comma) +
    theme_fig() +
    theme(legend.position = "none")
}

# =============================================================================
# 1. Original 2x3 layout
# =============================================================================
wa_a <- make_panel_a("WestAfrica") + ggtitle(SC_LABELS["WestAfrica"]) +
  theme(plot.title = element_text(size = 9, face = "bold"))
wa_b <- make_panel_b("WestAfrica") + ggtitle(SC_LABELS["WestAfrica"]) +
  theme(plot.title = element_text(size = 9, face = "bold"))
wa_c <- make_panel_c("WestAfrica") + ggtitle(SC_LABELS["WestAfrica"]) +
  theme(plot.title = element_text(size = 9, face = "bold"))

drc_a <- make_panel_a("DRC") + ggtitle(SC_LABELS["DRC"]) +
  theme(plot.title = element_text(size = 9, face = "bold"))
drc_b <- make_panel_b("DRC") + ggtitle(SC_LABELS["DRC"]) +
  theme(plot.title = element_text(size = 9, face = "bold"))
drc_c <- make_panel_c("DRC") + ggtitle(SC_LABELS["DRC"]) +
  theme(plot.title = element_text(size = 9, face = "bold"))

fig4_all <- (wa_a | wa_b | wa_c) /
  (drc_a | drc_b | drc_c) +
  plot_annotation(tag_levels = "a")

save_fig("figure_4", fig4_all, 14, 10)
message("Figure 4 saved.")

# =============================================================================
# 2. Alt layout
# =============================================================================
XLIM_A["DRC"] <- 2500

wa_a_alt <- make_panel_a("WestAfrica") +
  theme(legend.position   = c(0.70, 0.25),
        legend.background = element_rect(fill = "white", color = "grey90",
                                         linewidth = 0.3))
drc_a_alt <- make_panel_a("DRC")
wa_b_alt  <- make_panel_b("WestAfrica") +
  theme(legend.position   = c(0.70, 0.25),
        legend.background = element_rect(fill = "white", color = "grey90",
                                         linewidth = 0.3),
        legend.key.width  = unit(1.0, "cm"))
wa_c_alt  <- make_panel_c("WestAfrica")

fig4_alt <- (wa_a_alt | drc_a_alt) /
  (wa_b_alt  | wa_c_alt) +
  plot_annotation(tag_levels = "a")

save_fig("figure_4_alt", fig4_alt, 11.5, 8.5)
message("Figure 4 Alt saved.")

# =============================================================================
# 3. Alt2 layout: 3-column sketch style
# =============================================================================
XLIM_A["DRC"] <- 2500   # Restore original limit
wa_a_alt2 <- make_panel_a("WestAfrica") +
  labs(title = "West Africa Archetype", tag = "a") +
  theme(plot.title        = element_text(size = 9, face = "bold"),
        axis.title.y      = element_blank(),
        axis.title.x      = element_blank(), # Removed x-axis label for the top panel
        legend.position   = c(0.65, 0.2),
        legend.background = element_blank(),
        legend.key        = element_blank(),
        plot.margin       = margin(t = 5, r = 5, b = -4, l = 5)) # Reduced bottom margin

drc_a_alt2 <- make_panel_a("DRC") +
  labs(title = "DRC Archetype", tag = "b") +
  theme(plot.title     = element_text(size = 9, face = "bold"),
        axis.title.y   = element_blank(),
        legend.position = "none",
        plot.margin    = margin(t = -4, r = 5, b = 5, l = 5)) # Reduced top margin to narrow the gap

wa_b_alt2 <- make_panel_b("WestAfrica",
                          xlim_max  = 2.5,
                          x_breaks  = c(0, 1, 2, 2.5),
                          x_labels  = c("0", "1×", "2×", "2.5×")) +
  labs(title = NULL, tag = "c") +
  theme(legend.position   = c(0.65, 0.25),
        legend.background = element_blank(),
        legend.key        = element_blank(),
        legend.key.width  = unit(1.0, "cm"))

wa_c_alt2 <- make_panel_c("WestAfrica") +
  labs(title = NULL, tag = "d")

# Shared y-axis label
y_label <- wrap_elements(
  textGrob("HCW deaths averted", rot = 90,
           gp = gpar(fontsize = 9))
)

# Isolate the left column to lock the y-axis label centering and prevent flattening errors
left_col <- wrap_plots(y_label, (wa_a_alt2 / drc_a_alt2), widths = c(0.04, 1))

# Combine the isolated left column with panel B and panel C using your original structure
fig4_alt2 <- (left_col | wa_b_alt2 | wa_c_alt2) +
  plot_layout(widths = c(1.04, 1, 1)) &
  theme(plot.tag.position = "topleft")

save_fig("figure_4_alt2", fig4_alt2, 10, 4)
message("Figure 4 Alt2 saved.")
# =============================================================================
# 04_extract_figure4_numbers.R
#
# Pulls every number needed for the Figure 4 results paragraph.
# No text-generation step -- just the extracted values, printed with enough
# precision to verify by hand. Re-run any time the underlying CSVs change.
# =============================================================================

library(here)
library(dplyr)
library(tidyr)

options(pillar.sigfig = 6)   # show enough decimals to verify differences by hand

# ---- Load data (same CSVs as 03_plot_figure4.R) -----------------------------
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
#    (Figure 4a/b) -- absolute counts AND the unrounded difference, both
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
#    (Figure 4c, West Africa only, restricted to the plotted x-range)
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
#    (Figure 4d, West Africa)
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