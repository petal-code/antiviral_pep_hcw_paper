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

make_panel_b <- function(sc) {
  df <- panel_b_raw %>%
    filter(scenario == sc, dpc == 0, supply_ratio <= 3.3)
  
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
                       breaks = 0:3,
                       labels = c("0", "1×", "2×", "3×"),
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
#
#   Col 1 (left)   : West Africa Archetype (top) / DRC Archetype (bottom)
#                    shared "HCW deaths averted" y-label centred across both
#                    no individual y-axis titles on the panels
#   Col 2 (centre) : West Africa panel b — no title
#   Col 3 (right)  : West Africa panel c — no title
#   Legends: transparent background
# =============================================================================
XLIM_A["DRC"] <- 10000   # restore

wa_a_alt2 <- make_panel_a("WestAfrica") +
  labs(title = "West Africa Archetype") +
  theme(plot.title     = element_text(size = 9, face = "bold"),
        axis.title.y   = element_blank(),
        legend.position   = c(0.65, 0.25),
        legend.background = element_blank(),
        legend.key        = element_blank())

drc_a_alt2 <- make_panel_a("DRC") +
  labs(title = "DRC Archetype") +
  theme(plot.title     = element_text(size = 9, face = "bold"),
        axis.title.y   = element_blank(),
        legend.position = "none")

wa_b_alt2 <- make_panel_b("WestAfrica") +
  labs(title = NULL) +
  theme(legend.position   = c(0.65, 0.25),
        legend.background = element_blank(),
        legend.key        = element_blank(),
        legend.key.width  = unit(1.0, "cm"))

wa_c_alt2 <- make_panel_c("WestAfrica") +
  labs(title = NULL)

# Shared y-axis label centred across the two left panels
y_label <- wrap_elements(
  textGrob("HCW deaths averted", rot = 90,
           gp = gpar(fontsize = 9))
)

# Wrap y_label so patchwork skips it for tagging
fig4_alt2 <- (y_label | (wa_a_alt2 / drc_a_alt2) | wa_b_alt2 | wa_c_alt2) +
  plot_layout(widths = c(0.04, 1, 1, 1)) +
  plot_annotation(tag_levels = "a",
                  tag_prefix = "",
                  # skip the first element (y_label spacer) in tagging
  ) &
  theme(plot.tag.position = "topleft")

# patchwork tags the spacer too — override by manually setting tags
fig4_alt2[[1]] <- fig4_alt2[[1]] + labs(tag = "")

save_fig("figure_4_alt2", fig4_alt2, 14, 7)
message("Figure 4 Alt2 saved.")