# =============================================================================
# 02_plot_figure4.R
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
DPC_LABELS <- c("0" = "Prompt delivery (DPC 0)", "5" = "Delayed delivery (DPC 5)")
POLICY_LINETYPES <- c(B = "solid", A = "dashed")

XLIM_A <- c(WestAfrica = 30000, DRC = 2500)

SC_ORDER  <- c("WestAfrica", "DRC")
SC_LABELS <- c(WestAfrica = "West Africa (Worst)", DRC = "DRC (Middle, PlusPlus)")

save_fig <- function(filename_base, plot, width, height) {
  ggsave(file.path(OUT_DIR, paste0(filename_base, ".tiff")),
         plot, width = width, height = height, dpi = 400, units = "in")
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
                       name   = "Stockpile size (doses)") +
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
                       name   = "Dose supply relative to targeted demand") +
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

make_panel_c_dpc0_only <- function(sc) {
  df <- doses_df %>% filter(scenario == sc, dpc == 0)
  
  ribbon <- df %>%
    select(efficacy_scale, intrinsic_efficacy, policy, doses_per_death_med) %>%
    pivot_wider(id_cols = c(efficacy_scale, intrinsic_efficacy),
                names_from = policy, values_from = doses_per_death_med) %>%
    filter(!is.na(A) & !is.na(B))
  
  ggplot() +
    geom_ribbon(data = ribbon,
                aes(x = intrinsic_efficacy * 100, ymin = B, ymax = A),
                fill = DPC_COLORS[["0"]], alpha = 0.18) +
    geom_line(data = df %>% filter(policy == "B"),
              aes(x = intrinsic_efficacy * 100, y = doses_per_death_med),
              color = DPC_COLORS[["0"]], linetype = "solid", linewidth = 1.0) +
    geom_line(data = df %>% filter(policy == "A"),
              aes(x = intrinsic_efficacy * 100, y = doses_per_death_med),
              color = DPC_COLORS[["0"]], linetype = "dashed", linewidth = 1.0) +
    scale_x_continuous(breaks = seq(20, 90, by = 10),
                       labels = function(x) paste0(x, "%"),
                       name   = "Antiviral efficacy (%)") +
    scale_y_log10(name   = "doses / death averted (log)",
                  breaks = c(20, 50, 100, 300, 1000, 3000),
                  labels = scales::comma) +
    theme_fig() +
    theme(legend.position = "none")
}

# =============================================================================
# figure_4_alt3: 3-column sketch-style layout, panel d (doses/death) shows
# prompt delivery only (dpc 0); delayed delivery (dpc 5) lines removed.
# =============================================================================
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

wa_c_alt3 <- make_panel_c_dpc0_only("WestAfrica") +
  labs(title = NULL, tag = "d")

# Shared y-axis label
y_label <- wrap_elements(
  textGrob("HCW deaths averted", rot = 90,
           gp = gpar(fontsize = 9))
)

# Isolate the left column to lock the y-axis label centering and prevent flattening errors
left_col <- wrap_plots(y_label, (wa_a_alt2 / drc_a_alt2), widths = c(0.04, 1))

fig4_alt3 <- (left_col | wa_b_alt2 | wa_c_alt3) +
  plot_layout(widths = c(1.04, 1, 1)) &
  theme(plot.tag.position = "topleft")

save_fig("figure_4_alt3", fig4_alt3, 10, 4)
message("Figure 4 Alt3 saved.")