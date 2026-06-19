# =============================================================================
# 03_plot_figure4.R
#
# Plots Figure 4 panels a, b, c from pre-computed CSVs.
# Run 02_extract_figure4.R first.
#
# Panel a : HCW deaths averted vs stockpile size
#           Policy B only, DPC 0 and DPC 5
#           IQR band + vertical IQR tick at STOCKPILE_HIGHLIGHT
#
# Panel b : % HCW deaths averted vs supply/demand ratio
#           DPC 0 only, Policy A (broad) vs Policy B (targeted)
#
# Panel c : doses per death averted vs intrinsic efficacy (log y)
#           DPC 0 and DPC 5, Policy A (dashed) vs Policy B (solid)
#           Same DPC shaded between the two policy lines
# =============================================================================

library(here)
library(ggplot2)
library(dplyr)
library(patchwork)

source(here("analyses", "03_figure_template", "helper_functions_figure_1to4.R"))

OUT_DIR <- here("figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Single color palette used across all panels
DPC_COLORS <- c(
  "0" = "#1D9E75",   # teal  — dpc 0
  "5" = "#C0392B"    # red   — dpc 5
)
DPC_LABELS <- c("0" = "prompt delivery (dpc 0)",
                "5" = "delayed delivery (dpc 5)")

POLICY_LINETYPES <- c(B = "solid", A = "dashed")
POLICY_LABELS    <- c(B = "B lower (solid)", A = "A upper (dashed)")

save_fig <- function(filename_base, plot, width, height) {
  ggsave(file.path(OUT_DIR, paste0(filename_base, ".png")),
         plot, width = width, height = height, dpi = 400, units = "in")
  ggsave(file.path(OUT_DIR, paste0(filename_base, ".pdf")),
         plot, width = width, height = height, units = "in")
}

# =============================================================================
# Load data
# =============================================================================
stockpile_df <- read.csv(here("output_figgen", "figure4_stockpile_summary.csv"),
                         stringsAsFactors = FALSE) %>%
  mutate(dpc_chr = as.character(dpc))

doses_df <- read.csv(here("output_figgen", "figure4_doses_per_death.csv"),
                     stringsAsFactors = FALSE) %>%
  mutate(dpc_chr = as.character(dpc))

# =============================================================================
# Panel a: HCW deaths averted vs stockpile size (Policy B only)
# =============================================================================
panel_a_df <- stockpile_df %>%
  filter(policy == "B") %>%
  rename(stockpile = stockpile_doses)

# Wide format for between-DPC ribbon in panel a
panel_a_wide <- panel_a_df %>%
  select(stockpile, dpc_chr, deaths_averted_med) %>%
  tidyr::pivot_wider(names_from = dpc_chr, values_from = deaths_averted_med,
                     names_prefix = "dpc") %>%
  filter(!is.na(dpc0) & !is.na(dpc5))

panel_a <- ggplot() +
  # Ribbon between DPC 0 median and DPC 5 median — light DPC 5 color
  geom_ribbon(data = panel_a_wide,
              aes(x = stockpile, ymin = dpc5, ymax = dpc0),
              fill = "#C0392B", alpha = 0.12, color = NA) +
  # Median lines per DPC
  geom_line(data = panel_a_df,
            aes(x = stockpile, y = deaths_averted_med, color = dpc_chr),
            linewidth = 1.0) +
  scale_color_manual(values = DPC_COLORS, labels = DPC_LABELS, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02)),
                     limits = c(0, 30000),
                     name = "stockpile size (doses)",
                     labels = scales::comma) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                     limits = c(0, NA),
                     name = "HCW deaths averted") +
  theme_fig() +
  theme(legend.position = c(0.65, 0.25),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 8))

# =============================================================================
# Panel b: % deaths averted vs supply/demand ratio (DPC 0, Policy A vs B)
# =============================================================================
panel_b_df <- stockpile_df %>%
  filter(dpc == 0)

# x-axis limit: just past 3x supply ratio
panel_b_xlim <- max(panel_b_df$supply_ratio_med[panel_b_df$supply_ratio_med <= 3.3],
                    na.rm = TRUE) * 1.05

# Wide format for ribbon between A and B
panel_b_wide <- panel_b_df %>%
  select(supply_ratio_med, policy, pct_averted_med) %>%
  tidyr::pivot_wider(names_from = policy, values_from = pct_averted_med) %>%
  filter(!is.na(A) & !is.na(B))

panel_b <- ggplot(panel_b_df,
                  aes(x = supply_ratio_med, y = pct_averted_med,
                      linetype = policy)) +
  # Shade between A and B to show "value of targeting"
  geom_ribbon(data = panel_b_wide,
              aes(x = supply_ratio_med, ymin = A, ymax = B),
              inherit.aes = FALSE,
              fill = DPC_COLORS[["0"]], alpha = 0.12) +
  # Both policies in same DPC color, distinguished by linetype only
  geom_line(color = DPC_COLORS[["0"]], linewidth = 1.0) +
  # 1x reference line
  geom_vline(xintercept = 1, linetype = "dashed",
             color = "grey60", linewidth = 0.4) +
  scale_linetype_manual(values = POLICY_LINETYPES,
                        labels = c(B = "Policy B (targeted)",
                                   A = "Policy A (broad)"),
                        name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02)),
                     limits = c(0, panel_b_xlim),
                     breaks = 0:3,
                     labels = c("0", "1×", "2×", "3×"),
                     name = "supply (× targeted demand)") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                     limits = c(0, NA),
                     name = "% HCW deaths averted") +
  theme_fig() +
  theme(legend.position = c(0.65, 0.25),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.key.width = unit(1.5, "cm"),
        legend.text = element_text(size = 8)) +
  guides(linetype = guide_legend(
    override.aes = list(linewidth = 1.2, color = DPC_COLORS[["0"]])
  ))

# =============================================================================
# Panel c: doses per death averted vs intrinsic efficacy (log y)
#          DPC 0 (teal) and DPC 5 (amber/red)
#          Policy B = solid, Policy A = dashed; same DPC shaded between them
# =============================================================================
# Build wide data for ribbons (one row per dpc x efficacy_scale)
ribbon_df <- doses_df %>%
  select(dpc_chr, efficacy_scale, intrinsic_efficacy, policy,
         doses_per_death_med) %>%
  tidyr::pivot_wider(names_from = policy,
                     values_from = doses_per_death_med) %>%
  filter(!is.na(A) & !is.na(B))

panel_c <- ggplot() +
  # Shaded band between Policy B (solid, lower) and Policy A (dashed, higher)
  geom_ribbon(data = ribbon_df,
              aes(x = intrinsic_efficacy * 100,
                  ymin = B, ymax = A,
                  fill = dpc_chr),
              alpha = 0.18) +
  # Policy B lines (solid)
  geom_line(data = doses_df %>% filter(policy == "B"),
            aes(x = intrinsic_efficacy * 100,
                y = doses_per_death_med,
                color = dpc_chr),
            linetype = "solid", linewidth = 1.0) +
  # Policy A lines (dashed)
  geom_line(data = doses_df %>% filter(policy == "A"),
            aes(x = intrinsic_efficacy * 100,
                y = doses_per_death_med,
                color = dpc_chr),
            linetype = "dashed", linewidth = 1.0) +
  scale_color_manual(values = DPC_COLORS, guide = "none") +
  scale_fill_manual(values  = DPC_COLORS, guide = "none") +
  scale_x_continuous(breaks = seq(20, 90, by = 10),
                     labels = function(x) paste0(x, "%"),
                     name = "intrinsic efficacy (%)") +
  scale_y_log10(name = "doses / death averted (log)",
                breaks = c(20, 50, 100, 300, 1000, 3000),
                labels = scales::comma) +
  theme_fig() +
  theme(legend.position = "none")

# =============================================================================
# Combine and save
# =============================================================================
fig4_all <- (panel_a | panel_b | panel_c) +
  plot_annotation(tag_levels = list(c("a ", "b ", "c ")))

save_fig("figure_4", fig4_all, 10, 3)
message("Figure 4 saved.")

save_fig("figure_4_panel-a_stockpile",      panel_a, 5, 5)
save_fig("figure_4_panel-b_supply-demand",  panel_b, 5, 5)
save_fig("figure_4_panel-c_doses-per-death", panel_c, 5, 5)
message("Figure 4 split panels saved.")