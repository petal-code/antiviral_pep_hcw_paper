# =============================================================================
# 02_plot_figure_supp_delay.R
# Visualise the DPC-delay (efficacy drop-off delay) sensitivity sweep,
# including shift=0 (our default/main setting) as a black reference line/bar.
# Reads: output_figgen/figure_supp_delay_particle_summary.csv
#        data-processed/DPC_fixed_efficacy_varied_d50.rds
# =============================================================================
source(here::here("analyses", "03_figure_template", "helper_functions_figure_1to4.R"))
OUT_DIR <- here("figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

save_fig <- function(filename_base, plot, width, height) {
  ggsave(file.path(OUT_DIR, paste0(filename_base, ".png")),
         plot, width = width, height = height, dpi = 400, units = "in")
  ggsave(file.path(OUT_DIR, paste0(filename_base, ".pdf")),
         plot, width = width, height = height, units = "in")
}

SHIFTS <- c(0, 2, 4, 6, 8)  # 10 dropped -- efficacy already fully saturated by then

ARM_PANEL_ORDER <- c(
  "with_conflict_mid", "with_conflict_lo", "with_conflict_hi",
  "dpc_conflict_mid",  "dpc_conflict_lo",  "dpc_conflict_hi"
)

# Five easily-distinguishable colors, one per shift level (0 = our default/
# main setting, shown in black; ColorBrewer Dark2 for the rest)
DELAY_LEVELS <- c("0d (default)", paste0("+", setdiff(SHIFTS, 0), "d"))
DELAY_COLORS <- setNames(
  c("grey40", "#1b9e77", "#d95f02", "#7570b3", "#e7298a"),
  DELAY_LEVELS
)

# =============================================================================
# Panel 1: efficacy vs DPC curves, one line per shift level
# =============================================================================
curve_d50_dat <- readRDS(here("data-processed", "DPC_fixed_efficacy_varied_d50.rds"))

make_efficacy_fn_shifted <- function(efficacy_col, shift) {
  force(efficacy_col)
  force(shift)
  function(dpc) {
    approx(x = curve_d50_dat$dpc, y = curve_d50_dat[[efficacy_col]],
           xout = pmax(dpc - shift, 0), rule = 2)$y
  }
}

dpc_seq <- seq(0, max(curve_d50_dat$dpc) + max(SHIFTS), by = 0.1)

eff_curves_df <- do.call(rbind, lapply(SHIFTS, function(s) {
  eff_fn <- make_efficacy_fn_shifted("efficacy", s)
  label  <- if (s == 0) "0d (default)" else paste0("+", s, "d")
  data.frame(
    dpc          = dpc_seq,
    efficacy     = eff_fn(dpc_seq),
    shift_label  = factor(label, levels = DELAY_LEVELS)
  )
}))

panel_curves <- ggplot(eff_curves_df, aes(x = dpc, y = efficacy, color = shift_label)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = DELAY_COLORS, name = "Efficacy delay") +
  scale_y_continuous(limits = c(0, NA), labels = scales::percent) +
  labs(x = "Days post-exposure (DPC)", y = "Efficacy") +
  theme_fig() +
  theme(legend.position = "bottom")

# =============================================================================
# Panel 2: HCW deaths averted (%) by shift, faceted by arm
# =============================================================================
particle_df_delay <- read.csv(here("output_figgen", "figure_supp_delay_particle_summary.csv"),
                              stringsAsFactors = FALSE)

averted_summ <- particle_df_delay %>%
  filter(arm %in% ARM_PANEL_ORDER, shift %in% SHIFTS) %>%
  group_by(scenario, arm, shift) %>%
  summarise(
    median = median(pct_hcw_deaths_averted, na.rm = TRUE),
    lo95   = quantile(pct_hcw_deaths_averted, 0.025, na.rm = TRUE),
    hi95   = quantile(pct_hcw_deaths_averted, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    shift_label = factor(ifelse(shift == 0, "0d (default)", paste0("+", shift, "d")),
                         levels = DELAY_LEVELS),
    arm_label   = factor(arm, levels = ARM_PANEL_ORDER)
  )

panel_bars <- ggplot(averted_summ, aes(x = arm_label, y = median, fill = shift_label)) +
  geom_col(position = position_dodge(0.75), width = 0.7) +
  geom_errorbar(aes(ymin = lo95, ymax = hi95),
                position = position_dodge(0.75), width = 0.2, color = "black", linewidth = 0.4) +
  scale_fill_manual(values = DELAY_COLORS, name = "Efficacy delay") +
  scale_y_continuous(limits = c(0, 100), labels = function(x) paste0(x, "%")) +
  labs(x = NULL, y = "HCW deaths averted (%)") +
  theme_fig() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "bottom")

# =============================================================================
# Combine and save
# =============================================================================
fig_delay <- (panel_curves | panel_bars) +
  plot_layout(widths = c(1, 1.3)) +
  plot_annotation(tag_levels = "a")

save_fig("figure_supp_delay", fig_delay, 16, 6)
message("figure_supp_delay saved.")