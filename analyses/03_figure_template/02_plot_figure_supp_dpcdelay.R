# =============================================================================
# 02_plot_figure_supp_delay.R
# Visualise the DPC-delay (efficacy drop-off delay) sensitivity sweep,
# including shift=0 (our default/main setting) as a grey reference line/bar.
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

# Five easily-distinguishable colors, one per shift level (0 = our default/
# main setting, shown in grey; ColorBrewer Dark2 for the rest)
DELAY_LEVELS <- c("0d (default)", paste0("+", setdiff(SHIFTS, 0), "d"))
DELAY_COLORS <- setNames(
  c("grey40", "#1b9e77", "#d95f02", "#7570b3", "#e7298a"),
  DELAY_LEVELS
)

# =============================================================================
# Left panel: efficacy vs DPC curves, one line per shift level, with a very
# light ribbon spanning the optimistic (eighty_efficacy_hi) to pessimistic
# (eighty_efficacy_lo) bounds, shifted the same way as the median.
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
  mid_fn <- make_efficacy_fn_shifted("efficacy", s)
  lo_fn  <- make_efficacy_fn_shifted("eighty_efficacy_lo", s)  # pessimistic
  hi_fn  <- make_efficacy_fn_shifted("eighty_efficacy_hi", s)  # optimistic
  label  <- if (s == 0) "0d (default)" else paste0("+", s, "d")
  data.frame(
    dpc         = dpc_seq,
    mid         = mid_fn(dpc_seq),
    lo          = lo_fn(dpc_seq),
    hi          = hi_fn(dpc_seq),
    shift_label = factor(label, levels = DELAY_LEVELS)
  )
}))

panel_curves <- ggplot(eff_curves_df, aes(x = dpc, color = shift_label, fill = shift_label)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.10, color = NA) +
  geom_line(aes(y = mid), linewidth = 1.1) +
  scale_color_manual(values = DELAY_COLORS, name = NULL) +
  scale_fill_manual(values = DELAY_COLORS, name = NULL) +
  scale_y_continuous(limits = c(0, NA), labels = scales::percent) +
  labs(x = "Days post-exposure (DPC)", y = "Efficacy") +
  theme_fig() +
  theme(legend.position = "none")

# =============================================================================
# Right panel: HCW deaths averted (%) by arm, dodged/colored by shift.
#
# Arm order and tick layout match figure_supp_conflict_intensity: dpc_conflict
# (lo, mid, hi) then with_conflict (lo, mid, hi); line 1 of the x-axis ticks
# is the efficacy level, line 2 (scenario name, spanning each trio) is drawn
# separately below the axis in bold/larger text.
# =============================================================================
particle_df_delay <- read.csv(here("output_figgen", "figure_supp_delay_particle_summary.csv"),
                              stringsAsFactors = FALSE)

ARM_PANEL_ORDER <- c(
  "dpc_conflict_lo",  "dpc_conflict_mid",  "dpc_conflict_hi",
  "with_conflict_lo", "with_conflict_mid", "with_conflict_hi"
)

EFF_TICK_LABELS       <- c(lo = "Pessimistic", mid = "Central", hi = "Optimistic")
SCENARIO_GROUP_LABELS <- c(dpc_conflict = "Delayed dosing", with_conflict = "Delayed coverage + dosing")

x_tick_labels <- sapply(ARM_PANEL_ORDER, function(arm_name) {
  eff_key <- sub("^(with_conflict|dpc_conflict)_", "", arm_name)
  EFF_TICK_LABELS[[eff_key]]
})

GROUP_LABEL_Y <- -7  # below the y=0 axis line; needs coord_cartesian(clip="off")

group_label_df <- data.frame(
  x     = c(2, 5),
  label = c(SCENARIO_GROUP_LABELS[["dpc_conflict"]], SCENARIO_GROUP_LABELS[["with_conflict"]])
)

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
  geom_text(data = group_label_df, aes(x = x, y = GROUP_LABEL_Y, label = label),
            inherit.aes = FALSE, fontface = "bold", size = 4, vjust = 1) +
  scale_fill_manual(values = DELAY_COLORS, name = "Efficacy delay") +
  scale_x_discrete(labels = x_tick_labels) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0, 0.05))) +
  coord_cartesian(ylim = c(0, 100), clip = "off") +
  labs(x = NULL, y = "HCW deaths averted (%)") +
  theme_fig() +
  theme(
    legend.position       = c(0.5, 0.97),
    legend.justification  = c(0.5, 1),
    legend.direction      = "horizontal",
    legend.background     = element_blank(),
    legend.key            = element_blank(),
    legend.text           = element_text(size = 8),
    plot.margin           = margin(t = 5, r = 10, b = 32, l = 5)
  ) +
  guides(fill = guide_legend(nrow = 1))

# =============================================================================
# Combine and save
# =============================================================================
fig_delay <- (panel_curves | panel_bars) +
  plot_layout(widths = c(1, 1.3)) +
  plot_annotation(tag_levels = "a")

save_fig("figure_supp_delay", fig_delay, 16, 5)
message("figure_supp_delay saved.")