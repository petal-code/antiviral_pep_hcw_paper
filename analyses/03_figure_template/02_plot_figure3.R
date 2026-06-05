# =============================================================================
# 02_plot_figure3.R
# Coverage scenario comparison -- post-hoc applied to baseline runs
# =============================================================================
source(here::here("analyses", "03_figure_template", "helper_functions_figure_1to4.R"))
OUT_DIR <- here("figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

results <- load_results()

# =============================================================================
# Panels a, b, c: Coverage curve visualisations
# =============================================================================
make_coverage_plot <- function(cs) {
  spec  <- COVERAGE_SPECS[[cs]]
  df    <- data.frame(time = spec$times, coverage = spec$values)
  color <- COVERAGE_COLORS[cs]

  ggplot(df, aes(x = time, y = coverage)) +
    geom_line(color = color, linewidth = 1.2) +
    geom_point(color = color, size = 3) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    scale_x_continuous(breaks = c(0, 30, 60, 90)) +
    labs(x = "Days since outbreak start",
         y = "OBV coverage",
         title = COVERAGE_LABELS[match(cs, COVERAGE_LEVELS)],
         subtitle = "OBV coverage over time") +
    theme_fig()
}

ggsave(file.path(OUT_DIR, "figure_3_a.png"), make_coverage_plot("full"),
       width = 5, height = 4, dpi = 150)
ggsave(file.path(OUT_DIR, "figure_3_b.png"), make_coverage_plot("ramp_high"),
       width = 5, height = 4, dpi = 150)
ggsave(file.path(OUT_DIR, "figure_3_c.png"), make_coverage_plot("ramp_low"),
       width = 5, height = 4, dpi = 150)

# =============================================================================
# Build post-hoc OBV data: all efficacy x coverage combinations in one pass
# =============================================================================
message("Applying post-hoc OBV across efficacy x coverage grid...")

# Baseline rows needed once for the denominator
baseline_rows <- build_run_df_obv(results, "baseline")

obv_rows <- do.call(rbind, lapply(OBV_EFFICACY_LEVELS, function(eff) {
  do.call(rbind, lapply(COVERAGE_LEVELS, function(cov) {
    build_run_df_obv(results, eff, cov)
  }))
}))

all_rows <- rbind(baseline_rows, obv_rows)
pdf      <- make_particle_df(all_rows)

obv_df <- pdf %>%
  filter(arm %in% OBV_EFFICACY_LEVELS) %>%
  mutate(
    arm_label      = factor(OBV_EFFICACY_LABELS[match(arm, OBV_EFFICACY_LEVELS)],
                            levels = OBV_EFFICACY_LABELS),
    coverage_label = factor(COVERAGE_LABELS[match(coverage_scenario, COVERAGE_LEVELS)],
                            levels = COVERAGE_LABELS),
    scenario_label = factor(SCENARIO_LABELS[scenario], levels = SCENARIO_LABELS)
  )

# =============================================================================
# Panels d-f (HCW deaths averted) and g-i (HCW days lost averted)
# x = efficacy, colour = scenario, one panel per coverage scenario
# =============================================================================
sc_colors <- setNames(SCENARIO_COLORS, SCENARIO_LABELS)

make_box_plot <- function(cs, metric, y_label, title) {
  df <- obv_df %>%
    filter(coverage_scenario == cs, !is.na(.data[[metric]]))

  ggplot(df, aes(x = arm_label, y = .data[[metric]],
                 fill = scenario_label, color = scenario_label)) +
    geom_boxplot(outlier.size = 0.5, width = 0.6,
                 alpha = 0.6, position = position_dodge(0.75)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    scale_fill_manual(values  = sc_colors, name = NULL) +
    scale_color_manual(values = sc_colors, name = NULL) +
    labs(x = "OBV efficacy", y = y_label,
         title = title,
         subtitle = sprintf("%s | Boxplot across posterior particles",
                            COVERAGE_LABELS[match(cs, COVERAGE_LEVELS)])) +
    theme_fig() +
    theme(panel.grid.major.x = element_blank())
}

# Panels d, e, f -- HCW deaths averted
for (i in seq_along(COVERAGE_LEVELS)) {
  cs    <- COVERAGE_LEVELS[i]
  label <- letters[3 + i]   # d, e, f
  ggsave(
    file.path(OUT_DIR, sprintf("figure_3_%s.png", label)),
    make_box_plot(cs, "pct_hcw_deaths_averted",
                  "HCW deaths averted (%)",
                  "% HCW deaths averted"),
    width = 7, height = 5, dpi = 150
  )
}

# Panels g, h, i -- HCW days lost averted
for (i in seq_along(COVERAGE_LEVELS)) {
  cs    <- COVERAGE_LEVELS[i]
  label <- letters[6 + i]   # g, h, i
  ggsave(
    file.path(OUT_DIR, sprintf("figure_3_%s.png", label)),
    make_box_plot(cs, "pct_days_lost_averted",
                  "HCW days lost averted (%)",
                  "% HCW days lost averted"),
    width = 7, height = 5, dpi = 150
  )
}

message("Figure 3 panels saved: a-c (coverage curves), d-f (deaths averted), g-i (days lost averted)")
