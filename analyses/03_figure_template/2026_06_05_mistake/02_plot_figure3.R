# =============================================================================
# 02_plot_figure3.R
# =============================================================================
source(here::here("analyses", "03_figure_template", "helper_functions_figure_1to4.R"))
OUT_DIR <- here("figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Coverage specs (must match simulation settings)
COVERAGE_SPECS <- list(
  full      = list(times = c(0, 90), values = c(1.00, 1.00)),
  ramp_high = list(times = c(0, 30, 60, 90), values = c(0.20, 0.47, 0.73, 1.00)),
  ramp_low  = list(times = c(0, 30, 60, 90), values = c(0.20, 0.30, 0.40, 0.50))
)

# =============================================================================
# Panels a, b, c: Coverage curve time series
# =============================================================================
make_coverage_plot <- function(cs) {
  spec <- COVERAGE_SPECS[[cs]]
  df   <- data.frame(time = spec$times, coverage = spec$values)
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
# Load simulation results
# =============================================================================
results_full      <- load_results("full")
results_ramp_high <- load_results("ramp_high")
results_ramp_low  <- load_results("ramp_low")

run_full <- build_run_df(results_full)
base_df  <- run_full %>%
  filter(arm == "baseline") %>%
  group_by(scenario, particle_id) %>%
  summarise(
    baseline_hcw_deaths = mean(n_hcw_deaths),
    baseline_days_lost  = mean(hcw_days_lost),
    .groups = "drop"
  )

build_obv_df <- function(results) {
  build_run_df(results) %>%
    filter(arm %in% OBV_EFFICACY_LEVELS) %>%
    group_by(scenario, particle_id, arm, coverage_scenario) %>%
    summarise(
      n_hcw_deaths  = mean(n_hcw_deaths),
      hcw_days_lost = mean(hcw_days_lost),
      .groups = "drop"
    ) %>%
    left_join(base_df, by = c("scenario", "particle_id")) %>%
    mutate(
      pct_hcw_deaths_averted = 100 * (baseline_hcw_deaths - n_hcw_deaths) / baseline_hcw_deaths,
      pct_days_lost_averted  = 100 * (baseline_days_lost  - hcw_days_lost) / baseline_days_lost
    )
}

obv_df <- bind_rows(
  build_obv_df(results_full),
  build_obv_df(results_ramp_high),
  build_obv_df(results_ramp_low)
) %>%
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

# Panels d, e, f — HCW deaths averted
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

# Panels g, h, i — HCW days lost averted
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