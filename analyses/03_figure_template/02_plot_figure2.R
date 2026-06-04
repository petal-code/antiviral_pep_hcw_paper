# =============================================================================
# 02_plot_figure2.R
# =============================================================================
source(here::here("analyses", "03_figure_template", "helper_functions_figure_1to4.R"))
OUT_DIR <- here("figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

results <- load_results("full")
run_df  <- build_run_df(results)
pdf     <- make_particle_df(run_df)

obv_df <- pdf %>%
  filter(arm %in% OBV_EFFICACY_LEVELS) %>%
  mutate(
    arm_label      = factor(OBV_EFFICACY_LABELS[match(arm, OBV_EFFICACY_LEVELS)],
                            levels = OBV_EFFICACY_LABELS),
    scenario_label = SCENARIO_LABELS[scenario]
  )

make_summ <- function(df, metric, sc) {
  df %>%
    filter(scenario == sc) %>%
    group_by(arm_label) %>%
    summarise(
      median = median(.data[[metric]], na.rm = TRUE),
      lo95   = quantile(.data[[metric]], 0.025, na.rm = TRUE),
      hi95   = quantile(.data[[metric]], 0.975, na.rm = TRUE),
      .groups = "drop"
    )
}

make_bar_plot <- function(summ_df, sc, y_label, title) {
  obv_colors <- setNames(ARM_COLORS_OBV, OBV_EFFICACY_LABELS)
  
  ggplot(summ_df, aes(x = arm_label, y = median, fill = arm_label)) +
    geom_col(width = 0.6, alpha = 0.85) +
    geom_errorbar(
      aes(ymin = lo95, ymax = hi95, color = arm_label),
      width = 0, linewidth = 0.8
    ) +
    geom_point(aes(color = arm_label), size = 2) +
    scale_fill_manual(values  = obv_colors, guide = "none") +
    scale_color_manual(values = obv_colors, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(x = "OBV efficacy", y = y_label,
         title = title,
         subtitle = sprintf("%s | Bar: median | Line: 95%% CI across posterior particles",
                            SCENARIO_LABELS[sc])) +
    theme_fig() +
    theme(panel.grid.major.x = element_blank())
}

ggsave(file.path(OUT_DIR, "figure_2_a.png"),
       make_bar_plot(make_summ(obv_df, "pct_hcw_deaths_averted", "WestAfrica"),
                     "WestAfrica", "HCW deaths averted (%)",
                     "% HCW deaths averted by OBV efficacy"),
       width = 7, height = 5, dpi = 150)

ggsave(file.path(OUT_DIR, "figure_2_b.png"),
       make_bar_plot(make_summ(obv_df, "pct_days_lost_averted", "WestAfrica"),
                     "WestAfrica", "HCW days lost averted (%)",
                     "% HCW days lost averted by OBV efficacy"),
       width = 7, height = 5, dpi = 150)

ggsave(file.path(OUT_DIR, "figure_2_c.png"),
       make_bar_plot(make_summ(obv_df, "pct_hcw_deaths_averted", "DRC"),
                     "DRC", "HCW deaths averted (%)",
                     "% HCW deaths averted by OBV efficacy"),
       width = 7, height = 5, dpi = 150)

ggsave(file.path(OUT_DIR, "figure_2_d.png"),
       make_bar_plot(make_summ(obv_df, "pct_days_lost_averted", "DRC"),
                     "DRC", "HCW days lost averted (%)",
                     "% HCW days lost averted by OBV efficacy"),
       width = 7, height = 5, dpi = 150)

message("Figure 2 panels saved: a, b, c, d")