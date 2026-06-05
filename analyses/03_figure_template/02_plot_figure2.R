# =============================================================================
# 02_plot_figure2.R
# OBV efficacy comparison at full coverage -- post-hoc applied to baseline runs
# =============================================================================
source(here::here("analyses", "03_figure_template", "helper_functions_figure_1to4.R"))
OUT_DIR <- here("figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

results <- load_results()

# Build baseline + all efficacy arms at full coverage in a single pass
run_df <- do.call(rbind, c(
  list(build_run_df_obv(results, "baseline")),
  lapply(OBV_EFFICACY_LEVELS, function(eff) {
    build_run_df_obv(results, eff, "full")
  })
))

pdf <- make_particle_df(run_df) %>%
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
  light_col <- if (sc == "WestAfrica") "#fdd8a0" else "#b2e4d8"
  base_col  <- SCENARIO_COLORS[sc]
  dark_col  <- rgb(t(col2rgb(base_col) * 0.7), maxColorValue = 255)
  # 50–80%: gradient from light to exact scenario color; 90%: one shade darker
  obv_colors <- setNames(
    c(colorRampPalette(c(light_col, base_col))(4), dark_col),
    OBV_EFFICACY_LABELS
  )

  ggplot(summ_df, aes(x = arm_label, y = median, fill = arm_label)) +
    geom_col(aes(color = arm_label == "80%"), width = 0.8, alpha = 0.85, linewidth = 0.8) +
    scale_color_manual(values = c("TRUE" = "black", "FALSE" = NA), guide = "none") +
    geom_errorbar(
      aes(ymin = lo95, ymax = hi95),
      color = "black", width = 0.2, linewidth = 0.8
    ) +
    geom_point(color = "black", size = 2) +
    scale_fill_manual(values = obv_colors, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0)),
                       limits = c(0, 100)) +
    labs(x = "OBV efficacy", y = y_label,
         title = title,
         subtitle = sprintf("%s | Bar: median | Line: 95%% CI across posterior particles",
                            SCENARIO_LABELS[sc])) +
    theme_fig() +
    theme(panel.grid.major.x = element_blank())
}


# panels
fig2a <- make_bar_plot(make_summ(pdf, "pct_hcw_deaths_averted", "WestAfrica"),
                       "WestAfrica", "HCW deaths averted (%)",
                       "% HCW deaths averted by OBV efficacy")

fig2b <- make_bar_plot(make_summ(pdf, "pct_days_lost_averted", "WestAfrica"),
                       "WestAfrica", "HCW days lost averted (%)",
                       "% HCW days lost averted by OBV efficacy")

fig2c <- make_bar_plot(make_summ(pdf, "pct_hcw_deaths_averted", "DRC"),
                       "DRC", "HCW deaths averted (%)",
                       "% HCW deaths averted by OBV efficacy")

fig2d <- make_bar_plot(make_summ(pdf, "pct_days_lost_averted", "DRC"),
                       "DRC", "HCW days lost averted (%)",
                       "% HCW days lost averted by OBV efficacy")

ggsave(file.path(OUT_DIR, "figure_2_a.png"),
       fig2a,
       width = 7, height = 5, dpi = 150)

ggsave(file.path(OUT_DIR, "figure_2_b.png"),
       fig2b,
       width = 7, height = 5, dpi = 150)

ggsave(file.path(OUT_DIR, "figure_2_c.png"),
       fig2c,
       width = 7, height = 5, dpi = 150)

ggsave(file.path(OUT_DIR, "figure_2_d.png"),
       fig2d ,
       width = 7, height = 5, dpi = 150)

message("Figure 2 panels saved: a, b, c, d")

# Combine all panels
make_header <- function(label, angle = 0) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = label, fontface = "bold", size = 5,
             angle = angle) +
    theme_void()
}

strip_titles <- function(p) p + theme(plot.title = element_blank(),
                                       plot.subtitle = element_blank())

fig2_all <- (
  (make_header("West Africa") | make_header("DRC")) /
  (strip_titles(fig2a) | strip_titles(fig2c)) /
  (strip_titles(fig2b) | strip_titles(fig2d))
) +
  plot_layout(heights = c(0.08, 1, 1)) +
  plot_annotation(tag_levels = list(c("", "", "a", "c", "b", "d")))

ggsave(file.path(OUT_DIR, "figure_2_ALL.png"), fig2_all,
       width = 11, height = 8, dpi = 150, units = "in")

message("Figure 2 combined saved")
