# =============================================================================
# 02_plot_figure2.R
# OBV efficacy comparison at full coverage (50-90%)
# Reads pre-computed CSV from output_figgen/figure_2_run_summary.csv
# Run 02_extract_figure2.R first.
# =============================================================================
source(here::here("analyses", "03_figure_template", "helper_functions_figure_1to4.R"))
OUT_DIR <- here("figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FIG2_EFFICACY_LEVELS <- c("obv_50", "obv_60", "obv_70", "obv_80", "obv_90")
FIG2_EFFICACY_LABELS <- c("50%", "60%", "70%", "80%", "90%")

run_df <- read.csv(here("output_figgen", "figure_2_run_summary.csv"),
                   stringsAsFactors = FALSE)

pdf <- make_particle_df(run_df) %>%
  filter(arm %in% FIG2_EFFICACY_LEVELS) %>%
  mutate(
    arm_label      = factor(FIG2_EFFICACY_LABELS[match(arm, FIG2_EFFICACY_LEVELS)],
                            levels = FIG2_EFFICACY_LABELS),
    scenario_label = SCENARIO_LABELS[scenario]
  )

save_figure_data(pdf, "figure_2_particle_df.csv")

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

make_bar_plot <- function(summ_df, sc, y_label) {
  light_col  <- if (sc == "WestAfrica") "#fdd8a0" else "#b2e4d8"
  base_col   <- SCENARIO_COLORS[sc]
  dark_col   <- rgb(t(col2rgb(base_col) * 0.7), maxColorValue = 255)
  obv_colors <- setNames(
    colorRampPalette(c(light_col, dark_col))(length(FIG2_EFFICACY_LABELS)),
    FIG2_EFFICACY_LABELS
  )
  ggplot(summ_df, aes(x = arm_label, y = median, fill = arm_label)) +
    geom_col(aes(color = arm_label == "80%"), width = 0.8, alpha = 0.85, linewidth = 0.8) +
    scale_color_manual(values = c("TRUE" = "black", "FALSE" = NA), guide = "none") +
    geom_errorbar(aes(ymin = lo95, ymax = hi95), color = "black", width = 0.2, linewidth = 0.8) +
    geom_point(color = "black", size = 2) +
    scale_fill_manual(values = obv_colors, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0)), limits = c(0, 100),
                       labels = function(x) paste0(x, "%")) +
    labs(x = "Antiviral efficacy", y = y_label) +
    theme_fig()
}

make_header <- function(label) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = label, fontface = "bold", size = 4.5) +
    theme_void()
}

save_fig <- function(filename_base, plot, width, height) {
  ggsave(file.path(OUT_DIR, paste0(filename_base, ".png")),
         plot, width = width, height = height, dpi = 400, units = "in")
  ggsave(file.path(OUT_DIR, paste0(filename_base, ".pdf")),
         plot, width = width, height = height, units = "in")
}

fig2a <- make_bar_plot(make_summ(pdf, "pct_hcw_deaths_averted", "WestAfrica"), "WestAfrica", "HCW deaths averted")
fig2b <- make_bar_plot(make_summ(pdf, "pct_hcw_deaths_averted", "DRC"),        "DRC",        "HCW deaths averted")
fig2c <- make_bar_plot(make_summ(pdf, "pct_days_lost_averted",  "WestAfrica"), "WestAfrica", "HCW days lost averted")
fig2d <- make_bar_plot(make_summ(pdf, "pct_days_lost_averted",  "DRC"),        "DRC",        "HCW days lost averted")

fig2_all <- (
  (make_header("West Africa archetype") | make_header("DRC archetype")) /
    ((fig2a | fig2b) + plot_layout(axis_titles = "collect")) /
    ((fig2c | fig2d) + plot_layout(axis_titles = "collect"))
) +
  plot_layout(heights = c(0.08, 1, 1)) +
  plot_annotation(tag_levels = list(c("", "", "a ", "b ", "c ", "d ")))

save_fig("figure_2", fig2_all, 10, 6.5)
message("Figure 2 saved")