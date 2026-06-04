# =============================================================================
# 02_plot_figure1.R
# =============================================================================
source(here::here("analyses", "03_figure_template", "helper_functions_figure_1to4.R"))
OUT_DIR <- here("figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

results <- load_results("full")

# =============================================================================
# Panel A/C: Total infections histogram (baseline only)
# =============================================================================
run_df <- flatten_to_df(results)
hist_df <- run_df %>%
  filter(arm == "baseline") %>%
  group_by(scenario, particle_id) %>%
  summarise(n_infections = mean(n_infections), .groups = "drop")

make_hist <- function(sc) {
  df    <- filter(hist_df, scenario == sc)
  color <- SCENARIO_COLORS[sc]
  ggplot(df, aes(x = n_infections)) +
    geom_histogram(bins = 40, fill = color, alpha = 0.75, color = "white") +
    labs(x = "Total infections", y = "Posterior particles",
         title = SCENARIO_LABELS[sc],
         subtitle = "Baseline — distribution across posterior particles") +
    theme_fig()
}

# =============================================================================
# Panel B/D: Cumulative HCW deaths time series — baseline vs OBV 80%
# =============================================================================
ts_df <- build_weekly_ts(results, metric = "hcw_deaths",
                         arms = c("baseline", "obv_80"))
ts_colors <- c(baseline = "#555555", obv_80 = "#2166ac")
ts_labels <- c(baseline = "Baseline", obv_80 = "OBV 80%")

make_ts <- function(sc) {
  arms      <- c("baseline", "obv_80")
  ts_colors <- get_arm_colors(sc, arms)
  ts_labels <- ARM_LABELS[arms]
  
  df <- filter(ts_df, scenario == sc) %>%
    mutate(arm = factor(arm, levels = arms))
  ggplot(df, aes(x = week, color = arm, fill = arm)) +
    geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.12, color = NA) +
    geom_ribbon(aes(ymin = q25,  ymax = q75),  alpha = 0.25, color = NA) +
    geom_line(aes(y = q50), linewidth = 1.0) +
    scale_color_manual(values = ts_colors, labels = ts_labels, name = NULL) +
    scale_fill_manual( values = ts_colors, labels = ts_labels, name = NULL) +
    labs(x = "Days since outbreak start",
         y = "Cumulative HCW deaths",
         title = SCENARIO_LABELS[sc],
         subtitle = "Baseline vs OBV 80% | Line: median | Bands: 50% / 95% CI") +
    theme_fig()
}

# =============================================================================
# Save each panel independently
# =============================================================================
ggsave(file.path(OUT_DIR, "figure_1_a.png"), make_hist("WestAfrica"),
       width = 7, height = 5, dpi = 150)
ggsave(file.path(OUT_DIR, "figure_1_b.png"), make_ts("WestAfrica"),
       width = 7, height = 5, dpi = 150)
ggsave(file.path(OUT_DIR, "figure_1_c.png"), make_hist("DRC"),
       width = 7, height = 5, dpi = 150)
ggsave(file.path(OUT_DIR, "figure_1_d.png"), make_ts("DRC"),
       width = 7, height = 5, dpi = 150)

message("Figure 1 panels saved: a, b, c, d")