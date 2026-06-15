# =============================================================================
# 02_plot_figure4.R
# Efficacy x coverage heatmap
# Reads pre-computed CSV from output_figgen/figure_4_run_summary.csv
# Run 02_extract_figure4.R first.
# =============================================================================
source(here::here("analyses", "03_figure_template", "helper_functions_figure_1to4.R"))
OUT_DIR <- here("figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

run_df <- read.csv(here("output_figgen", "figure_4_run_summary.csv"),
                   stringsAsFactors = FALSE)

heatmap_df <- run_df %>%
  group_by(scenario, particle_id, arm, obv_efficacy, obv_coverage) %>%
  summarise(
    n_hcw_deaths       = mean(n_hcw_deaths),
    hcw_days_lost      = sum(hcw_days_lost),
    prevented_hcw      = sum(prevented_hcw),
    counterfactual_hcw = sum(counterfactual_hcw),
    baseline_days_lost = sum(baseline_days_lost),
    .groups = "drop"
  ) %>%
  mutate(
    pct_hcw_deaths_averted = ifelse(
      counterfactual_hcw > 0, 100 * prevented_hcw / counterfactual_hcw, NA_real_),
    pct_days_lost_averted = ifelse(
      !is.na(baseline_days_lost) & baseline_days_lost > 0,
      100 * (baseline_days_lost - hcw_days_lost) / baseline_days_lost, NA_real_)
  ) %>%
  group_by(scenario, obv_efficacy, obv_coverage) %>%
  summarise(
    median_deaths_averted    = median(pct_hcw_deaths_averted, na.rm = TRUE),
    median_days_lost_averted = median(pct_days_lost_averted,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    coverage_label = factor(paste0(round(obv_coverage * 100), "%"),
                            levels = paste0(c(10, 30, 50, 70, 90), "%")),
    efficacy_label = factor(paste0(round(obv_efficacy * 100), "%"),
                            levels = paste0(c(50, 60, 70, 80, 90), "%"))
  )

save_figure_data(heatmap_df, "figure_4_heatmap_df.csv")

make_heatmap <- function(sc, metric, fill_label, subtitle = NULL) {
  df       <- filter(heatmap_df, scenario == sc)
  sc_color <- SCENARIO_COLORS[[sc]]
  ggplot(df, aes(x = coverage_label, y = efficacy_label, fill = .data[[metric]])) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.1f%%", .data[[metric]])),
              size = 4, fontface = "bold", color = "grey20") +
    scale_fill_gradient(low = "white", high = sc_color, name = fill_label,
                        limits = c(0, 100), labels = function(x) paste0(x, "%"),
                        na.value = "grey90") +
    labs(x = "Antiviral coverage", y = "Antiviral efficacy", subtitle = subtitle) +
    theme_fig(base_size = 13) +
    theme(legend.position = "none", panel.grid = element_blank(),
          axis.line = element_blank())
}

make_header <- function(label) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = label, fontface = "bold", size = 4.5) +
    theme_void()
}

fig4a <- make_heatmap("WestAfrica", "median_deaths_averted",    "Deaths averted",    subtitle = "% HCW deaths averted")
fig4b <- make_heatmap("DRC",        "median_deaths_averted",    "Deaths averted",    subtitle = "% HCW deaths averted")
fig4c <- make_heatmap("WestAfrica", "median_days_lost_averted", "Days lost averted", subtitle = "% HCW days lost averted")
fig4d <- make_heatmap("DRC",        "median_days_lost_averted", "Days lost averted", subtitle = "% HCW days lost averted")

fig4_all <- (
  (make_header("West Africa archetype") | make_header("DRC archetype")) /
    ((fig4a + fig4b + fig4c + fig4d) + plot_layout(nrow = 2, axis_titles = "collect"))
) +
  plot_layout(guides = "collect", heights = c(0.1, 1)) +
  plot_annotation(tag_levels = list(c("", "", "a ", "b ", "c ", "d ")))

save_fig <- function(filename_base, plot, width, height) {
  ggsave(file.path(OUT_DIR, paste0(filename_base, ".png")),
         plot, width = width, height = height, dpi = 400, units = "in")
  ggsave(file.path(OUT_DIR, paste0(filename_base, ".pdf")),
         plot, width = width, height = height, units = "in")
}

save_fig("figure_4", fig4_all, 10, 6.5)

# Deaths-only variant
fig4_deaths <- (
  (make_header("West Africa archetype") | make_header("DRC archetype")) /
    ((fig4a | fig4b) + plot_layout(axis_titles = "collect"))
) +
  plot_layout(heights = c(0.1, 1)) +
  plot_annotation(tag_levels = list(c("", "", "a ", "b ")))
save_fig("figure_4_deaths-averted", fig4_deaths, 10, 4)

# Days-lost-only variant
fig4_days <- (
  (make_header("West Africa archetype") | make_header("DRC archetype")) /
    ((fig4c | fig4d) + plot_layout(axis_titles = "collect"))
) +
  plot_layout(heights = c(0.1, 1)) +
  plot_annotation(tag_levels = list(c("", "", "a ", "b ")))
save_fig("figure_4_days-averted", fig4_days, 10, 4)

message("Figure 4 saved")


############### aggregating number for the paper
run_df4 <- read.csv(here("output_figgen", "figure_4_run_summary.csv"))

particle_df4 <- run_df4 %>%
  group_by(scenario, particle_id, arm, obv_efficacy, obv_coverage) %>%
  summarise(
    prevented_hcw      = sum(prevented_hcw),
    counterfactual_hcw = sum(counterfactual_hcw),
    .groups = "drop"
  ) %>%
  mutate(
    pct_hcw_deaths_averted = ifelse(
      counterfactual_hcw > 0, 100 * prevented_hcw / counterfactual_hcw, NA_real_)
  )

particle_df4 %>%
  filter(obv_efficacy == 0.8, obv_coverage %in% c(0.1, 0.9)) %>%
  group_by(scenario, obv_coverage) %>%
  summarise(
    median_pct = median(pct_hcw_deaths_averted, na.rm = TRUE),
    lo_pct     = quantile(pct_hcw_deaths_averted, 0.025, na.rm = TRUE),
    hi_pct     = quantile(pct_hcw_deaths_averted, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~round(., 1))) %>%
  as.data.frame() %>% print()

particle_df4 %>%
  filter(obv_coverage == 0.9, obv_efficacy %in% c(0.5, 0.9)) %>%
  group_by(scenario, obv_efficacy) %>%
  summarise(
    median_pct = median(pct_hcw_deaths_averted, na.rm = TRUE),
    lo_pct     = quantile(pct_hcw_deaths_averted, 0.025, na.rm = TRUE),
    hi_pct     = quantile(pct_hcw_deaths_averted, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~round(., 1))) %>%
  as.data.frame() %>% print()

# Days lost averted by efficacy x coverage grid
particle_df4 %>%
  group_by(scenario, particle_id, arm, obv_efficacy, obv_coverage) %>%
  summarise(
    hcw_days_lost      = sum(hcw_days_lost,      na.rm = TRUE),
    baseline_days_lost = sum(baseline_days_lost, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    pct_days_lost_averted = ifelse(
      baseline_days_lost > 0,
      100 * (baseline_days_lost - hcw_days_lost) / baseline_days_lost,
      NA_real_
    )
  ) %>%
  # Key comparisons matching figure 4 text structure
  filter(
    (obv_efficacy == 0.8 & obv_coverage %in% c(0.1, 0.9)) |
      (obv_coverage == 0.9 & obv_efficacy %in% c(0.5, 0.9))
  ) %>%
  group_by(scenario, obv_efficacy, obv_coverage) %>%
  summarise(
    med = median(pct_days_lost_averted, na.rm = TRUE),
    lo  = quantile(pct_days_lost_averted, 0.025, na.rm = TRUE),
    hi  = quantile(pct_days_lost_averted, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~round(.x, 1))) %>%
  as.data.frame() %>% print()