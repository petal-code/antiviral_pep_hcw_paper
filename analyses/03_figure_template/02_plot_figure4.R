# =============================================================================
# 02_plot_figure4.R
# Efficacy x coverage heatmap -- fully post-hoc, no separate simulation needed
# =============================================================================
source(here::here(
  "analyses",
  "03_figure_template",
  "helper_functions_figure_1to4.R"
))
OUT_DIR <- here("figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Build panels ----
# Five scalar coverage levels (not ramp curves)
COVERAGE_GRID <- c(0.10, 0.30, 0.50, 0.70, 0.90)

message("Loading baseline results...")
results <- load_results()

message("Applying post-hoc OBV across efficacy x coverage grid...")

baseline_rows <- build_run_df_obv(results, "baseline")

grid_rows <- do.call(rbind, lapply(names(OBV_EFFICACY_VALUES), function(eff_name) {
  eff <- OBV_EFFICACY_VALUES[[eff_name]]
  do.call(rbind, lapply(COVERAGE_GRID, function(cov) {
    cov_spec  <- list(times = c(0, 1), values = c(cov, cov))
    cov_label <- sprintf("cov%02d", round(cov * 100))

    do.call(rbind, lapply(results, function(x) {
      run_seed <- x$particle_id * 1000L + x$rep
      obv      <- apply_obv_posthoc(x$tdf, eff, cov_spec, seed = run_seed)

      # Use prevented_flag from apply_obv_posthoc to keep prevented individuals
      # consistent with the prevented count.
      days_lost <- compute_hcw_days_lost(
        x$tdf, x$duration,
        obv_received = obv$obv_received,
        prevented    = obv$prevented_flag
      )

      data.frame(
        scenario           = x$scenario,
        particle_id        = x$particle_id,
        rep                = x$rep,
        arm                = eff_name,
        coverage_scenario  = cov_label,
        obv_efficacy       = eff,
        obv_coverage       = cov,
        n_infections       = x$n_infections,
        n_hcw_deaths       = x$n_hcw_deaths - obv$prevented_hcw,
        counterfactual_hcw = x$n_hcw_deaths,
        prevented_hcw      = obv$prevented_hcw,
        hcw_days_lost      = days_lost,
        stringsAsFactors   = FALSE
      )
    }))
  }))
}))

base_particle <- baseline_rows %>%
  group_by(scenario, particle_id) %>%
  summarise(baseline_days_lost = mean(hcw_days_lost), .groups = "drop")

heatmap_df <- grid_rows %>%
  group_by(scenario, particle_id, obv_efficacy, obv_coverage) %>%
  summarise(
    prevented_hcw      = sum(prevented_hcw),
    counterfactual_hcw = sum(counterfactual_hcw),
    hcw_days_lost      = mean(hcw_days_lost),
    .groups = "drop"
  ) %>%
  left_join(base_particle, by = c("scenario", "particle_id")) %>%
  mutate(
    pct_hcw_deaths_averted = ifelse(
      counterfactual_hcw > 0,
      100 * prevented_hcw / counterfactual_hcw,
      NA_real_
    ),
    pct_days_lost_averted = ifelse(
      !is.na(baseline_days_lost) & baseline_days_lost > 0,
      100 * (baseline_days_lost - hcw_days_lost) / baseline_days_lost,
      NA_real_
    )
  ) %>%
  group_by(scenario, obv_efficacy, obv_coverage) %>%
  summarise(
    median_deaths_averted    = median(pct_hcw_deaths_averted, na.rm = TRUE),
    median_days_lost_averted = median(pct_days_lost_averted,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    coverage_label = factor(paste0(round(obv_coverage * 100), "%"), levels = paste0(c(10, 30, 50, 70, 90), "%")),
    efficacy_label = factor(paste0(round(obv_efficacy * 100), "%"), levels = paste0(c(50, 60, 70, 80, 90), "%"))
  )

make_heatmap <- function(sc, metric, fill_label, palette = "YlOrRd") {
  df <- filter(heatmap_df, scenario == sc)

  ggplot(df, aes(x = coverage_label, y = efficacy_label, fill = .data[[metric]])) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(
      aes(label = sprintf("%.0f%%", .data[[metric]])),
      size = 4, fontface = "bold", color = "grey20"
    ) +
    scale_fill_distiller(
      palette   = palette,
      direction = 1,
      name      = fill_label,
      limits    = c(0, 100),
      labels    = function(x) paste0(x, "%"),
      na.value  = "grey90"
    ) +
    labs(x = "OBV coverage", y = "OBV efficacy") +
    theme_fig(base_size = 13) +
    theme(legend.position = "right", panel.grid = element_blank(), axis.ticks = element_blank())
}

# Combine panels ----
make_header <- function(label) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = label, fontface = "bold", size = 5) +
    theme_void()
}

fig4a <- make_heatmap("WestAfrica", "median_deaths_averted",    "Deaths averted",   palette = "YlOrRd")
fig4b <- make_heatmap("DRC",        "median_deaths_averted",    "Deaths averted",   palette = "YlOrRd")
fig4c <- make_heatmap("WestAfrica", "median_days_lost_averted", "Days lost averted", palette = "YlGnBu")
fig4d <- make_heatmap("DRC",        "median_days_lost_averted", "Days lost averted", palette = "YlGnBu")

fig4_all <- ((make_header("West Africa archetype") | make_header("DRC archetype")) /
               ((fig4a + fig4b + fig4c + fig4d) + plot_layout(nrow = 2, axis_titles = "collect"))) +
  plot_layout(guides = "collect", heights = c(0.08, 1)) +
  plot_annotation(tag_levels = list(c("", "", "a ", "b ", "c ", "d ")))

ggsave(
  file.path(OUT_DIR, "figure_4.png"),
  fig4_all, width = 14, height = 10, dpi = 150, units = "in"
)

message("Figure 4 saved")
