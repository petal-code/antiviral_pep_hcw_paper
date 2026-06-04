# =============================================================================
# 02_plot_figure4.R
# =============================================================================
source(here::here("analyses", "03_figure_template", "helper_functions_figure_1to4.R"))

OUT_DIR  <- here("figures")
GRID_DIR <- here("outputs", "simulation_fig4")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Load grid results
# =============================================================================
message("Loading grid results...")
files <- list.files(GRID_DIR, pattern = "\\.rds$", full.names = TRUE)
if (length(files) == 0) stop("No RDS files found in: ", GRID_DIR)
message(sprintf("  %d files found", length(files)))
results <- lapply(files, readRDS)

# =============================================================================
# Flatten to data frame
# =============================================================================
grid_df <- do.call(rbind, lapply(results, function(x) {
  days_lost <- compute_hcw_days_lost(x$tdf, x$duration)
  data.frame(
    scenario      = x$scenario,
    particle_id   = x$particle_id,
    obv_efficacy  = x$obv_efficacy,
    obv_coverage  = x$obv_coverage,
    rep           = x$rep,
    n_hcw_deaths  = x$n_hcw_deaths,
    hcw_days_lost = days_lost,
    stringsAsFactors = FALSE
  )
}))

# =============================================================================
# Attach baseline from simulation_fig1to3/full
# =============================================================================
message("Loading baseline...")
base_df <- build_run_df(load_results("full")) %>%
  filter(arm == "baseline") %>%
  group_by(scenario, particle_id) %>%
  summarise(
    baseline_hcw_deaths = mean(n_hcw_deaths),
    baseline_days_lost  = mean(hcw_days_lost),
    .groups = "drop"
  )

# =============================================================================
# Aggregate: mean over reps per particle, then median over particles
# =============================================================================
heatmap_df <- grid_df %>%
  group_by(scenario, particle_id, obv_efficacy, obv_coverage) %>%
  summarise(
    n_hcw_deaths  = mean(n_hcw_deaths),
    hcw_days_lost = mean(hcw_days_lost),
    .groups = "drop"
  ) %>%
  left_join(base_df, by = c("scenario", "particle_id")) %>%
  mutate(
    pct_hcw_deaths_averted = ifelse(
      baseline_hcw_deaths > 0,
      100 * (baseline_hcw_deaths - n_hcw_deaths) / baseline_hcw_deaths,
      NA_real_
    ),
    pct_days_lost_averted = ifelse(
      baseline_days_lost > 0,
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
    coverage_label = factor(paste0(round(obv_coverage * 100), "%"),
                            levels = paste0(c(10,30,50,70,90), "%")),
    efficacy_label = factor(paste0(round(obv_efficacy * 100), "%"),
                            levels = paste0(c(50,60,70,80,90), "%"))
  )

# =============================================================================
# Heatmap plot function
# =============================================================================
make_heatmap <- function(sc, metric, fill_label, title) {
  df <- filter(heatmap_df, scenario == sc)
  
  ggplot(df, aes(x = coverage_label, y = efficacy_label,
                 fill = .data[[metric]])) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.0f%%", .data[[metric]])),
              size = 4, fontface = "bold", color = "grey20") +
    scale_fill_distiller(
      palette   = "YlOrRd",
      direction = 1,
      name      = fill_label,
      limits    = c(0, 100),
      na.value  = "grey90"
    ) +
    labs(x = "OBV coverage", y = "OBV efficacy",
         title    = title,
         subtitle = sprintf("%s | Median across posterior particles",
                            SCENARIO_LABELS[sc]),
         caption  = "Stochastic branching process model (fiber)") +
    theme_fig(base_size = 13) +
    theme(legend.position = "right",
          panel.grid      = element_blank(),
          axis.ticks      = element_blank())
}

# =============================================================================
# Save four independent panels
# =============================================================================
panels <- list(
  list(sc = "WestAfrica", metric = "median_deaths_averted",
       fill_label = "Deaths averted (%)",
       title = "% HCW deaths averted — West Africa",
       file  = "figure_4_a.png"),
  list(sc = "WestAfrica", metric = "median_days_lost_averted",
       fill_label = "Days lost averted (%)",
       title = "% HCW days lost averted — West Africa",
       file  = "figure_4_b.png"),
  list(sc = "DRC", metric = "median_deaths_averted",
       fill_label = "Deaths averted (%)",
       title = "% HCW deaths averted — DRC",
       file  = "figure_4_c.png"),
  list(sc = "DRC", metric = "median_days_lost_averted",
       fill_label = "Days lost averted (%)",
       title = "% HCW days lost averted — DRC",
       file  = "figure_4_d.png")
)

for (panel in panels) {
  ggsave(file.path(OUT_DIR, panel$file),
         make_heatmap(panel$sc, panel$metric, panel$fill_label, panel$title),
         width = 7, height = 6, dpi = 150)
  message(sprintf("Saved: %s", panel$file))
}
