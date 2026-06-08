# =============================================================================
# 02_plot_figure4.R
# Efficacy x coverage heatmap -- fully post-hoc, no separate simulation needed
# =============================================================================
source(here::here("analyses", "03_figure_template", "helper_functions_figure_1to4.R"))
OUT_DIR <- here("figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Coverage grid for figure 4: five scalar levels (not the ramp curves)
COVERAGE_GRID <- c(0.10, 0.30, 0.50, 0.70, 0.90)

# =============================================================================
# Load baseline results
# =============================================================================
message("Loading baseline results...")
results <- load_results()

# =============================================================================
# Build efficacy x coverage grid post-hoc
#
# For figure 4 the coverage axis uses fixed scalar values (not ramp curves),
# so we construct scalar coverage specs on the fly.
# =============================================================================
message("Applying post-hoc OBV across efficacy x coverage grid...")

efficacy_grid <- OBV_EFFICACY_VALUES   # named vector: obv_50..obv_90

# Baseline rows
baseline_rows <- build_run_df_obv(results, "baseline")

grid_rows <- do.call(rbind, lapply(names(efficacy_grid), function(eff_name) {
  eff <- efficacy_grid[[eff_name]]
  do.call(rbind, lapply(COVERAGE_GRID, function(cov) {
    # Build a flat (scalar) coverage spec for this coverage level
    cov_spec <- list(times = c(0, 1), values = c(cov, cov))
    cov_label <- sprintf("cov%02d", round(cov * 100))
    
    do.call(rbind, lapply(results, function(x) {
      run_seed  <- x$particle_id * 1000L + x$rep
      obv       <- apply_obv_posthoc(x$tdf, eff, cov_spec, seed = run_seed)
      
      # Use prevented_flag returned directly from apply_obv_posthoc so that
      # the prevented individuals are consistent with the prevented count.
      days_lost <- compute_hcw_days_lost(x$tdf, x$duration,
                                         obv_received = obv$obv_received,
                                         prevented    = obv$prevented_flag)
      
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

# =============================================================================
# Aggregate to particle level (burden-weighted % averted)
# =============================================================================
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
    coverage_label = factor(paste0(round(obv_coverage * 100), "%"),
                            levels = paste0(c(10, 30, 50, 70, 90), "%")),
    efficacy_label = factor(paste0(round(obv_efficacy * 100), "%"),
                            levels = paste0(c(50, 60, 70, 80, 90), "%"))
  )

# =============================================================================
# Heatmap plot function
# =============================================================================
make_heatmap <- function(sc, metric, fill_label, title, palette = "YlOrRd") {
  df <- filter(heatmap_df, scenario == sc)

  ggplot(df, aes(x = coverage_label, y = efficacy_label,
                 fill = .data[[metric]])) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.0f%%", .data[[metric]])),
              size = 4, fontface = "bold", color = "grey20") +
    scale_fill_distiller(
      palette   = palette,
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
# Save four panels
# =============================================================================
fig4a <- make_heatmap("WestAfrica", "median_deaths_averted",
                      "Deaths averted (%)", "% HCW deaths averted -- West Africa",
                      palette = "YlOrRd")
fig4b <- make_heatmap("WestAfrica", "median_days_lost_averted",
                      "Days lost averted (%)", "% HCW days lost averted -- West Africa",
                      palette = "YlGnBu")
fig4c <- make_heatmap("DRC", "median_deaths_averted",
                      "Deaths averted (%)", "% HCW deaths averted -- DRC",
                      palette = "YlOrRd")
fig4d <- make_heatmap("DRC", "median_days_lost_averted",
                      "Days lost averted (%)", "% HCW days lost averted -- DRC",
                      palette = "YlGnBu")

ggsave(file.path(OUT_DIR, "figure_4_a.png"), fig4a, width = 7, height = 6, dpi = 150)
ggsave(file.path(OUT_DIR, "figure_4_b.png"), fig4b, width = 7, height = 6, dpi = 150)
ggsave(file.path(OUT_DIR, "figure_4_c.png"), fig4c, width = 7, height = 6, dpi = 150)
ggsave(file.path(OUT_DIR, "figure_4_d.png"), fig4d, width = 7, height = 6, dpi = 150)

message("Figure 4 individual panels saved")

# =============================================================================
# Composite Figure 4
# =============================================================================
make_header <- function(label) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = label, fontface = "bold", size = 5) +
    theme_void()
}

strip_titles <- function(p) p + theme(plot.title    = element_blank(),
                                       plot.subtitle = element_blank(),
                                       plot.caption  = element_blank())

fig4_all <- (
  (make_header("West Africa") | make_header("DRC")) /
  (strip_titles(fig4a)        | strip_titles(fig4c)) /
  (strip_titles(fig4b)        | strip_titles(fig4d))
) +
  plot_layout(guides = "collect", heights = c(0.08, 1, 1)) +
  plot_annotation(tag_levels = list(c("", "", "a", "c", "b", "d")))

ggsave(file.path(OUT_DIR, "figure_4_ALL.png"), fig4_all,
       width = 14, height = 10, dpi = 150, units = "in")

message("Figure 4 composite saved")