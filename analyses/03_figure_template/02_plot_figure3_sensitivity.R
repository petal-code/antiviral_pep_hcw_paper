# =============================================================================
# 02_plot_figure3_sensitivity.R
# Visualise Figure 3 sensitivity analysis results.
# Reads: output_figgen/figure_3_sens_run_summary.csv
# =============================================================================
source(here::here("analyses", "03_figure_template", "helper_functions_figure_1to4.R"))
OUT_DIR <- here("figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

sens_raw <- read.csv(here("output_figgen", "figure_3_sens_run_summary.csv"),
                     stringsAsFactors = FALSE)
save_fig <- function(filename_base, plot, width, height) {
  ggsave(file.path(OUT_DIR, paste0(filename_base, ".png")),
         plot, width = width, height = height, dpi = 400, units = "in")
  ggsave(file.path(OUT_DIR, paste0(filename_base, ".pdf")),
         plot, width = width, height = height, units = "in")
}
# =============================================================================
# Particle-level summaries
# =============================================================================
particle_df <- sens_raw %>%
  group_by(scenario, particle_id, arm, onset_day, max_cov, efficacy) %>%
  summarise(
    prevented_hcw      = sum(prevented_hcw,      na.rm = TRUE),
    counterfactual_hcw = sum(counterfactual_hcw, na.rm = TRUE),
    hcw_days_lost      = sum(hcw_days_lost,      na.rm = TRUE),
    baseline_days_lost = sum(baseline_days_lost, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    pct_hcw_deaths_averted = ifelse(
      counterfactual_hcw > 0,
      100 * prevented_hcw / counterfactual_hcw,
      NA_real_
    ),
    pct_days_lost_averted = ifelse(
      baseline_days_lost > 0,
      100 * (baseline_days_lost - hcw_days_lost) / baseline_days_lost,
      NA_real_
    )
  ) %>%
  filter(!is.na(pct_hcw_deaths_averted)) %>%
  mutate(
    onset_label    = factor(sprintf("Day %d", onset_day),
                            levels = sprintf("Day %d", seq(0, 100, 20))),
    cov_label      = factor(sprintf("%d%%", round(max_cov * 100)),
                            levels = sprintf("%d%%", seq(20, 100, 20))),
    scenario_label = factor(SCENARIO_LABELS[scenario], levels = SCENARIO_LABELS)
  )

# =============================================================================
# Heatmap: median % HCW deaths averted
# =============================================================================
heatmap_df <- particle_df %>%
  group_by(scenario, onset_day, max_cov, onset_label, cov_label) %>%
  summarise(
    median_pct_deaths = median(pct_hcw_deaths_averted, na.rm = TRUE),
    median_pct_days   = median(pct_days_lost_averted,  na.rm = TRUE),
    .groups = "drop"
  )

make_heatmap <- function(sc, metric = "median_pct_deaths",
                         fill_label = "HCW deaths\naverted") {
  df       <- filter(heatmap_df, scenario == sc)
  sc_color <- SCENARIO_COLORS[[sc]]
  title_map <- c(WestAfrica = "West Africa Archetype", DRC = "DRC Archetype")
  ggplot(df, aes(x = onset_label, y = cov_label, fill = .data[[metric]])) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.1f%%", .data[[metric]])),
              size = 3.5, fontface = "bold", color = "grey20") +
    scale_fill_gradient(low = "white", high = sc_color,
                        limits = c(0, 100),
                        labels = function(x) paste0(x, "%"),
                        name   = fill_label) +
    labs(x     = "Scale-up timing",
         y     = "Maximum coverage",
         title = title_map[sc]) +
    theme_fig() +
    theme(legend.position = "right",
          panel.grid      = element_blank(),
          axis.line       = element_blank(),
          plot.title      = element_text(size = 11, hjust = 0.5))
}

# Deaths averted heatmap
fig_heatmap_deaths <- (make_heatmap("WestAfrica") | make_heatmap("DRC")) +
  plot_annotation(tag_levels = list(c("a ", "b ")))
save_fig("figure_3_sens_heatmap_deaths", fig_heatmap_deaths, 12, 4.5)

# Days lost averted heatmap
fig_heatmap_days <- (
  make_heatmap("WestAfrica", "median_pct_days", "HCW days lost\naverted") |
    make_heatmap("DRC",        "median_pct_days", "HCW days lost\naverted")
) +
  plot_annotation(tag_levels = list(c("a ", "b ")))
save_fig("figure_3_sens_heatmap_days", fig_heatmap_days, 12, 4.5)

message("Figure 3 sensitivity heatmaps saved")

# =============================================================================
# Coverage schematic panel
# =============================================================================
make_coverage_schematic <- function() {
  onset   <- 100
  plateau <- 300
  max_cov <- 0.80
  
  # Manual curve: flat 0 until onset, then smooth S-shape to plateau
  t_flat    <- seq(0, onset, by = 1)
  t_ramp    <- seq(onset, plateau, by = 1)
  t_flat2   <- seq(plateau, 420, by = 1)
  
  cov_flat  <- rep(0, length(t_flat))
  cov_ramp  <- max_cov / (1 + exp(-0.04 * (t_ramp - (onset + plateau) / 2)))
  cov_ramp  <- cov_ramp - min(cov_ramp)  # force start at 0
  cov_ramp  <- cov_ramp / max(cov_ramp) * max_cov  # force end at max_cov
  cov_flat2 <- rep(max_cov, length(t_flat2))
  
  df_s <- data.frame(
    t        = c(t_flat, t_ramp[-1], t_flat2[-1]),
    coverage = c(cov_flat, cov_ramp[-1], cov_flat2[-1])
  )
  
  ggplot(df_s, aes(x = t, y = coverage)) +
    geom_line(color = "#CC3399", linewidth = 1.2) +
    # Scale-up timing arrow: x=0 to onset, sitting on y=0
    annotate("segment",
             x = 0, xend = onset, y = 0.02, yend = 0.02,
             arrow = arrow(ends = "both", length = unit(0.12, "cm"), type = "closed"),
             color = "grey40", linewidth = 0.5) +
    annotate("text", x = onset / 2, y = max_cov * 0.1,
             label = paste0("Scale-up\ntiming"), size = 2.6, color = "grey20", hjust = 0.5) +
    # Max coverage arrow: x=plateau, y=0 to max_cov
    annotate("segment",
             x = plateau, xend = plateau, y = 0, yend = max_cov,
             arrow = arrow(ends = "both", length = unit(0.12, "cm"), type = "closed"),
             color = "grey40", linewidth = 0.5) +
    annotate("text", x = plateau + 12, y = max_cov / 2,
             label = paste0("Maximum\ncoverage"), size = 2.6, color = "grey20",
             hjust = 0, vjust = 0.5) +
    scale_x_continuous(limits = c(0, 420), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1),   expand = c(0, 0.02)) +
    labs(x     = "Days since outbreak start",
         y     = "Antiviral coverage") +
    theme_fig() +
    theme(plot.margin = margin(10, 25, 10, 5),
          axis.text   = element_blank(),
          axis.ticks  = element_blank())
}

# =============================================================================
# Boxplot panels
# =============================================================================
make_boxplot <- function(sc, metric = "pct_hcw_deaths_averted",
                         y_label = "HCW deaths averted") {
  df        <- filter(particle_df, scenario == sc)
  sc_color  <- unname(SCENARIO_COLORS[sc])
  light_col <- if (sc == "WestAfrica") "#fdd8a0" else "#a8ddb5"
  fill_cols <- setNames(
    colorRampPalette(c(light_col, sc_color))(5),
    sprintf("%d%%", seq(20, 100, 20))
  )
  title_map <- c(WestAfrica = "West Africa Archetype", DRC = "DRC Archetype")
  ggplot(df, aes(x = onset_label, y = .data[[metric]], fill = cov_label)) +
    geom_boxplot(outlier.shape = NA, width = 0.7, color = "black",
                 linewidth = 0.3, position = position_dodge(0.8)) +
    scale_fill_manual(values = fill_cols, name = "Max coverage") +
    scale_y_continuous(limits = c(0, 100),
                       labels = function(x) paste0(x, "%")) +
    labs(x     = "Scale-up timing",
         y     = y_label,
         title = title_map[sc]) +
    theme_fig() +
    theme(legend.position = "bottom",
          plot.title      = element_text(size = 11, hjust = 0.5))
}

p_schematic <- make_coverage_schematic()

# Deaths averted boxplot
fig_boxplot_deaths <- (p_schematic | make_boxplot("WestAfrica") | make_boxplot("DRC")) +
  plot_layout(widths = c(1, 2, 2)) +
  plot_annotation(tag_levels = list(c("a ", "b ", "c "))) &
  theme(legend.position = "bottom")

save_fig("figure_3_sens_boxplot_deaths", fig_boxplot_deaths, 14, 4)
message("Figure 3 sensitivity boxplot (deaths averted) saved")

# Days lost averted boxplot
fig_boxplot_days <- (p_schematic |
                       make_boxplot("WestAfrica",
                                    metric  = "pct_days_lost_averted",
                                    y_label = "HCW days lost averted") |
                       make_boxplot("DRC",
                                    metric  = "pct_days_lost_averted",
                                    y_label = "HCW days lost averted")) +
  plot_layout(widths = c(1, 2, 2)) +
  plot_annotation(tag_levels = list(c("a ", "b ", "c "))) &
  theme(legend.position = "bottom")

save_fig("figure_3_sens_boxplot_days", fig_boxplot_days, 14, 4)
message("Figure 3 sensitivity boxplot (days lost averted) saved")

# =============================================================================
# Numbers for SI Section 3.8 sensitivity analysis text
# =============================================================================

# Immediate scale-up (onset=0) vs delayed (onset=100), across max coverage levels
sens_text <- particle_df %>%
  group_by(scenario, onset_day, max_cov) %>%
  summarise(
    lo_deaths  = quantile(pct_hcw_deaths_averted, 0.025, na.rm = TRUE),
    med_deaths = median(pct_hcw_deaths_averted, na.rm = TRUE),
    hi_deaths  = quantile(pct_hcw_deaths_averted, 0.975, na.rm = TRUE),
    lo_days    = quantile(pct_days_lost_averted, 0.025, na.rm = TRUE),
    med_days   = median(pct_days_lost_averted, na.rm = TRUE),
    hi_days    = quantile(pct_days_lost_averted, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~round(.x, 1)))

# Key comparisons:
# 1. Effect of onset timing at max_cov = 1.0 (100%)
sens_text %>% filter(max_cov == 1.0) %>%
  select(scenario, onset_day, med_deaths, lo_deaths, hi_deaths,
         med_days, lo_days, hi_days) %>%
  as.data.frame() %>% print()

# 2. Effect of max coverage at onset = 0
sens_text %>% filter(onset_day == 0) %>%
  select(scenario, max_cov, med_deaths, lo_deaths, hi_deaths,
         med_days, lo_days, hi_days) %>%
  as.data.frame() %>% print()

# 3. Worst case (onset=100, max_cov=0.2) vs best case (onset=0, max_cov=1.0)
sens_text %>%
  filter((onset_day == 100 & max_cov == 0.2) |
           (onset_day == 0   & max_cov == 1.0)) %>%
  select(scenario, onset_day, max_cov, med_deaths, lo_deaths, hi_deaths) %>%
  as.data.frame() %>% print()
# onset timing at max_cov = 1.0
sens_text %>% filter(max_cov == 1.0) %>%
  select(scenario, onset_day, med_deaths, lo_deaths, hi_deaths) %>%
  as.data.frame() %>% print()