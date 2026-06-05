# =============================================================================
# 02_plot_figure1.R
# =============================================================================
library(patchwork)
source(here::here("analyses", "03_figure_template", "helper_functions_figure_1to4.R"))
OUT_DIR <- here("figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

results <- load_results()

# =============================================================================
# Panel A/C: Weekly infection incidence -- bar (median) + error bar (95% CI)
# =============================================================================
ts_infections <- build_weekly_ts(results, metric = "infections",
                                 bin_width = 7,
                                 efficacy_name = "baseline") %>%
  mutate(week = week / 7)

make_infection_bar <- function(sc) {
  df    <- filter(ts_infections, scenario == sc)
  color <- unname(SCENARIO_COLORS[sc])
  
  ggplot(df, aes(x = week, y = q50)) +
    geom_col(fill = color, alpha = 0.70, width = 0.8) +
    geom_errorbar(aes(ymin = q025, ymax = q975),
                  width = 0.3, linewidth = 0.5, color = "grey30") +
    scale_x_continuous(breaks = seq(0, 35, by = 5), limits = c(0, 35),
                       expand = expansion(add = c(0.5, 0.5))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(x = "Weeks since outbreak start",
         y = "Weekly new infections") +
    theme_fig()
}

# =============================================================================
# Panel B/D: Cumulative HCW deaths -- baseline vs OBV 80% full coverage
# =============================================================================
ts_baseline <- build_weekly_ts(results, metric = "hcw_deaths",
                               bin_width = 7,
                               efficacy_name = "baseline")
ts_obv80    <- build_weekly_ts(results, metric = "hcw_deaths",
                               bin_width = 7,
                               efficacy_name = "obv_80",
                               coverage_name = "full")

ts_hcw_df <- bind_rows(
  mutate(ts_baseline, arm = "baseline"),
  mutate(ts_obv80,    arm = "obv_80")
) %>%
  mutate(week = week / 7)

make_ts <- function(sc) {
  arms         <- c("baseline", "obv_80")
  sc_color     <- unname(SCENARIO_COLORS[sc])
  ts_colors    <- setNames(c("grey50", sc_color), arms)
  ts_linetypes <- c(baseline = "solid", obv_80 = "dashed")
  ts_labels    <- c(baseline = "Without OBV", obv_80 = "With OBV (80% efficacy)")
  
  df <- filter(ts_hcw_df, scenario == sc) %>%
    mutate(arm = factor(arm, levels = arms))
  
  ggplot(df, aes(x = week, color = arm, fill = arm)) +
    geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.12, color = NA) +
    geom_ribbon(aes(ymin = q25,  ymax = q75),  alpha = 0.25, color = NA) +
    geom_line(aes(y = q50, linetype = arm), linewidth = 1.0) +
    scale_color_manual(values = ts_colors, labels = ts_labels, name = NULL) +
    scale_fill_manual( values = ts_colors, labels = ts_labels, name = NULL) +
    scale_linetype_manual(values = ts_linetypes, labels = ts_labels, name = NULL) +
    scale_x_continuous(breaks = seq(0, 35, by = 5), limits = c(0, 35),
                       expand = expansion(add = c(0.5, 0.5))) +
    labs(x = "Weeks since outbreak start",
         y = "Cumulative HCW deaths") +
    theme_fig() +
    theme(legend.key.width = unit(2, "cm"))
}

# =============================================================================
# Save individual panels
# =============================================================================
fig1a <- make_infection_bar("WestAfrica")
fig1b <- make_ts("WestAfrica")
fig1c <- make_infection_bar("DRC")
fig1d <- make_ts("DRC")

ggsave(file.path(OUT_DIR, "figure_1_a.png"), fig1a,
       width = 7, height = 5, dpi = 150)
ggsave(file.path(OUT_DIR, "figure_1_b.png"), fig1b,
       width = 7, height = 5, dpi = 150)
ggsave(file.path(OUT_DIR, "figure_1_c.png"), fig1c,
       width = 7, height = 5, dpi = 150)
ggsave(file.path(OUT_DIR, "figure_1_d.png"), fig1d,
       width = 7, height = 5, dpi = 150)

# =============================================================================
# Combined figure layouts
# =============================================================================
make_header <- function(label) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = label,
             fontface = "bold", size = 5) +
    theme_void()
}

# Side-by-side layout
fig1_all <- (
  (make_header("West Africa") | make_header("DRC")) /
    (fig1a | fig1c) /
    (fig1b | fig1d)
) +
  plot_layout(heights = c(0.12, 1, 2)) +
  plot_annotation(tag_levels = list(c("", "", "a ", "c ", "b ", "d ")))

ggsave(file.path(OUT_DIR, "figure_1_ALL.png"), fig1_all,
       width = 11, height = 6.5, dpi = 150, units = "in")

# Stacked layout
fig1_all_v2 <- (
  make_header("West Africa") / fig1a / fig1b /
    make_header("DRC")         / fig1c / fig1d
) +
  plot_layout(heights = c(0.12, 1, 2, 0.12, 1, 2)) +
  plot_annotation(tag_levels = list(c("", "a ", "b ", "", "c ", "d ")))

ggsave(file.path(OUT_DIR, "figure_1_ALL_v2.png"), fig1_all_v2,
       width = 6.5, height = 11, dpi = 150, units = "in")

message("Figure 1 saved")