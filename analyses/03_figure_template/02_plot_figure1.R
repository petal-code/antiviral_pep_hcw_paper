# =============================================================================
# 02_plot_figure1.R
# Epi curves of simulated outbreaks
# =============================================================================
source(here::here(
  "analyses",
  "03_figure_template",
  "helper_functions_figure_1to4.R"
))
OUT_DIR <- here("figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

results <- load_results()

# x_max_weeks <- function(sc) if (sc == "WestAfrica") 60 else 80
x_max_weeks <- function(sc) SCENARIO_X_MAX_DAYS[sc] / 7

# Build panels ----
# Weekly infections in entire population -- without OBV (v1 top panels)
ts_infections_allpop <- build_weekly_ts(
  results,
  metric = "infections",
  bin_width = 7,
  efficacy_name = "baseline"
) %>%
  mutate(week = week / 7) # days -> weeks

make_infection_bar <- function(sc, show_errorbars = TRUE) {
  x_max <- x_max_weeks(sc)
  df    <- filter(ts_infections_allpop, scenario == sc, week <= x_max)

  p <- ggplot(df, aes(x = week, y = q50)) +
    geom_col(fill = "grey50", width = 0.8, alpha = 0.5) +
    scale_x_continuous(breaks = seq(0, x_max, 5)) +
    labs(x = "Weeks since outbreak start", y = "Incident infections (all)") +
    theme_fig()
  if (show_errorbars) p <- p + geom_errorbar(aes(ymin = q25, ymax = q75), width = 0.25, linewidth = 0.5)
  p
}

# Weekly deaths in entire population -- without OBV (v2 top panels)
ts_deaths_allpop <- build_weekly_ts(
  results,
  metric = "deaths",
  bin_width = 7,
  efficacy_name = "baseline"
) %>%
  mutate(week = week / 7) # days -> weeks

make_death_bar <- function(sc, show_errorbars = TRUE) {
  x_max <- x_max_weeks(sc)
  df    <- filter(ts_deaths_allpop, scenario == sc, week <= x_max)

  p <- ggplot(df, aes(x = week, y = q50)) +
    geom_col(fill = "grey50", width = 0.8, alpha = 0.5) +
    scale_x_continuous(breaks = seq(0, x_max, 5)) +
    labs(x = "Weeks since outbreak start", y = "Incident deaths (all)") +
    theme_fig()
  if (show_errorbars) p <- p + geom_errorbar(aes(ymin = q25, ymax = q75), width = 0.25, linewidth = 0.5)
  p
}

# Weekly HCW deaths -- without OBV (v3 top panels)
ts_hcw_deaths_base <- build_weekly_ts(
  results,
  metric = "hcw_deaths_incidence",
  bin_width = 7,
  efficacy_name = "baseline"
) %>%
  mutate(week = week / 7, arm = "baseline") # days -> weeks

make_hcw_death_bar_baseline <- function(sc, show_errorbars = TRUE) {
  x_max <- x_max_weeks(sc)
  df    <- filter(ts_hcw_deaths_base, scenario == sc, week <= x_max)

  p <- ggplot(df, aes(x = week, y = q50)) +
    geom_col(fill = "grey50", width = 0.8, alpha = 0.5) +
    scale_x_continuous(limits = c(0, x_max), breaks = seq(0, x_max, 5)) +
    labs(x = "Weeks since outbreak start", y = "Incident HCW deaths") +
    theme_fig()
  if (show_errorbars) p <- p + geom_errorbar(aes(ymin = q25, ymax = q75), width = 0.25, linewidth = 0.5)
  p
}

# Weekly HCW deaths -- with and without OBV (v4 top panels)
ts_hcw_deaths_obv80 <- build_weekly_ts(
  results,
  metric = "hcw_deaths_incidence",
  bin_width = 7,
  efficacy_name = "obv_80",
  coverage_name = "full"
) %>%
  mutate(week = week / 7, arm = "obv_80") # days -> weeks

ts_hcw_inc <- bind_rows(ts_hcw_deaths_base, ts_hcw_deaths_obv80) %>%
  mutate(arm = factor(arm, levels = c("baseline", "obv_80")))

make_hcw_death_bar <- function(sc, show_errorbars = TRUE) {
  arms       <- c("baseline", "obv_80")
  sc_color   <- unname(SCENARIO_COLORS[sc])
  bar_colors <- setNames(c("grey50", sc_color), arms)
  bar_labels <- c(baseline = "Without OBV", obv_80 = "With OBV (80% efficacy, 100% coverage)")
  x_max      <- x_max_weeks(sc)
  df         <- filter(ts_hcw_inc, scenario == sc, week <= x_max)

  p <- ggplot(df, aes(x = week, y = q50, fill = arm)) +
    geom_col(width = 0.5, position = position_dodge(width = 0.5), alpha = 0.5) +
    scale_fill_manual(values = bar_colors, labels = bar_labels, name = NULL) +
    scale_x_continuous(limits = c(0, x_max), breaks = seq(0, x_max, 5)) +
    labs(x = "Weeks since outbreak start", y = "Incident HCW deaths") +
    theme_fig()
  if (show_errorbars) p <- p + geom_errorbar(
    aes(ymin = q25, ymax = q75),
    width = 0.25, linewidth = 0.25,
    position = position_dodge(width = 0.5)
  )
  p
}

# Cumulative HCW deaths -- with and without OBV (bottom panels, all version)
ts_baseline <- build_weekly_ts(
  results, metric = "hcw_deaths", bin_width = 7, efficacy_name = "baseline"
)
ts_obv80 <- build_weekly_ts(
  results, metric = "hcw_deaths", bin_width = 7, efficacy_name = "obv_80", coverage_name = "full"
)

ts_hcw_df <- bind_rows(
  mutate(ts_baseline, arm = "baseline"),
  mutate(ts_obv80,    arm = "obv_80")
) %>%
  mutate(week = week / 7) # days -> weeks

make_ts <- function(sc) {
  arms         <- c("baseline", "obv_80")
  sc_color     <- unname(SCENARIO_COLORS[sc])
  ts_colors    <- setNames(c("grey50", sc_color), arms)
  ts_linetypes <- c(baseline = "solid", obv_80 = "dashed")
  ts_labels    <- c(baseline = "Without OBV", obv_80 = "With OBV (80% efficacy, 100% coverage)")
  x_max        <- x_max_weeks(sc)
  df           <- filter(ts_hcw_df, scenario == sc, week <= x_max) %>%
    mutate(arm = factor(arm, levels = arms))

  ggplot(df, aes(x = week, color = arm, fill = arm)) +
    geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.25, color = NA) +
    geom_line(aes(y = q50, linetype = arm), linewidth = 1) +
    scale_color_manual(values = ts_colors, labels = ts_labels, name = NULL) +
    scale_fill_manual(values = ts_colors, labels = ts_labels, name = NULL) +
    scale_linetype_manual(values = ts_linetypes, labels = ts_labels, name = NULL) +
    scale_x_continuous(limits = c(0, x_max), breaks = seq(0, x_max, 5)) +
    labs(x = "Weeks since outbreak start", y = "Cumulative HCW deaths") +
    theme_fig() +
    theme(legend.key.width = unit(1, "cm"))
}

# Combine panels ----
# Column headers
make_header <- function(label) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = label, fontface = "bold", size = 4.5) +
    theme_void()
}

# Version 1: weekly deaths + cumulative HCW deaths --- likely will use this variant as the main fig
fig1_v2 <- ((make_header("West Africa archetype") | make_header("DRC archetype")) /
              ((make_death_bar("WestAfrica") | make_death_bar("DRC")) + plot_layout(axis_titles = "collect")) /
              ((make_ts("WestAfrica") | make_ts("DRC")) + plot_layout(axis_titles = "collect"))) +
  plot_layout(heights = c(0.2, 1, 2)) +
  plot_annotation(tag_levels = list(c("", "", "a ", "b ", "c ", "d ")))

ggsave(
  file.path(OUT_DIR, "figure_1_all-deaths-baseline-only.png"),
  fig1_v2, width = 10, height = 6.5, dpi = 400, units = "in"
)

# Version 2: weekly infections + cumulative HCW deaths
fig1_v1 <- ((make_header("West Africa archetype") | make_header("DRC archetype")) /
              ((make_infection_bar("WestAfrica") | make_infection_bar("DRC")) + plot_layout(axis_titles = "collect")) /
              ((make_ts("WestAfrica") | make_ts("DRC")) + plot_layout(axis_titles = "collect"))) +
  plot_layout(heights = c(0.2, 1, 2)) +
  plot_annotation(tag_levels = list(c("", "", "a ", "b ", "c ", "d ")))

ggsave(
  file.path(OUT_DIR, "figure_1_all-infections-baseline-only.png"),
  fig1_v1, width = 10, height = 6.5, dpi = 400, units = "in"
)

# Version 3: weekly HCW deaths (baseline only) + cumulative HCW deaths
fig1_v3 <- (
  (make_header("West Africa archetype") | make_header("DRC archetype")) /
    ((make_hcw_death_bar_baseline("WestAfrica") | make_hcw_death_bar_baseline("DRC")) + plot_layout(axis_titles = "collect")) /
    ((make_ts("WestAfrica") | make_ts("DRC"))) + plot_layout(axis_titles = "collect")) +
  plot_layout(heights = c(0.2, 1, 2)) +
  plot_annotation(tag_levels = list(c("", "", "a ", "b ", "c ", "d ")))

ggsave(
  file.path(OUT_DIR, "figure_1_HCW-deaths-baseline-only.png"),
  fig1_v3, width = 10, height = 6.5, dpi = 400, units = "in"
)

# Version 4: weekly HCW deaths (with/without OBV) + cumulative HCW deaths
fig1_v4 <- ((make_header("West Africa archetype") | make_header("DRC archetype")) /
              ((make_hcw_death_bar("WestAfrica") | make_hcw_death_bar("DRC")) + plot_layout(axis_titles = "collect")) /
              ((make_ts("WestAfrica") | make_ts("DRC")) + plot_layout(axis_titles = "collect"))) +
  plot_layout(heights = c(0.2, 1, 2)) +
  plot_annotation(tag_levels = list(c("", "", "a ", "c ", "b ", "d ")))

ggsave(
  file.path(OUT_DIR, "figure_1_HCW-deaths-baseline-obv.png"),
  fig1_v4, width = 10, height = 6.5, dpi = 400, units = "in"
)

# No-errorbar variants ----
fig1_v2_ne <- ((make_header("West Africa archetype") | make_header("DRC archetype")) /
                 ((make_death_bar("WestAfrica", FALSE) | make_death_bar("DRC", FALSE)) + plot_layout(axis_titles = "collect")) /
                 ((make_ts("WestAfrica") | make_ts("DRC")) + plot_layout(axis_titles = "collect"))) +
  plot_layout(heights = c(0.2, 1, 2)) +
  plot_annotation(tag_levels = list(c("", "", "a ", "b ", "c ", "d ")))

ggsave(
  file.path(OUT_DIR, "figure_1_all-deaths-baseline-only_no-errorbars.png"),
  fig1_v2_ne, width = 10, height = 6.5, dpi = 400, units = "in"
)

fig1_v1_ne <- ((make_header("West Africa archetype") | make_header("DRC archetype")) /
                 ((make_infection_bar("WestAfrica", FALSE) | make_infection_bar("DRC", FALSE)) + plot_layout(axis_titles = "collect")) /
                 ((make_ts("WestAfrica") | make_ts("DRC")) + plot_layout(axis_titles = "collect"))) +
  plot_layout(heights = c(0.2, 1, 2)) +
  plot_annotation(tag_levels = list(c("", "", "a ", "b ", "c ", "d ")))

ggsave(
  file.path(OUT_DIR, "figure_1_all-infections-baseline-only_no-errorbars.png"),
  fig1_v1_ne, width = 10, height = 6.5, dpi = 400, units = "in"
)

fig1_v3_ne <- (
  (make_header("West Africa archetype") | make_header("DRC archetype")) /
    ((make_hcw_death_bar_baseline("WestAfrica", FALSE) | make_hcw_death_bar_baseline("DRC", FALSE)) + plot_layout(axis_titles = "collect")) /
    ((make_ts("WestAfrica") | make_ts("DRC"))) + plot_layout(axis_titles = "collect")) +
  plot_layout(heights = c(0.2, 1, 2)) +
  plot_annotation(tag_levels = list(c("", "", "a ", "b ", "c ", "d ")))

ggsave(
  file.path(OUT_DIR, "figure_1_HCW-deaths-baseline-only_no-errorbars.png"),
  fig1_v3_ne, width = 10, height = 6.5, dpi = 400, units = "in"
)

fig1_v4_ne <- ((make_header("West Africa archetype") | make_header("DRC archetype")) /
                 ((make_hcw_death_bar("WestAfrica", FALSE) | make_hcw_death_bar("DRC", FALSE)) + plot_layout(axis_titles = "collect")) /
                 ((make_ts("WestAfrica") | make_ts("DRC")) + plot_layout(axis_titles = "collect"))) +
  plot_layout(heights = c(0.2, 1, 2)) +
  plot_annotation(tag_levels = list(c("", "", "a ", "c ", "b ", "d ")))

ggsave(
  file.path(OUT_DIR, "figure_1_HCW-deaths-baseline-obv_no-errorbars.png"),
  fig1_v4_ne, width = 10, height = 6.5, dpi = 400, units = "in"
)

message("Figure 1 variants saved")
