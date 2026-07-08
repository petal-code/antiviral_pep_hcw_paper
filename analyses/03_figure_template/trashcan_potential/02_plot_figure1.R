# =============================================================================
# 02_plot_figure1.R
# Epi curves of simulated outbreaks
# Reads pre-computed CSV from output_figgen/figure_1_weekly_ts.csv
# Run 01_extract_figure1.R first.
#
# NOTE: Only "figure_1_all-deaths-baseline-only" is generated. All other
# figure variants (infections, HCW-deaths-baseline-only, HCW-deaths-
# baseline-obv, no-errorbar loop) have been removed.
# =============================================================================
source(here::here("analyses", "03_figure_template", "helper_functions_figure_1to4.R"))
OUT_DIR <- here("figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Load pre-computed quantile data ----
ts_quantiles <- read.csv(here("output_figgen", "figure_1_weekly_ts.csv"),
                         stringsAsFactors = FALSE)

x_max_weeks <- function(sc) SCENARIO_X_MAX_DAYS[sc] / 7

get_ts <- function(sc, metric_name, arm_name) {
  x_max <- x_max_weeks(sc)
  ts_quantiles %>%
    filter(scenario == sc, metric == metric_name, arm == arm_name, week <= x_max)
}

# Panel functions ----
make_death_bar <- function(sc) {
  x_max <- x_max_weeks(sc)
  df    <- get_ts(sc, "deaths", "baseline")
  ggplot(df, aes(x = week, y = q50)) +
    geom_line(color = "grey50", linewidth = 1) +
    geom_ribbon(aes(ymin = q10, ymax = q90), fill = "grey50", alpha = 0.25, color = NA) +
    geom_line(aes(y = mean), color = "firebrick", linewidth = 1, linetype = "dashed") +
    scale_x_continuous(limits = c(0, x_max), breaks = seq(0, x_max, 5)) +
    labs(x = "Weeks since outbreak start", y = "Incident deaths (all)") +
    theme_fig()
}

make_ts <- function(sc) {
  arms         <- c("baseline", "obv")
  sc_color     <- unname(SCENARIO_COLORS[sc])
  ts_colors    <- setNames(c("grey50", sc_color), arms)
  ts_linetypes <- c(baseline = "solid", obv = "dashed")
  ts_labels    <- c(baseline = "Without antiviral", obv = "With antiviral (80% efficacy, 80% coverage)")
  x_max        <- x_max_weeks(sc)
  df <- bind_rows(
    get_ts(sc, "hcw_deaths", "baseline"),
    get_ts(sc, "hcw_deaths", "obv")
  ) %>% mutate(arm = factor(arm, levels = arms))
  ggplot(df, aes(x = week, color = arm, fill = arm)) +
    geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.25, color = NA) +
    geom_line(aes(y = q50, linetype = arm), linewidth = 1) +
    # geom_line(aes(y = q75, linetype = arm), linewidth = 1) +
    scale_color_manual(values = ts_colors, labels = ts_labels, name = NULL) +
    scale_fill_manual(values = ts_colors, labels = ts_labels, name = NULL) +
    scale_linetype_manual(values = ts_linetypes, labels = ts_labels, name = NULL) +
    scale_x_continuous(limits = c(0, x_max), breaks = seq(0, x_max, 5)) +
    labs(x = "Weeks since outbreak start", y = "Cumulative HCW deaths") +
    theme_fig() +
    theme(legend.key.width = unit(1, "cm"))
}

make_header <- function(label) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = label, fontface = "bold", size = 4.5) +
    theme_void()
}

# Save both PNG and PDF
save_fig <- function(filename_base, plot, width, height) {
  ggsave(file.path(OUT_DIR, paste0(filename_base, ".tiff")),
         plot, width = width, height = height, dpi = 400, units = "in")
  ggsave(file.path(OUT_DIR, paste0(filename_base, ".png")),
         plot, width = width, height = height, dpi = 400, units = "in")
  ggsave(file.path(OUT_DIR, paste0(filename_base, ".pdf")),
         plot, width = width, height = height, units = "in")
}

# Save figure ----
fig1_v2 <- ((make_header("West Africa archetype") | make_header("DRC archetype")) /
              ((make_death_bar("WestAfrica") | make_death_bar("DRC")) + plot_layout(axis_titles = "collect")) /
              ((make_ts("WestAfrica") | make_ts("DRC")) + plot_layout(axis_titles = "collect"))) +
  plot_layout(heights = c(0.2, 1, 2)) +
  plot_annotation(tag_levels = list(c("", "", "a ", "b ", "c ", "d ")))
save_fig("figure_1_all-deaths-baseline-only", fig1_v2, 10, 6.5)

message("Figure 1 (all-deaths-baseline-only) saved")

############### aggregating number for the paper
ts <- read.csv(here("output_figgen", "figure_1_weekly_ts.csv"))
ts %>%
  filter(metric == "hcw_deaths") %>%
  group_by(scenario, arm) %>%
  filter(week == max(week)) %>%
  select(scenario, arm, q50, q25, q75) %>%
  pivot_wider(names_from = arm, values_from = c(q50, q25, q75)) %>%
  mutate(
    averted_median = q50_baseline - q50_obv,
    pct_median     = 100 * averted_median / q50_baseline
  ) %>%
  select(scenario,
         baseline_median = q50_baseline, baseline_lo = q025_baseline, baseline_hi = q975_baseline,
         obv_median = q50_obv, obv_lo = q025_obv, obv_hi = q975_obv,
         averted_median, pct_median) %>%
  mutate(across(where(is.numeric), ~round(., 1))) %>%
  as.data.frame() %>%
  print()

particle_df <- read.csv(here("output_figgen", "figure_1_particle_cum_hcw_deaths.csv"))

particle_df %>%
  pivot_wider(names_from = arm, values_from = cum_hcw_deaths) %>%
  mutate(
    averted = baseline - obv,
    pct_averted = 100 * averted / baseline
  ) %>%
  group_by(scenario) %>%
  summarise(
    averted_median = median(averted),
    averted_lo     = quantile(averted, 0.025),
    averted_hi     = quantile(averted, 0.975),
    pct_median     = median(pct_averted),
    pct_lo         = quantile(pct_averted, 0.025),
    pct_hi         = quantile(pct_averted, 0.975)
  ) %>%
  mutate(across(where(is.numeric), ~round(., 1))) %>%
  as.data.frame() %>%
  print()