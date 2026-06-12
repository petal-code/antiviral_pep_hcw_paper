# =============================================================================
# 02_plot_figure3.R
# Coverage scenario comparison
# Reads pre-computed CSV from output_figgen/figure_3_run_summary.csv
# and output_figgen/figure_3_averted_weekly_ts.csv
# Run 02_extract_figure3.R first.
# =============================================================================
source(here::here("analyses", "03_figure_template", "helper_functions_figure_1to4.R"))
OUT_DIR <- here("figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FIG3_EFFICACY_LEVELS <- c("obv_50", "obv_60", "obv_70", "obv_80", "obv_90")
FIG3_EFFICACY_LABELS <- c("50%", "60%", "70%", "80%", "90%")

# =============================================================================
# Helper functions
# =============================================================================
make_col_header <- function(label) {
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

x_max_weeks <- function(sc) SCENARIO_X_MAX_DAYS[sc] / 7

# =============================================================================
# Load data and build particle-level data frame
# =============================================================================
run_df <- read.csv(here("output_figgen", "figure_3_run_summary.csv"),
                   stringsAsFactors = FALSE)

pdf <- make_particle_df(run_df) %>%
  filter(arm != "baseline") %>%
  mutate(
    coverage_name  = sub("__.*", "", arm),
    eff_name       = sub(".*__", "", arm),
    arm_label      = factor(FIG3_EFFICACY_LABELS[match(eff_name, FIG3_EFFICACY_LEVELS)],
                            levels = FIG3_EFFICACY_LABELS),
    coverage_label = factor(COVERAGE_LABELS[match(coverage_name, COVERAGE_LEVELS)],
                            levels = COVERAGE_LABELS),
    scenario_label = factor(SCENARIO_LABELS[scenario], levels = SCENARIO_LABELS)
  )

save_figure_data(pdf, "figure_3_particle_df.csv")

# =============================================================================
# Coverage curve panels
# =============================================================================
make_coverage_plot <- function(cs) {
  spec   <- COVERAGE_SPECS[[cs]]
  x_max  <- max(SCENARIO_X_MAX_DAYS)
  t_days <- seq(0, x_max, by = 1)
  df     <- data.frame(week = t_days / 7,
                       coverage = coverage_at_time(t_days, spec) * 100)
  ggplot(df, aes(x = week, y = coverage)) +
    geom_line(color = COVERAGE_COLORS[cs], linewidth = 1.2) +
    scale_y_continuous(limits = c(0, 100), labels = function(x) paste0(x, "%")) +
    scale_x_continuous(breaks = seq(0, x_max / 7, by = 13)) +
    labs(x = "Weeks since outbreak start", y = "Antiviral coverage") +
    theme_fig()
}

# =============================================================================
# Box plot panels (deaths averted / days lost averted)
# =============================================================================
make_box_plot <- function(cs, metric, y_label) {
  df <- pdf %>%
    filter(coverage_name == cs, !is.na(.data[[metric]])) %>%
    mutate(fill_group = paste(as.character(scenario_label),
                              as.character(arm_label), sep = "."))
  
  build_sc_fill <- function(sc_key) {
    sc_lbl    <- SCENARIO_LABELS[sc_key]
    base_col  <- SCENARIO_COLORS[sc_key]
    light_col <- if (sc_key == "WestAfrica") "#fdd8a0" else "#b2e4d8"
    dark_col  <- rgb(t(col2rgb(base_col) * 0.7), maxColorValue = 255)
    cols      <- colorRampPalette(c(light_col, dark_col))(length(FIG3_EFFICACY_LABELS))
    setNames(cols, paste(sc_lbl, FIG3_EFFICACY_LABELS, sep = "."))
  }
  
  fill_vals     <- c(build_sc_fill("WestAfrica"), build_sc_fill("DRC"))
  legend_breaks <- c(paste(SCENARIO_LABELS["WestAfrica"], "80%", sep = "."),
                     paste(SCENARIO_LABELS["DRC"],        "80%", sep = "."))
  
  ggplot(df, aes(x = arm_label, y = .data[[metric]], fill = fill_group)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, color = "black", linewidth = 0.4,
                 position = position_dodge(0.75)) +
    scale_fill_manual(values = fill_vals, breaks = rev(legend_breaks),
                      labels = c("DRC archetype", "West Africa archetype"), name = NULL) +
    scale_y_continuous(limits = c(0, 100), labels = function(x) paste0(x, "%")) +
    labs(x = "Antiviral efficacy", y = y_label) +
    theme_fig()
}

# =============================================================================
# Build and save: deaths averted / days lost averted figures
# =============================================================================
h1 <- make_col_header("Constant, Full Coverage")
h2 <- make_col_header("Ramp Up to High Coverage")
h3 <- make_col_header("Ramp Up to Medium Coverage")

p_a <- make_coverage_plot(COVERAGE_LEVELS[1])
p_b <- make_coverage_plot(COVERAGE_LEVELS[2])
p_c <- make_coverage_plot(COVERAGE_LEVELS[3])

p_d <- make_box_plot(COVERAGE_LEVELS[1], "pct_hcw_deaths_averted", "HCW deaths averted")
p_e <- make_box_plot(COVERAGE_LEVELS[2], "pct_hcw_deaths_averted", "HCW deaths averted")
p_f <- make_box_plot(COVERAGE_LEVELS[3], "pct_hcw_deaths_averted", "HCW deaths averted")
p_g <- make_box_plot(COVERAGE_LEVELS[1], "pct_days_lost_averted",  "HCW days lost averted")
p_h <- make_box_plot(COVERAGE_LEVELS[2], "pct_days_lost_averted",  "HCW days lost averted")
p_i <- make_box_plot(COVERAGE_LEVELS[3], "pct_days_lost_averted",  "HCW days lost averted")

figure_3_deaths <- (
  (h1 | h2 | h3) /
    ((p_a | p_b | p_c) + plot_layout(axis_titles = "collect")) /
    ((p_d | p_e | p_f) + plot_layout(axis_titles = "collect"))
) +
  plot_layout(guides = "collect", heights = c(0.2, 1, 3)) +
  plot_annotation(tag_levels = list(c("", "", "", "a ", "b ", "c ", "d ", "e ", "f "))) &
  theme(legend.position = "bottom")

save_fig("figure_3_deaths-averted", figure_3_deaths, 10, 6.5)

figure_3_days_lost <- (
  (h1 | h2 | h3) /
    ((p_a | p_b | p_c) + plot_layout(axis_titles = "collect")) /
    ((p_g | p_h | p_i) + plot_layout(axis_titles = "collect"))
) +
  plot_layout(guides = "collect", axes = "collect", heights = c(0.2, 1, 3)) +
  plot_annotation(tag_levels = list(c("", "", "", "a ", "b ", "c ", "d ", "e ", "f "))) &
  theme(legend.position = "bottom")

save_fig("figure_3_days-averted", figure_3_days_lost, 10, 6.5)

message("Figure 3 variants saved")

# =============================================================================
# Cumulative averted burden panels (2 scenarios x 3 coverage levels),
# with efficacy 50-90% as colored lines within each panel
# =============================================================================
averted_ts <- read.csv(here("output_figgen", "figure_3_averted_weekly_ts.csv"),
                       stringsAsFactors = FALSE)

EFF_COLORS <- setNames(
  colorRampPalette(c("#cfe8f7", "#08306b"))(length(FIG3_EFFICACY_LABELS)),
  FIG3_EFFICACY_LABELS
)

make_cumulative_panel <- function(metric_name, y_label, sc, cov) {
  x_max <- x_max_weeks(sc)
  
  df <- averted_ts %>%
    filter(metric == metric_name, scenario == sc, coverage_name == cov)
  
  cum_df <- summarise_cumulative_ts(df, extra_group_cols = c("coverage_name", "eff_name")) %>%
    filter(week <= x_max) %>%
    mutate(eff_label = factor(FIG3_EFFICACY_LABELS[match(eff_name, FIG3_EFFICACY_LEVELS)],
                              levels = FIG3_EFFICACY_LABELS))
  
  ggplot(cum_df, aes(x = week, color = eff_label, fill = eff_label)) +
    geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.12, color = NA) +
    geom_line(aes(y = q50), linewidth = 1) +
    scale_color_manual(values = EFF_COLORS, name = "Antiviral efficacy") +
    scale_fill_manual(values = EFF_COLORS, name = "Antiviral efficacy") +
    scale_x_continuous(limits = c(0, x_max), breaks = seq(0, x_max, 13)) +
    labs(x = "Weeks since outbreak start", y = y_label) +
    theme_fig()
}

build_cumulative_figure <- function(metric_name, y_label, filename_base) {
  p_wa_full      <- make_cumulative_panel(metric_name, y_label, "WestAfrica", "full")
  p_wa_ramp_high <- make_cumulative_panel(metric_name, y_label, "WestAfrica", "ramp_high")
  p_wa_ramp_low  <- make_cumulative_panel(metric_name, y_label, "WestAfrica", "ramp_low")
  p_drc_full      <- make_cumulative_panel(metric_name, y_label, "DRC", "full")
  p_drc_ramp_high <- make_cumulative_panel(metric_name, y_label, "DRC", "ramp_high")
  p_drc_ramp_low  <- make_cumulative_panel(metric_name, y_label, "DRC", "ramp_low")
  
  fig <- (
    (h1 | h2 | h3) /
      ((p_wa_full | p_wa_ramp_high | p_wa_ramp_low) + plot_layout(axis_titles = "collect")) /
      ((p_drc_full | p_drc_ramp_high | p_drc_ramp_low) + plot_layout(axis_titles = "collect"))
  ) +
    plot_layout(guides = "collect", heights = c(0.2, 1, 1)) +
    plot_annotation(tag_levels = list(c("", "", "", "a ", "b ", "c ", "d ", "e ", "f "))) &
    theme(legend.position = "bottom")
  
  save_fig(filename_base, fig, 10, 6.5)
}

build_cumulative_figure("averted_hcw_deaths_incidence",
                        "Cumulative HCW deaths averted",
                        "figure_3_cumulative-hcw-deaths-averted")

build_cumulative_figure("averted_hcw_days_lost_incidence",
                        "Cumulative HCW days lost averted",
                        "figure_3_cumulative-hcw-days-lost-averted")


# =============================================================================
# Weekly HCW deaths incidence overlay (80% efficacy):
# baseline (no OBV) + full / ramp_high / ramp_low, per scenario
# =============================================================================
weekly_80 <- read.csv(here("output_figgen", "figure_3_weekly_hcw_deaths_80.csv"),
                      stringsAsFactors = FALSE)

# Baseline is the same regardless of which arm it came from (OBV doesn't
# affect transmission); take it from the "full" arm only, label as "baseline".
weekly_80_clean <- weekly_80 %>%
  filter(
    (arm == "obv") |
      (arm == "baseline" & coverage_name == "full")
  ) %>%
  mutate(
    line_group = ifelse(arm == "baseline", "baseline", coverage_name)
  )

LINE_LEVELS <- c("baseline", "full", "ramp_high", "ramp_low")
LINE_LABELS <- c("No antiviral", "Full (100%)", "Ramp high (0%->80%)", "Ramp low (0%->50%)")
LINE_COLORS <- c(baseline = "grey40", full = "#1a9641",
                 ramp_high = "#fdae61", ramp_low = "#d7191c")

summarise_weekly <- function(df, sc) {
  x_max <- x_max_weeks(sc)
  df %>%
    filter(scenario == sc, week <= x_max) %>%
    group_by(line_group, week) %>%
    summarise(q50 = median(value), .groups = "drop") %>%
    mutate(line_group = factor(line_group, levels = LINE_LEVELS))
}

make_weekly_deaths_panel <- function(sc) {
  x_max <- x_max_weeks(sc)
  df <- summarise_weekly(weekly_80_clean, sc)
  ggplot(df, aes(x = week, y = q50, color = line_group)) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = LINE_COLORS,
                       breaks = LINE_LEVELS, labels = LINE_LABELS,
                       name = NULL) +
    scale_x_continuous(limits = c(0, x_max), breaks = seq(0, x_max, 13)) +
    labs(x = "Weeks since outbreak start", y = "Incident HCW deaths") +
    theme_fig()
}

p_d <- make_weekly_deaths_panel("WestAfrica")
p_e <- make_weekly_deaths_panel("DRC")

fig3_weekly_deaths <- (p_d | p_e) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = list(c("d ", "e "))) &
  theme(legend.position = "bottom")

save_fig("figure_3_weekly-hcw-deaths-80pct", fig3_weekly_deaths, 10, 4)

message("Figure 3 weekly HCW deaths overlay (80% efficacy) saved")

message("Figure 3 cumulative averted-burden panels saved")