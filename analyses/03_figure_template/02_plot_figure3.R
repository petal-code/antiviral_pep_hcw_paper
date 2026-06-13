# =============================================================================
# 02_plot_figure3.R
# Coverage scenario comparison
# Reads pre-computed CSV from output_figgen/figure_3_run_summary.csv
<<<<<<< HEAD
# and output_figgen/figure_3_weekly_hcw_deaths_80.csv
=======
>>>>>>> 873ecc12b709e76c5085cd6ebf2f57c289f1da8c
# Run 02_extract_figure3.R first.
# =============================================================================
source(here::here("analyses", "03_figure_template", "helper_functions_figure_1to4.R"))
OUT_DIR <- here("figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FIG3_EFFICACY_LEVELS <- c("obv_50", "obv_60", "obv_70", "obv_80", "obv_90")
FIG3_EFFICACY_LABELS <- c("50%", "60%", "70%", "80%", "90%")

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

# Coverage curve panels
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

<<<<<<< HEAD
x_max_weeks <- function(sc) SCENARIO_X_MAX_DAYS[sc] / 7

=======
>>>>>>> 873ecc12b709e76c5085cd6ebf2f57c289f1da8c
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

<<<<<<< HEAD
message("Figure 3 variants saved")

# =============================================================================
# Weekly HCW deaths incidence overlay (80% efficacy):
# baseline (no antiviral) + full / ramp_high / ramp_low, per scenario
# =============================================================================
weekly_80 <- read.csv(here("output_figgen", "figure_3_weekly_hcw_deaths_80.csv"),
                      stringsAsFactors = FALSE)

LINE_LEVELS <- c("baseline", "full", "ramp_high", "ramp_low")
LINE_LABELS <- c("No antiviral", "Full (100%)", "Ramp high (0%->80%)", "Ramp low (0%->50%)")
LINE_COLORS <- c(baseline = "grey40", full = "#1a9641",
                 ramp_high = "#fdae61", ramp_low = "#d7191c")

make_weekly_deaths_panel <- function(sc) {
  x_max <- x_max_weeks(sc)
  df <- weekly_80 %>%
    filter(scenario == sc, week <= x_max) %>%
    mutate(line_group = factor(line_group, levels = LINE_LEVELS))
  
  ggplot(df, aes(x = week, y = q50, color = line_group, fill = line_group)) +
    # geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.15, color = NA) +
    geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = LINE_COLORS,
                       breaks = LINE_LEVELS, labels = LINE_LABELS, name = NULL) +
    scale_fill_manual(values = LINE_COLORS,
                      breaks = LINE_LEVELS, labels = LINE_LABELS, name = NULL) +
    scale_x_continuous(limits = c(0, x_max), breaks = seq(0, x_max, 13)) +
    labs(x = "Weeks since outbreak start", y = "Incident HCW deaths") +
    theme_fig()
}

p_weekly_d <- make_weekly_deaths_panel("WestAfrica")
p_weekly_e <- make_weekly_deaths_panel("DRC")

fig3_weekly_deaths <- (p_weekly_d | p_weekly_e) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = list(c("d ", "e "))) &
  theme(legend.position = "bottom")

save_fig("figure_3_weekly-hcw-deaths-80pct", fig3_weekly_deaths, 10, 4)

message("Figure 3 weekly HCW deaths overlay (80% efficacy) saved")



# =============================================================================
# Combined figure: a-c coverage curves, d-e weekly HCW deaths overlay (80%),
# f-h deaths averted boxplots
# =============================================================================
p_f_full      <- make_box_plot(COVERAGE_LEVELS[1], "pct_hcw_deaths_averted", "HCW deaths averted")
p_f_ramp_high <- make_box_plot(COVERAGE_LEVELS[2], "pct_hcw_deaths_averted", "HCW deaths averted")
p_f_ramp_low  <- make_box_plot(COVERAGE_LEVELS[3], "pct_hcw_deaths_averted", "HCW deaths averted")

figure_3_combined <- (
  (h1 | h2 | h3) /
    ((p_a | p_b | p_c) + plot_layout(axis_titles = "collect")) /
    ((p_weekly_d | p_weekly_e) + plot_layout(axis_titles = "collect", widths = c(1, 1))) /
    ((p_f_full | p_f_ramp_high | p_f_ramp_low) + plot_layout(axis_titles = "collect"))
) +
  plot_layout(guides = "collect", heights = c(0.2, 1, 1, 1.5)) +
  plot_annotation(tag_levels = list(c("", "", "", "a ", "b ", "c ", "d ", "e ", "f ", "g ", "h "))) &
  theme(legend.position = "bottom")

save_fig("figure_3_combined", figure_3_combined, 10, 9)


message("Figure 3 combined (a-h) saved")





=======
message("Figure 3 variants saved")
>>>>>>> 873ecc12b709e76c5085cd6ebf2f57c289f1da8c
