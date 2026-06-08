# =============================================================================
# 02_plot_figure3.R
# Coverage scenario comparison -- post-hoc applied to baseline runs
# =============================================================================
source(here::here("analyses", "03_figure_template", "helper_functions_figure_1to4.R"))
OUT_DIR <- here("figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

results <- load_results()

# =============================================================================
# Panels a, b, c: Coverage curve visualisations
# =============================================================================
make_coverage_plot <- function(cs) {
  spec  <- COVERAGE_SPECS[[cs]]
  df    <- data.frame(time = spec$times, coverage = spec$values * 100)
  color <- COVERAGE_COLORS[cs]

  ggplot(df, aes(x = time, y = coverage)) +
    geom_line(color = color, linewidth = 1.2) +
    geom_point(color = color, size = 3) +
    scale_y_continuous(limits = c(0, 100)) +
    scale_x_continuous(breaks = c(0, 30, 60, 90)) +
    labs(x = "Days since outbreak start",
         y = "OBV coverage (%)") +
    theme_fig()
}

ggsave(file.path(OUT_DIR, "figure_3_a.png"), make_coverage_plot("full"),
       width = 5, height = 4, dpi = 150)
ggsave(file.path(OUT_DIR, "figure_3_b.png"), make_coverage_plot("ramp_high"),
       width = 5, height = 4, dpi = 150)
ggsave(file.path(OUT_DIR, "figure_3_c.png"), make_coverage_plot("ramp_low"),
       width = 5, height = 4, dpi = 150)

# =============================================================================
# Build post-hoc OBV data: all efficacy x coverage combinations in one pass
# =============================================================================
message("Applying post-hoc OBV across efficacy x coverage grid...")

# Baseline rows needed once for the denominator
baseline_rows <- build_run_df_obv(results, "baseline")

obv_rows <- do.call(rbind, lapply(OBV_EFFICACY_LEVELS, function(eff) {
  do.call(rbind, lapply(COVERAGE_LEVELS, function(cov) {
    build_run_df_obv(results, eff, cov)
  }))
}))

all_rows <- rbind(baseline_rows, obv_rows)
pdf      <- make_particle_df(all_rows)

obv_df <- pdf %>%
  filter(arm %in% OBV_EFFICACY_LEVELS) %>%
  mutate(
    arm_label      = factor(OBV_EFFICACY_LABELS[match(arm, OBV_EFFICACY_LEVELS)],
                            levels = OBV_EFFICACY_LABELS),
    coverage_label = factor(COVERAGE_LABELS[match(coverage_scenario, COVERAGE_LEVELS)],
                            levels = COVERAGE_LABELS),
    scenario_label = factor(SCENARIO_LABELS[scenario], levels = SCENARIO_LABELS)
  )

# =============================================================================
# Panels d-f (HCW deaths averted) and g-i (HCW days lost averted)
# x = efficacy, colour = scenario, one panel per coverage scenario
# =============================================================================
sc_colors <- setNames(SCENARIO_COLORS, SCENARIO_LABELS)

make_box_plot <- function(cs, metric, y_label) {
  df <- obv_df %>%
    filter(coverage_scenario == cs, !is.na(.data[[metric]])) %>%
    mutate(fill_group = paste(as.character(scenario_label),
                              as.character(arm_label), sep = "."))

  build_sc_fill <- function(sc_key) {
    sc_lbl    <- SCENARIO_LABELS[sc_key]
    base_col  <- SCENARIO_COLORS[sc_key]
    light_col <- if (sc_key == "WestAfrica") "#fdd8a0" else "#b2e4d8"
    dark_col  <- rgb(t(col2rgb(base_col) * 0.7), maxColorValue = 255)
    cols      <- c(colorRampPalette(c(light_col, base_col))(4), dark_col)
    setNames(cols, paste(sc_lbl, OBV_EFFICACY_LABELS, sep = "."))
  }

  fill_vals <- c(build_sc_fill("WestAfrica"), build_sc_fill("DRC"))
  legend_breaks <- c(
    paste(SCENARIO_LABELS["WestAfrica"], "80%", sep = "."),
    paste(SCENARIO_LABELS["DRC"],        "80%", sep = ".")
  )

  ggplot(df, aes(x = arm_label, y = .data[[metric]], fill = fill_group)) +
    geom_boxplot(aes(linewidth = arm_label == "80%"),
                 outlier.size = 0.5, width = 0.6, color = "black",
                 position = position_dodge(0.75)) +
    scale_fill_manual(values = fill_vals, breaks = legend_breaks,
                      labels = c("West Africa", "DRC"), name = NULL) +
    scale_linewidth_manual(values = c("TRUE" = 1.0, "FALSE" = 0.4), guide = "none") +
    scale_y_continuous(limits = c(0, 100)) +
    labs(x = "OBV efficacy", y = y_label) +
    theme_fig() +
    theme(panel.grid.major.x = element_blank())
}

# Panels d, e, f -- HCW deaths averted
for (i in seq_along(COVERAGE_LEVELS)) {
  cs    <- COVERAGE_LEVELS[i]
  label <- letters[3 + i]   # d, e, f
  ggsave(
    file.path(OUT_DIR, sprintf("figure_3_%s.png", label)),
    make_box_plot(cs, "pct_hcw_deaths_averted",
                  "HCW deaths averted (%)"),
    width = 7, height = 5, dpi = 150
  )
}

# =============================================================================
# Composite Figure 3: row 1 = coverage curves (a-c), row 2 = HCW deaths (d-f)
# =============================================================================
library(patchwork)

make_col_header <- function(label) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = label, fontface = "bold", size = 4.5) +
    theme_void()
}

h1 <- make_col_header("Constant, Full Coverage")
h2 <- make_col_header("Ramp Up to High Coverage")
h3 <- make_col_header("Ramp Up to Medium Coverage")

# Re-generate panel objects (individual PNGs already saved above)
p_a <- make_coverage_plot(COVERAGE_LEVELS[1])
p_b <- make_coverage_plot(COVERAGE_LEVELS[2])
p_c <- make_coverage_plot(COVERAGE_LEVELS[3])

p_d <- make_box_plot(COVERAGE_LEVELS[1], "pct_hcw_deaths_averted",
                     "HCW deaths averted (%)")
p_e <- make_box_plot(COVERAGE_LEVELS[2], "pct_hcw_deaths_averted",
                     "HCW deaths averted (%)")
p_f <- make_box_plot(COVERAGE_LEVELS[3], "pct_hcw_deaths_averted",
                     "HCW deaths averted (%)")

no_x  <- function(p) p + theme(axis.title.x = element_blank())
no_y  <- function(p) p + theme(axis.title.y = element_blank())
no_xy <- function(p) p + theme(axis.title.x = element_blank(), axis.title.y = element_blank())

figure_3_composite <- (
  (h1 | h2 | h3) /
  (no_x(p_a) | no_y(p_b) | no_xy(p_c)) /
  (no_x(p_d) | no_y(p_e) | no_xy(p_f))
) +
  plot_layout(guides = "collect", axes = "collect",
              heights = c(0.08, 1, 1)) +
  plot_annotation(tag_levels = list(c("", "", "", "a", "b", "c", "d", "e", "f"))) &
  theme(legend.position = "bottom")

ggsave(
  file.path(OUT_DIR, "figure_3_ALL.png"),
  figure_3_composite,
  width = 15, height = 11, dpi = 150
)

#message("Composite Figure 3 saved to ", file.path(OUT_DIR, "figure_3_composite.png"))
# Panels g, h, i -- HCW days lost averted
for (i in seq_along(COVERAGE_LEVELS)) {
  cs    <- COVERAGE_LEVELS[i]
  label <- letters[6 + i]   # g, h, i
  ggsave(
    file.path(OUT_DIR, sprintf("figure_3_%s.png", label)),
    make_box_plot(cs, "pct_days_lost_averted",
                  "HCW days lost averted (%)"),
    width = 7, height = 5, dpi = 150
  )
}


message("Figure 3 panels saved: a-c (coverage curves), d-f (deaths averted), g-i (days lost averted)")
