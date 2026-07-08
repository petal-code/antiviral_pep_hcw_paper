# =============================================================================
# 02_plot_figure3.R
# Coverage scenario comparison
# Reads pre-computed CSV from output_figgen/figure_3_run_summary.csv
# and output_figgen/figure_3_weekly_hcw_deaths_80.csv
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
  ggsave(file.path(OUT_DIR, paste0(filename_base, ".pdf")),
         plot, width = width, height = height, units = "in")
}

x_max_weeks <- function(sc) SCENARIO_X_MAX_DAYS[sc] / 7
h1 <- make_col_header("Scenario 1: Outbreak-ready /\npre-positioned PEP")
h2 <- make_col_header("Scenario 2: Deployment-ready /\nsupply-constrained PEP")
h3 <- make_col_header("Scenario 3: Evaluation-dependent /\nlimited-access PEP")

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
# Weekly HCW deaths incidence overlay (80% efficacy):
# baseline (no antiviral) + full / ramp_high / ramp_low, per scenario
# =============================================================================
weekly_80 <- read.csv(here("output_figgen", "figure_3_weekly_hcw_deaths_80.csv"),
                      stringsAsFactors = FALSE)

LINE_LEVELS <- c("baseline", "full", "ramp_high", "ramp_low")
LINE_LABELS <- c("No antiviral", "Scenario 1", "Scenario 2", "Scenario 3")
# LINE_LABELS <- c("No antiviral", "Full (100%)", "Ramp high (0%->80%)", "Ramp low (0%->50%)")
# LINE_COLORS <- c(baseline = "grey40", full = "#1a9641",
#                  ramp_high = "#fdae61", ramp_low = "#d7191c")
LINE_COLORS <- c(baseline = "grey40", full = "#5C995C", ramp_high = "#BCA05C", ramp_low = "#B3614C")  # sage / ochre / terracotta

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

# =============================================================================
# Figure 3 (redesign): coverage curves (A,B,C) on top; per-country weekly HCW-
# death incidence (D,E) on the LEFT; % HCW deaths averted vs efficacy +95% CrI
# (F,G) on the RIGHT.   AAABBBCCC / DDDDDDFFF / EEEEEEGGG
# =============================================================================
source(here::here("analyses", "03_figure_template", "helper_functions_figure_1to4.R"))
library(dplyr); library(ggplot2); library(patchwork)
OUT_DIR <- here::here("figures"); dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FIG3_EFFICACY_LEVELS <- c("obv_50", "obv_60", "obv_70", "obv_80", "obv_90")

# ---- coverage-scenario colours (Scenario 1 -> 3) -- pick ONE -----------------
# Country colours are orange (#d95f02) / teal (#1b9e77), so these AVOID those hues.
# Used consistently for the coverage curves, the D/E overlay lines and the F/G
# lines; baseline stays grey.
# COV_PAL <- c(full = "#2B3A67", ramp_high = "#496A81", ramp_low = "#66999B")  # A (default): blue -> purple -> red
# COV_PAL <- c(full = "#0072B2", ramp_high = "#CC79A7", ramp_low = "#E69F00") # B: Okabe-Ito (colourblind-safe)
# COV_PAL <- c(full = "#08519c", ramp_high = "#3182bd", ramp_low = "#9ecae1") # C: sequential blues (dark = most ready)
# COV_PAL <- c(full = "#1b7837", ramp_high = "#9970ab", ramp_low = "#762a83") # D: green -> purple (PRGn)
# COV_PAL <- c(full = "#1a9641", ramp_high = "#fdae61", ramp_low = "#d7191c")
# COV_PAL <- c(full = "#1F7039", ramp_high = "#CD945A", ramp_low = "#AD2628")
# COV_PAL <- c(full = "#428A4E", ramp_high = "#C2A347", ramp_low = "#B25243")  # muted forest / gold / brick
COV_PAL <- c(full = "#5C995C", ramp_high = "#BCA05C", ramp_low = "#B3614C")  # sage / ochre / terracotta

COVERAGE_COLORS <- COV_PAL
LINE_COLORS     <- c(baseline = "grey40", COV_PAL)

# ---- data -------------------------------------------------------------------
run_df <- read.csv(here::here("output_figgen", "figure_3_run_summary.csv"), stringsAsFactors = FALSE)
pdf <- make_particle_df(run_df) %>%
  filter(arm != "baseline") %>%
  mutate(coverage_name  = sub("__.*", "", arm),
         eff_name       = sub(".*__", "", arm),
         efficacy_pct   = as.numeric(sub("obv_", "", eff_name)),
         coverage_label = factor(COVERAGE_LABELS[match(coverage_name, COVERAGE_LEVELS)], levels = COVERAGE_LABELS))

weekly_80   <- read.csv(here::here("output_figgen", "figure_3_weekly_hcw_deaths_80.csv"), stringsAsFactors = FALSE)
x_max_weeks <- function(sc) SCENARIO_X_MAX_DAYS[sc] / 7

# ---- A,B,C: coverage curves (two-line scenario title) -----------------------
make_coverage_plot <- function(cs, title) {
  x_max <- max(SCENARIO_X_MAX_DAYS); t_days <- seq(0, x_max, by = 1)
  ggplot(data.frame(week = t_days / 7, coverage = coverage_at_time(t_days, COVERAGE_SPECS[[cs]]) * 100),
         aes(week, coverage)) +
    geom_line(color = COVERAGE_COLORS[cs], linewidth = 1.2) +
    scale_y_continuous(limits = c(0, 100), labels = function(x) paste0(x, "%")) +
    scale_x_continuous(breaks = seq(0, x_max / 7, by = 13)) +
    labs(x = "Weeks since outbreak start", y = "Antiviral coverage", title = title) +
    theme_fig() +
    theme(plot.title = element_text(face = "bold", size = 9, hjust = 0.5, lineheight = 0.95))
}

# ---- D,E: weekly incident HCW deaths overlay (80% efficacy); country title ---
LINE_LEVELS <- c("baseline", "full", "ramp_high", "ramp_low")
LINE_LABELS <- c("No antiviral", "Scenario 1", "Scenario 2", "Scenario 3")
make_weekly_deaths_panel <- function(sc, panel_title) {
  xm <- x_max_weeks(sc)
  df <- weekly_80 %>% filter(scenario == sc, week <= xm) %>%
    mutate(line_group = factor(line_group, levels = LINE_LEVELS))
  ggplot(df, aes(week, q50, color = line_group, fill = line_group)) +
    geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = LINE_COLORS, breaks = LINE_LEVELS, labels = LINE_LABELS, name = NULL) +
    scale_fill_manual(values = LINE_COLORS, breaks = LINE_LEVELS, labels = LINE_LABELS, name = NULL) +
    scale_x_continuous(limits = c(0, xm), breaks = seq(0, xm, 13)) +
    labs(x = "Weeks since outbreak start", y = "Incident HCW deaths", title = panel_title) +
    theme_fig() +
    theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5))
}

# ---- F,G: % HCW deaths averted vs efficacy (median + 95% CrI); NO title ------
make_averted_line_panel <- function(sc, metric = "pct_hcw_deaths_averted", y_label = "HCW deaths averted") {
  df <- pdf %>%
    filter(scenario == sc, eff_name %in% FIG3_EFFICACY_LEVELS, !is.na(.data[[metric]])) %>%
    group_by(coverage_name, coverage_label, efficacy_pct) %>%
    summarise(lo  = quantile(.data[[metric]], 0.025, names = FALSE),
              med = quantile(.data[[metric]], 0.5,   names = FALSE),
              hi  = quantile(.data[[metric]], 0.975, names = FALSE), .groups = "drop")
  ggplot(df, aes(efficacy_pct, med, color = coverage_label, fill = coverage_label)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1) + geom_point(size = 1.6) +
    scale_color_manual(values = setNames(COVERAGE_COLORS[COVERAGE_LEVELS], COVERAGE_LABELS), name = NULL) +
    scale_fill_manual(values  = setNames(COVERAGE_COLORS[COVERAGE_LEVELS], COVERAGE_LABELS), name = NULL) +
    scale_x_continuous(breaks = c(50, 60, 70, 80, 90), labels = function(x) paste0(x, "%")) +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    guides(color = "none", fill = "none") +              # F/G share D/E's legend
    labs(x = "Antiviral efficacy", y = y_label) +
    lims(y = c(0, NA)) +
    theme_fig()
}

# ---- country row order (top -> bottom). Per "respectively" = West Africa, DRC.
#      SWAP these two lines if you want DRC on top. ----------------------------
TOP_SC <- "WestAfrica"; TOP_LAB <- "West Africa Archetype"
BOT_SC <- "DRC";        BOT_LAB <- "DRC Archetype"

p_a <- make_coverage_plot("full",      "Scenario 1: Outbreak-ready /\npre-positioned PEP stockpiles")
p_b <- make_coverage_plot("ramp_high", "Scenario 2: Deployment-ready /\nsupply-constrained PEP")
p_c <- make_coverage_plot("ramp_low",  "Scenario 3: Evaluation-dependent /\nlimited-access PEP")
p_d <- make_weekly_deaths_panel(TOP_SC, TOP_LAB)   # D (top-left)
p_e <- make_weekly_deaths_panel(BOT_SC, BOT_LAB)   # E (bottom-left)
p_f <- make_averted_line_panel(TOP_SC)             # F (top-right, no title)
p_g <- make_averted_line_panel(BOT_SC)             # G (bottom-right, no title)

design <- "
AAABBBCCC
AAABBBCCC
DDDDDDFFF
DDDDDDFFF
DDDDDDFFF
EEEEEEGGG
EEEEEEGGG
EEEEEEGGG
"

figure_3_redesign <-
  p_a + p_b + p_c + p_d + p_e + p_f + p_g +
  plot_layout(design = design, guides = "collect") +     # row heights come from the design
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom",
        plot.tag = element_text(face = "bold"))

print(figure_3_redesign)
ggsave(file.path(OUT_DIR, "figure_3_redesign.pdf"), 
       figure_3_redesign, width = 8, height = 6.5)
# 6.5 * 8