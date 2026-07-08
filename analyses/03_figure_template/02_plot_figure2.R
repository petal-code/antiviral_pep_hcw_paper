# =============================================================================
# 02_plot_figure2.R  -- Figure-2 redesign (lines for weekly incidence,
#   dodged bars + error bars for averted outcomes)
#   D,E: weekly incidence = lines.
#   F,G: averted = DODGED bars + error bars.
#   Layout:  AAABBBCCC / DDDDDDFFF / EEEEEEGGG
# Prereqs: 02_extract_figure2.R has produced the two CSVs.
# =============================================================================
source(here::here("analyses", "03_figure_template", "helper_functions_figure_1to4.R"))
library(dplyr); library(ggplot2); library(patchwork)
OUT_DIR <- here::here("figures"); dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FIG2_EFFICACY_LEVELS <- c("obv_50", "obv_60", "obv_70", "obv_80", "obv_90")
EFF_NUM              <- c(50, 60, 70, 80, 90)

# ---- coverage colours (muted; distinct from country orange/teal) ------------
COV_PAL <- c(full = "#428A4E", ramp_high = "#C2A347", ramp_low = "#B25243")
COVERAGE_COLORS <- COV_PAL
LINE_COLORS <- c(baseline = "grey40", COV_PAL)
LINE_LEVELS <- c("baseline", "full", "ramp_high", "ramp_low")
LINE_LABELS <- c("No antiviral", "Scenario 1", "Scenario 2", "Scenario 3")

TOP_SC <- "WestAfrica"; TOP_LAB <- "West Africa Archetype"
BOT_SC <- "DRC";        BOT_LAB <- "DRC Archetype"

# ---- data -------------------------------------------------------------------
run_df <- read.csv(here::here("output_figgen", "figure_2_run_summary.csv"), stringsAsFactors = FALSE)
pdf <- make_particle_df(run_df) %>%
  filter(arm != "baseline") %>%
  mutate(coverage_name  = sub("__.*", "", arm),
         eff_name       = sub(".*__", "", arm),
         efficacy_pct   = as.numeric(sub("obv_", "", eff_name)),
         coverage_label = factor(COVERAGE_LABELS[match(coverage_name, COVERAGE_LEVELS)], levels = COVERAGE_LABELS))
weekly_80   <- read.csv(here::here("output_figgen", "figure_2_weekly_hcw_deaths_80.csv"), stringsAsFactors = FALSE)
x_max_weeks <- function(sc) SCENARIO_X_MAX_DAYS[sc] / 7
.cov_fill   <- setNames(COVERAGE_COLORS[COVERAGE_LEVELS], COVERAGE_LABELS)

# ---- A,B,C: coverage curves -------------------------------------------------
make_coverage_plot <- function(cs, title) {
  x_max <- max(SCENARIO_X_MAX_DAYS); t_days <- seq(0, x_max, by = 1)
  ggplot(data.frame(week = t_days/7, coverage = coverage_at_time(t_days, COVERAGE_SPECS[[cs]])*100), aes(week, coverage)) +
    geom_line(color = COVERAGE_COLORS[cs], linewidth = 1.2) +
    scale_y_continuous(limits = c(0,100), labels = function(x) paste0(x,"%")) +
    scale_x_continuous(breaks = seq(0, x_max/7, by = 13)) +
    labs(x = "Weeks since outbreak start", y = "Antiviral coverage", title = title) +
    theme_fig() + theme(plot.title = element_text(face="bold", size=9, hjust=0.5, lineheight=0.95))
}

# ---- D,E as LINES -----------------------------------------------------------
make_weekly_lines <- function(sc, panel_title) {
  xm <- x_max_weeks(sc)
  df <- weekly_80 %>% filter(scenario==sc, week<=xm) %>% mutate(line_group = factor(line_group, levels = LINE_LEVELS))
  ggplot(df, aes(week, q50, color=line_group, fill=line_group)) +
    geom_ribbon(aes(ymin=q25, ymax=q75), alpha=0.15, color=NA) + geom_line(linewidth=1) +
    scale_color_manual(values=LINE_COLORS, breaks=LINE_LEVELS, labels=LINE_LABELS, name=NULL) +
    scale_fill_manual(values=LINE_COLORS, breaks=LINE_LEVELS, labels=LINE_LABELS, name=NULL) +
    scale_x_continuous(limits=c(0,xm), breaks=seq(0,xm,13)) +
    labs(x="Weeks since outbreak start", y="Incident HCW deaths", title=panel_title) +
    theme_fig() + theme(plot.title=element_text(face="bold", size=11, hjust=0.5))
}

# ---- F,G shared summary -----------------------------------------------------
averted_summary <- function(sc, metric="pct_hcw_deaths_averted") {
  pdf %>% filter(scenario==sc, eff_name %in% FIG2_EFFICACY_LEVELS, !is.na(.data[[metric]])) %>%
    group_by(coverage_label, efficacy_pct) %>%
    summarise(lo=quantile(.data[[metric]],.025,names=FALSE),
              med=quantile(.data[[metric]],.5,names=FALSE),
              hi=quantile(.data[[metric]],.975,names=FALSE), .groups="drop") %>%
    mutate(coverage_label = factor(coverage_label, levels=COVERAGE_LABELS))
}

# ---- F,G as DODGED bars + error bars ----------------------------------------
make_averted_dodged <- function(sc, y_label="HCW deaths averted") {
  d <- averted_summary(sc); d$eff <- factor(d$efficacy_pct, levels=EFF_NUM, labels=paste0(EFF_NUM,"%"))
  pos <- position_dodge(width=0.8)
  ggplot(d, aes(eff, med, fill=coverage_label)) +
    geom_col(position=pos, width=0.7) +
    geom_errorbar(aes(ymin=lo, ymax=hi), position=pos, width=0.25, color="grey30", linewidth=0.4) +
    scale_fill_manual(values=.cov_fill, name=NULL) +
    scale_y_continuous(labels=function(x) paste0(x,"%")) +
    guides(fill="none") + labs(x="Antiviral efficacy", y=y_label) + theme_fig()
}

# ---- assembly ---------------------------------------------------------------
design <- "
AAAAAABBBBBBCCCCCC
AAAAAABBBBBBCCCCCC
DDDDDDDDDDFFFFFFFF
DDDDDDDDDDFFFFFFFF
DDDDDDDDDDFFFFFFFF
EEEEEEEEEEGGGGGGGG
EEEEEEEEEEGGGGGGGG
EEEEEEEEEEGGGGGGGG
"
TITLES <- c("Scenario 1: Outbreak-ready /\npre-positioned PEP stockpiles",
            "Scenario 2: Deployment-ready /\nsupply-constrained PEP",
            "Scenario 3: Evaluation-dependent /\nlimited-access PEP")

assemble <- function(de_fun, fg_fun) {
  p_a <- make_coverage_plot("full", TITLES[1]); p_b <- make_coverage_plot("ramp_high", TITLES[2])
  p_c <- make_coverage_plot("ramp_low", TITLES[3])
  p_d <- de_fun(TOP_SC, TOP_LAB); p_e <- de_fun(BOT_SC, BOT_LAB)
  p_f <- fg_fun(TOP_SC);          p_g <- fg_fun(BOT_SC)
  fig <- p_a + p_b + p_c + p_d + p_e + p_f + p_g +
    plot_layout(design = design, guides = "collect") +
    plot_annotation(tag_levels = "a") & theme(legend.position = "bottom",
                                              plot.tag = element_text(face = "bold"))
  print(fig); invisible(fig)
}

fig2 <- assemble(make_weekly_lines, make_averted_dodged)

ggsave(plot = fig2,
       filename = file.path(OUT_DIR, "figure_2_redesign_bars.pdf"),
       width = 11, height = 9)

ggsave(plot = fig2,
       filename = file.path(OUT_DIR, "figure_2_redesign_bars.tiff"), dpi = 400,
       width = 11, height = 9)
message("Saved figure_2_redesign_bars")
# =============================================================================
# Numbers for the Figure 2 text (80% efficacy, by coverage scenario)
# =============================================================================
fig2_text_numbers <- pdf %>%
  filter(eff_name == "obv_80") %>%
  group_by(scenario, coverage_name) %>%
  summarise(
    lo  = quantile(pct_hcw_deaths_averted, 0.025, na.rm = TRUE),
    med = quantile(pct_hcw_deaths_averted, 0.5,   na.rm = TRUE),
    hi  = quantile(pct_hcw_deaths_averted, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(c(lo, med, hi), ~round(.x, 1)))

print(fig2_text_numbers)

# Percent of Scenario 1 (full) benefit preserved under Scenario 2 / 3,
# based on median pct_hcw_deaths_averted
fig2_pct_preserved <- fig2_text_numbers %>%
  select(scenario, coverage_name, med) %>%
  tidyr::pivot_wider(names_from = coverage_name, values_from = med) %>%
  mutate(
    pct_preserved_ramp_high = round(100 * ramp_high / full, 1),
    pct_preserved_ramp_low  = round(100 * ramp_low  / full, 1)
  )

print(fig2_pct_preserved)