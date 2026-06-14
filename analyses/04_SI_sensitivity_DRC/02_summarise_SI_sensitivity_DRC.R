# =============================================================================
# 02_summarise_SI_sensitivity_DRC.R
#
# Turns the per-particle output of 01_run_SI_sensitivity_DRC.R into the SI
# deliverables for the two DRC-like sensitivity analyses:
#
#   (1) Transmissibility stress test (R0 / offspring means +10/+20/+30%)
#       Reports: no-PEP baseline HCW deaths, HCW deaths averted, % reduction.
#   (2) HCW-exposure upshift (hcw_risk_scalar +25/+50/+100%)
#       Reports: baseline HCW deaths, HCW deaths averted, HCW-days lost averted,
#                % reduction.
#
# Stress-level labels are derived from the scaling factors set in script 1, so
# they always match whatever factors are configured there.
#
# Posterior uncertainty is summarised as median + 95% credible interval
# (2.5/97.5% quantiles across particles) -- the same convention as the main
# figures (cf. 02_plot_figure2.R make_summ()). The as-fitted DRC archetype
# (x1.00) is included in each analysis as the reference level.
#
# Inputs : output_figgen/SI_sensitivity_DRC_particle_df.csv
# Outputs:
#   output_figgen/SI_sensitivity_DRC_summary_long.csv   (tidy: median/lo95/hi95)
#   output_figgen/SI_sensitivity_DRC_summary_table.csv  (formatted SI table)
#   figures/figure_S_DRC_transmissibility_stress_test.{pdf,png}
#   figures/figure_S_DRC_hcw_exposure_upshift.{pdf,png}
# =============================================================================

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

IN_CSV  <- here("output_figgen", "SI_sensitivity_DRC_particle_df.csv")
FIG_DIR <- here("figures")
OUT_DIR <- here("output_figgen")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

particle_df <- read.csv(IN_CSV, stringsAsFactors = FALSE)

DRC_COL   <- "#1b9e77"
EFF_LABEL <- function(eff) sprintf("%d%%", round(eff * 100))

# Light->dark green ramp, one shade per stress level (fitted = lightest).
level_palette <- function(n) colorRampPalette(c("#b2e4d8", DRC_COL, "#0b4f3c"))(n)

# -----------------------------------------------------------------------------
# Attach the shared reference (x1.00) cell to each analysis and label the levels.
# Labels are DERIVED from the scaling factor (factor 1 -> "+0% (fitted)", factor
# f -> "+<round((f-1)*100)>%") so they can never disagree with the factors set in
# 01_run_SI_sensitivity_DRC.R, whatever those happen to be.
# -----------------------------------------------------------------------------
level_label_of <- function(f) {
  ifelse(abs(f - 1) < 1e-9, "+0% (fitted)", sprintf("+%g%%", round((f - 1) * 100)))
}

make_analysis_df <- function(which_analysis, factor_col) {
  df <- particle_df %>%
    filter(analysis %in% c("reference", which_analysis)) %>%
    mutate(level_factor = .data[[factor_col]],
           level_label  = level_label_of(.data[[factor_col]]))
  ord <- sort(unique(df$level_factor))
  df$level_label <- factor(df$level_label, levels = level_label_of(ord))
  df$efficacy_label <- factor(EFF_LABEL(df$efficacy),
                              levels = EFF_LABEL(sort(unique(df$efficacy))))
  df
}

trans_df <- make_analysis_df("transmissibility", "r0_factor")
hcw_df   <- make_analysis_df("hcw_exposure",     "hcw_factor")

# -----------------------------------------------------------------------------
# Across-particle summary: median + 95% credible interval.
# -----------------------------------------------------------------------------
METRICS <- c("baseline_hcw_deaths", "pep_hcw_deaths", "hcw_deaths_averted",
             "pct_hcw_deaths_averted", "baseline_hcw_days_lost",
             "hcw_days_lost_averted", "pct_days_lost_averted")

summarise_levels <- function(df, analysis_name) {
  do.call(rbind, lapply(METRICS, function(m) {
    df %>%
      group_by(level_label, level_factor, efficacy, efficacy_label) %>%
      summarise(
        median = median(.data[[m]], na.rm = TRUE),
        lo95   = quantile(.data[[m]], 0.025, na.rm = TRUE),
        hi95   = quantile(.data[[m]], 0.975, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(analysis = analysis_name, metric = m)
  }))
}

summary_long <- bind_rows(
  summarise_levels(trans_df, "transmissibility"),
  summarise_levels(hcw_df,   "hcw_exposure")
) %>%
  select(analysis, level_label, level_factor, efficacy, efficacy_label,
         metric, median, lo95, hi95) %>%
  arrange(analysis, metric, level_factor, efficacy)

write.csv(summary_long, file.path(OUT_DIR, "SI_sensitivity_DRC_summary_long.csv"),
          row.names = FALSE)
message("Saved tidy summary: SI_sensitivity_DRC_summary_long.csv")

# -----------------------------------------------------------------------------
# Formatted SI table: one row per (analysis x level x efficacy), reporting
# median [95% CrI] for each requested quantity.
# -----------------------------------------------------------------------------
fmt_count <- function(m, l, h) sprintf("%.0f [%.0f, %.0f]", m, l, h)
fmt_pct   <- function(m, l, h) sprintf("%.1f [%.1f, %.1f]", m, l, h)

fmt_metric <- function(df_long, analysis_name, metric_name, formatter) {
  df_long %>%
    filter(analysis == analysis_name, metric == metric_name) %>%
    transmute(level_label, efficacy_label,
              value = formatter(median, lo95, hi95))
}

build_table <- function(analysis_name, include_days) {
  base <- fmt_metric(summary_long, analysis_name, "baseline_hcw_deaths",    fmt_count) %>%
    rename(`Baseline HCW deaths (no PEP)` = value)
  avert <- fmt_metric(summary_long, analysis_name, "hcw_deaths_averted",     fmt_count) %>%
    rename(`HCW deaths averted` = value)
  pct <- fmt_metric(summary_long, analysis_name, "pct_hcw_deaths_averted", fmt_pct) %>%
    rename(`HCW deaths averted (%)` = value)
  tab <- base %>% left_join(avert, by = c("level_label", "efficacy_label")) %>%
    left_join(pct, by = c("level_label", "efficacy_label"))
  if (include_days) {
    days <- fmt_metric(summary_long, analysis_name, "hcw_days_lost_averted", fmt_count) %>%
      rename(`HCW-days lost averted` = value)
    pdays <- fmt_metric(summary_long, analysis_name, "pct_days_lost_averted", fmt_pct) %>%
      rename(`HCW-days lost averted (%)` = value)
    tab <- tab %>% left_join(days, by = c("level_label", "efficacy_label")) %>%
      left_join(pdays, by = c("level_label", "efficacy_label"))
  }
  tab %>%
    mutate(analysis = analysis_name) %>%
    rename(`Stress level` = level_label, `Antiviral efficacy` = efficacy_label) %>%
    relocate(analysis)
}

summary_table <- bind_rows(
  build_table("transmissibility", include_days = TRUE),
  build_table("hcw_exposure",     include_days = TRUE)
)
write.csv(summary_table, file.path(OUT_DIR, "SI_sensitivity_DRC_summary_table.csv"),
          row.names = FALSE)
message("Saved formatted SI table: SI_sensitivity_DRC_summary_table.csv")

# Echo the headline arm (full coverage, 80% efficacy) to the console.
message("\n--- Headline arm: full coverage, 80% antiviral efficacy ---")
summary_table %>%
  filter(`Antiviral efficacy` == "80%") %>%
  as.data.frame() %>%
  print()

# -----------------------------------------------------------------------------
# Figures. Median points + 95% CrI error bars, one colour per stress level.
# -----------------------------------------------------------------------------
theme_si <- function() theme_classic(base_size = 11) +
  theme(legend.position = "top", panel.grid.major.y = element_line(colour = "grey92"))

plot_metric <- function(df_long, analysis_name, metric_name, y_label,
                        levels_in_order, pct_axis = FALSE) {
  d <- df_long %>% filter(analysis == analysis_name, metric == metric_name)
  d$level_label <- factor(d$level_label, levels = levels_in_order)
  pal <- setNames(level_palette(length(levels_in_order)), levels_in_order)
  p <- ggplot(d, aes(efficacy_label, median, colour = level_label, group = level_label)) +
    geom_line(position = position_dodge(width = 0.5), linewidth = 0.6, alpha = 0.6) +
    geom_errorbar(aes(ymin = lo95, ymax = hi95), width = 0.35,
                  position = position_dodge(width = 0.5), linewidth = 0.6) +
    geom_point(position = position_dodge(width = 0.5), size = 2) +
    scale_colour_manual(values = pal, name = NULL) +
    labs(x = "Antiviral efficacy", y = y_label) +
    theme_si()
  if (pct_axis)
    p <- p + scale_y_continuous(labels = function(x) paste0(x, "%"))
  p
}

# Baseline (no-PEP) HCW deaths by stress level -- collapse efficacy (the no-PEP
# baseline is independent of the PEP arm), summarising across particles x arms.
baseline_panel <- function(df, analysis_name, levels_in_order) {
  d <- df %>%
    group_by(level_label) %>%
    summarise(median = median(baseline_hcw_deaths, na.rm = TRUE),
              lo95   = quantile(baseline_hcw_deaths, 0.025, na.rm = TRUE),
              hi95   = quantile(baseline_hcw_deaths, 0.975, na.rm = TRUE),
              .groups = "drop")
  d$level_label <- factor(d$level_label, levels = levels_in_order)
  pal <- setNames(level_palette(length(levels_in_order)), levels_in_order)
  ggplot(d, aes(level_label, median, fill = level_label)) +
    geom_col(width = 0.7, alpha = 0.9) +
    geom_errorbar(aes(ymin = lo95, ymax = hi95), width = 0.2, linewidth = 0.6) +
    scale_fill_manual(values = pal, guide = "none") +
    labs(x = NULL, y = "Baseline HCW deaths (no PEP)") +
    theme_si()
}

# ---- Figure: transmissibility stress test ----
trans_levels <- levels(trans_df$level_label)
fig_trans <- (
  baseline_panel(trans_df, "transmissibility", trans_levels) +
    plot_metric(summary_long, "transmissibility", "hcw_deaths_averted",
                "HCW deaths averted (count)", trans_levels) +
    plot_metric(summary_long, "transmissibility", "pct_hcw_deaths_averted",
                "HCW deaths averted (%)", trans_levels, pct_axis = TRUE)
) +
  plot_layout(nrow = 1, guides = "collect") +
  plot_annotation(
    title = "DRC-like vaccine-free stress test: baseline transmissibility +10/+20/+30%",
    subtitle = "HCW-targeted PEP at 100% coverage. Median +/- 95% CrI across posterior particles.",
    tag_levels = "a"
  ) & theme(legend.position = "top")
ggsave(file.path(FIG_DIR, "figure_S_DRC_transmissibility_stress_test.pdf"),
       fig_trans, width = 12, height = 4.5)
ggsave(file.path(FIG_DIR, "figure_S_DRC_transmissibility_stress_test.png"),
       fig_trans, width = 12, height = 4.5, dpi = 320)
message("Saved figure: figure_S_DRC_transmissibility_stress_test.{pdf,png}")

# ---- Figure: HCW-exposure upshift ----
hcw_levels <- levels(hcw_df$level_label)
fig_hcw <- (
  baseline_panel(hcw_df, "hcw_exposure", hcw_levels) +
    plot_metric(summary_long, "hcw_exposure", "pct_hcw_deaths_averted",
                "HCW deaths averted (%)", hcw_levels, pct_axis = TRUE) +
    plot_metric(summary_long, "hcw_exposure", "pct_days_lost_averted",
                "HCW-days lost averted (%)", hcw_levels, pct_axis = TRUE)
) +
  plot_layout(nrow = 1, guides = "collect") +
  plot_annotation(
    title = "DRC-like HCW-exposure upshift: hcw_risk_scalar +25/+50/+100%",
    subtitle = "HCW-targeted PEP at 100% coverage. Median +/- 95% CrI across posterior particles.",
    tag_levels = "a"
  ) & theme(legend.position = "top")
ggsave(file.path(FIG_DIR, "figure_S_DRC_hcw_exposure_upshift.pdf"),
       fig_hcw, width = 12, height = 4.5)
ggsave(file.path(FIG_DIR, "figure_S_DRC_hcw_exposure_upshift.png"),
       fig_hcw, width = 12, height = 4.5, dpi = 320)
message("Saved figure: figure_S_DRC_hcw_exposure_upshift.{pdf,png}")

message("\nSI sensitivity summarising complete.")
