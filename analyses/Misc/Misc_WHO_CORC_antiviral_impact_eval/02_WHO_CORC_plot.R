# 02_WHO_CORC_plot.R  (analyses/Misc_WHO_CORC_antiviral_impact_eval)
# =============================================================================
# Obeldesivir (OBV PEP) impact analysis — plotting step.
#
# Reads the intermediate produced by 01_run_simulations.R and produces:
#
#   Epidemic curves (deaths over time AND HCW deaths over time), where each
#   parameter set's N_REPS replicates are summarised by their MEAN trajectory:
#     (i)  every parameter set's mean trajectory plotted individually;
#     (ii) the median trajectory with the 25% / 75% interval band (computed
#          across the parameter sets' mean trajectories).
#   Both with and without obeldesivir are overlaid for direct comparison.
#
#   A bar chart of the % of HCW deaths averted by obeldesivir (with the % of all
#   deaths averted shown alongside for context); bars are medians across the
#   parameter sets, error bars are the 25%/75% interval.
#
# Figures are written to outputs/.
# =============================================================================


# -----------------------------------------------------------------------------
# 0. Locate the repo root (works on any machine, from anywhere in the repo)
# -----------------------------------------------------------------------------
# here::here() finds the repo root by locating antiviral_pep_hcw_paper.Rproj, so there are
# no per-user paths to maintain. Requires the `here` package:
# install.packages("here").
REPO_ROOT    <- here::here()
ANALYSIS_DIR <- here::here("analyses", "Misc_WHO_CORC_antiviral_impact_eval")
OUTPUT_DIR   <- here::here("outputs", "misc", "WHO_CORC_outputs")
RESULTS_RDS  <- file.path(OUTPUT_DIR, "WHO_CORC_prelim_antiviral_simulation_results.rds")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)


# -----------------------------------------------------------------------------
# 1. LIBRARIES + DATA
# -----------------------------------------------------------------------------
library(ggplot2)
# patchwork is used only for the optional combined panel; load it if present.
have_patchwork <- requireNamespace("patchwork", quietly = TRUE)

if (!file.exists(RESULTS_RDS)) {
  stop("Intermediate not found: ", RESULTS_RDS,
       "\nRun 01_run_simulations.R first.", call. = FALSE)
}
res          <- readRDS(RESULTS_RDS)
trajectories <- res$trajectories
per_set      <- res$per_set_averted
cfg          <- res$config


# -----------------------------------------------------------------------------
# 2. PRESENTATION HELPERS (arm labels, palette, bin-unit label)
# -----------------------------------------------------------------------------
arm_levels <- c("no_obv", "obv")
arm_labels <- c("Without obeldesivir",
                sprintf("With obeldesivir (%.0f%% efficacy)", 100 * cfg$obv$efficacy))
arm_pal    <- setNames(c("#555555", "#1B9E77"), arm_labels)
to_arm     <- function(x) factor(arm_labels[match(x, arm_levels)], levels = arm_labels)

trajectories$arm_f <- to_arm(trajectories$arm)

bw   <- cfg$bin_width_days
unit <- if (bw == 7) "week" else if (bw == 1) "day" else sprintf("%d days", bw)
ylab_deaths <- sprintf("Deaths per %s", unit)
ylab_hcw    <- sprintf("HCW deaths per %s", unit)

subtitle_txt <- sprintf(
  "%s scenario | %d posterior parameter sets x %d replicates | obeldesivir: %.0f%% efficacy, %.0f%% coverage, %.0f%% adherence",
  cfg$scenario_id, cfg$n_sets, cfg$n_reps,
  100 * cfg$obv$efficacy, 100 * cfg$obv$coverage, 100 * cfg$obv$adherence
)

base_theme <- theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        legend.title    = element_blank(),
        plot.subtitle   = element_text(size = 8, colour = "grey30"),
        panel.grid.minor = element_blank())


# -----------------------------------------------------------------------------
# 3. EPIDEMIC CURVES
# -----------------------------------------------------------------------------
# (i) Individual per-set mean trajectories.
plot_individual <- function(ycol, ylab, title) {
  ggplot(trajectories,
         aes(x = time_days, y = .data[[ycol]],
             group = interaction(set_id, arm_f), colour = arm_f)) +
    geom_line(alpha = 0.25, linewidth = 0.3) +
    scale_colour_manual(values = arm_pal) +
    labs(x = "Time (days since outbreak seeding)", y = ylab,
         title = title, subtitle = subtitle_txt) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, linewidth = 1))) +
    base_theme
}

# (ii) Median trajectory with 25%/75% band, across the per-set mean trajectories.
summarise_band <- function(ycol) {
  f  <- function(v) quantile(v, c(0.25, 0.5, 0.75), names = FALSE)
  ag <- aggregate(list(q = trajectories[[ycol]]),
                  by  = list(arm_f = trajectories$arm_f, time_days = trajectories$time_days),
                  FUN = f)
  q  <- ag$q   # matrix: columns are the 25% / 50% / 75% quantiles
  data.frame(arm_f     = factor(ag$arm_f, levels = arm_labels),
             time_days = ag$time_days,
             lo = q[, 1], med = q[, 2], hi = q[, 3])
}
plot_band <- function(ycol, ylab, title) {
  band <- summarise_band(ycol)
  ggplot(band, aes(x = time_days, colour = arm_f, fill = arm_f)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.25, colour = NA) +
    geom_line(aes(y = med), linewidth = 0.9) +
    scale_colour_manual(values = arm_pal) +
    scale_fill_manual(values = arm_pal) +
    labs(x = "Time (days since outbreak seeding)", y = ylab,
         title = title, subtitle = subtitle_txt) +
    base_theme
}

p_deaths_ind  <- plot_individual("deaths_per_bin",     ylab_deaths,
                                 "Deaths over time: individual posterior trajectories")
p_deaths_band <- plot_band      ("deaths_per_bin",     ylab_deaths,
                                 "Deaths over time: median and 25-75% interval")
p_hcw_ind     <- plot_individual("hcw_deaths_per_bin", ylab_hcw,
                                 "HCW deaths over time: individual posterior trajectories")
p_hcw_band    <- plot_band      ("hcw_deaths_per_bin", ylab_hcw,
                                 "HCW deaths over time: median and 25-75% interval")

cowplot::plot_grid(p_deaths_ind, p_deaths_band, p_hcw_ind, p_hcw_band,
                   nrow = 2)

ggsave(file.path(OUTPUT_DIR, "antiviral_deaths_over_time_individual.png"),
       p_deaths_ind, width = 9, height = 5.5, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "antiviral_deaths_over_time_median_iqr.png"),
       p_deaths_band, width = 9, height = 5.5, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "antiviral_hcw_deaths_over_time_individual.png"),
       p_hcw_ind, width = 9, height = 5.5, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "antiviral_hcw_deaths_over_time_median_iqr.png"),
       p_hcw_band, width = 9, height = 5.5, dpi = 300)

# Combined 2x2 panel (deaths on top, HCW deaths on bottom; individual left,
# median/IQR right). Optional: only built if the patchwork package is installed.
if (have_patchwork) {
  library(patchwork)
  combined <- (p_deaths_ind + p_deaths_band) / (p_hcw_ind + p_hcw_band) +
    plot_layout(guides = "collect") &
    theme(legend.position = "top")
  ggsave(file.path(OUTPUT_DIR, "antiviral_epidemic_curves_combined.png"),
         combined, width = 15, height = 11, dpi = 300)
} else {
  message("Note: install the 'patchwork' package to also get the combined 2x2 panel.")
}


# -----------------------------------------------------------------------------
# 4. BAR CHART: % OF HCW DEATHS AVERTED (with % of all deaths for context)
# -----------------------------------------------------------------------------
qs <- function(x, p) quantile(x[is.finite(x)], p)
bar_df <- data.frame(
  metric = factor(c("HCW deaths", "All deaths"),
                  levels = c("HCW deaths", "All deaths")),
  median = c(median(per_set$pct_hcw_deaths_averted, na.rm = TRUE),
             median(per_set$pct_deaths_averted,     na.rm = TRUE)),
  lo     = c(qs(per_set$pct_hcw_deaths_averted, 0.25),
             qs(per_set$pct_deaths_averted,     0.25)),
  hi     = c(qs(per_set$pct_hcw_deaths_averted, 0.75),
             qs(per_set$pct_deaths_averted,     0.75))
)

p_bar <- ggplot(subset(bar_df, metric == "HCW deaths"), aes(x = metric, y = median, fill = metric)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.18, linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f%%", median)), vjust = -0.6, hjust = 1.05, size = 4) +
  scale_fill_manual(values = c("HCW deaths" = "#1B9E77", "All deaths" = "#999999"),
                    guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "% of HCW deaths averted by obeldesivir",
       title = "",
       subtitle = "") + 
         #paste0(paste(strwrap(subtitle_txt, width = 70), collapse = "\n"),
          #               "\nBars: median across parameter sets; error bars: 25-75% interval")) +
  theme_bw(base_size = 12) +
  theme(plot.subtitle = element_text(size = 8, colour = "grey30"))
        #panel.grid.major.x = element_blank(),
        # panel.grid.minor = element_blank())

ggsave(file.path(OUTPUT_DIR, "antiviral_pct_hcw_deaths_averted_bar.png"),
       p_bar, width = 8, height = 5.8, dpi = 300)


# -----------------------------------------------------------------------------
# 5. DONE
# -----------------------------------------------------------------------------
message("Figures written to: ", OUTPUT_DIR)
figs <- c("antiviral_deaths_over_time_individual.png",
          "antiviral_deaths_over_time_median_iqr.png",
          "antiviral_hcw_deaths_over_time_individual.png",
          "antiviral_hcw_deaths_over_time_median_iqr.png",
          "antiviral_pct_hcw_deaths_averted_bar.png")
if (have_patchwork) figs <- c(figs, "antiviral_epidemic_curves_combined.png")
for (f in figs) message("  - ", f)

p_hcw_update <- p_hcw_band +
  theme_bw() + 
  labs(x = "Time Since Outbreak Start (Days)", y = "Healthcare Worker\n Deaths Per Week",
       subtitle = "", title = "") +
  lims(x = c(0, 500)) +
  theme(legend.position = "bottom",
        legend.title = element_blank())

cowplot::plot_grid(p_hcw_update, p_bar, labels = c("A", "B"),
                   align = "h", axis = "b", rel_widths = c(3.5, 1))

