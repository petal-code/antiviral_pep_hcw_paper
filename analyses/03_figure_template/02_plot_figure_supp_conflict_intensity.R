# =============================================================================
# 02_plot_figure_supp_conflict_intensity.R
# Visualise the conflict-INTENSITY sensitivity sweep (baseline vs weak vs
# strong conflict impact on coverage & DPC).
# Reads: output_figgen/figure_supp_conflict_intensity_particle_summary.csv
#        data-processed/SDB_communityDeath_blended.rds
# =============================================================================
source(here::here("analyses", "03_figure_template", "helper_functions_figure_1to4.R"))
OUT_DIR <- here("figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

save_fig <- function(filename_base, plot, width, height) {
  ggsave(file.path(OUT_DIR, paste0(filename_base, ".png")),
         plot, width = width, height = height, dpi = 400, units = "in")
  ggsave(file.path(OUT_DIR, paste0(filename_base, ".pdf")),
         plot, width = width, height = height, units = "in")
}

ARM_NAMES <- c(
  "no_pep",
  "with_conflict_mid", "with_conflict_lo", "with_conflict_hi",
  "cov_conflict_mid",  "cov_conflict_lo",  "cov_conflict_hi",
  "dpc_conflict_mid",  "dpc_conflict_lo",  "dpc_conflict_hi",
  "optimistic_mid",    "optimistic_lo",    "optimistic_hi"
)

# baseline = existing scenario (black), weak = improved/mitigated (blue),
# strong = worsened (red)
COND_LEVELS <- c("baseline", "weak", "strong")
COND_LABELS <- c(baseline = "Baseline", weak = "Improved (weaker conflict impact)",
                 strong   = "Worsened (stronger conflict impact)")
COND_COLORS <- c(baseline = "grey40", weak = "#2166AC", strong = "#B2182B")

# Intensity conditions -- (coverage_max, dpc_max); baseline matches
# 01_analysis_figure3new_conflict_dpc.R, weak/strong match
# 01_analysis_figure3new_conflict_dpc_sensitivity_intensity.R
INTENSITY_PARAMS <- list(
  baseline = list(coverage_max = 80,  dpc_max = 4),
  weak     = list(coverage_max = 100, dpc_max = 1),
  strong   = list(coverage_max = 60,  dpc_max = 7)
)

# =============================================================================
# Left panel: coverage & DPC over time, for each condition
# =============================================================================
sdb <- readRDS(here("data-processed", "SDB_communityDeath_blended.rds"))

rescale_sdb_segment <- function(sdb_ref, day_out, orig_from, orig_to) {
  t        <- (day_out - day_out[1]) / (day_out[length(day_out)] - day_out[1])
  day_orig <- orig_from + t * (orig_to - orig_from)
  approx(sdb_ref$day, sdb_ref$value, xout = day_orig, rule = 2)$y
}

idx_150_325 <- sdb$day >= 150 & sdb$day <= 325
idx_325_400 <- sdb$day >  325 & sdb$day <= 400

sdb_tweaked <- sdb$value
sdb_tweaked[idx_150_325] <- rescale_sdb_segment(sdb, sdb$day[idx_150_325], 150, 200)
sdb_tweaked[idx_325_400] <- rescale_sdb_segment(sdb, sdb$day[idx_325_400], 200, 350)
sdb$value_tweaked <- sdb_tweaked

build_conflict_curves <- function(sdb, coverage_max, dpc_max) {
  sdb$coverage_conflict <- sdb$value_tweaked * coverage_max / max(sdb$value_tweaked)
  sdb$dpc_conflict      <- 1 + dpc_max * (1 - (sdb$value_tweaked / max(sdb$value_tweaked)))
  
  sub      <- sdb[sdb$day < 200, ]
  peak_row <- sub[which.max(sub$coverage_conflict), ]
  peak_day <- peak_row$day
  sdb$dpc_conflict[sdb$day <= peak_day] <- 1
  
  sdb
}

# Cap the DPC secondary axis at a round number that comfortably fits the
# strongest condition's peak DPC (~8 days), while the primary (coverage)
# axis always spans 0-100%.
DPC_AXIS_MAX <- 10
scale_factor <- DPC_AXIS_MAX / 100

curve_df <- do.call(rbind, lapply(COND_LEVELS, function(cond_name) {
  p    <- INTENSITY_PARAMS[[cond_name]]
  sdb_c <- build_conflict_curves(sdb, p$coverage_max, p$dpc_max)
  rbind(
    data.frame(day = sdb_c$day, value = sdb_c$coverage_conflict,
               curve_type = "Coverage", condition = cond_name),
    data.frame(day = sdb_c$day, value = sdb_c$dpc_conflict / scale_factor,
               curve_type = "DPC", condition = cond_name)
  )
})) %>%
  mutate(
    condition_label = factor(COND_LABELS[condition], levels = COND_LABELS[COND_LEVELS]),
    curve_type      = factor(curve_type, levels = c("Coverage", "DPC"))
  )

panel_curves <- ggplot(curve_df, aes(x = day, y = value, color = condition_label, linetype = curve_type)) +
  annotate("rect", xmin = 110, xmax = 300, ymin = -Inf, ymax = Inf,
           fill = "grey85", alpha = 0.6) +
  annotate("text", x = (110 + 300) / 2, y = 97,
           label = "conflict", size = 3.5, color = "grey30") +
  geom_line(linewidth = 1.0) +
  scale_color_manual(values = setNames(COND_COLORS[COND_LEVELS], COND_LABELS[COND_LEVELS]), name = NULL) +
  scale_linetype_manual(values = c(Coverage = "solid", DPC = "dashed"), name = NULL) +
  scale_y_continuous(
    name     = "Coverage (%)",
    limits   = c(0, 100),
    sec.axis = sec_axis(~ . * scale_factor, name = "DPC (days)")
  ) +
  labs(x = "Day") +
  theme_fig() +
  theme(legend.position = "bottom", legend.box = "vertical")

# =============================================================================
# Right panel: HCW deaths averted (%) by arm, dodged/colored by condition
# =============================================================================
particle_df_intensity <- read.csv(
  here("output_figgen", "figure_supp_conflict_intensity_particle_summary.csv"),
  stringsAsFactors = FALSE
)

ARM_PANEL_ORDER <- c(
  "with_conflict_mid", "with_conflict_lo", "with_conflict_hi",
  "dpc_conflict_mid",  "dpc_conflict_lo",  "dpc_conflict_hi"
)

averted_summ <- particle_df_intensity %>%
  filter(arm %in% ARM_PANEL_ORDER) %>%
  group_by(scenario, arm, condition) %>%
  summarise(
    median = median(pct_hcw_deaths_averted, na.rm = TRUE),
    lo95   = quantile(pct_hcw_deaths_averted, 0.025, na.rm = TRUE),
    hi95   = quantile(pct_hcw_deaths_averted, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    condition_label = factor(COND_LABELS[condition], levels = COND_LABELS[COND_LEVELS]),
    arm_label       = factor(arm, levels = ARM_PANEL_ORDER)
  )

panel_bars <- ggplot(averted_summ, aes(x = arm_label, y = median, fill = condition_label)) +
  geom_col(position = position_dodge(0.75), width = 0.7) +
  geom_errorbar(aes(ymin = lo95, ymax = hi95),
                position = position_dodge(0.75), width = 0.2, color = "black", linewidth = 0.4) +
  scale_fill_manual(values = setNames(COND_COLORS[COND_LEVELS], COND_LABELS[COND_LEVELS]), name = NULL) +
  scale_y_continuous(limits = c(0, 100), labels = function(x) paste0(x, "%")) +
  labs(x = NULL, y = "HCW deaths averted (%)") +
  theme_fig() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "bottom")

# =============================================================================
# Combine and save
# =============================================================================
fig_intensity <- (panel_curves | panel_bars) +
  plot_layout(widths = c(1, 1.3)) +
  plot_annotation(tag_levels = "a")

save_fig("figure_supp_conflict_intensity", fig_intensity, 14, 6)
message("figure_supp_conflict_intensity saved.")