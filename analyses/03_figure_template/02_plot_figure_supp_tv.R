# =============================================================================
# 02_plot_figure_supp_tv.R
# Visualise the time-varying-parameter sensitivity sweep.
# Reads: output_figgen/figure_supp_tv_particle_summary.csv
#        data-processed/final_six_scenario_values_original_approach.csv
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

SCENARIO_ID <- "Middle_DRC_ConflictSmoothed_PlusPlus"

# goodgood = strongest good, good = mild good, baseline = unperturbed,
# bad = mild bad, badbad = strongest bad
COND_LEVELS <- c("goodgood", "good", "baseline", "bad", "badbad")
COND_LABELS <- c(
  goodgood = "Good x2 (50%)",
  good     = "Good (25%)",
  baseline = "Baseline",
  bad      = "Bad (25%)",
  badbad   = "Bad x2 (50%)"
)
COND_COLORS <- c(
  goodgood = "#08519C",  # dark blue
  good     = "#6BAED6",  # light blue
  baseline = "black",
  bad      = "#FC9272",  # light red
  badbad   = "#A50026"   # dark red
)

# =============================================================================
# Left panel: sensitivity target curve -- prob_unsafe_funeral_comm over time,
# shown as one representative example of the five TV parameters perturbed
# =============================================================================
get_multiplier <- function(direction, type, tier) {
  increase <- if (tier == 2) 1.5 else 1.25
  decrease <- if (tier == 2) 0.5 else 0.75
  if (direction == "good_high") {
    if (type == "good") increase else decrease
  } else {
    if (type == "good") decrease else increase
  }
}

COND_PARAMS <- list(
  goodgood = list(type = "good", tier = 2),
  good     = list(type = "good", tier = 1),
  baseline = list(type = NA,     tier = NA),
  bad      = list(type = "bad",  tier = 1),
  badbad   = list(type = "bad",  tier = 2)
)

scenario_matrix <- read.csv(here("data-processed", "final_six_scenario_values_original_approach.csv"),
                            stringsAsFactors = FALSE)
scenario_rows <- scenario_matrix[scenario_matrix$scenario == SCENARIO_ID, ]

curve_df <- do.call(rbind, lapply(COND_LEVELS, function(cond_name) {
  p    <- COND_PARAMS[[cond_name]]
  mult <- if (cond_name == "baseline") 1 else get_multiplier("good_low", p$type, p$tier)
  vals <- pmin(scenario_rows$prob_unsafe_funeral_comm * mult, 1)  # capped param
  data.frame(
    relative_day = scenario_rows$relative_day,
    value         = vals,
    condition     = cond_name
  )
})) %>%
  mutate(condition_label = factor(COND_LABELS[condition], levels = COND_LABELS[COND_LEVELS]))

panel_target <- ggplot(curve_df, aes(x = relative_day, y = value, color = condition_label)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = setNames(COND_COLORS[COND_LEVELS], COND_LABELS[COND_LEVELS]), name = NULL) +
  scale_y_continuous(limits = c(0, NA), labels = scales::percent) +
  labs(x = "Relative day", y = "P(unsafe funeral, community)") +
  theme_fig() +
  theme(legend.position = "bottom")

# =============================================================================
# Load simulation-based particle summary
# =============================================================================
particle_df_tv <- read.csv(here("output_figgen", "figure_supp_tv_particle_summary.csv"),
                           stringsAsFactors = FALSE)

ARM_PANEL_ORDER <- c("with_conflict_mid", "with_conflict_lo", "with_conflict_hi")

summ_tv <- particle_df_tv %>%
  filter(arm %in% ARM_PANEL_ORDER) %>%
  group_by(scenario, arm, condition) %>%
  summarise(
    med_deaths    = median(n_hcw_deaths, na.rm = TRUE),
    lo_deaths     = quantile(n_hcw_deaths, 0.025, na.rm = TRUE),
    hi_deaths     = quantile(n_hcw_deaths, 0.975, na.rm = TRUE),
    med_averted   = median(pct_hcw_deaths_averted, na.rm = TRUE),
    lo_averted    = quantile(pct_hcw_deaths_averted, 0.025, na.rm = TRUE),
    hi_averted    = quantile(pct_hcw_deaths_averted, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    condition_label = factor(COND_LABELS[condition], levels = COND_LABELS[COND_LEVELS]),
    arm_label       = factor(arm, levels = ARM_PANEL_ORDER)
  )

# =============================================================================
# Middle panel: raw HCW deaths (absolute counts)
# =============================================================================
panel_deaths <- ggplot(summ_tv, aes(x = arm_label, y = med_deaths, fill = condition_label)) +
  geom_col(position = position_dodge(0.75), width = 0.7) +
  geom_errorbar(aes(ymin = lo_deaths, ymax = hi_deaths),
                position = position_dodge(0.75), width = 0.2, color = "black", linewidth = 0.4) +
  scale_fill_manual(values = setNames(COND_COLORS[COND_LEVELS], COND_LABELS[COND_LEVELS]), name = NULL) +
  labs(x = NULL, y = "HCW deaths") +
  theme_fig() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "bottom")

# =============================================================================
# Right panel: HCW deaths averted (%)
# =============================================================================
panel_averted <- ggplot(summ_tv, aes(x = arm_label, y = med_averted, fill = condition_label)) +
  geom_col(position = position_dodge(0.75), width = 0.7) +
  geom_errorbar(aes(ymin = lo_averted, ymax = hi_averted),
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
fig_tv <- (panel_target | panel_deaths | panel_averted) +
  plot_layout(widths = c(1, 1, 1), guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom")

save_fig("figure_supp_tv", fig_tv, 18, 6)
message("figure_supp_tv saved.")