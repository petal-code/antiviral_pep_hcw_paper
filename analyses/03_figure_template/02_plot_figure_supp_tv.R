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
  goodgood = "Good x2 (20%)",
  good     = "Good (10%)",
  baseline = "Baseline",
  bad      = "Bad (10%)",
  badbad   = "Bad x2 (20%)"
)
COND_COLORS <- c(
  goodgood = "#08519C",  # dark blue
  good     = "#6BAED6",  # light blue
  baseline = "grey40",
  bad      = "#FC9272",  # light red
  badbad   = "#A50026"   # dark red
)

get_multiplier <- function(direction, type, tier) {
  increase <- if (tier == 2) 1.2 else 1.1
  decrease <- if (tier == 2) 0.8 else 0.9
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

# =============================================================================
# Panel a: all five TV parameters, each as its own sub-panel (2x3 grid,
# last cell blank) showing the five condition curves
# =============================================================================
TV_PARAMS <- c("prob_hosp", "delay_hosp", "prob_unsafe_funeral_comm",
               "prob_unsafe_funeral_hosp", "prop_etu")
TV_LABELS <- c(
  prob_hosp                = "P(hospitalised)",
  delay_hosp                = "Hospitalisation delay",
  prob_unsafe_funeral_comm = "P(unsafe funeral, community)",
  prob_unsafe_funeral_hosp = "P(unsafe funeral, hospital)",
  prop_etu                  = "ETU proportion"
)
DIRECTION_MAP <- c(
  prob_hosp                = "good_high",
  delay_hosp                = "good_low",
  prob_unsafe_funeral_comm = "good_low",
  prob_unsafe_funeral_hosp = "good_low",
  prop_etu                  = "good_high"
)
CAPPED_PARAMS <- c("prob_hosp", "prob_unsafe_funeral_comm", "prob_unsafe_funeral_hosp", "prop_etu")

scenario_matrix <- read.csv(here("data-processed", "final_six_scenario_values_original_approach.csv"),
                            stringsAsFactors = FALSE)
scenario_rows <- scenario_matrix[scenario_matrix$scenario == SCENARIO_ID, ]

tv_curve_df <- do.call(rbind, lapply(TV_PARAMS, function(param) {
  do.call(rbind, lapply(COND_LEVELS, function(cond_name) {
    p    <- COND_PARAMS[[cond_name]]
    mult <- if (cond_name == "baseline") 1 else get_multiplier(DIRECTION_MAP[[param]], p$type, p$tier)
    vals <- scenario_rows[[param]] * mult
    if (param %in% CAPPED_PARAMS) vals <- pmin(vals, 1)
    data.frame(
      relative_day = scenario_rows$relative_day,
      value         = vals,
      condition     = cond_name,
      parameter     = param
    )
  }))
})) %>%
  mutate(
    condition_label = factor(COND_LABELS[condition], levels = COND_LABELS[COND_LEVELS]),
    parameter_label = factor(TV_LABELS[parameter],   levels = TV_LABELS[TV_PARAMS])
  )

panel_target <- ggplot(tv_curve_df, aes(x = relative_day, y = value, color = condition_label)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = setNames(COND_COLORS[COND_LEVELS], COND_LABELS[COND_LEVELS]), name = NULL) +
  facet_wrap(~ parameter_label, ncol = 5, scales = "free_y") +
  labs(x = "Relative day", y = "Parameter value") +
  theme_fig() +
  theme(legend.position = "none", strip.text = element_text(size = 9))

# =============================================================================
# Load simulation-based particle summary
# =============================================================================
particle_df_tv <- read.csv(here("output_figgen", "figure_supp_tv_particle_summary.csv"),
                           stringsAsFactors = FALSE)

# All arms here are "with_conflict" (no dpc_conflict counterpart in this
# sweep), so ticks just need the efficacy level, no second scenario line.
ARM_PANEL_ORDER <- c("with_conflict_lo", "with_conflict_mid", "with_conflict_hi")
EFF_TICK_LABELS <- c(lo = "Pessimistic", mid = "Central", hi = "Optimistic")
x_tick_labels   <- unname(EFF_TICK_LABELS[sub("^with_conflict_", "", ARM_PANEL_ORDER)])

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
# Panel b: raw HCW deaths (absolute counts)
# =============================================================================
panel_deaths <- ggplot(summ_tv, aes(x = arm_label, y = med_deaths, fill = condition_label)) +
  geom_col(position = position_dodge(0.75), width = 0.7) +
  geom_errorbar(aes(ymin = lo_deaths, ymax = hi_deaths),
                position = position_dodge(0.75), width = 0.2, color = "black", linewidth = 0.4) +
  scale_fill_manual(values = setNames(COND_COLORS[COND_LEVELS], COND_LABELS[COND_LEVELS]), name = NULL) +
  scale_x_discrete(labels = x_tick_labels) +
  labs(x = NULL, y = "HCW deaths") +
  theme_fig() +
  theme(legend.position = "none")

# =============================================================================
# Panel c: HCW deaths averted (%) -- legend lives here, inside the panel, top
# =============================================================================
panel_averted <- ggplot(summ_tv, aes(x = arm_label, y = med_averted, fill = condition_label)) +
  geom_col(position = position_dodge(0.75), width = 0.7) +
  geom_errorbar(aes(ymin = lo_averted, ymax = hi_averted),
                position = position_dodge(0.75), width = 0.2, color = "black", linewidth = 0.4) +
  scale_fill_manual(values = setNames(COND_COLORS[COND_LEVELS], COND_LABELS[COND_LEVELS]), name = NULL) +
  scale_x_discrete(labels = x_tick_labels) +
  scale_y_continuous(limits = c(0, 100), labels = function(x) paste0(x, "%")) +
  labs(x = NULL, y = "HCW deaths averted (%)") +
  theme_fig() +
  theme(
    legend.position       = c(0.5, 0.97),
    legend.justification  = c(0.5, 1),
    legend.direction      = "horizontal",
    legend.background     = element_blank(),
    legend.key            = element_blank(),
    legend.text           = element_text(size = 7)
  ) +
  guides(fill = guide_legend(nrow = 2))

# =============================================================================
# Combine and save
# 3x10 grid: top row (1x10) = TV parameter facets, spanning full width;
# bottom two rows (2x10) split evenly left/right between the deaths panel
# and the averted-% panel, each spanning both bottom rows.
# =============================================================================
tv_design <- "
AAAAAAAAAA
BBBBBCCCCC
BBBBBCCCCC
"

fig_tv <- panel_target + panel_deaths + panel_averted +
  plot_layout(design = tv_design) +
  plot_annotation(tag_levels = "a")

save_fig("figure_supp_tv", fig_tv, 14, 8)
message("figure_supp_tv saved.")