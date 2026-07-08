# =============================================================================
# 02_plot_figure3.R
#
# Combined figure: incident HCW deaths time series (left) + waterfall of
# cumulative HCW deaths by scenario with 95% CrI whiskers (middle/right) +
# HCW deaths averted (%) by scenario x efficacy arm (bottom).
# =============================================================================
source(here::here("analyses", "03_figure_template", "helper_functions_figure_1to4.R"))

OUT_DIR <- here("figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

save_fig <- function(filename_base, plot, width, height) {
  ggsave(file.path(OUT_DIR, paste0(filename_base, ".png")),
         plot, width = width, height = height, dpi = 400, units = "in")
  ggsave(file.path(OUT_DIR, paste0(filename_base, ".tiff")),
         plot, width = width, height = height, dpi = 400, units = "in")
  ggsave(file.path(OUT_DIR, paste0(filename_base, ".pdf")),
         plot, width = width, height = height, units = "in")
}

# =============================================================================
# Load data
# =============================================================================
ts_3       <- read.csv(here("output_figgen", "figure_3_weekly_ts.csv"))
particle_3 <- read.csv(here("output_figgen", "figure_3_particle_summary.csv"))

# =============================================================================
# Shared label/color constants
# =============================================================================
EFF_ARM_LABELS <- c(hi = "Optimistic", mid = "Central", lo = "Pessimistic")
EFF_ARM_ORDER  <- c("hi", "mid", "lo")

# Arms for the incidence time-series panel (mid efficacy, scenario comparison).
# "No PEP" points at no_pep_mid, the matched-seed counterfactual (tdf +
# prevented from the with_conflict_mid runs themselves).
ARM_LABELS_TS <- c(
  no_pep_mid        = "No PEP",
  optimistic_mid    = "Ideal (100% coverage, 0 delay)",
  dpc_conflict_mid  = "Delayed dosing",
  with_conflict_mid = "Delayed coverage + dosing"
)
ARM_COLORS_TS <- c(
  "No PEP"                          = "black",
  "Ideal (100% coverage, 0 delay)"  = "#1a9641",
  "Delayed dosing"                  = "#f58231",
  "Delayed coverage + dosing"       = "#d7191c"
)

# =============================================================================
# Left panel: incident HCW deaths time series (mid efficacy, scenario comparison)
# =============================================================================
make_ts_panel <- function(metric_name, y_label) {
  df <- ts_3 %>%
    filter(scenario == "DRC", arm %in% names(ARM_LABELS_TS), metric == metric_name) %>%
    mutate(
      day       = week * 7,
      arm_label = factor(ARM_LABELS_TS[arm], levels = ARM_LABELS_TS)
    )
  ggplot(df, aes(x = day, y = q50, color = arm_label, fill = arm_label)) +
    geom_ribbon(aes(ymin = pmax(q25, 0), ymax = q75), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = ARM_COLORS_TS, name = NULL) +
    scale_fill_manual(values = ARM_COLORS_TS, name = NULL) +
    scale_x_continuous(limits = c(0, 600), expand = c(0, 0)) +
    labs(x = "Days since outbreak start", y = y_label) +
    theme_fig() +
    theme(
      legend.position       = c(0.02, 0.98),
      legend.justification  = c(0, 1),
      legend.direction      = "vertical",
      legend.background     = element_blank(),
      legend.key            = element_blank(),
      legend.text           = element_text(size = 8),
      plot.margin           = margin(5, 10, 2, 5)
    ) +
    guides(color = guide_legend(ncol = 1), fill = guide_legend(ncol = 1))
}

panel_ts_incident <- make_ts_panel("hcw_deaths_incidence", "Mean weekly incident HCW deaths")

# =============================================================================
# Middle panel: waterfall of cumulative HCW deaths (Ideal -> Delayed dosing ->
# Delayed coverage + dosing -> No antiviral), with 95% CrI whiskers on each bar
# =============================================================================
y_max <- 140
WF_BW <- 0.28  # half-width of each bar

final_week <- max(ts_3$week[ts_3$scenario == "DRC" & ts_3$metric == "hcw_deaths"], na.rm = TRUE)

wf_vals <- ts_3 %>%
  filter(scenario == "DRC", metric == "hcw_deaths", week == final_week,
         arm %in% c("optimistic_mid", "dpc_conflict_mid", "with_conflict_mid", "no_pep_mid")) %>%
  mutate(stage = c(
    optimistic_mid    = "Ideal",
    dpc_conflict_mid  = "Delayed dosing",
    with_conflict_mid = "Delayed coverage + dosing",
    no_pep_mid        = "No antiviral"
  )[arm])

ideal_val  <- wf_vals$q50[wf_vals$stage == "Ideal"]
dpc_val    <- wf_vals$q50[wf_vals$stage == "Delayed dosing"]
both_val   <- wf_vals$q50[wf_vals$stage == "Delayed coverage + dosing"]
no_pep_val <- wf_vals$q50[wf_vals$stage == "No antiviral"]

wf_bars <- data.frame(
  x          = 1:4,
  ymin       = c(0,         ideal_val, dpc_val,  both_val),
  ymax       = c(ideal_val, dpc_val,   both_val, no_pep_val),
  fill_group = c("Ideal", "Delayed dosing", "Delayed coverage + dosing", "No antiviral")
)

WF_COLORS <- c(
  "Ideal"                     = "#1a9641",
  "Delayed dosing"            = "#f58231",
  "Delayed coverage + dosing" = "#d7191c",
  "No antiviral"              = "grey60"
)

wf_segs <- data.frame(
  x    = wf_bars$x[-nrow(wf_bars)] - WF_BW,
  xend = wf_bars$x[-1]              + WF_BW,
  y    = wf_bars$ymax[-nrow(wf_bars)]
)

wf_xlabels <- data.frame(
  x     = 1:4,
  label = c("Ideal", "Delayed\ndosing", "Delayed\ncoverage\n+ dosing", "No\nantiviral")
)

# 95% CrI whiskers at the top of each bar
wf_whiskers <- wf_bars %>%
  left_join(wf_vals %>% select(fill_group = stage, q50, q25, q75), by = "fill_group") %>%
  mutate(
    w_lo = ymax - (q50 - q25),
    w_hi = pmin(ymax + (q75 - q50), y_max)
  )

panel_waterfall <- ggplot() +
  geom_rect(data = wf_bars,
            aes(xmin = x - WF_BW, xmax = x + WF_BW, ymin = ymin, ymax = ymax, fill = fill_group),
            color = "grey50", linewidth = 0.25) +
  geom_segment(data = wf_segs,
               aes(x = x, xend = xend, y = y, yend = y),
               linetype = "dotted", color = "black", linewidth = 0.9) +
  geom_errorbar(data = wf_whiskers,
                aes(x = x, ymin = w_lo, ymax = w_hi),
                width = WF_BW * 0.8, color = "black", linewidth = 0.6) +
  scale_fill_manual(values = WF_COLORS, guide = "none") +
  scale_x_continuous(breaks = wf_xlabels$x, labels = wf_xlabels$label) +
  scale_y_continuous(limits = c(0, y_max), expand = c(0, 0)) +
  labs(x = NULL, y = "Cumulative HCW deaths") +
  theme_fig() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.ticks.x = element_blank())

# =============================================================================
# Bottom panel: HCW deaths averted (%) by scenario x efficacy arm
# =============================================================================
SCEN_ORDER  <- c("Ideal", "Delayed dosing", "Delayed coverage + dosing")
SCEN_COLORS <- c(
  "Ideal"                     = "#1a9641",
  "Delayed dosing"            = "#f58231",
  "Delayed coverage + dosing" = "#d7191c"
)

ARM_TO_SCEN <- c(
  optimistic_hi     = "Ideal",                     optimistic_mid    = "Ideal",                     optimistic_lo    = "Ideal",
  dpc_conflict_hi   = "Delayed dosing",             dpc_conflict_mid  = "Delayed dosing",             dpc_conflict_lo  = "Delayed dosing",
  with_conflict_hi  = "Delayed coverage + dosing",  with_conflict_mid = "Delayed coverage + dosing",  with_conflict_lo = "Delayed coverage + dosing"
)

averted_summary <- particle_3 %>%
  filter(arm %in% names(ARM_TO_SCEN)) %>%
  mutate(
    scenario_label = factor(ARM_TO_SCEN[arm], levels = SCEN_ORDER),
    eff_arm        = sub("^(optimistic|dpc_conflict|with_conflict)_", "", arm),
    eff_arm_label  = factor(EFF_ARM_LABELS[eff_arm], levels = EFF_ARM_LABELS[EFF_ARM_ORDER])
  ) %>%
  group_by(eff_arm_label, scenario_label) %>%
  summarise(
    median = median(pct_hcw_deaths_averted, na.rm = TRUE),
    lo     = quantile(pct_hcw_deaths_averted, 0.025, na.rm = TRUE),
    hi     = quantile(pct_hcw_deaths_averted, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

panel_averted <- ggplot(averted_summary,
                        aes(x = scenario_label, y = median,
                            fill = scenario_label, ymin = lo, ymax = hi)) +
  geom_col(position = position_dodge(0.7), width = 0.6, color = "grey70", linewidth = 0.3) +
  geom_errorbar(position = position_dodge(0.7), width = 0.2, linewidth = 0.4) +
  scale_fill_manual(values = SCEN_COLORS, name = NULL) +
  scale_y_continuous(limits = c(0, 100), labels = function(x) paste0(x, "%"), expand = c(0, 0)) +
  facet_wrap(~ eff_arm_label, nrow = 1) +
  labs(x = NULL, y = "HCW deaths averted (%)") +
  theme_fig() +
  theme(
    legend.position  = "none",
    axis.text.x      = element_text(angle = 25, hjust = 1),
    strip.background = element_blank(),
    strip.text       = element_text(size = 10, face = "plain")
  )

# =============================================================================
# Assemble and save
# =============================================================================
top_row <- wrap_plots(panel_ts_incident, panel_waterfall, widths = c(2, 1))
fig3    <- wrap_plots(top_row, panel_averted, ncol = 1, heights = c(2, 1)) +
  plot_annotation(tag_levels = "a")

save_fig("figure_3", fig3, 10, 7)
message("Figure 3 plotting complete.")