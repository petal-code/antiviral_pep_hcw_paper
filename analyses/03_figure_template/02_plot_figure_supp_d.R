# =============================================================================
# 02_plot_figure_supp_d.R
#
# Two standalone panels, neither derived from simulation output:
#   figure_supp_DPC      : antiviral efficacy vs. days post-exposure (DPC)
#   figure_supp_conflict : antiviral coverage & DPC trajectory over time,
#                          under the conflict scenario
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

# =============================================================================
# Load data
# =============================================================================
curve_d50_dat <- readRDS(here("data-processed", "DPC_fixed_efficacy_varied_d50.rds"))

# =============================================================================
# SDB: load and apply conflict tweak (150-400 day window)
# Shape-preserving x-axis rescale: curve form maintained, trough shifted
# from day 200 to day 325.
#   150-200 (original descent) stretched onto 150-325
#   200-350 (original recovery) squeezed onto 325-400
# Days outside 150-400 are untouched.
# =============================================================================
sdb <- readRDS(here("data-processed", "SDB_communityDeath_blended.rds"))

rescale_sdb_segment <- function(day_out, orig_from, orig_to) {
  t        <- (day_out - day_out[1]) / (day_out[length(day_out)] - day_out[1])
  day_orig <- orig_from + t * (orig_to - orig_from)
  approx(sdb$day, sdb$value, xout = day_orig, rule = 2)$y
}

sdb$value_tweaked <- sdb$value
idx_150_325 <- sdb$day >= 150 & sdb$day <= 325
idx_325_400 <- sdb$day >  325 & sdb$day <= 400
sdb$value_tweaked[idx_150_325] <- rescale_sdb_segment(sdb$day[idx_150_325], 150, 200)
sdb$value_tweaked[idx_325_400] <- rescale_sdb_segment(sdb$day[idx_325_400], 200, 350)

sdb$coverage_conflict <- sdb$value_tweaked * 80 / max(sdb$value_tweaked)
sdb$dpc_conflict      <- 1 + 6 * (1 - (sdb$value_tweaked / max(sdb$value_tweaked)))

sub      <- sdb[sdb$day < 200, ]
peak_row <- sub[which.max(sub$coverage_conflict), ]
peak_day <- peak_row$day

sdb$dpc_conflict[sdb$day <= peak_day] <- 1

# Flat version: coverage held constant at its peak once reached (used as
# the "Delayed dosing" scenario's coverage curve, unaffected by conflict).
sdb$coverage_flat <- sdb$coverage_conflict
sdb$coverage_flat[sdb$day > peak_day] <- peak_row$coverage_conflict

# Shared constants
scale_factor <- 7.5 / 80  # cap DPC secondary axis at 5.5 when coverage axis tops out at 80
cov_color    <- "#E08214"

COV_SCEN_COLORS <- c(
  "Ideal (100% coverage, 0 delay)" = "#1a9641",
  "Delayed dosing"                 = "#f58231",
  "Delayed coverage + dosing"      = "#d7191c"
)

# =============================================================================
# figure_supp_DPC: efficacy vs DPC
# Both the hi and lo efficacy bounds are drawn as black dashed lines (instead
# of green/red) to avoid the colors being mistaken for a meaningful encoding.
# The region between them is shaded light grey, and each dashed line is
# labeled directly (Optimistic / Pessimistic) near its right-hand end.
# =============================================================================
panel_dpc <- ggplot(curve_d50_dat, aes(x = dpc)) +
  geom_ribbon(aes(ymin = eighty_efficacy_lo, ymax = eighty_efficacy_hi),
              fill = "grey80", alpha = 0.4) +
  geom_line(aes(y = eighty_efficacy_hi), color = "black", linetype = "dashed", linewidth = 1) +
  geom_line(aes(y = efficacy),            color = "black", linewidth = 1.2) +
  geom_line(aes(y = eighty_efficacy_lo), color = "black", linetype = "dashed", linewidth = 1) +
  # Fixed label positions matching the requested locations: "Optimistic" up and
  # to the right along the upper dashed curve, "Pessimistic" lower-left along
  # the lower dashed curve's early decline. Adjust x/y here if the underlying
  # curve shape changes and the labels drift off the lines.
  annotate("text", x = 6.2,   y = 0.75, label = "Optimistic",  hjust = 0, vjust = 0, size = 3.2, color = "black") +
  annotate("text", x = 0.5, y = 0.35, label = "Pessimistic", hjust = 0, vjust = 0, size = 3.2, color = "black") +
  scale_y_continuous(limits = c(0, NA), labels = scales::percent) +
  labs(x = "Days post-exposure (DPC)", y = "Efficacy") +
  theme_fig()

save_fig("figure_supp_DPC", panel_dpc, 5, 4)

# =============================================================================
# figure_supp_conflict: coverage & DPC over time
# Conflict period shaded from day 110 to 300 (fixed)
#
# Lines are colored to match the "Delayed dosing" / "Delayed coverage +
# dosing" scenario colors used elsewhere:
#   - coverage_flat (coverage stays at its peak, unaffected by conflict)
#     is the coverage curve used by the "Delayed dosing" scenario -> orange
#   - coverage_conflict (coverage declines during conflict) is the
#     coverage curve used by the "Delayed coverage + dosing" scenario -> red
#   - dpc_conflict (DPC curve under conflict) is shared identically by
#     both scenarios, so it's kept neutral black/dashed rather than tied
#     to a single scenario color
# DPC secondary axis capped at 5.5.
# =============================================================================
panel_conflict <- ggplot(sdb, aes(x = day)) +
  annotate("rect", xmin = 110, xmax = 300, ymin = -Inf, ymax = Inf,
           fill = "grey85", alpha = 0.6) +
  annotate("text", x = (110 + 300) / 2, y = 78,
           label = "conflict", size = 3.5, color = "grey30") +
  geom_line(aes(y = coverage_flat),
            color = COV_SCEN_COLORS[["Delayed dosing"]], linetype = "solid", linewidth = 1.1) +
  geom_line(aes(y = coverage_conflict),
            color = COV_SCEN_COLORS[["Delayed coverage + dosing"]], linetype = "solid", linewidth = 1.1) +
  geom_line(aes(y = dpc_conflict / scale_factor),
            color = "black", linetype = "dashed", linewidth = 0.9) +
  scale_y_continuous(
    name     = "Coverage (%)",
    limits   = c(0, 80),
    sec.axis = sec_axis(~ . * scale_factor, name = "DPC (days)")
  ) +
  labs(x = "Day") +
  theme_fig() +
  theme(
    axis.title.y       = element_text(color = cov_color),
    axis.text.y        = element_text(color = cov_color),
    axis.title.y.right = element_text(color = "black"),
    axis.text.y.right  = element_text(color = "black")
  )

save_fig("figure_supp_conflict", panel_conflict, 6, 4)

message("figure_supp_DPC and figure_supp_conflict saved.")

# =============================================================================
# figure_3_panel_c_by_eff_incident
#
# Unlike the two panels above, this one IS derived from simulation output:
# weekly incident HCW deaths by scenario, faceted by efficacy arm (one panel
# per efficacy level, scenarios as lines within each panel). Uses the
# matched-seed no_pep_{eff} counterfactual for each efficacy level.
# Requires output_figgen/figure_3_weekly_ts.csv (see 02_extract_figure3.R).
# =============================================================================
ts_3 <- read.csv(here("output_figgen", "figure_3_weekly_ts.csv"))

EFF_ARM_LABELS <- c(hi = "Optimistic", mid = "Central", lo = "Pessimistic")
EFF_ARM_ORDER  <- c("hi", "mid", "lo")

EFF_SCENARIO_COLORS <- c(
  "No PEP"                          = "black",
  "Ideal (100% coverage, 0 delay)"  = "#1a9641",
  "Delayed dosing"                  = "#f58231",
  "Delayed coverage + dosing"       = "#d7191c"
)

make_panel_c_by_eff <- function(metric_name, y_label, eff, eff_title, show_legend = FALSE) {
  arm_map <- setNames(
    c("No PEP",
      "Ideal (100% coverage, 0 delay)",
      "Delayed dosing",
      "Delayed coverage + dosing"),
    c(paste0("no_pep_",       eff),
      paste0("optimistic_",    eff),
      paste0("dpc_conflict_",  eff),
      paste0("with_conflict_", eff))
  )
  df <- ts_3 %>%
    filter(scenario == "DRC", arm %in% names(arm_map), metric == metric_name) %>%
    mutate(
      day       = week * 7,
      arm_label = factor(arm_map[arm], levels = arm_map)
    )
  p <- ggplot(df, aes(x = day, y = q50, color = arm_label, fill = arm_label)) +
    geom_ribbon(aes(ymin = pmax(q25, 0), ymax = q75), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = EFF_SCENARIO_COLORS, name = NULL) +
    scale_fill_manual(values = EFF_SCENARIO_COLORS, name = NULL) +
    scale_x_continuous(limits = c(0, 600), expand = c(0, 0)) +
    labs(x = "Day", y = y_label,
         title = paste0("Delayed coverage + dosing analysis (", eff_title, " efficacy)")) +
    theme_fig()
  if (show_legend) {
    p <- p + theme(
      legend.position      = c(0.02, 0.98),
      legend.justification = c(0, 1),
      legend.direction     = "vertical",
      legend.background    = element_blank(),
      legend.key           = element_blank(),
      legend.text          = element_text(size = 8)
    ) + guides(color = guide_legend(ncol = 1),
               fill  = guide_legend(ncol = 1))
  } else {
    p <- p + theme(legend.position = "none")
  }
  p
}

panel_c_hi_inc  <- make_panel_c_by_eff("hcw_deaths_incidence", "Mean weekly incident HCW deaths", "hi",  "optimistic",  show_legend = TRUE)
panel_c_mid_inc <- make_panel_c_by_eff("hcw_deaths_incidence", "Mean weekly incident HCW deaths", "mid", "central",     show_legend = FALSE)
panel_c_lo_inc  <- make_panel_c_by_eff("hcw_deaths_incidence", "Mean weekly incident HCW deaths", "lo",  "pessimistic", show_legend = FALSE)

save_fig("figure_supp_inci_by_effi",
         (panel_c_hi_inc / panel_c_mid_inc / panel_c_lo_inc) +
           plot_annotation(tag_levels = "a"),
         6, 9)
