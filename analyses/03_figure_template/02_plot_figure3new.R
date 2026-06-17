# =============================================================================
# 02_plot_figure3new.R
# Visualise figure3_new: efficacy curve, coverage/DPC scenarios, weekly HCW
# deaths under each scenario, and decomposition of lost impact by efficacy arm.
#
# Panel a: efficacy vs DPC (data-driven, low/median/high curves)
# Panel b: coverage & DPC over time (with conflict / without conflict)
# Panel c: weekly HCW deaths (incident or cumulative; two separate files)
# Panel d: decomposition of HCW deaths averted % by efficacy arm
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

ts_3new       <- read.csv(here("output_figgen", "figure_3new_weekly_ts.csv"))
particle_3new <- read.csv(here("output_figgen", "figure_3new_particle_summary.csv"))

# =============================================================================
# Panel a: efficacy vs DPC, low/median/high curves (data-driven)
# =============================================================================
curve_d50_dat <- readRDS(here("data-processed", "DPC_fixed_efficacy_varied_d50.rds"))

panel_a <- ggplot(curve_d50_dat, aes(x = dpc)) +
  geom_line(aes(y = eighty_efficacy_hi), color = "#1a9641", linetype = "dashed", linewidth = 1) +
  geom_line(aes(y = efficacy),            color = "black",   linewidth = 1.2) +
  geom_line(aes(y = eighty_efficacy_lo), color = "#d7191c", linetype = "dotted", linewidth = 1) +
  scale_y_continuous(limits = c(0, NA), labels = scales::percent) +
  labs(x = "Days post-exposure (DPC)", y = "Efficacy",
       title = "Efficacy vs days post-exposure (DPC)") +
  theme_fig()

# =============================================================================
# Panel b: coverage & DPC over time (with conflict / without conflict)
# =============================================================================
sdb <- readRDS(here("data-processed", "SDB_communityDeath_blended.rds"))

sdb$coverage_conflict <- sdb$value * 80 / max(sdb$value)
sdb$dpc_conflict       <- 1 + 4 * (1 - sdb$value / max(sdb$value))

sub <- sdb[sdb$day < 200, ]
peak_row <- sub[which.max(sub$coverage_conflict), ]
peak_day <- peak_row$day

sdb$dpc_conflict[sdb$day <= peak_day] <- 1

sdb$coverage_noconflict <- sdb$coverage_conflict
sdb$coverage_noconflict[sdb$day > peak_day] <- peak_row$coverage_conflict
sdb$dpc_noconflict <- 1

scale_factor <- 5 / 80
cov_color <- "#E08214"
dpc_color <- "black"

conflict_start <- peak_day
conflict_end   <- sdb$day[which.min(sdb$dpc_conflict[sdb$day < 365])]

panel_b <- ggplot(sdb, aes(x = day)) +
  annotate("rect", xmin = conflict_start, xmax = conflict_end, ymin = -Inf, ymax = Inf,
           fill = "grey85", alpha = 0.6) +
  annotate("text", x = (conflict_start + conflict_end) / 2, y = 78,
           label = "conflict", size = 3.5, color = "grey30") +
  geom_line(aes(y = coverage_conflict), color = cov_color, linetype = "dashed", linewidth = 1.1) +
  geom_line(aes(y = coverage_noconflict), color = cov_color, linetype = "solid", linewidth = 1.1) +
  geom_line(aes(y = dpc_conflict / scale_factor), color = dpc_color, linetype = "dashed", linewidth = 0.9) +
  geom_line(aes(y = dpc_noconflict / scale_factor), color = dpc_color, linetype = "solid", linewidth = 0.9) +
  scale_y_continuous(
    name   = "Coverage",
    sec.axis = sec_axis(~ . * scale_factor, name = "DPC (days)")
  ) +
  labs(x = "Day", title = "Coverage & DPC over time") +
  theme_fig() +
  theme(
    axis.title.y       = element_text(color = cov_color),
    axis.text.y        = element_text(color = cov_color),
    axis.title.y.right  = element_text(color = dpc_color),
    axis.text.y.right   = element_text(color = dpc_color)
  )

# =============================================================================
# Panel c: weekly HCW deaths (mid efficacy arm), no PEP / with conflict /
# without conflict / optimistic. Two versions: incident and cumulative.
# =============================================================================
ARM_LABELS_C <- c(
  no_pep               = "No PEP",
  optimistic_mid       = "Ideal (100% coverage, 0 delay)",
  without_conflict_mid = "Without conflict",
  with_conflict_mid    = "With conflict"
)
ARM_COLORS_C <- c(
  no_pep               = "black",
  optimistic_mid       = "#1a9641",
  without_conflict_mid = "#E08214",
  with_conflict_mid    = "#d7191c"
)

make_panel_c <- function(sc, metric_name, y_label) {
  df <- ts_3new %>%
    filter(scenario == sc, arm %in% names(ARM_LABELS_C), metric == metric_name) %>%
    mutate(arm_label = factor(ARM_LABELS_C[arm], levels = ARM_LABELS_C))
  
  ggplot(df, aes(x = week, y = q50, color = arm_label)) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = setNames(ARM_COLORS_C, ARM_LABELS_C), name = NULL) +
    labs(x = "Week", y = y_label,
         title = sprintf("%s (mid efficacy)", SCENARIO_LABELS[sc])) +
    theme_fig()
}

panel_c_incident   <- make_panel_c("DRC", "hcw_deaths_incidence", "Weekly incident HCW deaths")
panel_c_cumulative <- make_panel_c("DRC", "hcw_deaths",            "Cumulative HCW deaths")

# =============================================================================
# Panel d: decomposition of HCW deaths averted % by efficacy arm
# =============================================================================
EFF_ARM_LABELS <- c(lo = "Pessimistic", mid = "Central", hi = "Optimistic")
EFF_ARM_ORDER  <- c("lo", "mid", "hi")

decomp_df <- particle_3new %>%
  filter(arm %in% c("optimistic_mid", "optimistic_lo", "optimistic_hi",
                    "without_conflict_mid", "without_conflict_lo", "without_conflict_hi",
                    "with_conflict_mid", "with_conflict_lo", "with_conflict_hi")) %>%
  mutate(
    eff_arm  = sub("^(optimistic|without_conflict|with_conflict)_", "", arm),
    cov_scen = sub("_(mid|lo|hi)$", "", arm)
  ) %>%
  group_by(scenario, eff_arm, cov_scen) %>%
  summarise(median_pct = median(pct_hcw_deaths_averted, na.rm = TRUE), .groups = "drop")

decomp_wide <- decomp_df %>%
  pivot_wider(names_from = cov_scen, values_from = median_pct) %>%
  mutate(
    realised = with_conflict,
    dpc_loss = pmax(without_conflict - with_conflict, 0),
    cov_loss = pmax(optimistic - without_conflict, 0)
  ) %>%
  select(scenario, eff_arm, realised, dpc_loss, cov_loss) %>%
  pivot_longer(cols = c(realised, dpc_loss, cov_loss),
               names_to = "component", values_to = "value") %>%
  mutate(
    component = factor(component, levels = c("realised", "dpc_loss", "cov_loss"),
                       labels = c("Realised", "Conflict (dpc)", "Coverage (cov)")),
    eff_arm_label = factor(EFF_ARM_LABELS[eff_arm],
                           levels = EFF_ARM_LABELS[EFF_ARM_ORDER])
  )

DECOMP_COLORS <- c("Realised" = "#1a9641", "Conflict (dpc)" = "#d7191c", "Coverage (cov)" = "#fdae61")

make_panel_d <- function(sc) {
  df <- filter(decomp_wide, scenario == sc)
  ggplot(df, aes(x = eff_arm_label, y = value, fill = component)) +
    geom_col(position = "stack", width = 0.6, color = "black", linewidth = 0.3) +
    scale_fill_manual(values = DECOMP_COLORS, name = NULL) +
    labs(x = "Antiviral efficacy", y = "HCW deaths averted (%)",
         title = sprintf("Decomposition by DPC (%s)", SCENARIO_LABELS[sc])) +
    theme_fig()
}

panel_d <- make_panel_d("DRC")

# =============================================================================
# Save individual panels
# =============================================================================
save_fig("figure_3new_panel_a_efficacy", panel_a, 5, 4)
save_fig("figure_3new_panel_b_coverage_dpc", panel_b, 6, 4)
save_fig("figure_3new_panel_c_incident", panel_c_incident, 6, 4)
save_fig("figure_3new_panel_c_cumulative", panel_c_cumulative, 6, 4)
save_fig("figure_3new_panel_d_decomposition", panel_d, 5, 4)

# Combined layout (a top-left, b top-right, c bottom-left, d bottom-right)
fig_combined <- (panel_a | panel_b) / (panel_c_cumulative | panel_d) +
  plot_annotation(tag_levels = "a")

save_fig("figure_3new_combined", fig_combined, 11, 8)

message("Figure 3new plotting complete.")