# =============================================================================
# 02_plot_figure3new.R
# =============================================================================
source(here::here("analyses", "03_figure_template", "helper_functions_figure_1to4.R"))
source(here::here("functions", "setup_model_parameters.R"))
library(fiber)

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
ts_3new         <- read.csv(here("output_figgen", "figure_3new_weekly_ts.csv"))
particle_3new   <- read.csv(here("output_figgen", "figure_3new_particle_summary.csv"))
pep_uptake_3new <- read.csv(here("output_figgen", "figure_3new_pep_uptake_summary.csv"))
period_raw      <- read.csv(here("output_figgen", "figure_3new_period_summary.csv"))
curve_d50_dat   <- readRDS(here("data-processed", "DPC_fixed_efficacy_varied_d50.rds"))

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
# sdb$dpc_conflict      <- 1 + 9 * (1 - (sdb$value_tweaked / max(sdb$value_tweaked))^2)
sdb$dpc_conflict      <- 1 + 4 * (1 - (sdb$value_tweaked / max(sdb$value_tweaked)))

sub      <- sdb[sdb$day < 200, ]
peak_row <- sub[which.max(sub$coverage_conflict), ]
peak_day <- peak_row$day

sdb$dpc_conflict[sdb$day <= peak_day] <- 1

# Flat versions: unaffected dimension held constant at peak.
# NOTE: kept here for reference / potential reuse elsewhere, but no longer
# plotted in panel b -- panel b now only shows the conflict-scenario curves
# that are actually used in panel c.
sdb$coverage_flat <- sdb$coverage_conflict
sdb$coverage_flat[sdb$day > peak_day] <- peak_row$coverage_conflict
sdb$dpc_flat <- 1

# Shared constants
scale_factor <- 5.5 / 80  # cap DPC secondary axis at 5.5 when coverage axis tops out at 80
cov_color    <- "#E08214"
dpc_color    <- "black"
day_max_tv   <- max(ts_3new$week[ts_3new$scenario == "DRC"], na.rm = TRUE) * 7

# =============================================================================
# Shared label/color constants
# =============================================================================
EFF_ARM_LABELS <- c(hi = "Optimistic", mid = "Central", lo = "Pessimistic")
EFF_ARM_ORDER  <- c("hi", "mid", "lo")

COV_SCEN_LABELS <- c(
  optimistic    = "Ideal (100% coverage, 0 delay)",
  dpc_conflict  = "DPC impacted",
  with_conflict = "Both impacted"
)
COV_SCEN_ORDER  <- c("optimistic", "dpc_conflict", "with_conflict")
COV_SCEN_COLORS <- c(
  "Ideal (100% coverage, 0 delay)" = "#1a9641",
  "DPC impacted"                   = "#f58231",
  "Both impacted"                  = "#d7191c"
)

# Arms for panel c (mid efficacy, scenario comparison).
# "No PEP" now points at no_pep_mid, the matched-seed counterfactual
# (tdf + prevented from the with_conflict_mid runs themselves), rather
# than the separately-simulated no_pep arc. See 02_extract_figure3new.R.
ARM_LABELS_C <- c(
  no_pep_mid        = "No PEP",
  optimistic_mid    = "Ideal (100% coverage, 0 delay)",
  dpc_conflict_mid  = "DPC impacted",
  with_conflict_mid = "Both impacted"
)
ARM_COLORS_C <- c(
  "No PEP"                          = "black",
  "Ideal (100% coverage, 0 delay)"  = "#1a9641",
  "DPC impacted"                    = "#f58231",
  "Both impacted"                   = "#d7191c"
)

# Arms for panel bottom (three arms only) -- same matched-seed no_pep swap
ARM_LABELS_BOTTOM <- c(
  no_pep_mid         = "No PEP",
  optimistic_mid     = "Ideal (100% coverage, 0 delay)",
  with_conflict_mid  = "Both impacted"
)
ARM_COLORS_BOTTOM <- c(
  "No PEP"                          = "black",
  "Ideal (100% coverage, 0 delay)"  = "#1a9641",
  "Both impacted"                   = "#d7191c"
)

# Arms for efficacy comparison panel (both impacted, all three efficacy levels)
# NOTE: this panel still uses the single separately-simulated "no_pep" arc
# as its reference, since there is no single matched-seed no_pep series that
# applies across all three efficacy levels at once.
ARM_LABELS_EFF <- c(
  no_pep            = "No PEP",
  with_conflict_hi  = "Optimistic efficacy",
  with_conflict_mid = "Central efficacy",
  with_conflict_lo  = "Pessimistic efficacy"
)
ARM_COLORS_EFF <- c(
  "No PEP"               = "black",
  "Optimistic efficacy"  = "#1a9641",
  "Central efficacy"     = "#f58231",
  "Pessimistic efficacy" = "#d7191c"
)

# Colors for efficacy-faceted panels (scenario as lines within each panel)
EFF_SCENARIO_COLORS <- c(
  "No PEP"                          = "black",
  "Ideal (100% coverage, 0 delay)"  = "#1a9641",
  "DPC impacted"                    = "#f58231",
  "Both impacted"                   = "#d7191c"
)

# Period alpha levels for early/late boxplot distinction
PERIOD_ALPHA <- c("early (day 0-200)" = 0.35, "late (day 201+)" = 1.0)
PERIOD_LABELS <- c("Early (day 0-200)", "Late (day 201+)")

# =============================================================================
# Pre-process pep_uptake_3new: parse arm labels, drop cov_conflict
# =============================================================================
pep_uptake_3new <- pep_uptake_3new %>%
  mutate(
    eff_arm        = sub("^(optimistic|with_conflict|cov_conflict|dpc_conflict)_", "", arm),
    cov_scen       = sub("_(mid|lo|hi)$", "", arm),
    eff_arm_label  = factor(EFF_ARM_LABELS[eff_arm],   levels = EFF_ARM_LABELS[EFF_ARM_ORDER]),
    cov_scen_label = factor(COV_SCEN_LABELS[cov_scen], levels = COV_SCEN_LABELS[COV_SCEN_ORDER])
  ) %>%
  filter(!is.na(cov_scen_label))

# =============================================================================
# Pre-process period_raw:
#   1. Extract no_pep deaths as common denominator per particle x rep x period
#   2. Join onto PEP arms and compute pct_averted
#   3. Drop cov_conflict, add factor labels
# =============================================================================
no_pep_denom <- period_raw %>%
  filter(arm == "no_pep") %>%
  select(particle_id, rep, period, no_pep_deaths = n_hcw_deaths)

# period_3new <- period_raw %>%
#   filter(arm != "no_pep") %>%
#   left_join(no_pep_denom, by = c("particle_id", "rep", "period")) %>%
#   mutate(
#     pct_averted    = ifelse(
#       !is.na(no_pep_deaths) & no_pep_deaths > 0,
#       100 * (no_pep_deaths - n_hcw_deaths) / no_pep_deaths,
#       NA_real_
#     ),
#     eff_arm        = sub("^(optimistic|with_conflict|cov_conflict|dpc_conflict)_", "", arm),
#     cov_scen       = sub("_(mid|lo|hi)$", "", arm),
#     eff_arm_label  = factor(EFF_ARM_LABELS[eff_arm],   levels = EFF_ARM_LABELS[EFF_ARM_ORDER]),
#     cov_scen_label = factor(COV_SCEN_LABELS[cov_scen], levels = COV_SCEN_LABELS[COV_SCEN_ORDER]),
#     period         = factor(period, levels = c("early (day 0-200)", "late (day 201+)"))
#   ) %>%
#   filter(!is.na(cov_scen_label))
period_3new <- period_raw %>%
  filter(arm != "no_pep") %>%
  mutate(
    eff_arm        = sub("^(optimistic|with_conflict|cov_conflict|dpc_conflict)_", "", arm),
    cov_scen       = sub("_(mid|lo|hi)$", "", arm),
    eff_arm_label  = factor(EFF_ARM_LABELS[eff_arm],   levels = EFF_ARM_LABELS[EFF_ARM_ORDER]),
    cov_scen_label = factor(COV_SCEN_LABELS[cov_scen], levels = COV_SCEN_LABELS[COV_SCEN_ORDER]),
    period         = factor(period, levels = c("early (day 0-200)", "late (day 201+)"))
  ) %>%
  filter(!is.na(cov_scen_label))

# =============================================================================
# Panel a: efficacy vs DPC
# Both the hi and lo efficacy bounds are drawn as black dashed lines (instead
# of green/red) to avoid the colors being mistaken for a meaningful encoding.
# The region between them is shaded light grey, and each dashed line is
# labeled directly (Optimistic / Pessimistic) near its right-hand end.
# =============================================================================
panel_a <- ggplot(curve_d50_dat, aes(x = dpc)) +
  geom_ribbon(aes(ymin = eighty_efficacy_lo, ymax = eighty_efficacy_hi),
              fill = "grey80", alpha = 0.4) +
  geom_line(aes(y = eighty_efficacy_hi), color = "black", linetype = "dashed", linewidth = 1) +
  geom_line(aes(y = efficacy),            color = "black", linewidth = 1.2) +
  geom_line(aes(y = eighty_efficacy_lo), color = "black", linetype = "dashed", linewidth = 1) +
  # Fixed label positions matching the requested locations: "Optimistic" up and
  # to the right along the upper dashed curve, "Pessimistic" lower-left along
  # the lower dashed curve's early decline. Adjust x/y here if the underlying
  # curve shape changes and the labels drift off the lines.
  annotate("text", x = 8,   y = 0.75, label = "Optimistic",  hjust = 0, vjust = 0, size = 3.2, color = "black") +
  annotate("text", x = 1.5, y = 0.35, label = "Pessimistic", hjust = 0, vjust = 0, size = 3.2, color = "black") +
  scale_y_continuous(limits = c(0, NA), labels = scales::percent) +
  labs(x = "Days post-exposure (DPC)", y = "Efficacy") +
  theme_fig()

# =============================================================================
# Panel b: coverage & DPC over time
# Conflict period shaded from day 110 to 300 (fixed)
#
# Only the DPC/coverage curves actually used by panel c's scenarios are
# shown, and lines are colored to match panel c's scenario colors so the
# two panels read together directly:
#   - coverage_flat (coverage stays at its peak, unaffected by conflict)
#     is the coverage curve used by the "DPC impacted" scenario -> orange
#   - coverage_conflict (coverage declines during conflict) is the
#     coverage curve used by the "Both impacted" scenario -> red
#   - dpc_conflict (DPC curve under conflict) is shared identically by
#     both the "DPC impacted" and "Both impacted" scenarios, so it's kept
#     neutral black/dashed rather than tied to a single scenario color
# dpc_flat (the near-flat DPC curve used by the "Ideal" scenario) is
# intentionally NOT shown -- it's a flat horizontal line that adds no
# information here.
# DPC secondary axis capped at 5.5.
# =============================================================================
panel_b <- ggplot(sdb, aes(x = day)) +
  annotate("rect", xmin = 110, xmax = 300, ymin = -Inf, ymax = Inf,
           fill = "grey85", alpha = 0.6) +
  annotate("text", x = (110 + 300) / 2, y = 78,
           label = "conflict", size = 3.5, color = "grey30") +
  geom_line(aes(y = coverage_flat),
            color = COV_SCEN_COLORS[["DPC impacted"]], linetype = "solid", linewidth = 1.1) +
  geom_line(aes(y = coverage_conflict),
            color = COV_SCEN_COLORS[["Both impacted"]], linetype = "solid", linewidth = 1.1) +
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

# =============================================================================
# Panel c: HCW deaths time series (mid efficacy, scenario comparison)
# x-axis fixed to 420 days to match panel b
# =============================================================================
make_panel_c <- function(metric_name, y_label) {
  df <- ts_3new %>%
    filter(scenario == "DRC", arm %in% names(ARM_LABELS_C), metric == metric_name) %>%
    mutate(
      day       = week * 7,
      arm_label = factor(ARM_LABELS_C[arm], levels = ARM_LABELS_C)
    )
  ggplot(df, aes(x = day, y = q50, color = arm_label, fill = arm_label)) +
    geom_ribbon(aes(ymin = pmax(q25, 0), ymax = q75), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = ARM_COLORS_C, name = NULL) +
    scale_fill_manual(values = ARM_COLORS_C, name = NULL) +
    scale_x_continuous(limits = c(0, 420), expand = c(0, 0)) +
    labs(x = "Day", y = y_label) +
    theme_fig() +
    theme(legend.position = c(0.3, 0.85))
}

panel_c_incident   <- make_panel_c("hcw_deaths_incidence", "Mean weekly incident HCW deaths")
panel_c_cumulative <- make_panel_c("hcw_deaths",           "Mean cumulative HCW deaths")

# =============================================================================
# Panel c (efficacy variant): both impacted, all three efficacy levels
# =============================================================================
make_panel_c_eff <- function(metric_name, y_label) {
  df <- ts_3new %>%
    filter(scenario == "DRC", arm %in% names(ARM_LABELS_EFF), metric == metric_name) %>%
    mutate(
      day       = week * 7,
      arm_label = factor(ARM_LABELS_EFF[arm], levels = ARM_LABELS_EFF)
    )
  ggplot(df, aes(x = day, y = q50, color = arm_label, fill = arm_label)) +
    geom_ribbon(aes(ymin = pmax(q25, 0), ymax = q75), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = ARM_COLORS_EFF, name = NULL) +
    scale_fill_manual(values = ARM_COLORS_EFF, name = NULL) +
    scale_x_continuous(limits = c(0, 420), expand = c(0, 0)) +
    labs(x = "Day", y = y_label,
         title = "DRC archetype (both impacted, by efficacy)") +
    theme_fig()
}

panel_c_eff_incident   <- make_panel_c_eff("hcw_deaths_incidence", "Mean weekly incident HCW deaths")
panel_c_eff_cumulative <- make_panel_c_eff("hcw_deaths",           "Mean cumulative HCW deaths")

# =============================================================================
# Panel c (efficacy-faceted): one panel per efficacy level, scenarios as lines
# Uses the matched-seed no_pep_{eff} counterfactual for each efficacy level.
# =============================================================================
make_panel_c_by_eff <- function(metric_name, y_label, eff, eff_title) {
  arm_map <- setNames(
    c("No PEP",
      "Ideal (100% coverage, 0 delay)",
      "DPC impacted",
      "Both impacted"),
    c(paste0("no_pep_",       eff),
      paste0("optimistic_",    eff),
      paste0("dpc_conflict_",  eff),
      paste0("with_conflict_", eff))
  )
  df <- ts_3new %>%
    filter(scenario == "DRC", arm %in% names(arm_map), metric == metric_name) %>%
    mutate(
      day       = week * 7,
      arm_label = factor(arm_map[arm], levels = arm_map)
    )
  ggplot(df, aes(x = day, y = q50, color = arm_label, fill = arm_label)) +
    geom_ribbon(aes(ymin = pmax(q25, 0), ymax = q75), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = EFF_SCENARIO_COLORS, name = NULL) +
    scale_fill_manual(values = EFF_SCENARIO_COLORS, name = NULL) +
    scale_x_continuous(limits = c(0, 420), expand = c(0, 0)) +
    labs(x = "Day", y = y_label,
         title = paste0("DRC archetype (", eff_title, " efficacy)")) +
    theme_fig()
}

panel_c_hi_cum  <- make_panel_c_by_eff("hcw_deaths", "Mean cumulative HCW deaths", "hi",  "optimistic")
panel_c_mid_cum <- make_panel_c_by_eff("hcw_deaths", "Mean cumulative HCW deaths", "mid", "central")
panel_c_lo_cum  <- make_panel_c_by_eff("hcw_deaths", "Mean cumulative HCW deaths", "lo",  "pessimistic")

panel_c_hi_inc  <- make_panel_c_by_eff("hcw_deaths_incidence", "Mean weekly incident HCW deaths", "hi",  "optimistic")
panel_c_mid_inc <- make_panel_c_by_eff("hcw_deaths_incidence", "Mean weekly incident HCW deaths", "mid", "central")
panel_c_lo_inc  <- make_panel_c_by_eff("hcw_deaths_incidence", "Mean weekly incident HCW deaths", "lo",  "pessimistic")

# =============================================================================
# Panel d: decomposition of HCW deaths averted % by efficacy arm, shown as
# a waterfall/bridge chart (full Ideal bar -> DPC loss step -> Coverage loss
# step -> full Realised bar), per efficacy arm.
#
# Decomposition steps (computed at particle level before taking medians):
#   realised = with_conflict                       (both coverage and DPC impacted)
#   cov_loss = dpc_conflict - with_conflict        (additional loss from imperfect coverage)
#   dpc_loss = optimistic   - dpc_conflict         (additional loss from DPC delay)
#
# NOTE: this panel still uses pct_hcw_deaths_averted from particle_3new,
# which is computed against the separately-simulated no_pep arc (it is
# not affected by the matched-seed no_pep_{eff} swap used in panels b/c).
# =============================================================================
DECOMP_ARMS <- c(
  "optimistic_mid",    "optimistic_lo",    "optimistic_hi",
  "dpc_conflict_mid",  "dpc_conflict_lo",   "dpc_conflict_hi",
  "with_conflict_mid", "with_conflict_lo",  "with_conflict_hi"
)

decomp_particle <- particle_3new %>%
  filter(arm %in% DECOMP_ARMS) %>%
  mutate(
    eff_arm  = sub("^(optimistic|dpc_conflict|with_conflict)_", "", arm),
    cov_scen = sub("_(mid|lo|hi)$", "", arm)
  ) %>%
  select(particle_id, eff_arm, cov_scen, pct_hcw_deaths_averted) %>%
  pivot_wider(names_from = cov_scen, values_from = pct_hcw_deaths_averted) %>%
  mutate(
    realised = with_conflict,
    cov_loss = pmax(dpc_conflict - with_conflict, 0),
    dpc_loss = pmax(optimistic   - dpc_conflict,  0)
  )

decomp_summary <- decomp_particle %>%
  group_by(eff_arm) %>%
  summarise(
    realised_med = median(realised, na.rm = TRUE),
    cov_loss_med = median(cov_loss, na.rm = TRUE),
    dpc_loss_med = median(dpc_loss, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    scenario      = "DRC",
    eff_arm_label = factor(EFF_ARM_LABELS[eff_arm],
                           levels = EFF_ARM_LABELS[EFF_ARM_ORDER])
  )

DECOMP_COLORS <- c(
  "Realised"      = COV_SCEN_COLORS[["Ideal (100% coverage, 0 delay)"]],
  "Coverage loss" = COV_SCEN_COLORS[["Both impacted"]],
  "DPC loss"      = COV_SCEN_COLORS[["DPC impacted"]]
)

# =============================================================================
# Waterfall layout for panel d
#
# For each efficacy arm, four bars are drawn left to right:
#   1. "Ideal"          full bar, 0 -> optimistic value             (green)
#   2. "DPC loss"        floating bar, dpc_conflict -> optimistic    (orange)
#   3. "Coverage loss"   floating bar, realised -> dpc_conflict      (red)
#   4. "Realised"        full bar, 0 -> realised value               (green)
# Dotted horizontal segments bridge the top/bottom edges between
# consecutive bars, matching a standard waterfall/bridge chart. Bar borders
# are a thin grey (rather than black) so the dotted bridge line -- drawn as
# a separate, later layer -- reads clearly on top instead of disappearing
# into a heavy bar outline. Bar spacing is wider than bar width so each
# bridge segment has enough room between bars to actually be visible; each
# segment runs from one bar's edge to the very edge of the next bar.
# =============================================================================
WF_BAR_WIDTH   <- 0.7
WF_BAR_SPACING <- 1.4  # distance between consecutive bar centers within a group (> WF_BAR_WIDTH so the gap, and the bridge line in it, is visible)
WF_GROUP_GAP   <- 1.2  # extra horizontal gap inserted between efficacy-arm groups

build_waterfall_df <- function(summary_df, eff_arm_order = EFF_ARM_ORDER) {
  summary_df <- summary_df %>%
    mutate(
      optimistic_med   = realised_med + cov_loss_med + dpc_loss_med,
      dpc_conflict_med = realised_med + cov_loss_med
    )
  
  bars <- list()
  segs <- list()
  group_centers <- numeric(0)
  group_width   <- 3 * WF_BAR_SPACING
  
  for (i in seq_along(eff_arm_order)) {
    ea  <- eff_arm_order[i]
    row <- summary_df %>% filter(eff_arm == ea)
    if (nrow(row) == 0) next
    
    base_x <- (i - 1) * (group_width + WF_GROUP_GAP)
    x1 <- base_x; x2 <- base_x + WF_BAR_SPACING
    x3 <- base_x + 2 * WF_BAR_SPACING; x4 <- base_x + 3 * WF_BAR_SPACING
    group_centers <- c(group_centers, (x1 + x4) / 2)
    
    bars[[length(bars) + 1]] <- data.frame(
      eff_arm = ea, stage = "Ideal",
      xmin = x1 - WF_BAR_WIDTH / 2, xmax = x1 + WF_BAR_WIDTH / 2,
      ymin = 0, ymax = row$optimistic_med,
      fill_group = "Realised"
    )
    bars[[length(bars) + 1]] <- data.frame(
      eff_arm = ea, stage = "DPC loss",
      xmin = x2 - WF_BAR_WIDTH / 2, xmax = x2 + WF_BAR_WIDTH / 2,
      ymin = row$dpc_conflict_med, ymax = row$optimistic_med,
      fill_group = "DPC loss"
    )
    bars[[length(bars) + 1]] <- data.frame(
      eff_arm = ea, stage = "Coverage loss",
      xmin = x3 - WF_BAR_WIDTH / 2, xmax = x3 + WF_BAR_WIDTH / 2,
      ymin = row$realised_med, ymax = row$dpc_conflict_med,
      fill_group = "Coverage loss"
    )
    bars[[length(bars) + 1]] <- data.frame(
      eff_arm = ea, stage = "Realised",
      xmin = x4 - WF_BAR_WIDTH / 2, xmax = x4 + WF_BAR_WIDTH / 2,
      ymin = 0, ymax = row$realised_med,
      fill_group = "Realised"
    )
    
    # each segment spans exactly from one bar's edge to the adjacent bar's edge
    segs[[length(segs) + 1]] <- data.frame(
      x = x1 + WF_BAR_WIDTH / 2, xend = x2 - WF_BAR_WIDTH / 2,
      y = row$optimistic_med,   yend = row$optimistic_med
    )
    segs[[length(segs) + 1]] <- data.frame(
      x = x2 + WF_BAR_WIDTH / 2, xend = x3 - WF_BAR_WIDTH / 2,
      y = row$dpc_conflict_med, yend = row$dpc_conflict_med
    )
    segs[[length(segs) + 1]] <- data.frame(
      x = x3 + WF_BAR_WIDTH / 2, xend = x4 - WF_BAR_WIDTH / 2,
      y = row$realised_med,     yend = row$realised_med
    )
  }
  
  list(
    bars          = do.call(rbind, bars),
    segs          = do.call(rbind, segs),
    group_centers = group_centers
  )
}

make_panel_d <- function(sc) {
  df <- filter(decomp_summary, scenario == sc)
  wf <- build_waterfall_df(df)
  
  axis_labels_df <- data.frame(
    x     = wf$group_centers,
    label = EFF_ARM_LABELS[EFF_ARM_ORDER]
  )
  
  ggplot() +
    geom_rect(data = wf$bars,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_group),
              color = "grey70", linewidth = 0.3) +
    geom_segment(data = wf$segs,
                 aes(x = x, xend = xend, y = y, yend = yend),
                 linetype = "dotted", color = "black", linewidth = 0.7) +
    scale_fill_manual(values = DECOMP_COLORS, name = NULL) +
    scale_x_continuous(breaks = axis_labels_df$x, labels = axis_labels_df$label) +
    scale_y_continuous(limits = c(0, 100), labels = function(x) paste0(x, "%")) +
    labs(x = "Antiviral efficacy", y = "HCW deaths averted (%)") +
    theme_fig() +
    theme(axis.ticks.x = element_blank())
}

panel_d <- make_panel_d("DRC")


# =============================================================================
# Panel e: PEP uptake % by coverage scenario x efficacy arm
# =============================================================================
uptake_summary <- pep_uptake_3new %>%
  group_by(eff_arm_label, cov_scen_label) %>%
  summarise(
    median_pct = median(pct_received, na.rm = TRUE),
    lo_pct     = quantile(pct_received, 0.025, na.rm = TRUE),
    hi_pct     = quantile(pct_received, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

panel_e <- ggplot(uptake_summary, aes(x = eff_arm_label, y = median_pct, fill = cov_scen_label)) +
  geom_col(position = position_dodge(0.7), width = 0.6, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = lo_pct, ymax = hi_pct),
                position = position_dodge(0.7), width = 0.2, linewidth = 0.4) +
  scale_fill_manual(values = COV_SCEN_COLORS, name = NULL) +
  scale_y_continuous(limits = c(0, 100), labels = function(x) paste0(x, "%")) +
  labs(x = "Antiviral efficacy", y = "HCW receiving PEP (%)",
       title = "PEP uptake by scenario") +
  theme_fig()

# =============================================================================
# Panel f: DPC distribution among infected PEP recipients
# =============================================================================
panel_f <- ggplot(
  pep_uptake_3new %>% filter(!is.na(dpc_mean)),
  aes(x = cov_scen_label, y = dpc_median, fill = cov_scen_label)
) +
  geom_boxplot(outlier.shape = NA, width = 0.6, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = COV_SCEN_COLORS, name = NULL, guide = "none") +
  facet_wrap(~ eff_arm_label, nrow = 1) +
  labs(x = NULL, y = "DPC among infected PEP recipients (days)",
       title = "DPC distribution by scenario") +
  theme_fig() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# =============================================================================
# Period-stratified panels
# Early (day 0-200): lighter fill; Late (day 201+): full fill
# Each x tick shows two boxplots (early vs late) side by side
# =============================================================================
make_period_boxplot <- function(data, y_var, y_label, title_text) {
  ggplot(data,
         aes(x = cov_scen_label, y = .data[[y_var]],
             fill = cov_scen_label, alpha = period)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, color = "black", linewidth = 0.3,
                 position = position_dodge(0.75)) +
    scale_fill_manual(values = COV_SCEN_COLORS, name = NULL, guide = "none") +
    scale_alpha_manual(values = PERIOD_ALPHA, name = NULL, labels = PERIOD_LABELS) +
    facet_wrap(~ eff_arm_label, nrow = 1) +
    labs(x = NULL, y = y_label, title = title_text) +
    theme_fig() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

panel_dpc_period <- make_period_boxplot(
  period_3new %>% filter(!is.na(dpc_median)),
  y_var      = "dpc_median",
  y_label    = "DPC among PEP recipients (days)",
  title_text = "DPC distribution by period, efficacy, and scenario"
)

panel_deaths_period <- make_period_boxplot(
  period_3new,
  y_var      = "n_hcw_deaths",
  y_label    = "HCW deaths",
  title_text = "HCW deaths by period, efficacy, and scenario"
)

panel_averted_period <- make_period_boxplot(
  period_3new %>% filter(!is.na(pct_averted)),
  y_var      = "pct_averted",
  y_label    = "HCW deaths averted (%)",
  title_text = "HCW deaths averted by period, efficacy, and scenario"
) +
  scale_y_continuous(labels = function(x) paste0(x, "%"))

# =============================================================================
# Panel top: time-varying parameters + SDB (normalised)
# =============================================================================
mp <- make_model_parameters(
  scenario_id     = "Middle_DRC_ConflictSmoothed_PlusPlus",
  scenario_matrix = read_scenario_matrix(here("data-processed",
                                              "final_six_scenario_values_original_approach.csv")),
  overrides       = list(seeding_cases = 25L, check_final_size = 10000)
)

TV_LABELS <- c(
  prob_hospitalised_genPop     = "P(hospitalised), general population",
  prob_hospitalised_hcw        = "P(hospitalised), HCW",
  hospitalisation_delay_factor = "Hospitalisation delay factor",
  p_unsafe_funeral_comm_genPop = "P(unsafe funeral, community), general population",
  p_unsafe_funeral_comm_hcw   = "P(unsafe funeral, community), HCW",
  p_unsafe_funeral_hosp_genPop = "P(unsafe funeral, hospital), general population",
  p_unsafe_funeral_hosp_hcw   = "P(unsafe funeral, hospital), HCW",
  prop_etu                     = "ETU proportion",
  ppe_coverage_hcw             = "PPE coverage, HCW"
)

TV_COLORS <- c(
  "P(hospitalised), general population"               = "#e6194b",
  "P(hospitalised), HCW"                               = "#3cb44b",
  "Hospitalisation delay factor"                       = "#4363d8",
  "P(unsafe funeral, community), general population"   = "#f58231",
  "P(unsafe funeral, community), HCW"                  = "#911eb4",
  "P(unsafe funeral, hospital), general population"    = "#42d4f4",
  "P(unsafe funeral, hospital), HCW"                   = "#f032e6",
  "ETU proportion"                                      = "#9A6324",
  "PPE coverage, HCW"                                   = "#000000",
  "SDB blended success (conflict driver)"               = "#bcf60c"
)

t_seq_tv <- seq(0, day_max_tv, by = 1)

tv_long <- do.call(rbind, lapply(names(mp$tv_args), function(nm) {
  fn      <- mp$tv_args[[nm]]
  raw_val <- fn(t_seq_tv)
  max_val <- max(raw_val, na.rm = TRUE)
  data.frame(
    day        = t_seq_tv,
    value_norm = if (max_val > 0) raw_val / max_val else raw_val,
    parameter  = TV_LABELS[nm],
    stringsAsFactors = FALSE
  )
}))

sdb_tv <- data.frame(
  day        = sdb$day,
  value_norm = sdb$value_tweaked / max(sdb$value_tweaked, na.rm = TRUE),
  parameter  = "SDB blended success (conflict driver)",
  stringsAsFactors = FALSE
)

panel_top <- ggplot(bind_rows(tv_long, sdb_tv),
                    aes(x = day, y = value_norm, color = parameter)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = TV_COLORS) +
  scale_x_continuous(limits = c(0, day_max_tv), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Day", y = "Normalised value (0-1)",
       title = "Time-varying parameters (DRC archetype, conflict)") +
  theme_fig() +
  theme(legend.position = "right", legend.text = element_text(size = 7)) +
  guides(color = guide_legend(ncol = 1))

# Panel bottom: weekly HCW deaths (three arms)
panel_bottom <- ts_3new %>%
  filter(scenario == "DRC", arm %in% names(ARM_LABELS_BOTTOM),
         metric == "hcw_deaths_incidence") %>%
  mutate(
    day       = week * 7,
    arm_label = factor(ARM_LABELS_BOTTOM[arm], levels = ARM_LABELS_BOTTOM)
  ) %>%
  ggplot(aes(x = day, y = q50, color = arm_label, fill = arm_label)) +
  geom_ribbon(aes(ymin = pmax(q25, 0), ymax = q75), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = ARM_COLORS_BOTTOM, name = NULL) +
  scale_fill_manual(values = ARM_COLORS_BOTTOM, name = NULL) +
  scale_x_continuous(limits = c(0, day_max_tv), expand = c(0, 0)) +
  labs(x = "Day", y = "Mean weekly incident HCW deaths") +
  theme_fig()

# =============================================================================
# Save all panels
# =============================================================================
save_fig("figure_3new_panel_a_efficacy",         panel_a,            5,   4)
save_fig("figure_3new_panel_b_coverage_dpc",     panel_b,            6,   4)
save_fig("figure_3new_panel_c_incident",         panel_c_incident,   6,   4)
save_fig("figure_3new_panel_c_cumulative",       panel_c_cumulative, 6,   4)
save_fig("figure_3new_panel_d_decomposition",    panel_d,            5,   4)
save_fig("figure_3new_panel_e_uptake",           panel_e,            7,   4.5)
save_fig("figure_3new_panel_f_dpc_distribution", panel_f,            10,  4.5)

save_fig("figure_3new_panel_c_eff_incident",     panel_c_eff_incident,   6, 4)
save_fig("figure_3new_panel_c_eff_cumulative",   panel_c_eff_cumulative, 6, 4)

save_fig("figure_3new_panel_c_by_eff_cumulative",
         panel_c_hi_cum | panel_c_mid_cum | panel_c_lo_cum, 15, 4)
save_fig("figure_3new_panel_c_by_eff_incident",
         panel_c_hi_inc | panel_c_mid_inc | panel_c_lo_inc, 15, 4)

save_fig("figure_3new_dpc_by_period",     panel_dpc_period,     12, 6)
save_fig("figure_3new_deaths_by_period",  panel_deaths_period,  12, 6)
save_fig("figure_3new_averted_by_period", panel_averted_period, 12, 6)

save_fig("figure_3new_tv_params_and_deaths",
         panel_top / panel_bottom +
           plot_layout(heights = c(1, 1)) +
           plot_annotation(tag_levels = "a"),
         7, 7)

save_fig("figure_3new_combined",
         (panel_a | panel_b) / (panel_c_cumulative | panel_d) +
           plot_annotation(tag_levels = "a"),
         11, 8)

# Alternate combined version: panel c shows incident (rather than cumulative)
# HCW deaths. Appended after the original combined figure.
save_fig("figure_3new_combined_incident",
         (panel_a | panel_b) / (panel_c_incident | panel_d) +
           plot_annotation(tag_levels = "a"),
         11, 8)

message("Figure 3new plotting complete.")

panel_recip_deaths_period <- make_period_boxplot(
  period_3new,
  y_var      = "n_pep_died",
  y_label    = "HCW deaths among PEP recipients",
  title_text = "HCW deaths among PEP recipients by period, efficacy, and scenario"
)

panel_recip_averted_period <- make_period_boxplot(
  period_3new %>% filter(!is.na(pct_averted_recipients)),
  y_var      = "pct_averted_recipients",
  y_label    = "Deaths averted among PEP recipients (%)",
  title_text = "Deaths averted among PEP recipients by period, efficacy, and scenario"
) +
  scale_y_continuous(limits = c(0, 100),
                     labels = function(x) paste0(x, "%"))

save_fig("figure_3new_recip_deaths_by_period",  panel_recip_deaths_period,  12, 6)
save_fig("figure_3new_recip_averted_by_period", panel_recip_averted_period, 12, 6)

# =============================================================================
# SI export: efficacy curve (old panel a) + coverage/DPC trajectory (old panel b)
# =============================================================================
save_fig("figure_3new_SIexport",
         (panel_a | panel_b) +
           plot_annotation(tag_levels = "a"),
         11, 4)

# =============================================================================
# figure_3new_combined_final
#
# Layout: (new_panel_a | new_panel_b) / new_panel_c
#   widths  c(2, 1) -- panel a takes 2/3, panel b takes 1/3
#   heights c(2, 1) -- top row takes 2/3, bottom row takes 1/3
#
# new_panel_a (= panel_c_cumulative): cumulative HCW deaths time series,
#   central efficacy, four scenarios.
#
# new_panel_b: stacked bar chart, central efficacy, four scenarios.
#   Each bar's total height = No PEP cumulative deaths (fixed denominator),
#   so the grey top section = deaths averted, coloured bottom = deaths that
#   still occurred. Scenarios left to right: No PEP, Both impacted,
#   DPC impacted, Ideal.
#
# new_panel_c: % HCW deaths averted by efficacy assumption (Optimistic /
#   Central / Pessimistic), three scenarios side by side per facet
#   (Ideal / DPC impacted / Both impacted), bars = median,
#   error bars = 2.5-97.5 percentile range.
# =============================================================================

# --- new_panel_b data: cumulative deaths at end of follow-up, mid efficacy ---
final_week_b <- max(
  ts_3new$week[ts_3new$scenario == "DRC" & ts_3new$metric == "hcw_deaths"],
  na.rm = TRUE
)

PANEL_B_ARMS <- c(
  no_pep_mid        = "No PEP",
  with_conflict_mid = "Both impacted",
  dpc_conflict_mid  = "DPC impacted",
  optimistic_mid    = "Ideal"
)

PANEL_B_COLORS_DIED <- c(
  "No PEP"        = "grey60",
  "Both impacted" = "#d7191c",
  "DPC impacted"  = "#f58231",
  "Ideal"         = "#1a9641"
)

cum_end <- ts_3new %>%
  filter(scenario == "DRC", metric == "hcw_deaths", week == final_week_b,
         arm %in% names(PANEL_B_ARMS)) %>%
  mutate(arm_label = factor(PANEL_B_ARMS[arm], levels = PANEL_B_ARMS))

# No PEP mean deaths = fixed total bar height for all bars
no_pep_total <- cum_end$q50[cum_end$arm == "no_pep_mid"]

# Stack order: averted (grey) on bottom, died (scenario color) on top.
# This makes the "floor" of the colored section a visual baseline for
# comparing how much was averted across scenarios.
panel_b_bars <- cum_end %>%
  mutate(
    died    = q50,
    averted = no_pep_total - q50
  ) %>%
  select(arm_label, died, averted) %>%
  tidyr::pivot_longer(cols = c(averted, died),
                      names_to = "segment", values_to = "value") %>%
  mutate(segment = factor(segment, levels = c("averted", "died")))

# All averted segments and No PEP died segment use the same grey so the
# "averted" region reads as one unified neutral backdrop across all bars.
new_panel_b <- ggplot(panel_b_bars,
                      aes(x = arm_label, y = value, fill = interaction(segment, arm_label))) +
  geom_col(width = 1.0, color = "grey50", linewidth = 0.2) +
  scale_fill_manual(
    values = c(
      "averted.No PEP"        = "grey75",
      "averted.Both impacted" = "grey75",
      "averted.DPC impacted"  = "grey75",
      "averted.Ideal"         = "grey75",
      "died.No PEP"           = "grey75",  # No PEP: same grey throughout
      "died.Both impacted"    = "#d7191c",
      "died.DPC impacted"     = "#f58231",
      "died.Ideal"            = "#1a9641"
    ),
    guide = "none"
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Cumulative HCW deaths") +
  theme_fig() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

# --- new_panel_c: % averted across all three efficacy assumptions -----------
PANEL_C_SCEN_ORDER  <- c("Ideal", "DPC impacted", "Both impacted")
PANEL_C_SCEN_COLORS <- c(
  "Ideal"         = "#1a9641",
  "DPC impacted"  = "#f58231",
  "Both impacted" = "#d7191c"
)

PANEL_C_ARM_MAP <- c(
  optimistic_hi    = "Ideal",
  optimistic_mid   = "Ideal",
  optimistic_lo    = "Ideal",
  dpc_conflict_hi  = "DPC impacted",
  dpc_conflict_mid = "DPC impacted",
  dpc_conflict_lo  = "DPC impacted",
  with_conflict_hi = "Both impacted",
  with_conflict_mid = "Both impacted",
  with_conflict_lo = "Both impacted"
)

panel_c_averted <- particle_3new %>%
  filter(arm %in% names(PANEL_C_ARM_MAP)) %>%
  mutate(
    scenario_label = factor(PANEL_C_ARM_MAP[arm], levels = PANEL_C_SCEN_ORDER),
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

new_panel_c <- ggplot(panel_c_averted,
                      aes(x = scenario_label, y = median,
                          fill = scenario_label, ymin = lo, ymax = hi)) +
  geom_col(position = position_dodge(0.7), width = 0.6,
           color = "grey70", linewidth = 0.3) +
  geom_errorbar(position = position_dodge(0.7), width = 0.2, linewidth = 0.4) +
  scale_fill_manual(values = PANEL_C_SCEN_COLORS, name = NULL) +
  scale_y_continuous(limits = c(0, 100), labels = function(x) paste0(x, "%"),
                     expand = c(0, 0)) +
  facet_wrap(~ eff_arm_label, nrow = 1) +
  labs(x = NULL, y = "HCW deaths averted (%)") +
  theme_fig() +
  theme(
    legend.position  = "none",
    axis.text.x      = element_text(angle = 25, hjust = 1),
    strip.background = element_blank(),
    strip.text       = element_text(size = 10, face = "plain")
  )

# --- assemble and save ------------------------------------------------------
# Shared y upper limit for panels a and b: take the max of No PEP q50
# across all time points so both panels share the same scale.
y_max_ab <- max(
  ts_3new$q75[ts_3new$scenario == "DRC" &
                ts_3new$metric   == "hcw_deaths" &
                ts_3new$arm      == "no_pep_mid"],
  na.rm = TRUE
) * 1.0  # 5% headroom

# Rebuild panel a with matching y range, corrected x label, and tighter margin
panel_a_final <- make_panel_c("hcw_deaths", "Mean cumulative HCW deaths") +
  scale_y_continuous(limits = c(0, y_max_ab), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 600), expand = c(0, 0)) +
  labs(x = "Days since outbreak start") +
  theme(plot.margin = margin(5, 10, 2, 5))

# Rebuild panel b with the same y upper limit
new_panel_b_final <- new_panel_b +
  scale_y_continuous(limits = c(0, y_max_ab), expand = c(0, 0))

# Use wrap_plots with explicit area layout to enforce 2:1 width ratio between
# panel a and panel b in the top row. The simple (a | b) / c patchwork
# syntax ignores widths when the bottom panel spans both columns.
# plot_annotation tag_levels applied after wrapping so labels a/b/c appear.
top_row <- wrap_plots(panel_a_final, new_panel_b_final, widths = c(2, 1))

final_fig <- wrap_plots(top_row, new_panel_c, ncol = 1, heights = c(2, 1)) +
  plot_annotation(tag_levels = "a")

save_fig("figure_3new_combined_final", final_fig, 14, 9)

