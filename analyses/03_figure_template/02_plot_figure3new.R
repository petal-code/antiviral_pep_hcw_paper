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
  annotate("text", x = 6.6,   y = 0.75, label = "Optimistic",  hjust = 0, vjust = 0, size = 3.2, color = "black") +
  annotate("text", x = 0.8, y = 0.35, label = "Pessimistic", hjust = 0, vjust = 0, size = 3.2, color = "black") +
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
  ggplot(df, aes(x = day, y = mean_val, color = arm_label, fill = arm_label)) +
    geom_ribbon(aes(ymin = pmax(ci_lo, 0), ymax = ci_hi), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = ARM_COLORS_C, name = NULL) +
    scale_fill_manual(values = ARM_COLORS_C, name = NULL) +
    scale_x_continuous(limits = c(0, 420), expand = c(0, 0)) +
    labs(x = "Day", y = y_label) +
    theme_fig()
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
  ggplot(df, aes(x = day, y = mean_val, color = arm_label, fill = arm_label)) +
    geom_ribbon(aes(ymin = pmax(ci_lo, 0), ymax = ci_hi), alpha = 0.15, color = NA) +
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
  ggplot(df, aes(x = day, y = mean_val, color = arm_label, fill = arm_label)) +
    geom_ribbon(aes(ymin = pmax(ci_lo, 0), ymax = ci_hi), alpha = 0.15, color = NA) +
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
  ggplot(aes(x = day, y = mean_val, color = arm_label, fill = arm_label)) +
  geom_ribbon(aes(ymin = pmax(ci_lo, 0), ymax = ci_hi), alpha = 0.15, color = NA) +
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
# 03_text_numbers_figure3new.R
#
# Pulls the specific summary numbers needed to fill in the XX placeholders
# in the Figure 3 results paragraph:
#
#   1. Cumulative HCW deaths (mean, 95% CI) at the end of follow-up, mid
#      efficacy, for three arms:
#        - achievable-impact / Ideal      (optimistic_mid)
#        - delay only, coverage preserved (dpc_conflict_mid)
#        - delay + reduced coverage       (with_conflict_mid)
#
#   2. % of the achievable HCW deaths averted that was realised under
#      conflict disruption (i.e. with_conflict's pct_hcw_deaths_averted as
#      a fraction of optimistic's pct_hcw_deaths_averted), computed per
#      particle and summarised as median + 95% CI, for each of the
#      optimistic/central/pessimistic delay-efficacy assumptions.
#
# Reads from the CSVs already produced by 02_extract_figure3new.R. Prints a
# flat list of the six numbers at the end, ready to copy/paste.
# =============================================================================
library(dplyr)
library(here)

ts_3new       <- read.csv(here("output_figgen", "figure_3new_weekly_ts.csv"))
particle_3new <- read.csv(here("output_figgen", "figure_3new_particle_summary.csv"))

EFF_ARM_LABELS <- c(hi = "Optimistic", mid = "Central", lo = "Pessimistic")
EFF_ARM_ORDER  <- c("hi", "mid", "lo")

fmt_n <- function(mean_val, lo, hi, digits = 1) {
  sprintf("%s (%s\u2013%s)",
          formatC(mean_val, format = "f", digits = digits),
          formatC(lo,        format = "f", digits = digits),
          formatC(hi,        format = "f", digits = digits))
}

fmt_pct <- function(med, lo, hi, digits = 1) {
  sprintf("%s%% (%s\u2013%s%%)",
          formatC(med, format = "f", digits = digits),
          formatC(lo,  format = "f", digits = digits),
          formatC(hi,  format = "f", digits = digits))
}

# =============================================================================
# 1. Cumulative HCW deaths at end of follow-up, mid efficacy, three arms
# =============================================================================
final_week <- max(
  ts_3new$week[ts_3new$scenario == "DRC" & ts_3new$metric == "hcw_deaths"],
  na.rm = TRUE
)

cum_deaths <- ts_3new %>%
  filter(scenario == "DRC", metric == "hcw_deaths", week == final_week,
         arm %in% c("optimistic_mid", "dpc_conflict_mid", "with_conflict_mid")) %>%
  select(arm, mean_val, ci_lo, ci_hi)

deaths_ideal <- cum_deaths %>% filter(arm == "optimistic_mid")
deaths_dpc   <- cum_deaths %>% filter(arm == "dpc_conflict_mid")
deaths_both  <- cum_deaths %>% filter(arm == "with_conflict_mid")

stopifnot(nrow(deaths_ideal) == 1, nrow(deaths_dpc) == 1, nrow(deaths_both) == 1)

txt_deaths_ideal <- fmt_n(deaths_ideal$mean_val, deaths_ideal$ci_lo, deaths_ideal$ci_hi)
txt_deaths_dpc   <- fmt_n(deaths_dpc$mean_val,   deaths_dpc$ci_lo,   deaths_dpc$ci_hi)
txt_deaths_both  <- fmt_n(deaths_both$mean_val,  deaths_both$ci_lo,  deaths_both$ci_hi)


# =============================================================================
# 2. % of achievable HCW deaths averted realised under conflict disruption,
#    by delay-efficacy assumption
# =============================================================================
realised_ratio <- particle_3new %>%
  filter(arm %in% c(
    "optimistic_hi", "optimistic_mid", "optimistic_lo",
    "with_conflict_hi", "with_conflict_mid", "with_conflict_lo"
  )) %>%
  mutate(
    eff_arm  = sub("^(optimistic|with_conflict)_", "", arm),
    cov_scen = sub("_(mid|lo|hi)$", "", arm)
  ) %>%
  select(particle_id, eff_arm, cov_scen, pct_hcw_deaths_averted) %>%
  tidyr::pivot_wider(names_from = cov_scen, values_from = pct_hcw_deaths_averted) %>%
  mutate(
    realised_pct_of_achievable = 100 * with_conflict / optimistic
  ) %>%
  filter(is.finite(realised_pct_of_achievable))

realised_summary <- realised_ratio %>%
  group_by(eff_arm) %>%
  summarise(
    median = median(realised_pct_of_achievable, na.rm = TRUE),
    ci_lo  = quantile(realised_pct_of_achievable, 0.025, na.rm = TRUE),
    ci_hi  = quantile(realised_pct_of_achievable, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(eff_arm_label = EFF_ARM_LABELS[eff_arm]) %>%
  arrange(match(eff_arm, EFF_ARM_ORDER))


txt_realised_hi  <- fmt_pct(
  realised_summary$median[realised_summary$eff_arm == "hi"],
  realised_summary$ci_lo[realised_summary$eff_arm == "hi"],
  realised_summary$ci_hi[realised_summary$eff_arm == "hi"]
)
txt_realised_mid <- fmt_pct(
  realised_summary$median[realised_summary$eff_arm == "mid"],
  realised_summary$ci_lo[realised_summary$eff_arm == "mid"],
  realised_summary$ci_hi[realised_summary$eff_arm == "mid"]
)
txt_realised_lo  <- fmt_pct(
  realised_summary$median[realised_summary$eff_arm == "lo"],
  realised_summary$ci_lo[realised_summary$eff_arm == "lo"],
  realised_summary$ci_hi[realised_summary$eff_arm == "lo"]
)

# =============================================================================
# 3. Final number list (copy/paste these)
# =============================================================================
cat("=== NUMBERS ===\n")
cat("Cumulative HCW deaths, Ideal (achievable):       ", txt_deaths_ideal, "\n")
cat("Cumulative HCW deaths, delay only:                ", txt_deaths_dpc, "\n")
cat("Cumulative HCW deaths, delay + reduced coverage:  ", txt_deaths_both, "\n")
cat("Realised % of achievable, Optimistic efficacy:    ", txt_realised_hi, "\n")
cat("Realised % of achievable, Central efficacy:       ", txt_realised_mid, "\n")
cat("Realised % of achievable, Pessimistic efficacy:   ", txt_realised_lo, "\n")