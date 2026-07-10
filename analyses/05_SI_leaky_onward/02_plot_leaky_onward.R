# =============================================================================
# 02_plot_leaky_onward.R
#
# Reads the per-run metrics from 01_analysis_leaky_onward.R and produces a single
# 2x2 (a-d) HCW-deaths figure for the leaky-transmission sensitivity. Aggregation
# order: median over the 10 replicates within each posterior draw, then median +
# 95% interval across the draws. No re-simulation -- everything is derived from
# _intermediate/leaky_onward_per_run.rds.
#
# Panels (HCW deaths throughout; residual transmissibility r on the x-axis, where
# r = 0 is fully transmission-blocking and r = 1 is no effect on transmission):
#   a. HCW deaths averted by OBV.
#   b. Additional HCW deaths arising from residual transmission (untreated onward
#      cases; the leaky arm alone, no differencing).
#   c. HCW deaths averted as a % of all HCW deaths in the no-OBV world.
#   d. Additional HCW deaths as a % of all HCW deaths in the no-OBV world.
#
# Requires patchwork (2x2 assembly) in addition to ggplot2/dplyr/tidyr.
# =============================================================================

library(here)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

source(here("analyses", "05_SI_leaky_onward", "leaky_onward_helpers.R"))

INT_DIR <- here("analyses", "05_SI_leaky_onward", "_intermediate")
FIG_DIR <- here("figures")
CSV_DIR <- here("output_figgen")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CSV_DIR, recursive = TRUE, showWarnings = FALSE)

per_run <- readRDS(file.path(INT_DIR, "leaky_onward_per_run.rds"))

SC_LABELS <- c(DRC = "DRC", WestAfrica = "West Africa")
per_run$scenario_label <- ifelse(per_run$scenario %in% names(SC_LABELS),
                                 SC_LABELS[per_run$scenario], per_run$scenario)

# Full aggregated table (reference / other metrics).
summ <- summarise_leaky_onward(per_run[, c("scenario", "particle", "rep", "r", "metric", "value")])
write.csv(summ, file.path(CSV_DIR, "leaky_onward_summary.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# Per-run wide table: r-independent scalars joined onto the per-r metrics, so
# every percentage is formed PER RUN before aggregation.
# -----------------------------------------------------------------------------
scalars <- per_run %>%
  filter(is.na(r)) %>%
  select(scenario, scenario_label, particle, rep, metric, value) %>%
  pivot_wider(names_from = metric, values_from = value)

per_r <- per_run %>%
  filter(!is.na(r),
         metric %in% c("net_deaths_averted_hcw", "leaky_deaths_accruing_hcw")) %>%
  select(scenario, particle, rep, r, metric, value) %>%
  pivot_wider(names_from = metric, values_from = value)

run_wide <- per_r %>%
  left_join(scalars, by = c("scenario", "particle", "rep")) %>%
  mutate(
    total_hcw_deaths_no_obv = tdf_deaths_hcw + no_obv_deaths_hcw,   # whole no-OBV HCW toll
    hcw_averted    = net_deaths_averted_hcw,
    hcw_additional = leaky_deaths_accruing_hcw,                     # untreated onward HCW deaths
    pct_hcw_averted    = ifelse(total_hcw_deaths_no_obv > 0,
                                100 * net_deaths_averted_hcw / total_hcw_deaths_no_obv, NA_real_),
    pct_hcw_additional = ifelse(total_hcw_deaths_no_obv > 0,
                                100 * leaky_deaths_accruing_hcw / total_hcw_deaths_no_obv, NA_real_)
  )

# Median over reps within draw, then median + 95% across draws.
agg_col <- function(df, col) {
  df %>%
    group_by(scenario, scenario_label, particle, r) %>%
    summarise(v = median(.data[[col]], na.rm = TRUE), .groups = "drop") %>%
    group_by(scenario, scenario_label, r) %>%
    summarise(median = median(v, na.rm = TRUE),
              lo95   = quantile(v, 0.025, na.rm = TRUE, names = FALSE),
              hi95   = quantile(v, 0.975, na.rm = TRUE, names = FALSE),
              .groups = "drop")
}

agg_averted        <- agg_col(run_wide, "hcw_averted")
agg_additional     <- agg_col(run_wide, "hcw_additional")
agg_pct_averted    <- agg_col(run_wide, "pct_hcw_averted")
agg_pct_additional <- agg_col(run_wide, "pct_hcw_additional")

# Reference (dashed): the currently-reported effect = directly-treated HCWs only.
ref_reported <- scalars %>%
  group_by(scenario, scenario_label) %>%
  summarise(value = median(reported_index_deaths_hcw, na.rm = TRUE), .groups = "drop")

# Analysis-1 table (kept): current index-only vs true index+downstream averted.
a1 <- scalars %>%
  group_by(scenario, scenario_label) %>%
  summarise(reported_index_hcw = median(reported_index_deaths_hcw, na.rm = TRUE),
            total_averted_hcw   = median(no_obv_deaths_hcw, na.rm = TRUE), .groups = "drop") %>%
  mutate(downstream_hcw = total_averted_hcw - reported_index_hcw,
         fold_hcw       = total_averted_hcw / reported_index_hcw)
write.csv(a1, file.path(CSV_DIR, "leaky_onward_analysis1_table.csv"), row.names = FALSE)
print(a1)

# -----------------------------------------------------------------------------
# 2x2 figure
# -----------------------------------------------------------------------------
# Archetype colours (house palette, matching helper_functions_figure_1to4.R).
SCEN_COLS <- c("DRC" = "#1b9e77", "West Africa" = "#d95f02")

panel <- function(dat, ylab, pct = FALSE, show_x = FALSE, ref = NULL) {
  g <- ggplot(dat, aes(r, median, colour = scenario_label, fill = scenario_label)) +
    geom_ribbon(aes(ymin = lo95, ymax = hi95), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.85) +
    facet_wrap(~ scenario_label, scales = "free_y") +
    scale_colour_manual(values = SCEN_COLS) +
    scale_fill_manual(values = SCEN_COLS) +
    scale_x_continuous(breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0.01, 0.02))) +
    labs(x = if (show_x) "Residual transmissibility of treated cases, r" else NULL, y = ylab) +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          strip.background  = element_rect(fill = "grey92", colour = NA),
          strip.text        = element_text(face = "bold", size = 9),
          axis.title        = element_text(size = 9),
          legend.position   = "none")   # archetype is already the facet strip
  if (!is.null(ref))
    g <- g + geom_hline(data = ref, aes(yintercept = value),
                        linetype = "dashed", colour = "grey40", linewidth = 0.35)
  if (pct) g <- g + scale_y_continuous(limits = c(0, NA), labels = function(x) paste0(x, "%"),
                                       expand = expansion(mult = c(0, 0.08)))
  else     g <- g + scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08)))
  g
}

pa <- panel(agg_averted,        "HCW deaths averted by OBV", ref = ref_reported)
pb <- panel(agg_additional,     "Additional HCW deaths from residual transmission")
pc <- panel(agg_pct_averted,    "HCW deaths averted\n(% of HCW deaths without OBV)",
            pct = TRUE, show_x = TRUE)
pd <- panel(agg_pct_additional, "Additional HCW deaths\n(% of HCW deaths without OBV)",
            pct = TRUE, show_x = TRUE)

fig <- (pa | pb) / (pc | pd) +
  plot_annotation(
    tag_levels = "a",
    caption = paste0(
      "r = residual transmissibility of treated cases (0 = onward transmission fully blocked, the current model; ",
      "1 = transmission unaffected). Line = median, band = 95% interval across posterior draws. ",
      "Dashed line (a): currently reported effect (directly-treated HCWs only)."
    )
  ) &
  theme(plot.tag     = element_text(face = "bold", size = 12),
        plot.caption = element_text(hjust = 0, size = 7.5, colour = "grey30"))

ggsave(file.path(FIG_DIR, "fig_leaky_onward_hcw_2x2.pdf"), fig, width = 9, height = 6.4)
ggsave(file.path(FIG_DIR, "fig_leaky_onward_hcw_2x2.png"), fig, width = 9, height = 6.4, dpi = 300)

# Tidy CSV of exactly what is plotted.
write.csv(bind_rows(
  agg_averted        %>% mutate(panel = "a_hcw_deaths_averted"),
  agg_additional     %>% mutate(panel = "b_hcw_deaths_additional"),
  agg_pct_averted    %>% mutate(panel = "c_pct_hcw_deaths_averted"),
  agg_pct_additional %>% mutate(panel = "d_pct_hcw_deaths_additional")
), file.path(CSV_DIR, "leaky_onward_hcw_figure.csv"), row.names = FALSE)

message("Wrote fig_leaky_onward_hcw_2x2.{pdf,png} to ", FIG_DIR)
