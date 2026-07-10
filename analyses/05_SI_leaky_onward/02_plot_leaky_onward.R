# =============================================================================
# 02_plot_leaky_onward.R
#
# Reads the per-run metrics from 01_analysis_leaky_onward.R and produces the
# DEATHS-focused, HCW-headline reviewer figures. Aggregation order (as specified):
# median over the 10 replicates within each posterior draw, then median + 95%
# interval across the draws. No re-simulation here -- everything is derived from
# _intermediate/leaky_onward_per_run.rds.
#
# Figures (HCW deaths throughout):
#   1. fig_leaky_onward_hcw_averted        : HCW deaths averted by OBV vs r.
#   2. fig_leaky_onward_hcw_leaked         : additional HCW deaths that leak
#                                            through vs r (absolute).
#   3. fig_leaky_onward_hcw_pct_averted    : % of all HCW deaths (no-OBV world)
#                                            averted by OBV vs r.
#   4. fig_leaky_onward_hcw_leaked_frac    : additional (leaked) HCW deaths as a
#                                            % of all HCW deaths (no-OBV world).
# Plus leaky_onward_summary.csv and leaky_onward_analysis1_table.csv.
#
# On the "two simulations" behind the averted line: for each base run we replay
# the blocked HCWs (a) with no drug (A1 -> no_obv_deaths_hcw, the full benefit,
# r-independent) and (b) with the leaky drug (A2 -> leaky_deaths_accruing_hcw,
# per r). Averted(r) = A1 - A2(r). A1 is shared across r within a run, and the
# leaked HCW amount is small, so the averted curve is dominated by A1 and is not
# a noisy point-by-point difference (that was only an issue for TOTAL deaths,
# which we no longer plot). Figures 2 and 4 use A2 alone -- no differencing.
# =============================================================================

library(here)
library(ggplot2)
library(dplyr)
library(tidyr)

source(here("analyses", "05_SI_leaky_onward", "leaky_onward_helpers.R"))

INT_DIR <- here("analyses", "05_SI_leaky_onward", "_intermediate")
FIG_DIR <- here("figures")
CSV_DIR <- here("output_figgen")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CSV_DIR, recursive = TRUE, showWarnings = FALSE)

per_run <- readRDS(file.path(INT_DIR, "leaky_onward_per_run.rds"))

SC_LABELS <- c(DRC = "DRC (Middle)", WestAfrica = "West Africa (Worst)")
per_run$scenario_label <- ifelse(per_run$scenario %in% names(SC_LABELS),
                                 SC_LABELS[per_run$scenario], per_run$scenario)

# Full aggregated table (kept for reference / other metrics).
summ <- summarise_leaky_onward(per_run[, c("scenario", "particle", "rep", "r", "metric", "value")])
write.csv(summ, file.path(CSV_DIR, "leaky_onward_summary.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# Build a per-run wide table: r-independent scalars joined onto the per-r metrics,
# so every derived quantity (percentages, fractions) is formed PER RUN before any
# aggregation -- the statistically correct order.
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
    total_hcw_deaths_no_obv = tdf_deaths_hcw + no_obv_deaths_hcw,  # whole no-OBV epidemic, HCW
    hcw_averted             = net_deaths_averted_hcw,
    hcw_leaked              = leaky_deaths_accruing_hcw,
    pct_hcw_averted         = ifelse(total_hcw_deaths_no_obv > 0,
                                     100 * net_deaths_averted_hcw / total_hcw_deaths_no_obv, NA_real_),
    pct_hcw_leaked          = ifelse(total_hcw_deaths_no_obv > 0,
                                     100 * leaky_deaths_accruing_hcw / total_hcw_deaths_no_obv, NA_real_)
  )

# Aggregate one column: median over reps within draw, then median + 95% across draws.
agg_col <- function(df, col) {
  by_draw <- df %>%
    group_by(scenario, scenario_label, particle, r) %>%
    summarise(v = median(.data[[col]], na.rm = TRUE), .groups = "drop")
  by_draw %>%
    group_by(scenario, scenario_label, r) %>%
    summarise(median = median(v, na.rm = TRUE),
              lo95   = quantile(v, 0.025, na.rm = TRUE, names = FALSE),
              hi95   = quantile(v, 0.975, na.rm = TRUE, names = FALSE),
              .groups = "drop")
}

agg_averted    <- agg_col(run_wide, "hcw_averted")
agg_leaked     <- agg_col(run_wide, "hcw_leaked")
agg_pct_avert  <- agg_col(run_wide, "pct_hcw_averted")
agg_pct_leaked <- agg_col(run_wide, "pct_hcw_leaked")

# Reference line for Fig 1: the current index-only reported number (per scenario).
ref_reported <- scalars %>%
  group_by(scenario, scenario_label) %>%
  summarise(value = median(reported_index_deaths_hcw, na.rm = TRUE), .groups = "drop")

# -----------------------------------------------------------------------------
# Analysis 1 table (kept): current index-only vs true index+downstream averted.
# -----------------------------------------------------------------------------
a1 <- scalars %>%
  group_by(scenario, scenario_label) %>%
  summarise(reported_index_hcw = median(reported_index_deaths_hcw, na.rm = TRUE),
            total_averted_hcw   = median(no_obv_deaths_hcw, na.rm = TRUE),
            reported_index_total = median(reported_index_deaths_total, na.rm = TRUE),
            total_averted_total  = median(no_obv_deaths_total, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(downstream_hcw   = total_averted_hcw - reported_index_hcw,
         fold_hcw         = total_averted_hcw / reported_index_hcw)
write.csv(a1, file.path(CSV_DIR, "leaky_onward_analysis1_table.csv"), row.names = FALSE)
print(a1)

# -----------------------------------------------------------------------------
# Shared aesthetics
# -----------------------------------------------------------------------------
BLUE   <- "#2166AC"
ORANGE <- "#D6604D"
base_theme <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.text       = element_text(face = "bold"),
        plot.title       = element_text(face = "bold"),
        plot.subtitle    = element_text(colour = "grey30"),
        plot.title.position = "plot")
x_scale <- scale_x_continuous(breaks = seq(0, 1, 0.2), expand = expansion(mult = c(0.01, 0.02)))
xlab_r  <- "Residual transmissibility of a treated person  (r:  0 = current model → 1 = no effect on transmission)"

save_fig <- function(p, name, h = 4.2) {
  ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), p, width = 9, height = h)
  ggsave(file.path(FIG_DIR, paste0(name, ".png")), p, width = 9, height = h, dpi = 220)
}

# -----------------------------------------------------------------------------
# Figure 1: HCW deaths averted by OBV vs r
# -----------------------------------------------------------------------------
p1 <- ggplot(agg_averted, aes(r, median)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), fill = BLUE, alpha = 0.15) +
  geom_hline(data = ref_reported, aes(yintercept = value),
             linetype = "dashed", colour = "grey45", linewidth = 0.4) +
  geom_text(data = ref_reported, aes(x = 0.02, y = value, label = "currently reported (directly-treated only)"),
            hjust = 0, vjust = -0.6, size = 2.9, colour = "grey45") +
  geom_line(colour = BLUE, linewidth = 1) +
  facet_wrap(~ scenario_label, scales = "free_y") +
  x_scale +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
  labs(x = xlab_r, y = "HCW deaths averted by OBV",
       title = "OBV's averted HCW deaths hold up as onward transmission becomes leaky",
       subtitle = "Median (95% interval) across posterior draws; treated HCWs keep OBV's protection against death") +
  base_theme
save_fig(p1, "fig_leaky_onward_hcw_averted")

# -----------------------------------------------------------------------------
# Figure 2: additional HCW deaths that leak through vs r (absolute; A2 only)
# -----------------------------------------------------------------------------
p2 <- ggplot(agg_leaked, aes(r, median)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), fill = ORANGE, alpha = 0.15) +
  geom_line(colour = ORANGE, linewidth = 1) +
  facet_wrap(~ scenario_label, scales = "free_y") +
  x_scale +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
  labs(x = xlab_r, y = "Additional HCW deaths that leak through",
       title = "Additional HCW deaths a leaky drug allows through",
       subtitle = "Deaths among untreated onward HCWs only; the drug keeps protecting the downstream HCWs it reaches") +
  base_theme
save_fig(p2, "fig_leaky_onward_hcw_leaked")

# -----------------------------------------------------------------------------
# Figure 3: % of all HCW deaths (no-OBV world) averted by OBV vs r
# -----------------------------------------------------------------------------
p3 <- ggplot(agg_pct_avert, aes(r, median)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), fill = BLUE, alpha = 0.15) +
  geom_line(colour = BLUE, linewidth = 1) +
  facet_wrap(~ scenario_label) +
  x_scale +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08)),
                     labels = function(x) paste0(x, "%")) +
  labs(x = xlab_r, y = "HCW deaths averted (% of all HCW deaths without OBV)",
       title = "Share of HCW deaths OBV averts is stable across the leakiness range",
       subtitle = "Averted HCW deaths as a percentage of all HCW deaths in the no-OBV world") +
  base_theme
save_fig(p3, "fig_leaky_onward_hcw_pct_averted")

# -----------------------------------------------------------------------------
# Figure 4: leaked HCW deaths as % of all HCW deaths (no-OBV world)
# -----------------------------------------------------------------------------
p4 <- ggplot(agg_pct_leaked, aes(r, median)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), fill = ORANGE, alpha = 0.15) +
  geom_line(colour = ORANGE, linewidth = 1) +
  facet_wrap(~ scenario_label) +
  x_scale +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08)),
                     labels = function(x) paste0(x, "%")) +
  labs(x = xlab_r, y = "Leaked HCW deaths (% of all HCW deaths without OBV)",
       title = "Even with no effect on transmission, the leaked HCW deaths are a small fraction",
       subtitle = "Additional untreated HCW deaths as a percentage of all HCW deaths in the no-OBV world") +
  base_theme
save_fig(p4, "fig_leaky_onward_hcw_leaked_frac")

# Tidy CSVs of exactly what is plotted.
write.csv(bind_rows(
  agg_averted    %>% mutate(quantity = "hcw_deaths_averted"),
  agg_leaked     %>% mutate(quantity = "hcw_deaths_leaked"),
  agg_pct_avert  %>% mutate(quantity = "pct_hcw_deaths_averted"),
  agg_pct_leaked %>% mutate(quantity = "pct_hcw_deaths_leaked")
), file.path(CSV_DIR, "leaky_onward_hcw_figures.csv"), row.names = FALSE)

message("Wrote 4 HCW figures to ", FIG_DIR, " and CSVs to ", CSV_DIR)
