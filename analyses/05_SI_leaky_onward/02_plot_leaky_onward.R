# =============================================================================
# 02_plot_leaky_onward.R
#
# Reads the per-run metrics from 01_analysis_leaky_onward.R, aggregates them
# (median over the 10 replicates within each posterior draw, then median + 95%
# interval across the 200 draws), and produces the reviewer-response outputs.
# DEATHS are the focus throughout, with HCW deaths as the headline.
#
# Outputs (to output_figgen/ and figures/):
#   - leaky_onward_summary.csv          : the full aggregated table.
#   - leaky_onward_analysis1_table.csv  : Analysis 1 (conservatism of the current
#                                          index-only reporting).
#   - fig_leaky_onward_hcw_deaths.pdf/png : Analysis 2 headline -- HCW deaths
#                                          averted vs how leaky the drug is.
#   - fig_leaky_onward_accruing.pdf/png   : companion -- HCW deaths that leak
#                                          through, and those still averted
#                                          downstream by re-treatment.
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
summ    <- summarise_leaky_onward(per_run)
write.csv(summ, file.path(CSV_DIR, "leaky_onward_summary.csv"), row.names = FALSE)

SC_LABELS <- c(DRC = "DRC (Middle)", WestAfrica = "West Africa (Worst)")
lab_sc <- function(x) ifelse(x %in% names(SC_LABELS), SC_LABELS[x], x)

# Convenience: pull one r-independent metric's median per scenario.
scalar_of <- function(metric) {
  summ %>% filter(is.na(r), metric == !!metric) %>%
    select(scenario, value = median)
}

# =============================================================================
# Analysis 1 -- the current index-only reporting is conservative
#
# Reported (index-only) vs true averted (index + downstream, from the no-OBV
# counterfactual), for total and HCW deaths. The ratio shows how much larger the
# true averted burden is than what the paper currently credits OBV with.
# =============================================================================
a1 <- bind_rows(
  scalar_of("reported_index_deaths_total") %>% mutate(who = "Total", quantity = "reported_index"),
  scalar_of("reported_index_deaths_hcw")   %>% mutate(who = "HCW",   quantity = "reported_index"),
  scalar_of("no_obv_deaths_total")         %>% mutate(who = "Total", quantity = "total_averted"),
  scalar_of("no_obv_deaths_hcw")           %>% mutate(who = "HCW",   quantity = "total_averted")
) %>%
  pivot_wider(names_from = quantity, values_from = value) %>%
  mutate(downstream_averted = total_averted - reported_index,
         fold_vs_reported   = total_averted / reported_index,
         scenario_label     = lab_sc(scenario)) %>%
  arrange(scenario, who)

write.csv(a1, file.path(CSV_DIR, "leaky_onward_analysis1_table.csv"), row.names = FALSE)
message("Analysis 1 (deaths averted: reported index-only vs true index+downstream):")
print(a1)

# =============================================================================
# Analysis 2 -- HCW deaths averted vs how leaky OBV is (the headline)
#
# net_deaths_averted_hcw(r) = no-OBV HCW deaths - HCW deaths that still occur
# under a leaky drug (keeping OBV's death protection on everyone it treats).
# Shown against two references: the current reported (index-only) number, and
# the full sterilizing benefit (r = 0). The message: the HCW-death benefit holds
# up across the whole leakiness range and never falls below what we report.
# =============================================================================
make_net_panel <- function(who = c("hcw", "total")) {
  who <- match.arg(who)
  net_metric <- paste0("net_deaths_averted_", who)
  rep_metric <- paste0("reported_index_deaths_", who)

  net <- summ %>% filter(!is.na(r), metric == net_metric) %>%
    mutate(scenario_label = lab_sc(scenario))
  ref_reported <- scalar_of(rep_metric) %>% mutate(scenario_label = lab_sc(scenario))

  ggplot(net, aes(r, median)) +
    geom_ribbon(aes(ymin = lo95, ymax = hi95), alpha = 0.18, fill = "#2c7fb8") +
    geom_line(colour = "#2c7fb8", linewidth = 0.9) +
    geom_hline(data = ref_reported, aes(yintercept = value),
               linetype = "dashed", colour = "grey35") +
    geom_text(data = ref_reported, aes(x = 0.02, y = value,
              label = "current reported (index-only)"),
              hjust = 0, vjust = -0.5, size = 2.8, colour = "grey35") +
    facet_wrap(~ scenario_label, scales = "free_y") +
    scale_x_continuous(breaks = seq(0, 1, 0.2)) +
    labs(
      x = "Residual transmissibility of a treated person  (r:  0 = current model, 1 = no effect on transmission)",
      y = if (who == "hcw") "HCW deaths averted by OBV" else "Total deaths averted by OBV",
      title = if (who == "hcw")
        "OBV's averted HCW deaths are robust to leaky onward transmission"
      else "OBV's averted total deaths vs leaky onward transmission",
      subtitle = "Median (95% interval) across posterior draws; treated people keep OBV's protection against death"
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))
}

p_hcw <- make_net_panel("hcw")
ggsave(file.path(FIG_DIR, "fig_leaky_onward_hcw_deaths.pdf"), p_hcw, width = 9, height = 4.2)
ggsave(file.path(FIG_DIR, "fig_leaky_onward_hcw_deaths.png"), p_hcw, width = 9, height = 4.2, dpi = 200)

p_tot <- make_net_panel("total")
ggsave(file.path(FIG_DIR, "fig_leaky_onward_total_deaths.pdf"), p_tot, width = 9, height = 4.2)

# =============================================================================
# Companion -- where the leaked deaths go
#
# As OBV gets leakier, some HCW deaths leak through (the coverage/efficacy gaps),
# but the drug keeps saving downstream HCWs it re-treats. Plotting both shows the
# leaked deaths stay small precisely because of that continued downstream
# protection.
# =============================================================================
acc <- summ %>%
  filter(!is.na(r), metric %in% c("leaky_deaths_accruing_hcw",
                                  "leaky_deaths_averted_downstream_hcw")) %>%
  mutate(scenario_label = lab_sc(scenario),
         series = recode(metric,
                         leaky_deaths_accruing_hcw           = "HCW deaths that leak through",
                         leaky_deaths_averted_downstream_hcw = "HCW deaths still averted downstream (re-treatment)"))

p_acc <- ggplot(acc, aes(r, median, colour = series, fill = series)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ scenario_label, scales = "free_y") +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  scale_colour_manual(values = c("HCW deaths that leak through" = "#d95f0e",
                                 "HCW deaths still averted downstream (re-treatment)" = "#31a354")) +
  scale_fill_manual(values = c("HCW deaths that leak through" = "#d95f0e",
                               "HCW deaths still averted downstream (re-treatment)" = "#31a354")) +
  labs(x = "Residual transmissibility of a treated person (r)",
       y = "HCW deaths", colour = NULL, fill = NULL,
       title = "Where the leaked HCW deaths go",
       subtitle = "The drug keeps treating downstream HCWs, so the deaths that leak through stay small") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(file.path(FIG_DIR, "fig_leaky_onward_accruing.pdf"), p_acc, width = 9, height = 4.6)
ggsave(file.path(FIG_DIR, "fig_leaky_onward_accruing.png"), p_acc, width = 9, height = 4.6, dpi = 200)

message("Wrote figures to ", FIG_DIR, " and CSVs to ", CSV_DIR)
