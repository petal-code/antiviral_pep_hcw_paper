# =============================================================================
# 02_plot_figure5.R
# Dose efficiency: doses per death averted under Policy A vs Policy B
# Reads pre-computed CSVs from output_figgen/
# Run 02_extract_figure5.R first.
# =============================================================================
source(here::here("analyses", "03_figure_template", "helper_functions_figure_1to4.R"))
OUT_DIR <- here("figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OBV_EFFICACY_A   <- 0.80
DOSES_PER_COURSE <- 20L

FIG5_EFFICACY_LEVELS <- c("obv_20", "obv_30", "obv_40",
                          "obv_50", "obv_60", "obv_70", "obv_80")

dose_df       <- read.csv(here("output_figgen", "figure_5_dose_summary.csv"),
                          stringsAsFactors = FALSE)
ppe_particles <- read.csv(here("output_figgen", "figure_5_ppe_by_particle.csv"),
                          stringsAsFactors = FALSE)

particle_df <- dose_df %>%
  group_by(scenario, particle_id, eff_name, efficacy) %>%
  summarise(doses_B         = mean(doses_B),
            n_prevented_hcw = mean(n_prevented_hcw),
            .groups = "drop") %>%
  left_join(ppe_particles, by = c("scenario", "particle_id")) %>%
  mutate(
    averted   = n_prevented_hcw * efficacy,
    base_pool = ifelse((1 - ppe_efficacy) > 0,
                       doses_B / (1 - ppe_efficacy), doses_B),
    dda_B = ifelse(averted > 0, doses_B   * DOSES_PER_COURSE / averted, NA_real_),
    dda_A = ifelse(averted > 0, base_pool * DOSES_PER_COURSE / averted, NA_real_),
    scenario_label = factor(SCENARIO_LABELS[scenario], levels = SCENARIO_LABELS)
  )

save_figure_data(particle_df, "figure_5_particle_df.csv")

# Color scheme
sc_fills <- c(
  setNames(c("#fdd0a2", unname(SCENARIO_COLORS["WestAfrica"])),
           paste(SCENARIO_LABELS["WestAfrica"], c("A", "B"), sep = ".")),
  setNames(c("#a8ddb5", unname(SCENARIO_COLORS["DRC"])),
           paste(SCENARIO_LABELS["DRC"],        c("A", "B"), sep = "."))
)

# Panel a: DDA boxplot at 80% efficacy
panel_a_df <- particle_df %>%
  filter(eff_name == "obv_80") %>%
  select(scenario, particle_id, scenario_label, A = dda_A, B = dda_B) %>%
  tidyr::pivot_longer(cols = c(A, B), names_to = "policy", values_to = "value") %>%
  mutate(group = paste(as.character(scenario_label), policy, sep = ".")) %>%
  filter(!is.na(value))

panel_a <- ggplot(panel_a_df, aes(x = scenario_label, y = value, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, color = "black",
               linewidth = 0.25, position = position_dodge(0.75)) +
  scale_fill_manual(
    values = sc_fills,
    breaks = c(paste(SCENARIO_LABELS["WestAfrica"], "A", sep = "."),
               paste(SCENARIO_LABELS["WestAfrica"], "B", sep = "."),
               paste(SCENARIO_LABELS["DRC"],        "A", sep = "."),
               paste(SCENARIO_LABELS["DRC"],        "B", sep = ".")),
    labels = c("West Africa — Policy A", "West Africa — Policy B",
               "DRC — Policy A",         "DRC — Policy B"),
    name = NULL
  ) +
  scale_y_continuous(
    limits = c(0, 350),
    breaks = seq(0, 350, by = 25),
    expand = expansion(mult = c(0, 0.05))
  ) +
  coord_cartesian(ylim = c(0, 225)) +
  scale_x_discrete(labels = c("West Africa (Worst)"    = "West Africa\narchetype",
                              "DRC (Middle, PlusPlus)" = "DRC\narchetype")) +
  labs(x = NULL, y = sprintf("Doses per death averted\n(1 course = %d doses)",
                             DOSES_PER_COURSE)) +
  theme_fig() +
  theme(legend.position = "bottom", legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size = 8)) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

# Panel b: DDA vs OBV efficacy (20-80%)
sweep_df <- particle_df %>%
  group_by(scenario, scenario_label, efficacy) %>%
  summarise(med_dda_A = median(dda_A, na.rm = TRUE),
            med_dda_B = median(dda_B, na.rm = TRUE),
            q25_dda_B = quantile(dda_B, 0.25, na.rm = TRUE),
            q75_dda_B = quantile(dda_B, 0.75, na.rm = TRUE),
            .groups = "drop")

policy_a_lines <- particle_df %>%
  filter(eff_name == "obv_80") %>%
  group_by(scenario, scenario_label) %>%
  summarise(med_dda_A = median(dda_A, na.rm = TRUE), .groups = "drop")

panel_b <- ggplot() +
  geom_ribbon(data = sweep_df,
              aes(x = efficacy * 100, ymin = q25_dda_B, ymax = q75_dda_B,
                  fill = scenario_label), alpha = 0.2) +
  geom_line(data = sweep_df,
            aes(x = efficacy * 100, y = med_dda_B, color = scenario_label),
            linewidth = 1.0) +
  geom_hline(data = policy_a_lines,
             aes(yintercept = med_dda_A, color = scenario_label),
             linetype = "dashed", linewidth = 1.0) +
  scale_color_manual(values = setNames(SCENARIO_COLORS, SCENARIO_LABELS), name = NULL) +
  scale_fill_manual( values = setNames(SCENARIO_COLORS, SCENARIO_LABELS), name = NULL) +
  scale_x_continuous(breaks = seq(20, 80, by = 10),
                     labels = function(x) paste0(x, "%"),
                     name   = "Antiviral efficacy (Policy B)") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                     name   = sprintf("Doses per death averted\n(1 course = %d doses)",
                                      DOSES_PER_COURSE)) +
  theme_fig() +
  theme(legend.position  = "none",
        panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
        panel.grid.minor = element_line(color = "grey95", linewidth = 0.2))

fig5_all <- (panel_a | panel_b) +
  plot_layout(widths = c(2, 3)) +
  plot_annotation(tag_levels = list(c("a ", "b ")))

save_fig <- function(filename_base, plot, width, height) {
  ggsave(file.path(OUT_DIR, paste0(filename_base, ".png")),
         plot, width = width, height = height, dpi = 400, units = "in")
  ggsave(file.path(OUT_DIR, paste0(filename_base, ".pdf")),
         plot, width = width, height = height, units = "in")
}

save_fig("figure_5", fig5_all, 12, 5)
message("Figure 5 saved")

# =============================================================================
# Split versions: panel a and panel b saved as separate figures
# =============================================================================
# fig5_panel_a <- panel_a +
#   plot_annotation(tag_levels = list(c("a ")))
# save_fig("figure_5_panel-a_dose-efficiency-boxplot", fig5_panel_a, 5, 5)
# 
# fig5_panel_b <- panel_b +
#   plot_annotation(tag_levels = list(c("b ")))
# save_fig("figure_5_panel-b_dose-efficiency-vs-efficacy", fig5_panel_b, 7, 5)

# fig5_panel_a <- panel_a
<<<<<<< HEAD
save_fig("figure_5_panel-a_dose-efficiency-boxplot", panel_a, 4, 4)

# fig5_panel_b <- panel_b
save_fig("figure_5_panel-b_dose-efficiency-vs-efficacy", panel_b, 6, 4)
=======
save_fig("figure_5_panel-a_dose-efficiency-boxplot", panel_a, 5, 5)

# fig5_panel_b <- panel_b
save_fig("figure_5_panel-b_dose-efficiency-vs-efficacy", panel_b, 7, 5)
>>>>>>> 873ecc12b709e76c5085cd6ebf2f57c289f1da8c

message("Figure 5 split panels saved")
# =============================================================================
# Panel a, alternative split: Policy A boxplots in one panel,
# Policy B boxplots in another panel (separate y-scales)
# =============================================================================
panel_a_A <- panel_a_df %>%
  filter(policy == "A") %>%
  ggplot(aes(x = scenario_label, y = value, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, color = "black",
               linewidth = 0.25, position = position_dodge(0.75)) +
  scale_fill_manual(
    values = sc_fills,
    breaks = c(paste(SCENARIO_LABELS["WestAfrica"], "A", sep = "."),
               paste(SCENARIO_LABELS["DRC"],        "A", sep = ".")),
    guide = "none"
  ) +
  scale_y_continuous(limits = c(0, 225), breaks = seq(0, 225, by = 25),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_x_discrete(labels = c("West Africa (Worst)"    = "West Africa\narchetype",
                              "DRC (Middle, PlusPlus)" = "DRC\narchetype")) +
  labs(x = NULL, y = sprintf("Doses per death averted\n(1 course = %d doses)",
                             DOSES_PER_COURSE),
       title = "Policy A") +
  theme_fig() +
  theme(legend.position = "none",
        plot.title = element_text(size = 11, hjust = 0.5))

panel_a_B <- panel_a_df %>%
  filter(policy == "B") %>%
  ggplot(aes(x = scenario_label, y = value, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, color = "black",
               linewidth = 0.25, position = position_dodge(0.75)) +
  scale_fill_manual(
    values = sc_fills,
    breaks = c(paste(SCENARIO_LABELS["WestAfrica"], "B", sep = "."),
               paste(SCENARIO_LABELS["DRC"],        "B", sep = ".")),
    guide = "none"
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_x_discrete(labels = c("West Africa (Worst)"    = "West Africa\narchetype",
                              "DRC (Middle, PlusPlus)" = "DRC\narchetype")) +
  labs(x = NULL, y = sprintf("Doses per death averted\n(1 course = %d doses)",
                             DOSES_PER_COURSE),
       title = "Policy B") +
  theme_fig() +
  theme(legend.position = "none",
        plot.title = element_text(size = 11, hjust = 0.5))

fig5_panel_a_AB <- (panel_a_A | panel_a_B) +
  plot_annotation(tag_levels = list(c("a ", "b ")))

save_fig("figure_5_panel-a_dose-efficiency-boxplot_A-vs-B", fig5_panel_a_AB, 8, 5)

message("Figure 5 panel a A-vs-B split saved")

