# =============================================================================
# 02_plot_figure5.R
# Dose efficiency: deaths averted per dose (DDA) under Policy A vs Policy B
#
# Policy A: give OBV to ALL HCW (back-calculated pre-PPE/ETU pool) at 80%
# Policy B: give OBV only to hospital-acquired HCW at varying efficacy (0-80%)
#
# Panel a: DDA boxplot by scenario, Policy A vs B (x-axis = scenario)
# Panel b: DDA curve — Policy A horizontal dashed, Policy B curve vs efficacy
# =============================================================================
source(here::here("analyses", "03_figure_template", "helper_functions_figure_1to4.R"))
OUT_DIR <- here("figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OBV_EFFICACY_A <- 0.80

message("Loading baseline results...")
results <- load_results()

# =============================================================================
# Per-run metrics
# =============================================================================
run_df <- do.call(rbind, lapply(results, function(x) {
  cases  <- x$tdf
  is_hcw <- !is.na(cases$class) & cases$class == "HCW"
  hcw    <- cases[is_hcw, ]
  
  n_all_hcw  <- nrow(hcw)
  died_all   <- !is.na(hcw$outcome) & hcw$outcome
  n_died_all <- sum(died_all)
  
  hosp_acq    <- !is.na(hcw$infection_location) & hcw$infection_location == "hospital"
  n_hosp      <- sum(hosp_acq)
  died_hosp   <- hosp_acq & died_all
  n_died_hosp <- sum(died_hosp)
  
  data.frame(
    scenario     = x$scenario,
    particle_id  = x$particle_id,
    rep          = x$rep,
    ppe_efficacy = if (!is.null(x$ppe_efficacy)) x$ppe_efficacy else NA_real_,
    etu_efficacy = if (!is.null(x$etu_efficacy)) x$etu_efficacy else NA_real_,
    n_all_hcw    = n_all_hcw,
    n_died_all   = n_died_all,
    n_hosp       = n_hosp,
    n_died_hosp  = n_died_hosp,
    stringsAsFactors = FALSE
  )
}))

# Aggregate over reps per particle
particle_df <- run_df %>%
  group_by(scenario, particle_id) %>%
  summarise(
    n_all_hcw    = mean(n_all_hcw),
    n_died_all   = mean(n_died_all),
    n_hosp       = mean(n_hosp),
    n_died_hosp  = mean(n_died_hosp),
    ppe_efficacy = first(ppe_efficacy),
    etu_efficacy = first(etu_efficacy),
    .groups = "drop"
  ) %>%
  mutate(
    # Base pool: back-calculate using PPE only (ETU affects transmission, not infection)
    base_pool = ifelse(
      (1 - ppe_efficacy) > 0,
      n_hosp / (1 - ppe_efficacy),
      n_hosp
    ),
    # averted same for both policies: only hospital-acquired HCW deaths can be averted
    averted_fixed = n_died_hosp * OBV_EFFICACY_A,
    # Policy A DDA: doses = base pool, averted = hospital-acquired deaths
    dda_A = ifelse(averted_fixed > 0, base_pool / averted_fixed, NA_real_),
    scenario_label = factor(SCENARIO_LABELS[scenario], levels = SCENARIO_LABELS)
  )

# =============================================================================
# Sweep Policy B efficacy
# =============================================================================
EFF_SWEEP_GRID <- seq(0, 100, by = 2)

sweep_df <- do.call(rbind, lapply(EFF_SWEEP_GRID, function(pct) {
  eff_b <- OBV_EFFICACY_A * pct / 100
  particle_df %>%
    mutate(
      eff_pct_of_A = pct,
      averted_B    = n_died_hosp * eff_b,
      dda_B        = ifelse(averted_B > 0, n_hosp / averted_B, NA_real_)
    ) %>%
    group_by(scenario, scenario_label, eff_pct_of_A) %>%
    summarise(
      med_dda_A = median(dda_A,  na.rm = TRUE),
      med_dda_B = median(dda_B,  na.rm = TRUE),
      q25_dda_B = quantile(dda_B, 0.25, na.rm = TRUE),
      q75_dda_B = quantile(dda_B, 0.75, na.rm = TRUE),
      .groups = "drop"
    )
}))

policy_a_lines <- particle_df %>%
  group_by(scenario, scenario_label) %>%
  summarise(med_dda_A = median(dda_A, na.rm = TRUE), .groups = "drop")

# =============================================================================
# Color scheme
# =============================================================================
sc_fills <- c(
  setNames(c("#fdd0a2", unname(SCENARIO_COLORS["WestAfrica"])),
           paste(SCENARIO_LABELS["WestAfrica"], c("A", "B"), sep = ".")),
  setNames(c("#a8ddb5", unname(SCENARIO_COLORS["DRC"])),
           paste(SCENARIO_LABELS["DRC"],        c("A", "B"), sep = "."))
)

# =============================================================================
# Panel a: DDA boxplot — x = scenario, fill = Policy A/B shade
# =============================================================================
panel_a_df <- particle_df %>%
  mutate(
    dda_B = ifelse(n_died_hosp * OBV_EFFICACY_A > 0,
                   n_hosp / (n_died_hosp * OBV_EFFICACY_A), NA_real_)
  ) %>%
  select(scenario, particle_id, scenario_label, A = dda_A, B = dda_B) %>%
  tidyr::pivot_longer(cols = c(A, B),
                      names_to = "policy", values_to = "value") %>%
  mutate(group = paste(as.character(scenario_label), policy, sep = ".")) %>%
  filter(!is.na(value))

panel_a <- ggplot(panel_a_df, aes(x = scenario_label, y = value, fill = group)) +
  geom_boxplot(outlier.size = 0.3, width = 0.6, color = "black",
               linewidth = 0.25, position = position_dodge(0.75)) +
  scale_fill_manual(
    values = sc_fills,
    breaks = c(paste(SCENARIO_LABELS["WestAfrica"], "A", sep = "."),
               paste(SCENARIO_LABELS["WestAfrica"], "B", sep = "."),
               paste(SCENARIO_LABELS["DRC"],        "A", sep = "."),
               paste(SCENARIO_LABELS["DRC"],        "B", sep = ".")),
    labels = c("West Africa — Policy A", "West Africa — Policy B",
               "DRC — Policy A",         "DRC — Policy B"),
    name   = NULL
  ) +
  # scale_y_log10(
  #   # expand = expansion(mult = c(0.05, 0.05)),
  #   limits = c(0.1, 100),
  #   breaks = c(0.1, 1, 10, 100),
  #   labels = c("0.1", "1", "10", "100")
  # ) +
  scale_y_log10(
    limits = c(2, 16),                   
    breaks = c(2, 4, 8, 16),          
    labels = c("2", "4", "8", "16") 
  ) +
  labs(x = NULL, y = "Doses per death averted (log scale)") +
  scale_x_discrete(labels = c(
    "West Africa (Worst)" = "West Africa\narchetype",
    "DRC (Middle, PlusPlus)" = "DRC\narchetype"
  )) +
  theme_fig() +
  theme(legend.position  = "bottom",
        legend.key.size  = unit(0.4, "cm"),
        legend.text      = element_text(size = 8)) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

# =============================================================================
# Panel b: DDA curve — x = Policy B efficacy % of Policy A
# =============================================================================
panel_b <- ggplot() +
  geom_ribbon(
    data  = sweep_df %>% filter(eff_pct_of_A >= 10),
    aes(x = eff_pct_of_A, ymin = q25_dda_B, ymax = q75_dda_B,
        fill = scenario_label),
    alpha = 0.2
  ) +
  geom_line(
    data = sweep_df %>% filter(eff_pct_of_A >= 10),
    aes(x = eff_pct_of_A, y = med_dda_B, color = scenario_label),
    linewidth = 1.0
  ) +
  geom_hline(
    data     = policy_a_lines,
    aes(yintercept = med_dda_A, color = scenario_label),
    linetype = "dashed", linewidth = 1.0
  ) +
  scale_color_manual(values = setNames(SCENARIO_COLORS, SCENARIO_LABELS), name = NULL) +
  scale_fill_manual( values = setNames(SCENARIO_COLORS, SCENARIO_LABELS), name = NULL) +
  scale_x_continuous(
    limits = c(10, 100),
    breaks = seq(10, 100, by = 10),
    labels = function(x) paste0(x, "%"),
    name   = sprintf("Policy B OBV efficacy (%% of Policy A's %.0f%%)",
                     OBV_EFFICACY_A * 100)
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05)),
    name   = "Doses per death averted"
  ) +
  theme_fig() +
  theme(legend.position  = "none",
        panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
        panel.grid.minor = element_line(color = "grey95", linewidth = 0.2))

# =============================================================================
# Combine and save
# =============================================================================
library(patchwork)

fig5_all <- (panel_a | panel_b) +
  plot_layout(widths = c(2, 3)) 
  # plot_annotation(
  #   tag_levels = list(c("a ", "b ")),
  #   caption    = sprintf(
  #     "Policy A: all HCW (back-calculated pre-PPE pool), OBV efficacy = %.0f%% | Policy B: hospital-acquired HCW only, efficacy varies",
  #     OBV_EFFICACY_A * 100
  #   )
  

ggsave(
  file.path(OUT_DIR, "figure_5.png"),
  fig5_all, width = 12, height = 5, dpi = 400, units = "in"
)

message("Figure 5 saved")
