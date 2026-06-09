# ============================================================================
# west_africa_checking.R
# ----------------------------------------------------------------------------
# Diagnostic script (NOT part of the main pipeline output).
#
# Fits the West Africa data with Model A TWICE -- once WITH the targeted tweak
# priors (and the WA_UFC_01 down-weight) and once WITHOUT any of them -- and
# overlays the two fitted curves parameter-by-parameter. The gap between the two
# curves for a given parameter is exactly "what the tweaks are doing", which is
# the thing we want to see and judge.
#
# Inputs : data-processed/wa_anchors.csv
# Stan   : stan-models/modelA_partialpool_estimateQ_withTweaks.stan
#          stan-models/modelA_partialpool_estimateQ_noTweaks.stan
# Outputs: data-processed/wa_checking_curves.csv
#          data-processed/wa_checking_overlay.png
# ============================================================================

suppressPackageStartupMessages({
  library(cmdstanr)
  library(dplyr); library(tidyr); library(stringr)
  library(readr); library(tibble); library(ggplot2)
})

source("helpers.R")
set.seed(123)

wa_anchors <- read_csv("data-processed/wa_anchors.csv", show_col_types = FALSE)

# ---- Shared per-parameter metadata (identical to script 01) ----------------
domain_meta <- tribble(
  ~parameter,                ~abs_min, ~abs_max,
  "delay_hosp",              0.5,      10.0,
  "p_hosp",                  0.0,      1.0,
  "p_ETU",                   0.0,      1.0,
  "latent_IPC",              0.0,      1.0,
  "p_unsafe_funeral_comm",   0.0,      1.0,
  "p_unsafe_funeral_hosp",   0.0,      0.10
)

param_meta <- wa_anchors %>%
  group_by(parameter) %>%
  summarise(lb_prior_mean = first(lower_bound), ub_prior_mean = first(upper_bound),
            direction = first(direction), .groups = "drop") %>%
  left_join(domain_meta, by = "parameter") %>%
  mutate(param_id = match(parameter, PARAM_LEVELS),
         increases = if_else(direction == "up", 1L, 0L)) %>%
  arrange(param_id) %>%
  build_param_support()

J <- nrow(param_meta)
max_day  <- max(wa_anchors$relative_day, na.rm = TRUE)
tau_pred <- seq(0, 1, length.out = 100)

# ---- Tweak configuration (same values as script 01) ------------------------
tweak_param <- "p_unsafe_funeral_comm"
j_tweak <- param_meta$param_id[param_meta$parameter == tweak_param]
anchor_to_downweight <- "WA_UFC_01"
anchor_downweight_sd_mult <- 1.75

hyper <- list(
  mu_t50_raw_prior_mean = 0.0,  mu_t50_raw_prior_sd = 1.25, sigma_t50_prior_sd = 0.50,
  mu_log_k_prior_mean = log(8), mu_log_k_prior_sd = 0.80,   sigma_log_k_prior_sd = 0.40,
  sigma_frac_prior_meanlog = log(0.12), sigma_frac_prior_sdlog = 0.60
)

mod_with <- cmdstan_model("stan-models/modelA_partialpool_estimateQ_withTweaks.stan")
mod_no   <- cmdstan_model("stan-models/modelA_partialpool_estimateQ_noTweaks.stan")

# ---- Fit Model A with or without the tweaks --------------------------------
fit_wa <- function(tweaks_on) {

  fit_df <- wa_anchors %>%
    left_join(select(param_meta, parameter, param_id, increases), by = "parameter") %>%
    mutate(tau = relative_day / max_day, y_obs = value_used,
           weight = if_else(is.na(weight), 1, weight),
           obs_sd_mult = 1 / pmax(weight, 0.25)) %>%
    arrange(param_id, tau)

  # The down-weight is part of the tweak package, so it only applies WITH tweaks.
  if (tweaks_on) {
    fit_df <- fit_df %>%
      mutate(obs_sd_mult = if_else(
        anchor_id == anchor_to_downweight & parameter == tweak_param,
        obs_sd_mult * anchor_downweight_sd_mult, obs_sd_mult))
  }

  base_data <- c(list(
    N = nrow(fit_df), J = J,
    param_id = fit_df$param_id, tau = fit_df$tau,
    y_obs = fit_df$y_obs, obs_sd_mult = fit_df$obs_sd_mult,
    abs_min = param_meta$abs_min, abs_max = param_meta$abs_max,
    lower_floor = param_meta$lower_floor, lower_cap = param_meta$lower_cap,
    upper_cap = param_meta$upper_cap,
    lb_prior_mean = param_meta$lb_prior_mean, lb_prior_sd = param_meta$lb_prior_sd,
    ub_prior_mean = param_meta$ub_prior_mean, ub_prior_sd = param_meta$ub_prior_sd,
    increases = param_meta$increases,
    N_pred = length(tau_pred), tau_pred = tau_pred
  ), hyper)

  if (tweaks_on) {
    zero <- rep(0L, J)
    tweak_data <- list(
      use_t50_tweak = `[<-`(zero, j_tweak, 1L),
      t50_tweak_mean = `[<-`(rep(0.5, J), j_tweak, 0.72),
      t50_tweak_sd   = `[<-`(rep(1.0, J), j_tweak, 0.07),
      use_logk_tweak = `[<-`(zero, j_tweak, 1L),
      logk_tweak_mean = `[<-`(rep(log(8), J), j_tweak, log(11)),
      logk_tweak_sd   = `[<-`(rep(1.0, J), j_tweak, 0.35),
      use_upper_tweak = `[<-`(zero, j_tweak, 1L),
      upper_tweak_mean = `[<-`(param_meta$ub_prior_mean, j_tweak, 0.95),
      upper_tweak_sd   = `[<-`(rep(1.0, J), j_tweak, 0.04),
      use_lower_tweak = zero,
      lower_tweak_mean = param_meta$lb_prior_mean,
      lower_tweak_sd   = rep(1.0, J)
    )
    stan_data <- c(base_data, tweak_data)
    mod <- mod_with
  } else {
    stan_data <- base_data
    mod <- mod_no
  }

  fit <- mod$sample(data = stan_data, seed = 123, chains = 4, parallel_chains = 4,
                    iter_warmup = 1500, iter_sampling = 1500,
                    adapt_delta = 0.97, max_treedepth = 13, refresh = 0)

  fit$summary(variables = "theta_pred") %>%
    mutate(param_id = as.integer(str_match(variable, "\\[([0-9]+),([0-9]+)\\]")[, 2]),
           grid_id  = as.integer(str_match(variable, "\\[([0-9]+),([0-9]+)\\]")[, 3]),
           tau = tau_pred[grid_id], relative_day = tau * max_day,
           variant = if_else(tweaks_on, "with tweaks", "no tweaks")) %>%
    left_join(select(param_meta, param_id, parameter), by = "param_id") %>%
    select(variant, parameter, relative_day, mean, q5, q95)
}

curves <- bind_rows(fit_wa(TRUE), fit_wa(FALSE)) %>%
  mutate(panel = factor(PANEL_LOOKUP[parameter], levels = unname(PANEL_LOOKUP)))

write_csv(curves, "data-processed/wa_checking_curves.csv")

# ---- Overlay plot: with-tweaks vs no-tweaks, per parameter -----------------
anchor_points <- wa_anchors %>%
  mutate(panel = factor(PANEL_LOOKUP[parameter], levels = unname(PANEL_LOOKUP)))

p <- ggplot(curves, aes(relative_day, mean, colour = variant, fill = variant)) +
  geom_ribbon(aes(ymin = q5, ymax = q95), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) +
  geom_point(data = anchor_points, aes(relative_day, value_used),
             inherit.aes = FALSE, colour = "grey25", size = 1.8) +
  facet_wrap(~ panel, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = c("with tweaks" = "#d62728", "no tweaks" = "#1f77b4")) +
  scale_fill_manual(values   = c("with tweaks" = "#d62728", "no tweaks" = "#1f77b4")) +
  labs(title = "West Africa Model A: with vs without targeted tweaks",
       subtitle = "Grey points are the literature anchors; ribbons are 90% intervals",
       x = "Relative outbreak day", y = NULL, colour = NULL, fill = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold"))

ggsave("data-processed/wa_checking_overlay.png", p, width = 12, height = 9, dpi = 150)
message("west_africa_checking.R complete. Wrote wa_checking_curves.csv and wa_checking_overlay.png")
