# ============================================================================
# west_africa_checking.R   (DIAGNOSTIC -- not part of the main pipeline output)
# ----------------------------------------------------------------------------
# PURPOSE
#   Fit the West Africa data with Model A TWICE -- once WITH the targeted tweak
#   priors (and the WA_UFC_01 down-weight), once WITHOUT any of them -- and
#   overlay the two fitted curves parameter-by-parameter. The gap between the two
#   curves for a parameter is exactly "what the tweaks are doing", which is the
#   thing to look at and judge. Nothing here feeds 03 and nothing is saved; the
#   overlay figure is printed to the active graphics device only.
#
# Inputs : data-processed/WestAfrica_QCurve_PreppedData.rds   (uses $anchors)
# Stan   : stan-models/modelA_partialpool_estimateQ_withTweaks.stan
#          stan-models/modelA_partialpool_estimateQ_noTweaks.stan
# Outputs: none (the overlay plot is displayed, not saved)
# ============================================================================

suppressPackageStartupMessages({
  library(cmdstanr)
  library(dplyr); library(tidyr); library(stringr)
  library(readr); library(tibble); library(ggplot2)
})

source(here::here("analyses", "01_latent_response_parameter_estimation", "helpers.R"))
set.seed(123)

# The cleaned West Africa anchors (same input as the production fit in 01).
wa_anchors <- readRDS(file.path(DIR_PROCESSED, "WestAfrica_QCurve_PreppedData.rds"))$anchors

# ---- Shared per-parameter metadata (identical to script 01) ----------------
# Hard admissible domains, then collapse to one row per parameter and add the
# endpoint support/priors. (See 01 section 2 for the param_meta column meanings.)
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
# These mirror 01 exactly so the "with tweaks" fit here matches the production one.
tweak_param <- "p_unsafe_funeral_comm"
j_tweak <- param_meta$param_id[param_meta$parameter == tweak_param]   # row to switch on
anchor_to_downweight <- "WA_UFC_01"
anchor_downweight_sd_mult <- 1.75

# Hyperpriors for the partial pooling and observation scale (see 01 section 5).
hyper <- list(
  mu_t50_raw_prior_mean = 0.0,  mu_t50_raw_prior_sd = 1.25, sigma_t50_prior_sd = 0.50,
  mu_log_k_prior_mean = log(8), mu_log_k_prior_sd = 0.80,   sigma_log_k_prior_sd = 0.40,
  sigma_frac_prior_meanlog = log(0.12), sigma_frac_prior_sdlog = 0.60
)

# The with-tweaks model has the tweak-prior hooks; the no-tweaks model omits them
# entirely. Compile both once.
mod_with <- cmdstan_model(file.path(DIR_STAN, "modelA_partialpool_estimateQ_withTweaks.stan"))
mod_no   <- cmdstan_model(file.path(DIR_STAN, "modelA_partialpool_estimateQ_noTweaks.stan"))

# ---- Fit Model A with or without the tweaks --------------------------------
# Returns a tidy table of the fitted parameter curves tagged by `variant`.
fit_wa <- function(tweaks_on) {

  # Per-observation fit table (param_id, tau, y_obs, obs_sd_mult) -- as in 01.
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

  # Data common to both models (everything except the tweak hooks).
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
    # Build the tweak-prior vectors. `[<-`(x, i, v) is base R for "return x with
    # element i replaced by v" -- here: all-OFF vectors with the targeted
    # parameter (j_tweak) switched on to the same values used in 01.
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
    stan_data <- base_data        # no tweak hooks at all
    mod <- mod_no
  }

  fit <- mod$sample(data = stan_data, seed = 123, chains = 4, parallel_chains = 4,
                    iter_warmup = 1500, iter_sampling = 1500,
                    adapt_delta = 0.97, max_treedepth = 13, refresh = 0)

  # Tidy the fitted parameter curves (theta_pred[j,m]) and tag the variant.
  fit$summary(variables = "theta_pred") %>%
    mutate(param_id = as.integer(str_match(variable, "\\[([0-9]+),([0-9]+)\\]")[, 2]),
           grid_id  = as.integer(str_match(variable, "\\[([0-9]+),([0-9]+)\\]")[, 3]),
           tau = tau_pred[grid_id], relative_day = tau * max_day,
           variant = if_else(tweaks_on, "with tweaks", "no tweaks")) %>%
    left_join(select(param_meta, param_id, parameter), by = "param_id") %>%
    select(variant, parameter, relative_day, mean, q5, q95)
}

# Fit both ways and stack the results; add a human-readable panel label.
curves <- bind_rows(fit_wa(TRUE), fit_wa(FALSE)) %>%
  mutate(panel = factor(PANEL_LOOKUP[parameter], levels = unname(PANEL_LOOKUP)))

# ---- Overlay plot: with-tweaks vs no-tweaks, per parameter -----------------
# One facet per parameter; the literature anchors are shown as grey points so the
# two fits can be judged against the data they were fit to.
anchor_points <- wa_anchors %>%
  mutate(panel = factor(PANEL_LOOKUP[parameter], levels = unname(PANEL_LOOKUP)))

p <- ggplot(curves, aes(relative_day, mean, colour = variant, fill = variant)) +
  geom_ribbon(aes(ymin = q5, ymax = q95), alpha = 0.15, colour = NA) +   # 90% intervals
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
p
