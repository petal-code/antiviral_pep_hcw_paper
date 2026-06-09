# ============================================================================
# DRC_no_conflict_checking.R
# ----------------------------------------------------------------------------
# Diagnostic script (NOT part of the main pipeline output).
#
# Explores how to model the DRC NO-CONFLICT scenario, by fitting it three ways
# and overlaying the per-parameter curves:
#
#   (1) MODEL B  -- the original approach. The empirical SDB no-conflict success
#       curve is supplied as the FIXED shared Q; only the endpoints are estimated.
#
#   (2) MODEL A, SDB at full weight -- estimate the Q shape (partially pooled
#       across parameters) using ALL the data: the literature anchors for every
#       parameter PLUS the SDB no-conflict points as observations of community
#       unsafe funerals. Because the SDB points are far more numerous than the
#       sparse anchors, they dominate the community-unsafe-funeral fit and, via
#       the partial pooling, pull the shared shape toward the SDB curve.
#
#   (3) MODEL A, SDB down-weighted -- the same, but the SDB points are
#       collectively down-weighted to about SDB_TARGET_EQUIV_ANCHORS anchors'
#       worth of information, so they inform the shape without numerically
#       swamping the literature anchors.
#
# The point of (2) vs (3) is to SEE how much the dense SDB series dominates the
# pooled estimate -- the concern that motivated keeping no-conflict out of the
# main pipeline for now. The script also prints the estimated pooling SDs
# (sigma_t50, sigma_log_k) for the two Model A fits: small values mean strong
# pooling (the SDB shape is imposed on all parameters), large values mean each
# parameter keeps its own shape.
#
# Inputs : data-processed/DRC_QCurve/DRC_QCurve_PreppedData.rds   (uses $anchors and $no_conflict_qseries)
# Stan   : stan-models/modelB_fixedQ_boundsOnly.stan
#          stan-models/modelA_partialpool_estimateQ_noTweaks.stan
# Outputs: data-processed/drc_no_conflict_checking_curves.csv
#          data-processed/drc_no_conflict_checking_overlay.png
# ============================================================================

suppressPackageStartupMessages({
  library(cmdstanr)
  library(dplyr); library(tidyr); library(stringr)
  library(readr); library(tibble); library(ggplot2)
})

source(here::here("analyses", "01_latent_response_parameter_estimation", "helpers.R"))
set.seed(123)

# How many literature anchors' worth of information the whole SDB series should
# carry in the down-weighted Model A fit.
SDB_TARGET_EQUIV_ANCHORS <- 3

# ----------------------------------------------------------------------------
# Data and shared metadata
# ----------------------------------------------------------------------------
# Pull what we need out of the bundled DRC prep object (from 00).
drc_prep    <- readRDS(file.path(DIR_PROCESSED, "DRC_QCurve/DRC_QCurve_PreppedData.rds"))
qseries     <- drc_prep$no_conflict_qseries   # the truncated no-conflict Q series
drc_anchors <- drc_prep$anchors               # cleaned DRC literature anchors

# The no-conflict horizon ends where its (truncated) Q series ends; clip anchors
# to the same window.
max_day     <- max(qseries$relative_day, na.rm = TRUE)
drc_anchors <- drc_anchors %>% filter(relative_day <= max_day)

domain_meta <- tribble(
  ~parameter,                ~abs_min, ~abs_max,
  "delay_hosp",              0.5,      12.0,
  "p_hosp",                  0.0,      1.0,
  "p_ETU",                   0.0,      1.0,
  "latent_IPC",              0.0,      1.0,
  "p_unsafe_funeral_comm",   0.0,      1.0,
  "p_unsafe_funeral_hosp",   0.0,      0.010
)

param_meta <- drc_anchors %>%
  group_by(parameter) %>%
  summarise(lb_prior_mean = first(lower_bound), ub_prior_mean = first(upper_bound),
            direction = first(direction), .groups = "drop") %>%
  left_join(domain_meta, by = "parameter") %>%
  mutate(param_id = match(parameter, PARAM_LEVELS),
         increases = if_else(direction == "up", 1L, 0L)) %>%
  arrange(param_id) %>%
  build_param_support()

J <- nrow(param_meta)
tau_pred <- seq(0, 1, length.out = 100)

hyper <- list(
  mu_t50_raw_prior_mean = 0.0,  mu_t50_raw_prior_sd = 1.25, sigma_t50_prior_sd = 0.50,
  mu_log_k_prior_mean = log(8), mu_log_k_prior_sd = 0.80,   sigma_log_k_prior_sd = 0.40,
  sigma_frac_prior_meanlog = log(0.12), sigma_frac_prior_sdlog = 0.60
)
SIGMA_FRAC_PRIOR_MEANLOG <- log(0.12); SIGMA_FRAC_PRIOR_SDLOG <- 0.60

mod_B <- cmdstan_model(file.path(DIR_STAN, "modelB_fixedQ_boundsOnly.stan"))
mod_A <- cmdstan_model(file.path(DIR_STAN, "modelA_partialpool_estimateQ_noTweaks.stan"))

# Helper: pull a tidy per-parameter theta curve out of a fit.
tidy_theta <- function(fit, variant) {
  fit$summary(variables = "theta_pred") %>%
    mutate(param_id = as.integer(str_match(variable, "\\[([0-9]+),([0-9]+)\\]")[, 2]),
           grid_id  = as.integer(str_match(variable, "\\[([0-9]+),([0-9]+)\\]")[, 3]),
           tau = tau_pred[grid_id], relative_day = tau * max_day, variant = variant) %>%
    left_join(select(param_meta, param_id, parameter), by = "param_id") %>%
    select(variant, parameter, relative_day, mean)
}

# ============================================================================
# (1) MODEL B -- SDB no-conflict curve supplied as the fixed shared Q
# ============================================================================
q_at_day <- make_interp(qseries$relative_day, qseries$q_value)
success_max <- max(qseries$success_smoothed, na.rm = TRUE)

fitB_df <- drc_anchors %>%
  left_join(select(param_meta, parameter, param_id, increases), by = "parameter") %>%
  mutate(q_obs = clip01(q_at_day(relative_day)), y_obs = value_used,
         weight = if_else(is.na(weight), 1, weight), obs_sd_mult = 1 / pmax(weight, 0.25)) %>%
  arrange(param_id)

tau_grid <- sort(unique(c(seq(0, 1, length.out = 250), qseries$tau_q)))
q_pred   <- clip01(make_interp(qseries$tau_q, qseries$q_value)(tau_grid))

# Endpoint tweaks: same community/hospital unsafe-funeral anchoring as script 02.
zero <- rep(0L, J)
twB <- list(use_upper = zero, upper_mean = param_meta$ub_prior_mean, upper_sd = rep(1, J),
            use_lower = zero, lower_mean = param_meta$lb_prior_mean, lower_sd = rep(1, J))
j_ufc <- param_meta$param_id[param_meta$parameter == "p_unsafe_funeral_comm"]
j_ufh <- param_meta$param_id[param_meta$parameter == "p_unsafe_funeral_hosp"]
twB$use_upper[j_ufc] <- 1L; twB$upper_mean[j_ufc] <- 1.00;                twB$upper_sd[j_ufc] <- 0.02
twB$use_lower[j_ufc] <- 1L; twB$lower_mean[j_ufc] <- clip01(1 - success_max); twB$lower_sd[j_ufc] <- 0.02
twB$use_upper[j_ufh] <- 1L; twB$upper_mean[j_ufh] <- 0.010;  twB$upper_sd[j_ufh] <- 0.003
twB$use_lower[j_ufh] <- 1L; twB$lower_mean[j_ufh] <- 0.0005; twB$lower_sd[j_ufh] <- 0.001

dataB <- list(
  N = nrow(fitB_df), J = J, N_pred = length(tau_grid),
  param_id = fitB_df$param_id, q_obs = fitB_df$q_obs, y_obs = fitB_df$y_obs,
  obs_sd_mult = fitB_df$obs_sd_mult,
  abs_min = param_meta$abs_min, abs_max = param_meta$abs_max,
  lower_floor = param_meta$lower_floor, lower_cap = param_meta$lower_cap, upper_cap = param_meta$upper_cap,
  lb_prior_mean = param_meta$lb_prior_mean, lb_prior_sd = param_meta$lb_prior_sd,
  ub_prior_mean = param_meta$ub_prior_mean, ub_prior_sd = param_meta$ub_prior_sd,
  increases = param_meta$increases,
  sigma_frac_prior_meanlog = SIGMA_FRAC_PRIOR_MEANLOG, sigma_frac_prior_sdlog = SIGMA_FRAC_PRIOR_SDLOG,
  use_upper_tweak = twB$use_upper, upper_tweak_mean = twB$upper_mean, upper_tweak_sd = twB$upper_sd,
  use_lower_tweak = twB$use_lower, lower_tweak_mean = twB$lower_mean, lower_tweak_sd = twB$lower_sd,
  q_pred = q_pred
)
fitB <- mod_B$sample(data = dataB, seed = 123, chains = 4, parallel_chains = 4,
                     iter_warmup = 1500, iter_sampling = 1500,
                     adapt_delta = 0.98, max_treedepth = 13, refresh = 0)

# Model B predicts on tau_grid (not tau_pred); tidy it directly and apply the
# absolute community-unsafe-funeral override.
ufc_on_tau <- make_interp(qseries$tau_q, qseries$unsafe_funeral_comm_proxy)
curveB <- fitB$summary(variables = "theta_pred") %>%
  mutate(param_id = as.integer(str_match(variable, "\\[([0-9]+),([0-9]+)\\]")[, 2]),
         grid_id  = as.integer(str_match(variable, "\\[([0-9]+),([0-9]+)\\]")[, 3]),
         tau = tau_grid[grid_id], relative_day = tau * max_day, variant = "Model B (fixed SDB Q)") %>%
  left_join(select(param_meta, param_id, parameter), by = "param_id") %>%
  mutate(mean = if_else(parameter == "p_unsafe_funeral_comm", clip01(ufc_on_tau(tau)), mean)) %>%
  select(variant, parameter, relative_day, mean)

# ============================================================================
# (2)/(3) MODEL A -- estimate Q from anchors + SDB points (full / down-weighted)
# ============================================================================
# SDB no-conflict observations of community unsafe funerals: y = 1 - success, at
# the real data-bearing support points.
sdb_obs <- qseries %>%
  filter(n_eligible_sum > 0) %>%
  transmute(parameter = "p_unsafe_funeral_comm",
            relative_day, value_used = unsafe_funeral_comm_proxy)
n_sdb <- nrow(sdb_obs)

fit_modelA <- function(sdb_downweight) {

  # SDB observation weight: full = 1 each; down-weighted so the whole SDB series
  # carries ~SDB_TARGET_EQUIV_ANCHORS anchors' worth of weight.
  sdb_weight <- if (sdb_downweight) SDB_TARGET_EQUIV_ANCHORS / n_sdb else 1

  anchor_rows <- drc_anchors %>%
    transmute(parameter, relative_day, y_obs = value_used,
              obs_sd_mult = 1 / pmax(if_else(is.na(weight), 1, weight), 0.25))
  sdb_rows <- sdb_obs %>%
    transmute(parameter, relative_day, y_obs = value_used,
              obs_sd_mult = 1 / sdb_weight)   # no 0.25 floor: allow strong down-weighting

  # The Model A fit table stacks the sparse literature anchors (all parameters)
  # with the many SDB observations (community unsafe funerals only). Columns:
  # param_id, tau, y_obs, obs_sd_mult -- same contract as the 01 fit table.
  fitA_df <- bind_rows(anchor_rows, sdb_rows) %>%
    left_join(select(param_meta, parameter, param_id, increases), by = "parameter") %>%
    mutate(tau = relative_day / max_day) %>%
    arrange(param_id, tau)

  dataA <- c(list(
    N = nrow(fitA_df), J = J,
    param_id = fitA_df$param_id, tau = fitA_df$tau,
    y_obs = fitA_df$y_obs, obs_sd_mult = fitA_df$obs_sd_mult,
    abs_min = param_meta$abs_min, abs_max = param_meta$abs_max,
    lower_floor = param_meta$lower_floor, lower_cap = param_meta$lower_cap, upper_cap = param_meta$upper_cap,
    lb_prior_mean = param_meta$lb_prior_mean, lb_prior_sd = param_meta$lb_prior_sd,
    ub_prior_mean = param_meta$ub_prior_mean, ub_prior_sd = param_meta$ub_prior_sd,
    increases = param_meta$increases,
    N_pred = length(tau_pred), tau_pred = tau_pred
  ), hyper)

  fit <- mod_A$sample(data = dataA, seed = 123, chains = 4, parallel_chains = 4,
                      iter_warmup = 1500, iter_sampling = 1500,
                      adapt_delta = 0.97, max_treedepth = 13, refresh = 0)

  variant <- if (sdb_downweight) "Model A (SDB down-weighted)" else "Model A (SDB full weight)"
  # Report the pooling SDs: how strongly the SDB-dominated UFC shape is imposed
  # on the other parameters.
  pooling <- fit$summary(variables = c("sigma_t50", "sigma_log_k"))
  cat("\nPooling strength -- ", variant, ":\n", sep = ""); print(pooling[, c("variable", "mean", "q5", "q95")])

  tidy_theta(fit, variant)
}

curveA_full <- fit_modelA(sdb_downweight = FALSE)
curveA_dw   <- fit_modelA(sdb_downweight = TRUE)

# ============================================================================
# Combine and plot
# ============================================================================
curves <- bind_rows(curveB, curveA_full, curveA_dw) %>%
  mutate(panel = factor(PANEL_LOOKUP[parameter], levels = unname(PANEL_LOOKUP)))
write_csv(curves, file.path(DIR_PROCESSED, "drc_no_conflict_checking_curves.csv"))

# Overlay the SDB community-unsafe-funeral data points for reference.
sdb_points <- sdb_obs %>%
  mutate(parameter = "p_unsafe_funeral_comm",
         panel = factor(PANEL_LOOKUP[parameter], levels = unname(PANEL_LOOKUP)))

p <- ggplot(curves, aes(relative_day, mean, colour = variant)) +
  geom_line(linewidth = 0.9) +
  geom_point(data = sdb_points, aes(relative_day, value_used),
             inherit.aes = FALSE, colour = "grey50", size = 1.2, alpha = 0.7) +
  facet_wrap(~ panel, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = c(
    "Model B (fixed SDB Q)"        = "#000000",
    "Model A (SDB full weight)"    = "#d62728",
    "Model A (SDB down-weighted)"  = "#1f77b4")) +
  labs(title = "DRC no-conflict: Model B vs Model A (SDB full / down-weighted)",
       subtitle = "Grey points: SDB community unsafe-funeral data (1 - success)",
       x = "Relative outbreak day", y = NULL, colour = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold"))

ggsave(file.path(DIR_PROCESSED, "drc_no_conflict_checking_overlay.png"), p, width = 12, height = 9, dpi = 150)
message("DRC_no_conflict_checking.R complete. Wrote curves CSV and overlay PNG.")
