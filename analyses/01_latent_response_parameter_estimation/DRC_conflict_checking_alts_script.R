# ============================================================================
# DRC_conflict_checking_alts_script.R   (DIAGNOSTIC / EXPLORATION -- not pipeline)
# ----------------------------------------------------------------------------
# PURPOSE
#   Prototype an ALTERNATIVE to Model B for the DRC conflict scenarios: a Stan
#   model whose functional form admits BOTH partial pooling AND the
#   non-monotonic (rise -> dip -> recover) shape that the conflict creates.
#   We fit it to the DRC "conflict" and "conflict++" SDB data and overlay the
#   resulting curves, conflict vs conflict++, one panel per parameter.
#
# THE FUNCTIONAL FORM
#   For response parameter j, the response-quality curve is decomposed as
#
#       Q_j(tau) = S_j(tau) * (1 - depth_j * Pulse(tau))
#
#     * S_j(tau)   : a normalised monotone logistic (the underlying response
#                    improvement). Its midpoint t50_j and steepness k_j are
#                    PARTIALLY POOLED across the six parameters (as in Model A).
#     * Pulse(tau) : an asymmetric Gaussian bump in [0,1] peaking at t_dip, with
#                    separate rise width (w_on) and recovery width (w_off) so a
#                    sharp hit can recover slowly. SHARED across parameters
#                    (one conflict event hits the whole response apparatus).
#     * depth_j    : how hard the conflict knocks parameter j down, in [0,1],
#                    PARTIALLY POOLED across parameters (logit-normal hierarchy).
#
#   The non-monotonicity lives entirely in the shared Pulse; the pooled pieces
#   (the baseline shape and the per-parameter depth) stay well-defined.
#
# WHAT IDENTIFIES THE DIP (and why this is the interesting bit)
#   The conflict dip is only visible in the DENSE safe-and-dignified-burial (SDB)
#   success series. We therefore feed that series in as the dense observations of
#   community unsafe funerals (y = 1 - success); the other five parameters bring
#   only their sparse literature anchors. So the shared Pulse is pinned almost
#   entirely by the one data-dense parameter, and the partial pooling then
#   propagates that conflict signal to the data-poor parameters. The overlay
#   makes the "how much does the SDB series bleed into everything else" question
#   (the SDB-domination concern) directly visible.
#
# Inputs : data-processed/DRC_QCurve_PreppedData.rds (uses $anchors,
#          $conflict_qseries, $conflict_plusplus_qseries)
# Output : none -- two comparison plots are printed to the graphics device only.
#
# NOTE: this is an EXPLORATION. The model is intentionally richer (and slower)
# than the production Models A/B; the Stan code lives inline (compiled via
# write_stan_file) rather than in stan-models/, to keep it out of the pipeline.
# ============================================================================

suppressPackageStartupMessages({
  library(cmdstanr)
  library(dplyr); library(tidyr); library(stringr)
  library(readr); library(tibble); library(ggplot2)
})

source(here::here("analyses", "01_latent_response_parameter_estimation", "helpers.R"))
set.seed(123)

# ----------------------------------------------------------------------------
# 1. The alternative Stan model (inline)
# ----------------------------------------------------------------------------
# Pooled, non-monotonic Q. See the header for the full description. The blocks
# mirror the production models (endpoint transforms + Jacobians, robust
# Student-t likelihood); the new pieces are the shared Pulse and the pooled
# per-parameter dip depth.
stan_code <- '
functions {
  // Normalised monotone logistic: S(0)=0, S(1)=1 over the response window.
  real s_norm(real tau, real t50, real k) {
    real q0 = inv_logit(k * (0 - t50));
    real q1 = inv_logit(k * (1 - t50));
    return (inv_logit(k * (tau - t50)) - q0) / (q1 - q0);
  }
  // Asymmetric Gaussian pulse in [0,1], peaking at 1 at t_dip. Rise width w_on
  // (before the peak) and recovery width w_off (after), so the conflict can hit
  // sharply and recover slowly when w_off > w_on.
  real pulse(real tau, real t_dip, real w_on, real w_off) {
    real w = tau < t_dip ? w_on : w_off;
    return exp(-0.5 * square((tau - t_dip) / w));
  }
}
data {
  int<lower=1> N;                              // observations (anchors + dense SDB)
  int<lower=1> J;                              // parameters (6)
  array[N] int<lower=1, upper=J> param_id;
  vector<lower=0, upper=1>[N] tau;
  vector[N] y_obs;
  vector<lower=0>[N] obs_sd_mult;
  array[J] int<lower=0, upper=1> increases;

  // Per-parameter endpoint support / priors (identical construction to A and B).
  vector[J] lower_floor; vector[J] lower_cap; vector[J] upper_cap;
  vector[J] lb_prior_mean; vector[J] lb_prior_sd;
  vector[J] ub_prior_mean; vector[J] ub_prior_sd;

  // Hyperpriors: pooled baseline shape.
  real mu_t50_raw_prior_mean; real<lower=0> mu_t50_raw_prior_sd; real<lower=0> sigma_t50_prior_sd;
  real mu_log_k_prior_mean;   real<lower=0> mu_log_k_prior_sd;   real<lower=0> sigma_log_k_prior_sd;
  // Hyperpriors: pooled conflict depth (on the logit scale).
  real mu_depth_prior_mean;   real<lower=0> mu_depth_prior_sd;   real<lower=0> sigma_depth_prior_sd;
  // Priors: shared pulse timing / widths.
  real t_dip_prior_mean;      real<lower=0> t_dip_prior_sd;
  real log_w_on_prior_mean;   real<lower=0> log_w_on_prior_sd;
  real log_w_off_prior_mean;  real<lower=0> log_w_off_prior_sd;
  // Prior: observation noise fraction.
  real sigma_frac_prior_meanlog; real<lower=0> sigma_frac_prior_sdlog;

  int<lower=1> N_pred; vector<lower=0, upper=1>[N_pred] tau_pred;
}
parameters {
  // Pooled baseline shape (non-centred).
  real mu_t50_raw; real<lower=0> sigma_t50; vector[J] t50_std;
  real mu_log_k;   real<lower=0> sigma_log_k; vector[J] log_k_std;
  // Pooled conflict depth (non-centred, logit scale).
  real mu_depth_raw; real<lower=0> sigma_depth; vector[J] depth_std;
  // Shared pulse.
  real<lower=0, upper=1> t_dip;
  real log_w_on; real log_w_off;
  // Per-parameter endpoints + observation scale.
  vector[J] lower_raw; vector[J] upper_gap_raw;
  vector<lower=0>[J] sigma_frac;
}
transformed parameters {
  vector<lower=0, upper=1>[J] t50_param;
  vector<lower=0>[J] k_param;
  vector<lower=0, upper=1>[J] depth_param;
  vector[J] lower_est; vector[J] upper_est; vector[J] span_est;
  real w_on  = exp(log_w_on);
  real w_off = exp(log_w_off);

  for (j in 1:J) {
    real s_l; real s_u; real lower_j; real upper_j;
    t50_param[j]   = inv_logit(mu_t50_raw + sigma_t50 * t50_std[j]);
    k_param[j]     = exp(mu_log_k + sigma_log_k * log_k_std[j]);
    depth_param[j] = inv_logit(mu_depth_raw + sigma_depth * depth_std[j]);

    s_l = inv_logit(lower_raw[j]);
    lower_j = lower_floor[j] + (lower_cap[j] - lower_floor[j]) * s_l;
    s_u = inv_logit(upper_gap_raw[j]);
    upper_j = lower_j + (upper_cap[j] - lower_j) * s_u;
    lower_est[j] = lower_j; upper_est[j] = upper_j; span_est[j] = upper_j - lower_j;
  }
}
model {
  // Pooled baseline shape.
  mu_t50_raw ~ normal(mu_t50_raw_prior_mean, mu_t50_raw_prior_sd);
  sigma_t50  ~ normal(0, sigma_t50_prior_sd);
  t50_std    ~ normal(0, 1);
  mu_log_k    ~ normal(mu_log_k_prior_mean, mu_log_k_prior_sd);
  sigma_log_k ~ normal(0, sigma_log_k_prior_sd);
  log_k_std   ~ normal(0, 1);
  // Pooled conflict depth.
  mu_depth_raw ~ normal(mu_depth_prior_mean, mu_depth_prior_sd);
  sigma_depth  ~ normal(0, sigma_depth_prior_sd);
  depth_std    ~ normal(0, 1);
  // Shared pulse.
  t_dip     ~ normal(t_dip_prior_mean, t_dip_prior_sd);
  log_w_on  ~ normal(log_w_on_prior_mean, log_w_on_prior_sd);
  log_w_off ~ normal(log_w_off_prior_mean, log_w_off_prior_sd);
  // Observation scale.
  sigma_frac ~ lognormal(sigma_frac_prior_meanlog, sigma_frac_prior_sdlog);

  // Endpoint priors with the usual manual Jacobian corrections.
  for (j in 1:J) {
    real s_l = inv_logit(lower_raw[j]);
    real s_u = inv_logit(upper_gap_raw[j]);
    target += normal_lpdf(lower_est[j] | lb_prior_mean[j], lb_prior_sd[j]);
    target += log(lower_cap[j] - lower_floor[j]) + log(s_l) + log1m(s_l);
    target += normal_lpdf(upper_est[j] | ub_prior_mean[j], ub_prior_sd[j]);
    target += log(upper_cap[j] - lower_est[j]) + log(s_u) + log1m(s_u);
  }

  // Observation model.
  for (n in 1:N) {
    int j = param_id[n];
    real Q = s_norm(tau[n], t50_param[j], k_param[j])
             * (1 - depth_param[j] * pulse(tau[n], t_dip, w_on, w_off));
    real mu = increases[j] == 1 ? lower_est[j] + span_est[j] * Q
                                : upper_est[j] - span_est[j] * Q;
    y_obs[n] ~ student_t(4, mu, sigma_frac[j] * span_est[j] * obs_sd_mult[n]);
  }
}
generated quantities {
  matrix[J, N_pred] Q_pred;       // the non-monotonic response-quality curve
  matrix[J, N_pred] theta_pred;   // the parameter value (native units)
  for (j in 1:J) {
    for (m in 1:N_pred) {
      real Q = s_norm(tau_pred[m], t50_param[j], k_param[j])
               * (1 - depth_param[j] * pulse(tau_pred[m], t_dip, w_on, w_off));
      Q_pred[j, m] = Q;
      theta_pred[j, m] = increases[j] == 1 ? lower_est[j] + span_est[j] * Q
                                           : upper_est[j] - span_est[j] * Q;
    }
  }
}
'
mod <- cmdstan_model(write_stan_file(stan_code))

# ----------------------------------------------------------------------------
# 2. Inputs and shared metadata
# ----------------------------------------------------------------------------
drc_prep    <- readRDS(file.path(DIR_PROCESSED, "DRC_QCurve_PreppedData.rds"))
drc_anchors <- drc_prep$anchors

# DRC admissible domains (same as Model B in script 02).
domain_meta <- tribble(
  ~parameter,                ~abs_min, ~abs_max,
  "delay_hosp",              0.5,      12.0,
  "p_hosp",                  0.0,      1.0,
  "p_ETU",                   0.0,      1.0,
  "latent_IPC",              0.0,      1.0,
  "p_unsafe_funeral_comm",   0.0,      1.0,
  "p_unsafe_funeral_hosp",   0.0,      0.010
)

# One row per parameter, with endpoint support/priors (see 01 for the columns).
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
stopifnot(J == length(PARAM_LEVELS))

tau_pred <- seq(0, 1, length.out = 200)   # fine grid so the dip is well resolved

# Priors passed to Stan (kept here, explicit, so they are easy to inspect/tune).
hyper <- list(
  mu_t50_raw_prior_mean = 0.0,    mu_t50_raw_prior_sd = 1.25, sigma_t50_prior_sd  = 0.50,
  mu_log_k_prior_mean   = log(8), mu_log_k_prior_sd   = 0.80, sigma_log_k_prior_sd = 0.40,
  mu_depth_prior_mean   = 0.0,    mu_depth_prior_sd   = 1.50, sigma_depth_prior_sd = 1.00,
  t_dip_prior_mean      = 0.35,   t_dip_prior_sd      = 0.25,
  log_w_on_prior_mean   = log(0.08), log_w_on_prior_sd  = 0.70,
  log_w_off_prior_mean  = log(0.20), log_w_off_prior_sd = 0.70,
  sigma_frac_prior_meanlog = log(0.12), sigma_frac_prior_sdlog = 0.60
)

# ----------------------------------------------------------------------------
# 3. Fit one DRC scenario to the alternative model
# ----------------------------------------------------------------------------
# Assemble observations: sparse literature anchors for every parameter, PLUS the
# dense SDB community-burial series (1 - success) as community-unsafe-funeral
# observations. The SDB points are weighted by sample size (more eligible burials
# -> tighter), rescaled so the median SDB point carries an anchor-like weight.
fit_alt <- function(qseries, label) {

  message("\n==== Fitting alternative (pooled, non-monotonic) model: ", label, " ====")
  max_day <- max(c(drc_anchors$relative_day, qseries$relative_day), na.rm = TRUE)

  anchor_obs <- drc_anchors %>%
    transmute(parameter, relative_day, y_obs = value_used, src = "anchor",
              obs_sd_mult = 1 / pmax(if_else(is.na(weight), 1, weight), 0.25))

  sdb_rows <- qseries %>% filter(n_eligible_sum > 0)
  med_n <- median(sdb_rows$n_eligible_sum)
  sdb_obs <- sdb_rows %>%
    transmute(parameter = "p_unsafe_funeral_comm", relative_day,
              y_obs = unsafe_funeral_comm_proxy, src = "sdb",
              obs_sd_mult = sqrt(med_n / pmax(n_eligible_sum, 1)))

  obs <- bind_rows(anchor_obs, sdb_obs) %>%
    left_join(select(param_meta, parameter, param_id, increases), by = "parameter") %>%
    mutate(tau = relative_day / max_day) %>%
    arrange(param_id, tau)

  stan_data <- c(list(
    N = nrow(obs), J = J,
    param_id = obs$param_id, tau = obs$tau,
    y_obs = obs$y_obs, obs_sd_mult = obs$obs_sd_mult,
    increases = param_meta$increases,
    lower_floor = param_meta$lower_floor, lower_cap = param_meta$lower_cap,
    upper_cap = param_meta$upper_cap,
    lb_prior_mean = param_meta$lb_prior_mean, lb_prior_sd = param_meta$lb_prior_sd,
    ub_prior_mean = param_meta$ub_prior_mean, ub_prior_sd = param_meta$ub_prior_sd,
    N_pred = length(tau_pred), tau_pred = tau_pred
  ), hyper)

  fit <- mod$sample(
    data = stan_data, seed = 123,
    chains = 4, parallel_chains = 4,
    iter_warmup = 1000, iter_sampling = 1000,
    adapt_delta = 0.95, max_treedepth = 12, refresh = 250
  )
  cat("\nDiagnostics (", label, "):\n", sep = ""); print(fit$diagnostic_summary())
  # Pooling + conflict-shape posterior summary (the quantities of interest).
  cat("\nShape/conflict posterior (", label, "):\n", sep = "")
  print(fit$summary(variables = c("sigma_t50", "sigma_log_k", "sigma_depth",
                                  "t_dip", "w_on", "w_off", "depth_param"))[,
        c("variable", "mean", "q5", "q95")])

  # Tidy a [J, N_pred] matrix variable into a long per-parameter curve.
  tidy_mat <- function(varname) {
    fit$summary(variables = varname) %>%
      mutate(param_id = as.integer(str_match(variable, "\\[([0-9]+),([0-9]+)\\]")[, 2]),
             grid_id  = as.integer(str_match(variable, "\\[([0-9]+),([0-9]+)\\]")[, 3]),
             tau = tau_pred[grid_id], relative_day = tau * max_day, scenario = label) %>%
      left_join(select(param_meta, param_id, parameter), by = "param_id") %>%
      select(scenario, parameter, relative_day, mean, q5, q95)
  }

  list(theta = tidy_mat("theta_pred"),
       Q     = tidy_mat("Q_pred"),
       obs   = obs %>% mutate(scenario = label))
}

# ----------------------------------------------------------------------------
# 4. Fit conflict and conflict++ and overlay them per parameter
# ----------------------------------------------------------------------------
fit_conflict <- fit_alt(drc_prep$conflict_qseries,          "conflict")
fit_pluspl   <- fit_alt(drc_prep$conflict_plusplus_qseries, "conflict++")

panel_lvls <- unname(PANEL_LOOKUP)
add_panel  <- function(df) mutate(df, panel = factor(PANEL_LOOKUP[parameter], levels = panel_lvls))

theta_all <- bind_rows(fit_conflict$theta, fit_pluspl$theta) %>% add_panel()
Q_all     <- bind_rows(fit_conflict$Q,     fit_pluspl$Q)     %>% add_panel()

# Data points for the overlay: literature anchors (grey, shared) and the dense
# SDB community points (coloured by scenario, since the ++ collapse changes them).
anchor_pts <- fit_conflict$obs %>% filter(src == "anchor") %>% add_panel()
sdb_pts    <- bind_rows(fit_conflict$obs, fit_pluspl$obs) %>% filter(src == "sdb") %>% add_panel()

scen_cols <- c("conflict" = "#1f77b4", "conflict++" = "#d62728")

# --- Plot 1: the parameter curves (native units), conflict vs conflict++ ---
p_theta <- ggplot(theta_all, aes(relative_day, mean, colour = scenario, fill = scenario)) +
  geom_ribbon(aes(ymin = q5, ymax = q95), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.9) +
  geom_point(data = anchor_pts, aes(relative_day, y_obs),
             inherit.aes = FALSE, colour = "grey30", size = 1.8) +
  geom_point(data = sdb_pts, aes(relative_day, y_obs, colour = scenario),
             inherit.aes = FALSE, size = 0.8, alpha = 0.5) +
  facet_wrap(~ panel, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = scen_cols) + scale_fill_manual(values = scen_cols) +
  labs(title = "DRC alternative model: pooled + non-monotonic fitted parameter curves",
       subtitle = "conflict vs conflict++ overlaid; grey = literature anchors, faint points = dense SDB community data",
       x = "Relative outbreak day", y = NULL, colour = NULL, fill = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom", strip.text = element_text(face = "bold"))

# --- Plot 2: the latent response-quality Q_j(tau), where the dip is clearest ---
p_Q <- ggplot(Q_all, aes(relative_day, mean, colour = scenario, fill = scenario)) +
  geom_ribbon(aes(ymin = q5, ymax = q95), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ panel, scales = "fixed", ncol = 2) +
  scale_colour_manual(values = scen_cols) + scale_fill_manual(values = scen_cols) +
  labs(title = "DRC alternative model: latent response-quality Q_j(tau)",
       subtitle = "The shared conflict 'pulse' produces the rise-dip-recover shape; depth is pooled across parameters",
       x = "Relative outbreak day", y = "Q (0 = worst, 1 = best)", colour = NULL, fill = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom", strip.text = element_text(face = "bold"))

print(p_theta)   # display only; not saved
print(p_Q)       # display only; not saved

message("\nDRC_conflict_checking_alts_script.R complete (plots printed, nothing saved).")
