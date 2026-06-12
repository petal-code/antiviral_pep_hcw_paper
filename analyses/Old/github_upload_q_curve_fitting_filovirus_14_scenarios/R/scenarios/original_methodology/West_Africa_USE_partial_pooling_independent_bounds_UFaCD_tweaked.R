# ============================================================
# West Africa worst-case:
# parameter-specific finite-window Q_j(t) with partial pooling
# on shape parameters only
#
# SMALL TARGETED TWEAK ONLY FOR:
#   p_unsafe_funeral_comm
#
# All other parameters remain unchanged from the current model.
# ============================================================

suppressPackageStartupMessages({
  library(cmdstanr)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(tibble)
  library(readr)
})

# ------------------------------------------------------------
# User settings
# ------------------------------------------------------------
workbook_path <- find_input_file(c("filovirus_three_scenario_curve_inputs_bestcase_recreate.xlsx"), "original methodology curve workbook")
sheet_name    <- "Worst_WestAfrica"
stan_file     <- "west_africa_partial_pool_normalised_Q_estimated_bounds_tweaked_ufc.stan"
out_prefix    <- "west_africa_partial_pool_normalisedQ_estimated_bounds_tweaked_ufc"

# ------------------------------------------------------------
# Hard admissible domains for each parameter
# ------------------------------------------------------------
domain_meta <- tibble::tribble(
  ~parameter,                 ~abs_min, ~abs_max,
  "delay_hosp",               0.5,      10.0,
  "p_hosp",                   0.0,      1.0,
  "p_ETU",                    0.0,      1.0,
  "latent_IPC",               0.0,      1.0,
  "p_unsafe_funeral_comm",    0.0,      1.0,
  "p_unsafe_funeral_hosp",    0.0,      0.10
)

# ------------------------------------------------------------
# Prior / support settings for lower and upper bounds
# ------------------------------------------------------------
support_half_width_mult <- 0.75

bound_sd_frac_of_initial_span <- 0.15
bound_sd_frac_of_domain       <- 0.03

# ------------------------------------------------------------
# Priors for parameter-specific observation scale:
# sigma_y_j = sigma_frac_j * (upper_j - lower_j)
# ------------------------------------------------------------
sigma_frac_prior_meanlog <- log(0.12)
sigma_frac_prior_sdlog   <- 0.60

# ------------------------------------------------------------
# Hyperpriors for partial pooling of Q_j shape parameters
# ------------------------------------------------------------
mu_t50_raw_prior_mean <- 0.0
mu_t50_raw_prior_sd   <- 1.25
sigma_t50_prior_sd    <- 0.50

mu_log_k_prior_mean   <- log(8)
mu_log_k_prior_sd     <- 0.80
sigma_log_k_prior_sd  <- 0.40

# ------------------------------------------------------------
# TARGETED TWEAK FOR p_unsafe_funeral_comm
#
# Interpretation:
#   - slightly later decline
#   - slightly steeper decline
#   - slightly higher upper bound
#
# These are extra priors ON TOP OF the existing model.
# Set any use_* flag to FALSE to turn that piece off.
# ------------------------------------------------------------
tweak_param <- "p_unsafe_funeral_comm"

use_t50_tweak   <- TRUE
t50_tweak_mean  <- 0.72     # later midpoint
t50_tweak_sd    <- 0.07

use_logk_tweak  <- TRUE
logk_tweak_mean <- log(11)  # steeper decline than current baseline
logk_tweak_sd   <- 0.35

use_upper_tweak <- TRUE
upper_tweak_mean <- 0.95    # keeps early unsafe funeral risk near ceiling
upper_tweak_sd   <- 0.04

# Optional extra downweight on one discordant anchor
use_anchor_downweight      <- TRUE
anchor_to_downweight       <- "WA_UFC_01"
anchor_downweight_sd_mult  <- 1.75   # >1 means looser fit to that anchor

# ------------------------------------------------------------
# Write Stan model to disk
# ------------------------------------------------------------
stan_code <- '
functions {
  real q_raw(real tau, real t50, real k) {
    return inv_logit(k * (tau - t50));
  }

  real q_norm(real tau, real t50, real k) {
    real q0 = q_raw(0, t50, k);
    real q1 = q_raw(1, t50, k);
    return (q_raw(tau, t50, k) - q0) / (q1 - q0);
  }
}

data {
  int<lower=1> N;
  int<lower=1> J;

  array[N] int<lower=1, upper=J> param_id;

  vector<lower=0, upper=1>[N] tau;
  vector[N] y_obs;
  vector<lower=0>[N] obs_sd_mult;

  vector[J] abs_min;
  vector[J] abs_max;

  vector[J] lower_floor;
  vector[J] lower_cap;
  vector[J] upper_cap;

  vector[J] lb_prior_mean;
  vector[J] lb_prior_sd;
  vector[J] ub_prior_mean;
  vector[J] ub_prior_sd;

  array[J] int<lower=0, upper=1> increases;

  real sigma_frac_prior_meanlog;
  real<lower=0> sigma_frac_prior_sdlog;

  real mu_t50_raw_prior_mean;
  real<lower=0> mu_t50_raw_prior_sd;
  real<lower=0> sigma_t50_prior_sd;

  real mu_log_k_prior_mean;
  real<lower=0> mu_log_k_prior_sd;
  real<lower=0> sigma_log_k_prior_sd;

  // Targeted extra priors for selected parameters
  array[J] int<lower=0, upper=1> use_t50_tweak;
  vector<lower=0, upper=1>[J] t50_tweak_mean;
  vector<lower=0>[J] t50_tweak_sd;

  array[J] int<lower=0, upper=1> use_logk_tweak;
  vector[J] logk_tweak_mean;
  vector<lower=0>[J] logk_tweak_sd;

  array[J] int<lower=0, upper=1> use_upper_tweak;
  vector[J] upper_tweak_mean;
  vector<lower=0>[J] upper_tweak_sd;

  int<lower=1> N_pred;
  vector<lower=0, upper=1>[N_pred] tau_pred;
}

parameters {
  // Partial pooling for parameter-specific Q_j shapes
  real mu_t50_raw;
  real<lower=0> sigma_t50;
  vector[J] t50_std;

  real mu_log_k;
  real<lower=0> sigma_log_k;
  vector[J] log_k_std;

  // Unpooled parameter-specific lower/upper bounds
  vector[J] lower_raw;
  vector[J] upper_gap_raw;

  // Unpooled parameter-specific observation scale
  vector<lower=0>[J] sigma_frac;
}

transformed parameters {
  vector<lower=0, upper=1>[J] t50_param;
  vector<lower=0>[J] k_param;
  vector[J] lower_est;
  vector[J] upper_est;
  vector[J] span_est;

  for (j in 1:J) {
    real s_l;
    real s_u;
    real lower_j;
    real upper_j;

    t50_param[j] = inv_logit(mu_t50_raw + sigma_t50 * t50_std[j]);
    k_param[j]   = exp(mu_log_k + sigma_log_k * log_k_std[j]);

    s_l = inv_logit(lower_raw[j]);
    lower_j = lower_floor[j] + (lower_cap[j] - lower_floor[j]) * s_l;

    s_u = inv_logit(upper_gap_raw[j]);
    upper_j = lower_j + (upper_cap[j] - lower_j) * s_u;

    lower_est[j] = lower_j;
    upper_est[j] = upper_j;
    span_est[j]  = upper_j - lower_j;
  }
}

model {
  // ----------------------------
  // Hyperpriors for partial pooling
  // ----------------------------
  mu_t50_raw ~ normal(mu_t50_raw_prior_mean, mu_t50_raw_prior_sd);
  sigma_t50  ~ normal(0, sigma_t50_prior_sd);
  t50_std    ~ normal(0, 1);

  mu_log_k   ~ normal(mu_log_k_prior_mean, mu_log_k_prior_sd);
  sigma_log_k ~ normal(0, sigma_log_k_prior_sd);
  log_k_std  ~ normal(0, 1);

  // ----------------------------
  // Unpooled priors for sigma_frac
  // ----------------------------
  sigma_frac ~ lognormal(sigma_frac_prior_meanlog, sigma_frac_prior_sdlog);

  // ----------------------------
  // Unpooled priors for lower/upper bounds
  // ----------------------------
  for (j in 1:J) {
    real s_l = inv_logit(lower_raw[j]);
    real s_u = inv_logit(upper_gap_raw[j]);

    target += normal_lpdf(lower_est[j] | lb_prior_mean[j], lb_prior_sd[j]);

    // Jacobian for lower_raw -> lower_est
    target += log(lower_cap[j] - lower_floor[j])
              + log(s_l)
              + log1m(s_l);

    target += normal_lpdf(upper_est[j] | ub_prior_mean[j], ub_prior_sd[j]);

    // Jacobian for upper_gap_raw -> upper_est, conditional on lower_est
    target += log(upper_cap[j] - lower_est[j])
              + log(s_u)
              + log1m(s_u);
  }

  // ----------------------------
  // EXTRA TARGETED PRIORS
  // ----------------------------
  for (j in 1:J) {
    if (use_t50_tweak[j] == 1) {
      target += normal_lpdf(t50_param[j] | t50_tweak_mean[j], t50_tweak_sd[j]);
    }
    if (use_logk_tweak[j] == 1) {
      target += normal_lpdf(log(k_param[j]) | logk_tweak_mean[j], logk_tweak_sd[j]);
    }
    if (use_upper_tweak[j] == 1) {
      target += normal_lpdf(upper_est[j] | upper_tweak_mean[j], upper_tweak_sd[j]);
    }
  }

  // ----------------------------
  // Observation model on original scale
  // ----------------------------
  for (n in 1:N) {
    int j = param_id[n];
    real q = q_norm(tau[n], t50_param[j], k_param[j]);
    real mu;

    if (increases[j] == 1) {
      mu = lower_est[j] + span_est[j] * q;
    } else {
      mu = upper_est[j] - span_est[j] * q;
    }

    y_obs[n] ~ student_t(
      4,
      mu,
      sigma_frac[j] * span_est[j] * obs_sd_mult[n]
    );
  }
}

generated quantities {
  matrix[J, N_pred] Q_pred;
  matrix[J, N_pred] theta_pred;
  vector[N] log_lik;

  for (j in 1:J) {
    for (m in 1:N_pred) {
      real q = q_norm(tau_pred[m], t50_param[j], k_param[j]);
      Q_pred[j, m] = q;

      if (increases[j] == 1) {
        theta_pred[j, m] = lower_est[j] + span_est[j] * q;
      } else {
        theta_pred[j, m] = upper_est[j] - span_est[j] * q;
      }
    }
  }

  for (n in 1:N) {
    int j = param_id[n];
    real q = q_norm(tau[n], t50_param[j], k_param[j]);
    real mu;

    if (increases[j] == 1) {
      mu = lower_est[j] + span_est[j] * q;
    } else {
      mu = upper_est[j] - span_est[j] * q;
    }

    log_lik[n] = student_t_lpdf(
      y_obs[n] |
      4,
      mu,
      sigma_frac[j] * span_est[j] * obs_sd_mult[n]
    );
  }
}
'

writeLines(stan_code, con = stan_file)

# ------------------------------------------------------------
# Read workbook
# ------------------------------------------------------------
anchors <- read_excel(
  workbook_path,
  sheet = sheet_name,
  skip = 2
)

# ------------------------------------------------------------
# Basic checks
# ------------------------------------------------------------
required_cols <- c(
  "anchor_id", "parameter", "relative_day", "value_used",
  "fit_role", "include_in_fit", "weight", "direction",
  "lower_bound", "upper_bound"
)

missing_cols <- setdiff(required_cols, names(anchors))
if (length(missing_cols) > 0) {
  stop(
    "These required columns are missing from the workbook sheet: ",
    paste(missing_cols, collapse = ", ")
  )
}

# ------------------------------------------------------------
# Clean workbook columns
# ------------------------------------------------------------
anchors <- anchors %>%
  mutate(
    parameter      = as.character(parameter),
    fit_role       = trimws(as.character(fit_role)),
    include_in_fit = toupper(trimws(as.character(include_in_fit))),
    direction      = tolower(trimws(as.character(direction))),
    anchor_id      = as.character(anchor_id),
    relative_day   = as.numeric(relative_day),
    value_used     = as.numeric(value_used),
    lower_bound    = as.numeric(lower_bound),
    upper_bound    = as.numeric(upper_bound),
    weight         = as.numeric(weight)
  )

# ------------------------------------------------------------
# Parameter order and labels
# ------------------------------------------------------------
param_levels <- c(
  "delay_hosp",
  "p_hosp",
  "p_ETU",
  "latent_IPC",
  "p_unsafe_funeral_comm",
  "p_unsafe_funeral_hosp"
)

panel_lookup <- c(
  "delay_hosp"             = "Delay to hospitalisation",
  "p_hosp"                 = "Probability of hospitalisation",
  "p_ETU"                  = "Proportion in ETU / ETC",
  "latent_IPC"             = "Latent IPC index",
  "p_unsafe_funeral_comm"  = "Unsafe funeral after community death",
  "p_unsafe_funeral_hosp"  = "Unsafe funeral after hospital death"
)

panel_order <- c(
  "Delay to hospitalisation",
  "Probability of hospitalisation",
  "Proportion in ETU / ETC",
  "Latent IPC index",
  "Unsafe funeral after community death",
  "Unsafe funeral after hospital death"
)

# ------------------------------------------------------------
# Rows used in the fit
# ------------------------------------------------------------
fit_df <- anchors %>%
  filter(parameter %in% param_levels) %>%
  filter(include_in_fit == "YES") %>%
  filter(!is.na(value_used)) %>%
  filter(!is.na(relative_day))

if (nrow(fit_df) == 0) {
  stop("No usable rows found after filtering include_in_fit == 'YES'.")
}

# ------------------------------------------------------------
# Parameter metadata
# Workbook lower/upper become prior centres
# ------------------------------------------------------------
param_meta <- fit_df %>%
  group_by(parameter) %>%
  summarise(
    lb_prior_mean = first(lower_bound),
    ub_prior_mean = first(upper_bound),
    direction     = first(direction),
    .groups = "drop"
  ) %>%
  left_join(domain_meta, by = "parameter") %>%
  mutate(
    param_id   = match(parameter, param_levels),
    increases  = if_else(direction == "up", 1L, 0L),
    span0      = ub_prior_mean - lb_prior_mean,
    domain_w   = abs_max - abs_min
  ) %>%
  arrange(param_id)

if (!all(param_levels %in% param_meta$parameter)) {
  missing_params <- setdiff(param_levels, param_meta$parameter)
  stop(
    "These fitted parameters were not found in the workbook sheet: ",
    paste(missing_params, collapse = ", ")
  )
}

if (any(is.na(param_meta$abs_min)) || any(is.na(param_meta$abs_max))) {
  stop("At least one parameter is missing abs_min / abs_max in domain_meta.")
}

if (any(param_meta$ub_prior_mean <= param_meta$lb_prior_mean)) {
  stop("At least one parameter has workbook upper_bound <= lower_bound.")
}

# ------------------------------------------------------------
# Construct admissible support for lower/upper estimates
# ------------------------------------------------------------
param_meta <- param_meta %>%
  mutate(
    lower_floor = pmax(abs_min, lb_prior_mean - support_half_width_mult * span0),

    lower_cap_raw = pmin(
      lb_prior_mean + support_half_width_mult * span0,
      ub_prior_mean - 0.05 * span0,
      abs_max - 0.05 * domain_w
    ),

    lower_cap = pmax(lower_floor + 1e-3, lower_cap_raw),

    upper_cap_raw = pmin(
      abs_max,
      ub_prior_mean + support_half_width_mult * span0
    ),

    upper_cap = pmax(lower_cap + 1e-3, upper_cap_raw),

    lb_prior_sd = pmax(
      bound_sd_frac_of_initial_span * span0,
      bound_sd_frac_of_domain * domain_w,
      1e-3
    ),

    ub_prior_sd = pmax(
      bound_sd_frac_of_initial_span * span0,
      bound_sd_frac_of_domain * domain_w,
      1e-3
    )
  )

if (any(param_meta$lower_cap <= param_meta$lower_floor)) {
  stop("At least one parameter has lower_cap <= lower_floor.")
}
if (any(param_meta$upper_cap <= param_meta$lower_cap)) {
  stop("At least one parameter has upper_cap <= lower_cap.")
}

# ------------------------------------------------------------
# Join row-level ids and directions
# ------------------------------------------------------------
fit_df <- fit_df %>%
  left_join(
    param_meta %>%
      select(parameter, param_id, increases),
    by = "parameter"
  )

# ------------------------------------------------------------
# Relative time
# ------------------------------------------------------------
max_day <- max(fit_df$relative_day, na.rm = TRUE)

fit_df <- fit_df %>%
  mutate(
    tau = relative_day / max_day,
    y_obs = value_used,
    weight = if_else(is.na(weight), 1.0, weight),
    obs_sd_mult = 1 / pmax(weight, 0.25)
  ) %>%
  arrange(param_id, tau)

# ------------------------------------------------------------
# OPTIONAL extra anchor downweight for WA_UFC_01
# ------------------------------------------------------------
if (use_anchor_downweight) {
  fit_df <- fit_df %>%
    mutate(
      obs_sd_mult = if_else(
        anchor_id == anchor_to_downweight & parameter == tweak_param,
        obs_sd_mult * anchor_downweight_sd_mult,
        obs_sd_mult
      )
    )
}

if (any(is.na(fit_df$y_obs))) stop("Some y_obs values are NA.")
if (any(is.na(fit_df$obs_sd_mult))) stop("Some obs_sd_mult values are NA.")

cat("\nRows per parameter used in Stan fit:\n")
print(table(fit_df$parameter))

cat("\nFit roles used:\n")
print(table(fit_df$fit_role, useNA = "ifany"))

# ------------------------------------------------------------
# Build targeted tweak vectors
# ------------------------------------------------------------
J <- nrow(param_meta)

use_t50_tweak_vec   <- rep(0L, J)
t50_tweak_mean_vec  <- rep(0.5, J)
t50_tweak_sd_vec    <- rep(1.0, J)

use_logk_tweak_vec  <- rep(0L, J)
logk_tweak_mean_vec <- rep(log(8), J)
logk_tweak_sd_vec   <- rep(1.0, J)

use_upper_tweak_vec <- rep(0L, J)
upper_tweak_mean_vec <- param_meta$ub_prior_mean
upper_tweak_sd_vec   <- rep(1.0, J)

j_tweak <- param_meta$param_id[param_meta$parameter == tweak_param]

if (length(j_tweak) != 1) {
  stop("Could not uniquely identify tweak parameter in param_meta.")
}

if (use_t50_tweak) {
  use_t50_tweak_vec[j_tweak]  <- 1L
  t50_tweak_mean_vec[j_tweak] <- t50_tweak_mean
  t50_tweak_sd_vec[j_tweak]   <- t50_tweak_sd
}

if (use_logk_tweak) {
  use_logk_tweak_vec[j_tweak]  <- 1L
  logk_tweak_mean_vec[j_tweak] <- logk_tweak_mean
  logk_tweak_sd_vec[j_tweak]   <- logk_tweak_sd
}

if (use_upper_tweak) {
  use_upper_tweak_vec[j_tweak]   <- 1L
  upper_tweak_mean_vec[j_tweak]  <- upper_tweak_mean
  upper_tweak_sd_vec[j_tweak]    <- upper_tweak_sd
}

# ------------------------------------------------------------
# Build Stan data
# ------------------------------------------------------------
tau_pred <- seq(0, 1, length.out = 100)

stan_data <- list(
  N = nrow(fit_df),
  J = nrow(param_meta),
  param_id = fit_df$param_id,
  tau = fit_df$tau,
  y_obs = fit_df$y_obs,
  obs_sd_mult = fit_df$obs_sd_mult,
  abs_min = param_meta$abs_min,
  abs_max = param_meta$abs_max,
  lower_floor = param_meta$lower_floor,
  lower_cap = param_meta$lower_cap,
  upper_cap = param_meta$upper_cap,
  lb_prior_mean = param_meta$lb_prior_mean,
  lb_prior_sd = param_meta$lb_prior_sd,
  ub_prior_mean = param_meta$ub_prior_mean,
  ub_prior_sd = param_meta$ub_prior_sd,
  increases = param_meta$increases,
  sigma_frac_prior_meanlog = sigma_frac_prior_meanlog,
  sigma_frac_prior_sdlog = sigma_frac_prior_sdlog,
  mu_t50_raw_prior_mean = mu_t50_raw_prior_mean,
  mu_t50_raw_prior_sd = mu_t50_raw_prior_sd,
  sigma_t50_prior_sd = sigma_t50_prior_sd,
  mu_log_k_prior_mean = mu_log_k_prior_mean,
  mu_log_k_prior_sd = mu_log_k_prior_sd,
  sigma_log_k_prior_sd = sigma_log_k_prior_sd,
  use_t50_tweak = use_t50_tweak_vec,
  t50_tweak_mean = t50_tweak_mean_vec,
  t50_tweak_sd = t50_tweak_sd_vec,
  use_logk_tweak = use_logk_tweak_vec,
  logk_tweak_mean = logk_tweak_mean_vec,
  logk_tweak_sd = logk_tweak_sd_vec,
  use_upper_tweak = use_upper_tweak_vec,
  upper_tweak_mean = upper_tweak_mean_vec,
  upper_tweak_sd = upper_tweak_sd_vec,
  N_pred = length(tau_pred),
  tau_pred = tau_pred
)

# ------------------------------------------------------------
# Compile and sample
# ------------------------------------------------------------
mod <- cmdstan_model(stan_file)

fit <- mod$sample(
  data = stan_data,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1500,
  iter_sampling = 1500,
  adapt_delta = 0.97,
  max_treedepth = 13,
  refresh = 200
)

# ------------------------------------------------------------
# Diagnostics
# ------------------------------------------------------------
cat("\nPosterior summaries for partial-pooled Q shape and unpooled bounds:\n")
print(
  fit$summary(
    variables = c(
      "mu_t50_raw", "sigma_t50",
      "mu_log_k", "sigma_log_k",
      "t50_param", "k_param",
      "lower_est", "upper_est", "sigma_frac"
    )
  )
)

cat("\nDiagnostic summary:\n")
print(fit$diagnostic_summary())

# ------------------------------------------------------------
# Extract parameter-specific Q summaries
# ------------------------------------------------------------
q_summ <- fit$summary(variables = "Q_pred") %>%
  mutate(
    param_id = as.integer(str_match(variable, "Q_pred\\[([0-9]+),([0-9]+)\\]")[, 2]),
    grid_id  = as.integer(str_match(variable, "Q_pred\\[([0-9]+),([0-9]+)\\]")[, 3])
  ) %>%
  left_join(
    param_meta %>% select(param_id, parameter),
    by = "param_id"
  ) %>%
  mutate(
    tau = tau_pred[grid_id],
    relative_day = tau * max_day,
    panel_title = factor(panel_lookup[parameter], levels = panel_order)
  ) %>%
  select(parameter, param_id, grid_id, tau, relative_day, mean, median, sd, q5, q95, panel_title)

# ------------------------------------------------------------
# Extract fitted parameter curves
# ------------------------------------------------------------
curve_summ <- fit$summary(variables = "theta_pred") %>%
  mutate(
    param_id = as.integer(str_match(variable, "theta_pred\\[([0-9]+),([0-9]+)\\]")[, 2]),
    grid_id  = as.integer(str_match(variable, "theta_pred\\[([0-9]+),([0-9]+)\\]")[, 3])
  ) %>%
  left_join(
    param_meta %>% select(param_id, parameter),
    by = "param_id"
  ) %>%
  mutate(
    tau = tau_pred[grid_id],
    relative_day = tau * max_day,
    panel_title = factor(panel_lookup[parameter], levels = panel_order)
  ) %>%
  select(parameter, param_id, grid_id, tau, relative_day, mean, median, sd, q5, q95, panel_title)

# ------------------------------------------------------------
# Extract lower/upper summaries
# ------------------------------------------------------------
lower_summ <- fit$summary(variables = "lower_est") %>%
  mutate(
    param_id = as.integer(str_match(variable, "lower_est\\[([0-9]+)\\]")[, 2]),
    bound = "lower"
  )

upper_summ <- fit$summary(variables = "upper_est") %>%
  mutate(
    param_id = as.integer(str_match(variable, "upper_est\\[([0-9]+)\\]")[, 2]),
    bound = "upper"
  )

bound_summ <- bind_rows(lower_summ, upper_summ) %>%
  left_join(
    param_meta %>% select(param_id, parameter),
    by = "param_id"
  ) %>%
  select(parameter, bound, mean, median, sd, q5, q95, rhat, ess_bulk, ess_tail)

# ------------------------------------------------------------
# Plot points
# ------------------------------------------------------------
plot_points <- fit_df %>%
  mutate(
    panel_title = factor(panel_lookup[parameter], levels = panel_order)
  )

# ------------------------------------------------------------
# Main six-panel plot
# ------------------------------------------------------------
p_main <- ggplot(curve_summ, aes(x = relative_day, y = mean)) +
  geom_ribbon(
    aes(ymin = q5, ymax = q95),
    fill = "#c7dcec",
    alpha = 0.9
  ) +
  geom_line(
    linewidth = 0.9,
    colour = "#1f77b4"
  ) +
  geom_point(
    data = plot_points,
    aes(x = relative_day, y = y_obs),
    inherit.aes = FALSE,
    colour = "#ff7f0e",
    size = 2.2
  ) +
  geom_text(
    data = plot_points,
    aes(x = relative_day, y = y_obs, label = anchor_id),
    inherit.aes = FALSE,
    hjust = 0,
    nudge_x = 2,
    size = 2.5,
    colour = "grey25"
  ) +
  facet_wrap(~ panel_title, scales = "free_y", ncol = 2) +
  theme_bw(base_size = 11) +
  labs(
    title = "West Africa: partial-pooled finite-window Q with targeted UFC-community tweak",
    subtitle = "Only p_unsafe_funeral_comm receives extra informative priors; all other parameters unchanged",
    x = "Relative outbreak day",
    y = NULL
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.minor = element_blank()
  )

print(p_main)

ggsave(
  filename = paste0(out_prefix, "_curves.png"),
  plot = p_main,
  width = 12,
  height = 10,
  dpi = 180
)

# ------------------------------------------------------------
# Latent Q_j plot
# ------------------------------------------------------------
p_q <- ggplot(q_summ, aes(x = relative_day, y = mean)) +
  geom_ribbon(aes(ymin = q5, ymax = q95), fill = "#dbe8f4", alpha = 0.9) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ panel_title, scales = "free_y", ncol = 2) +
  theme_bw(base_size = 11) +
  labs(
    title = "West Africa: parameter-specific finite-window Q_j(t)",
    subtitle = "Each Q_j starts at 0 and ends at 1; only UFC-community receives a targeted extra prior",
    x = "Relative outbreak day",
    y = "Q_j(t)"
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.minor = element_blank()
  )

print(p_q)

ggsave(
  filename = paste0(out_prefix, "_Qplot.png"),
  plot = p_q,
  width = 12,
  height = 10,
  dpi = 180
)

# ------------------------------------------------------------
# BP-ready matrix
# ------------------------------------------------------------
bp_matrix <- curve_summ %>%
  select(parameter, relative_day, tau, mean) %>%
  pivot_wider(names_from = parameter, values_from = mean) %>%
  mutate(
    prob_hosp = p_hosp,
    delay_hosp = delay_hosp,
    prob_unsafe_funeral_comm = p_unsafe_funeral_comm,
    prob_unsafe_funeral_hosp = p_unsafe_funeral_hosp,
    prob_unsafe_funeral_etu = 0,
    prop_etu = p_ETU,
    ipc_helper = latent_IPC
  ) %>%
  select(
    relative_day,
    tau,
    prob_hosp,
    delay_hosp,
    prob_unsafe_funeral_comm,
    prob_unsafe_funeral_hosp,
    prob_unsafe_funeral_etu,
    prop_etu,
    ipc_helper
  )

# ------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------
write_csv(curve_summ, paste0(out_prefix, "_curve_summaries.csv"))
write_csv(q_summ, paste0(out_prefix, "_Q_summaries.csv"))
write_csv(bound_summ, paste0(out_prefix, "_bound_summaries.csv"))
write_csv(bp_matrix, paste0(out_prefix, "_bp_input_matrix.csv"))
write_csv(plot_points, paste0(out_prefix, "_anchor_rows_used.csv"))
write_csv(param_meta, paste0(out_prefix, "_prior_setup.csv"))

tweak_summary <- tibble(
  tweak_param = tweak_param,
  use_t50_tweak = use_t50_tweak,
  t50_tweak_mean = t50_tweak_mean,
  t50_tweak_sd = t50_tweak_sd,
  use_logk_tweak = use_logk_tweak,
  logk_tweak_mean = logk_tweak_mean,
  logk_tweak_sd = logk_tweak_sd,
  use_upper_tweak = use_upper_tweak,
  upper_tweak_mean = upper_tweak_mean,
  upper_tweak_sd = upper_tweak_sd,
  use_anchor_downweight = use_anchor_downweight,
  anchor_to_downweight = anchor_to_downweight,
  anchor_downweight_sd_mult = anchor_downweight_sd_mult
)

write_csv(tweak_summary, paste0(out_prefix, "_tweak_settings.csv"))

cat("\nDone.\n")
cat("Files written:\n")
cat(" - ", stan_file, "\n", sep = "")
cat(" - ", paste0(out_prefix, "_curves.png"), "\n", sep = "")
cat(" - ", paste0(out_prefix, "_Qplot.png"), "\n", sep = "")
cat(" - ", paste0(out_prefix, "_curve_summaries.csv"), "\n", sep = "")
cat(" - ", paste0(out_prefix, "_Q_summaries.csv"), "\n", sep = "")
cat(" - ", paste0(out_prefix, "_bound_summaries.csv"), "\n", sep = "")
cat(" - ", paste0(out_prefix, "_bp_input_matrix.csv"), "\n", sep = "")
cat(" - ", paste0(out_prefix, "_anchor_rows_used.csv"), "\n", sep = "")
cat(" - ", paste0(out_prefix, "_prior_setup.csv"), "\n", sep = "")
cat(" - ", paste0(out_prefix, "_tweak_settings.csv"), "\n", sep = "")
