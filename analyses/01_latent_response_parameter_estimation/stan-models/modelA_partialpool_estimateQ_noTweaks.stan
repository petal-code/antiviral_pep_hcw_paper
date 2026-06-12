// ============================================================================
// MODEL A (NO tweak priors) -- the clean baseline
// ----------------------------------------------------------------------------
// This is identical to modelA_partialpool_estimateQ_withTweaks.stan EXCEPT that
// it contains NONE of the optional "tweak" priors: no use_*_tweak data inputs,
// and no extra informative priors in the model block. Everything that the data
// say about each parameter's curve shape and endpoints comes only from the
// literature anchors, the workbook-derived endpoint priors, and the partial
// pooling across parameters.
//
// It exists so we can fit the West Africa data with and without the targeted
// tweaks and overlay the two fits parameter-by-parameter (see
// west_africa_checking.R). The difference between the two fits is exactly "what
// the tweaks are doing".
//
// See the with-tweaks file for the full description of the generative model and
// the partial-pooling structure; the comments below are kept brief to avoid
// duplication.
// ============================================================================

functions {
  // Raw logistic, and the finite-window normalisation making Q(0)=0, Q(1)=1.
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

  int<lower=1> N_pred;
  vector<lower=0, upper=1>[N_pred] tau_pred;
}

parameters {
  real mu_t50_raw;
  real<lower=0> sigma_t50;
  vector[J] t50_std;

  real mu_log_k;
  real<lower=0> sigma_log_k;
  vector[J] log_k_std;

  vector[J] lower_raw;
  vector[J] upper_gap_raw;

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
  // Hyperpriors for the partially-pooled shape.
  mu_t50_raw  ~ normal(mu_t50_raw_prior_mean, mu_t50_raw_prior_sd);
  sigma_t50   ~ normal(0, sigma_t50_prior_sd);
  t50_std     ~ normal(0, 1);

  mu_log_k    ~ normal(mu_log_k_prior_mean, mu_log_k_prior_sd);
  sigma_log_k ~ normal(0, sigma_log_k_prior_sd);
  log_k_std   ~ normal(0, 1);

  // Observation scale.
  sigma_frac ~ lognormal(sigma_frac_prior_meanlog, sigma_frac_prior_sdlog);

  // Endpoint priors with manual Jacobian corrections (see with-tweaks file).
  for (j in 1:J) {
    real s_l = inv_logit(lower_raw[j]);
    real s_u = inv_logit(upper_gap_raw[j]);

    target += normal_lpdf(lower_est[j] | lb_prior_mean[j], lb_prior_sd[j]);
    target += log(lower_cap[j] - lower_floor[j]) + log(s_l) + log1m(s_l);

    target += normal_lpdf(upper_est[j] | ub_prior_mean[j], ub_prior_sd[j]);
    target += log(upper_cap[j] - lower_est[j]) + log(s_u) + log1m(s_u);
  }

  // Observation model: robust Student-t(4), noise proportional to the span.
  for (n in 1:N) {
    int j = param_id[n];
    real q = q_norm(tau[n], t50_param[j], k_param[j]);
    real mu;

    if (increases[j] == 1) {
      mu = lower_est[j] + span_est[j] * q;
    } else {
      mu = upper_est[j] - span_est[j] * q;
    }

    y_obs[n] ~ student_t(4, mu, sigma_frac[j] * span_est[j] * obs_sd_mult[n]);
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
    log_lik[n] = student_t_lpdf(y_obs[n] | 4, mu, sigma_frac[j] * span_est[j] * obs_sd_mult[n]);
  }
}
