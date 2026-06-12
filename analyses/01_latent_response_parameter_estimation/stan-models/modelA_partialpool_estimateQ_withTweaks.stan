// ============================================================================
// MODEL A (WITH optional targeted "tweak" priors)
// ----------------------------------------------------------------------------
// Partial-pooled, finite-window logistic Q-curve model that ESTIMATES BOTH:
//   (1) the shape of each parameter's response-quality curve  Q_j(tau), and
//   (2) the per-parameter magnitude endpoints (lower_j, upper_j).
//
// This is the "West Africa" model. It is used when we want the data (a handful
// of literature anchor points per parameter) to tell us how each response
// parameter evolves over the outbreak, WITHOUT supplying an external
// response-quality curve. (Contrast with Model B, where the Q curve is supplied
// as fixed data from the DRC safe-and-dignified-burial time series.)
//
// ---------------------------------------------------------------------------
// THE GENERATIVE STORY
// ---------------------------------------------------------------------------
// Outbreak progress is measured on a normalised clock tau in [0, 1]
// (tau = relative_day / max_anchor_day). Each response parameter j (e.g.
// probability of hospitalisation, delay to hospitalisation, ...) moves
// monotonically from a "worst-response" endpoint to a "best-response" endpoint
// as the response matures. The shape of that movement is a logistic curve,
// normalised so it starts exactly at 0 and ends exactly at 1 over the window:
//
//     Q_j(tau) = q_norm(tau ; t50_j, k_j)
//
// where t50_j is the midpoint (when the parameter is halfway through its change)
// and k_j is the steepness. The observed value of parameter j is then a linear
// interpolation between its two magnitude endpoints:
//
//     increasing parameter:  mu = lower_j + (upper_j - lower_j) * Q_j(tau)
//     decreasing parameter:  mu = upper_j - (upper_j - lower_j) * Q_j(tau)
//
// ---------------------------------------------------------------------------
// PARTIAL POOLING (this is across PARAMETERS, not across data sources)
// ---------------------------------------------------------------------------
// The shape parameters (t50_j, k_j) are PARTIALLY POOLED across the J response
// parameters via a hierarchical (non-centred) prior. Intuitively: we expect the
// different response parameters to improve on broadly similar timescales, so we
// shrink each parameter's (t50_j, k_j) toward a shared group mean. The amount of
// shrinkage is governed by sigma_t50 and sigma_log_k, which are themselves
// estimated. The magnitude endpoints (lower_j, upper_j) are NOT pooled - they
// live on different scales (days vs probabilities) and are estimated
// independently with informative, workbook-derived priors.
//
// ---------------------------------------------------------------------------
// "TWEAK" PRIORS (the part to scrutinise)
// ---------------------------------------------------------------------------
// On top of the model above, this version exposes OPTIONAL per-parameter
// "tweak" priors. Each is an extra informative normal prior placed on a single
// parameter's t50, log(k), upper endpoint, or lower endpoint. They let a chosen
// parameter "escape" the partial-pooling shrinkage when there is external
// knowledge about its shape or level. They are switched on per parameter via the
// use_*_tweak[] integer flags (0 = off, 1 = on); when all are 0 this model is
// numerically identical to the no-tweaks version. These priors are researcher
// choices and should each be justifiable from evidence, not from curve
// aesthetics - hence the separate no-tweaks model for comparison.
// ============================================================================

functions {
  // Raw (un-normalised) logistic response-quality curve.
  real q_raw(real tau, real t50, real k) {
    return inv_logit(k * (tau - t50));
  }

  // Finite-window normalisation: rescale the logistic so that it equals exactly
  // 0 at tau = 0 and exactly 1 at tau = 1. This makes Q a clean 0->1 progress
  // index over the observed outbreak window, regardless of (t50, k).
  real q_norm(real tau, real t50, real k) {
    real q0 = q_raw(0, t50, k);
    real q1 = q_raw(1, t50, k);
    return (q_raw(tau, t50, k) - q0) / (q1 - q0);
  }
}

data {
  // ---- Sizes -------------------------------------------------------------
  int<lower=1> N;                              // number of anchor observations
  int<lower=1> J;                              // number of response parameters

  // ---- Observations ------------------------------------------------------
  array[N] int<lower=1, upper=J> param_id;     // which parameter each obs belongs to
  vector<lower=0, upper=1>[N] tau;             // normalised outbreak time of each obs
  vector[N] y_obs;                             // observed (literature anchor) value
  vector<lower=0>[N] obs_sd_mult;              // per-obs noise multiplier (1/weight)

  // ---- Hard admissible support per parameter -----------------------------
  vector[J] abs_min;                           // floor of what is physically sensible
  vector[J] abs_max;                           // ceiling of what is physically sensible

  // ---- Transform limits for the estimated endpoints ----------------------
  vector[J] lower_floor;                       // lower_j is constrained to
  vector[J] lower_cap;                         //   [lower_floor, lower_cap]
  vector[J] upper_cap;                         // upper_j is constrained to [lower_j, upper_cap]

  // ---- Priors on the endpoints (centres from the evidence workbook) ------
  vector[J] lb_prior_mean;
  vector[J] lb_prior_sd;
  vector[J] ub_prior_mean;
  vector[J] ub_prior_sd;

  // ---- Direction of each parameter ---------------------------------------
  array[J] int<lower=0, upper=1> increases;    // 1 = improves upward, 0 = downward

  // ---- Prior for the per-parameter observation scale ---------------------
  // sigma_y_j = sigma_frac_j * (upper_j - lower_j); sigma_frac_j ~ lognormal(.)
  real sigma_frac_prior_meanlog;
  real<lower=0> sigma_frac_prior_sdlog;

  // ---- Hyperpriors for the partially-pooled curve shape ------------------
  real mu_t50_raw_prior_mean;
  real<lower=0> mu_t50_raw_prior_sd;
  real<lower=0> sigma_t50_prior_sd;

  real mu_log_k_prior_mean;
  real<lower=0> mu_log_k_prior_sd;
  real<lower=0> sigma_log_k_prior_sd;

  // ---- OPTIONAL targeted tweak priors (per parameter) --------------------
  // Each block: a 0/1 switch, a prior mean, and a prior sd, for every parameter.
  // Only entries with use_*_tweak[j] == 1 contribute to the posterior.
  array[J] int<lower=0, upper=1> use_t50_tweak;   // extra prior on t50_j
  vector<lower=0, upper=1>[J] t50_tweak_mean;
  vector<lower=0>[J] t50_tweak_sd;

  array[J] int<lower=0, upper=1> use_logk_tweak;  // extra prior on log(k_j)
  vector[J] logk_tweak_mean;
  vector<lower=0>[J] logk_tweak_sd;

  array[J] int<lower=0, upper=1> use_upper_tweak; // extra prior on upper_j
  vector[J] upper_tweak_mean;
  vector<lower=0>[J] upper_tweak_sd;

  array[J] int<lower=0, upper=1> use_lower_tweak; // extra prior on lower_j
  vector[J] lower_tweak_mean;
  vector<lower=0>[J] lower_tweak_sd;

  // ---- Prediction grid ---------------------------------------------------
  int<lower=1> N_pred;
  vector<lower=0, upper=1>[N_pred] tau_pred;
}

parameters {
  // Partially-pooled curve-shape parameters (non-centred parameterisation).
  real mu_t50_raw;                 // group-mean midpoint (on the logit scale)
  real<lower=0> sigma_t50;         // between-parameter sd of the midpoint
  vector[J] t50_std;               // standardised per-parameter midpoint offsets

  real mu_log_k;                   // group-mean log-steepness
  real<lower=0> sigma_log_k;       // between-parameter sd of log-steepness
  vector[J] log_k_std;             // standardised per-parameter steepness offsets

  // Independent (un-pooled) per-parameter magnitude endpoints, in unconstrained
  // space; transformed to the admissible window in `transformed parameters`.
  vector[J] lower_raw;
  vector[J] upper_gap_raw;

  // Independent per-parameter observation scale (as a fraction of the span).
  vector<lower=0>[J] sigma_frac;
}

transformed parameters {
  vector<lower=0, upper=1>[J] t50_param;   // per-parameter midpoint in (0,1)
  vector<lower=0>[J] k_param;              // per-parameter steepness (>0)
  vector[J] lower_est;                     // per-parameter lower magnitude endpoint
  vector[J] upper_est;                     // per-parameter upper magnitude endpoint
  vector[J] span_est;                      // upper_est - lower_est

  for (j in 1:J) {
    real s_l;
    real s_u;
    real lower_j;
    real upper_j;

    // Shape: pull each parameter's t50/k from the shared hierarchical mean.
    t50_param[j] = inv_logit(mu_t50_raw + sigma_t50 * t50_std[j]);
    k_param[j]   = exp(mu_log_k + sigma_log_k * log_k_std[j]);

    // Lower endpoint: squashed into [lower_floor, lower_cap].
    s_l = inv_logit(lower_raw[j]);
    lower_j = lower_floor[j] + (lower_cap[j] - lower_floor[j]) * s_l;

    // Upper endpoint: squashed into [lower_j, upper_cap] (always above lower_j).
    s_u = inv_logit(upper_gap_raw[j]);
    upper_j = lower_j + (upper_cap[j] - lower_j) * s_u;

    lower_est[j] = lower_j;
    upper_est[j] = upper_j;
    span_est[j]  = upper_j - lower_j;
  }
}

model {
  // ---- Hyperpriors for the partially-pooled shape ------------------------
  mu_t50_raw  ~ normal(mu_t50_raw_prior_mean, mu_t50_raw_prior_sd);
  sigma_t50   ~ normal(0, sigma_t50_prior_sd);     // half-normal (sigma_t50 >= 0)
  t50_std     ~ normal(0, 1);

  mu_log_k    ~ normal(mu_log_k_prior_mean, mu_log_k_prior_sd);
  sigma_log_k ~ normal(0, sigma_log_k_prior_sd);   // half-normal
  log_k_std   ~ normal(0, 1);

  // ---- Prior for the observation scale -----------------------------------
  sigma_frac ~ lognormal(sigma_frac_prior_meanlog, sigma_frac_prior_sdlog);

  // ---- Priors on the endpoints, with manual Jacobian corrections ---------
  // We placed normal priors directly on the CONSTRAINED endpoints (lower_est,
  // upper_est), but sampled the UNCONSTRAINED lower_raw / upper_gap_raw. The
  // change-of-variables therefore needs an explicit log-Jacobian term for each.
  for (j in 1:J) {
    real s_l = inv_logit(lower_raw[j]);
    real s_u = inv_logit(upper_gap_raw[j]);

    target += normal_lpdf(lower_est[j] | lb_prior_mean[j], lb_prior_sd[j]);
    // d(lower_est)/d(lower_raw) = (lower_cap - lower_floor) * s_l * (1 - s_l)
    target += log(lower_cap[j] - lower_floor[j]) + log(s_l) + log1m(s_l);

    target += normal_lpdf(upper_est[j] | ub_prior_mean[j], ub_prior_sd[j]);
    // d(upper_est)/d(upper_gap_raw) = (upper_cap - lower_est) * s_u * (1 - s_u)
    target += log(upper_cap[j] - lower_est[j]) + log(s_u) + log1m(s_u);
  }

  // ---- OPTIONAL targeted tweak priors ------------------------------------
  // Extra informative priors on selected parameters. Each is a researcher choice.
  for (j in 1:J) {
    if (use_t50_tweak[j] == 1) {
      target += normal_lpdf(t50_param[j]      | t50_tweak_mean[j],   t50_tweak_sd[j]);
    }
    if (use_logk_tweak[j] == 1) {
      target += normal_lpdf(log(k_param[j])   | logk_tweak_mean[j],  logk_tweak_sd[j]);
    }
    if (use_upper_tweak[j] == 1) {
      target += normal_lpdf(upper_est[j]      | upper_tweak_mean[j], upper_tweak_sd[j]);
    }
    if (use_lower_tweak[j] == 1) {
      target += normal_lpdf(lower_est[j]      | lower_tweak_mean[j], lower_tweak_sd[j]);
    }
  }

  // ---- Observation model -------------------------------------------------
  // Student-t(4) likelihood: robust to the occasional discordant literature
  // anchor. The noise scale is proportional to each parameter's own span.
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
  // Posterior predictive curves on the regular prediction grid:
  //   Q_pred[j, m]      = the normalised response-quality curve for parameter j
  //   theta_pred[j, m]  = the parameter's value (native units) at that grid point
  matrix[J, N_pred] Q_pred;
  matrix[J, N_pred] theta_pred;
  vector[N] log_lik;          // pointwise log-likelihood (for LOO / model checks)

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
