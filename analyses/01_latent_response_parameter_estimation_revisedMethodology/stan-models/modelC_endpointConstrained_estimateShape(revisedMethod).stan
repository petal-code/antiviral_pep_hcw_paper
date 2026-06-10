// ============================================================================
// MODEL C  (REVISED METHODOLOGY: endpoint-constrained partial pooling)
// ----------------------------------------------------------------------------
// Partial-pooled, finite-window logistic Q-curve model that estimates ONLY the
// SHAPE of each parameter's response-quality curve Q_j(tau). Unlike Model A
// (original methodology), the two per-parameter magnitude endpoints are NOT
// estimated here - they are supplied as fixed DATA (theta_start, theta_end),
// locked beforehand to early/late literature-window extrema (see
// lock_endpoints() in helpers(revisedMethod).R).
//
// This is the "West Africa" revised-methodology model. It sits between Model A
// and Model B:
//   * like Model A, it ESTIMATES the curve shape (t50, k), partially pooled
//     across parameters;
//   * like Model B, it treats part of the problem as fixed - but here it is the
//     ENDPOINTS that are fixed (Model B fixes the shared Q and estimates the
//     endpoints; Model C does the opposite).
//
// ---------------------------------------------------------------------------
// THE GENERATIVE STORY
// ---------------------------------------------------------------------------
// Outbreak progress is measured on a normalised clock tau in [0, 1]
// (tau = relative_day / endpoint_day). Each response parameter j moves from its
// fixed "worst-response" endpoint theta_start_j to its fixed "best-response"
// endpoint theta_end_j as the response matures, along a normalised logistic:
//
//     Q_j(tau) = q_norm(tau ; t50_j, k_j)            (runs exactly 0 -> 1)
//     mu_j(tau) = theta_start_j + (theta_end_j - theta_start_j) * Q_j(tau)
//
// Direction is already baked into the endpoints by the caller: an increasing
// parameter has theta_start < theta_end, a decreasing one has theta_start >
// theta_end. So the mean is the same simple interpolation for every parameter -
// no separate increasing/decreasing branch is needed (contrast Model A).
//
// ---------------------------------------------------------------------------
// PARTIAL POOLING (across PARAMETERS, not across data sources)
// ---------------------------------------------------------------------------
// The shape parameters (t50_j, k_j) are PARTIALLY POOLED across the J response
// parameters via a hierarchical (non-centred) prior, exactly as in Model A: we
// expect the parameters to improve on broadly similar timescales, so each
// parameter's (t50_j, k_j) is shrunk toward a shared group mean, with the
// shrinkage strength (sigma_t50, sigma_log_k) itself estimated. There are no
// magnitude endpoints to pool here - they are fixed data.
//
// ---------------------------------------------------------------------------
// LIKELIHOOD
// ---------------------------------------------------------------------------
// Student-t(4): robust to the occasional discordant literature anchor. The
// observation scale obs_sd is supplied per observation (built in 01 from each
// parameter's endpoint span and the anchor weight), so there is no separate
// per-parameter sigma to estimate.
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
  vector<lower=1e-9>[N] obs_sd;                // per-obs observation scale (built in 01)

  // ---- FIXED magnitude endpoints (the locked literature extrema) ---------
  // Direction is encoded here: increasing -> start < end, decreasing -> start > end.
  vector[J] theta_start;                       // value at Q = 0 (worst response)
  vector[J] theta_end;                         // value at Q = 1 (best response)

  // ---- Hyperpriors for the partially-pooled curve shape ------------------
  real mu_t50_raw_prior_mean;
  real<lower=0> mu_t50_raw_prior_sd;
  real<lower=0> sigma_t50_prior_sd;

  real mu_log_k_prior_mean;
  real<lower=0> mu_log_k_prior_sd;
  real<lower=0> sigma_log_k_prior_sd;

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
}

transformed parameters {
  vector<lower=0, upper=1>[J] t50_param;   // per-parameter midpoint in (0,1)
  vector<lower=0>[J] k_param;              // per-parameter steepness (>0)

  for (j in 1:J) {
    // Shape: pull each parameter's t50/k from the shared hierarchical mean.
    t50_param[j] = inv_logit(mu_t50_raw + sigma_t50 * t50_std[j]);
    k_param[j]   = exp(mu_log_k + sigma_log_k * log_k_std[j]);
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

  // ---- Observation model -------------------------------------------------
  // Student-t(4) likelihood with the fixed-endpoint interpolation as the mean.
  for (n in 1:N) {
    int j = param_id[n];
    real q = q_norm(tau[n], t50_param[j], k_param[j]);
    real mu = theta_start[j] + (theta_end[j] - theta_start[j]) * q;
    y_obs[n] ~ student_t(4, mu, obs_sd[n]);
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
      theta_pred[j, m] = theta_start[j] + (theta_end[j] - theta_start[j]) * q;
    }
  }

  for (n in 1:N) {
    int j = param_id[n];
    real q = q_norm(tau[n], t50_param[j], k_param[j]);
    real mu = theta_start[j] + (theta_end[j] - theta_start[j]) * q;
    log_lik[n] = student_t_lpdf(y_obs[n] | 4, mu, obs_sd[n]);
  }
}
