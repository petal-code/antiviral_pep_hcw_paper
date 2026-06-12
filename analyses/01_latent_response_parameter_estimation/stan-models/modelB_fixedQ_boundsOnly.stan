// ============================================================================
// MODEL B (fixed Q supplied as data; estimate magnitude endpoints only)
// ----------------------------------------------------------------------------
// This is the "DRC" model. Unlike Model A, it does NOT estimate the shape of the
// response-quality curve. Instead, a single shared Q curve is supplied as DATA
// (q_obs at the anchor times, q_pred on the prediction grid). That curve comes
// from the empirical DRC safe-and-dignified-burial (SDB / Warsame) success time
// series, reconstructed and smoothed in 00_DataPreparation_and_Cleaning.R.
//
// WHY NO CURVE ESTIMATION, AND WHY NO PARTIAL POOLING?
//   The DRC "conflict" scenarios are defined by the jagged, conflict-interrupted
//   shape of the empirical SDB curve. A smooth logistic (Model A) would erase
//   exactly the feature we care about, so here the empirical Q shape is taken as
//   fixed and shared by every parameter. Because there is no shape to estimate,
//   there is nothing to pool: the only free quantities are each parameter's two
//   magnitude endpoints, which live on different scales and are deliberately
//   estimated INDEPENDENTLY (with informative, workbook-derived priors). Pooling
//   them would be wrong (it would shrink a delay-in-days toward a probability).
//
// THE MODEL
//   The shared Q in [0,1] says how far along the response is at each time. Each
//   parameter j rides that shared clock between its own endpoints:
//       increasing parameter:  mu = lower_j + (upper_j - lower_j) * Q
//       decreasing parameter:  mu = upper_j - (upper_j - lower_j) * Q
//   so an anchor near Q=0 informs the worst-response endpoint and an anchor near
//   Q=1 informs the best-response endpoint.
//
// TWEAK PRIORS
//   Optional informative priors on the upper and lower endpoints are exposed
//   (use_upper_tweak[], use_lower_tweak[]). The DRC fit uses these to (a) anchor
//   community unsafe funerals to begin at 1 when Q=0, and to its absolute
//   Warsame floor (1 - max success) at Q=1, and (b) strongly regularise the
//   near-zero hospital unsafe-funeral parameter. All off when the flags are 0.
// ============================================================================

data {
  // ---- Sizes -------------------------------------------------------------
  int<lower=1> N;                              // number of anchor observations
  int<lower=1> J;                              // number of response parameters
  int<lower=1> N_pred;                         // size of the prediction grid

  // ---- Observations ------------------------------------------------------
  array[N] int<lower=1, upper=J> param_id;     // which parameter each obs belongs to
  vector<lower=0, upper=1>[N] q_obs;           // SHARED Q evaluated at each obs time (DATA)
  vector[N] y_obs;                             // observed (literature anchor) value
  vector<lower=0>[N] obs_sd_mult;              // per-obs noise multiplier (1/weight)

  // ---- Hard admissible support and transform limits ----------------------
  vector[J] abs_min;
  vector[J] abs_max;
  vector[J] lower_floor;                       // lower_j in [lower_floor, lower_cap]
  vector[J] lower_cap;
  vector[J] upper_cap;                         // upper_j in [lower_j, upper_cap]

  // ---- Priors on the endpoints (centres from the evidence workbook) ------
  vector[J] lb_prior_mean;
  vector[J] lb_prior_sd;
  vector[J] ub_prior_mean;
  vector[J] ub_prior_sd;

  // ---- Direction of each parameter ---------------------------------------
  array[J] int<lower=0, upper=1> increases;    // 1 = improves upward, 0 = downward

  // ---- Prior for the per-parameter observation scale ---------------------
  real sigma_frac_prior_meanlog;
  real<lower=0> sigma_frac_prior_sdlog;

  // ---- OPTIONAL targeted endpoint tweak priors ---------------------------
  array[J] int<lower=0, upper=1> use_upper_tweak;
  vector[J] upper_tweak_mean;
  vector<lower=0>[J] upper_tweak_sd;

  array[J] int<lower=0, upper=1> use_lower_tweak;
  vector[J] lower_tweak_mean;
  vector<lower=0>[J] lower_tweak_sd;

  // ---- Shared Q evaluated on the prediction grid (DATA) ------------------
  vector<lower=0, upper=1>[N_pred] q_pred;
}

parameters {
  // Only the per-parameter magnitude endpoints and observation scale are free.
  // No shape parameters, no hierarchy.
  vector[J] lower_raw;
  vector[J] upper_gap_raw;
  vector<lower=0>[J] sigma_frac;
}

transformed parameters {
  vector[J] lower_est;
  vector[J] upper_est;
  vector[J] span_est;

  for (j in 1:J) {
    real s_l;
    real s_u;
    real lower_j;
    real upper_j;

    // Lower endpoint squashed into [lower_floor, lower_cap].
    s_l = inv_logit(lower_raw[j]);
    lower_j = lower_floor[j] + (lower_cap[j] - lower_floor[j]) * s_l;

    // Upper endpoint squashed into [lower_j, upper_cap].
    s_u = inv_logit(upper_gap_raw[j]);
    upper_j = lower_j + (upper_cap[j] - lower_j) * s_u;

    lower_est[j] = lower_j;
    upper_est[j] = upper_j;
    span_est[j]  = upper_j - lower_j;
  }
}

model {
  // Observation scale.
  sigma_frac ~ lognormal(sigma_frac_prior_meanlog, sigma_frac_prior_sdlog);

  // Endpoint priors with manual Jacobian corrections (same logic as Model A).
  for (j in 1:J) {
    real s_l = inv_logit(lower_raw[j]);
    real s_u = inv_logit(upper_gap_raw[j]);

    target += normal_lpdf(lower_est[j] | lb_prior_mean[j], lb_prior_sd[j]);
    target += log(lower_cap[j] - lower_floor[j]) + log(s_l) + log1m(s_l);

    target += normal_lpdf(upper_est[j] | ub_prior_mean[j], ub_prior_sd[j]);
    target += log(upper_cap[j] - lower_est[j]) + log(s_u) + log1m(s_u);
  }

  // Optional targeted endpoint tweaks.
  for (j in 1:J) {
    if (use_upper_tweak[j] == 1) {
      target += normal_lpdf(upper_est[j] | upper_tweak_mean[j], upper_tweak_sd[j]);
    }
    if (use_lower_tweak[j] == 1) {
      target += normal_lpdf(lower_est[j] | lower_tweak_mean[j], lower_tweak_sd[j]);
    }
  }

  // Observation model: the shared, FIXED q_obs drives every parameter; only the
  // endpoints (and hence the span) are learned.
  for (n in 1:N) {
    int j = param_id[n];
    real mu;

    if (increases[j] == 1) {
      mu = lower_est[j] + span_est[j] * q_obs[n];
    } else {
      mu = upper_est[j] - span_est[j] * q_obs[n];
    }

    y_obs[n] ~ student_t(4, mu, sigma_frac[j] * span_est[j] * obs_sd_mult[n]);
  }
}

generated quantities {
  // Predicted curves. Q_pred is the SAME supplied shared curve for every
  // parameter (carried through for convenience); theta_pred is each parameter's
  // value implied by riding that shared curve between its fitted endpoints.
  matrix[J, N_pred] Q_pred;
  matrix[J, N_pred] theta_pred;
  vector[N] log_lik;

  for (j in 1:J) {
    for (m in 1:N_pred) {
      real q = q_pred[m];
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
    real mu;
    if (increases[j] == 1) {
      mu = lower_est[j] + span_est[j] * q_obs[n];
    } else {
      mu = upper_est[j] - span_est[j] * q_obs[n];
    }
    log_lik[n] = student_t_lpdf(y_obs[n] | 4, mu, sigma_frac[j] * span_est[j] * obs_sd_mult[n]);
  }
}
