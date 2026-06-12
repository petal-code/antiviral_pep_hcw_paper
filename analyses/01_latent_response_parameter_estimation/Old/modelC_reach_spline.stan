// ============================================================================
// MODEL C (SPLINE) -- identical to modelC_reach_rw.stan EXCEPT the latent reach
// curve is a penalised B-spline (P-spline) on logit(community-death proportion)
// instead of a random walk. Same two-clock structure, same supplied SDB curve,
// same endpoint estimation from anchors.
//
//   logit p_comm(t) = B(t) * beta      with a 2nd-order-difference penalty on beta
//                                       (smoothness controlled by sigma_beta).
//
// The basis matrix B (G x K) is built in the harness (splines::bs on the grid
// days) and supplied as data. Everything downstream (R_prob, R_clock, endpoints,
// anchor likelihood, generated quantities) is the same as the RW model.
// ============================================================================

data {
  int<lower=2> G;
  // (no grid_dt -- smoothness lives in the penalty, not a per-step scaling)

  int<lower=1> Mc;
  array[Mc] int<lower=1, upper=G> cd_grid;
  array[Mc] int<lower=0> cd_n;
  array[Mc] int<lower=1> cd_N;

  real<lower=0, upper=1> c_hi;
  real<lower=0, upper=1> c_lo;

  vector<lower=0, upper=1>[G] s_grid;
  real<lower=0, upper=1> s_ref;

  // ---- spline basis ------------------------------------------------------
  int<lower=3> K;                              // number of basis functions
  matrix[G, K] B;                              // B-spline basis on the grid days

  int<lower=1> J;
  vector[J] abs_min;        vector[J] abs_max;
  vector[J] lower_floor;    vector[J] lower_cap;   vector[J] upper_cap;
  vector[J] lb_prior_mean;  vector[J] lb_prior_sd;
  vector[J] ub_prior_mean;  vector[J] ub_prior_sd;
  array[J] int<lower=0, upper=1> increases;
  array[J] int<lower=1, upper=2> clock_type;

  int<lower=1> N;
  array[N] int<lower=1, upper=J> param_id;
  vector[N] y_obs;
  vector<lower=0>[N] obs_sd_mult;
  array[N] int<lower=1, upper=G> anchor_grid;

  real sigma_frac_prior_meanlog;  real<lower=0> sigma_frac_prior_sdlog;
  real<lower=0> spline_sigma_prior_sd;         // half-normal scale for the P-spline penalty sd

  array[J] int<lower=0, upper=1> use_upper_tweak; vector[J] upper_tweak_mean; vector<lower=0>[J] upper_tweak_sd;
  array[J] int<lower=0, upper=1> use_lower_tweak; vector[J] lower_tweak_mean; vector<lower=0>[J] lower_tweak_sd;
}

parameters {
  vector[K] beta;                              // spline coefficients
  real<lower=0> sigma_beta;                    // P-spline smoothness sd (2nd diffs)

  vector[J] lower_raw;
  vector[J] upper_gap_raw;
  vector<lower=0>[J] sigma_frac;
}

transformed parameters {
  vector[G] pc_logit = B * beta;
  vector<lower=0, upper=1>[G] p_comm;
  vector<lower=0, upper=1>[G] R_prob;
  vector<lower=0, upper=1>[G] R_clock;
  vector[J] lower_est;
  vector[J] upper_est;
  vector[J] span_est;

  for (g in 1:G) {
    p_comm[g]  = inv_logit(pc_logit[g]);
    R_prob[g]  = 1 - p_comm[g];
    R_clock[g] = fmin(1, fmax(0, (c_hi - p_comm[g]) / (c_hi - c_lo)));
  }

  for (j in 1:J) {
    real s_l = inv_logit(lower_raw[j]);
    real lower_j = lower_floor[j] + (lower_cap[j] - lower_floor[j]) * s_l;
    real s_u = inv_logit(upper_gap_raw[j]);
    real upper_j = lower_j + (upper_cap[j] - lower_j) * s_u;
    lower_est[j] = lower_j;
    upper_est[j] = upper_j;
    span_est[j]  = upper_j - lower_j;
  }
}

model {
  // ---- P-spline: weak priors on the first two coefs, 2nd-order-difference
  //      random-walk penalty on the rest (smoothness ~ sigma_beta) ----------
  beta[1] ~ normal(0, 2);
  beta[2] ~ normal(0, 2);
  sigma_beta ~ normal(0, spline_sigma_prior_sd);
  for (k in 3:K) beta[k] ~ normal(2 * beta[k - 1] - beta[k - 2], sigma_beta);

  // ---- community-death binomial likelihood ------------------------------
  for (i in 1:Mc) cd_n[i] ~ binomial(cd_N[i], p_comm[cd_grid[i]]);

  // ---- endpoint priors + Jacobians (Model B) ----------------------------
  sigma_frac ~ lognormal(sigma_frac_prior_meanlog, sigma_frac_prior_sdlog);
  for (j in 1:J) {
    real s_l = inv_logit(lower_raw[j]);
    real s_u = inv_logit(upper_gap_raw[j]);
    target += normal_lpdf(lower_est[j] | lb_prior_mean[j], lb_prior_sd[j]);
    target += log(lower_cap[j] - lower_floor[j]) + log(s_l) + log1m(s_l);
    target += normal_lpdf(upper_est[j] | ub_prior_mean[j], ub_prior_sd[j]);
    target += log(upper_cap[j] - lower_est[j]) + log(s_u) + log1m(s_u);
    if (use_upper_tweak[j] == 1) target += normal_lpdf(upper_est[j] | upper_tweak_mean[j], upper_tweak_sd[j]);
    if (use_lower_tweak[j] == 1) target += normal_lpdf(lower_est[j] | lower_tweak_mean[j], lower_tweak_sd[j]);
  }

  // ---- anchor likelihood -------------------------------------------------
  for (n in 1:N) {
    int j = param_id[n];
    real q = (clock_type[j] == 1)
               ? R_clock[anchor_grid[n]]
               : fmin(1, s_grid[anchor_grid[n]] / s_ref);
    real mu = (increases[j] == 1)
                ? lower_est[j] + span_est[j] * q
                : upper_est[j] - span_est[j] * q;
    y_obs[n] ~ student_t(4, mu, sigma_frac[j] * span_est[j] * obs_sd_mult[n]);
  }
}

generated quantities {
  matrix[J, G] theta_pred;
  vector[G] reach_prob   = R_prob;
  vector[G] reach_clock  = R_clock;
  vector[G] comm_death_p = p_comm;
  vector[G] uf_comm_pred;

  for (j in 1:J) {
    for (g in 1:G) {
      real q = (clock_type[j] == 1) ? R_clock[g] : fmin(1, s_grid[g] / s_ref);
      theta_pred[j, g] = (increases[j] == 1)
                           ? lower_est[j] + span_est[j] * q
                           : upper_est[j] - span_est[j] * q;
    }
  }
  for (g in 1:G) uf_comm_pred[g] = fmin(1, fmax(0, 1 - R_prob[g] * s_grid[g]));
}
