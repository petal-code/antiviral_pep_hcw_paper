// ============================================================================
// MODEL D -- BOTH response clocks SUPPLIED EMPIRICALLY (reach + SDB); only the
// endpoint magnitudes are estimated from the literature anchors (the Model B
// machinery). There is NO latent curve and NO smoothing model: the reach curve
// R(t) is built directly from the community-death data OUTSIDE Stan (exactly as
// the SDB success curve s(t) already is), and handed in as data. Stan therefore
// only solves the original-methodology endpoint problem (lower_j, upper_j from
// anchors) against two fixed clocks.
// ----------------------------------------------------------------------------
// Contrast with Model C: Model C estimated R(t) latently (random walk / spline
// on logit community-death proportion, with a binomial likelihood). That latent
// curve is smoothed by its prior and so irons out the early community-death
// shoulder. Model D removes the latent layer entirely -- R(t) tracks the data as
// tightly as the empirical construction allows.
//
// TWO SUPPLIED CLOCKS (both DATA):
//   * REACH  R_clock(t) = (c_hi - p_comm(t)) / (c_hi - c_lo)   in [0,1]  (scaling)
//            R_prob(t)  = 1 - p_comm(t)                                  (burial term)
//   * SDB    Sclk(t)    = s(t) / s_ref                                   (scaling)
//
// PARAMETER -> CLOCK (clock_type[]):
//   delay_hosp, p_hosp, p_ETU            -> reach clock  R_clock(t)
//   latent_IPC (toggle), uf_hosp         -> sdb   clock  Sclk(t)
//   p_unsafe_funeral_comm                -> DETERMINISTIC 1 - R_prob(t)*s(t)
//                                           (handled in generated quantities)
// ============================================================================

data {
  int<lower=2> G;                              // number of grid points

  // ---- supplied EMPIRICAL clocks on the grid (DATA, not estimated) --------
  vector<lower=0, upper=1>[G] R_clock;         // normalised reach clock (endpoint scaling)
  vector<lower=0, upper=1>[G] R_prob;          // reach probability (for the 1 - R*s term)
  vector<lower=0, upper=1>[G] s_grid;          // empirical SDB success on the grid
  real<lower=0, upper=1> s_ref;                // SDB clock normaliser (Sclk = s/s_ref)

  // ---- endpoint-parameter metadata (J params; uf_comm excluded) ----------
  int<lower=1> J;
  vector[J] abs_min;        vector[J] abs_max;
  vector[J] lower_floor;    vector[J] lower_cap;   vector[J] upper_cap;
  vector[J] lb_prior_mean;  vector[J] lb_prior_sd;
  vector[J] ub_prior_mean;  vector[J] ub_prior_sd;
  array[J] int<lower=0, upper=1> increases;
  array[J] int<lower=1, upper=2> clock_type;   // 1 = reach (R_clock), 2 = sdb (Sclk)

  // ---- anchors -----------------------------------------------------------
  int<lower=1> N;
  array[N] int<lower=1, upper=J> param_id;
  vector[N] y_obs;
  vector<lower=0>[N] obs_sd_mult;
  array[N] int<lower=1, upper=G> anchor_grid;  // grid index at each anchor's day

  // ---- priors ------------------------------------------------------------
  real sigma_frac_prior_meanlog;  real<lower=0> sigma_frac_prior_sdlog;

  // ---- optional endpoint tweaks (as Model B) -----------------------------
  array[J] int<lower=0, upper=1> use_upper_tweak; vector[J] upper_tweak_mean; vector<lower=0>[J] upper_tweak_sd;
  array[J] int<lower=0, upper=1> use_lower_tweak; vector[J] lower_tweak_mean; vector<lower=0>[J] lower_tweak_sd;
}

parameters {
  vector[J] lower_raw;
  vector[J] upper_gap_raw;
  vector<lower=0>[J] sigma_frac;
}

transformed parameters {
  vector[J] lower_est;
  vector[J] upper_est;
  vector[J] span_est;

  // Endpoint transforms (identical to Model B / Model C).
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

  // ---- anchor likelihood: clock value depends on the parameter's clock ---
  for (n in 1:N) {
    int j = param_id[n];
    real q = (clock_type[j] == 1)
               ? R_clock[anchor_grid[n]]                       // supplied reach clock
               : fmin(1, s_grid[anchor_grid[n]] / s_ref);      // supplied SDB clock
    real mu = (increases[j] == 1)
                ? lower_est[j] + span_est[j] * q
                : upper_est[j] - span_est[j] * q;
    y_obs[n] ~ student_t(4, mu, sigma_frac[j] * span_est[j] * obs_sd_mult[n]);
  }
}

generated quantities {
  // Per-endpoint-parameter native-unit curve on the grid (the harness assembles
  // these + the deterministic uf_comm into PARAM_LEVELS order). The reach and SDB
  // clocks are data, so they are simply plotted from R; only theta_pred and the
  // deterministic uf_comm depend on the estimated endpoints.
  matrix[J, G] theta_pred;
  vector[G] uf_comm_pred;                       // 1 - R_prob * s  (community unsafe funerals)

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
