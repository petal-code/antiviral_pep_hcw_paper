// ============================================================================
// MODEL C (RW) -- latent REACH curve estimated; SDB curve supplied; endpoints
// estimated from anchors (original-methodology machinery, reused from Model B).
// ----------------------------------------------------------------------------
// TWO RESPONSE CLOCKS (vs Model B's single SDB clock):
//   * REACH  R(t) -- LATENT, estimated here from the community-death counts via a
//                    random walk on logit(community-death proportion). Drives the
//                    care-access parameters (and optionally latent_IPC).
//   * SDB    s(t) -- SUPPLIED as DATA (the empirical Warsame success curve), as in
//                    Model B. Drives latent_IPC (default) / p_unsafe_funeral_hosp.
//
// PARAMETER -> CLOCK (set per parameter via clock_type[]):
//   delay_hosp, p_hosp, p_ETU            -> reach clock  R_clock(t)   (estimated)
//   latent_IPC (default), uf_hosp        -> sdb   clock  Sclk(t)=s/s_ref (data)
//   p_unsafe_funeral_comm                -> DETERMINISTIC 1 - R_prob(t)*s(t)
//                                           (NOT endpoint-scaled; handled in R from
//                                            reach_prob + s_grid -- excluded from J)
//
// Two reach quantities, deliberately distinct:
//   R_prob(t)  = 1 - p_comm(t)                       actual reach PROBABILITY
//                                                    (used in the 1 - R*s burial term)
//   R_clock(t) = (c_hi - p_comm(t)) / (c_hi - c_lo)  normalised [0,1] CLOCK
//                                                    (used for endpoint scaling; 4a)
//
// Endpoints (lower_j, upper_j) for the J endpoint-parameters are estimated from
// the literature anchors exactly as Model B (same transforms, Jacobians, tweaks).
// ============================================================================

data {
  // ---- latent reach grid (regular spacing assumed) -----------------------
  int<lower=2> G;                              // number of reach-grid points
  real<lower=0> grid_dt;                       // grid spacing in days (RW scaling)

  // ---- community-death observations (binomial) ---------------------------
  int<lower=1> Mc;
  array[Mc] int<lower=1, upper=G> cd_grid;     // nearest reach-grid index per obs
  array[Mc] int<lower=0> cd_n;                 // community deaths (numerator)
  array[Mc] int<lower=1> cd_N;                 // all deaths (denominator)

  // ---- reach-clock normalisation (4a) ------------------------------------
  real<lower=0, upper=1> c_hi;                 // worst community-death prop -> R_clock=0
  real<lower=0, upper=1> c_lo;                 // best  community-death prop -> R_clock=1

  // ---- supplied SDB curve (DATA) on the grid -----------------------------
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
  array[N] int<lower=1, upper=G> anchor_grid;  // reach-grid index at each anchor's day

  // ---- priors ------------------------------------------------------------
  real sigma_frac_prior_meanlog;  real<lower=0> sigma_frac_prior_sdlog;
  real<lower=0> rw_sigma_prior_sd;             // half-normal scale for the RW step sd

  // ---- optional endpoint tweaks (as Model B) -----------------------------
  array[J] int<lower=0, upper=1> use_upper_tweak; vector[J] upper_tweak_mean; vector<lower=0>[J] upper_tweak_sd;
  array[J] int<lower=0, upper=1> use_lower_tweak; vector[J] lower_tweak_mean; vector<lower=0>[J] lower_tweak_sd;
}

parameters {
  real pc0_logit;                              // logit community-death prop at grid pt 1
  vector[G - 1] rw_z;                          // standardised RW innovations
  real<lower=0> rw_sigma;                      // RW step sd (per sqrt(day))

  vector[J] lower_raw;
  vector[J] upper_gap_raw;
  vector<lower=0>[J] sigma_frac;
}

transformed parameters {
  vector[G] pc_logit;
  vector<lower=0, upper=1>[G] p_comm;
  vector<lower=0, upper=1>[G] R_prob;          // reach probability (for the burial term)
  vector<lower=0, upper=1>[G] R_clock;         // normalised reach clock (for scaling)
  vector[J] lower_est;
  vector[J] upper_est;
  vector[J] span_est;

  // Random walk on the logit community-death proportion.
  pc_logit[1] = pc0_logit;
  for (g in 2:G) pc_logit[g] = pc_logit[g - 1] + rw_sigma * sqrt(grid_dt) * rw_z[g - 1];

  for (g in 1:G) {
    p_comm[g]  = inv_logit(pc_logit[g]);
    R_prob[g]  = 1 - p_comm[g];
    R_clock[g] = fmin(1, fmax(0, (c_hi - p_comm[g]) / (c_hi - c_lo)));
  }

  // Endpoint transforms (identical to Model B).
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
  // ---- latent reach (RW) priors -----------------------------------------
  pc0_logit ~ normal(0, 1.5);
  rw_z      ~ std_normal();
  rw_sigma  ~ normal(0, rw_sigma_prior_sd);    // half-normal (rw_sigma constrained > 0)

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

  // ---- anchor likelihood: clock value depends on the parameter's clock ---
  for (n in 1:N) {
    int j = param_id[n];
    real q = (clock_type[j] == 1)
               ? R_clock[anchor_grid[n]]                       // latent reach clock
               : fmin(1, s_grid[anchor_grid[n]] / s_ref);      // supplied SDB clock
    real mu = (increases[j] == 1)
                ? lower_est[j] + span_est[j] * q
                : upper_est[j] - span_est[j] * q;
    y_obs[n] ~ student_t(4, mu, sigma_frac[j] * span_est[j] * obs_sd_mult[n]);
  }
}

generated quantities {
  // Per-endpoint-parameter native-unit curve on the grid (the harness assembles
  // these + the deterministic uf_comm = 1 - R_prob .* s_grid into PARAM_LEVELS order).
  matrix[J, G] theta_pred;
  vector[G] reach_prob   = R_prob;
  vector[G] reach_clock  = R_clock;
  vector[G] comm_death_p = p_comm;
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
