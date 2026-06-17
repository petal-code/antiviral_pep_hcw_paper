// ============================================================================
// ODV NHP delay-to-initiation efficacy: FULL Bayesian piecewise-exponential model
// ----------------------------------------------------------------------------
// This is the Stan counterpart to the profiled-likelihood fit in
// analyses/odv_nhp_delay_efficacy/01_fit_odv_delay_efficacy.R.
//
// It implements the SAME survival model, but "completely": the K piecewise-
// constant baseline interval hazards are estimated as free parameters jointly
// with the two delay-efficacy shape parameters, instead of being profiled out
// analytically and summarised with a 2-parameter Laplace approximation. The
// posterior is sampled by HMC/NUTS, so baseline-hazard uncertainty, parameter
// correlation, boundary effects and posterior skew are all propagated into the
// fitted efficacy curve.
//
// LIKELIHOOD
//   The animal-level survival data are split into (animal x interval) rows with
//   a time-at-risk `exposure` and a death indicator `event`. For a piecewise-
//   constant hazard this is exactly a Poisson likelihood on the split rows:
//       event[n] ~ Poisson( exposure[n] * hr[n] * lambda[interval[n]] )
//   where lambda[k] is the baseline hazard in interval k and hr[n] is the
//   delay-specific hazard ratio (= 1 - efficacy) for post-initiation rows in a
//   treated arm, and 1 for vehicle / pre-initiation rows. (This equals the
//   piecewise-exponential survival likelihood up to a parameter-free constant.)
//
// DELAY-EFFICACY CURVE (identical to script 01)
//   efficacy(d) = E0 * { g(d) - g(d_zero) } / { g(0) - g(d_zero) }
//   g(d)        = 1 / (1 + exp(k * (d - d50)))  =  inv_logit(k * (d50 - d))
//   - E0     : efficacy at d = 0 (immediate treatment), on the hazard scale
//   - d50    : location of the logistic decline
//   - k      : steepness (FIXED at k_fixed, or estimated if fit_k == 1)
//   - d_zero : day by which efficacy is constrained to zero
//   The scaling forces efficacy(0) = E0 and efficacy(d_zero) = 0.
//
// PRIORS
//   - E0           : flat over [0,1] (no sampling statement => uniform), as in
//                    the bounds-only optimisation of script 01.
//   - d50          : the decline LOCATION is not identified by initiation days
//                    1-4 alone (the likelihood only sees efficacy at the
//                    observed days, and for d50 beyond ~5-6 the curve is flat
//                    across them). It is therefore given a weakly-informative
//                    normal prior centred in the observed window and capped at
//                    d_zero. Set d50_prior_sd very large to recover a flat prior.
//   - lambda[k]    : by default a flat (improper) prior on [0, inf). Stan adds
//                    the constraint's Jacobian automatically, so omitting a
//                    sampling statement gives a uniform prior on the constrained
//                    (>= 0) scale -- the most uninformative choice. The posterior
//                    is still proper because every retained interval has positive
//                    total exposure. Set use_hazard_prior = 1 for a very diffuse
//                    but guaranteed-proper exponential prior instead.
//   - k (if fitted): weakly-informative lognormal centred on the fixed default.
// ============================================================================

data {
  // ---- Sizes -------------------------------------------------------------
  int<lower=1> N;                                 // split (animal x interval) rows
  int<lower=1> K;                                 // number of baseline-hazard intervals
  int<lower=1> D;                                 // number of observed ODV initiation days

  // ---- Split survival data (one row per animal per risk interval) --------
  array[N] int<lower=1, upper=K> interval_idx;    // baseline-hazard interval for the row
  array[N] int<lower=0>          event;           // death indicator (0/1) in the interval
  vector<lower=0>[N]             exposure;        // time at risk in the interval
  array[N] int<lower=0, upper=1> post_treatment;  // 1 if the row is post-initiation (treated arm)
  array[N] int<lower=1, upper=D> dpc_idx;         // index into dpc_obs (set to 1 when unused)

  // ---- Curve settings ----------------------------------------------------
  vector<lower=0>[D] dpc_obs;                     // observed initiation days (e.g. 1,2,3,4)
  real<lower=0>      d_zero;                       // day by which efficacy is forced to 0
  real<lower=0>      eps_hr;                        // numerical floor on the hazard ratio

  // ---- d50 prior (identifies the otherwise non-identified decline location) --
  real<lower=0> d50_prior_mean;                    // centre of the weakly-informative d50 prior
  real<lower=0> d50_prior_sd;                      // sd (set very large to recover a flat prior)

  // ---- k: fixed or fitted ------------------------------------------------
  int<lower=0, upper=1> fit_k;                     // 0 = fix k at k_fixed, 1 = estimate k
  real<lower=0>         k_fixed;                   // value used when fit_k == 0
  real                  k_prior_logmean;           // lognormal prior (used when fit_k == 1)
  real<lower=0>         k_prior_logsd;

  // ---- Baseline-hazard prior switch --------------------------------------
  int<lower=0, upper=1> use_hazard_prior;          // 0 = flat improper, 1 = diffuse exponential
  real<lower=0>         hazard_prior_rate;         // exponential rate (used when == 1)

  // ---- Efficacy curve grid for posterior summaries -----------------------
  int<lower=1>       n_grid;
  vector<lower=0>[n_grid] grid_dpc;
}

parameters {
  real<lower=0, upper=1>      E0;                  // efficacy at day 0 (hazard scale)
  real<lower=0, upper=d_zero> d50;                 // location of the logistic decline (capped at d_zero)
  vector<lower=0>[K]      lambda;                  // baseline interval hazards (jointly estimated)
  array[fit_k] real<lower=0> k_param;              // length 1 if estimating k, else length 0
}

transformed parameters {
  real k;                                          // steepness: fitted or fixed
  vector[D] eff;                                   // efficacy at each observed dpc

  // Resolve the steepness once: fitted value if requested, otherwise the fixed one.
  if (fit_k == 1) {
    k = k_param[1];
  } else {
    k = k_fixed;
  }

  // Efficacy at each observed initiation day, on the hazard scale.
  {
    real g0 = inv_logit(k * (d50 - 0.0));         // g(0)
    real gz = inv_logit(k * (d50 - d_zero));      // g(d_zero)
    for (j in 1:D) {
      real gd = inv_logit(k * (d50 - dpc_obs[j]));
      eff[j] = E0 * (gd - gz) / (g0 - gz);
    }
  }
}

model {
  // ---- Priors ------------------------------------------------------------
  // E0: flat over [0,1] (no statement => uniform on support).
  // d50: weakly-informative normal (truncated to [0, d_zero] by the bounds),
  //      because the decline location is not identified by initiation days 1-4.
  // lambda: flat improper by default; optional diffuse proper alternative.
  d50 ~ normal(d50_prior_mean, d50_prior_sd);
  if (use_hazard_prior == 1) {
    lambda ~ exponential(hazard_prior_rate);
  }
  if (fit_k == 1) {
    k_param[1] ~ lognormal(k_prior_logmean, k_prior_logsd);
  }

  // ---- Likelihood: piecewise-exponential == Poisson on split rows --------
  {
    vector[N] mu;
    for (n in 1:N) {
      real hr = 1;                                 // vehicle / pre-treatment rows
      if (post_treatment[n] == 1) {                // treated, post-initiation rows
        hr = fmax(eps_hr, 1 - eff[dpc_idx[n]]);    // hazard ratio = 1 - efficacy
      }
      mu[n] = exposure[n] * hr * lambda[interval_idx[n]];
    }
    event ~ poisson(mu);
  }
}

generated quantities {
  // Posterior efficacy curve on the supplied grid (forced to 0 at/after d_zero).
  // This replaces the Laplace / MVN-draw uncertainty step of the R script: take
  // pointwise quantiles of efficacy_grid across draws for the credible ribbon.
  vector[n_grid] efficacy_grid;

  // Per-row log-likelihood (for LOO / model comparison, matching repo convention).
  vector[N] log_lik;

  {
    real g0 = inv_logit(k * (d50 - 0.0));
    real gz = inv_logit(k * (d50 - d_zero));
    for (g in 1:n_grid) {
      real e;
      if (grid_dpc[g] >= d_zero) {
        e = 0;
      } else {
        real gd = inv_logit(k * (d50 - grid_dpc[g]));
        e = E0 * (gd - gz) / (g0 - gz);
      }
      efficacy_grid[g] = fmin(1, fmax(0, e));       // keep within [0,1]
    }
  }

  for (n in 1:N) {
    real hr = 1;
    real mu_n;
    if (post_treatment[n] == 1) {
      hr = fmax(eps_hr, 1 - eff[dpc_idx[n]]);
    }
    mu_n = exposure[n] * hr * lambda[interval_idx[n]];
    log_lik[n] = poisson_lpmf(event[n] | mu_n);
  }
}
