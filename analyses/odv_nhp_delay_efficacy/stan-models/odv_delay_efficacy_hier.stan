// ============================================================================
// ODV NHP delay efficacy -- HIERARCHICAL (partial-pooling) baseline hazards
// ----------------------------------------------------------------------------
// Variant of odv_delay_efficacy.stan. The base model gives each of the K
// baseline interval hazards its own flat (improper) prior; marginalising K such
// independent, diffusely-priored nuisance hazards biases the efficacy upward and
// makes it over-precise (the incidental-parameters / Neyman-Scott problem).
//
// Here the K log-hazards are PARTIALLY POOLED:
//     log_lambda_k ~ Normal(mu_log, tau_log)            (non-centred)
// with hyperpriors on mu_log (overall level) and tau_log (between-interval sd).
// Pooling shrinks the sparse intervals toward a common level, which cuts the
// effective number of free nuisance parameters and hence the bias -- an
// alternative to simply coarsening the baseline (script 04).
//
// Everything else -- the scaled-logistic delay-efficacy curve, the Poisson /
// piecewise-exponential likelihood, the weakly-informative d50 prior and cap,
// and the optional steepness k -- is identical to the base model.
// ============================================================================

data {
  int<lower=1> N;
  int<lower=1> K;
  int<lower=1> D;

  array[N] int<lower=1, upper=K> interval_idx;
  array[N] int<lower=0>          event;
  vector<lower=0>[N]             exposure;
  array[N] int<lower=0, upper=1> post_treatment;
  array[N] int<lower=1, upper=D> dpc_idx;

  vector<lower=0>[D] dpc_obs;
  real<lower=0>      d_zero;
  real<lower=0>      eps_hr;

  real<lower=0> d50_prior_mean;
  real<lower=0> d50_prior_sd;

  int<lower=0, upper=1> fit_k;
  real<lower=0>         k_fixed;
  real                  k_prior_logmean;
  real<lower=0>         k_prior_logsd;

  // ---- Hierarchical hyperpriors on the log baseline hazards ----
  real          mu_prior_mean;     // prior mean for mu_log (overall log-hazard)
  real<lower=0> mu_prior_sd;
  real<lower=0> tau_prior_sd;      // half-normal scale for tau_log (pooling sd)

  int<lower=1>            n_grid;
  vector<lower=0>[n_grid] grid_dpc;
}

parameters {
  real<lower=0, upper=1>      E0;
  real<lower=0, upper=d_zero> d50;
  array[fit_k] real<lower=0>  k_param;

  real          mu_log;            // mean log baseline hazard
  real<lower=0> tau_log;           // between-interval sd of the log hazard
  vector[K]     z;                 // non-centred standard-normal offsets
}

transformed parameters {
  real k;
  vector[D] eff;
  vector[K] log_lambda = mu_log + tau_log * z;     // partially pooled log-hazards
  vector[K] lambda     = exp(log_lambda);

  if (fit_k == 1) {
    k = k_param[1];
  } else {
    k = k_fixed;
  }

  {
    real g0 = inv_logit(k * (d50 - 0.0));
    real gz = inv_logit(k * (d50 - d_zero));
    for (j in 1:D) {
      real gd = inv_logit(k * (d50 - dpc_obs[j]));
      eff[j] = E0 * (gd - gz) / (g0 - gz);
    }
  }
}

model {
  // ---- Hierarchical prior on the baseline log-hazards (non-centred) ----
  z       ~ std_normal();                        // => log_lambda ~ N(mu_log, tau_log)
  mu_log  ~ normal(mu_prior_mean, mu_prior_sd);
  tau_log ~ normal(0, tau_prior_sd);             // half-normal (tau_log >= 0)

  // ---- Curve priors (E0 flat over [0,1]; d50 weakly-informative + capped) ----
  d50 ~ normal(d50_prior_mean, d50_prior_sd);
  if (fit_k == 1) {
    k_param[1] ~ lognormal(k_prior_logmean, k_prior_logsd);
  }

  // ---- Likelihood: piecewise-exponential == Poisson on split rows ----
  {
    vector[N] mu;
    for (n in 1:N) {
      real hr = 1;
      if (post_treatment[n] == 1) {
        hr = fmax(eps_hr, 1 - eff[dpc_idx[n]]);
      }
      mu[n] = exposure[n] * hr * lambda[interval_idx[n]];
    }
    event ~ poisson(mu);
  }
}

generated quantities {
  vector[n_grid] efficacy_grid;
  vector[N]      log_lik;

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
      efficacy_grid[g] = fmin(1, fmax(0, e));
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
