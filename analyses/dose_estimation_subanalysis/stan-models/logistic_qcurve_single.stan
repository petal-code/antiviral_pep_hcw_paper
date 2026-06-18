// ============================================================================
// logistic_qcurve_single.stan
// ----------------------------------------------------------------------------
// A SINGLE-curve logistic "Q-curve" fit for the dose-estimation subanalysis.
//
// SAME FAMILY OF MODEL as Model A (the West Africa
// modelA_partialpool_estimateQ_*.stan): a logistic response-quality curve
// scaled between a lower endpoint L and an upper endpoint U,
//
//     mu(day) = L + (U - L) * inv_logit( k * (day - t50) )
//
// fit with a robust Student-t(4) observation model. It differs from Model A in
// three deliberate ways, all driven by what this subanalysis needs:
//
//   1. ONE curve only. There is a single response-quality series here (the dose
//      coverage proportions), so there is no shape to partially pool -- the
//      hierarchy in Model A collapses to a single (t50, k).
//
//   2. RAW relative-day logistic, NOT renormalised to the observation window.
//      Model A uses q_norm(), which rescales the logistic so it hits exactly 0
//      at tau = 0 and exactly 1 at tau = 1 (the edges of the observed window).
//      That is the right choice when the magnitude is carried by literature
//      endpoints and you never look outside the window. Here we explicitly
//      PROJECT THE CURVE FORWARD (up to a year past the data), so we keep the
//      logistic in raw day units: the estimated growth rate k then governs the
//      extrapolation directly, exactly "assuming the same growth rate as
//      estimated".
//
//   3. VERY informative endpoint priors pinning L ~ 0 and U ~ 1 (i.e. the curve
//      runs from a 0% floor to a 100% ceiling). Because the endpoints are
//      bounded to [0, 1] with prior means at the boundaries, these act as tight
//      half-normals that hold the asymptotes at 0% / 100% while the data inform
//      only the SHAPE (midpoint t50 and steepness k). This encodes the modelling
//      assumption that dose coverage ultimately reaches 100% and the observed
//      values are partway up that 0 -> 100% curve.
// ============================================================================

functions {
  // Logistic curve in raw day units, scaled between endpoints L and U.
  real logistic_curve(real day, real L, real U, real t50, real k) {
    return L + (U - L) * inv_logit(k * (day - t50));
  }
}

data {
  // ---- Observations -------------------------------------------------------
  int<lower=1> N;                       // number of (day, proportion) points
  vector[N] day;                        // relative day of each observation
  vector<lower=0, upper=1>[N] y;        // observed proportion in [0, 1]

  // ---- VERY informative endpoint priors (the curve's min / max level) -----
  real lower_prior_mean;                // prior centre for L (default 0)
  real<lower=0> lower_prior_sd;         // tight  -> L pinned near 0
  real upper_prior_mean;                // prior centre for U (default 1)
  real<lower=0> upper_prior_sd;         // tight  -> U pinned near 1

  // ---- Weakly informative shape priors ------------------------------------
  real t50_prior_mean;                  // prior centre for the midpoint (days)
  real<lower=0> t50_prior_sd;
  real logk_prior_mean;                 // prior centre for log(growth rate)
  real<lower=0> logk_prior_sd;

  // ---- Observation-noise prior (half-normal scale) ------------------------
  real<lower=0> sigma_prior_sd;

  // ---- Prediction grid (incl. the forward extrapolation days) -------------
  int<lower=1> N_pred;
  vector[N_pred] day_pred;
}

parameters {
  real<lower=0, upper=1> L;             // lower endpoint (min level)
  real<lower=0, upper=1> U;             // upper endpoint (max level)
  real t50;                             // logistic midpoint (days)
  real log_k;                           // log growth rate (k > 0)
  real<lower=1e-3> sigma;               // observation noise scale; floored a hair
                                        // above 0 so it can never underflow to an
                                        // exact 0 student_t scale during warmup
}

transformed parameters {
  real<lower=0> k = exp(log_k);
}

model {
  // Endpoint priors. With the [0,1] bounds and prior means at 0 / 1, these are
  // tight half-normals that hold the asymptotes at the 0% floor and 100% ceiling.
  L     ~ normal(lower_prior_mean, lower_prior_sd);
  U     ~ normal(upper_prior_mean, upper_prior_sd);

  // Shape priors (weakly informative; the data drive these).
  t50   ~ normal(t50_prior_mean, t50_prior_sd);
  log_k ~ normal(logk_prior_mean, logk_prior_sd);

  // Observation scale (half-normal via the <lower=0> constraint).
  sigma ~ normal(0, sigma_prior_sd);

  // Robust Student-t(4) likelihood: tolerant of the slight non-monotonicity in
  // the data (the last point can dip below the previous one) and of the tension
  // between the observed plateau and the pinned U = 1 ceiling.
  for (n in 1:N)
    y[n] ~ student_t(4, logistic_curve(day[n], L, U, t50, k), sigma);
}

generated quantities {
  // The fitted / extrapolated Q curve on the (possibly forward-projected) grid.
  vector[N_pred] q_pred;
  vector[N] log_lik;                    // pointwise log-lik (LOO / model checks)

  for (m in 1:N_pred)
    q_pred[m] = logistic_curve(day_pred[m], L, U, t50, k);

  for (n in 1:N)
    log_lik[n] = student_t_lpdf(y[n] | 4, logistic_curve(day[n], L, U, t50, k), sigma);
}
