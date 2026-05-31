## =============================================================================
## calculate_model_approx_r0.R
##
## SHARED single-type (genPop-dominant) R0 approximation at t = 0 for the
## revamped fiber NPI parameterisation, plus the inversion that turns a target
## (R0, prop_funeral) into offspring means. Used by BOTH ABC approaches:
##   * HCW-risk fit: calls solve_offspring_means_for_R0() once (efficacies fixed).
##   * NPI-efficacy fit: caches compute_R0_invariants() once, then recomputes the
##     cheap closed-form D / F per particle as etu_efficacy varies.
##
## Conditioning convention: prob_death_comm / prob_hospitalised_* are used
## DIRECTLY as P(event | symptomatic) (no division by prob_symptomatic), exactly
## matching fiber's complete_offspring_info() (whose Step 0 note removed the old
## prob/prob_symptomatic conversions).
##
## Design notes:
##   * hospital_quarantine_efficacy(0) is the prop_etu(0)-weighted mixture of
##     the fixed scalars etu_efficacy and general_hospital_quarantine_efficacy
##     (no more anchored-floor ipc_helper shape).
##   * Because etu_efficacy is ABC-FITTED, the direct multiplier D depends on
##     the particle. So D (and F) can no longer be precomputed once. Instead we
##     precompute the EFFICACY-INDEPENDENT invariants ONCE
##     (compute_R0_invariants(): Q_g, death-location probs, and the resolved
##     t = 0 prop_etu / unsafe-funeral inputs), then recompute the cheap
##     closed-form D and F per particle from those invariants + the particle's
##     efficacies (D_from_invariants(), F_from_invariants()).
##
## Source setup_model_parameters.R first; fiber must be loaded.
## =============================================================================


## Null-coalesce: returns b if a is NULL.
`%||%` <- function(a, b) if (is.null(a)) b else a


## Resolve a parameter that may be a scalar OR a function(t), at t = 0.
.at_t0 <- function(x) {
  if (is.null(x))     return(NULL)
  if (is.function(x)) return(x(0))
  x
}


## hospital_quarantine_efficacy(0) under the revamped coverage-mixture model.
.hq_eff0_from_parts <- function(prop_etu_0, etu_efficacy,
                                general_hospital_quarantine_efficacy) {
  prop_etu_0 * etu_efficacy +
    (1 - prop_etu_0) * general_hospital_quarantine_efficacy
}


## Derive hospital_quarantine_efficacy at t = 0 from an args list.
.hospital_quarantine_efficacy_t0 <- function(args) {
  if (!is.null(args$hospital_quarantine_efficacy)) {
    return(.at_t0(args$hospital_quarantine_efficacy))
  }
  prop_etu_0 <- .at_t0(args$prop_etu)
  etu_eff    <- args$etu_efficacy
  gen_eff    <- args$general_hospital_quarantine_efficacy
  if (is.null(prop_etu_0) || is.null(etu_eff) || is.null(gen_eff)) {
    stop(
      "Cannot derive hospital_quarantine_efficacy at t = 0. Supply either ",
      "`hospital_quarantine_efficacy`, or all of `prop_etu`, `etu_efficacy`, ",
      "and `general_hospital_quarantine_efficacy` in args.",
      call. = FALSE
    )
  }
  .hq_eff0_from_parts(prop_etu_0, etu_eff, gen_eff)
}


## -----------------------------------------------------------------------------
## Monte Carlo Q_g at t = 0, mirroring complete_offspring_info() logic.
## EFFICACY-INDEPENDENT: depends only on natural history + prob_hosp / delay at
## t = 0. Returns Q_g and the realised hospitalisation / death-location probs.
## -----------------------------------------------------------------------------
compute_Q_genPop_from_args <- function(args, n = 50000, seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  prob_hosp_g <- .at_t0(args$prob_hospitalised_genPop)
  hdf         <- .at_t0(args$hospitalisation_delay_factor) %||% 1.0
  if (is.null(prob_hosp_g)) {
    stop("`args$prob_hospitalised_genPop` is missing.", call. = FALSE)
  }

  ## 1. Incubation period
  T_incub <- args$incubation_period(n)

  ## 2. Symptomatic?
  symptomatic <- as.logical(rbinom(n, 1, args$prob_symptomatic))

  ## 3. Would-be community outcome (death/recovery, ignoring hospitalisation).
  ## Conditional-CFR convention (matches complete_offspring_info()): prob_death_comm
  ## IS P(die | symptomatic, community) directly — no division by prob_symptomatic.
  would_die_comm <- symptomatic &
    as.logical(rbinom(n, 1, args$prob_death_comm))

  T_comm_out <- T_incub
  if (any(would_die_comm)) {
    T_comm_out[would_die_comm] <- T_comm_out[would_die_comm] +
      args$onset_to_death(sum(would_die_comm))
  }
  if (any(!would_die_comm)) {
    T_comm_out[!would_die_comm] <- T_comm_out[!would_die_comm] +
      args$onset_to_recovery(sum(!would_die_comm))
  }

  ## 4. Potential hospitalisation. prob_hosp_g IS P(hospitalised | symptomatic)
  ## directly (conditional-CFR convention; matches complete_offspring_info()).
  potentially_hosp <- symptomatic &
    as.logical(rbinom(n, 1, prob_hosp_g))

  T_hosp <- rep(NA_real_, n)
  if (any(potentially_hosp)) {
    T_hosp[potentially_hosp] <- T_incub[potentially_hosp] +
      args$onset_to_hospitalisation(sum(potentially_hosp)) * hdf
  }

  ## 5. Realised hospitalisation: T_hosp must beat community outcome
  realised_hosp <- potentially_hosp & !is.na(T_hosp) & (T_hosp < T_comm_out)

  ## 6. Outcome status & time for realised-hospitalised parents
  second_chance_death <- if (args$prob_death_comm > 0) {
    args$prob_death_hosp / args$prob_death_comm
  } else 0
  dies_in_hosp <- logical(n)
  idx <- which(realised_hosp & would_die_comm)
  if (length(idx) > 0) {
    dies_in_hosp[idx] <- as.logical(rbinom(length(idx), 1, second_chance_death))
  }
  outcome_death <- would_die_comm
  outcome_death[realised_hosp] <- dies_in_hosp[realised_hosp]

  T_out <- T_comm_out
  idx <- which(realised_hosp & dies_in_hosp)
  if (length(idx) > 0) {
    T_out[idx] <- T_hosp[idx] + args$hospitalisation_to_death(length(idx))
  }
  idx <- which(realised_hosp & !dies_in_hosp)
  if (length(idx) > 0) {
    T_out[idx] <- T_hosp[idx] + args$hospitalisation_to_recovery(length(idx))
  }

  ## 7. Post-admission generation-time mass fraction (0 for non-hospitalised)
  q_h <- numeric(n)
  if (any(realised_hosp)) {
    F_out  <- pgamma(T_out[realised_hosp],
                     shape = args$Tg_shape_genPop, rate = args$Tg_rate_genPop)
    F_hosp <- pgamma(T_hosp[realised_hosp],
                     shape = args$Tg_shape_genPop, rate = args$Tg_rate_genPop)
    valid <- F_out > .Machine$double.eps
    q_h_vals <- numeric(sum(realised_hosp))
    q_h_vals[valid] <- (F_out[valid] - F_hosp[valid]) / F_out[valid]
    q_h[realised_hosp] <- q_h_vals
  }

  list(
    Q_g             = mean(q_h),
    p_realised_hosp = mean(realised_hosp),
    p_die_comm      = mean(outcome_death & !realised_hosp),
    p_die_hosp      = mean(outcome_death &  realised_hosp)
  )
}


## -----------------------------------------------------------------------------
## Efficacy-independent invariants for the R0 inversion (compute ONCE).
## -----------------------------------------------------------------------------
## Bundles the Monte-Carlo quantities (Q_g, death-location probs) with the
## resolved t = 0 scenario inputs (prop_etu(0) and the genPop unsafe-funeral
## probabilities). NONE of these depend on the conditional efficacies, so they
## can be cached per scenario and reused for every ABC particle.
compute_R0_invariants <- function(args, n = 50000, seed = NULL) {
  ps <- compute_Q_genPop_from_args(args = args, n = n, seed = seed)

  prop_etu_0 <- .at_t0(args$prop_etu)
  p_uf_cgp_0 <- .at_t0(args$p_unsafe_funeral_comm_genPop)
  p_uf_hgp_0 <- .at_t0(args$p_unsafe_funeral_hosp_genPop)

  if (is.null(prop_etu_0)) stop("`args$prop_etu` is required.", call. = FALSE)
  if (is.null(p_uf_cgp_0)) stop("`args$p_unsafe_funeral_comm_genPop` is required.", call. = FALSE)
  if (is.null(p_uf_hgp_0)) stop("`args$p_unsafe_funeral_hosp_genPop` is required.", call. = FALSE)

  list(
    Q_g             = ps$Q_g,
    p_realised_hosp = ps$p_realised_hosp,
    p_die_comm      = ps$p_die_comm,
    p_die_hosp      = ps$p_die_hosp,
    prop_etu_0      = prop_etu_0,
    p_uf_cgp_0      = p_uf_cgp_0,
    p_uf_hgp_0      = p_uf_hgp_0
  )
}


## Direct-transmission multiplier D from invariants + this particle's hospital
## quarantine efficacies. D = 1 - hospital_quarantine_efficacy(0) * Q_g.
## (ppe_efficacy does NOT enter D: PPE only thins HCW recipients, and this is
## the genPop-dominant single-type approximation.)
D_from_invariants <- function(inv, etu_efficacy,
                              general_hospital_quarantine_efficacy) {
  hq0 <- .hq_eff0_from_parts(
    inv$prop_etu_0, etu_efficacy, general_hospital_quarantine_efficacy
  )
  1 - hq0 * inv$Q_g
}


## Funeral-transmission multiplier F from invariants + safe-funeral efficacy.
F_from_invariants <- function(inv, safe_funeral_efficacy) {
  inv$p_die_comm * (1 - safe_funeral_efficacy * (1 - inv$p_uf_cgp_0)) +
    inv$p_die_hosp * (1 - safe_funeral_efficacy * (1 - inv$p_uf_hgp_0))
}


## Invert (R0, prop_funeral) -> offspring means given multipliers D, F.
##   R0_direct  = mn_offspring_genPop  * D
##   R0_funeral = mn_offspring_funeral * F
solve_offspring_means <- function(R0, prop_funeral, D, F_fun) {
  if (!is.finite(D) || D <= 0) {
    stop("Direct-transmission multiplier D = ", signif(D, 4),
         " is non-positive; cannot invert for mn_offspring_genPop.",
         call. = FALSE)
  }
  mn_genPop <- (1 - prop_funeral) * R0 / D

  if (prop_funeral == 0) {
    mn_funeral <- 0
  } else if (!is.finite(F_fun) || F_fun <= 0) {
    stop("Funeral-transmission multiplier F = ", signif(F_fun, 4),
         " is non-positive but a positive funeral share was requested.",
         call. = FALSE)
  } else {
    mn_funeral <- prop_funeral * R0 / F_fun
  }

  list(mn_genPop = mn_genPop, mn_funeral = mn_funeral)
}


## -----------------------------------------------------------------------------
## Forward single-type R0 at t = 0 (diagnostic / sanity check)
## -----------------------------------------------------------------------------
R0_single_type_from_args <- function(args, n = 50000, seed = NULL) {

  p_uf_cgp_0 <- .at_t0(args$p_unsafe_funeral_comm_genPop)
  p_uf_hgp_0 <- .at_t0(args$p_unsafe_funeral_hosp_genPop)
  hq_eff_0   <- .hospital_quarantine_efficacy_t0(args)
  if (is.null(p_uf_cgp_0)) stop("`args$p_unsafe_funeral_comm_genPop` is required.", call. = FALSE)
  if (is.null(p_uf_hgp_0)) stop("`args$p_unsafe_funeral_hosp_genPop` is required.", call. = FALSE)

  ps <- compute_Q_genPop_from_args(args = args, n = n, seed = seed)

  R0_direct <- args$mn_offspring_genPop * (1 - hq_eff_0 * ps$Q_g)

  R0_funeral_comm <- ps$p_die_comm * args$mn_offspring_funeral *
    (1 - args$safe_funeral_efficacy * (1 - p_uf_cgp_0))
  R0_funeral_hosp <- ps$p_die_hosp * args$mn_offspring_funeral *
    (1 - args$safe_funeral_efficacy * (1 - p_uf_hgp_0))
  R0_funeral <- R0_funeral_comm + R0_funeral_hosp

  list(
    R0                              = R0_direct + R0_funeral,
    R0_direct                       = R0_direct,
    R0_funeral                      = R0_funeral,
    Q_g                             = ps$Q_g,
    p_realised_hosp                 = ps$p_realised_hosp,
    p_die_comm                      = ps$p_die_comm,
    p_die_hosp                      = ps$p_die_hosp,
    hospital_quarantine_efficacy_t0 = hq_eff_0,
    p_unsafe_funeral_comm_genPop_t0 = p_uf_cgp_0,
    p_unsafe_funeral_hosp_genPop_t0 = p_uf_hgp_0
  )
}


## -----------------------------------------------------------------------------
## Convenience inversion at given efficacies (main-process sanity check).
## -----------------------------------------------------------------------------
## Mirrors the old solve_offspring_means_for_R0() interface: computes the
## invariants from `args`, then D and F at the efficacies already in `args`.
## During ABC the per-particle path uses compute_R0_invariants() +
## D_from_invariants()/F_from_invariants() directly (see
## abc_calibration_functions_npi.R).
solve_offspring_means_for_R0 <- function(R0,
                                         args,
                                         proportion_transmission_from_funerals,
                                         n    = 50000,
                                         seed = NULL) {
  if (!is.numeric(R0) || length(R0) != 1L || is.na(R0) || R0 <= 0) {
    stop("`R0` must be a single positive number.", call. = FALSE)
  }
  pi_f <- proportion_transmission_from_funerals
  if (!is.numeric(pi_f) || length(pi_f) != 1L || is.na(pi_f) || pi_f < 0 || pi_f > 1) {
    stop("`proportion_transmission_from_funerals` must be a single number in [0, 1].",
         call. = FALSE)
  }

  inv <- compute_R0_invariants(args = args, n = n, seed = seed)
  D <- D_from_invariants(inv, args$etu_efficacy,
                         args$general_hospital_quarantine_efficacy)
  F_fun <- F_from_invariants(inv, args$safe_funeral_efficacy)
  means <- solve_offspring_means(R0, pi_f, D, F_fun)

  list(
    mn_offspring_genPop_required    = means$mn_genPop,
    mn_offspring_funeral_required   = means$mn_funeral,
    R0_target                       = R0,
    proportion_transmission_from_funerals = pi_f,
    D_direct_multiplier             = D,
    F_funeral_multiplier            = F_fun,
    invariants                      = inv,
    Q_g                             = inv$Q_g,
    p_die_comm                      = inv$p_die_comm,
    p_die_hosp                      = inv$p_die_hosp,
    hospital_quarantine_efficacy_t0 = .hq_eff0_from_parts(
      inv$prop_etu_0, args$etu_efficacy,
      args$general_hospital_quarantine_efficacy
    )
  )
}
