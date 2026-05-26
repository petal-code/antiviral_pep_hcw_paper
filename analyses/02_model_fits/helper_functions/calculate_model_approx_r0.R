## =============================================================================
## R0_single_type_from_args.R
##
## Single-type (genPop-dominant) R0 approximation for a fiber simulation,
## evaluated at the very beginning of the outbreak (t = 0). All time-varying
## inputs (prob_hospitalised_genPop, hospitalisation_delay_factor, the unsafe
## funeral curves, and the prop_etu / ipc_helper inputs that feed into
## hospital_quarantine_efficacy) are read at t = 0.
##
## Usage:
##     source("setup_model_parameters.R")
##     scenario_matrix <- read_scenario_matrix("final_four_scenario_values.csv")
##     mp <- make_model_parameters(
##       scenario_id     = "Worst_WestAfrica",
##       scenario_matrix = scenario_matrix
##     )
##
##     R0_single_type_from_args(mp$args, n = 50000, seed = 1)
## =============================================================================


## Null-coalesce: returns b if a is NULL.
`%||%` <- function(a, b) if (is.null(a)) b else a


## Resolve a parameter that may be a scalar OR a function(t), at t = 0.
.at_t0 <- function(x) {
  if (is.null(x))     return(NULL)
  if (is.function(x)) return(x(0))
  x
}


## Derive hospital_quarantine_efficacy at t = 0 from the prop_etu / ipc_helper /
## etu_efficacy_baseline triplet (matching the updated model's internal logic).
.hospital_quarantine_efficacy_t0 <- function(args) {
  if (!is.null(args$hospital_quarantine_efficacy)) {
    return(.at_t0(args$hospital_quarantine_efficacy))
  }
  prop_etu_0 <- .at_t0(args$prop_etu)
  ipc_0      <- .at_t0(args$ipc_helper)
  etu_b      <- args$etu_efficacy_baseline
  if (is.null(prop_etu_0) || is.null(ipc_0) || is.null(etu_b)) {
    stop(
      "Cannot derive hospital_quarantine_efficacy at t = 0. Supply either ",
      "`hospital_quarantine_efficacy`, or all of `prop_etu`, `ipc_helper`, ",
      "and `etu_efficacy_baseline` in args.",
      call. = FALSE
    )
  }
  etu_efficacy_0 <- etu_b + (1 - etu_b) * ipc_0
  prop_etu_0 * etu_efficacy_0 + (1 - prop_etu_0) * ipc_0
}


## -----------------------------------------------------------------------------
## Monte Carlo Q_g at t = 0, mirroring complete_offspring_info() logic.
## Returns Q_g (expected post-admission GT fraction averaged over parents),
## together with the realised hospitalisation and death-location probabilities.
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

  ## 3. Would-be community outcome (death/recovery, ignoring hospitalisation)
  prob_death_given_symp_comm <- if (args$prob_symptomatic > 0) {
    args$prob_death_comm / args$prob_symptomatic
  } else 0
  would_die_comm <- symptomatic &
    as.logical(rbinom(n, 1, prob_death_given_symp_comm))

  T_comm_out <- T_incub
  if (any(would_die_comm)) {
    T_comm_out[would_die_comm] <- T_comm_out[would_die_comm] +
      args$onset_to_death(sum(would_die_comm))
  }
  if (any(!would_die_comm)) {
    T_comm_out[!would_die_comm] <- T_comm_out[!would_die_comm] +
      args$onset_to_recovery(sum(!would_die_comm))
  }

  ## 4. Potential hospitalisation
  prob_hosp_given_symp <- if (args$prob_symptomatic > 0) {
    prob_hosp_g / args$prob_symptomatic
  } else 0
  potentially_hosp <- symptomatic &
    as.logical(rbinom(n, 1, prob_hosp_given_symp))

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
## Main: single-type R0 at t = 0
## -----------------------------------------------------------------------------
R0_single_type_from_args <- function(args, n = 50000, seed = NULL) {

  ## Resolve time-varying scenario components at t = 0
  p_uf_cgp_0 <- .at_t0(args$p_unsafe_funeral_comm_genPop)
  p_uf_hgp_0 <- .at_t0(args$p_unsafe_funeral_hosp_genPop)
  hq_eff_0   <- .hospital_quarantine_efficacy_t0(args)
  if (is.null(p_uf_cgp_0)) stop("`args$p_unsafe_funeral_comm_genPop` is required.", call. = FALSE)
  if (is.null(p_uf_hgp_0)) stop("`args$p_unsafe_funeral_hosp_genPop` is required.", call. = FALSE)

  ## Q_g and per-parent death-location probabilities at t = 0
  ps <- compute_Q_genPop_from_args(args = args, n = n, seed = seed)

  ## Direct (community + hospital) contribution
  R0_direct <- args$mn_offspring_genPop * (1 - hq_eff_0 * ps$Q_g)

  ## Funeral contributions, split by where the parent died
  R0_funeral_comm <- ps$p_die_comm * args$mn_offspring_funeral *
    (1 - args$safe_funeral_efficacy * (1 - p_uf_cgp_0))
  R0_funeral_hosp <- ps$p_die_hosp * args$mn_offspring_funeral *
    (1 - args$safe_funeral_efficacy * (1 - p_uf_hgp_0))
  R0_funeral <- R0_funeral_comm + R0_funeral_hosp

  list(
    R0                              = R0_direct + R0_funeral,
    R0_direct                       = R0_direct,
    R0_funeral                      = R0_funeral,
    R0_funeral_comm                 = R0_funeral_comm,
    R0_funeral_hosp                 = R0_funeral_hosp,
    Q_g                             = ps$Q_g,
    p_realised_hosp                 = ps$p_realised_hosp,
    p_die_comm                      = ps$p_die_comm,
    p_die_hosp                      = ps$p_die_hosp,
    ## Resolved time-varying inputs at t = 0, for sanity-checking:
    prob_hospitalised_genPop_t0     = .at_t0(args$prob_hospitalised_genPop),
    hospital_quarantine_efficacy_t0 = hq_eff_0,
    p_unsafe_funeral_comm_genPop_t0 = p_uf_cgp_0,
    p_unsafe_funeral_hosp_genPop_t0 = p_uf_hgp_0
  )
}

## =============================================================================
## solve_offspring_means_for_R0.R
##
## Inverse of R0_single_type_from_args(): given a target R0 and a target share
## of transmission via funerals, return the values of mn_offspring_genPop and
## mn_offspring_funeral that produce them under the single-type approximation.
##
## Depends on helpers from R0_single_type_from_args.R:
##   .at_t0, .hospital_quarantine_efficacy_t0, compute_Q_genPop_from_args, %||%
## Source that file first.
##
## Logic. Under the single-type approximation:
##     R0_direct  = mn_offspring_genPop  * D
##     R0_funeral = mn_offspring_funeral * F
## with
##     D = 1 - hospital_quarantine_efficacy(0) * Q_g
##     F = p_die_comm * (1 - sfe * (1 - p_uf_cgp(0)))
##       + p_die_hosp * (1 - sfe * (1 - p_uf_hgp(0)))
## Both D and F depend only on parent attributes (Q_g, death-location
## probabilities) and resolved time-varying inputs — not on the means.
## Hence for a target R0 and funeral share pi:
##     mn_offspring_genPop  = (1 - pi) * R0 / D
##     mn_offspring_funeral =       pi * R0 / F
## =============================================================================

solve_offspring_means_for_R0 <- function(R0,
                                         args,
                                         proportion_transmission_from_funerals,
                                         n    = 50000,
                                         seed = NULL) {

  ## --- Validate inputs ----------------------------------------------------
  if (!is.numeric(R0) || length(R0) != 1L || is.na(R0) || R0 <= 0) {
    stop("`R0` must be a single positive number.", call. = FALSE)
  }
  pi_f <- proportion_transmission_from_funerals
  if (!is.numeric(pi_f) || length(pi_f) != 1L || is.na(pi_f) || pi_f < 0 || pi_f > 1) {
    stop("`proportion_transmission_from_funerals` must be a single number in [0, 1].",
         call. = FALSE)
  }

  ## --- Resolve time-varying inputs at t = 0 -------------------------------
  p_uf_cgp_0 <- .at_t0(args$p_unsafe_funeral_comm_genPop)
  p_uf_hgp_0 <- .at_t0(args$p_unsafe_funeral_hosp_genPop)
  hq_eff_0   <- .hospital_quarantine_efficacy_t0(args)
  if (is.null(p_uf_cgp_0)) stop("`args$p_unsafe_funeral_comm_genPop` is required.", call. = FALSE)
  if (is.null(p_uf_hgp_0)) stop("`args$p_unsafe_funeral_hosp_genPop` is required.", call. = FALSE)

  ## --- Q_g and per-parent death-location probabilities at t = 0 -----------
  ps <- compute_Q_genPop_from_args(args = args, n = n, seed = seed)

  ## --- Compute the D and F multipliers ------------------------------------
  D <- 1 - hq_eff_0 * ps$Q_g
  F <- ps$p_die_comm * (1 - args$safe_funeral_efficacy * (1 - p_uf_cgp_0)) +
    ps$p_die_hosp * (1 - args$safe_funeral_efficacy * (1 - p_uf_hgp_0))

  if (!is.finite(D) || D <= 0) {
    stop(
      "Direct-transmission multiplier D = ", signif(D, 4), " is non-positive. ",
      "This usually means hospital_quarantine_efficacy(0) * Q_g >= 1; the direct ",
      "channel is fully shut and the target R0 cannot be produced through genPop offspring alone.",
      call. = FALSE
    )
  }

  ## --- Target contributions -----------------------------------------------
  R0_direct_target  <- (1 - pi_f) * R0
  R0_funeral_target <-       pi_f * R0

  ## --- Solve for the means ------------------------------------------------
  mn_offspring_genPop_required <- R0_direct_target / D

  if (R0_funeral_target == 0) {
    mn_offspring_funeral_required <- 0
  } else if (!is.finite(F) || F <= 0) {
    stop(
      "Funeral-transmission multiplier F = ", signif(F, 4), " is non-positive ",
      "(under current params, deaths or unsafe-funeral conditions are too rare), ",
      "but a positive funeral share was requested. Reduce ",
      "`proportion_transmission_from_funerals` to 0, or revisit prob_death_*, ",
      "safe_funeral_efficacy, or p_unsafe_funeral_*.",
      call. = FALSE
    )
  } else {
    mn_offspring_funeral_required <- R0_funeral_target / F
  }

  ## --- Output -------------------------------------------------------------
  list(
    mn_offspring_genPop_required    = mn_offspring_genPop_required,
    mn_offspring_funeral_required   = mn_offspring_funeral_required,

    ## Implied R0 decomposition (exact under the formula, MC-noise free):
    R0_target                       = R0,
    R0_direct_target                = R0_direct_target,
    R0_funeral_target               = R0_funeral_target,
    proportion_transmission_from_funerals = pi_f,

    ## Multipliers used in the inversion:
    D_direct_multiplier             = D,
    F_funeral_multiplier            = F,

    ## Underlying MC quantities (for sanity-checking):
    Q_g                             = ps$Q_g,
    p_realised_hosp                 = ps$p_realised_hosp,
    p_die_comm                      = ps$p_die_comm,
    p_die_hosp                      = ps$p_die_hosp,

    ## Resolved time-varying inputs at t = 0:
    hospital_quarantine_efficacy_t0 = hq_eff_0,
    p_unsafe_funeral_comm_genPop_t0 = p_uf_cgp_0,
    p_unsafe_funeral_hosp_genPop_t0 = p_uf_hgp_0
  )
}


## =============================================================================
## Example usage
## =============================================================================
# source("setup_model_parameters.R")
# scenario_matrix <- read_scenario_matrix("final_four_scenario_values.csv")
# mp <- make_model_parameters(
#   scenario_id     = "Worst_WestAfrica",
#   scenario_matrix = scenario_matrix
# )
#
# sol <- solve_offspring_means_for_R0(
#   R0                                    = 1.5,
#   args                                  = mp$args,
#   proportion_transmission_from_funerals = 0.45,
#   n                                     = 50000,
#   seed                                  = 273
# )
# sol$mn_offspring_genPop_required
# sol$mn_offspring_funeral_required
#
# ## Sanity check: plug back in and recompute R0
# args2 <- mp$args
# args2$mn_offspring_genPop  <- sol$mn_offspring_genPop_required
# args2$mn_offspring_funeral <- sol$mn_offspring_funeral_required
# R0_single_type_from_args(args2, n = 50000, seed = 273)$R0      # should be ~1.5

