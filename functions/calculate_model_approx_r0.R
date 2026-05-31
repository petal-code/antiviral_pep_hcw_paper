## =============================================================================
## calculate_model_approx_R0.R
##
## Single-type (genPop-dominant) R0 approximation for a fiber simulation,
## evaluated at the very beginning of the outbreak (t = 0). All time-varying
## inputs (prob_hospitalised_genPop, hospitalisation_delay_factor, the unsafe
## funeral curves, and prop_etu — the only time-varying input that feeds
## hospital_quarantine_efficacy in the new engine) are read at t = 0.
##
## ----------------------------------------------------------------------------
## CHANGED: hospital_quarantine_efficacy(t) formulation
##
## Old (pre-refactor):
##     etu_eff(t)       = etu_efficacy_baseline + (1 - etu_efficacy_baseline) * ipc_helper(t)
##     hq_eff(t)        = prop_etu(t) * etu_eff(t) + (1 - prop_etu(t)) * ipc_helper(t)
##
## New (this file matches the current fiber engine):
##     hq_eff(t)        = prop_etu(t) * etu_efficacy
##                      + (1 - prop_etu(t)) * general_hospital_quarantine_efficacy
##
## Both efficacies are fixed scalars in [0, 1]; only prop_etu(t) is time-varying.
## There is no longer an anchored-floor / IPC-lift on the ETU efficacy, and
## ipc_helper does not enter hq at all.
##
## ----------------------------------------------------------------------------
## R0 math (unchanged; PPE is omitted because the single-type approximation
## treats every offspring of a genPop parent as genPop, so receiver PPE on
## HCWs does not enter):
##
##     R0 ≈ mn_offspring_genPop  · (1 - hq_eff(0) * Q_g)
##        + mn_offspring_funeral · [ p_die_comm * (1 - sfe * (1 - p_uf_cgp(0)))
##                                 + p_die_hosp * (1 - sfe * (1 - p_uf_hgp(0))) ]
##
## where Q_g is the expected post-admission GT mass fraction averaged over
## genPop parents (non-hospitalised parents contribute 0), estimated by Monte
## Carlo using the same per-parent logic as complete_offspring_info().
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


## Derive hospital_quarantine_efficacy at t = 0 from the new triplet:
##   prop_etu(t) (may be time-varying), etu_efficacy (scalar), and
##   general_hospital_quarantine_efficacy (scalar).
##
## A direct `hospital_quarantine_efficacy` field (scalar or function) is still
## honoured if supplied, as a back-door for testing specific hq values without
## going through the prop_etu / efficacy decomposition. The production fiber
## engine does not provide one and instead always derives from the triplet.
.hospital_quarantine_efficacy_t0 <- function(args) {
  if (!is.null(args$hospital_quarantine_efficacy)) {
    return(.at_t0(args$hospital_quarantine_efficacy))
  }
  prop_etu_0  <- .at_t0(args$prop_etu)
  etu_eff     <- args$etu_efficacy                          # fixed scalar
  general_eff <- args$general_hospital_quarantine_efficacy  # fixed scalar
  if (is.null(prop_etu_0) || is.null(etu_eff) || is.null(general_eff)) {
    stop(
      "Cannot derive hospital_quarantine_efficacy at t = 0. Supply either ",
      "`hospital_quarantine_efficacy` directly, or all of `prop_etu`, ",
      "`etu_efficacy`, and `general_hospital_quarantine_efficacy` in args.",
      call. = FALSE
    )
  }
  prop_etu_0 * etu_eff + (1 - prop_etu_0) * general_eff
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
  
  ## 3. Would-be community outcome (death/recovery, ignoring hospitalisation).
  ## Under the conditional-CFR ("Option B") interpretation used by the new
  ## engine, prob_death_comm is P(die | symp, community) directly.
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
  
  ## 4. Potential hospitalisation (conditional on symptomatic). prob_hosp_g
  ## is P(hosp | symp) directly under the new conditional-CFR convention.
  potentially_hosp <- symptomatic &
    as.logical(rbinom(n, 1, prob_hosp_g))
  
  T_hosp <- rep(NA_real_, n)
  if (any(potentially_hosp)) {
    T_hosp[potentially_hosp] <- T_incub[potentially_hosp] +
      args$onset_to_hospitalisation(sum(potentially_hosp)) * hdf
  }
  
  ## 5. Realised hospitalisation: T_hosp must beat community outcome
  realised_hosp <- potentially_hosp & !is.na(T_hosp) & (T_hosp < T_comm_out)
  
  ## 6. Outcome status & time for realised-hospitalised parents.
  ## Second-chance survival uses prob_death_hosp / prob_death_comm so the
  ## realised hospital CFR equals prob_death_hosp (the same identity the new
  ## complete_offspring_info() relies on; sanity-checked upstream by fiber).
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
## solve_offspring_means_for_R0
##
## Inverse of R0_single_type_from_args(): given a target R0 and a target share
## of transmission via funerals, return the (mn_offspring_genPop,
## mn_offspring_funeral) that produce them. D and F multipliers depend only on
## parent attributes and t=0 inputs (not on the means), so this is a direct
## algebraic solve — no root-finding needed.
##
##     mn_offspring_genPop  = (1 - pi) * R0 / D
##     mn_offspring_funeral =       pi * R0 / F
## with
##     D = 1 - hq_eff(0) * Q_g
##     F = p_die_comm * (1 - sfe * (1 - p_uf_cgp(0)))
##       + p_die_hosp * (1 - sfe * (1 - p_uf_hgp(0)))
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
  
  list(
    mn_offspring_genPop_required          = mn_offspring_genPop_required,
    mn_offspring_funeral_required         = mn_offspring_funeral_required,
    
    R0_target                             = R0,
    R0_direct_target                      = R0_direct_target,
    R0_funeral_target                     = R0_funeral_target,
    proportion_transmission_from_funerals = pi_f,
    
    D_direct_multiplier                   = D,
    F_funeral_multiplier                  = F,
    
    Q_g                                   = ps$Q_g,
    p_realised_hosp                       = ps$p_realised_hosp,
    p_die_comm                            = ps$p_die_comm,
    p_die_hosp                            = ps$p_die_hosp,
    
    hospital_quarantine_efficacy_t0       = hq_eff_0,
    p_unsafe_funeral_comm_genPop_t0       = p_uf_cgp_0,
    p_unsafe_funeral_hosp_genPop_t0       = p_uf_hgp_0
  )
}