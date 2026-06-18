## =============================================================================
## calculate_model_approx_rt.R
##
## Time-resolved single-type (genPop-dominant) reproduction-number curves for the
## fiber NPI parameterisation, computed the SAME analytical (Monte-Carlo over
## synthetic parents) way as calculate_model_approx_r0.R -- NO full simulation.
##
## Emits TWO curves over a calendar-time grid, from a SINGLE Monte-Carlo pass:
##
##   * R_inst(t)  -- instantaneous / frozen-conditions reproduction number
##                   (Cori sense): the per-generation multiplier if response
##                   conditions stayed at their time-t values. Every time-varying
##                   curve is resolved AT t. Empirical counterpart in fiber:
##                   compute_reproduction_number(...)$R_instantaneous.
##
##   * R_case(t)  -- cohort / case reproduction number: the expected offspring of
##                   the cohort INFECTED at t, integrating over the conditions
##                   each member actually experiences across its own natural
##                   history. Time-varying curves are resolved at t + (the
##                   event's own delay). Empirical counterpart in fiber:
##                   compute_reproduction_number(...)$R_case.
##
## Both are CONTROL reproduction numbers (susceptible_deplete = FALSE): no S/N
## factor. They coincide only while the response curves are flat over a
## generation interval; otherwise they differ by a forward GI convolution.
##
## Relationship to the t = 0 R0 approximation (sanity property, see tests):
##   With the same args, R_inst(0) reduces to R0_single_type_from_args(args)$R0
##   up to Monte-Carlo / draw-order differences (this file draws FULL-LENGTH
##   natural-history vectors so it can hold them fixed -- common random numbers --
##   across the whole grid; the R0 file draws per-subset, so the two match to MC
##   tolerance, NOT bit-for-bit). R_case(0) does NOT generally equal R0: the day-0
##   cohort already lives through the early ramp.
##
## SELF-CONTAINED. This file defines its own internal helpers and does NOT source
## or reach into calculate_model_approx_r0.R (whose deterministic seed/n inversion
## the ABC pipeline depends on must stay byte-stable). Same input contract as the
## R0 helpers: a single `args` list (merged make_base_args() + the time-varying
## curves from build_time_varying_args()). `fiber` need not be loaded to RUN these
## (the natural-history samplers and curves are already closures inside args), but
## the curves are normally fiber make_time_varying() objects exposing breakpoints
## via attr(curve, "times").
##
## Conditioning convention (preserved from the R0 code): prob_death_comm,
## prob_death_hosp, prob_hospitalised_genPop are P(event | symptomatic) DIRECTLY
## (no division by prob_symptomatic). Upstream constraint
## prob_death_hosp <= prob_death_comm keeps the second-chance ratio in [0, 1].
##
## SCOPE / NON-GOALS: single-type genPop-dominant only (the HCW/PPE pathway is
## collapsed, so ppe_coverage_hcw * ppe_efficacy does NOT enter these curves);
## depletion off; OBV PEP gate off. See the README block at the bottom of the
## file's roxygen for the multi-type / depletion extensions that are out of scope.
## =============================================================================


## Null-coalesce: returns b if a is NULL.
`%||%` <- function(a, b) if (is.null(a)) b else a


## -----------------------------------------------------------------------------
## Internal helpers (own copies; do NOT depend on calculate_model_approx_r0.R)
## -----------------------------------------------------------------------------

## Resolve a parameter that may be a scalar OR a function(t) at time(s) `t`.
## Vectorised over `t`: a curve is called once on the whole length-n time vector;
## a scalar is recycled to length(t). Returns NULL for NULL input.
.resolve_at <- function(x, t) {
  if (is.null(x))     return(NULL)
  if (is.function(x)) return(x(t))
  rep(x, length(t))
}

## Post-admission hospital quarantine efficacy hq at time(s) `t`.
## Mirrors the R0 code: use a supplied hospital_quarantine_efficacy curve if
## present, otherwise derive the prop_etu-weighted mixture of the two fixed
## scalars. Vectorised over `t`.
.hq_at <- function(args, t) {
  if (!is.null(args$hospital_quarantine_efficacy)) {
    return(.resolve_at(args$hospital_quarantine_efficacy, t))
  }
  pe <- .resolve_at(args$prop_etu, t)
  if (is.null(pe)) {
    stop("Cannot resolve hospital_quarantine_efficacy: supply either ",
         "`hospital_quarantine_efficacy` or `prop_etu` (+ etu_efficacy and ",
         "general_hospital_quarantine_efficacy) in args.", call. = FALSE)
  }
  pe * args$etu_efficacy + (1 - pe) * args$general_hospital_quarantine_efficacy
}

## Collect the union of breakpoints (attr(curve, "times")) across the model's
## time-varying inputs, used to build a default evaluation grid. Returns a sorted
## unique numeric vector (possibly empty if every input is scalar).
.collect_breakpoints <- function(args) {
  curve_names <- c(
    "prop_etu", "hospital_quarantine_efficacy",
    "prob_hospitalised_genPop", "hospitalisation_delay_factor",
    "p_unsafe_funeral_comm_genPop", "p_unsafe_funeral_hosp_genPop",
    "mn_offspring_genPop", "mn_offspring_funeral"
  )
  ts <- numeric(0)
  for (nm in curve_names) {
    x <- args[[nm]]
    if (is.function(x)) {
      tt <- attr(x, "times")
      if (!is.null(tt)) ts <- c(ts, as.numeric(tt))
    }
  }
  sort(unique(ts))
}


## -----------------------------------------------------------------------------
## (1) Draw the common-random-number (CRN) natural-history basis ONCE.
## -----------------------------------------------------------------------------
#' Draw the fixed natural-history variate basis for the R(t) Monte Carlo.
#'
#' Draws FULL-LENGTH vectors so every per-parent quantity can be held fixed
#' across the whole time grid and across both modes (common random numbers): the
#' shape of the curve is then driven by parameter movement, not by independent
#' resampling noise at each grid point, and one `seed` makes the entire two-curve
#' grid deterministic. Binary decisions are stored as held uniforms and
#' re-thresholded per t against the (possibly time-varying) probabilities.
#'
#' @param args model argument list (samplers `incubation_period`,
#'   `onset_to_hospitalisation`, `onset_to_death`, `onset_to_recovery`,
#'   `hospitalisation_to_death`, `hospitalisation_to_recovery` must be present).
#' @param n number of synthetic parents.
#' @param seed optional integer seed.
#' @return a named list of length-n variate / uniform vectors.
draw_nh_variates <- function(args, n = 50000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  need <- c("incubation_period", "onset_to_hospitalisation", "onset_to_death",
            "onset_to_recovery", "hospitalisation_to_death",
            "hospitalisation_to_recovery")
  miss <- need[!vapply(need, function(k) is.function(args[[k]]), logical(1))]
  if (length(miss) > 0) {
    stop("args is missing natural-history sampler(s): ",
         paste(miss, collapse = ", "), ".", call. = FALSE)
  }

  list(
    n               = n,
    T_incub         = args$incubation_period(n),
    ## per-parent uniforms for the four binary decisions
    U_symp          = runif(n),
    U_die_comm      = runif(n),
    U_pot_hosp      = runif(n),
    U_second_chance = runif(n),
    ## RAW (unscaled) onset->hospitalisation; scaled per t by the delay factor
    raw_onset_hosp  = args$onset_to_hospitalisation(n),
    ## full-length outcome-timing draws (indexed, never re-drawn)
    onset_death     = args$onset_to_death(n),
    onset_recovery  = args$onset_to_recovery(n),
    hosp_death      = args$hospitalisation_to_death(n),
    hosp_recovery   = args$hospitalisation_to_recovery(n)
  )
}


## -----------------------------------------------------------------------------
## (2) Evaluate the single-type R at ONE time `t`, ONE mode.
## -----------------------------------------------------------------------------
#' Single-type reproduction number at one calendar time, one mode.
#'
#' @param variates output of [draw_nh_variates()].
#' @param args model argument list.
#' @param t scalar calendar time.
#' @param mode "case" (forward: curves at t + event delay) or "instantaneous"
#'   (frozen: curves at t).
#' @param hq_fidelity "admission" (v1, default: hq at the parent's admission time
#'   t + T_hosp, held constant over the post-admission window) or "weighted" (v2:
#'   GT-mass-weighted integral of hq over the post-admission window).
#' @param n_subpoints number of sub-intervals for hq_fidelity = "weighted".
#' @return a named list: R, R_direct, R_funeral, Q_g, p_die_comm, p_die_hosp,
#'   hq_eff, prob_hosp_mean, D.
Rt_single_type_at <- function(variates, args, t,
                              mode        = c("case", "instantaneous"),
                              hq_fidelity = c("admission", "weighted"),
                              n_subpoints = 6L) {
  mode        <- match.arg(mode)
  hq_fidelity <- match.arg(hq_fidelity)

  n       <- variates$n
  T_incub <- variates$T_incub

  ## --- onset clock: where prob_hosp / delay are resolved -------------------
  ## case  -> the parent's own onset, t + T_incub
  ## inst  -> frozen at t
  onset_clock <- if (mode == "case") t + T_incub else rep(t, n)

  ## --- natural history (mirrors compute_Q_genPop_from_args, with CRN) -------
  symptomatic    <- variates$U_symp     < args$prob_symptomatic
  would_die_comm <- symptomatic & (variates$U_die_comm < args$prob_death_comm)

  prob_hosp_vec  <- .resolve_at(args$prob_hospitalised_genPop, onset_clock)
  if (is.null(prob_hosp_vec)) {
    stop("`args$prob_hospitalised_genPop` is missing.", call. = FALSE)
  }
  hdf_vec <- if (is.null(args$hospitalisation_delay_factor)) rep(1.0, n) else
    .resolve_at(args$hospitalisation_delay_factor, onset_clock)

  potentially_hosp <- symptomatic & (variates$U_pot_hosp < prob_hosp_vec)

  T_comm_out <- T_incub + ifelse(would_die_comm,
                                 variates$onset_death, variates$onset_recovery)
  T_hosp     <- T_incub + variates$raw_onset_hosp * hdf_vec
  realised_hosp <- potentially_hosp & (T_hosp < T_comm_out)

  ## second-chance survival once hospitalised (conditional-CFR convention)
  sc <- if (args$prob_death_comm > 0) args$prob_death_hosp / args$prob_death_comm else 0
  dies_in_hosp <- realised_hosp & would_die_comm & (variates$U_second_chance < sc)

  outcome_death <- would_die_comm
  outcome_death[realised_hosp] <- dies_in_hosp[realised_hosp]

  ## outcome (death/recovery) calendar time, relative to infection
  T_out <- T_comm_out
  idx_hd <- realised_hosp & dies_in_hosp
  idx_hr <- realised_hosp & !dies_in_hosp
  T_out[idx_hd] <- T_hosp[idx_hd] + variates$hosp_death[idx_hd]
  T_out[idx_hr] <- T_hosp[idx_hr] + variates$hosp_recovery[idx_hr]

  ## --- post-admission GT-mass fraction q_h ---------------------------------
  q_h <- numeric(n)
  if (any(realised_hosp)) {
    F_out  <- pgamma(T_out[realised_hosp],
                     shape = args$Tg_shape_genPop, rate = args$Tg_rate_genPop)
    F_hosp <- pgamma(T_hosp[realised_hosp],
                     shape = args$Tg_shape_genPop, rate = args$Tg_rate_genPop)
    valid  <- F_out > .Machine$double.eps
    qv <- numeric(sum(realised_hosp))
    qv[valid] <- (F_out[valid] - F_hosp[valid]) / F_out[valid]
    q_h[realised_hosp] <- qv
  }
  Q_g <- mean(q_h)

  ## --- DIRECT term: R_direct = mn_genPop * (1 - mean(hq * q_h)) -------------
  if (hq_fidelity == "admission") {
    ## hq at the parent's admission time, held constant post-admission (v1).
    admission_clock <- if (mode == "case") t + T_hosp else rep(t, n)
    hq_vec  <- .hq_at(args, admission_clock)
    removed <- hq_vec * q_h                      # 0 wherever q_h == 0
  } else {
    ## GT-mass-weighted hq integral over [T_hosp, T_out] (v2).
    removed <- .removed_weighted(args, t, mode, T_hosp, T_out, q_h,
                                 realised_hosp, n_subpoints,
                                 args$Tg_shape_genPop, args$Tg_rate_genPop)
  }
  D        <- 1 - mean(removed)
  mn_g     <- .resolve_at(args$mn_offspring_genPop, t)   # parent infection time = t
  R_direct <- mn_g * D

  ## effective (GT-mass-weighted) hq actually applied -- diagnostic only
  hq_eff <- if (sum(q_h) > 0) sum(removed) / sum(q_h) else NA_real_

  ## --- FUNERAL term: folded into the MC expectation ------------------------
  ## resolve unsafe-funeral probs + funeral mean at the DEATH clock.
  death_clock <- if (mode == "case") t + T_out else rep(t, n)
  p_uf_c <- .resolve_at(args$p_unsafe_funeral_comm_genPop, death_clock)
  p_uf_h <- .resolve_at(args$p_unsafe_funeral_hosp_genPop, death_clock)
  if (is.null(p_uf_c) || is.null(p_uf_h)) {
    stop("`args$p_unsafe_funeral_comm_genPop` and ",
         "`args$p_unsafe_funeral_hosp_genPop` are required.", call. = FALSE)
  }
  mn_f <- .resolve_at(args$mn_offspring_funeral, death_clock)
  sfe  <- args$safe_funeral_efficacy

  die_comm <- outcome_death & !realised_hosp
  die_hosp <- outcome_death &  realised_hosp

  fun_contrib <- numeric(n)
  fun_contrib[die_comm] <- mn_f[die_comm] * (1 - sfe * (1 - p_uf_c[die_comm]))
  fun_contrib[die_hosp] <- mn_f[die_hosp] * (1 - sfe * (1 - p_uf_h[die_hosp]))
  R_funeral <- mean(fun_contrib)

  list(
    R              = R_direct + R_funeral,
    R_direct       = R_direct,
    R_funeral      = R_funeral,
    Q_g            = Q_g,
    p_die_comm     = mean(die_comm),
    p_die_hosp     = mean(die_hosp),
    hq_eff         = hq_eff,
    prob_hosp_mean = mean(prob_hosp_vec),
    D              = D
  )
}


## GT-mass-weighted post-admission hq integral (hq_fidelity = "weighted", v2).
## For each realised-hospitalised parent, replace hq * q_h by
##   ( int_{T_hosp}^{T_out} hq(t + s) dF_Gamma(s) ) / F_Gamma(T_out),
## approximated with `K` equal sub-intervals (midpoint hq, exact GT-mass weights
## from the genPop generation-time CDF). With a flat hq this collapses to
## hq * q_h, so "weighted" and "admission" agree when prop_etu is locally flat.
.removed_weighted <- function(args, t, mode, T_hosp, T_out, q_h,
                              realised_hosp, K, shape, rate) {
  n <- length(T_hosp)
  removed <- numeric(n)
  rh <- which(realised_hosp & q_h > 0)
  if (length(rh) == 0L) return(removed)

  Th <- T_hosp[rh]; To <- T_out[rh]
  Fo <- pgamma(To, shape = shape, rate = rate)
  acc <- numeric(length(rh))
  for (k in seq_len(K)) {
    lo  <- Th + (k - 1) / K * (To - Th)
    hi  <- Th +  k      / K * (To - Th)
    mid <- 0.5 * (lo + hi)
    dF  <- pgamma(hi, shape = shape, rate = rate) -
           pgamma(lo, shape = shape, rate = rate)
    clock <- if (mode == "case") t + mid else rep(t, length(mid))
    acc <- acc + .hq_at(args, clock) * dF
  }
  removed[rh] <- ifelse(Fo > .Machine$double.eps, acc / Fo, 0)
  removed
}


## -----------------------------------------------------------------------------
## (3) Main entry: BOTH curves over a calendar-time grid, one MC pass.
## -----------------------------------------------------------------------------
#' Time-resolved single-type reproduction-number curves (R_inst and R_case).
#'
#' Draws the natural-history basis ONCE ([draw_nh_variates()]) then evaluates,
#' for every grid time and both modes, the single-type analytical reproduction
#' number ([Rt_single_type_at()]). Returns one row per grid time with the two
#' curves (and per-mode diagnostics) side by side.
#'
#' @param args model argument list (merged make_base_args() + the time-varying
#'   curves from build_time_varying_args()).
#' @param times numeric grid of calendar times. If NULL (default), built from the
#'   union of the time-varying curves' breakpoints (attr(curve, "times")),
#'   spanning floor(min)..ceiling(max) in steps of `by`. NOTE: curves resolved
#'   beyond their breakpoint range are clamped to the endpoint value
#'   (make_time_varying uses approxfun(rule = 2)); R_case lookups at t + delay can
#'   reach past the scenario span and will use the endpoint value.
#' @param n number of synthetic parents.
#' @param seed optional integer seed (one seed -> the entire two-curve grid is
#'   deterministic).
#' @param hq_fidelity "admission" (default) or "weighted" (see [Rt_single_type_at()]).
#' @param by grid step in days when `times` is NULL (default 1 = daily). To
#'   overlay against fiber's compute_reproduction_number() match its `bin_width`
#'   (default 7) downstream rather than coarsening this grid.
#' @param n_subpoints sub-intervals for hq_fidelity = "weighted".
#' @param long if TRUE, also/instead return a tidy long frame
#'   (time, mode, quantity, value) for the headline quantities.
#' @return a wide data.frame (default) with columns: time, R_inst,
#'   R_inst_direct, R_inst_funeral, R_case, R_case_direct, R_case_funeral,
#'   Q_g_inst, Q_g_case, hq_inst, hq_case, p_die_comm_case, p_die_hosp_case,
#'   prob_hosp_case. If long = TRUE, a long data.frame instead.
#' @examples
#' \dontrun{
#'   library(fiber)
#'   source("functions/setup_model_parameters.R")
#'   source("functions/calculate_model_approx_rt.R")
#'   sm  <- read_scenario_matrix("data-processed/final_six_scenario_values_original_approach.csv")
#'   mp  <- make_model_parameters("Middle_DRC_ConflictSmoothed_PlusPlus", sm)
#'   rtc <- Rt_curve_single_type(mp$args, n = 50000, seed = 1)
#'   with(rtc, plot(time, R_case, type = "l"));  with(rtc, lines(time, R_inst, lty = 2))
#'   abline(h = 1, col = "grey")
#' }
Rt_curve_single_type <- function(args, times = NULL, n = 50000, seed = NULL,
                                 hq_fidelity = "admission", by = 1,
                                 n_subpoints = 6L, long = FALSE) {

  if (is.null(args$mn_offspring_genPop) || is.null(args$mn_offspring_funeral)) {
    stop("args must carry mn_offspring_genPop and mn_offspring_funeral ",
         "(set by the R0 inversion / build_abc_model_args*()).", call. = FALSE)
  }

  if (is.null(times)) {
    bp <- .collect_breakpoints(args)
    times <- if (length(bp) == 0L) 0 else
      seq(floor(min(bp)), ceiling(max(bp)), by = by)
  }
  times <- sort(unique(as.numeric(times)))

  variates <- draw_nh_variates(args, n = n, seed = seed)

  rows <- lapply(times, function(t) {
    ins <- Rt_single_type_at(variates, args, t, mode = "instantaneous",
                             hq_fidelity = hq_fidelity, n_subpoints = n_subpoints)
    cas <- Rt_single_type_at(variates, args, t, mode = "case",
                             hq_fidelity = hq_fidelity, n_subpoints = n_subpoints)
    data.frame(
      time            = t,
      R_inst          = ins$R,
      R_inst_direct   = ins$R_direct,
      R_inst_funeral  = ins$R_funeral,
      R_case          = cas$R,
      R_case_direct   = cas$R_direct,
      R_case_funeral  = cas$R_funeral,
      Q_g_inst        = ins$Q_g,
      Q_g_case        = cas$Q_g,
      hq_inst         = ins$hq_eff,
      hq_case         = cas$hq_eff,
      p_die_comm_case = cas$p_die_comm,
      p_die_hosp_case = cas$p_die_hosp,
      prob_hosp_case  = cas$prob_hosp_mean
    )
  })
  wide <- do.call(rbind, rows)

  if (!long) return(wide)

  ## tidy long frame for the headline quantities (R, R_direct, R_funeral, Q_g, hq)
  mk <- function(mode, R, Rd, Rf, Qg, hq) {
    data.frame(
      time = wide$time, mode = mode,
      quantity = rep(c("R", "R_direct", "R_funeral", "Q_g", "hq"),
                     each = nrow(wide)),
      value = c(R, Rd, Rf, Qg, hq)
    )
  }
  rbind(
    mk("instantaneous", wide$R_inst, wide$R_inst_direct, wide$R_inst_funeral,
       wide$Q_g_inst, wide$hq_inst),
    mk("case", wide$R_case, wide$R_case_direct, wide$R_case_funeral,
       wide$Q_g_case, wide$hq_case)
  )
}
