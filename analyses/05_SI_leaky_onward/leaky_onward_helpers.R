# leaky_onward_helpers.R
# -----------------------------------------------------------------------------
# Helpers for the "leaky obeldesivir" sensitivity analysis (reviewer response).
# Pure functions only (fiber + base R) so they can be unit-tested in isolation
# and shipped to future/PSOCK workers. No file IO, no plotting here.
#
# THE QUESTION
#   The main model assumes obeldesivir (OBV), when it works, fully prevents a
#   health-worker (HCW) infection -- the treated person is removed from the
#   transmission chain and passes nothing on. A reviewer notes this is optimistic:
#   a real antiviral might reduce, not abolish, onward transmission. This analysis
#   asks what happens to OBV's estimated DEATH impact (the paper's main outcome,
#   with a focus on HCW deaths) if the drug is "leaky" to any degree.
#
# THE MECHANISM (see estimate_leaky_onward)
#   For a finished WITH-OBV run, `out$prevented_completed` is the set of
#   infections OBV prevented (the "treated & effective" people). We forward-
#   simulate those people instead of deleting them, letting a treated person
#   transmit at a fraction r in [0, 1] of normal ("residual transmissibility").
#   r = 0 reproduces the current model's transmission (treated pass nothing on);
#   r = 1 means OBV does nothing to transmission (the no-OBV world). The drug
#   keeps acting downstream at the SAME coverage, so onward HCWs it reaches are
#   treated too. Validity rests on fiber's no-depletion regime (documented in
#   estimate_leaky_onward), under which each prevented infection's subtree is an
#   independent branching process, so this post-hoc replay is unbiased.
#
# DEATHS, AND WHY THEY ARE BUFFERED (the point of the analysis)
#   We keep OBV's protection against DEATH on everyone it treats -- treated
#   people can be (leakily) infected and transmit, but do not die. This is a
#   COUNTING choice on the simulated tree: a case is a death only if it died AND
#   was not treated-effectively (`outcome & !engaged`). So the extra deaths a
#   leaky drug causes come only from the onward cases the drug MISSES (the
#   coverage / efficacy gaps), because the onward HCWs it reaches are kept alive.
#
# WHAT WE COMPUTE, per base run (one posterior draw x one replicate):
#   Two forward simulations of the SAME `prevented_completed`:
#     A1  no-OBV counterfactual : OBV off, r = 1 -> every prevented lineage plays
#         out with no drug at all. Deaths here = D_noOBV = what OBV averts in full
#         (index people + their whole downstream cascade). r-independent.
#     A2  leaky-OBV, per r      : OBV on (drug keeps treating downstream), treated
#         transmit at r, treated kept alive. Deaths that still occur =
#         D_leaky(r) = deaths among the UNTREATED onward cases only.
#   Because fiber has no depletion, the rest of the epidemic (the realised tdf)
#   is unchanged whether we add these lineages back or not, so:
#         deaths averted by OBV at leakiness r  =  D_noOBV - D_leaky(r).
#   This sits ABOVE the number the paper currently reports (index-only,
#   `sum(prevented_completed$outcome)`) for every r, and only erodes down toward
#   it as r -> 1 -- never below. That is the robustness result, for deaths.
#
#   Everything is tracked for TOTAL and HCW deaths; HCW deaths are the focus.
# -----------------------------------------------------------------------------


# Count deaths on a simulated leaky/no-OBV tree returned by
# estimate_leaky_onward(return_trees = TRUE). A tree row is one infection
# with logical `outcome` (TRUE = died in the natural-history model, i.e. the
# would-be death BEFORE OBV's death protection is applied), logical `engaged`
# (TRUE = treated-and-effective, i.e. OBV reached them), character `class`
# ("HCW"/"genPop"), and integer `generation` (seeds = 1, descendants > 1).
#
# Returns a one-row data.frame of death counts split every way we need:
#   *_all       : all would-be deaths (ignores treatment) -- used for the no-OBV
#                 counterfactual, where nobody is protected.
#   *_untreated : would-be deaths among NOT-treated cases -- the deaths that
#                 actually occur once OBV's death protection is applied (leaky arm).
#   *_treated_desc : would-be deaths among treated DESCENDANTS (generation > 1) --
#                 deaths OBV averts downstream by re-treating the onward cases.
# `_hcw` variants restrict to class == "HCW".
count_tree_deaths <- function(tree) {
  zero <- data.frame(
    deaths_all_total = 0L,          deaths_all_hcw = 0L,
    deaths_untreated_total = 0L,    deaths_untreated_hcw = 0L,
    deaths_treated_desc_total = 0L, deaths_treated_desc_hcw = 0L,
    n_infections = 0L, n_descendants = 0L
  )
  if (is.null(tree) || nrow(tree) == 0L) return(zero)

  o    <- tree$outcome %in% TRUE               # would-be death (NA -> FALSE)
  eng  <- tree$engaged %in% TRUE               # treated-and-effective
  hcw  <- !is.na(tree$class) & tree$class == "HCW"
  desc <- !is.na(tree$generation) & tree$generation > 1L

  data.frame(
    deaths_all_total          = sum(o),
    deaths_all_hcw            = sum(o & hcw),
    deaths_untreated_total    = sum(o & !eng),
    deaths_untreated_hcw      = sum(o & !eng & hcw),
    deaths_treated_desc_total = sum(o & eng & desc),
    deaths_treated_desc_hcw   = sum(o & eng & desc & hcw),
    n_infections              = nrow(tree),
    n_descendants             = sum(desc)
  )
}


# Realised deaths in a base run's transmission tree (the sterilizing WITH-OBV
# epidemic). Restricted to realised rows (non-NA infection time). These are the
# deaths the current pipeline already counts in the WITH-OBV arm.
count_tdf_deaths <- function(tdf) {
  tdf <- tdf[!is.na(tdf$time_infection_absolute), , drop = FALSE]
  o   <- tdf$outcome %in% TRUE
  hcw <- !is.na(tdf$class) & tdf$class == "HCW"
  data.frame(tdf_deaths_total = sum(o), tdf_deaths_hcw = sum(o & hcw))
}


# Core per-run computation. `out` is one fiber::branching_process_main() result
# from a WITH-OBV run (must carry `prevented_completed` and `sim_info$params`,
# i.e. the leaky-branch fiber). Returns a tidy long data.frame with columns
# (r, metric, value); r is NA for the r-independent (A1 / baseline) metrics.
#
#   seed          : base RNG seed for this run's forward sims (A1 and each A2 get
#                   deterministic offsets, so the whole thing is reproducible).
#   r_grid        : residual-transmissibility grid for the leaky arm.
#   max_tree_size : per-forward-sim cap (passed to estimate_leaky_onward, which
#                   warns rather than truncating silently). Choose >= the base
#                   run's check_final_size so the no-OBV cascade is not clipped.
compute_leaky_onward_metrics <- function(out,
                                         seed,
                                         r_grid = seq(0, 1, by = 0.1),
                                         max_tree_size = 1e5) {
  pc     <- out$prevented_completed
  params <- out$sim_info$params
  if (is.null(params)) {
    stop("out$sim_info$params is missing -- this run needs the leaky-branch fiber.",
         call. = FALSE)
  }

  emit <- function(r, metric, value) {
    data.frame(r = r, metric = metric, value = as.numeric(value),
               stringsAsFactors = FALSE)
  }

  rows <- list()

  # --- base-run (sterilizing WITH-OBV) realised deaths; r-independent ---------
  td <- count_tdf_deaths(out$tdf)
  rows[[length(rows) + 1L]] <- emit(NA, "tdf_deaths_total", td$tdf_deaths_total)
  rows[[length(rows) + 1L]] <- emit(NA, "tdf_deaths_hcw",   td$tdf_deaths_hcw)

  # Number and would-be deaths of the prevented (treated-effective) index cases:
  # exactly what the current pipeline reports as "averted" (index-only).
  n_prev   <- if (is.null(pc)) 0L else nrow(pc)
  p_total  <- if (is.null(pc) || n_prev == 0L) 0L else sum(pc$outcome %in% TRUE)
  p_hcw    <- if (is.null(pc) || n_prev == 0L) 0L else
                sum((pc$outcome %in% TRUE) & (pc$class == "HCW"))
  rows[[length(rows) + 1L]] <- emit(NA, "n_prevented",                n_prev)
  rows[[length(rows) + 1L]] <- emit(NA, "reported_index_deaths_total", p_total)
  rows[[length(rows) + 1L]] <- emit(NA, "reported_index_deaths_hcw",   p_hcw)

  # Nothing prevented -> nothing to forward-simulate; everything else is zero and
  # the leaky arm is flat at zero. Emit zeros so empty runs drop out cleanly.
  if (n_prev == 0L) {
    for (m in c("no_obv_deaths_total", "no_obv_deaths_hcw",
                "downstream_deaths_averted_total", "downstream_deaths_averted_hcw")) {
      rows[[length(rows) + 1L]] <- emit(NA, m, 0)
    }
    for (r in r_grid) {
      for (m in c("leaky_deaths_accruing_total", "leaky_deaths_accruing_hcw",
                  "leaky_deaths_averted_downstream_total", "leaky_deaths_averted_downstream_hcw",
                  "net_deaths_averted_total", "net_deaths_averted_hcw",
                  "onward_infections")) {
        rows[[length(rows) + 1L]] <- emit(r, m, 0)
      }
    }
    return(do.call(rbind, rows))
  }

  # --- A1: no-OBV counterfactual (OBV off, full transmission) -----------------
  # Deaths here = D_noOBV: the index people AND their entire downstream cascade,
  # all dying at the natural CFR (no drug anywhere). r-independent.
  a1 <- estimate_leaky_onward(
    prevented_completed        = pc,
    params                     = utils::modifyList(params, list(obv_pep_enabled = FALSE)),
    residual_transmissibility  = 1,
    n_replicates               = 1L,
    max_tree_size              = max_tree_size,
    seed                       = seed,
    return_trees               = TRUE
  )
  d_noobv <- count_tree_deaths(a1$trees[[1]])
  # D_noOBV counts every would-be death (no protection): use *_all.
  no_obv_total <- d_noobv$deaths_all_total
  no_obv_hcw   <- d_noobv$deaths_all_hcw
  rows[[length(rows) + 1L]] <- emit(NA, "no_obv_deaths_total", no_obv_total)
  rows[[length(rows) + 1L]] <- emit(NA, "no_obv_deaths_hcw",   no_obv_hcw)
  # Downstream-only averted = full averted minus index-only reported.
  rows[[length(rows) + 1L]] <- emit(NA, "downstream_deaths_averted_total", no_obv_total - p_total)
  rows[[length(rows) + 1L]] <- emit(NA, "downstream_deaths_averted_hcw",   no_obv_hcw   - p_hcw)

  # --- A2: leaky-OBV per r (OBV on, treated kept alive) -----------------------
  for (ri in seq_along(r_grid)) {
    r <- r_grid[ri]
    a2 <- estimate_leaky_onward(
      prevented_completed       = pc,
      params                    = params,                 # OBV stays on downstream
      residual_transmissibility = r,
      n_replicates              = 1L,
      max_tree_size             = max_tree_size,
      seed                      = seed + 1000L * ri,       # distinct, reproducible
      return_trees              = TRUE
    )
    d <- count_tree_deaths(a2$trees[[1]])

    # Deaths that STILL occur under a leaky drug = untreated would-be deaths
    # (treated people, index and downstream, are kept alive).
    accruing_total <- d$deaths_untreated_total
    accruing_hcw   <- d$deaths_untreated_hcw

    # Deaths OBV averts by re-treating the onward cases (the buffer): would-be
    # deaths among treated DESCENDANTS.
    averted_ds_total <- d$deaths_treated_desc_total
    averted_ds_hcw   <- d$deaths_treated_desc_hcw

    # Deaths averted by OBV at this leakiness = no-OBV deaths minus the deaths
    # that still occur. Stays >= the index-only reported number for all r.
    net_total <- no_obv_total - accruing_total
    net_hcw   <- no_obv_hcw   - accruing_hcw

    rows[[length(rows) + 1L]] <- emit(r, "leaky_deaths_accruing_total",           accruing_total)
    rows[[length(rows) + 1L]] <- emit(r, "leaky_deaths_accruing_hcw",             accruing_hcw)
    rows[[length(rows) + 1L]] <- emit(r, "leaky_deaths_averted_downstream_total", averted_ds_total)
    rows[[length(rows) + 1L]] <- emit(r, "leaky_deaths_averted_downstream_hcw",   averted_ds_hcw)
    rows[[length(rows) + 1L]] <- emit(r, "net_deaths_averted_total",              net_total)
    rows[[length(rows) + 1L]] <- emit(r, "net_deaths_averted_hcw",                net_hcw)
    rows[[length(rows) + 1L]] <- emit(r, "onward_infections",                     d$n_descendants)
  }

  do.call(rbind, rows)
}


# Aggregate the per-run long frame across the ensemble, in the order the user
# specified: FIRST take the median over the stochastic replicates within each
# posterior draw, THEN summarise across draws (median + 95% interval). Works on a
# data.frame with columns (scenario, particle, rep, r, metric, value); `r` may be
# NA for r-independent metrics (kept as its own group via a stable label).
summarise_leaky_onward <- function(per_run) {
  # r as a grouping label that keeps NA distinct (base R aggregate drops NA keys).
  per_run$r_lab <- ifelse(is.na(per_run$r), "NA", format(per_run$r))

  # Step 1: median over replicates within (scenario, particle, r, metric).
  by_draw <- stats::aggregate(
    value ~ scenario + particle + r_lab + metric,
    data = per_run, FUN = median
  )

  # Step 2: across draws -> median + central 95% interval.
  q <- function(x, p) stats::quantile(x, probs = p, names = FALSE, type = 7)
  agg <- do.call(rbind, lapply(
    split(by_draw, list(by_draw$scenario, by_draw$r_lab, by_draw$metric), drop = TRUE),
    function(g) {
      data.frame(
        scenario = g$scenario[1],
        r        = if (g$r_lab[1] == "NA") NA_real_ else as.numeric(g$r_lab[1]),
        metric   = g$metric[1],
        n_draws  = nrow(g),
        median   = median(g$value),
        lo95     = q(g$value, 0.025),
        hi95     = q(g$value, 0.975),
        mean     = mean(g$value),
        stringsAsFactors = FALSE
      )
    }
  ))
  rownames(agg) <- NULL
  agg[order(agg$scenario, agg$metric, agg$r), ]
}
