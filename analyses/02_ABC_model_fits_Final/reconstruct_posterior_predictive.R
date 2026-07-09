# reconstruct_posterior_predictive.R
# =============================================================================
# Reviewer-response plotting for the FINAL decoupled ABC-SMC fits.
#
# Drives off the saved .rds (each file IS the EasyABC `result`, with
# $param / $stats / $weights) plus its sidecar *_metadata.json (priors,
# summary_stats, observed_summaries, fit_params). Nothing is hardcoded per
# scenario beyond the list of fits to plot.
#
# We work from FOUR fits -- NS4 and NS5-with-takeoff (NOT the duration NS5),
# for each of the two archetypes (West Africa, DRC). A SEPARATE set of plots
# is produced for every fit so a final choice can be made downstream.
#
# For each fit it produces / prints:
#   (6) Posterior-predictive checks  -> gg_fit_ratio() | gg_pp_hist()
#         - sim/observed ratio pointranges AND histograms for EVERY summary
#           stat, now INCLUDING hcw_fraction (previously dropped).
#         - plus a quantitative goodness-of-fit table: posterior-predictive
#           median, 95% predictive interval, 95% coverage flag, and a
#           posterior-predictive p-value per summary stat.
#   (2) Prior vs posterior parameter distributions -> gg_prior_post()
#         - small number of bins; flat prior behind, posterior in front.
#         - plus a parameter-estimate table (weighted median + 95% CrI),
#           which answers "what were the fitted estimates" (incl. R0).
#   (4) Parameter correlation -> gg_pairs_lower()
#         - lower-triangle-only ggplot pairs plot + weighted correlation matrix.
#
# NOTE: written against the data structures in outputs/02_ABC_model_fits_Final/.
#       Not executed in-repo (no R here) -- run and report any errors to iterate.
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(ggplot2); library(jsonlite); library(cowplot)
})
`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- the four fits we work from --------------------------------------------
# id / archetype used for labelling + colour; rds resolved under outputs/.
OUT_SUB <- c("outputs", "02_ABC_model_fits_Final")
FITS <- list(
  WestAfrica_NS4 = list(
    archetype = "WestAfrica", label = "", #"NS4 (base 4)",
    rds = "fiber_ABC_SMC_Worst_WestAfrica_Decoupled_20260608_162044_check_NP5_NS4_NBREPS_30_NBSIMUL_472.rds"),
  # WestAfrica_NS5 = list(
  #   archetype = "WestAfrica", label = "NS5 (base 4 + takeoff)",
  #   rds = "fiber_ABC_SMC_Worst_WestAfrica_Decoupled_20260609_094145_check_NP5_NS5_NBREPS_30_NBSIMUL_472.rds"),
  DRC_NS4 = list(
    archetype = "DRC", label = "", #"NS4 (base 4)",
    rds = "fiber_ABC_SMC_Middle_DRC_ConflictSmoothed_PlusPlus_Decoupled_20260607_215621_check_NP5_NS4_NBREPS_30_NBSIMUL_590.rds")#,
  # DRC_NS5 = list(
  #   archetype = "DRC", label = "NS5 (base 4 + takeoff)",
  #   rds = "fiber_ABC_SMC_Middle_DRC_ConflictSmoothed_PlusPlus_Decoupled_20260609_201905_check_NP5_NS5_NBREPS_30_NBSIMUL_590.rds")
)

# ---- knobs ------------------------------------------------------------------
COUNTRY_COLS <- c(WestAfrica = "#d95f02", DRC = "#1b9e77")   # match repo archetype palette
NAME_FULL    <- c(WestAfrica = "West Africa", DRC = "DRC")

N_POST          <- 100000L   # posterior-predictive / posterior resample size
PP_SEED         <- 1L
PP_DROP_STATS   <- "takeoff" # stats to exclude from the posterior-predictive plots/table
                             #   (takeoff is a 0/1 indicator -> degenerate as a PP check;
                             #    still used in the fit, just not shown here)
PRIOR_POST_BINS <- 8L       # "small number of bins" for the prior-vs-posterior plots
PAIRS_N         <- 2000L    # points drawn in the pairs plot
SAVE_PLOTS      <- FALSE     # TRUE = write one file per fit to OUT_DIR
OUT_DIR         <- here(OUT_SUB[1], OUT_SUB[2], "reviewer_response")

PRETTY <- c(takeoff          = "Take-off (0/1)",
            log_n_deaths     = "log(Total Deaths)",
            log_n_hcw_deaths = "log(HCW deaths)",
            hcw_fraction     = "HCW fraction",
            d_p05_p95        = "death-date 5-95% span (days)",
            log_peak_height  = "log(Peak Height)")
pretty_of <- function(s) ifelse(s %in% names(PRETTY), PRETTY[s], s)

PARAM_PRETTY <- c(R0 = "R0", prop_funeral = "Unsafe funeral prop.",
                  etu_efficacy = "ETU efficacy", ppe_efficacy = "PPE efficacy",
                  hcw_risk_scalar = "HCW risk scalar")

# simulate n draws from one prior spec c(dist, a, b)
sim_prior_row <- function(spec, n) {
  d <- spec[1]; a <- as.numeric(spec[2]); b <- as.numeric(spec[3])
  switch(d, unif = runif(n, a, b), normal = rnorm(n, a, b),
         lognormal = rlnorm(n, a, b), stop("unhandled prior distribution: ", d))
}

# =============================================================================
# Load one fit: rds (param/stats/weights) + sidecar metadata json -------------
# =============================================================================
load_fit <- function(key) {
  sc      <- FITS[[key]]
  rds     <- here(OUT_SUB[1], OUT_SUB[2], sc$rds)
  metafile <- sub("\\.rds$", "_metadata.json", rds)
  stopifnot(file.exists(rds), file.exists(metafile))

  res <- readRDS(rds)
  md  <- jsonlite::fromJSON(metafile)

  stat_names  <- md$summary_stats
  fit_params  <- md$fit_params
  observed    <- setNames(as.numeric(md$observed_summaries), stat_names)
  priors      <- md$priors
  if (is.list(priors)) priors <- do.call(rbind, lapply(priors, unlist))  # -> char matrix n x 3

  stats <- as.data.frame(res$stats);  colnames(stats) <- stat_names
  param <- as.data.frame(res$param);  colnames(param) <- fit_params
  weights <- as.numeric(res$weights)          # res$weights can be an n x 1 matrix -> flatten

  # ONE weighted resample -> aligned posterior stats + params (same particles)
  set.seed(PP_SEED)
  idx <- sample(seq_len(nrow(stats)), size = N_POST, replace = TRUE, prob = weights)

  list(key = key, name = paste0(NAME_FULL[[sc$archetype]], " - ", sc$label),
       archetype = sc$archetype, col = COUNTRY_COLS[[sc$archetype]] %||% "#666666",
       stat_names = stat_names, fit_params = fit_params, observed = observed,
       priors = priors, weights = weights, param = param,
       post = stats[idx, , drop = FALSE], post_param = param[idx, , drop = FALSE])
}

# =============================================================================
# (6) posterior-predictive checks --------------------------------------------
# =============================================================================

# (A) simulated / observed fit-ratio -- one row per stat: a SMOOTHED ridgeline
# density of the posterior-predictive ratio, with the median + 95% pointrange
# drawn on top. `smooth` = density bandwidth multiplier (higher = smoother);
# `ridge_height` = how tall each density is in row units (<1 avoids overlap).
gg_fit_ratio <- function(fit, stats, col, smooth = 1.5, ridge_height = 0.85) {
  n     <- length(stats)
  ypos  <- setNames(rev(seq_len(n)), stats)     # first stat on top (highest y)
  labs  <- pretty_of(stats)

  # median + 95% CI of the ratio, per stat (the point + line)
  rdf <- do.call(rbind, lapply(stats, function(s) {
    q <- quantile(fit$post[[s]] / fit$observed[[s]], c(0.025, 0.5, 0.975), names = FALSE)
    data.frame(stat = s, y = ypos[[s]], lo = q[1], med = q[2], hi = q[3])
  }))

  # smoothed ridgeline density of the ratio, per stat (tails trimmed at 0.5/99.5%
  # so extreme draws don't squash the display; each row normalised to same height)
  dens <- do.call(rbind, lapply(stats, function(s) {
    r   <- fit$post[[s]] / fit$observed[[s]]
    rng <- quantile(r, c(0.005, 0.995), names = FALSE)
    d   <- density(r, adjust = smooth, from = rng[1], to = rng[2])
    data.frame(stat = s, y0 = ypos[[s]], x = d$x,
               h = d$y / max(d$y) * ridge_height)
  }))

  ggplot() +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
    geom_ribbon(data = dens, aes(x = x, ymin = y0, ymax = y0 + h, group = stat),
                fill = col, colour = col, alpha = 0.25, linewidth = 0.3) +
    geom_pointrange(data = rdf, aes(x = med, y = y, xmin = lo, xmax = hi),
                    colour = col, linewidth = 0.9, size = 0.5) +
    scale_y_continuous(breaks = unname(ypos[stats]), labels = labs,
                       expand = expansion(add = c(0.3, 0.15))) +   # ribbon already spans ridge_height; small top margin only
    labs(x = "Simulated / Observed", y = NULL) +
    theme_bw(base_size = 10)
}

# (B) posterior-predictive histograms, faceted; observed (black) + posterior
# median/95% (grey solid/dashed) overlaid.
gg_pp_hist <- function(fit, stats, fill_col, ncol = 2L) {
  post <- fit$post[stats]
  long <- utils::stack(post); names(long) <- c("value", "stat")
  long$stat <- factor(pretty_of(as.character(long$stat)), levels = pretty_of(stats))

  qdf <- do.call(rbind, lapply(stats, function(s) {
    q <- quantile(post[[s]], c(0.025, 0.5, 0.975), names = FALSE)
    data.frame(stat = pretty_of(s), value = q, kind = c("ci", "med", "ci"))
  }))
  qdf$stat <- factor(qdf$stat, levels = pretty_of(stats))
  odf <- data.frame(stat = factor(pretty_of(stats), levels = pretty_of(stats)),
                    value = as.numeric(fit$observed[stats]))

  ggplot(long, aes(value)) +
    geom_histogram(bins = 12, fill = fill_col, colour = "white", alpha = 0.85) +
    geom_vline(data = qdf, aes(xintercept = value, linetype = kind),
               colour = "grey25", linewidth = 0.5) +
    geom_vline(data = odf, aes(xintercept = value), colour = "black", linewidth = 1) +
    scale_linetype_manual(values = c(med = "solid", ci = "dashed"), guide = "none") +
    facet_wrap(~ stat, scales = "free", ncol = ncol) +
    labs(x = "Summary statistic value", y = "Frequency") +
    theme_bw(base_size = 10)
}

# NB: the combined posterior-predictive figure is assembled in the driver
# (pp_row), which excludes PP_DROP_STATS (takeoff) but keeps hcw_fraction.

# quantitative goodness-of-fit / coverage table (comment 6, "anything else")
gof_table <- function(fit) {
  do.call(rbind, lapply(setdiff(fit$stat_names, PP_DROP_STATS), function(s) {
    x   <- fit$post[[s]]; obs <- fit$observed[[s]]
    q   <- quantile(x, c(0.025, 0.5, 0.975), names = FALSE)
    ppp <- mean(x >= obs)                                   # posterior-predictive p-value
    data.frame(fit = fit$name, summary_stat = s,
               observed = round(obs, 4), pp_median = round(q[2], 4),
               pp_lo95 = round(q[1], 4), pp_hi95 = round(q[3], 4),
               obs_in_95 = obs >= q[1] && obs <= q[3],
               pp_pvalue = round(ppp, 3),
               row.names = NULL)
  }))
}

# =============================================================================
# (2) prior vs posterior parameter distributions ------------------------------
# =============================================================================
gg_prior_post <- function(fit, n_bins = PRIOR_POST_BINS) {
  recs <- list()
  for (pj in seq_along(fit$fit_params)) {
    pname  <- fit$fit_params[pj]
    a <- as.numeric(fit$priors[pj, 2]); b <- as.numeric(fit$priors[pj, 3])
    breaks <- seq(a, b, length.out = n_bins + 1)            # breaks aligned to prior bounds -> flat prior
    prior_draws <- sim_prior_row(fit$priors[pj, ], nrow(fit$post_param))
    post_draws  <- pmin(pmax(fit$post_param[[pname]], a), b)
    # divide bin counts by the number of draws -> probability mass per bin
    # (each distribution sums to 1 across its bins; prior & posterior use the
    # same n draws so the two are directly comparable).
    prc <- hist(prior_draws, breaks = breaks, plot = FALSE, include.lowest = TRUE)$counts / length(prior_draws)
    poc <- hist(post_draws,  breaks = breaks, plot = FALSE, include.lowest = TRUE)$counts / length(post_draws)
    base <- data.frame(parameter = unname(PARAM_PRETTY[pname] %||% pname),
                       xmin = head(breaks, -1), xmax = tail(breaks, -1))
    recs[[length(recs) + 1]] <- cbind(base, source = "Prior",     count = prc)
    recs[[length(recs) + 1]] <- cbind(base, source = "Posterior", count = poc)
  }
  d <- do.call(rbind, recs)
  d$parameter <- factor(d$parameter, levels = unname(PARAM_PRETTY[fit$fit_params]))
  d$source    <- factor(d$source, levels = c("Prior", "Posterior"))   # Prior drawn first = behind
  ggplot(d, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = count, fill = source)) +
    geom_rect(alpha = 0.6, colour = NA) +
    facet_wrap(~ parameter, scales = "free", nrow = 1) +
    scale_fill_manual(values = c(Prior = "grey70", Posterior = fit$col)) +
    labs(x = "Parameter value", y = "Probability mass",
         fill = NULL) +                                     # no title; legend handled by driver
    theme_bw(base_size = 10)
}

# weighted median + 95% CrI per parameter (answers "fitted estimates", incl R0)
param_estimate_table <- function(fit) {
  do.call(rbind, lapply(fit$fit_params, function(p) {
    q <- quantile(fit$post_param[[p]], c(0.025, 0.5, 0.975), names = FALSE)
    data.frame(fit = fit$name, parameter = p,
               median = round(q[2], 3), lo95 = round(q[1], 3), hi95 = round(q[3], 3),
               row.names = NULL)
  }))
}

# =============================================================================
# (4) parameter correlation ---------------------------------------------------
# =============================================================================
# LOWER-TRIANGLE-ONLY pairs plot, built from ggplot panels arranged with cowplot
# (no GGally, no duplicated upper triangle). Diagonal = variable name; lower
# cells = posterior scatter with the weighted correlation printed; upper = blank.
# Returns list(plot, cor) where cor is the exact weighted correlation matrix.
gg_pairs_lower <- function(fit) {
  ps <- fit$fit_params; np <- length(ps)
  dp <- fit$post_param[sample(nrow(fit$post_param), min(PAIRS_N, nrow(fit$post_param))),
                       ps, drop = FALSE]
  cormat <- cov.wt(fit$param[ps], wt = fit$weights, cor = TRUE)$cor   # exact weighted cor
  labs <- unname(PARAM_PRETTY[ps]); labs[is.na(labs)] <- ps[is.na(labs)]

  # reduced (np-1) x (np-1) grid: rows = params 2..np (y), cols = params 1..np-1 (x).
  # dropping the fully-empty top row and right column makes the populated triangle
  # fill the whole plotting space (bottom row + left column reach every edge).
  nr <- np - 1L
  cells <- vector("list", nr * nr)
  for (r in seq_len(nr)) for (cc in seq_len(nr)) {
    i <- r + 1L; j <- cc                                    # y-param i, x-param j
    idx <- (r - 1L) * nr + cc
    if (cc <= r) {                                          # lower triangle (j < i): scatter
      d <- data.frame(x = dp[[ps[j]]], y = dp[[ps[i]]])
      cells[[idx]] <- ggplot(d, aes(x, y)) +
        geom_point(alpha = 0.15, size = 0.5, colour = fit$col) +
        annotate("text", x = Inf, y = Inf, label = sprintf("r = %.2f", cormat[i, j]),
                 hjust = 1.1, vjust = 1.4, size = 3.4, colour = "grey20") +
        labs(x = if (i == np) labs[j] else NULL,           # x labels only on bottom row
             y = if (j == 1)  labs[i] else NULL) +          # y labels only on first column
        theme_bw(base_size = 10) +
        theme(axis.title = element_text(size = 9),
              axis.text  = element_text(size = 7))
    } else {
      cells[[idx]] <- NULL                                  # upper triangle: blank
    }
  }
  list(plot = cowplot::plot_grid(plotlist = cells, ncol = nr, align = "hv"),
       cor  = cormat)                                       # no title (per-fit title added by driver if needed)
}

# =============================================================================
# DRIVER -- combined West Africa (A) / DRC (B) figures + tables ---------------
#   - posterior-predictive checks : A,B (West Africa) over C,D (DRC)
#   - prior vs posterior          : A (West Africa) over B (DRC), shared legend
#   - parameter pairs             : A (West Africa) over B (DRC)
# Picks one fit per archetype from the (uncommented) FITS list.
# =============================================================================
if (SAVE_PLOTS) dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# load all (uncommented) fits; pick one per archetype -> A = West Africa, B = DRC
fits <- lapply(names(FITS), load_fit); names(fits) <- names(FITS)
by_arch <- function(a) { for (f in fits) if (f$archetype == a) return(f); NULL }
wa  <- by_arch("WestAfrica")
drc <- by_arch("DRC")
stopifnot(!is.null(wa), !is.null(drc))

# ---- (6) combined posterior-predictive checks: A,B (West Africa) / C,D (DRC) --
pp_row <- function(fit, labels) {
  stats <- setdiff(fit$stat_names, PP_DROP_STATS)
  cowplot::plot_grid(gg_fit_ratio(fit, stats, fit$col),
                     gg_pp_hist(fit, stats, fit$col, ncol = 2L),
                     labels = labels, rel_widths = c(0.8, 1.2), nrow = 1)
}
combined_pp <- cowplot::plot_grid(pp_row(wa,  c("A", "B")),
                                  pp_row(drc, c("C", "D")), ncol = 1)
print(combined_pp)

# ---- (2) combined prior vs posterior: A (West Africa) / B (DRC), shared legend right
pvp_wa  <- gg_prior_post(wa)
pvp_drc <- gg_prior_post(drc)
pvp_leg  <- cowplot::get_legend(pvp_wa + theme(legend.position = "right"))
pvp_body <- cowplot::plot_grid(pvp_wa  + theme(legend.position = "none"),
                               pvp_drc + theme(legend.position = "none"),
                               ncol = 1, labels = c("A", "B"))
combined_pvp <- cowplot::plot_grid(pvp_body, pvp_leg, ncol = 2, rel_widths = c(1, 0.13))
print(combined_pvp)

# ---- (4) combined pairs: A (West Africa) / B (DRC) stacked --------------------
prs_wa  <- gg_pairs_lower(wa)
prs_drc <- gg_pairs_lower(drc)
# add a thin blank strip above each grid so the A/B labels sit ABOVE the top-row
# y-axis titles instead of overlapping them.
pad_top <- function(g) cowplot::plot_grid(NULL, g, ncol = 1, rel_heights = c(0.04, 1))
combined_pairs <- cowplot::plot_grid(
  pad_top(prs_wa$plot), pad_top(prs_drc$plot), ncol = 1,
  labels = c("A", "B"), label_size = 16, label_x = 0.005, label_y = 1,
  hjust = 0, vjust = 1)
print(combined_pairs)

# supplementary pairs figure: 7 wide x 11 high
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(OUT_DIR, "Figure_S4_PairsPlot.pdf"), combined_pairs, width = 7, height = 11)

# ---- tables over all active fits --------------------------------------------
gof_all <- lapply(fits, gof_table)
est_all <- lapply(fits, param_estimate_table)
cor_all <- list(WestAfrica = prs_wa$cor, DRC = prs_drc$cor)

if (SAVE_PLOTS) {
  ggsave(file.path(OUT_DIR, "pp_checks_WA_DRC.pdf"),          combined_pp,    width = 10, height = 8)
  ggsave(file.path(OUT_DIR, "prior_vs_posterior_WA_DRC.pdf"), combined_pvp,   width = 12, height = 6)
}

# ---- combined tables (printed; optionally written) --------------------------
gof_table_all <- do.call(rbind, gof_all);  rownames(gof_table_all) <- NULL
est_table_all <- do.call(rbind, est_all);  rownames(est_table_all) <- NULL

cat("\n===== Goodness-of-fit / coverage (posterior-predictive) =====\n")
print(gof_table_all)
cat("\n===== Fitted parameter estimates (weighted median [95% CrI]) =====\n")
print(est_table_all)
cat("\n===== Posterior parameter correlation matrices =====\n")
for (k in names(cor_all)) { cat("\n--", k, "--\n"); print(round(cor_all[[k]], 2)) }

if (SAVE_PLOTS) {
  write.csv(gof_table_all, file.path(OUT_DIR, "gof_coverage_table.csv"), row.names = FALSE)
  write.csv(est_table_all, file.path(OUT_DIR, "parameter_estimates_table.csv"), row.names = FALSE)
}
