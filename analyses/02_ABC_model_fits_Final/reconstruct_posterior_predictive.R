# reconstruct_posterior_predictive.R
# =============================================================================
# Read the two FINAL decoupled ABC-SMC fits and draw, FOR EACH COUNTRY, ONE figure
# combining (A) the simulated/observed fit-ratio plot on the left and (B) the 2x2
# of posterior-predictive summary-stat histograms on the right. DRC and West
# Africa are drawn as SEPARATE figures, each in its country colour.
#
# Driven straight off the saved .RDS (each file IS the `result`, with $param/
# $stats/$weights and a self-describing $run_metadata, from which the observed
# targets + stat names are read -- nothing hardcoded per scenario).
#
# The plotting functions RETURN ggplot objects; the driver arranges them with
# cowplot::plot_grid (labels A / B).
# =============================================================================

suppressPackageStartupMessages({ library(here); library(ggplot2) })
`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- the two fits (as supplied) --------------------------------------------
SCENARIOS <- list(
  WestAfrica = list(
    id           = "Worst_WestAfrica",
    rds          = here("outputs", "02_ABC_model_fits_Final",
                        "fiber_ABC_SMC_Worst_WestAfrica_Decoupled_20260608_162044_check_NP5_NS4_NBREPS_30_NBSIMUL_472.RDS"),
    scenario_csv = here("data-processed", "final_six_scenario_values_original_approach.csv"),
    param_names  = c("R0", "prop_funeral", "etu_efficacy", "ppe_efficacy", "hcw_risk_scalar"),
    check_final_size = 40000
  ),
  DRC = list(
    id           = "Middle_DRC_ConflictSmoothed_PlusPlus",
    rds          = here("outputs", "02_ABC_model_fits_Final",
                        "fiber_ABC_SMC_Middle_DRC_ConflictSmoothed_PlusPlus_Decoupled_20260607_215621_check_NP5_NS4_NBREPS_30_NBSIMUL_590.RDS"),
    scenario_csv = here("data-processed", "final_six_scenario_values_original_approach.csv"),
    param_names  = c("R0", "prop_funeral", "etu_efficacy", "ppe_efficacy", "hcw_risk_scalar"),
    check_final_size = 10000
  )
)

# ---- knobs ------------------------------------------------------------------
# >>> COUNTRY COLOURS <<< no per-country palette exists in the repo (it colours by
# arm: Okabe-Ito #D55E00 / #0072B2). Defaulting to that pair, one per country --
# SET THESE to whatever you've used for DRC / West Africa elsewhere.
COUNTRY_COLS <- c(WestAfrica = "#D55E00", DRC = "#0072B2")
NAME_FULL    <- c(WestAfrica = "West Africa", DRC = "DRC")

# Summaries to drop so the histograms form a clean 2x2 (hcw_fraction is the
# redundant one -- "deliberately redundant with the two count logs"). Set to
# character(0) to show all 5 (the grid then becomes 3x2).
DROP_STATS <- "hcw_fraction"

PRETTY <- c(log_n_deaths     = "log(deaths)",
            log_n_hcw_deaths = "log(HCW deaths)",
            hcw_fraction     = "HCW fraction",
            d_p05_p95        = "death-date 5-95% span (days)",
            log_peak_height  = "log(peak height)")
pretty_of <- function(s) ifelse(s %in% names(PRETTY), PRETTY[s], s)

N_POST     <- 10000L     # posterior-predictive resample size
PP_SEED    <- 1L
SAVE_PLOTS <- FALSE      # TRUE = also write one PNG per country to OUT_DIR
OUT_DIR    <- here("outputs", "02_ABC_model_fits_Final")

# =============================================================================
# Load + posterior-predictive resample (one fit) -----------------------------
load_fit <- function(sc) {
  stopifnot(file.exists(sc$rds))
  res <- readRDS(sc$rds)
  md  <- res$run_metadata %||% list()
  stat_names <- md$summary_stats %||% colnames(res$stats)
  observed   <- md$observed_summaries
  if (is.null(stat_names) || is.null(observed))
    stop("No summary_stats/observed_summaries in $run_metadata for ", sc$id, ".")
  stats <- as.data.frame(res$stats); colnames(stats) <- stat_names
  set.seed(PP_SEED)
  idx <- sample(seq_len(nrow(stats)), size = N_POST, replace = TRUE, prob = res$weights)
  list(res = res, md = md, stat_names = stat_names, observed = observed,
       post = stats[idx, , drop = FALSE])
}

# (A) simulated / observed fit-ratio -- horizontal pointrange, one row per stat.
gg_fit_ratio <- function(fit, stats, col) {
  rdf <- do.call(rbind, lapply(stats, function(s) {
    q <- quantile(fit$post[[s]] / fit$observed[[s]], c(0.025, 0.5, 0.975), names = FALSE)
    data.frame(stat = pretty_of(s), lo = q[1], med = q[2], hi = q[3])
  }))
  rdf$stat <- factor(rdf$stat, levels = rev(pretty_of(stats)))   # first stat on top
  ggplot(rdf, aes(med, stat)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
    geom_pointrange(aes(xmin = lo, xmax = hi), colour = col, linewidth = 0.9, size = 0.5) +
    labs(x = "Simulated / Observed", y = NULL,
         caption = "ratio on the fitted scale (log for counts)") +
    theme_bw(base_size = 10) +
    theme(plot.caption = element_text(size = 7, colour = "grey40"))
}

# (B) posterior-predictive histograms, faceted 2x2; bars in the country colour,
# observed target (black) and posterior median/95% (grey solid/dashed) overlaid.
gg_pp_hist <- function(fit, stats, fill_col, ncol = 2L) {
  post <- fit$post[stats]
  long <- utils::stack(post)                       # -> values, ind
  names(long) <- c("value", "stat")
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
    labs(x = NULL, y = "Draws",
         subtitle = "black = observed; grey = posterior median (solid) / 95% (dashed)") +
    theme_bw(base_size = 10) +
    theme(strip.text = element_text(face = "bold"),
          plot.subtitle = element_text(size = 7, colour = "grey30"))
}

# build the combined A | B figure for one country
country_figure <- function(name, fit) {
  stats <- setdiff(names(fit$observed), DROP_STATS)   # 4 -> clean 2x2
  col   <- COUNTRY_COLS[[name]] %||% "#0072B2"
  body  <- cowplot::plot_grid(
    gg_fit_ratio(fit, stats, col),
    gg_pp_hist(fit, stats, col, ncol = 2L),
    labels = c("A", "B"), rel_widths = c(0.8, 1.2), nrow = 1)
  title <- cowplot::ggdraw() +
    cowplot::draw_label(sprintf("%s -- posterior-predictive checks", NAME_FULL[[name]] %||% name),
                        fontface = "bold", x = 0.01, hjust = 0, size = 13)
  cowplot::plot_grid(title, body, ncol = 1, rel_heights = c(0.08, 1))
}

# =============================================================================
# DRIVER -- one SEPARATE figure per country (returned in `figs`) --------------
figs <- list()
for (name in names(SCENARIOS)) {
  message("==== ", name, " ====")
  fit  <- load_fit(SCENARIOS[[name]])
  figs[[name]] <- country_figure(name, fit)
  print(figs[[name]])
  if (SAVE_PLOTS) {
    f <- file.path(OUT_DIR, sprintf("posterior_predictive_%s.png", name))
    ggsave(f, figs[[name]], width = 11, height = 5, dpi = 200)
    message("  wrote ", f)
  }
}
