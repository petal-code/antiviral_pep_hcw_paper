# ============================================================================
# DRC_conflict_checking_RW_script.R   (DIAGNOSTIC / EXPLORATION -- not pipeline)
# ----------------------------------------------------------------------------
# PURPOSE
#   A second alternative to Model B for the DRC conflict Q (companion to
#   DRC_conflict_checking_alts_script.R). Here the response-quality curve is a
#   FLEXIBLE latent function fit DIRECTLY to the safe-and-dignified-burial (SDB)
#   success COUNTS with a binomial likelihood -- so the uncertainty is real
#   (wide where few burials were observed, tight where many) rather than the
#   zero uncertainty of the deterministic rolling-mean curve that Model B uses.
#
# THE MODEL
#   For each weekly bin t, let n_t = eligible burials and y_t = successful ones.
#       y_t ~ Binomial(n_t, p_t),   logit(p_t) = f_t
#   and f_t follows a random walk over the (ordered) weekly bins:
#       RW1 (order 1):  f_t ~ Normal(f_{t-1}, sigma)            -- flexible level
#       RW2 (order 2):  f_t ~ Normal(2 f_{t-1} - f_{t-2}, sigma) -- smoother (slope)
#   sigma controls smoothness (penalises week-to-week jumps). This is the cheap,
#   stable, discrete cousin of a Gaussian process; it makes NO shape assumption
#   (so it handles the non-monotonic rise-dip-recover) but, unlike a GP, lives on
#   the data grid, so it does not mean-revert / fan out where there are no data.
#
#   The relative response-quality index used elsewhere is then
#       Q_t = p_t / max(p)         (max-scaled per posterior draw),
#   so we can put it on the same footing as the Model B conflict Q.
#
# CONFLICT vs CONFLICT++
#   conflict   : the real SDB counts.
#   conflict++ : the same bins, but the successful counts are forced to 0 over
#                days 200-300 (the count-space analogue of the pipeline's forced
#                response collapse). n is unchanged: "those weeks happened, but
#                every burial failed".
#
# WHAT TO LOOK FOR
#   * The WIDTH of the Q / p ribbons -- this is the whole point. Compare against
#     the Model B fixed curve (dashed): how much genuine uncertainty was being
#     thrown away by treating the rolling mean as known?
#   * Whether the flexible fit just chases SDB noise (RW1) or stays sensibly
#     smooth (RW2) -- flip RW_ORDER to see.
#   * It is still a SINGLE data source (burials). This does not address the
#     SDB-domination question; it only quantifies the SDB curve honestly.
#
# CAVEATS
#   * Counts are POOLED across provinces (volume-weighted), which differs from
#     the pipeline's unweighted province-mean-of-proportions -- a deliberate,
#     more standard choice for a count likelihood.
#   * Models only the OBSERVED window; extending Q to day 730 would still need
#     the pipeline's hold-flat convention past the last SDB bin.
#   * The binomial naturally down-weights small-n bins, so the ad-hoc early-spike
#     suppression used in 00 is not applied here.
#
# Inputs : data-processed/DRC_QCurve/DRC_QCurve_PreppedData.rds
#          (uses $province_weekly_qc for the counts, and $conflict_qseries /
#           $conflict_plusplus_qseries for the Model B reference curves)
# Output : none -- two comparison plots are printed to the graphics device only.
# ============================================================================

suppressPackageStartupMessages({
  library(cmdstanr)
  library(dplyr); library(tidyr); library(stringr)
  library(readr); library(tibble); library(ggplot2)
})

source(here::here("analyses", "01_latent_response_parameter_estimation", "helpers.R"))
set.seed(123)

# Random-walk order: 1 = flexible level (wigglier), 2 = penalise curvature (smoother).
RW_ORDER   <- 2L
# Exponential prior rate on the RW innovation sd (smoothness knob). RW2 wants a
# tighter innovation scale than RW1, hence a larger rate (smaller mean).
SIGMA_RATE <- if (RW_ORDER == 1L) 2 else 5

PLUSPLUS_WINDOW <- c(200, 300)   # day window over which conflict++ forces success -> 0

# ----------------------------------------------------------------------------
# 1. The random-walk binomial model (inline Stan)
# ----------------------------------------------------------------------------
stan_code <- '
data {
  int<lower=1> T;                 // number of weekly bins
  int<lower=1, upper=2> rw_order; // 1 = RW on level, 2 = RW on slope
  array[T] int<lower=0> n;        // eligible burials per bin
  array[T] int<lower=0> y;        // successful burials per bin (y <= n)
  real<lower=0> sigma_rate;       // exponential prior rate on the RW innovation sd
}
parameters {
  real f_init1;                   // logit success at bin 1
  real f_init2;                   // logit success at bin 2 (used only when rw_order == 2)
  vector[T] z;                    // standardised RW innovations (leading entries unused)
  real<lower=0> sigma;            // RW innovation scale (smoothness)
}
transformed parameters {
  vector[T] f;                    // latent logit success at each bin (non-centred RW)
  if (rw_order == 1) {
    f[1] = f_init1;
    for (t in 2:T) f[t] = f[t-1] + sigma * z[t];
  } else {
    f[1] = f_init1;
    f[2] = f_init2;
    for (t in 3:T) f[t] = 2 * f[t-1] - f[t-2] + sigma * z[t];
  }
}
model {
  f_init1 ~ normal(0, 3);         // broad prior on the initial level (logit scale)
  f_init2 ~ normal(0, 3);
  z       ~ normal(0, 1);         // standardised innovations (z[1], and z[2] for RW2, just sample the prior)
  sigma   ~ exponential(sigma_rate);
  // Honest count likelihood: a bin with few eligible burials informs f only weakly.
  y ~ binomial_logit(n, f);
}
generated quantities {
  vector[T] p = inv_logit(f);     // success probability per bin (absolute scale)
  vector[T] Q = p / max(p);       // relative response-quality index (max-scaled), per draw
}
'
mod <- cmdstan_model(write_stan_file(stan_code))

# ----------------------------------------------------------------------------
# 2. Build the weekly success counts (pooled across provinces)
# ----------------------------------------------------------------------------
drc_prep <- readRDS(file.path(DIR_PROCESSED, "DRC_QCurve/DRC_QCurve_PreppedData.rds"))

# counts: one row per weekly bin, with y = successful and n = eligible burials,
# pooled across the two provinces. (This is the data the binomial sees.)
counts <- drc_prep$province_weekly_qc %>%
  group_by(relative_day) %>%
  summarise(y = sum(n_successful_response), n = sum(n_eligible), .groups = "drop") %>%
  filter(n > 0) %>%
  arrange(relative_day)

# conflict   : the real counts.
# conflict++ : successful counts forced to 0 inside the collapse window (n kept).
counts_conflict <- counts
counts_pluspl   <- counts %>%
  mutate(y = if_else(relative_day >= PLUSPLUS_WINDOW[1] & relative_day <= PLUSPLUS_WINDOW[2],
                     0L, y))

# Model B reference curves (the deterministic fixed-Q that this model is an
# alternative to), for overlay.
ref <- bind_rows(
  transmute(drc_prep$conflict_qseries,          relative_day, q_value, success_smoothed, scenario = "conflict"),
  transmute(drc_prep$conflict_plusplus_qseries, relative_day, q_value, success_smoothed, scenario = "conflict++")
)

# ----------------------------------------------------------------------------
# 3. Fit one scenario
# ----------------------------------------------------------------------------
# Returns tidy posterior curves for p (success probability) and Q (relative
# index), one row per bin, tagged with the scenario label.
fit_rw <- function(counts_df, label) {
  message("\n==== Fitting RW", RW_ORDER, " binomial model: ", label, " ====")

  stan_data <- list(
    T = nrow(counts_df), rw_order = RW_ORDER,
    n = counts_df$n, y = counts_df$y, sigma_rate = SIGMA_RATE
  )

  fit <- mod$sample(
    data = stan_data, seed = 123,
    chains = 4, parallel_chains = 4,
    iter_warmup = 1000, iter_sampling = 1000,
    adapt_delta = 0.95, max_treedepth = 12, refresh = 250
  )
  cat("\nDiagnostics (", label, "):\n", sep = ""); print(fit$diagnostic_summary())
  cat("\nRW innovation sd sigma (", label, "):\n", sep = "")
  print(fit$summary(variables = "sigma")[, c("variable", "mean", "q5", "q95")])

  # Pull a per-bin vector variable into a tidy curve keyed by relative_day.
  pull_curve <- function(varname) {
    fit$summary(variables = varname) %>%
      mutate(t = as.integer(str_match(variable, "\\[([0-9]+)\\]")[, 2]),
             relative_day = counts_df$relative_day[t],
             scenario = label) %>%
      select(scenario, relative_day, mean, q5, q95)
  }
  list(p = pull_curve("p"), Q = pull_curve("Q"))
}

# ----------------------------------------------------------------------------
# 4. Fit conflict and conflict++ and overlay
# ----------------------------------------------------------------------------
fit_conflict <- fit_rw(counts_conflict, "conflict")
fit_pluspl   <- fit_rw(counts_pluspl,   "conflict++")

p_all <- bind_rows(fit_conflict$p, fit_pluspl$p)
Q_all <- bind_rows(fit_conflict$Q, fit_pluspl$Q)

# Raw pooled proportions (the data), sized by sample size, for the p plot.
data_pts <- bind_rows(
  transmute(counts_conflict, relative_day, prop = y / n, n, scenario = "conflict"),
  transmute(counts_pluspl,   relative_day, prop = y / n, n, scenario = "conflict++")
)

scen_cols <- c("conflict" = "#1f77b4", "conflict++" = "#d62728")

# --- Plot 1: success probability p(t), with the count data and Model B reference ---
p_success <- ggplot(p_all, aes(relative_day, mean, colour = scenario, fill = scenario)) +
  geom_ribbon(aes(ymin = q5, ymax = q95), alpha = 0.15, colour = NA) +   # 90% credible band
  geom_line(linewidth = 0.9) +
  geom_point(data = data_pts, aes(relative_day, prop, size = n),
             inherit.aes = FALSE, colour = "grey40", alpha = 0.45) +
  geom_line(data = ref, aes(relative_day, success_smoothed, colour = scenario),
            inherit.aes = FALSE, linetype = "dashed", linewidth = 0.7) +
  scale_size_area(max_size = 4) +
  scale_colour_manual(values = scen_cols) + scale_fill_manual(values = scen_cols) +
  labs(title = paste0("DRC SDB success probability: RW", RW_ORDER, " binomial fit"),
       subtitle = "Solid + band = posterior mean & 90% interval; points = weekly data (size = eligible burials); dashed = Model B fixed curve",
       x = "Relative outbreak day", y = "P(successful burial)", colour = NULL, fill = NULL, size = "eligible") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")

# --- Plot 2: the relative response-quality index Q(t), vs the Model B fixed Q ---
p_Q <- ggplot(Q_all, aes(relative_day, mean, colour = scenario, fill = scenario)) +
  geom_ribbon(aes(ymin = q5, ymax = q95), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) +
  geom_line(data = ref, aes(relative_day, q_value, colour = scenario),
            inherit.aes = FALSE, linetype = "dashed", linewidth = 0.7) +
  scale_colour_manual(values = scen_cols) + scale_fill_manual(values = scen_cols) +
  labs(title = paste0("DRC response-quality Q(t): RW", RW_ORDER, " binomial fit vs Model B fixed Q"),
       subtitle = "Solid + band = flexible fit with propagated uncertainty; dashed = Model B fixed (zero-uncertainty) Q",
       x = "Relative outbreak day", y = "Q (max-scaled, 0-1)", colour = NULL, fill = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")

print(p_success)   # display only; not saved
print(p_Q)         # display only; not saved

message("\nDRC_conflict_checking_RW_script.R complete (RW", RW_ORDER,
        " fits done; plots printed, nothing saved).")
