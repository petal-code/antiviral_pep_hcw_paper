# ============================================================================
# 02_DRC_QCurve_Fitting_Original.R
# ----------------------------------------------------------------------------
# WHAT THIS SCRIPT DOES
#   Fits the DRC "conflict" and "conflict++" response-parameter curves using
#   MODEL B, and saves them for the combine step (03).
#
# THE MODEL IN ONE PARAGRAPH (and how it differs from Model A)
#   Unlike Model A (West Africa), Model B does NOT estimate the shape of the
#   response-quality curve Q. Instead a single shared Q(t) curve is supplied as
#   DATA -- it is the empirical fraction of safe-and-dignified burials that were
#   successful over time (reconstructed in 00 from the Warsame line-list). Every
#   parameter rides that same fixed curve between its own two magnitude
#   endpoints, and ONLY those endpoints are estimated. Two consequences:
#     * the jagged, conflict-interrupted shape of the empirical curve is
#       preserved exactly (a smooth fit would erase the very signal we want); and
#     * there is no curve shape to pool, so (correctly) there is no partial
#       pooling here -- each parameter's endpoints are estimated independently.
#
#   ONE SPECIAL CASE: community unsafe funerals are reported on the ABSOLUTE
#   Warsame scale (1 - success), not the relative Q scale, so that parameter's
#   fitted curve is overridden deterministically at the end (see step 4).
#
#   "conflict++" is identical to "conflict" except the input Q series already has
#   a forced collapse (success -> 0) baked in over days 200-300 (done in 00), so
#   every parameter collapses with it -- no extra logic is needed here.
#
# Inputs : data-processed/DRC_QCurve_PreppedData.rds   (from 00; uses $anchors,
#          $conflict_qseries, $conflict_plusplus_qseries)
# Stan   : stan-models/modelB_fixedQ_boundsOnly.stan
# Output : data-processed/drc_conflict_fit.rds
#          data-processed/drc_conflict_plusplus_fit.rds
#
# NOTE on the no-conflict scenario: by design it is NOT fit here. Its data are
# still prepared in 00, and it is explored separately in DRC_no_conflict_checking.R
# (Model B vs Model A). Its horizon is used by 03 to time-stretch the West Africa
# conflict scenarios.
# ============================================================================

suppressPackageStartupMessages({
  library(cmdstanr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(tibble)
  library(ggplot2)   # display-only plots of the fitted curves at the end
})

source(here::here("analyses", "01_latent_response_parameter_estimation", "helpers.R"))

set.seed(123)

# ----------------------------------------------------------------------------
# 1. Read the bundled DRC inputs
# ----------------------------------------------------------------------------
# DRC_QCurve_PreppedData.rds (from 00) bundles the cleaned anchors, the three Q series, the
# durations and a QC table. Here we need the anchors now and the two conflict Q
# series further down.
drc_prep    <- readRDS(file.path(DIR_PROCESSED, "DRC_QCurve/DRC_QCurve_PreppedData.rds"))
drc_anchors <- drc_prep$anchors

# ----------------------------------------------------------------------------
# 2. Per-parameter metadata (DRC admissible domains)
# ----------------------------------------------------------------------------
# Hard admissible domain [abs_min, abs_max] per parameter. These differ slightly
# from West Africa: the delay-to-hospitalisation ceiling is 12 days (vs 10) and
# hospital unsafe funerals are capped much tighter at 0.010 (vs 0.10).
domain_meta <- tribble(
  ~parameter,                ~abs_min, ~abs_max,
  "delay_hosp",              0.5,      12.0,
  "p_hosp",                  0.0,      1.0,
  "p_ETU",                   0.0,      1.0,
  "latent_IPC",              0.0,      1.0,
  "p_unsafe_funeral_comm",   0.0,      1.0,
  "p_unsafe_funeral_hosp",   0.0,      0.010
)

# Collapse anchors to one row per parameter and add the endpoint support/priors.
# (param_meta has the same columns as in 01 -- see that script's section 2 for
# the full column-by-column description.)
param_meta <- drc_anchors %>%
  group_by(parameter) %>%
  summarise(
    lb_prior_mean = first(lower_bound),   # prior centre for the LOWER endpoint
    ub_prior_mean = first(upper_bound),   # prior centre for the UPPER endpoint
    direction     = first(direction),
    .groups = "drop"
  ) %>%
  left_join(domain_meta, by = "parameter") %>%
  mutate(
    param_id  = match(parameter, PARAM_LEVELS),
    increases = if_else(direction == "up", 1L, 0L)
  ) %>%
  arrange(param_id) %>%
  build_param_support()

J <- nrow(param_meta)

# Lognormal prior on each parameter's observation-noise fraction (sigma_y_j is a
# fraction of that parameter's own span). Shared across parameters and scenarios.
SIGMA_FRAC_PRIOR_MEANLOG <- log(0.12)
SIGMA_FRAC_PRIOR_SDLOG   <- 0.60

# Compile the Stan model once and reuse it for both scenario fits.
mod <- cmdstan_model(file.path(DIR_STAN, "modelB_fixedQ_boundsOnly.stan"))

# ----------------------------------------------------------------------------
# 3. Fit one DRC scenario (called twice: conflict, then conflict++)
# ----------------------------------------------------------------------------
# Given a Q series (the shared empirical curve) and the DRC anchors, this:
#   1. maps each anchor onto the shared Q (q_obs = Q at the anchor's day),
#   2. builds a dense prediction grid of the shared Q,
#   3. fits Model B (estimating only the endpoints), and
#   4. returns the per-parameter native-unit curves, with the community
#      unsafe-funeral curve replaced by its ABSOLUTE Warsame proxy (1 - success).
fit_drc_scenario <- function(qseries, label) {

  message("\n==== Fitting DRC scenario: ", label, " ====")

  # The shared horizon spans both the literature anchors and the SDB Q series.
  max_day <- max(c(drc_anchors$relative_day, qseries$relative_day), na.rm = TRUE)

  # The largest empirical success value sets the ABSOLUTE floor of community
  # unsafe funerals (1 - max success); used in the endpoint tweak prior below.
  success_max <- max(qseries$success_smoothed, na.rm = TRUE)

  # ---- 1. Map each anchor onto the shared Q --------------------------------
  # Each anchor sits at some outbreak day; q_obs is the shared response-quality
  # index Q read off the SDB curve at that day. So an anchor early in the
  # response (Q near 0) informs that parameter's worst-state endpoint, and a late
  # anchor (Q near 1) informs its best-state endpoint.
  q_at_day <- make_interp(qseries$relative_day, qseries$q_value)
  # fit_df columns: param_id (which parameter), q_obs (shared Q at the anchor's
  # day, in [0,1]), y_obs (observed value), obs_sd_mult (= 1/weight; looser when
  # larger).
  fit_df <- drc_anchors %>%
    left_join(select(param_meta, parameter, param_id, increases), by = "parameter") %>%
    mutate(
      q_obs       = clip01(q_at_day(relative_day)),
      y_obs       = value_used,
      weight      = if_else(is.na(weight), 1, weight),
      obs_sd_mult = 1 / pmax(weight, 0.25)
    ) %>%
    arrange(param_id)

  # ---- 2. Dense prediction grid of the shared Q ----------------------------
  # A 250-point regular grid in normalised time, UNIONED with the exact empirical
  # support points (tau_q) so the predicted curve passes through the real data
  # points (including the empirical maximum where Q = 1). q_pred is the shared Q
  # interpolated onto that grid.
  tau_grid <- sort(unique(c(seq(0, 1, length.out = 250), qseries$tau_q)))
  q_on_tau <- make_interp(qseries$tau_q, qseries$q_value)
  q_pred   <- clip01(q_on_tau(tau_grid))

  # ---- 3. Endpoint tweak priors (Model B exposes per-parameter upper/lower) -
  # Start with all tweaks OFF, then anchor the two unsafe-funeral parameters:
  #   community: high end (at Q=0) pinned to 1; low end (at Q=1) to the absolute
  #              Warsame floor 1 - max(success).
  #   hospital : strongly regularised to a near-zero band.
  zero_flags <- rep(0L, J)
  tw <- list(
    use_upper = zero_flags, upper_mean = param_meta$ub_prior_mean, upper_sd = rep(1.0, J),
    use_lower = zero_flags, lower_mean = param_meta$lb_prior_mean, lower_sd = rep(1.0, J)
  )
  # Small helper: switch on the upper and/or lower tweak for one named parameter.
  set_tweak <- function(tw, p, upper = NULL, lower = NULL) {
    j <- param_meta$param_id[param_meta$parameter == p]
    if (!is.null(upper)) { tw$use_upper[j] <- 1L; tw$upper_mean[j] <- upper$mean; tw$upper_sd[j] <- upper$sd }
    if (!is.null(lower)) { tw$use_lower[j] <- 1L; tw$lower_mean[j] <- lower$mean; tw$lower_sd[j] <- lower$sd }
    tw
  }
  tw <- set_tweak(tw, "p_unsafe_funeral_comm",
                  upper = list(mean = 1.00,                    sd = 0.02),
                  lower = list(mean = clip01(1 - success_max), sd = 0.02))
  tw <- set_tweak(tw, "p_unsafe_funeral_hosp",
                  upper = list(mean = 0.010,  sd = 0.003),
                  lower = list(mean = 0.0005, sd = 0.001))

  # Assemble the Stan data list (must match the model's `data` block exactly).
  # Note q_obs and q_pred are DATA here -- the curve shape is given, not fit.
  stan_data <- list(
    N = nrow(fit_df), J = J, N_pred = length(tau_grid),
    param_id = fit_df$param_id, q_obs = fit_df$q_obs,         # shared Q at each anchor
    y_obs = fit_df$y_obs, obs_sd_mult = fit_df$obs_sd_mult,
    abs_min = param_meta$abs_min, abs_max = param_meta$abs_max,
    lower_floor = param_meta$lower_floor, lower_cap = param_meta$lower_cap,
    upper_cap = param_meta$upper_cap,
    lb_prior_mean = param_meta$lb_prior_mean, lb_prior_sd = param_meta$lb_prior_sd,
    ub_prior_mean = param_meta$ub_prior_mean, ub_prior_sd = param_meta$ub_prior_sd,
    increases = param_meta$increases,
    sigma_frac_prior_meanlog = SIGMA_FRAC_PRIOR_MEANLOG,
    sigma_frac_prior_sdlog   = SIGMA_FRAC_PRIOR_SDLOG,
    use_upper_tweak = tw$use_upper, upper_tweak_mean = tw$upper_mean, upper_tweak_sd = tw$upper_sd,
    use_lower_tweak = tw$use_lower, lower_tweak_mean = tw$lower_mean, lower_tweak_sd = tw$lower_sd,
    q_pred = q_pred                                            # shared Q on the prediction grid
  )

  # Draw posterior samples (adapt_delta slightly higher than the West Africa fit).
  fit <- mod$sample(
    data = stan_data, seed = 123,
    chains = 4, parallel_chains = 4,
    iter_warmup = 1500, iter_sampling = 1500,
    adapt_delta = 0.98, max_treedepth = 13, refresh = 200
  )
  cat("\nSampler diagnostics (", label, "):\n", sep = ""); print(fit$diagnostic_summary())

  # ---- 4. Per-parameter native-unit curves ---------------------------------
  # theta_pred[j, m] is parameter j's value at grid point m. Parse the [j,m]
  # indices out of the variable names and attach parameter / tau / day.
  curve_summ <- fit$summary(variables = "theta_pred") %>%
    mutate(
      param_id = as.integer(str_match(variable, "\\[([0-9]+),([0-9]+)\\]")[, 2]),
      grid_id  = as.integer(str_match(variable, "\\[([0-9]+),([0-9]+)\\]")[, 3]),
      tau          = tau_grid[grid_id],
      relative_day = tau * max_day
    ) %>%
    left_join(select(param_meta, param_id, parameter), by = "param_id") %>%
    select(parameter, param_id, tau, relative_day, mean, q5, q95)

  # Deterministic override for community unsafe funerals: replace the fitted
  # (relative-Q) curve with the ABSOLUTE Warsame proxy 1 - success, interpolated
  # onto the grid. (Its floor is 1 - max(success), not 0 -- see 00.)
  ufc_on_tau <- make_interp(qseries$tau_q, qseries$unsafe_funeral_comm_proxy)
  curve_summ <- curve_summ %>%
    mutate(mean = if_else(parameter == "p_unsafe_funeral_comm",
                          clip01(ufc_on_tau(tau)), mean))

  # The shared Q on the prediction grid; 03 uses this as the scenario's q_value.
  q_grid <- tibble(tau = tau_grid, relative_day = tau_grid * max_day, q_value = q_pred)

  # Return everything 03 needs for this scenario:
  #   curve_summ : per-parameter native-unit curve (parameter, tau, relative_day, mean, q5, q95)
  #   q_grid     : the shared Q over the grid (tau, relative_day, q_value)
  #   param_meta : per-parameter metadata; max_day : the scenario horizon in days
  list(curve_summ = curve_summ, q_grid = q_grid, param_meta = param_meta, max_day = max_day)
}

# ----------------------------------------------------------------------------
# 4. Run both conflict scenarios and save
# ----------------------------------------------------------------------------
# The two Q series differ only in that the "++" one has the forced collapse baked
# in (done in 00); the SAME Model B is fit to each.
drc_conflict_qseries          <- drc_prep$conflict_qseries
drc_conflict_plusplus_qseries <- drc_prep$conflict_plusplus_qseries

drc_conflict_fit          <- fit_drc_scenario(drc_conflict_qseries,          "drc_conflict")
drc_conflict_plusplus_fit <- fit_drc_scenario(drc_conflict_plusplus_qseries, "drc_conflict_plusplus")

saveRDS(drc_conflict_fit,          file.path(DIR_PROCESSED, "drc_conflict_fit.rds"))
saveRDS(drc_conflict_plusplus_fit, file.path(DIR_PROCESSED, "drc_conflict_plusplus_fit.rds"))

# ----------------------------------------------------------------------------
# 5. Plot each fitted scenario with its data on top  (display only)
# ----------------------------------------------------------------------------
# For a fitted scenario, draw one facet per parameter: the posterior-MEAN curve
# (blue) with its 90% interval, the literature anchors the model was fit to
# (orange points), and -- in the community unsafe-funeral panel only -- the SDB
# data points (1 - success) that actually drive that curve (grey). The plots are
# printed to the active graphics device and are deliberately NOT saved.
plot_drc_fit <- function(fit, qseries, label) {

  # The community unsafe-funeral mean was deterministically overridden to the
  # absolute Warsame proxy, so its fitted q5/q95 are now stale; collapse the
  # ribbon there to avoid drawing a misleading interval in that one panel.
  curve_plot_df <- fit$curve_summ %>%
    mutate(q5  = if_else(parameter == "p_unsafe_funeral_comm", mean, q5),
           q95 = if_else(parameter == "p_unsafe_funeral_comm", mean, q95),
           panel = factor(PANEL_LOOKUP[parameter], levels = unname(PANEL_LOOKUP)))

  # Literature anchors the model was fit to (one set, shared across scenarios).
  anchor_plot_df <- drc_anchors %>%
    mutate(panel = factor(PANEL_LOOKUP[parameter], levels = unname(PANEL_LOOKUP)))

  # SDB community unsafe-funeral data points (1 - success) behind the UFC curve.
  sdb_plot_df <- qseries %>%
    filter(n_eligible_sum > 0) %>%
    transmute(relative_day,
              value = unsafe_funeral_comm_proxy,
              panel = factor(PANEL_LOOKUP[["p_unsafe_funeral_comm"]], levels = unname(PANEL_LOOKUP)))

  ggplot(curve_plot_df, aes(relative_day, mean)) +
    geom_ribbon(aes(ymin = q5, ymax = q95), fill = "#1f77b4", alpha = 0.20) +
    geom_line(colour = "#1f77b4", linewidth = 0.9) +
    geom_point(data = sdb_plot_df, aes(relative_day, value),
               inherit.aes = FALSE, colour = "grey55", size = 1, alpha = 0.7) +
    geom_point(data = anchor_plot_df, aes(relative_day, value_used),
               inherit.aes = FALSE, colour = "#ff7f0e", size = 2) +
    facet_wrap(~ panel, scales = "free_y", ncol = 2) +
    labs(title = paste0("DRC ", label, ": fitted parameter curves vs data"),
         subtitle = "Blue = posterior mean + 90% interval; orange = literature anchors; grey = SDB community proxy",
         x = "Relative outbreak day", y = NULL) +
    theme_bw(base_size = 11) +
    theme(strip.text = element_text(face = "bold"))
}

print(plot_drc_fit(drc_conflict_fit,          drc_conflict_qseries,          "conflict"))     # display only
print(plot_drc_fit(drc_conflict_plusplus_fit, drc_conflict_plusplus_qseries, "conflict++"))   # display only

message("\n02_DRC_QCurve_Fitting_Original.R complete. Saved DRC conflict + conflict++ fits.")
