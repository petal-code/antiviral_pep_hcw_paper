# ============================================================
# DRC middle scenario: conflict-including specification
#
# Charlie's requested conflict scenario:
#   - DO NOT bridge conflict.
#   - Use the monthly average successful-response trajectory, then smooth that
#     black average line before constructing Q.
#   - In the p_unsafe_funeral_comm panel, show the empirical Warsame-derived
#     unsafe-funeral proxy on the absolute scale, i.e. 1 - successful_response(t).
#   - Use the smoothed black average SUCCESS line as the empirical backbone,
#     augment it with a day-0 anchor at Q = 0, and scale it by its own maximum
#     among finite-date Warsame support points only so that Q is a relative 0-1 response index.
#
#         Q_conflict(t) = success_avg_smoothed(t) / max(success_avg_smoothed), with Q(0) = 0
#
#     This is deliberately NOT min-max scaling: late empirical dips are not
#     redefined as Q = 0 unless successful response truly returns to zero.
#
#     Important orientation fix: Q_conflict is an improvement / response-capacity
#     signal, so it must be based on successful response, NOT 1 - successful response.
#     The 1 - transformation is used only when plotting the p_unsafe_funeral_comm
#     reference line, because that panel is on the unsafe-funeral-risk scale.
#     For this visual reference, use 1 - successful_response(t), not 1 - Q_conflict(t),
#     so the best empirical SDB response bottoms out at 1 - max(successful_response)
#     rather than being artificially forced to 0.
#
#   - Use this shared empirical smoothed Q_conflict(t) as the common shape for
#     ALL DRC parameters, so they all inherit the same broad conflict interruptions.
#   - Fit parameter-specific lower/upper bounds (independent, prior-anchored)
#     against the workbook anchors.
#
# Important:
#   - This is NOT the same as the earlier DRC scripts where the Warsame data
#     were inserted as direct y-observations for p_unsafe_funeral_comm.
#   - Here, the Warsame monthly average defines the SHARED Q shape directly.
#   - The workbook anchors are then fitted against that shared Q shape.
#
# Default evidence handling:
#   - Keeps the curated workbook anchors in the fit.
#   - Removes DRC_UFC_01 by default because it is a summary of the Warsame
#     evidence stream and would otherwise double count it.
#
# ============================================================

suppressPackageStartupMessages({
  library(cmdstanr)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(tibble)
  library(readr)
})

# ------------------------------------------------------------
# User paths
# ------------------------------------------------------------
curve_workbook <- find_input_file(c("filovirus_three_scenario_curve_inputs_bestcase_recreate.xlsx"), "original methodology curve workbook")
curve_sheet    <- "Middle_DRC_2018_2019"

sdb_workbook <- find_input_file(c("evd_drc_sdb_performance_datasets_pub.xlsx", "evd_drc_sdb_performance_datasets_pub(1).xlsx", "evd_drc_sdb_performance_datasets_pub (1).xlsx"), "Warsame SDB workbook")
if (!file.exists(sdb_workbook)) {
  sdb_workbook_alt <- "evd_drc_sdb_performance_datasets_pub(1).xlsx"
  if (file.exists(sdb_workbook_alt)) {
    sdb_workbook <- sdb_workbook_alt
  }
}
if (!file.exists(sdb_workbook)) {
  sdb_workbook_alt2 <- "evd_drc_sdb_performance_datasets_pub (1).xlsx"
  if (file.exists(sdb_workbook_alt2)) {
    sdb_workbook <- sdb_workbook_alt2
  }
}

stan_file  <- "drc_conflict_warsame_weekly_Q_relative_absolute_UFC_rolling4_sparse_greydots.stan"
out_prefix <- "drc_conflict_plusplus_original_method_directCollapse"

# ------------------------------------------------------------
# Settings for reconstructing the Warsame orange line
# ------------------------------------------------------------
provinces_to_average <- c("North Kivu", "Ituri")

# Orange line denominator:
# Warsame Fig. 3 defines a successful response to include not only literal
# outcome_lshtm == "success", but also cases where an SDB was not ultimately
# needed after testing / assessment (the body could be released appropriately).
# Therefore the orange line numerator is:
#   success + "sdb not needed"
# and the denominator is:
#   success + "sdb not needed" + failure
# False alerts and unclear outcomes are excluded from the empirical performance
# denominator.
success_denominator_mode <- "successful_response_vs_failure"

# Primary choice remains the unweighted province mean, i.e. average the
# province-specific black lines rather than weighting by alert volume.
average_method <- "unweighted_province_mean"  # or "count_weighted"

# Temporal aggregation
# Use epidemiological weeks to explicitly recreate the step-like Warsame curve,
# rather than monthly smoothing/aggregation.
aggregation_unit <- "epi_week"  # "epi_week", "monthly", or "four_week"
bin_width_days <- 28            # only used when aggregation_unit == "four_week"

# Primary version uses all SDB responses, consistent with the published line.
sdb_origin_filter <- "all"     # "all" or "community"

# Early success spikes can still be tiny-n artefacts.
initial_small_n_action <- "zero_success"  # "zero_success", "drop", or "none"
initial_spike_max_day <- 75
initial_spike_min_eligible <- 10
initial_spike_success_threshold <- 0.50

# Minimum eligible SDB count required for a bin to contribute to the
# shared Q_conflict shape. With weekly reconstruction, keep this low so the
# plotted line explicitly follows the Warsame curve.
min_sdb_eligible_for_q <- 1

# Smooth the empirical weekly Warsame reconstruction using a 4-week rolling
# average. This keeps the empirical week-by-week reconstruction but removes
# some of the visual/parameter wiggle from individual noisy weeks.
#
# The default is centred because this is a retrospective reconstruction, not
# an online/real-time response rule. Switch to "trailing" if you want the
# curve to use only information available up to that week.
q_smoothing_method <- "rolling4"   # "none", "rolling4", "loess", or "spline"
rolling_window_n <- 4
rolling_align <- "center"          # "center" or "trailing"
rolling_use_eligible_weights <- FALSE  # FALSE = simple rolling mean of reconstructed points

# Plotting only: show a sparse set of grey rolling-average reference dots
# in the p_unsafe_funeral_comm panel. The fitted curve still uses the full
# weekly rolling-average trajectory; this just avoids visually overcrowding
# the panel with one dot for every epidemiological week. For weekly data and
# rolling_window_n = 4, this shows approximately one reference dot per 4 weeks.
plot_sparse_rolling_ufc_dots <- TRUE
plot_rolling_reference_every_n <- rolling_window_n
plot_always_include_start_end_ufc_dots <- TRUE

# ------------------------------------------------------------
# DRC with conflict ++ sensitivity
# ------------------------------------------------------------
# This is the old/original-methodology version of the DRC++ scenario.
# We retain the original shared empirical Q framework, but impose a temporary
# complete response collapse across the Warsame rolling support points between
# relative day 200 and 300.
#
# Mechanistically:
#   - Successful SDB response is set to 0 during the plateau.
#   - Therefore p_unsafe_funeral_comm = 1 during the plateau.
#   - Since Q_conflict is based on successful response, Q drops to 0 during
#     the plateau.
#   - All other parameters move accordingly through the original shared-Q
#     mapping: increasing response parameters are pulled down; decreasing
#     adverse parameters are pulled up.
plusplus_window_day <- c(200, 300)
plusplus_success_value <- 0
plusplus_unsafe_funeral_value <- 1

q_loess_span <- 0.45
q_spline_df  <- 5

# ------------------------------------------------------------
# Workbook anchor handling
# ------------------------------------------------------------
keep_workbook_ufc_anchors <- TRUE

# DRC_UFC_00 conflicts with the requirement that p_unsafe_funeral_comm
# should start at 1, and DRC_UFC_01 is a Warsame summary that would double
# count the line-list evidence. Both are therefore reference-only by default.
remove_superseded_ufc_anchor_ids <- c("DRC_UFC_00", "DRC_UFC_01")
plot_removed_ufc_anchors_as_reference <- TRUE

# ------------------------------------------------------------
# Parameter setup
# ------------------------------------------------------------
param_levels <- c(
  "delay_hosp",
  "p_hosp",
  "p_ETU",
  "latent_IPC",
  "p_unsafe_funeral_comm",
  "p_unsafe_funeral_hosp"
)

panel_lookup <- c(
  "delay_hosp"             = "Delay to hospitalisation",
  "p_hosp"                 = "Probability of hospitalisation",
  "p_ETU"                  = "Proportion in ETU / ETC",
  "latent_IPC"             = "Latent IPC index",
  "p_unsafe_funeral_comm"  = "Unsafe funeral after community death",
  "p_unsafe_funeral_hosp"  = "Unsafe funeral after hospital death"
)

panel_order <- c(
  "Delay to hospitalisation",
  "Probability of hospitalisation",
  "Proportion in ETU / ETC",
  "Latent IPC index",
  "Unsafe funeral after community death",
  "Unsafe funeral after hospital death"
)

# Hard admissible supports
domain_meta <- tibble::tribble(
  ~parameter,                 ~abs_min, ~abs_max,
  "delay_hosp",               0.5,      12.0,
  "p_hosp",                   0.0,      1.0,
  "p_ETU",                    0.0,      1.0,
  "latent_IPC",               0.0,      1.0,
  "p_unsafe_funeral_comm",    0.0,      1.0,
  "p_unsafe_funeral_hosp",    0.0,      0.010
)

# Lower/upper support around workbook prior centres
support_half_width_mult <- 0.75
bound_sd_frac_of_initial_span <- 0.15
bound_sd_frac_of_domain       <- 0.03

# Observation scale prior
sigma_frac_prior_meanlog <- log(0.12)
sigma_frac_prior_sdlog   <- 0.60

# Force p_unsafe_funeral_comm to begin at 1 when the shared Q starts at 0.
use_ufc_comm_start1_tweak <- TRUE
ufc_comm_upper_tweak_mean <- 1.00
ufc_comm_upper_tweak_sd   <- 0.02

# Important: Q_conflict is relative and max-scaled, but unsafe funerals are
# absolute. If the best empirical successful-response value is, for example, 0.875,
# the minimum absolute unsafe-funeral proxy should be 1 - 0.875 = 0.125,
# NOT zero. This tweak anchors the lower level of p_unsafe_funeral_comm to
# 1 - max(successful response).
use_ufc_comm_absolute_warsame_tweak <- TRUE
ufc_comm_lower_warsame_sd <- 0.02

# Strong regularisation for p_unsafe_funeral_hosp only
use_ufh_tweaks <- TRUE
ufh_upper_tweak_mean <- 0.010
ufh_upper_tweak_sd   <- 0.003
ufh_lower_tweak_mean <- 0.0005
ufh_lower_tweak_sd   <- 0.001

# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------
clip01 <- function(x) pmin(1, pmax(0, x))

rolling_mean_ordered <- function(y, k = 4, align = "center", weights = NULL) {
  if (!align %in% c("center", "trailing")) {
    stop("rolling_align must be either 'center' or 'trailing'.")
  }
  if (length(y) == 0) return(numeric(0))
  if (k < 1) stop("Rolling window k must be >= 1.")

  n <- length(y)
  out <- rep(NA_real_, n)

  if (is.null(weights)) {
    weights <- rep(1, n)
  }

  for (i in seq_len(n)) {
    if (align == "trailing") {
      idx <- seq.int(max(1, i - k + 1), i)
    } else {
      # Even-width centred windows are necessarily slightly asymmetric.
      # For k = 4 this uses one point before, the current point, and two after.
      left_n <- floor((k - 1) / 2)
      right_n <- (k - 1) - left_n
      idx <- seq.int(max(1, i - left_n), min(n, i + right_n))
    }

    yy <- y[idx]
    ww <- weights[idx]
    ok <- is.finite(yy) & is.finite(ww) & ww > 0

    if (any(ok)) {
      out[i] <- stats::weighted.mean(yy[ok], w = ww[ok])
    }
  }

  out
}


smooth_shared_q_source <- function(df, x_col, y_col,
                                   method = "loess",
                                   rolling_window_n = 4,
                                   rolling_align = "center",
                                   rolling_use_eligible_weights = FALSE,
                                   loess_span = 0.45,
                                   spline_df = 5) {
  out <- df %>% arrange(.data[[x_col]])
  x <- out[[x_col]]
  y <- out[[y_col]]

  if (method == "none") {
    out$success_avg_smoothed <- y
    return(out)
  }

  if (method == "rolling4") {
    out$success_avg_smoothed <- clip01(
      rolling_mean_ordered(
        y = y,
        k = rolling_window_n,
        align = rolling_align,
        weights = if (rolling_use_eligible_weights) out$n_eligible_sum else NULL
      )
    )
    out$success_avg_smoothed[!is.finite(out$success_avg_smoothed)] <- y[!is.finite(out$success_avg_smoothed)]
    return(out)
  }

  if (length(unique(x)) < 4) {
    warning("Too few unique x values to smooth robustly; returning unsmoothed Warsame average.")
    out$success_avg_smoothed <- y
    return(out)
  }

  if (method == "loess") {
    fit <- loess(
      formula = y ~ x,
      span = loess_span,
      degree = 2,
      family = "gaussian",
      control = loess.control(surface = "direct")
    )
    yhat <- as.numeric(predict(fit, newdata = data.frame(x = x)))
  } else if (method == "spline") {
    fit <- smooth.spline(x = x, y = y, df = spline_df)
    yhat <- as.numeric(predict(fit, x = x)$y)
  } else {
    stop("q_smoothing_method must be one of: 'none', 'rolling4', 'loess', or 'spline'.")
  }

  yhat[!is.finite(yhat)] <- y[!is.finite(yhat)]
  out$success_avg_smoothed <- clip01(yhat)
  out
}

# ------------------------------------------------------------
# Stan model: shared empirical Q_conflict(t), independent lower/upper bounds
# ------------------------------------------------------------
stan_code <- '
data {
  int<lower=1> N;
  int<lower=1> J;
  int<lower=1> N_pred;

  array[N] int<lower=1, upper=J> param_id;

  vector<lower=0, upper=1>[N] q_obs;
  vector[N] y_obs;
  vector<lower=0>[N] obs_sd_mult;

  vector[J] abs_min;
  vector[J] abs_max;
  vector[J] lower_floor;
  vector[J] lower_cap;
  vector[J] upper_cap;
  vector[J] lb_prior_mean;
  vector[J] lb_prior_sd;
  vector[J] ub_prior_mean;
  vector[J] ub_prior_sd;
  array[J] int<lower=0, upper=1> increases;

  real sigma_frac_prior_meanlog;
  real<lower=0> sigma_frac_prior_sdlog;

  array[J] int<lower=0, upper=1> use_upper_tweak;
  vector[J] upper_tweak_mean;
  vector<lower=0>[J] upper_tweak_sd;

  array[J] int<lower=0, upper=1> use_lower_tweak;
  vector[J] lower_tweak_mean;
  vector<lower=0>[J] lower_tweak_sd;

  vector<lower=0, upper=1>[N_pred] q_pred;
}

parameters {
  vector[J] lower_raw;
  vector[J] upper_gap_raw;
  vector<lower=0>[J] sigma_frac;
}

transformed parameters {
  vector[J] lower_est;
  vector[J] upper_est;
  vector[J] span_est;

  for (j in 1:J) {
    real s_l;
    real s_u;
    real lower_j;
    real upper_j;

    s_l = inv_logit(lower_raw[j]);
    lower_j = lower_floor[j] + (lower_cap[j] - lower_floor[j]) * s_l;

    s_u = inv_logit(upper_gap_raw[j]);
    upper_j = lower_j + (upper_cap[j] - lower_j) * s_u;

    lower_est[j] = lower_j;
    upper_est[j] = upper_j;
    span_est[j]  = upper_j - lower_j;
  }
}

model {
  sigma_frac ~ lognormal(sigma_frac_prior_meanlog, sigma_frac_prior_sdlog);

  for (j in 1:J) {
    real s_l = inv_logit(lower_raw[j]);
    real s_u = inv_logit(upper_gap_raw[j]);

    target += normal_lpdf(lower_est[j] | lb_prior_mean[j], lb_prior_sd[j]);
    target += log(lower_cap[j] - lower_floor[j]) + log(s_l) + log1m(s_l);

    target += normal_lpdf(upper_est[j] | ub_prior_mean[j], ub_prior_sd[j]);
    target += log(upper_cap[j] - lower_est[j]) + log(s_u) + log1m(s_u);
  }

  for (j in 1:J) {
    if (use_upper_tweak[j] == 1) {
      target += normal_lpdf(upper_est[j] | upper_tweak_mean[j], upper_tweak_sd[j]);
    }
    if (use_lower_tweak[j] == 1) {
      target += normal_lpdf(lower_est[j] | lower_tweak_mean[j], lower_tweak_sd[j]);
    }
  }

  for (n in 1:N) {
    int j = param_id[n];
    real mu;

    if (increases[j] == 1) {
      mu = lower_est[j] + span_est[j] * q_obs[n];
    } else {
      mu = upper_est[j] - span_est[j] * q_obs[n];
    }

    y_obs[n] ~ student_t(4, mu, sigma_frac[j] * span_est[j] * obs_sd_mult[n]);
  }
}

generated quantities {
  matrix[J, N_pred] Q_pred;
  matrix[J, N_pred] theta_pred;
  vector[N] log_lik;

  for (j in 1:J) {
    for (m in 1:N_pred) {
      real q = q_pred[m];
      Q_pred[j, m] = q;

      if (increases[j] == 1) {
        theta_pred[j, m] = lower_est[j] + span_est[j] * q;
      } else {
        theta_pred[j, m] = upper_est[j] - span_est[j] * q;
      }
    }
  }

  for (n in 1:N) {
    int j = param_id[n];
    real mu;

    if (increases[j] == 1) {
      mu = lower_est[j] + span_est[j] * q_obs[n];
    } else {
      mu = upper_est[j] - span_est[j] * q_obs[n];
    }

    log_lik[n] = student_t_lpdf(
      y_obs[n] | 4, mu, sigma_frac[j] * span_est[j] * obs_sd_mult[n]
    );
  }
}
'
writeLines(stan_code, con = stan_file)

# ------------------------------------------------------------
# Read the scenario evidence workbook
# ------------------------------------------------------------
anchors_raw <- read_excel(curve_workbook, sheet = curve_sheet, skip = 2)

required_cols <- c(
  "anchor_id", "parameter", "relative_day", "value_used",
  "fit_role", "include_in_fit", "weight", "direction",
  "lower_bound", "upper_bound"
)

missing_cols <- setdiff(required_cols, names(anchors_raw))
if (length(missing_cols) > 0) {
  stop("Missing required workbook columns: ", paste(missing_cols, collapse = ", "))
}

anchors <- anchors_raw %>%
  mutate(
    anchor_id      = as.character(anchor_id),
    parameter      = as.character(parameter),
    fit_role       = trimws(as.character(fit_role)),
    include_in_fit = toupper(trimws(as.character(include_in_fit))),
    direction      = tolower(trimws(as.character(direction))),
    relative_day   = as.numeric(relative_day),
    value_used     = as.numeric(value_used),
    lower_bound    = as.numeric(lower_bound),
    upper_bound    = as.numeric(upper_bound),
    weight         = as.numeric(weight),
    source_type    = "workbook_anchor"
  )

# ------------------------------------------------------------
# Build monthly Warsame orange-line values by province, with NO bridging
# and NO no-conflict plateau. This is the conflict-including version.
# ------------------------------------------------------------
admin_units <- read_excel(sdb_workbook, sheet = "admin_units")
sdb_dataset <- read_excel(sdb_workbook, sheet = "sdb_dataset")

admin_units <- admin_units %>%
  mutate(
    hz = as.character(hz),
    province = as.character(province)
  )

sdb_dataset <- sdb_dataset %>%
  mutate(
    hz = as.character(hz),
    outcome_lshtm = tolower(trimws(as.character(outcome_lshtm))),
    origin_cat = tolower(trimws(as.character(origin_cat))),
    date = as.Date(date)
  )

province_lookup <- admin_units %>%
  distinct(hz, province)

sdb_prov <- sdb_dataset %>%
  left_join(province_lookup, by = "hz") %>%
  filter(province %in% provinces_to_average) %>%
  # Important: some line-list rows have epi_year / epi_week but no actual date.
  # If these undated rows are retained, they can form NA-date weekly bins with
  # prop_success = 1 and therefore become the empirical maximum used to scale Q,
  # while never appearing on the date-based plots or prediction grid. That was
  # why the plotted Q curve topped out below 1 even though the internal support
  # check could still pass. Drop undated rows before reconstructing the Warsame
  # curve so that the maximum used for Q is the maximum of the plotted orange line.
  filter(!is.na(date))

if (sdb_origin_filter == "community") {
  sdb_prov <- sdb_prov %>%
    filter(origin_cat == "community")
}

if (nrow(sdb_prov) == 0) {
  stop("No SDB line-list rows found after filtering to requested provinces / origin setting.")
}

sdb_start_date <- min(sdb_prov$date, na.rm = TRUE)

if (aggregation_unit == "epi_week") {
  sdb_prov <- sdb_prov %>%
    mutate(
      # Monday-start week, matching the epi-week grouping used in the workbook.
      # This keeps bin dates identical across provinces before averaging.
      bin_start = date - as.integer(format(date, "%u")) + 1,
      bin_mid_date = bin_start + 3,
      bin_label = paste0(epi_year, "-W", sprintf("%02d", as.integer(epi_week)))
    )
} else if (aggregation_unit == "monthly") {
  sdb_prov <- sdb_prov %>%
    mutate(
      bin_start = as.Date(format(date, "%Y-%m-01")),
      bin_mid_date = bin_start + 14,
      bin_label = format(bin_start, "%Y-%m")
    )
} else if (aggregation_unit == "four_week") {
  sdb_prov <- sdb_prov %>%
    mutate(
      days_from_start = as.numeric(date - sdb_start_date),
      bin_index = floor(days_from_start / bin_width_days),
      bin_start = sdb_start_date + bin_index * bin_width_days,
      bin_mid_date = bin_start + floor(bin_width_days / 2),
      bin_label = paste0("bin_", bin_index)
    )
} else {
  stop("aggregation_unit must be one of: 'epi_week', 'monthly', or 'four_week'.")
}

prov_binned <- sdb_prov %>%
  group_by(province, bin_label, bin_start, bin_mid_date) %>%
  summarise(
    n_success = sum(outcome_lshtm == "success", na.rm = TRUE),
    n_successful_response = sum(outcome_lshtm %in% c("success", "sdb not needed"), na.rm = TRUE),
    n_failure = sum(outcome_lshtm == "failure", na.rm = TRUE),
    n_not_needed = sum(outcome_lshtm == "sdb not needed", na.rm = TRUE),
    n_false_alert = sum(outcome_lshtm == "false alert", na.rm = TRUE),
    n_unclear = sum(outcome_lshtm == "unclear", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    n_eligible = case_when(
      success_denominator_mode == "successful_response_vs_failure" ~ n_successful_response + n_failure,
      # Retained only as a sensitivity option. Because n_successful_response
      # already includes "sdb not needed", do not add n_not_needed again here.
      success_denominator_mode == "successful_response_all_nonmissing" ~ n_successful_response + n_failure + n_false_alert,
      TRUE ~ NA_real_
    ),
    prop_success_raw = if_else(n_eligible > 0, n_successful_response / n_eligible, NA_real_),
    relative_day = as.numeric(bin_mid_date - sdb_start_date),
    small_n_initial_spike = initial_small_n_action != "none" &
      relative_day <= initial_spike_max_day &
      n_eligible < initial_spike_min_eligible &
      prop_success_raw >= initial_spike_success_threshold,
    prop_success_clean = case_when(
      initial_small_n_action == "zero_success" & small_n_initial_spike ~ 0,
      initial_small_n_action == "drop" & small_n_initial_spike ~ NA_real_,
      TRUE ~ prop_success_raw
    )
  ) %>%
  filter(!is.na(prop_success_clean)) %>%
  arrange(province, bin_mid_date)

sdb_avg <- prov_binned %>%
  group_by(bin_label, bin_start, bin_mid_date) %>%
  summarise(
    relative_day = as.numeric(first(bin_mid_date) - sdb_start_date),
    n_provinces_with_data = n_distinct(province),
    n_eligible_sum = sum(n_eligible, na.rm = TRUE),
    prop_success_unweighted = mean(prop_success_clean, na.rm = TRUE),
    prop_success_count_weighted = if_else(
      sum(n_eligible, na.rm = TRUE) > 0,
      weighted.mean(prop_success_clean, w = pmax(n_eligible, 0), na.rm = TRUE),
      NA_real_
    ),
    any_initial_spike_adjusted = any(small_n_initial_spike),
    .groups = "drop"
  ) %>%
  arrange(bin_mid_date) %>%
  mutate(
    prop_success_avg = case_when(
      average_method == "unweighted_province_mean" ~ prop_success_unweighted,
      average_method == "count_weighted" ~ prop_success_count_weighted,
      TRUE ~ NA_real_
    )
  )

if (any(is.na(sdb_avg$prop_success_avg))) {
  stop("Some averaged SDB success values are NA; check average_method / denominator settings.")
}

# Keep bins sufficiently data-rich to define the conflict Q shape.
q_source <- sdb_avg %>%
  filter(n_eligible_sum >= min_sdb_eligible_for_q) %>%
  filter(is.finite(relative_day), !is.na(bin_mid_date), is.finite(prop_success_avg)) %>%
  arrange(bin_mid_date)

if (nrow(q_source) < 3) {
  stop("Too few monthly SDB bins remain after min_sdb_eligible_for_q filter.")
}

q_source <- smooth_shared_q_source(
  df = q_source,
  x_col = "relative_day",
  y_col = "prop_success_avg",
  method = q_smoothing_method,
  rolling_window_n = rolling_window_n,
  rolling_align = rolling_align,
  rolling_use_eligible_weights = rolling_use_eligible_weights,
  loess_span = q_loess_span,
  spline_df = q_spline_df
)

# ------------------------------------------------------------
# DRC++ temporary response-collapse plateau
# ------------------------------------------------------------
# The ++ scenario sets the Warsame rolling successful-response signal to zero
# from the first Warsame support point at/after day 200 through the last support
# point at/before day 300. This affects the shared Q_conflict(t) itself and
# therefore all response parameters under the original methodology.
q_source_before_plusplus <- q_source

plusplus_support_points <- q_source %>%
  dplyr::filter(is.finite(relative_day),
                relative_day >= plusplus_window_day[1],
                relative_day <= plusplus_window_day[2]) %>%
  dplyr::arrange(relative_day)

if (nrow(plusplus_support_points) == 0) {
  warning("No Warsame Q support points found inside the DRC++ window; no ++ plateau applied.")
  plusplus_plateau_start_day <- NA_real_
  plusplus_plateau_end_day <- NA_real_
} else {
  plusplus_plateau_start_day <- min(plusplus_support_points$relative_day, na.rm = TRUE)
  plusplus_plateau_end_day <- max(plusplus_support_points$relative_day, na.rm = TRUE)
}

q_source <- q_source %>%
  dplyr::mutate(
    is_drc_plusplus_collapse_point =
      is.finite(relative_day) &
      is.finite(plusplus_plateau_start_day) &
      relative_day >= plusplus_plateau_start_day &
      relative_day <= plusplus_plateau_end_day,
    prop_success_avg_original = prop_success_avg,
    success_avg_smoothed_original = success_avg_smoothed,
    prop_success_avg = dplyr::if_else(
      is_drc_plusplus_collapse_point,
      plusplus_success_value,
      prop_success_avg
    ),
    success_avg_smoothed = dplyr::if_else(
      is_drc_plusplus_collapse_point,
      plusplus_success_value,
      success_avg_smoothed
    )
  )

plusplus_points_forced <- q_source %>%
  dplyr::filter(is_drc_plusplus_collapse_point) %>%
  dplyr::mutate(
    forced_success_avg_smoothed = success_avg_smoothed,
    forced_unsafe_funeral_comm_proxy_abs = plusplus_unsafe_funeral_value
  )

readr::write_csv(
  plusplus_points_forced,
  paste0(out_prefix, "_plusplus_points_forced_to_collapse.csv")
)

# ------------------------------------------------------------
# Convert the empirical black average SUCCESS line into the shared Q_conflict shape.
# Do not invert this with 1 - success here. Q is the response-improvement signal.
#
# Important: Q_conflict is a RELATIVE 0-1 response index.
#
# We therefore scale the empirical smoothed successful-response trajectory by its
# own maximum only:
#
#   Q_conflict(t) = success_avg_smoothed(t) / max(success_avg_smoothed)
#
# after augmenting the series with Q(0) = 0.
#
# This is NOT min-max rescaling. We do not subtract the empirical minimum,
# because doing so made late low SDB values become Q ~= 0 and implied that
# the response had returned to outbreak-day-zero failure. Here, the best
# achieved empirical response is Q = 1, while later dips remain partial dips
# unless the empirical success proportion truly returns to zero.
#
# However, because Q is relative, the unsafe-funeral probability should NOT
# be plotted as 1 - Q. It should be plotted on the absolute Warsame scale:
#
#   p_unsafe_funeral_comm proxy = 1 - success_avg_smoothed
#
# Thus if the best achieved successful-response value is 0.875, the minimum
# unsafe-funeral proxy is 0.125, not 0.
# ------------------------------------------------------------
# Global x-axis horizon is shared between workbook anchors and the SDB series.
max_day <- max(
  max(anchors$relative_day, na.rm = TRUE),
  max(q_source$relative_day, na.rm = TRUE)
)

if (!is.finite(max_day) || max_day <= 0) {
  stop("Could not construct a valid global DRC time horizon.")
}

# Add the forced outbreak-start anchor on the ABSOLUTE successful-response scale.
# Then compute Q from that same empirical successful-response line. This keeps
# Q relative (0-1) while keeping the unsafe-funeral curve absolute.
q_source <- bind_rows(
  tibble(
    bin_label = "START",
    bin_start = sdb_start_date,
    bin_mid_date = sdb_start_date,
    relative_day = 0,
    n_provinces_with_data = 0,
    n_eligible_sum = 0,
    prop_success_unweighted = NA_real_,
    prop_success_count_weighted = NA_real_,
    any_initial_spike_adjusted = FALSE,
    prop_success_avg = NA_real_,
    success_avg_smoothed = 0,
    is_drc_plusplus_collapse_point = FALSE,
    prop_success_avg_original = NA_real_,
    success_avg_smoothed_original = NA_real_
  ),
  q_source
) %>%
  arrange(relative_day, bin_mid_date) %>%
  distinct(relative_day, .keep_all = TRUE)

q_scale_max <- max(q_source$success_avg_smoothed, na.rm = TRUE)
if (!is.finite(q_scale_max) || q_scale_max <= 0) {
  stop("Could not construct Q_conflict: maximum successful-response value is not positive.")
}

q_source <- q_source %>%
  mutate(
    tau_q = relative_day / max_day,
    # Q is relative: the best empirical successful-response value achieved in
    # this reconstruction is Q = 1. There is no min-subtraction and no later
    # prediction-grid re-scaling.
    q_conflict_shape = clip01(success_avg_smoothed / q_scale_max),
    # The community unsafe-funeral proxy is absolute, not relative.
    # Therefore it is 1 - the Warsame successful-response curve, not 1 - Q.
    unsafe_funeral_comm_proxy_abs = clip01(1 - success_avg_smoothed),
    is_drc_plusplus_collapse_point = dplyr::coalesce(is_drc_plusplus_collapse_point, FALSE)
  )

if (abs(min(q_source$q_conflict_shape, na.rm = TRUE) - 0) > 1e-8 ||
    abs(max(q_source$q_conflict_shape, na.rm = TRUE) - 1) > 1e-8) {
  stop("Q_conflict support points should span exactly 0 to 1 after empirical max-scaling.")
}

q_final_max <- max(q_source$q_conflict_shape, na.rm = TRUE)

# ------------------------------------------------------------
# Workbook anchor handling
# ------------------------------------------------------------
anchors_for_fit <- anchors %>%
  filter(parameter %in% param_levels) %>%
  filter(!(anchor_id %in% remove_superseded_ufc_anchor_ids)) %>%
  filter(!(parameter == "p_unsafe_funeral_comm" & !keep_workbook_ufc_anchors))

removed_reference_points <- anchors %>%
  filter(anchor_id %in% remove_superseded_ufc_anchor_ids) %>%
  filter(parameter %in% param_levels) %>%
  filter(!is.na(value_used), !is.na(relative_day))

# ------------------------------------------------------------
# QC plots for the conflict-derived shared Q shape
# ------------------------------------------------------------
avg_success_raw_plot <- sdb_avg %>%
  mutate(
    province = "Average raw",
    value = prop_success_avg
  )

avg_success_smooth_plot <- q_source %>%
  transmute(
    bin_mid_date,
    value = success_avg_smoothed
  )

p_sdb_qc <- ggplot() +
  geom_line(
    data = prov_binned,
    aes(x = bin_mid_date, y = prop_success_clean, colour = province),
    linewidth = 0.8,
    alpha = 0.8
  ) +
  geom_point(
    data = prov_binned,
    aes(x = bin_mid_date, y = prop_success_clean, colour = province, size = n_eligible),
    alpha = 0.6
  ) +
  geom_line(
    data = avg_success_raw_plot,
    aes(x = bin_mid_date, y = value),
    linewidth = 1.0,
    colour = "black",
    alpha = 0.5,
    linetype = 3
  ) +
  geom_line(
    data = avg_success_smooth_plot,
    aes(x = bin_mid_date, y = value),
    linewidth = 1.5,
    colour = "black"
  ) +
  geom_point(
    data = avg_success_raw_plot,
    aes(x = bin_mid_date, y = value),
    size = 2.0,
    colour = "black",
    alpha = 0.5
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw(base_size = 11) +
  labs(
    title = "DRC Warsame SDB orange-line reconstruction: weekly successful-response curve",
    subtitle = paste0("Aggregation: ", aggregation_unit,
                      "; province average: ", average_method,
                      "; smoothing: ", q_smoothing_method),
    x = NULL,
    y = "Successful SDB proportion",
    colour = "Province",
    size = "Eligible SDBs"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

print(p_sdb_qc)

ggsave(
  filename = paste0(out_prefix, "_warsame_orange_line_qc.png"),
  plot = p_sdb_qc,
  width = 10,
  height = 6,
  dpi = 180
)

p_q_conflict_qc <- ggplot(q_source, aes(x = bin_mid_date, y = q_conflict_shape)) +
  geom_line(linewidth = 1.2, colour = "steelblue4") +
  geom_point(aes(size = n_eligible_sum), colour = "steelblue4", alpha = 0.85) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_bw(base_size = 11) +
  labs(
    title = "DRC shared conflict Q shape",
    subtitle = paste0("Empirical averaged successful-response line, where success includes SDB not needed, empirically max-scaled using the observed maximum; aggregation=", aggregation_unit, "; smoothing=", q_smoothing_method),
    x = NULL,
    y = "Q_conflict(t)",
    size = "Eligible SDBs"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

print(p_q_conflict_qc)

ggsave(
  filename = paste0(out_prefix, "_shared_conflict_Q_qc.png"),
  plot = p_q_conflict_qc,
  width = 10,
  height = 6,
  dpi = 180
)

write_csv(prov_binned, paste0(out_prefix, "_province_binned_sdb_success.csv"))
write_csv(sdb_avg, paste0(out_prefix, "_averaged_conflict_binned_sdb_success.csv"))
write_csv(q_source, paste0(out_prefix, "_shared_conflict_Q_points.csv"))
write_csv(q_source_before_plusplus, paste0(out_prefix, "_shared_conflict_Q_points_before_plusplus.csv"))

# ------------------------------------------------------------
# Parameter metadata from workbook priors only
# ------------------------------------------------------------
param_meta <- anchors %>%
  filter(parameter %in% param_levels) %>%
  filter(include_in_fit == "YES") %>%
  group_by(parameter) %>%
  summarise(
    lb_prior_mean = first(lower_bound),
    ub_prior_mean = first(upper_bound),
    direction     = first(direction),
    .groups = "drop"
  ) %>%
  left_join(domain_meta, by = "parameter") %>%
  mutate(
    param_id   = match(parameter, param_levels),
    increases  = if_else(direction == "up", 1L, 0L),
    span0      = ub_prior_mean - lb_prior_mean,
    domain_w   = abs_max - abs_min
  ) %>%
  arrange(param_id)

if (!all(param_levels %in% param_meta$parameter)) {
  missing_params <- setdiff(param_levels, param_meta$parameter)
  stop("Missing fitted parameters in workbook: ", paste(missing_params, collapse = ", "))
}

if (any(param_meta$ub_prior_mean <= param_meta$lb_prior_mean)) {
  stop("At least one parameter has workbook upper_bound <= lower_bound.")
}

param_meta <- param_meta %>%
  mutate(
    lower_floor = pmax(abs_min, lb_prior_mean - support_half_width_mult * span0),
    lower_cap_raw = pmin(
      lb_prior_mean + support_half_width_mult * span0,
      ub_prior_mean - 0.05 * span0,
      abs_max - 0.05 * domain_w
    ),
    lower_cap = pmax(lower_floor + 1e-4, lower_cap_raw),
    upper_cap_raw = pmin(abs_max, ub_prior_mean + support_half_width_mult * span0),
    upper_cap = pmax(lower_cap + 1e-4, upper_cap_raw),
    lb_prior_sd = pmax(
      bound_sd_frac_of_initial_span * span0,
      bound_sd_frac_of_domain * domain_w,
      1e-4
    ),
    ub_prior_sd = pmax(
      bound_sd_frac_of_initial_span * span0,
      bound_sd_frac_of_domain * domain_w,
      1e-4
    )
  )

if (any(param_meta$lower_cap <= param_meta$lower_floor)) {
  stop("At least one parameter has lower_cap <= lower_floor.")
}
if (any(param_meta$upper_cap <= param_meta$lower_cap)) {
  stop("At least one parameter has upper_cap <= lower_cap.")
}

cat("\nParameter prior setup:\n")
print(
  param_meta %>%
    select(parameter, lb_prior_mean, ub_prior_mean, abs_min, abs_max,
           lower_floor, lower_cap, upper_cap, lb_prior_sd, ub_prior_sd)
)

# ------------------------------------------------------------
# Fit rows: workbook anchors only, mapped onto shared Q_conflict(t)
# ------------------------------------------------------------
fit_df <- anchors_for_fit %>%
  filter(parameter %in% param_levels) %>%
  filter(include_in_fit == "YES") %>%
  filter(!is.na(value_used)) %>%
  filter(!is.na(relative_day)) %>%
  left_join(param_meta %>% select(parameter, param_id, increases), by = "parameter") %>%
  mutate(
    tau = relative_day / max_day,
    q_obs = approx(
      x = q_source$tau_q,
      y = q_source$q_conflict_shape,
      xout = tau,
      method = "linear",
      rule = 2
    )$y,
    y_obs = value_used,
    weight = if_else(is.na(weight), 1.0, weight),
    obs_sd_mult = 1 / pmax(weight, 0.25),
    panel_title = factor(panel_lookup[parameter], levels = panel_order)
  ) %>%
  arrange(param_id, tau)

removed_reference_points <- removed_reference_points %>%
  left_join(param_meta %>% select(parameter, param_id, increases), by = "parameter") %>%
  mutate(
    tau = relative_day / max_day,
    q_obs = approx(
      x = q_source$tau_q,
      y = q_source$q_conflict_shape,
      xout = tau,
      method = "linear",
      rule = 2
    )$y,
    y_obs = value_used,
    panel_title = factor(panel_lookup[parameter], levels = panel_order),
    source_type = "removed_reference_anchor"
  )

if (any(is.na(fit_df$y_obs))) stop("Some y_obs values are NA.")
if (any(is.na(fit_df$obs_sd_mult))) stop("Some obs_sd_mult values are NA.")
if (any(is.na(fit_df$q_obs))) stop("Some q_obs values are NA.")

cat("\nRows per parameter used in Stan fit:\n")
print(table(fit_df$parameter))

cat("\nRows by source type:\n")
print(table(fit_df$source_type, useNA = "ifany"))

# ------------------------------------------------------------
# Targeted tweak vectors for p_unsafe_funeral_hosp only
# ------------------------------------------------------------
J <- nrow(param_meta)

use_upper_tweak_vec  <- rep(0L, J)
upper_tweak_mean_vec <- param_meta$ub_prior_mean
upper_tweak_sd_vec   <- rep(1.0, J)

use_lower_tweak_vec  <- rep(0L, J)
lower_tweak_mean_vec <- param_meta$lb_prior_mean
lower_tweak_sd_vec   <- rep(1.0, J)

if (use_ufh_tweaks) {
  j_ufh <- param_meta$param_id[param_meta$parameter == "p_unsafe_funeral_hosp"]
  if (length(j_ufh) != 1) stop("Could not uniquely identify p_unsafe_funeral_hosp.")

  use_upper_tweak_vec[j_ufh]  <- 1L
  upper_tweak_mean_vec[j_ufh] <- ufh_upper_tweak_mean
  upper_tweak_sd_vec[j_ufh]   <- ufh_upper_tweak_sd

  use_lower_tweak_vec[j_ufh]  <- 1L
  lower_tweak_mean_vec[j_ufh] <- ufh_lower_tweak_mean
  lower_tweak_sd_vec[j_ufh]   <- ufh_lower_tweak_sd
}

if (use_ufc_comm_start1_tweak || use_ufc_comm_absolute_warsame_tweak) {
  j_ufc <- param_meta$param_id[param_meta$parameter == "p_unsafe_funeral_comm"]
  if (length(j_ufc) != 1) stop("Could not uniquely identify p_unsafe_funeral_comm.")

  if (use_ufc_comm_start1_tweak) {
    use_upper_tweak_vec[j_ufc]  <- 1L
    upper_tweak_mean_vec[j_ufc] <- ufc_comm_upper_tweak_mean
    upper_tweak_sd_vec[j_ufc]   <- ufc_comm_upper_tweak_sd
  }

  if (use_ufc_comm_absolute_warsame_tweak) {
    # Since p_unsafe_funeral_comm decreases with Q and Q is scaled as
    # success / max(success), anchoring the lower level to 1 - max(success)
    # makes the fitted p_unsafe_funeral_comm trajectory approximate the
    # absolute empirical proxy 1 - success, rather than 1 - relative_Q.
    use_lower_tweak_vec[j_ufc]  <- 1L
    lower_tweak_mean_vec[j_ufc] <- clip01(1 - q_scale_max)
    lower_tweak_sd_vec[j_ufc]   <- ufc_comm_lower_warsame_sd
  }
}

# ------------------------------------------------------------
# Build Stan data
# ------------------------------------------------------------
# Build a prediction grid that includes BOTH a regular time grid and the exact
# empirical Q support points. This matters because the weekly Warsame curve can
# have narrow peaks. If the prediction grid misses the exact support day where
# Q = 1, the plotted Q curve can appear to top out below 1 even though the
# empirical support points themselves have been correctly scaled.
#
# The safest approach is to construct the grid as a data frame and carry the
# exact support-point Q values through directly, rather than relying only on
# interpolation and then checking whether interpolation happened to hit the
# maximum exactly.
q_pred_regular <- tibble(
  tau = seq(0, 1, length.out = 250)
) %>%
  mutate(
    q_pred = approx(
      x = q_source$tau_q,
      y = q_source$q_conflict_shape,
      xout = tau,
      method = "linear",
      rule = 2
    )$y
  )

# Include the exact empirical support points in the plotting/prediction grid.
# This means the grid naturally contains Q = 1 at the empirical maximum.
# Do not add artificial anchors or force max(q_pred) to 1 afterwards.
q_pred_support <- q_source %>%
  transmute(
    tau = tau_q,
    q_pred = q_conflict_shape
  )

q_pred_df <- bind_rows(
  q_pred_regular,
  q_pred_support
) %>%
  filter(is.finite(tau), is.finite(q_pred)) %>%
  mutate(
    tau = clip01(tau),
    q_pred = clip01(q_pred)
  ) %>%
  group_by(tau) %>%
  summarise(q_pred = max(q_pred, na.rm = TRUE), .groups = "drop") %>%
  arrange(tau)

tau_pred <- q_pred_df$tau
q_pred_vec <- q_pred_df$q_pred

# For p_unsafe_funeral_comm, do NOT reconstruct the absolute curve from the
# relative Q scale. The Q curve is deliberately scaled to 0-1 for the shared
# response index, but the community-unsafe-funeral curve must be on the
# absolute Warsame scale:
#
#   p_unsafe_funeral_comm(t) = 1 - successful_response(t)
#
# Using q_pred_vec here would reintroduce the exact muddle we are trying to
# avoid: the best observed Q value is 1 by construction, but the best observed
# successful-response proportion is less than 1. Therefore interpolate the
# absolute unsafe-funeral proxy directly from q_source.
unsafe_funeral_comm_pred_vec <- approx(
  x = q_source$tau_q,
  y = q_source$unsafe_funeral_comm_proxy_abs,
  xout = tau_pred,
  method = "linear",
  rule = 2
)$y %>% clip01()

if (abs(min(q_pred_vec, na.rm = TRUE) - 0) > 1e-8 ||
    abs(max(q_pred_vec, na.rm = TRUE) - 1) > 1e-8) {
  stop("The prediction-grid Q values still do not span 0 to 1 after adding exact finite-date support points. Check q_source$q_conflict_shape and missing dates.")
}

if (min(unsafe_funeral_comm_pred_vec, na.rm = TRUE) < -1e-8) {
  stop("The absolute Warsame p_unsafe_funeral_comm curve unexpectedly went below 0.")
}

stan_data <- list(
  N = nrow(fit_df),
  J = nrow(param_meta),
  N_pred = length(tau_pred),
  param_id = fit_df$param_id,
  q_obs = fit_df$q_obs,
  y_obs = fit_df$y_obs,
  obs_sd_mult = fit_df$obs_sd_mult,
  abs_min = param_meta$abs_min,
  abs_max = param_meta$abs_max,
  lower_floor = param_meta$lower_floor,
  lower_cap = param_meta$lower_cap,
  upper_cap = param_meta$upper_cap,
  lb_prior_mean = param_meta$lb_prior_mean,
  lb_prior_sd = param_meta$lb_prior_sd,
  ub_prior_mean = param_meta$ub_prior_mean,
  ub_prior_sd = param_meta$ub_prior_sd,
  increases = param_meta$increases,
  sigma_frac_prior_meanlog = sigma_frac_prior_meanlog,
  sigma_frac_prior_sdlog = sigma_frac_prior_sdlog,
  use_upper_tweak = use_upper_tweak_vec,
  upper_tweak_mean = upper_tweak_mean_vec,
  upper_tweak_sd = upper_tweak_sd_vec,
  use_lower_tweak = use_lower_tweak_vec,
  lower_tweak_mean = lower_tweak_mean_vec,
  lower_tweak_sd = lower_tweak_sd_vec,
  q_pred = q_pred_vec
)

# ------------------------------------------------------------
# Compile and sample
# ------------------------------------------------------------
mod <- cmdstan_model(stan_file)

fit <- mod$sample(
  data = stan_data,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1500,
  iter_sampling = 1500,
  adapt_delta = 0.98,
  max_treedepth = 13,
  refresh = 200
)

# ------------------------------------------------------------
# Diagnostics
# ------------------------------------------------------------
cat("\nPosterior summaries:\n")
print(
  fit$summary(
    variables = c(
      "lower_est", "upper_est", "sigma_frac"
    )
  )
)

cat("\nDiagnostic summary:\n")
print(fit$diagnostic_summary())

# ------------------------------------------------------------
# Extract shared Q and parameter curves
# ------------------------------------------------------------
q_shape_summ <- fit$summary(variables = "Q_pred") %>%
  mutate(
    param_id = as.integer(str_match(variable, "Q_pred\\[([0-9]+),([0-9]+)\\]")[, 2]),
    grid_id  = as.integer(str_match(variable, "Q_pred\\[([0-9]+),([0-9]+)\\]")[, 3])
  ) %>%
  left_join(param_meta %>% select(param_id, parameter), by = "param_id") %>%
  mutate(
    tau = tau_pred[grid_id],
    relative_day = tau * max_day,
    panel_title = factor(panel_lookup[parameter], levels = panel_order)
  ) %>%
  select(parameter, param_id, grid_id, tau, relative_day, mean, median, sd, q5, q95, panel_title)

curve_summ <- fit$summary(variables = "theta_pred") %>%
  mutate(
    param_id = as.integer(str_match(variable, "theta_pred\\[([0-9]+),([0-9]+)\\]")[, 2]),
    grid_id  = as.integer(str_match(variable, "theta_pred\\[([0-9]+),([0-9]+)\\]")[, 3])
  ) %>%
  left_join(param_meta %>% select(param_id, parameter), by = "param_id") %>%
  mutate(
    tau = tau_pred[grid_id],
    relative_day = tau * max_day,
    panel_title = factor(panel_lookup[parameter], levels = panel_order)
  ) %>%
  select(parameter, param_id, grid_id, tau, relative_day, mean, median, sd, q5, q95, panel_title)

# ------------------------------------------------------------
# Important final orientation/magnitude fix for community unsafe funerals
# ------------------------------------------------------------
# Q_conflict is deliberately relative:
#
#   Q_conflict(t) = successful_response(t) / max(successful_response)
#
# This is the right scale for the shared response shape used by hospitalisation,
# ETU and IPC. However, p_unsafe_funeral_comm is intended to be on the
# ABSOLUTE Warsame empirical scale:
#
#   p_unsafe_funeral_comm(t) = 1 - successful_response(t)
#
# Therefore it must NOT be read directly from theta_pred fitted against Q,
# and it must also NOT be reconstructed as 1 - Q_conflict. Q_conflict is
# scaled to 0-1, but the Warsame successful-response trajectory is on an
# absolute probability scale and its maximum is less than 1.
#
# The deterministic override below therefore uses the directly interpolated
# absolute Warsame proxy:
#
#   unsafe_funeral_comm_pred_vec = 1 - successful_response(t)
#
# This is the same quantity shown by the grey dashed/circle overlay in the
# p_unsafe_funeral_comm panel.
j_ufc_abs <- param_meta$param_id[param_meta$parameter == "p_unsafe_funeral_comm"]
if (length(j_ufc_abs) != 1) stop("Could not uniquely identify p_unsafe_funeral_comm for absolute override.")

ufc_abs_curve_summ <- tibble(
  parameter = "p_unsafe_funeral_comm",
  param_id = j_ufc_abs,
  grid_id = seq_along(q_pred_vec),
  tau = tau_pred,
  relative_day = tau_pred * max_day,
  mean = unsafe_funeral_comm_pred_vec,
  median = mean,
  sd = 0,
  q5 = mean,
  q95 = mean,
  panel_title = factor("Unsafe funeral after community death", levels = panel_order)
)

curve_summ <- curve_summ %>%
  filter(parameter != "p_unsafe_funeral_comm") %>%
  bind_rows(ufc_abs_curve_summ) %>%
  arrange(param_id, grid_id)

lower_summ <- fit$summary(variables = "lower_est") %>%
  mutate(
    param_id = as.integer(str_match(variable, "lower_est\\[([0-9]+)\\]")[, 2]),
    bound = "lower"
  )

upper_summ <- fit$summary(variables = "upper_est") %>%
  mutate(
    param_id = as.integer(str_match(variable, "upper_est\\[([0-9]+)\\]")[, 2]),
    bound = "upper"
  )

bound_summ <- bind_rows(lower_summ, upper_summ) %>%
  left_join(param_meta %>% select(param_id, parameter), by = "param_id") %>%
  select(parameter, bound, mean, median, sd, q5, q95, rhat, ess_bulk, ess_tail)

# ------------------------------------------------------------
# Plots
# ------------------------------------------------------------
plot_points_workbook <- fit_df %>%
  filter(source_type == "workbook_anchor")

plot_points_removed <- removed_reference_points

# Empirical shared conflict reference on the p_unsafe_funeral_comm scale.
#
# Important: Q_conflict is relative and max-scaled. The unsafe-funeral panel
# should therefore show the absolute empirical proxy 1 - successful_response,
# not 1 - Q_conflict. This prevents the curve hitting zero merely because
# Q has been scaled to 1 at the best achieved response.
plot_points_conflict_ufc_proxy <- q_source %>%
  arrange(relative_day) %>%
  mutate(
    rolling_ref_row = row_number(),
    keep_sparse_ref = case_when(
      is_drc_plusplus_collapse_point ~ TRUE,
      !plot_sparse_rolling_ufc_dots ~ TRUE,
      plot_always_include_start_end_ufc_dots & rolling_ref_row %in% c(1L, n()) ~ TRUE,
      TRUE ~ ((rolling_ref_row - 1L) %% plot_rolling_reference_every_n) == 0L
    )
  ) %>%
  filter(keep_sparse_ref) %>%
  transmute(
    relative_day = relative_day,
    y_obs = unsafe_funeral_comm_proxy_abs,
    parameter = "p_unsafe_funeral_comm",
    panel_title = factor("Unsafe funeral after community death", levels = panel_order),
    is_drc_plusplus_collapse_point = is_drc_plusplus_collapse_point
  )

# Empirical shared Q reference retained for optional checks, but not plotted
# in the panelised Q diagnostic to keep the final figure uncluttered.
plot_points_conflict_q_shape <- q_source %>%
  transmute(
    relative_day = relative_day,
    q_shape = q_conflict_shape
  ) %>%
  tidyr::crossing(
    panel_title = factor(panel_order, levels = panel_order)
  )

p_main <- ggplot(curve_summ, aes(x = relative_day, y = mean)) +
  geom_ribbon(aes(ymin = q5, ymax = q95), fill = "#c7dcec", alpha = 0.9) +
  geom_line(linewidth = 0.9, colour = "#1f77b4") +
  geom_point(
    data = plot_points_workbook,
    aes(x = relative_day, y = y_obs),
    inherit.aes = FALSE,
    colour = "#ff7f0e",
    size = 2.2
  ) +
  geom_text(
    data = plot_points_workbook,
    aes(x = relative_day, y = y_obs, label = anchor_id),
    inherit.aes = FALSE,
    hjust = 0,
    nudge_x = 2,
    size = 2.5,
    colour = "grey25"
  ) +
  geom_point(
    data = plot_points_conflict_ufc_proxy,
    aes(x = relative_day, y = y_obs),
    inherit.aes = FALSE,
    colour = "grey45",
    size = 2.0,
    alpha = 0.9
  ) +
  geom_point(
    data = plot_points_conflict_ufc_proxy %>% filter(is_drc_plusplus_collapse_point),
    aes(x = relative_day, y = y_obs),
    inherit.aes = FALSE,
    colour = "red3",
    size = 2.4,
    alpha = 0.95
  ) +
  geom_point(
    data = plot_points_removed,
    aes(x = relative_day, y = y_obs),
    inherit.aes = FALSE,
    shape = 4,
    colour = "grey20",
    size = 2.3,
    stroke = 0.9
  ) +
  facet_wrap(~ panel_title, scales = "free_y", ncol = 2) +
  theme_bw(base_size = 11) +
  labs(
    title = "DRC with conflict ++: original-methodology direct-collapse sensitivity",
    subtitle = paste0("Original shared empirical Q framework retained; Warsame successful response forced to zero over the ++ plateau; p_unsafe_funeral_comm is absolute 1 - successful_response; red points mark forced collapse; aggregation=", aggregation_unit, "; smoothing=", q_smoothing_method),
    x = "Relative outbreak day",
    y = NULL
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.minor = element_blank()
  )

print(p_main)

ggsave(
  filename = paste0(out_prefix, "_curves.png"),
  plot = p_main,
  width = 12,
  height = 10,
  dpi = 180
)

p_q <- ggplot(q_shape_summ, aes(x = relative_day, y = mean)) +
  geom_ribbon(aes(ymin = q5, ymax = q95), fill = "#dbe8f4", alpha = 0.9) +
  geom_line(linewidth = 0.9) +
  scale_y_continuous(limits = c(0, 1)) +
  facet_wrap(~ panel_title, scales = "fixed", ncol = 2) +
  theme_bw(base_size = 11) +
  labs(
    title = "DRC with conflict ++: shared empirical Q shape with temporary collapse plateau",
    subtitle = paste0("Q is successful-response / max(successful-response), with successful response forced to zero over the DRC++ plateau; aggregation=", aggregation_unit, "; smoothing=", q_smoothing_method),
    x = "Relative outbreak day",
    y = "Shared Q_conflict(t)"
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.minor = element_blank()
  )

print(p_q)

ggsave(
  filename = paste0(out_prefix, "_Qplot.png"),
  plot = p_q,
  width = 12,
  height = 10,
  dpi = 180
)

# ------------------------------------------------------------
# BP-ready matrix
# ------------------------------------------------------------
q_value_df <- tibble(
  relative_day = tau_pred * max_day,
  tau = tau_pred,
  q_value = q_pred_vec
)

bp_matrix <- curve_summ %>%
  select(parameter, relative_day, tau, mean) %>%
  pivot_wider(names_from = parameter, values_from = mean) %>%
  mutate(
    prob_hosp = p_hosp,
    delay_hosp = delay_hosp,
    prob_unsafe_funeral_comm = p_unsafe_funeral_comm,
    prob_unsafe_funeral_hosp = p_unsafe_funeral_hosp,
    prob_unsafe_funeral_etu = 0,
    prop_etu = p_ETU,
    ipc_helper = latent_IPC
  ) %>%
  select(
    relative_day,
    tau,
    prob_hosp,
    delay_hosp,
    prob_unsafe_funeral_comm,
    prob_unsafe_funeral_hosp,
    prob_unsafe_funeral_etu,
    prop_etu,
    ipc_helper
  ) %>%
  left_join(q_value_df, by = c("relative_day", "tau"))

# ------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------
write_csv(curve_summ, paste0(out_prefix, "_curve_summaries.csv"))
write_csv(q_shape_summ, paste0(out_prefix, "_Q_summaries.csv"))
write_csv(bound_summ, paste0(out_prefix, "_bound_summaries.csv"))
write_csv(bp_matrix, paste0(out_prefix, "_bp_input_matrix.csv"))
write_csv(plot_points_workbook, paste0(out_prefix, "_anchor_rows_used.csv"))
write_csv(param_meta, paste0(out_prefix, "_prior_setup.csv"))

settings_summary <- tibble(
  setting = c(
    "provinces_to_average",
    "average_method",
    "aggregation_unit",
    "sdb_origin_filter",
    "initial_small_n_action",
    "min_sdb_eligible_for_q",
    "q_smoothing_method",
    "rolling_window_n",
    "rolling_align",
    "rolling_use_eligible_weights",
    "plot_sparse_rolling_ufc_dots",
    "plot_rolling_reference_every_n",
    "q_loess_span",
    "q_spline_df",
    "q_scale_max_abs_success",
    "q_final_max_before_defensive_rescale",
    "q_scaling",
    "plusplus_window_day",
    "plusplus_plateau_start_day",
    "plusplus_plateau_end_day",
    "plusplus_success_value",
    "ufc_comm_min_absolute_proxy",
    "ufc_comm_proxy_scale",
    "remove_superseded_ufc_anchor_ids"
  ),
  value = c(
    paste(provinces_to_average, collapse = ", "),
    average_method,
    aggregation_unit,
    sdb_origin_filter,
    initial_small_n_action,
    as.character(min_sdb_eligible_for_q),
    q_smoothing_method,
    as.character(rolling_window_n),
    rolling_align,
    as.character(rolling_use_eligible_weights),
    as.character(plot_sparse_rolling_ufc_dots),
    as.character(plot_rolling_reference_every_n),
    as.character(q_loess_span),
    as.character(q_spline_df),
    as.character(q_scale_max),
    as.character(q_final_max),
    "success_avg_smoothed / max(success_avg_smoothed) using finite-date Warsame support points only; no subtraction of empirical minimum and no artificial prediction-grid forcing",
    paste(plusplus_window_day, collapse = " to "),
    as.character(plusplus_plateau_start_day),
    as.character(plusplus_plateau_end_day),
    as.character(plusplus_success_value),
    as.character(clip01(1 - q_scale_max)),
    "p_unsafe_funeral_comm curve is directly interpolated from absolute 1 - success_avg_smoothed, not reconstructed from Q_conflict",
    paste(remove_superseded_ufc_anchor_ids, collapse = ", ")
  )
)

write_csv(settings_summary, paste0(out_prefix, "_settings_summary.csv"))

cat("\nDone.\n")
cat("Files written:\n")
cat(" - ", stan_file, "\n", sep = "")
cat(" - ", paste0(out_prefix, "_warsame_orange_line_qc.png"), "\n", sep = "")
cat(" - ", paste0(out_prefix, "_shared_conflict_Q_qc.png"), "\n", sep = "")
cat(" - ", paste0(out_prefix, "_curve_summaries.csv"), "\n", sep = "")
cat(" - ", paste0(out_prefix, "_Q_summaries.csv"), "\n", sep = "")
cat(" - ", paste0(out_prefix, "_bound_summaries.csv"), "\n", sep = "")
cat(" - ", paste0(out_prefix, "_bp_input_matrix.csv"), "\n", sep = "")
cat(" - ", paste0(out_prefix, "_anchor_rows_used.csv"), "\n", sep = "")
cat(" - ", paste0(out_prefix, "_prior_setup.csv"), "\n", sep = "")
cat(" - ", paste0(out_prefix, "_province_binned_sdb_success.csv"), "\n", sep = "")
cat(" - ", paste0(out_prefix, "_averaged_conflict_binned_sdb_success.csv"), "\n", sep = "")
cat(" - ", paste0(out_prefix, "_shared_conflict_Q_points.csv"), "\n", sep = "")
cat(" - ", paste0(out_prefix, "_shared_conflict_Q_points_before_plusplus.csv"), "\n", sep = "")
cat(" - ", paste0(out_prefix, "_plusplus_points_forced_to_collapse.csv"), "\n", sep = "")
cat(" - ", paste0(out_prefix, "_settings_summary.csv"), "\n", sep = "")
