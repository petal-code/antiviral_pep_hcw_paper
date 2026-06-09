# ============================================================
# DRC non-conflict middle scenario: Warsame-only shared-Q specification,
# original/full DRC time horizon, original no-conflict plateau
#
# This is the sensitivity Jacob intended as "script 1":
#   - Keep the same time horizon as the original DRC non-conflict scenario.
#   - Keep the same no-conflict bridge/plateau used in the original DRC
#     non-conflict script: each province improves to its data-rich monthly
#     maximum and then plateaus; the two province-specific plateau curves are
#     averaged to make the black no-conflict Warsame line.
#   - Change only the source of the shared Q shape: instead of estimating Q
#     from all workbook parameters, define Q entirely from that bridged
#     no-conflict Warsame black line.
#
# In other words, this script should NOT use the raw conflict-interrupted
# weekly Warsame trajectory. It uses the original monthly no-conflict plateau
# trajectory, but lets that trajectory provide the shared Q shape for all
# parameters.
#
#   Q_no_conflict(t) = successful_SDB_black_line(t) / max(successful_SDB_black_line), with Q(0) = 0
#   p_unsafe_funeral_comm(t) = 1 - successful_SDB_black_line(t)
#
# The UFC curve is on the absolute Warsame scale and is NOT 1 - Q.
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

stan_file  <- "drc_non_conflict_warsame_only_Q_minus_conflict_time_original_plateau.stan"
out_prefix <- "drc_non_conflict_warsame_only_Q_minus_conflict_time_original_plateau"

# ------------------------------------------------------------
# Settings for reconstructing the Warsame orange line
# ------------------------------------------------------------
provinces_to_average <- c("North Kivu", "Ituri")

# Orange line denominator.
# To match the original DRC non-conflict script and reproduce the same no-conflict
# black-line plateau, use the original published-curve reconstruction:
#   success / (success + failure)
# i.e. do not add "sdb not needed" to the numerator for this specific
# backwards-compatible sensitivity. If you want the corrected successful-response
# definition instead, set this to "successful_response_vs_failure", but note that
# the plateau will then differ from the original non-conflict plot.
success_denominator_mode <- "success_vs_failure"

# Primary choice remains the unweighted province mean, i.e. average the
# province-specific black lines rather than weighting by alert volume.
average_method <- "unweighted_province_mean"  # or "count_weighted"

# Temporal aggregation.
# For this original-timeline non-conflict sensitivity, use the same monthly
# aggregation as the original no-conflict script. This avoids carrying the raw
# week-to-week conflict-interrupted Warsame wiggles into the non-conflict Q.
aggregation_unit <- "monthly"  # "epi_week", "monthly", or "four_week"
bin_width_days <- 28            # only used when aggregation_unit == "four_week"

# Primary version uses all SDB responses, consistent with the published line.
sdb_origin_filter <- "all"     # "all" or "community"

# Early success spikes can still be tiny-n artefacts.
initial_small_n_action <- "zero_success"  # "zero_success", "drop", or "none"
initial_spike_max_day <- 75
initial_spike_min_eligible <- 10
initial_spike_success_threshold <- 0.50

# No-conflict construction: for each province, once the reconstructed
# successful-SDB no-conflict trajectory reaches its province-specific maximum, later
# deterioration is interpreted as conflict-related and replaced by a plateau.
use_no_conflict_plateau <- TRUE
plateau_max_min_eligible <- 25
plateau_after_first_max <- TRUE

# Keep the original/full DRC time horizon in this version. The paired
# "minus conflict time" scripts set this to TRUE and truncate at the first
# data-rich maximum of the bridged Warsame black line.
remove_conflict_time_from_non_conflict <- TRUE
non_conflict_end_min_eligible <- 25
non_conflict_end_tolerance <- 1e-8

# Minimum eligible SDB count required for a bin to contribute to the
# shared Q_conflict shape. With weekly reconstruction, keep this low so the
# plotted line explicitly follows the Warsame curve.
min_sdb_eligible_for_q <- 1

# Do not apply an additional rolling average in this script: the intended shape
# is the original monthly no-conflict plateau black line itself.
q_smoothing_method <- "none"   # "none", "rolling4", "loess", or "spline"
rolling_window_n <- 4
rolling_align <- "center"          # "center" or "trailing"
rolling_use_eligible_weights <- FALSE  # FALSE = simple rolling mean of reconstructed points

# Plotting only: show a sparse set of grey rolling-average reference dots
# in the p_unsafe_funeral_comm panel. The fitted curve still uses the full
# weekly rolling-average trajectory; this just avoids visually overcrowding
# the panel with one dot for every epidemiological week. For weekly data and
# rolling_window_n = 4, this shows approximately one reference dot per 4 weeks.
plot_sparse_rolling_ufc_dots <- TRUE
plot_rolling_reference_every_n <- 1
plot_always_include_start_end_ufc_dots <- TRUE

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



apply_no_conflict_plateau <- function(
    df,
    y_col,
    date_col,
    eligible_col,
    min_eligible_for_max = 25,
    plateau_after_first_max = TRUE
) {
  out <- df[order(df[[date_col]]), , drop = FALSE]

  y <- as.numeric(out[[y_col]])
  eligible <- as.numeric(out[[eligible_col]])

  valid_idx <- which(!is.na(y) & !is.na(eligible) & eligible >= min_eligible_for_max)
  if (length(valid_idx) == 0) {
    out[[paste0(y_col, "_no_conflict")]] <- y
    out$plateau_applied <- FALSE
    out$plateau_onset_date <- as.Date(NA)
    out$plateau_max_success <- NA_real_
    return(out)
  }

  max_success <- max(y[valid_idx], na.rm = TRUE)
  if (!is.finite(max_success)) {
    out[[paste0(y_col, "_no_conflict")]] <- y
    out$plateau_applied <- FALSE
    out$plateau_onset_date <- as.Date(NA)
    out$plateau_max_success <- NA_real_
    return(out)
  }

  max_idx <- valid_idx[which(y[valid_idx] >= max_success - 1e-8)]
  plateau_start_idx <- if (plateau_after_first_max) min(max_idx) else max(max_idx)

  y_no_conflict <- y
  y_no_conflict[seq.int(plateau_start_idx, length(y_no_conflict))] <- max_success

  out[[paste0(y_col, "_no_conflict")]] <- pmin(1, pmax(0, y_no_conflict))
  out$plateau_applied <- seq_len(nrow(out)) >= plateau_start_idx
  out$plateau_onset_date <- out[[date_col]][plateau_start_idx]
  out$plateau_max_success <- max_success
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
# Build monthly Warsame successful-SDB values by province, apply the
# no-conflict plateau rule province-by-province, then average provinces.
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

prov_binned_unbridged <- sdb_prov %>%
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
    # Keep both definitions available, but default to the original non-conflict
    # reconstruction so the time horizon and plateau match the earlier script.
    n_success_metric = case_when(
      success_denominator_mode == "success_vs_failure" ~ as.numeric(n_success),
      success_denominator_mode == "successful_response_vs_failure" ~ as.numeric(n_successful_response),
      success_denominator_mode == "successful_response_all_nonmissing" ~ as.numeric(n_successful_response),
      TRUE ~ NA_real_
    ),
    n_eligible = case_when(
      success_denominator_mode == "success_vs_failure" ~ n_success + n_failure,
      success_denominator_mode == "successful_response_vs_failure" ~ n_successful_response + n_failure,
      # Retained only as a sensitivity option. Because n_successful_response
      # already includes "sdb not needed", do not add n_not_needed again here.
      success_denominator_mode == "successful_response_all_nonmissing" ~ n_successful_response + n_failure + n_false_alert,
      TRUE ~ NA_real_
    ),
    prop_success_raw = if_else(n_eligible > 0, n_success_metric / n_eligible, NA_real_),
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

prov_binned_list <- split(prov_binned_unbridged, prov_binned_unbridged$province)

prov_binned <- bind_rows(lapply(prov_binned_list, function(df_prov) {
  this_province <- unique(df_prov$province)
  if (length(this_province) != 1) stop("Province split failed.")

  if (use_no_conflict_plateau) {
    out <- apply_no_conflict_plateau(
      df = df_prov,
      y_col = "prop_success_clean",
      date_col = "bin_mid_date",
      eligible_col = "n_eligible",
      min_eligible_for_max = plateau_max_min_eligible,
      plateau_after_first_max = plateau_after_first_max
    ) %>%
      mutate(
        prop_success_no_conflict = prop_success_clean_no_conflict,
        prop_success_bridged = prop_success_no_conflict,
        bridged_by_conflict_rule = plateau_applied
      )
  } else {
    out <- df_prov %>%
      mutate(
        prop_success_no_conflict = prop_success_clean,
        prop_success_bridged = prop_success_no_conflict,
        plateau_applied = FALSE,
        plateau_onset_date = as.Date(NA),
        plateau_max_success = NA_real_,
        bridged_by_conflict_rule = FALSE
      )
  }

  out
})) %>%
  arrange(bin_mid_date, province)

sdb_avg <- prov_binned %>%
  group_by(bin_label, bin_start, bin_mid_date) %>%
  summarise(
    relative_day = as.numeric(first(bin_mid_date) - sdb_start_date),
    n_provinces_with_data = n_distinct(province),
    n_eligible_sum = sum(n_eligible, na.rm = TRUE),
    prop_success_unweighted = mean(prop_success_no_conflict, na.rm = TRUE),
    prop_success_count_weighted = if_else(
      sum(n_eligible, na.rm = TRUE) > 0,
      weighted.mean(prop_success_no_conflict, w = pmax(n_eligible, 0), na.rm = TRUE),
      NA_real_
    ),
    any_initial_spike_adjusted = any(small_n_initial_spike),
    any_plateau_applied = any(plateau_applied),
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

# Keep bins sufficiently data-rich to define the no-conflict Q shape.
q_source <- sdb_avg %>%
  filter(n_eligible_sum >= min_sdb_eligible_for_q) %>%
  filter(is.finite(relative_day), !is.na(bin_mid_date), is.finite(prop_success_avg)) %>%
  arrange(bin_mid_date)

if (nrow(q_source) < 3) {
  stop("Too few SDB bins remain after min_sdb_eligible_for_q filter.")
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
# Remove conflict-time tail from the no-conflict scenario
# ------------------------------------------------------------
non_conflict_end_day <- max(q_source$relative_day, na.rm = TRUE)
non_conflict_end_date <- max(q_source$bin_mid_date, na.rm = TRUE)

if (remove_conflict_time_from_non_conflict) {
  end_candidates <- q_source %>%
    filter(
      is.finite(relative_day),
      is.finite(success_avg_smoothed),
      is.finite(n_eligible_sum),
      n_eligible_sum >= non_conflict_end_min_eligible
    )

  if (nrow(end_candidates) < 1) {
    stop("Cannot remove conflict time: no data-rich no-conflict Warsame bins available to define endpoint.")
  }

  end_success <- max(end_candidates$success_avg_smoothed, na.rm = TRUE)
  non_conflict_end_row <- end_candidates %>%
    filter(success_avg_smoothed >= end_success - non_conflict_end_tolerance) %>%
    arrange(relative_day) %>%
    slice(1)

  non_conflict_end_day <- non_conflict_end_row$relative_day[[1]]
  non_conflict_end_date <- non_conflict_end_row$bin_mid_date[[1]]

  if (!is.finite(non_conflict_end_day) || non_conflict_end_day <= 0) {
    stop("Invalid no-conflict endpoint after conflict-time removal.")
  }

  q_source <- q_source %>%
    filter(relative_day <= non_conflict_end_day)

  sdb_avg <- sdb_avg %>%
    filter(relative_day <= non_conflict_end_day)

  prov_binned <- prov_binned %>%
    filter(relative_day <= non_conflict_end_day)

  message(
    "Removed conflict time from no-conflict Warsame-only DRC scenario: endpoint = day ",
    round(non_conflict_end_day, 1),
    " (", non_conflict_end_date, ")."
  )
}

# ------------------------------------------------------------
# Convert the empirical black average SUCCESS line into the shared Q_conflict shape.
# Do not invert this with 1 - success here. Q is the response-improvement signal.
#
# Important: Q_conflict is a RELATIVE 0-1 response index.
#
# We therefore scale the empirical smoothed successful-SDB no-conflict trajectory by its
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
max_day <- if (remove_conflict_time_from_non_conflict) {
  non_conflict_end_day
} else {
  max(
    max(anchors$relative_day, na.rm = TRUE),
    max(q_source$relative_day, na.rm = TRUE)
  )
}

if (!is.finite(max_day) || max_day <= 0) {
  stop("Could not construct a valid global DRC time horizon.")
}

# Add the forced outbreak-start anchor on the ABSOLUTE successful-response scale.
# Then compute Q from that same empirical successful-SDB no-conflict black line. This keeps
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
    success_avg_smoothed = 0
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
    unsafe_funeral_comm_proxy_abs = clip01(1 - success_avg_smoothed)
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

if (remove_conflict_time_from_non_conflict) {
  anchors_for_fit <- anchors_for_fit %>%
    filter(relative_day <= max_day)

  removed_reference_points <- removed_reference_points %>%
    filter(relative_day <= max_day)
}

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
    title = "DRC Warsame SDB orange-line reconstruction: monthly no-conflict successful-SDB curve",
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
    title = "DRC no-conflict Warsame-only shared Q shape",
    subtitle = paste0("Empirical averaged successful-SDB no-conflict black line, from the original monthly no-conflict black line, empirically max-scaled using the observed maximum; aggregation=", aggregation_unit, "; smoothing=", q_smoothing_method),
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
  filename = paste0(out_prefix, "_shared_no_conflict_Q_qc.png"),
  plot = p_q_conflict_qc,
  width = 10,
  height = 6,
  dpi = 180
)

write_csv(prov_binned, paste0(out_prefix, "_province_binned_sdb_success.csv"))
write_csv(sdb_avg, paste0(out_prefix, "_averaged_no_conflict_binned_sdb_success.csv"))
write_csv(q_source, paste0(out_prefix, "_shared_no_conflict_Q_points.csv"))

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
    dplyr::select(parameter, lb_prior_mean, ub_prior_mean, abs_min, abs_max,
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
  left_join(param_meta %>% dplyr::select(parameter, param_id, increases), by = "parameter") %>%
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
  left_join(param_meta %>% dplyr::select(parameter, param_id, increases), by = "parameter") %>%
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
  left_join(param_meta %>% dplyr::select(param_id, parameter), by = "param_id") %>%
  mutate(
    tau = tau_pred[grid_id],
    relative_day = tau * max_day,
    panel_title = factor(panel_lookup[parameter], levels = panel_order)
  ) %>%
  dplyr::select(parameter, param_id, grid_id, tau, relative_day, mean, median, sd, q5, q95, panel_title)

curve_summ <- fit$summary(variables = "theta_pred") %>%
  mutate(
    param_id = as.integer(str_match(variable, "theta_pred\\[([0-9]+),([0-9]+)\\]")[, 2]),
    grid_id  = as.integer(str_match(variable, "theta_pred\\[([0-9]+),([0-9]+)\\]")[, 3])
  ) %>%
  left_join(param_meta %>% dplyr::select(param_id, parameter), by = "param_id") %>%
  mutate(
    tau = tau_pred[grid_id],
    relative_day = tau * max_day,
    panel_title = factor(panel_lookup[parameter], levels = panel_order)
  ) %>%
  dplyr::select(parameter, param_id, grid_id, tau, relative_day, mean, median, sd, q5, q95, panel_title)

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
# scaled to 0-1, but the Warsame successful-SDB no-conflict trajectory is on an
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
  left_join(param_meta %>% dplyr::select(param_id, parameter), by = "param_id") %>%
  dplyr::select(parameter, bound, mean, median, sd, q5, q95, rhat, ess_bulk, ess_tail)

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
    panel_title = factor("Unsafe funeral after community death", levels = panel_order)
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
    title = "DRC no-conflict Warsame-only scenario: original monthly plateau Q, minus conflict time, minus conflict time",
    subtitle = paste0("Q is relative for response parameters; p_unsafe_funeral_comm is absolute 1 - successful_response; grey circles show sparse monthly no-conflict absolute UFC proxy; aggregation=", aggregation_unit, "; smoothing=", q_smoothing_method, "; DRC_UFC_00 and DRC_UFC_01 reference-only"),
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
    title = "DRC no-conflict scenario: shared Warsame-only Q shape applied to all response parameters",
    subtitle = paste0("Q is the empirical averaged successful-SDB no-conflict trajectory, from the original monthly no-conflict black line, empirically max-scaled after adding Q(0)=0; aggregation=", aggregation_unit, "; smoothing=", q_smoothing_method),
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
bp_matrix <- curve_summ %>%
  dplyr::select(parameter, relative_day, tau, mean) %>%
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
  dplyr::select(
    relative_day,
    tau,
    prob_hosp,
    delay_hosp,
    prob_unsafe_funeral_comm,
    prob_unsafe_funeral_hosp,
    prob_unsafe_funeral_etu,
    prop_etu,
    ipc_helper
  )

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
    "use_no_conflict_plateau",
    "remove_conflict_time_from_non_conflict",
    "non_conflict_end_day",
    "non_conflict_end_date",
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
    as.character(use_no_conflict_plateau),
    as.character(remove_conflict_time_from_non_conflict),
    as.character(non_conflict_end_day),
    as.character(non_conflict_end_date),
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
cat(" - ", paste0(out_prefix, "_shared_no_conflict_Q_qc.png"), "\n", sep = "")
cat(" - ", paste0(out_prefix, "_curve_summaries.csv"), "\n", sep = "")
cat(" - ", paste0(out_prefix, "_Q_summaries.csv"), "\n", sep = "")
cat(" - ", paste0(out_prefix, "_bound_summaries.csv"), "\n", sep = "")
cat(" - ", paste0(out_prefix, "_bp_input_matrix.csv"), "\n", sep = "")
cat(" - ", paste0(out_prefix, "_anchor_rows_used.csv"), "\n", sep = "")
cat(" - ", paste0(out_prefix, "_prior_setup.csv"), "\n", sep = "")
cat(" - ", paste0(out_prefix, "_province_binned_sdb_success.csv"), "\n", sep = "")
cat(" - ", paste0(out_prefix, "_averaged_no_conflict_binned_sdb_success.csv"), "\n", sep = "")
cat(" - ", paste0(out_prefix, "_shared_no_conflict_Q_points.csv"), "\n", sep = "")
cat(" - ", paste0(out_prefix, "_settings_summary.csv"), "\n", sep = "")
