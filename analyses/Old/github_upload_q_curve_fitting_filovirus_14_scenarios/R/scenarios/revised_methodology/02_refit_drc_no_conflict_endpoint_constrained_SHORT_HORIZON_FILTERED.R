# ============================================================
# DRC no-conflict scenario
# Endpoint-constrained partial-pooling refit / matrix rebuild
#
# Key change versus earlier scripts:
#   * Parameter start/end magnitudes are locked to early/late
#     literature-supported extrema where available; curated
#     scenario ranges are used only as fallbacks.
#   * The partial-pooling model estimates only the response-shape
#     parameters (t50 and slope k), not lower/upper magnitudes.
#   * Orange points are plotted with vertical uncertainty bars.
#     If row-level low/high values are not provided, this script
#     uses a conservative plotting/likelihood default based on a
#     fraction of the endpoint span.
#
# Inputs expected in working directory:
#   - filovirus_model_parameter_table_4_scenarios_with_etu_baseline(2).xlsx
#
# Outputs:
#   - worst_west_africa_endpoint_constrained_matrix.csv
#   - worst_west_africa_endpoint_constrained_curve_summaries.csv
#   - worst_west_africa_endpoint_constrained_anchor_table_used.csv
#   - worst_west_africa_endpoint_constrained_endpoint_constrained_plot.png/pdf
# ============================================================

suppressPackageStartupMessages({
  library(cmdstanr)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(posterior)
})

# -----------------------
# User settings
# -----------------------
parameter_workbook <- "filovirus_model_parameter_table_4_scenarios_with_etu_baseline.xlsx"
source_sheet <- "DRC no conflict"
out_prefix <- "drc_no_conflict_endpoint_constrained_zero_plateau"
stan_file <- paste0(out_prefix, "_endpoint_constrained_partial_pool.stan")

scenario_name <- "Middle_DRC_NoConflict"
# Shortened no-conflict horizon: remove conflict time from the historical DRC timeline.
# Anchors beyond this horizon are excluded from fitting and plotting because they
# refer to the conflict-affected historical period, not the counterfactual
# shortened no-conflict scenario.
scenario_duration_days <- 193
plot_x_max_days <- 200
fit_anchor_max_day <- scenario_duration_days
n_pred <- 100

# Optional: if you manually curate row-specific literature uncertainty,
# save a CSV with columns: anchor_id,value_low,value_high.
manual_anchor_bounds_csv <- "manual_anchor_bounds.csv"

# Default observation uncertainty when row-specific ranges are not available.
# This affects fitting, not the fixed endpoints.
default_obs_sd_frac_of_span <- 0.08
min_obs_sd_abs <- 1e-4

# -----------------------
# Mappings
# -----------------------
param_map <- c(
  "prob_hospitalised_genPop; prob_hospitalised_hcw" = "prob_hosp",
  "hospitalisation_delay_factor" = "delay_hosp",
  "p_unsafe_funeral_comm_genPop; p_unsafe_funeral_comm_hcw" = "prob_unsafe_funeral_comm",
  "p_unsafe_funeral_hosp_genPop; p_unsafe_funeral_hosp_hcw" = "prob_unsafe_funeral_hosp",
  "prob_unsafe_funeral_etu" = "prob_unsafe_funeral_etu",
  "prop_etu" = "prop_etu",
  "ipc_helper; ppe_efficacy_hcw" = "ipc_helper"
)

raw_param_map <- c(
  "delay_hosp" = "delay_hosp",
  "p_hosp" = "prob_hosp",
  "p_ETU" = "prop_etu",
  "latent_IPC" = "ipc_helper",
  "p_unsafe_funeral_comm" = "prob_unsafe_funeral_comm",
  "p_unsafe_funeral_hosp" = "prob_unsafe_funeral_hosp"
)

direction_map <- c(
  "prob_hosp" = "increasing",
  "delay_hosp" = "decreasing",
  "prob_unsafe_funeral_comm" = "decreasing",
  "prob_unsafe_funeral_hosp" = "decreasing",
  "prob_unsafe_funeral_etu" = "flat",
  "prop_etu" = "increasing",
  "ipc_helper" = "increasing"
)

panel_lookup <- c(
  "delay_hosp" = "Delay to hospitalisation",
  "prob_hosp" = "Probability of hospitalisation",
  "prop_etu" = "Proportion in ETU / ETC",
  "ipc_helper" = "Latent IPC / PPE index",
  "prob_unsafe_funeral_comm" = "Unsafe funeral after community death",
  "prob_unsafe_funeral_hosp" = "Unsafe funeral after hospital death",
  "prob_unsafe_funeral_etu" = "Unsafe funeral after ETU death"
)

matrix_cols <- c(
  "prob_hosp", "delay_hosp", "prob_unsafe_funeral_comm",
  "prob_unsafe_funeral_hosp", "prob_unsafe_funeral_etu",
  "prop_etu", "ipc_helper"
)

# -----------------------
# Helpers
# -----------------------
clip01 <- function(x) pmin(1, pmax(0, x))

extract_range <- function(x) {
  x <- as.character(x)
  x <- stringr::str_replace_all(x, "[−—]", "–")
  if (stringr::str_detect(x, "–")) {
    parts <- stringr::str_split_fixed(x, "–", 2)
  } else if (stringr::str_detect(x, "\\d\\s*-\\s*\\d")) {
    parts <- stringr::str_split_fixed(x, "\\s*-\\s*", 2)
  } else {
    nums <- stringr::str_extract_all(x, "[-+]?\\d*\\.?\\d+(?:[eE][-+]?\\d+)?")[[1]]
    if (length(nums) == 0) return(c(NA_real_, NA_real_))
    if (length(nums) == 1) return(rep(as.numeric(nums[1]), 2))
    return(as.numeric(nums[1:2]))
  }
  out <- rep(NA_real_, 2)
  for (i in 1:2) {
    nums <- stringr::str_extract_all(parts[i], "[-+]?\\d*\\.?\\d+(?:[eE][-+]?\\d+)?")[[1]]
    if (length(nums) > 0) out[i] <- as.numeric(nums[1])
  }
  out
}

parse_anchor_value <- function(x) {
  x <- as.character(x)
  day <- as.numeric(stringr::str_match(x, "relative\\s+day\\s+([0-9.]+)\\s*:")[,2])
  val <- as.numeric(stringr::str_match(x, ":\\s*([0-9.]+)")[,2])
  tibble(relative_day = day, value_used = val)
}

parse_from_description <- function(x, key) {
  pat <- paste0(key, ":\\s*([^;]+)")
  out <- stringr::str_match(as.character(x), pat)[,2]
  stringr::str_squish(out)
}

make_endpoint_values <- function(parameter, lo, hi) {
  direction <- direction_map[[parameter]]
  if (is.na(lo) || is.na(hi)) stop("Missing endpoint for ", parameter)
  if (direction == "increasing") return(c(start = lo, end = hi))
  if (direction == "decreasing") return(c(start = hi, end = lo))
  return(c(start = lo, end = lo))
}

q_raw <- function(tau, t50, k) plogis(k * (tau - t50))
q_norm <- function(tau, t50, k) {
  q0 <- q_raw(0, t50, k)
  q1 <- q_raw(1, t50, k)
  (q_raw(tau, t50, k) - q0) / (q1 - q0)
}

# -----------------------
# Read parameter table
# -----------------------
raw <- readxl::read_excel(parameter_workbook, sheet = source_sheet, col_names = FALSE)
names(raw)[1:7] <- c("Section","Parameter","Symbol","Description","Value_range","Reference","URL")

summary_ranges <- raw %>%
  filter(.data$Section == "Time-varying response") %>%
  filter(.data$Parameter %in% names(param_map)) %>%
  mutate(
    parameter = unname(param_map[.data$Parameter]),
    range_pair = lapply(.data$Value_range, extract_range),
    endpoint_low = vapply(range_pair, function(z) z[1], numeric(1)),
    endpoint_high = vapply(range_pair, function(z) z[2], numeric(1))
  ) %>%
  select(parameter, endpoint_low, endpoint_high, Value_range, Reference, URL)

anchors <- raw %>%
  filter(.data$Section == "Time-varying curve anchors") %>%
  filter(!is.na(.data$Parameter), .data$Parameter != "Orange-point fitted anchors") %>%
  mutate(
    anchor_id = as.character(.data$Parameter),
    raw_parameter = parse_from_description(.data$Description, "Parameter"),
    parameter = unname(raw_param_map[raw_parameter]),
    fit_role = parse_from_description(.data$Description, "fit role"),
    evidence_tier = parse_from_description(.data$Description, "evidence tier")
  ) %>%
  bind_cols(parse_anchor_value(.$Value_range)) %>%
  filter(.data$parameter %in% matrix_cols) %>%
  filter(is.finite(.data$relative_day), is.finite(.data$value_used)) %>%
  mutate(
    value_low = value_used,
    value_high = value_used
  )

if (file.exists(manual_anchor_bounds_csv)) {
  manual_bounds <- readr::read_csv(manual_anchor_bounds_csv, show_col_types = FALSE) %>%
    select(anchor_id, value_low, value_high)
  anchors <- anchors %>%
    select(-value_low, -value_high) %>%
    left_join(manual_bounds, by = "anchor_id") %>%
    mutate(
      value_low = if_else(is.na(value_low), value_used, as.numeric(value_low)),
      value_high = if_else(is.na(value_high), value_used, as.numeric(value_high))
    )
}

anchors <- anchors %>%
  left_join(
    summary_ranges %>% rename(endpoint_value_range = Value_range, endpoint_reference = Reference, endpoint_url = URL),
    by = "parameter"
  ) %>%
  mutate(
    tau = clip01(relative_day / scenario_duration_days),
    direction = unname(direction_map[parameter])
  )

if (nrow(anchors) == 0) stop("No usable anchors found for sheet: ", source_sheet)

# Critical shortened-timeline fix:
# Drop all literature anchors beyond the counterfactual no-conflict horizon.
# These points belong to the historical conflict-affected tail and should not
# elongate fitted no-conflict curves or trigger terminal-zero plateaus.
anchors_pre_horizon_filter <- anchors
anchors <- anchors %>%
  filter(relative_day <= fit_anchor_max_day + 1e-8)

post_horizon_anchors <- anchors_pre_horizon_filter %>%
  anti_join(anchors %>% select(anchor_id), by = "anchor_id")

readr::write_csv(
  post_horizon_anchors %>%
    select(any_of(c(
      "anchor_id", "parameter", "relative_day", "value_used",
      "fit_role", "evidence_tier", "endpoint_reference", "endpoint_url"
    ))),
  paste0(out_prefix, "_post_horizon_anchors_excluded.csv")
)

if (nrow(anchors) == 0) {
  stop(
    "All anchors were beyond the shortened DRC no-conflict horizon. ",
    "Check scenario_duration_days / fit_anchor_max_day."
  )
}

# Build endpoint table.
# IMPORTANT: endpoints are locked using pragmatic endpoint-window extrema,
# not earliest/latest anchors and not global scenario summary ranges.
#
# Rationale:
#   - increasing parameters should start at an early observed minimum and end at
#     a late observed maximum;
#   - decreasing parameters should start at an early observed maximum and end at
#     a late observed minimum;
#   - if a terminal zero is observed for unsafe-funeral parameters, the endpoint
#     is forced to zero even if the zero anchor occurs before the final model day;
#   - if no credible early/late anchor exists, fall back to the curated scenario
#     range for now. Those rows are flagged in the endpoint table so they can be
#     revisited as model-predicted endpoints if needed.
#
# This should retain the sensible partial-pooled shapes from the earlier plots,
# while making the curve begin/end at the intended literature-supported values.
early_endpoint_window_day <- 50
# For the shortened DRC no-conflict scenario, define the late endpoint window
# relative to the shortened scenario end, not the original conflict timeline.
late_endpoint_window_width <- 50
late_endpoint_start_day   <- max(0, scenario_duration_days - late_endpoint_window_width)
zero_endpoint_parameters <- c(
  "prob_unsafe_funeral_comm",
  "prob_unsafe_funeral_hosp",
  "prob_unsafe_funeral_etu"
)
zero_tolerance <- 1e-12

eligible_endpoint_anchors <- anchors %>%
  filter(!is.na(parameter), is.finite(value_used), is.finite(relative_day)) %>%
  filter(is.na(fit_role) | !fit_role %in% c("reference_only", "plot_only"))

choose_endpoint_row <- function(df, direction, window, endpoint_type) {
  if (nrow(window) == 0) return(NULL)

  # endpoint_type is "start" or "end". The value selection depends on direction.
  if (direction == "increasing") {
    idx <- if (endpoint_type == "start") which.min(window$value_used) else which.max(window$value_used)
  } else if (direction == "decreasing") {
    idx <- if (endpoint_type == "start") which.max(window$value_used) else which.min(window$value_used)
  } else {
    # Flat parameter: use the earliest row in the window.
    idx <- which.min(window$relative_day)
  }
  window[idx[1], , drop = FALSE]
}

endpoint_anchor_values <- eligible_endpoint_anchors %>%
  group_by(parameter) %>%
  group_modify(function(.x, .y) {
    p <- .y$parameter[[1]]
    direction <- unname(direction_map[p])

    early_window <- .x %>% filter(relative_day <= early_endpoint_window_day)
    late_window  <- .x %>% filter(relative_day >= late_endpoint_start_day)

    start_row <- choose_endpoint_row(.x, direction, early_window, "start")

    # Unsafe-funeral parameters are special: if an explicit zero appears anywhere,
    # force the final endpoint to zero. This prevents the curve ending above zero
    # simply because the zero anchor is not inside the late endpoint window.
    zero_rows <- if (p %in% zero_endpoint_parameters) {
      .x %>% filter(abs(value_used) <= zero_tolerance)
    } else {
      .x[0, , drop = FALSE]
    }
    if (nrow(zero_rows) > 0) {
      end_row <- zero_rows %>% arrange(desc(relative_day)) %>% slice(1)
      end_source <- "terminal_zero_anchor_anywhere"
    } else {
      end_row <- choose_endpoint_row(.x, direction, late_window, "end")
      end_source <- if (!is.null(end_row)) "late_window_extreme_anchor" else NA_character_
    }

    tibble(
      start_anchor_day = if (!is.null(start_row)) start_row$relative_day[[1]] else NA_real_,
      end_anchor_day   = if (!is.null(end_row))   end_row$relative_day[[1]]   else NA_real_,
      start_anchor_value = if (!is.null(start_row)) start_row$value_used[[1]] else NA_real_,
      end_anchor_value   = if (!is.null(end_row))   end_row$value_used[[1]]   else NA_real_,
      start_anchor_low   = if (!is.null(start_row)) start_row$value_low[[1]] else NA_real_,
      start_anchor_high  = if (!is.null(start_row)) start_row$value_high[[1]] else NA_real_,
      end_anchor_low     = if (!is.null(end_row))   end_row$value_low[[1]]   else NA_real_,
      end_anchor_high    = if (!is.null(end_row))   end_row$value_high[[1]]  else NA_real_,
      start_endpoint_source = if (!is.null(start_row)) "early_window_extreme_anchor" else NA_character_,
      end_endpoint_source   = end_source
    )
  }) %>%
  ungroup()

endpoint_table <- summary_ranges %>%
  filter(parameter %in% matrix_cols) %>%
  left_join(endpoint_anchor_values, by = "parameter") %>%
  mutate(
    fallback_start_value = mapply(function(p, lo, hi) make_endpoint_values(p, lo, hi)[["start"]],
                                  parameter, endpoint_low, endpoint_high),
    fallback_end_value = mapply(function(p, lo, hi) make_endpoint_values(p, lo, hi)[["end"]],
                                parameter, endpoint_low, endpoint_high),
    start_value = if_else(is.finite(start_anchor_value), start_anchor_value, fallback_start_value),
    end_value   = if_else(is.finite(end_anchor_value),   end_anchor_value,   fallback_end_value),
    start_endpoint_source = if_else(is.na(start_endpoint_source), "summary_range_fallback_no_early_anchor", start_endpoint_source),
    end_endpoint_source   = if_else(is.na(end_endpoint_source),   "summary_range_fallback_no_late_anchor",  end_endpoint_source),
    endpoint_source = paste(start_endpoint_source, end_endpoint_source, sep = " / "),
    direction = unname(direction_map[parameter]),
    # For ordinary endpoints, Q reaches its endpoint at the scenario end.
    # For terminal-zero unsafe-funeral anchors, Q reaches the zero anchor at
    # that anchor day and the parameter remains at zero thereafter.
    endpoint_day_for_tau = if_else(
      end_endpoint_source == "terminal_zero_anchor_anywhere" & is.finite(end_anchor_day),
      pmax(end_anchor_day, 1),
      as.numeric(scenario_duration_days)
    )
  )

# Audit endpoint choices explicitly; this is the first file to check if a panel
# does not start/end where expected.
readr::write_csv(endpoint_table, paste0(out_prefix, "_endpoint_table_used.csv"))

# Recompute anchor tau using the endpoint timing actually used for each parameter.
# This matters for terminal-zero unsafe-funeral parameters: if the literature
# anchor reaches zero before the overall scenario end, the fitted shape must
# reach Q=1 at that zero-anchor day, then stay flat at zero afterwards.
anchors <- anchors %>%
  select(-any_of("endpoint_day_for_tau")) %>%
  left_join(endpoint_table %>% select(parameter, endpoint_day_for_tau), by = "parameter") %>%
  mutate(
    endpoint_day_for_tau = if_else(is.finite(endpoint_day_for_tau), endpoint_day_for_tau, as.numeric(scenario_duration_days)),
    tau = clip01(relative_day / endpoint_day_for_tau)
  )

# Add flat ETU funeral if not present in anchors but present in summary.
anchors_fit <- anchors %>%
  filter(parameter != "prob_unsafe_funeral_etu")

param_levels <- endpoint_table %>%
  filter(parameter != "prob_unsafe_funeral_etu") %>%
  pull(parameter)

J <- length(param_levels)
param_id_lookup <- setNames(seq_along(param_levels), param_levels)

stan_data <- anchors_fit %>%
  mutate(
    param_id = unname(param_id_lookup[parameter])
  ) %>%
  filter(!is.na(param_id))

endpoint_fit <- endpoint_table %>%
  filter(parameter %in% param_levels) %>%
  arrange(match(parameter, param_levels))

span <- abs(endpoint_fit$end_value - endpoint_fit$start_value)
obs_sd <- pmax(
  abs(stan_data$value_high - stan_data$value_low) / 4,
  default_obs_sd_frac_of_span * span[stan_data$param_id],
  min_obs_sd_abs
)

write_csv(
  anchors %>% select(any_of(c("anchor_id", "parameter", "relative_day", "tau", "value_used", "value_low", "value_high",
                            "fit_role", "evidence_tier", "Reference", "URL", "endpoint_reference", "endpoint_url"))),
  paste0(out_prefix, "_anchor_table_used.csv")
)

# -----------------------
# Stan model: endpoints fixed, Q shape partially pooled
# -----------------------
stan_code <- '
functions {
  real q_raw(real tau, real t50, real k) {
    return inv_logit(k * (tau - t50));
  }
  real q_norm(real tau, real t50, real k) {
    real q0 = q_raw(0, t50, k);
    real q1 = q_raw(1, t50, k);
    return (q_raw(tau, t50, k) - q0) / (q1 - q0);
  }
}
data {
  int<lower=1> N;
  int<lower=1> J;
  array[N] int<lower=1, upper=J> param_id;
  vector<lower=0, upper=1>[N] tau;
  vector[N] y_obs;
  vector<lower=1e-9>[N] obs_sd;
  vector[J] theta_start;
  vector[J] theta_end;
  real mu_t50_raw_prior_mean;
  real<lower=0> mu_t50_raw_prior_sd;
  real mu_log_k_prior_mean;
  real<lower=0> mu_log_k_prior_sd;
  int<lower=1> N_pred;
  vector<lower=0, upper=1>[N_pred] tau_pred;
}
parameters {
  real mu_t50_raw;
  real<lower=0> sigma_t50;
  vector[J] t50_std;

  real mu_log_k;
  real<lower=0> sigma_log_k;
  vector[J] log_k_std;
}
transformed parameters {
  vector<lower=0, upper=1>[J] t50_param;
  vector<lower=0>[J] k_param;
  for (j in 1:J) {
    t50_param[j] = inv_logit(mu_t50_raw + sigma_t50 * t50_std[j]);
    k_param[j] = exp(mu_log_k + sigma_log_k * log_k_std[j]);
  }
}
model {
  mu_t50_raw ~ normal(mu_t50_raw_prior_mean, mu_t50_raw_prior_sd);
  sigma_t50 ~ normal(0, 0.50);
  t50_std ~ normal(0, 1);

  mu_log_k ~ normal(mu_log_k_prior_mean, mu_log_k_prior_sd);
  sigma_log_k ~ normal(0, 0.40);
  log_k_std ~ normal(0, 1);

  for (n in 1:N) {
    int j = param_id[n];
    real q = q_norm(tau[n], t50_param[j], k_param[j]);
    real mu = theta_start[j] + (theta_end[j] - theta_start[j]) * q;
    y_obs[n] ~ student_t(4, mu, obs_sd[n]);
  }
}
generated quantities {
  matrix[J, N_pred] Q_pred;
  matrix[J, N_pred] theta_pred;
  for (j in 1:J) {
    for (m in 1:N_pred) {
      real q = q_norm(tau_pred[m], t50_param[j], k_param[j]);
      Q_pred[j, m] = q;
      theta_pred[j, m] = theta_start[j] + (theta_end[j] - theta_start[j]) * q;
    }
  }
}
'
writeLines(stan_code, stan_file)

data_list <- list(
  N = nrow(stan_data),
  J = J,
  param_id = stan_data$param_id,
  tau = stan_data$tau,
  y_obs = stan_data$value_used,
  obs_sd = obs_sd,
  theta_start = endpoint_fit$start_value,
  theta_end = endpoint_fit$end_value,
  mu_t50_raw_prior_mean = 0,
  mu_t50_raw_prior_sd = 1.25,
  mu_log_k_prior_mean = log(8),
  mu_log_k_prior_sd = 0.80,
  N_pred = n_pred,
  tau_pred = seq(0, 1, length.out = n_pred)
)

mod <- cmdstanr::cmdstan_model(stan_file)
fit <- mod$sample(
  data = data_list,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  seed = 123
)

# -----------------------
# Summarise predictions
# -----------------------
q_summary <- fit$summary(variables = "Q_pred")
theta_summary <- fit$summary(variables = "theta_pred")

parse_idx <- function(x) {
  m <- stringr::str_match(x, "\\[(\\d+),(\\d+)\\]")
  tibble(param_id = as.integer(m[,2]), grid_id = as.integer(m[,3]))
}

curve_summ <- bind_cols(
  parse_idx(theta_summary$variable),
  theta_summary %>% select(mean, median, q5, q95)
) %>%
  mutate(
    parameter = param_levels[param_id],
    tau = data_list$tau_pred[grid_id],
    endpoint_day_for_tau = endpoint_fit$endpoint_day_for_tau[param_id],
    relative_day = tau * endpoint_day_for_tau,
    panel_title = unname(panel_lookup[parameter])
  )

q_summ <- bind_cols(
  parse_idx(q_summary$variable),
  q_summary %>% select(mean, median, q5, q95)
) %>%
  mutate(
    parameter = param_levels[param_id],
    tau = data_list$tau_pred[grid_id],
    endpoint_day_for_tau = endpoint_fit$endpoint_day_for_tau[param_id],
    relative_day = tau * endpoint_day_for_tau,
    panel_title = unname(panel_lookup[parameter])
  )

# Add fixed ETU funeral curve if present.
if ("prob_unsafe_funeral_etu" %in% endpoint_table$parameter) {
  etu_row <- endpoint_table %>% filter(parameter == "prob_unsafe_funeral_etu") %>% slice(1)
  fixed_etu <- tibble(
    param_id = max(curve_summ$param_id) + 1L,
    grid_id = seq_len(n_pred),
    mean = etu_row$start_value,
    median = etu_row$start_value,
    q5 = etu_row$start_value,
    q95 = etu_row$start_value,
    parameter = "prob_unsafe_funeral_etu",
    tau = data_list$tau_pred,
    relative_day = tau * scenario_duration_days,
    panel_title = unname(panel_lookup[parameter])
  )
  curve_summ <- bind_rows(curve_summ, fixed_etu)
}

# If any parameter reaches its endpoint before the scenario end, explicitly
# append a plateau row at the final scenario day. This makes the plot and the
# exported matrix hold the endpoint value rather than stretching the decline to
# the final model time point.
plateau_rows <- curve_summ %>%
  group_by(parameter) %>%
  filter(max(relative_day, na.rm = TRUE) < scenario_duration_days - 1e-8) %>%
  arrange(relative_day, .by_group = TRUE) %>%
  slice_tail(n = 1) %>%
  ungroup() %>%
  mutate(
    grid_id = max(curve_summ$grid_id, na.rm = TRUE) + row_number(),
    tau = 1,
    relative_day = as.numeric(scenario_duration_days)
  )

if (nrow(plateau_rows) > 0) {
  curve_summ <- bind_rows(curve_summ, plateau_rows) %>%
    arrange(parameter, relative_day)
}

readr::write_csv(curve_summ, paste0(out_prefix, "_curve_summaries.csv"))
readr::write_csv(q_summ, paste0(out_prefix, "_Q_summaries.csv"))

# Matrix output: one row per time point, one column per BP time-varying parameter.
matrix_out <- tibble(
  scenario = scenario_name,
  time_index = seq_len(n_pred),
  tau = data_list$tau_pred,
  relative_day = tau * scenario_duration_days
)

# Use mean fitted Q averaged across parameters as a diagnostic q_value on the
# common matrix grid. For parameters with early terminal endpoints, values after
# their endpoint are carried forward by interpolation below.
q_diag <- q_summ %>%
  group_by(grid_id) %>%
  summarise(q_value = mean(mean, na.rm = TRUE), .groups = "drop") %>%
  arrange(grid_id)

matrix_out$q_value <- q_diag$q_value[seq_len(n_pred)]

matrix_days <- matrix_out$relative_day
for (p in matrix_cols) {
  if (p %in% curve_summ$parameter) {
    cs <- curve_summ %>%
      filter(parameter == p) %>%
      arrange(relative_day) %>%
      group_by(relative_day) %>%
      summarise(mean = mean(mean, na.rm = TRUE), .groups = "drop")
    matrix_out[[p]] <- stats::approx(
      x = cs$relative_day,
      y = cs$mean,
      xout = matrix_days,
      method = "linear",
      rule = 2,
      ties = mean
    )$y
  } else {
    # Missing parameter: fall back to fixed zero if ETU funeral, otherwise NA.
    matrix_out[[p]] <- ifelse(p == "prob_unsafe_funeral_etu", 0, NA_real_)
  }
}

readr::write_csv(matrix_out, paste0(out_prefix, "_matrix.csv"))

# -----------------------
# Plot
# -----------------------
plot_parameters <- setdiff(matrix_cols, "prob_unsafe_funeral_etu")

plot_df <- curve_summ %>%
  filter(parameter %in% plot_parameters) %>%
  mutate(panel_title = factor(panel_title, levels = panel_lookup[plot_parameters]))

anchor_plot <- anchors %>%
  filter(parameter %in% plot_parameters) %>%
  filter(relative_day <= fit_anchor_max_day + 1e-8) %>%
  mutate(
    panel_title = unname(panel_lookup[parameter]),
    panel_title = factor(panel_title, levels = panel_lookup[plot_parameters])
  )

p <- ggplot() +
  geom_ribbon(
    data = plot_df,
    aes(x = relative_day, ymin = q5, ymax = q95),
    alpha = 0.18
  ) +
  geom_line(
    data = plot_df,
    aes(x = relative_day, y = mean),
    linewidth = 0.9
  ) +
  geom_errorbar(
    data = anchor_plot,
    aes(x = relative_day, ymin = value_low, ymax = value_high),
    width = 0,
    colour = "orange3",
    linewidth = 0.55
  ) +
  geom_point(
    data = anchor_plot,
    aes(x = relative_day, y = value_used),
    colour = "orange3",
    size = 2.1
  ) +
  geom_text(
    data = anchor_plot,
    aes(x = relative_day, y = value_used, label = anchor_id),
    colour = "orange4",
    size = 2.4,
    hjust = -0.05,
    vjust = -0.6,
    check_overlap = TRUE
  ) +
  facet_wrap(~panel_title, scales = "free_y", ncol = 2) +
  labs(
    title = paste0(scenario_name, ": endpoint-constrained time-varying parameter curves"),
    subtitle = "Endpoints locked to early/late literature extrema; terminal-zero anchors plateau at zero; shape partially pooled",
    x = "Relative day",
    y = "Parameter value"
  ) +
  coord_cartesian(xlim = c(0, plot_x_max_days), clip = "off") +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave(paste0(out_prefix, "_endpoint_constrained_plot.png"), p, width = 11, height = 8, dpi = 300)
ggsave(paste0(out_prefix, "_endpoint_constrained_plot.pdf"), p, width = 11, height = 8)

print(fit$summary(variables = c("mu_t50_raw", "sigma_t50", "mu_log_k", "sigma_log_k")))
message("Done: ", out_prefix)
