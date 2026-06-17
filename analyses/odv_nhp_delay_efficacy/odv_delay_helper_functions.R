# ============================================================================
# odv_delay_helper_functions.R
# ----------------------------------------------------------------------------
# Shared data-prep and post-fit helpers for the ODV NHP delay-efficacy Stan
# fits (scripts 04 / 05). The deterministic data processing mirrors
# 01_fit_odv_delay_efficacy.R / 02_fit_odv_delay_efficacy_stan.R exactly; it is
# factored out here so the baseline-coarsening (04) and partial-pooling (05)
# variants do not duplicate ~200 lines each.
# ============================================================================

suppressPackageStartupMessages({
  library(survival)
})

parse_death_times <- function(x) {
  x <- trimws(as.character(x))
  if (is.na(x) || x == "") return(numeric(0))
  as.numeric(strsplit(x, ";", fixed = TRUE)[[1]])
}

clamp01 <- function(x) pmin(1, pmax(0, x))

km_surv_at <- function(dat, t_eval) {
  sf <- survfit(Surv(time, event) ~ 1, data = dat)
  as.numeric(summary(sf, times = t_eval, extend = TRUE)$surv[1])
}

# Read + validate the raw arm-level CSV and check internal consistency.
read_odv_raw <- function(raw_path) {
  if (!file.exists(raw_path)) stop("Could not find raw data file: ", raw_path)
  raw_dat <- read.csv(raw_path, stringsAsFactors = FALSE)
  required_cols <- c("arm", "dpc", "n_total", "n_survivors_reported",
                     "death_times", "censor_day")
  miss <- setdiff(required_cols, names(raw_dat))
  if (length(miss) > 0)
    stop("Raw data file is missing required columns: ", paste(miss, collapse = ", "))
  raw_dat$n_deaths_extracted <- vapply(
    raw_dat$death_times, function(x) length(parse_death_times(x)), integer(1))
  raw_dat$n_survivors_extracted <- raw_dat$n_total - raw_dat$n_deaths_extracted
  if (any(raw_dat$n_survivors_extracted != raw_dat$n_survivors_reported))
    stop("Extracted death times do not match reported survivor counts in raw CSV.")
  raw_dat
}

# Expand arm-level rows to one row per animal.
expand_individual <- function(raw_dat) {
  make_arm_df <- function(row) {
    dts <- parse_death_times(row$death_times)
    n_total <- as.integer(row$n_total)
    n_cens <- n_total - length(dts)
    if (n_cens < 0) stop("More death times than animals in arm: ", row$arm)
    data.frame(
      subject_id = sprintf("%s_%02d", row$arm, seq_len(n_total)),
      arm = row$arm, dpc = as.integer(row$dpc),
      time = c(dts, rep(as.numeric(row$censor_day), n_cens)),
      event = c(rep(1L, length(dts)), rep(0L, n_cens)),
      stringsAsFactors = FALSE)
  }
  ind <- do.call(rbind, lapply(seq_len(nrow(raw_dat)),
                               function(i) make_arm_df(raw_dat[i, ])))
  ind$time <- as.numeric(ind$time)
  ind$event <- as.integer(ind$event)
  ind
}

# Empirical hazard-scale efficacy points (1 - KM hazard ratio) + bootstrap CIs.
compute_empirical_points <- function(individual_dat, observed_dpc, settings) {
  veh <- individual_dat[individual_dat$arm == "vehicle", ]
  odv <- individual_dat[individual_dat$arm != "vehicle", ]
  H_veh <- -log(km_surv_at(veh, settings$t_censor_plot))

  ep <- do.call(rbind, lapply(observed_dpc, function(d) {
    dd <- odv[odv$dpc == d, ]
    S_d <- km_surv_at(dd, settings$t_censor_plot)
    data.frame(dpc = d, n_total = nrow(dd), n_survivors = sum(dd$event == 0L),
               survival = S_d, efficacy_hazard_scale = 1 - (-log(S_d)) / H_veh,
               stringsAsFactors = FALSE)
  }))

  boot <- matrix(NA_real_, nrow = settings$B_emp, ncol = length(observed_dpc))
  for (b in seq_len(settings$B_emp)) {
    vb <- veh[sample(seq_len(nrow(veh)), replace = TRUE), ]
    Hb <- -log(km_surv_at(vb, settings$t_censor_plot))
    if (!is.finite(Hb) || Hb <= 0) next
    for (j in seq_along(observed_dpc)) {
      dd <- odv[odv$dpc == observed_dpc[j], ]
      db <- dd[sample(seq_len(nrow(dd)), replace = TRUE), ]
      Hd <- -log(km_surv_at(db, settings$t_censor_plot))
      if (is.finite(Hd)) boot[b, j] <- 1 - Hd / Hb
    }
  }
  ep$efficacy_lo <- NA_real_; ep$efficacy_hi <- NA_real_
  for (j in seq_along(observed_dpc)) {
    vals <- boot[, j]; vals <- vals[is.finite(vals)]
    if (length(vals) >= 50) {
      ep$efficacy_lo[ep$dpc == observed_dpc[j]] <- quantile(vals, 0.025)
      ep$efficacy_hi[ep$dpc == observed_dpc[j]] <- quantile(vals, 0.975)
    }
  }
  ep$efficacy_hazard_scale <- clamp01(ep$efficacy_hazard_scale)
  ep$efficacy_lo <- clamp01(ep$efficacy_lo)
  ep$efficacy_hi <- clamp01(ep$efficacy_hi)
  ep
}

# Interior interval breakpoints (excluding 0 and t_fit_end).
# FINE: a breakpoint at every observed death time (matches scripts 01/02).
make_fine_cuts <- function(individual_dat, observed_dpc, settings) {
  ad <- sort(unique(individual_dat$time[individual_dat$event == 1L]))
  ad <- ad[ad <= settings$t_fit_end]
  br <- sort(unique(c(0, ad, observed_dpc, settings$t_fit_end)))
  br[-c(1, length(br))]
}

# COARSE: keep the dpc breakpoints (so treatment timing stays exact) plus a
# small user-supplied set of death-region cuts.
make_coarse_cuts <- function(observed_dpc, settings, death_cuts) {
  br <- sort(unique(c(0, observed_dpc, death_cuts, settings$t_fit_end)))
  br[-c(1, length(br))]
}

# Split animal-level data into (animal x interval) rows given interior cuts.
build_split_data <- function(individual_dat, settings, cuts) {
  fd <- individual_dat
  fd$time_fit  <- pmin(fd$time, settings$t_fit_end)
  fd$event_fit <- as.integer(fd$event == 1L & fd$time <= settings$t_fit_end)
  sd <- survSplit(Surv(time_fit, event_fit) ~ ., data = fd,
                  cut = cuts, episode = "interval")
  sd$tstart   <- as.numeric(sd$tstart)
  sd$tstop    <- as.numeric(sd$time_fit)
  sd$exposure <- sd$tstop - sd$tstart
  sd <- sd[sd$exposure > 0, ]
  sd$interval     <- factor(sd$interval)
  sd$interval_idx <- as.integer(sd$interval)
  sd$post_treatment <- as.integer(sd$dpc > 0 & sd$tstart >= sd$dpc)
  list(split_dat = sd, K = length(levels(sd$interval)))
}

# Stan data fields common to BOTH the base and hierarchical models. Each driver
# appends the fields specific to its model (hazard-prior switch vs hyperpriors).
build_common_stan_data <- function(split_dat, observed_dpc, settings, K, curve_grid) {
  dpc_idx <- match(split_dat$dpc, observed_dpc)
  dpc_idx[is.na(dpc_idx)] <- 1L
  list(
    N = nrow(split_dat), K = K, D = length(observed_dpc),
    interval_idx   = as.integer(split_dat$interval_idx),
    event          = as.integer(split_dat$event_fit),
    exposure       = as.numeric(split_dat$exposure),
    post_treatment = as.integer(split_dat$post_treatment),
    dpc_idx        = as.integer(dpc_idx),
    dpc_obs        = as.numeric(observed_dpc),
    d_zero         = settings$dpc_zero,
    eps_hr         = settings$eps_hr,
    fit_k           = as.integer(settings$fit_k),
    k_fixed         = settings$k_fixed,
    k_prior_logmean = settings$k_prior_logmean,
    k_prior_logsd   = settings$k_prior_logsd,
    d50_prior_mean  = settings$d50_prior_mean,
    d50_prior_sd    = settings$d50_prior_sd,
    n_grid   = length(curve_grid),
    grid_dpc = as.numeric(curve_grid)
  )
}

# MAP (mode) check, robust to older CmdStan where jacobian= is unavailable.
run_map <- function(mod, stan_data, seed, vars = c("E0", "d50", "k")) {
  tryCatch({
    opt <- tryCatch(
      mod$optimize(data = stan_data, jacobian = FALSE, seed = seed),
      error = function(e) mod$optimize(data = stan_data, seed = seed))
    as.data.frame(opt$summary(vars))
  }, error = function(e) {
    message("optimize() failed (", conditionMessage(e), "); skipping MAP check.")
    NULL
  })
}

# Post-fit summaries: parameter table, pointwise fitted curve, compact fit row.
summarise_stan_fit <- function(fit, settings, individual_dat, empirical_points,
                               curve_grid, diag_summ,
                               summary_vars = c("E0", "d50", "k")) {
  param_summ <- fit$summary(summary_vars)

  grid_draws <- fit$draws(variables = "efficacy_grid", format = "draws_matrix")
  gq <- apply(grid_draws, 2, quantile, probs = c(0.025, 0.5, 0.975), names = FALSE)
  curve_dat <- data.frame(dpc = curve_grid, efficacy = gq[2, ],
                          efficacy_lo = gq[1, ], efficacy_hi = gq[3, ])

  eff_summ <- fit$summary("eff")
  rmse_emp <- sqrt(mean((eff_summ$median - empirical_points$efficacy_hazard_scale)^2))

  gp <- function(nm, col) param_summ[[col]][param_summ$variable == nm]
  fit_summary <- data.frame(
    E0_median = gp("E0", "median"), E0_lo = gp("E0", "q5"), E0_hi = gp("E0", "q95"),
    d50_median = gp("d50", "median"), d50_lo = gp("d50", "q5"), d50_hi = gp("d50", "q95"),
    k_median = gp("k", "median"), k_estimated = as.logical(settings$fit_k),
    dpc_zero = settings$dpc_zero, rmse_empirical_points = rmse_emp,
    n_animals = nrow(individual_dat), n_events = sum(individual_dat$event),
    n_divergent = sum(diag_summ$num_divergent),
    max_rhat = max(param_summ$rhat, na.rm = TRUE),
    min_ess_bulk = min(param_summ$ess_bulk, na.rm = TRUE),
    stringsAsFactors = FALSE)

  list(param_summary = as.data.frame(param_summ),
       fitted_curve = curve_dat, fit_summary = fit_summary)
}
