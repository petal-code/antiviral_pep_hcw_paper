# ============================================================================
# ODV NHP delay efficacy -- PARTIAL-POOLING (hierarchical) Stan fit
# ============================================================================
#
# Why this script
# ---------------
# An alternative to coarsening the baseline (script 04) for the same problem:
# the fine 12-interval baseline in script 02 has many independent flat-prior
# nuisance hazards, and marginalising them biases efficacy upward and makes it
# over-precise (the incidental-parameters / Neyman-Scott problem). Here we KEEP
# the fine 12-interval baseline but place a hierarchical prior on the log
# hazards,
#     log_lambda_k ~ Normal(mu_log, tau_log),
# so the sparse intervals are shrunk together (partial pooling). This reduces the
# effective number of free nuisance parameters and hence the bias, while keeping
# the flexible baseline shape.
#
# Uses stan-models/odv_delay_efficacy_hier.stan.
# Output: data-processed/odv_nhp_delay/odv_ebov_rhesus_delay_efficacy_fit_stan_partialpool.rds
#
# Run from anywhere in the repo:
#   Rscript analyses/odv_nhp_delay_efficacy/05_fit_odv_delay_efficacy_stan_partialpool.R
# ============================================================================

suppressPackageStartupMessages({
  library(cmdstanr)
})
source(here::here("analyses", "odv_nhp_delay_efficacy", "odv_delay_helper_functions.R"))

# ---- Paths ----
raw_path  <- here::here("data-raw", "odv_nhp_delay", "odv_ebov_rhesus_delay_survival.csv")
out_dir   <- here::here("data-processed", "odv_nhp_delay")
out_path  <- file.path(out_dir, "odv_ebov_rhesus_delay_efficacy_fit_stan_partialpool.rds")
stan_file <- here::here("analyses", "odv_nhp_delay_efficacy",
                        "stan-models", "odv_delay_efficacy_hier.stan")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Settings (shared values match script 02; fine baseline kept) ----
settings <- list(
  t_censor_plot = 28, t_fit_end = 15, dpc_zero = 15, eps_hr = 1e-6,
  B_emp = 500, seed = 1,
  fit_k = FALSE, k_fixed = 1, k_prior_logmean = log(1), k_prior_logsd = 0.5,
  d50_prior_mean = 4, d50_prior_sd = 2,

  # --- Hierarchical hyperpriors on the log baseline hazards ---
  # mu_log: overall log-hazard level (data-informed; vague prior centred on a
  #         daily hazard of ~0.1). tau_log: between-interval sd, half-normal --
  #         smaller tau => more pooling => more bias reduction. Tunable.
  mu_prior_mean = log(0.1), mu_prior_sd = 2, tau_prior_sd = 1,

  chains = 4, iter_warmup = 1500, iter_sampling = 1500,
  adapt_delta = 0.99, max_treedepth = 12   # higher adapt_delta for the hierarchy
)
set.seed(settings$seed)

# ---- Data prep (shared helpers; FINE breakpoints as in script 02) ----
raw_dat        <- read_odv_raw(raw_path)
individual_dat <- expand_individual(raw_dat)
observed_dpc   <- sort(unique(individual_dat$dpc[individual_dat$dpc > 0]))
empirical_points <- compute_empirical_points(individual_dat, observed_dpc, settings)

cuts <- make_fine_cuts(individual_dat, observed_dpc, settings)
sp   <- build_split_data(individual_dat, settings, cuts)
split_dat <- sp$split_dat
K <- sp$K
message("Fine baseline with partial pooling: K = ", K, " intervals")

# ---- Stan data (hierarchical model: common fields + hyperpriors) ----
curve_grid <- seq(0, settings$dpc_zero, by = 0.02)
stan_data <- c(
  build_common_stan_data(split_dat, observed_dpc, settings, K, curve_grid),
  list(mu_prior_mean = settings$mu_prior_mean,
       mu_prior_sd   = settings$mu_prior_sd,
       tau_prior_sd  = settings$tau_prior_sd)
)

# ---- Fit ----
mod <- cmdstan_model(stan_file)
fit <- mod$sample(
  data = stan_data, seed = settings$seed,
  chains = settings$chains, parallel_chains = settings$chains,
  iter_warmup = settings$iter_warmup, iter_sampling = settings$iter_sampling,
  adapt_delta = settings$adapt_delta, max_treedepth = settings$max_treedepth,
  refresh = 200
)
cat("\nSampler diagnostics:\n"); diag_summ <- fit$diagnostic_summary(); print(diag_summ)
map_estimate <- run_map(mod, stan_data, settings$seed)
cat("\nMAP (mode):\n"); print(map_estimate)

# ---- Summaries + save (also summarise the pooling hyperparameters) ----
ss <- summarise_stan_fit(fit, settings, individual_dat, empirical_points,
                         curve_grid, diag_summ,
                         summary_vars = c("E0", "d50", "k", "mu_log", "tau_log"))
cat("\nPooling hyperparameters (mu_log = overall level, tau_log = between-interval sd):\n")
print(ss$param_summary[ss$param_summary$variable %in% c("mu_log", "tau_log"), ])

output <- list(
  metadata = list(
    analysis = "ODV NHP delay-efficacy (Stan, hierarchical/partial-pooled baseline)",
    method = paste(
      "Piecewise-exponential survival model fit in Stan with a FINE", K,
      "-interval baseline whose log-hazards are partially pooled",
      "(log_lambda_k ~ Normal(mu_log, tau_log)), to reduce the incidental-parameter",
      "(Neyman-Scott) inflation of efficacy seen with independent flat-prior hazards."),
    raw_path = raw_path, output_path = out_path, stan_model = stan_file,
    created_by = "analyses/odv_nhp_delay_efficacy/05_fit_odv_delay_efficacy_stan_partialpool.R",
    settings = settings, K = K
  ),
  raw_data = raw_dat, individual_survival_data = individual_dat,
  split_survival_data = split_dat, empirical_points = empirical_points,
  stan_data = stan_data, fitted_curve = ss$fitted_curve,
  param_summary = ss$param_summary, fit_summary = ss$fit_summary,
  diagnostics = diag_summ, map_estimate = map_estimate,
  draws = as.data.frame(fit$draws(c("E0", "d50", "k", "mu_log", "tau_log"), format = "df"))
)
saveRDS(output, out_path)
message("Saved partial-pooling fit to: ", out_path)
print(ss$fit_summary)
