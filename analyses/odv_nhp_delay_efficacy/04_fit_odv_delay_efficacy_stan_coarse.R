# ============================================================================
# ODV NHP delay efficacy -- COARSENED-BASELINE Stan fit
# ============================================================================
#
# Why this script
# ---------------
# 02_fit_odv_delay_efficacy_stan.R uses a fine piecewise-constant baseline
# hazard (one interval per observed death time => K = 12 for only 16 deaths).
# Marginalising that many independent, flat-prior baseline hazards biases the
# fitted efficacy upward and makes it over-precise (the incidental-parameters /
# Neyman-Scott problem; see the analysis notes). This script reduces the number
# of nuisance hazards by COARSENING the baseline to ~6 intervals, which moves the
# posterior back onto the profiled/MLE answer with honest (wider) intervals.
#
# Note: this reuses the SAME Stan model as 02 (stan-models/odv_delay_efficacy.stan).
# The baseline granularity is a property of the split data fed in, not the model,
# so only the breakpoints change here. The baseline hazards keep their flat prior.
#
# Output: data-processed/odv_nhp_delay/odv_ebov_rhesus_delay_efficacy_fit_stan_coarse.rds
#
# Run from anywhere in the repo:
#   Rscript analyses/odv_nhp_delay_efficacy/04_fit_odv_delay_efficacy_stan_coarse.R
# ============================================================================

suppressPackageStartupMessages({
  library(cmdstanr)
})
source(here::here("analyses", "odv_nhp_delay_efficacy", "odv_delay_helper_functions.R"))

# ---- Paths ----
raw_path  <- here::here("data-raw", "odv_nhp_delay", "odv_ebov_rhesus_delay_survival.csv")
out_dir   <- here::here("data-processed", "odv_nhp_delay")
out_path  <- file.path(out_dir, "odv_ebov_rhesus_delay_efficacy_fit_stan_coarse.rds")
stan_file <- here::here("analyses", "odv_nhp_delay_efficacy",
                        "stan-models", "odv_delay_efficacy.stan")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Settings (shared values match script 02) ----
settings <- list(
  t_censor_plot = 28, t_fit_end = 15, dpc_zero = 15, eps_hr = 1e-6,
  B_emp = 500, seed = 1,
  fit_k = FALSE, k_fixed = 1, k_prior_logmean = log(1), k_prior_logsd = 0.5,
  use_hazard_prior = 0L, hazard_prior_rate = 1e-3,   # flat baseline hazards (as in 02)
  d50_prior_mean = 4, d50_prior_sd = 2,

  # --- COARSE baseline: keep the dpc breakpoints {1,2,3,4} so treatment timing
  # stays exact, plus a single death-region cut at day 9. This gives 6 intervals:
  #   [0,1] [1,2] [2,3] [3,4] [4,9] [9,15]   (~8 deaths either side of day 9)
  coarse_death_cuts = c(9),

  chains = 4, iter_warmup = 1000, iter_sampling = 1000,
  adapt_delta = 0.95, max_treedepth = 12
)
set.seed(settings$seed)

# ---- Data prep (shared helpers; coarse breakpoints) ----
raw_dat        <- read_odv_raw(raw_path)
individual_dat <- expand_individual(raw_dat)
observed_dpc   <- sort(unique(individual_dat$dpc[individual_dat$dpc > 0]))
empirical_points <- compute_empirical_points(individual_dat, observed_dpc, settings)

cuts <- make_coarse_cuts(observed_dpc, settings, settings$coarse_death_cuts)
sp   <- build_split_data(individual_dat, settings, cuts)
split_dat <- sp$split_dat
K <- sp$K
message("Coarsened baseline: K = ", K, " intervals")

# ---- Stan data (base model: common fields + flat-hazard-prior switch) ----
curve_grid <- seq(0, settings$dpc_zero, by = 0.02)
stan_data <- c(
  build_common_stan_data(split_dat, observed_dpc, settings, K, curve_grid),
  list(use_hazard_prior  = as.integer(settings$use_hazard_prior),
       hazard_prior_rate = settings$hazard_prior_rate)
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

# ---- Summaries + save ----
ss <- summarise_stan_fit(fit, settings, individual_dat, empirical_points,
                         curve_grid, diag_summ)

output <- list(
  metadata = list(
    analysis = "ODV NHP delay-efficacy (Stan, coarsened 6-interval baseline)",
    method = paste(
      "Piecewise-exponential survival model fit in Stan with a COARSE",
      "(", K, "-interval) baseline hazard, to reduce the incidental-parameter",
      "(Neyman-Scott) inflation of efficacy seen with the fine 12-interval baseline."),
    raw_path = raw_path, output_path = out_path, stan_model = stan_file,
    created_by = "analyses/odv_nhp_delay_efficacy/04_fit_odv_delay_efficacy_stan_coarse.R",
    settings = settings, K = K
  ),
  raw_data = raw_dat, individual_survival_data = individual_dat,
  split_survival_data = split_dat, empirical_points = empirical_points,
  stan_data = stan_data, fitted_curve = ss$fitted_curve,
  param_summary = ss$param_summary, fit_summary = ss$fit_summary,
  diagnostics = diag_summ, map_estimate = map_estimate,
  draws = as.data.frame(fit$draws(c("E0", "d50", "k"), format = "df"))
)
saveRDS(output, out_path)
message("Saved coarse-baseline fit to: ", out_path)
print(ss$fit_summary)
