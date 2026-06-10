# ============================================================================
# 04_DRC_QCurve_ModelC_ReachEstimation.R   (TRIAL / SCAFFOLD -- untested)
# ----------------------------------------------------------------------------
# Model C: estimate a LATENT REACH curve R(t) from the DRC community-death
# proportion (binomial), keep the empirical SDB success curve s(t) SUPPLIED (as
# Model B), and re-wire the parameter -> clock map per the reach/SDB decomposition:
#
#   delay_hosp, p_hosp, p_ETU      -> reach clock  R_clock(t)         (estimated)
#   latent_IPC                     -> SDB clock (default) OR reach     (toggle)
#   p_unsafe_funeral_hosp          -> SDB clock                       (near-zero)
#   p_unsafe_funeral_comm          -> DETERMINISTIC 1 - R_prob(t)*s(t) (no endpoints)
#
# Endpoints for the 5 endpoint-parameters are estimated from the literature
# anchors with the SAME machinery as Model B (original methodology).
#
# Fits four variants for comparison: {RW, spline} x {IPC on SDB, IPC on reach}.
#
# Prereqs (run first): 00_DataPreparation_and_Cleaning.R (-> DRC_QCurve_PreppedData.rds)
#                      00b_DataPreparation_CommunityDeaths.R (-> DRC_CommunityDeaths_Prepped.rds)
#
# >>> This is a first scaffold to TRIAL the model: it has not been run here (no
#     CmdStan in the build env). Expect to debug grid sizes / priors / sampler
#     settings on first run. Key knobs are flagged below.
# ============================================================================

suppressPackageStartupMessages({
  library(cmdstanr); library(dplyr); library(tidyr); library(stringr)
  library(tibble); library(ggplot2); library(splines)
})
source(here::here("analyses", "01_latent_response_parameter_estimation", "helpers.R"))
set.seed(123)

# ---- knobs (EDIT) ----------------------------------------------------------
GRID_DT       <- 14L          # reach-grid spacing (days); G = horizon/GRID_DT + 1
SPLINE_DF     <- 14L          # B-spline basis size (spline model)
RW_SIGMA_PSD  <- 0.30         # half-normal scale on the RW step sd (per sqrt-day, logit)
SPL_SIGMA_PSD <- 0.50         # half-normal scale on the P-spline 2nd-diff sd
SAMP <- list(chains = 4, parallel_chains = 4, iter_warmup = 1500,
             iter_sampling = 1500, adapt_delta = 0.98, max_treedepth = 13, refresh = 200)
VARIANTS <- expand.grid(model = c("rw", "spline"),
                        ipc_clock = c("sdb", "reach"), stringsAsFactors = FALSE)

# ---- 1. Inputs -------------------------------------------------------------
drc_prep <- readRDS(file.path(DIR_PROCESSED, "DRC_QCurve", "DRC_QCurve_PreppedData.rds"))
anchors  <- drc_prep$anchors
qseries  <- drc_prep$conflict_qseries          # SUPPLIED SDB curve (conflict scenario)
cd_prep  <- readRDS(file.path(DIR_PROCESSED, "DRC_QCurve", "DRC_CommunityDeaths_Prepped.rds"))
cd_obs   <- cd_prep$obs

# Endpoint parameters = all six EXCEPT community unsafe funerals (now deterministic).
ENDPOINT_PARAMS <- setdiff(PARAM_LEVELS, "p_unsafe_funeral_comm")
# Which clock each endpoint parameter rides (clock_type: 1 = reach, 2 = SDB).
clock_for <- function(ipc_clock) c(
  delay_hosp            = 1L,
  p_hosp                = 1L,
  p_ETU                 = 1L,
  latent_IPC            = if (ipc_clock == "reach") 1L else 2L,
  p_unsafe_funeral_hosp = 2L
)[ENDPOINT_PARAMS]

# Hard admissible domains (mirrors 02_DRC_QCurve_Fitting_Original.R).
domain_meta <- tribble(
  ~parameter,                ~abs_min, ~abs_max,
  "delay_hosp",              0.5,      12.0,
  "p_hosp",                  0.0,      1.0,
  "p_ETU",                   0.0,      1.0,
  "latent_IPC",              0.0,      1.0,
  "p_unsafe_funeral_hosp",   0.0,      0.010
)

param_meta <- anchors %>%
  filter(parameter %in% ENDPOINT_PARAMS) %>%
  group_by(parameter) %>%
  summarise(lb_prior_mean = first(lower_bound),
            ub_prior_mean = first(upper_bound),
            direction     = first(direction), .groups = "drop") %>%
  left_join(domain_meta, by = "parameter") %>%
  mutate(param_id  = match(parameter, ENDPOINT_PARAMS),
         increases = if_else(direction == "up", 1L, 0L)) %>%
  arrange(param_id) %>%
  build_param_support()
J <- nrow(param_meta)

# ---- 2. Reach grid + supplied SDB curve on the grid ------------------------
grid_day <- seq(0L, HORIZON_DAYS, by = GRID_DT)
G        <- length(grid_day)
nearest  <- function(day) vapply(day, function(d) which.min(abs(grid_day - d)), integer(1))

# Supplied SDB ACTUAL success on the grid (success, not the normalised q_value):
#   s_ref = max success  =>  Sclk = s/s_ref = q_value (matches Model B),
#   and the burial term uses the actual success s.
s_success_grid <- clip01(make_interp(qseries$relative_day, qseries$success_smoothed)(grid_day))
s_ref          <- max(qseries$success_smoothed, na.rm = TRUE)

# ---- 3. Anchor -> grid (endpoint params only) ------------------------------
fit_anchors <- anchors %>%
  filter(parameter %in% ENDPOINT_PARAMS) %>%
  left_join(select(param_meta, parameter, param_id), by = "parameter") %>%
  mutate(weight      = if_else(is.na(weight), 1, weight),
         obs_sd_mult = 1 / pmax(weight, 0.25),
         anchor_grid = nearest(relative_day)) %>%
  arrange(param_id)

# ---- 4. Endpoint tweaks (only the near-zero hospital unsafe funeral) --------
zero <- rep(0L, J)
tw <- list(use_upper = zero, upper_mean = param_meta$ub_prior_mean, upper_sd = rep(1.0, J),
           use_lower = zero, lower_mean = param_meta$lb_prior_mean, lower_sd = rep(1.0, J))
jh <- match("p_unsafe_funeral_hosp", param_meta$parameter)
if (!is.na(jh)) {
  tw$use_upper[jh] <- 1L; tw$upper_mean[jh] <- 0.010;  tw$upper_sd[jh] <- 0.003
  tw$use_lower[jh] <- 1L; tw$lower_mean[jh] <- 0.0005; tw$lower_sd[jh]  <- 0.001
}

# ---- 5. Common Stan data + per-variant additions ---------------------------
base_data <- list(
  G = G, Mc = nrow(cd_obs),
  cd_grid = nearest(cd_obs$relative_day), cd_n = cd_obs$n_comm, cd_N = cd_obs$N_deaths,
  c_hi = cd_prep$c_hi, c_lo = cd_prep$c_lo,
  s_grid = s_success_grid, s_ref = s_ref,
  J = J,
  abs_min = param_meta$abs_min, abs_max = param_meta$abs_max,
  lower_floor = param_meta$lower_floor, lower_cap = param_meta$lower_cap, upper_cap = param_meta$upper_cap,
  lb_prior_mean = param_meta$lb_prior_mean, lb_prior_sd = param_meta$lb_prior_sd,
  ub_prior_mean = param_meta$ub_prior_mean, ub_prior_sd = param_meta$ub_prior_sd,
  increases = param_meta$increases,
  N = nrow(fit_anchors), param_id = fit_anchors$param_id,
  y_obs = fit_anchors$value_used, obs_sd_mult = fit_anchors$obs_sd_mult,
  anchor_grid = fit_anchors$anchor_grid,
  sigma_frac_prior_meanlog = log(0.12), sigma_frac_prior_sdlog = 0.60,
  use_upper_tweak = tw$use_upper, upper_tweak_mean = tw$upper_mean, upper_tweak_sd = tw$upper_sd,
  use_lower_tweak = tw$use_lower, lower_tweak_mean = tw$lower_mean, lower_tweak_sd = tw$lower_sd
)

mods <- list(
  rw     = cmdstan_model(file.path(DIR_STAN, "modelC_reach_rw.stan")),
  spline = cmdstan_model(file.path(DIR_STAN, "modelC_reach_spline.stan"))
)
B_spline <- bs(grid_day, df = SPLINE_DF, intercept = TRUE)   # G x K basis (spline model)

# ---- 6. Fit one variant ----------------------------------------------------
fit_variant <- function(model, ipc_clock) {
  message("\n==== Model C [", model, " | IPC on ", ipc_clock, "] ====")
  sd <- base_data
  sd$clock_type <- unname(clock_for(ipc_clock))
  if (model == "rw") {
    sd$grid_dt <- GRID_DT; sd$rw_sigma_prior_sd <- RW_SIGMA_PSD
  } else {
    sd$K <- ncol(B_spline); sd$B <- unclass(B_spline); sd$spline_sigma_prior_sd <- SPL_SIGMA_PSD
  }
  fit <- mods[[model]]$sample(data = sd, seed = 123,
    chains = SAMP$chains, parallel_chains = SAMP$parallel_chains,
    iter_warmup = SAMP$iter_warmup, iter_sampling = SAMP$iter_sampling,
    adapt_delta = SAMP$adapt_delta, max_treedepth = SAMP$max_treedepth, refresh = SAMP$refresh)
  print(fit$diagnostic_summary())

  # Assemble the 6 native-unit parameter curves on the grid:
  th <- fit$summary("theta_pred") %>%
    mutate(j = as.integer(str_match(variable, "\\[([0-9]+),([0-9]+)\\]")[, 2]),
           g = as.integer(str_match(variable, "\\[([0-9]+),([0-9]+)\\]")[, 3]),
           parameter = ENDPOINT_PARAMS[j], relative_day = grid_day[g]) %>%
    select(parameter, relative_day, mean, q5, q95)
  ufc <- fit$summary("uf_comm_pred") %>%
    mutate(g = as.integer(str_match(variable, "\\[([0-9]+)\\]")[, 2]),
           parameter = "p_unsafe_funeral_comm", relative_day = grid_day[g]) %>%
    select(parameter, relative_day, mean, q5, q95)
  reach <- fit$summary(c("reach_clock", "reach_prob", "comm_death_p")) %>%
    mutate(metric = str_remove(variable, "\\[.*"),
           g = as.integer(str_match(variable, "\\[([0-9]+)\\]")[, 2]),
           relative_day = grid_day[g]) %>%
    select(metric, relative_day, mean, q5, q95)

  list(model = model, ipc_clock = ipc_clock,
       curves = bind_rows(th, ufc), reach = reach, fit = fit)
}

# ---- 7. Fit all four variants ----------------------------------------------
fits <- lapply(seq_len(nrow(VARIANTS)), function(i)
  fit_variant(VARIANTS$model[i], VARIANTS$ipc_clock[i]))
names(fits) <- sprintf("%s_IPC-%s", VARIANTS$model, VARIANTS$ipc_clock)

saveRDS(fits, file.path(DIR_PROCESSED, "DRC_QCurve", "DRC_QCurve_ModelC_fits.rds"))
message("Saved DRC_QCurve_ModelC_fits.rds")

# ---- 8. Overlay the four variants' parameter curves (display only) ---------
curve_df <- bind_rows(lapply(names(fits), function(nm)
  mutate(fits[[nm]]$curves, variant = nm))) %>%
  mutate(panel = factor(PANEL_LOOKUP[parameter], levels = unname(PANEL_LOOKUP)))

print(
  ggplot(curve_df, aes(relative_day, mean, colour = variant, fill = variant)) +
    geom_ribbon(aes(ymin = q5, ymax = q95), colour = NA, alpha = 0.12) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~ panel, scales = "free_y", ncol = 2) +
    labs(title = "DRC Model C: parameter curves, four variants",
         subtitle = "reach estimated (RW vs spline) x latent_IPC on SDB vs reach clock",
         x = "Relative outbreak day", y = NULL) +
    theme_bw(base_size = 11) + theme(legend.position = "top")
)

# ---- 9. Reach posterior vs the community-death data (the key check) ---------
reach_df <- bind_rows(lapply(names(fits), function(nm)
  filter(fits[[nm]]$reach, metric == "comm_death_p") %>% mutate(variant = nm)))
print(
  ggplot(reach_df, aes(relative_day, mean, colour = variant, fill = variant)) +
    geom_ribbon(aes(ymin = q5, ymax = q95), colour = NA, alpha = 0.12) +
    geom_line(linewidth = 0.8) +
    geom_point(data = cd_obs, aes(relative_day, n_comm / N_deaths),
               inherit.aes = FALSE, colour = "black", size = 1.6, alpha = 0.7) +
    labs(title = "Latent community-death proportion vs data (posterior-predictive check)",
         x = "Relative outbreak day", y = "Community-death proportion") +
    theme_bw(base_size = 11) + theme(legend.position = "top")
)
