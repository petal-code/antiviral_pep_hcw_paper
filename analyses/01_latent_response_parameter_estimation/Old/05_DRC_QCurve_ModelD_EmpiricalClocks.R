# ============================================================================
# 05_DRC_QCurve_ModelD_EmpiricalClocks.R   (TRIAL / SCAFFOLD -- untested)
# ----------------------------------------------------------------------------
# Model D: BOTH response clocks are SUPPLIED EMPIRICALLY -- no latent curve.
#
# Where Model C (script 04) estimated the reach curve R(t) latently (a random
# walk / spline on the community-death proportion, smoothed by its prior and so
# unable to follow the early community-death shoulder), Model D builds R(t)
# DIRECTLY from the community-death data, exactly as the SDB success curve s(t)
# is already supplied. Stan then only solves the original-methodology endpoint
# problem (estimate lower_j / upper_j from the literature anchors) against the
# two fixed empirical clocks. The reach curve tracks the data as tightly as the
# empirical construction allows.
#
#   delay_hosp, p_hosp, p_ETU      -> reach clock  R_clock(t)   (EMPIRICAL, supplied)
#   latent_IPC                     -> SDB clock (default) OR reach (toggle)
#   p_unsafe_funeral_hosp          -> SDB clock                  (near-zero)
#   p_unsafe_funeral_comm          -> DETERMINISTIC 1 - R_prob(t)*s(t)
#
# Fits FOUR variants: {IPC on SDB, IPC on reach} x {conflict, conflict++ s(t)}.
# (No RW/spline dimension -- the reach curve is not fitted.)
#
# Prereqs (run first): 00_DataPreparation_and_Cleaning.R (-> DRC_QCurve_PreppedData.rds)
#                      00b_DataPreparation_CommunityDeaths.R (-> DRC_CommunityDeaths_Prepped.rds)
#
# >>> First scaffold to TRIAL the model: not run here (no CmdStan in the build
#     env). Expect to debug priors / grid / sampler on first run. Key knobs below.
# ============================================================================

suppressPackageStartupMessages({
  library(cmdstanr); library(dplyr); library(tidyr); library(stringr)
  library(tibble); library(ggplot2)
})
source(here::here("analyses", "01_latent_response_parameter_estimation", "helpers.R"))
set.seed(123)

# ---- knobs (EDIT) ----------------------------------------------------------
# REACH_SMOOTH_K: optional light centred rolling mean over the community-death
# cell series. 1 = NONE (the empirical curve passes through every pooled cell and
# fully captures the early shoulder). Raise it ONLY if the raw curve is too
# jagged to use -- it then begins to wash the shoulder back out, defeating the
# point of going empirical. Because there is no binomial weighting any more, small
# -N cells now have FULL leverage: their spikes are raw sampling noise. The
# denominator-sized points in section 9 show which cells those are.
GRID_DT        <- 14L          # grid spacing (days); G = horizon/GRID_DT + 1
REACH_SMOOTH_K <- 1L           # rolling-mean window over cells (1 = no smoothing)
SAMP <- list(chains = 4, parallel_chains = 4, iter_warmup = 1000,
             iter_sampling = 1500, adapt_delta = 0.95, max_treedepth = 12, refresh = 200)

SCENARIOS <- c("conflict", "conflict_plusplus")
VARIANTS  <- expand.grid(ipc_clock = c("sdb", "reach"),
                         scenario  = SCENARIOS, stringsAsFactors = FALSE)

# ---- 1. Inputs -------------------------------------------------------------
drc_prep <- readRDS(file.path(DIR_PROCESSED, "DRC_QCurve", "DRC_QCurve_PreppedData.rds"))
anchors  <- drc_prep$anchors
cd_prep  <- readRDS(file.path(DIR_PROCESSED, "DRC_QCurve", "DRC_CommunityDeaths_Prepped.rds"))
cd_obs   <- cd_prep$obs

qseries_for <- function(scenario) {
  switch(scenario,
    conflict          = drc_prep$conflict_qseries,
    conflict_plusplus = drc_prep$conflict_plusplus_qseries,
    stop("unknown scenario: ", scenario))
}

ENDPOINT_PARAMS <- setdiff(PARAM_LEVELS, "p_unsafe_funeral_comm")
clock_for <- function(ipc_clock) c(
  delay_hosp            = 1L,
  p_hosp                = 1L,
  p_ETU                 = 1L,
  latent_IPC            = if (ipc_clock == "reach") 1L else 2L,
  p_unsafe_funeral_hosp = 2L
)[ENDPOINT_PARAMS]

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

# ---- 2. Grid + EMPIRICAL reach curve from the community-death data ----------
grid_day <- seq(0L, HORIZON_DAYS, by = GRID_DT)
G        <- length(grid_day)
nearest  <- function(day) vapply(day, function(d) which.min(abs(grid_day - d)), integer(1))

# Pool the community-death observations into grid cells (denominator-weighted:
# pooled proportion = sum(n_comm) / sum(N_deaths) within the cell), then linearly
# interpolate across cells onto the full grid. This is the EMPIRICAL reach driver:
# no latent model, no Bayesian smoothing. rule = 2 (in make_interp) holds the
# curve flat before the first / after the last observed cell.
cd_cell <- cd_obs %>%
  mutate(gi = nearest(relative_day)) %>%
  group_by(gi) %>%
  summarise(n = sum(n_comm), N = sum(N_deaths), .groups = "drop") %>%
  arrange(gi) %>%
  mutate(p = n / N, day = grid_day[gi])
p_cell <- if (REACH_SMOOTH_K > 1L) rolling_mean_centered(cd_cell$p, k = REACH_SMOOTH_K) else cd_cell$p

p_comm_emp  <- clip01(make_interp(cd_cell$day, p_cell)(grid_day))
c_hi        <- cd_prep$c_hi
c_lo        <- cd_prep$c_lo
R_prob_emp  <- clip01(1 - p_comm_emp)
R_clock_emp <- clip01((c_hi - p_comm_emp) / (c_hi - c_lo))

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

# ---- 5. Common Stan data (reach clock is fixed; s(t) injected per-fit) ------
base_data <- list(
  G = G,
  R_clock = R_clock_emp, R_prob = R_prob_emp,    # SUPPLIED empirical reach clocks
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

mod <- cmdstan_model(file.path(DIR_STAN, "modelD_empiricalClocks.stan"))

sdb_on_grid <- function(scenario) {
  qs <- qseries_for(scenario)
  list(s_grid = clip01(make_interp(qs$relative_day, qs$success_smoothed)(grid_day)),
       s_ref  = max(qs$success_smoothed, na.rm = TRUE))
}

# ---- 6. Fit one variant ----------------------------------------------------
fit_variant <- function(ipc_clock, scenario) {
  message("\n==== Model D [IPC on ", ipc_clock, " | ", scenario, "] ====")
  sdb <- sdb_on_grid(scenario)
  sd  <- base_data
  sd$s_grid     <- sdb$s_grid
  sd$s_ref      <- sdb$s_ref
  sd$clock_type <- unname(clock_for(ipc_clock))
  fit <- mod$sample(data = sd, seed = 123,
    chains = SAMP$chains, parallel_chains = SAMP$parallel_chains,
    iter_warmup = SAMP$iter_warmup, iter_sampling = SAMP$iter_sampling,
    adapt_delta = SAMP$adapt_delta, max_treedepth = SAMP$max_treedepth, refresh = SAMP$refresh)
  print(fit$diagnostic_summary())

  th <- fit$summary("theta_pred") %>%
    mutate(j = as.integer(str_match(variable, "\\[([0-9]+),([0-9]+)\\]")[, 2]),
           g = as.integer(str_match(variable, "\\[([0-9]+),([0-9]+)\\]")[, 3]),
           parameter = ENDPOINT_PARAMS[j], relative_day = grid_day[g]) %>%
    select(parameter, relative_day, mean, q5, q95)
  ufc <- fit$summary("uf_comm_pred") %>%
    mutate(g = as.integer(str_match(variable, "\\[([0-9]+)\\]")[, 2]),
           parameter = "p_unsafe_funeral_comm", relative_day = grid_day[g]) %>%
    select(parameter, relative_day, mean, q5, q95)

  list(ipc_clock = ipc_clock, scenario = scenario,
       curve_type = sprintf("IPC-%s", ipc_clock),
       curves = bind_rows(th, ufc), fit = fit)
}

# ---- 7. Fit all four variants ----------------------------------------------
fits <- lapply(seq_len(nrow(VARIANTS)), function(i)
  fit_variant(VARIANTS$ipc_clock[i], VARIANTS$scenario[i]))
names(fits) <- sprintf("IPC-%s_%s", VARIANTS$ipc_clock, VARIANTS$scenario)

saveRDS(list(fits = fits, p_comm_emp = p_comm_emp, R_clock_emp = R_clock_emp,
             R_prob_emp = R_prob_emp, grid_day = grid_day),
        file.path(DIR_PROCESSED, "DRC_QCurve", "DRC_QCurve_ModelD_fits.rds"))
message("Saved DRC_QCurve_ModelD_fits.rds")

# ---- 8. Parameter curves vs the DATA (display only) ------------------------
# colour = IPC-clock; linetype = scenario; orange = literature anchors;
# grey = SDB community proxy (as in 02 / 04).
curve_df <- bind_rows(lapply(names(fits), function(nm) {
  f <- fits[[nm]]
  mutate(f$curves, variant = nm, curve_type = f$curve_type, scenario = f$scenario)
})) %>%
  mutate(panel    = factor(PANEL_LOOKUP[parameter], levels = unname(PANEL_LOOKUP)),
         scenario = factor(scenario, levels = SCENARIOS))

anchor_plot_df <- anchors %>%
  filter(parameter %in% ENDPOINT_PARAMS) %>%
  transmute(relative_day, value_used,
            panel = factor(PANEL_LOOKUP[parameter], levels = unname(PANEL_LOOKUP)))

ufc_proxy_df <- qseries_for("conflict") %>%
  filter(n_eligible_sum > 0) %>%
  transmute(relative_day, value = unsafe_funeral_comm_proxy,
            panel = factor(PANEL_LOOKUP[["p_unsafe_funeral_comm"]], levels = unname(PANEL_LOOKUP)))

print(
  ggplot(curve_df, aes(relative_day, mean, colour = curve_type,
                       linetype = scenario, group = variant)) +
    geom_line(linewidth = 0.8) +
    geom_point(data = ufc_proxy_df, aes(relative_day, value),
               inherit.aes = FALSE, colour = "grey55", size = 1, alpha = 0.7) +
    geom_point(data = anchor_plot_df, aes(relative_day, value_used),
               inherit.aes = FALSE, colour = "#ff7f0e", size = 2) +
    facet_wrap(~ panel, scales = "free_y", ncol = 2) +
    labs(title = "DRC Model D: parameter curves vs data (empirical reach + SDB clocks)",
         subtitle = paste("colour = IPC-clock; linetype = scenario;",
                          "orange = literature anchors; grey = SDB community proxy"),
         x = "Relative outbreak day", y = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position = "top", strip.text = element_text(face = "bold"))
)

# ---- 9. The EMPIRICAL reach curve vs the community-death data ---------------
# This is the whole point: R(t) is built directly from the data, so the curve
# passes through it (no posterior interval -- it is not estimated). Point AREA =
# denominator (N_deaths), so it is visible which cells the curve is following.
emp_df <- tibble(relative_day = grid_day, p_comm = p_comm_emp)
print(
  ggplot(emp_df, aes(relative_day, p_comm)) +
    geom_line(colour = "#1f77b4", linewidth = 1) +
    geom_point(data = cd_obs, aes(relative_day, n_comm / N_deaths, size = N_deaths),
               inherit.aes = FALSE, shape = 21, fill = "grey25",
               colour = "white", stroke = 0.3, alpha = 0.85) +
    scale_size_area(max_size = 7, name = "deaths (N)") +
    labs(title = "Empirical community-death proportion fed into Model D (no latent fit)",
         subtitle = sprintf("blue = empirical curve through the data (REACH_SMOOTH_K = %d); point area = denominator",
                            REACH_SMOOTH_K),
         x = "Relative outbreak day", y = "Community-death proportion") +
    theme_bw(base_size = 11)
)
