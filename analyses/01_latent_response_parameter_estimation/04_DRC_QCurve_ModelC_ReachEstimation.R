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
# Fits EIGHT variants for comparison:
#   {RW, spline} x {IPC on SDB, IPC on reach} x {conflict, conflict++ SDB curve}.
# Reach is estimated from the community-death data ONLY, so it is invariant to
# the IPC-clock toggle and to the scenario (which only swaps the SUPPLIED s(t));
# the scenario therefore moves only the SDB-clock parameter curves, not reach.
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
# NB: RW_SIGMA / SPL_SIGMA / SPLINE_DF control how WIGGLY the latent reach curve
# is allowed to be. They were loosened (RW 0.30->0.60, spline 0.50->1.00, df
# 14->20) so the curve has the freedom to rise into the early community-death
# "shoulder" around day ~150-230 instead of being ironed flat by the smoothness
# prior. Whether it SHOULD chase those points depends on their denominators -
# see the denominator-sized points in section 9.
GRID_DT       <- 14L          # reach-grid spacing (days); G = horizon/GRID_DT + 1
SPLINE_DF     <- 20L          # B-spline basis size (spline model)
RW_SIGMA_PSD  <- 0.60         # half-normal scale on the RW step sd (per sqrt-day, logit)
SPL_SIGMA_PSD <- 1.00         # half-normal scale on the P-spline 2nd-diff sd
SAMP <- list(chains = 4, parallel_chains = 4, iter_warmup = 1500,
             iter_sampling = 1500, adapt_delta = 0.98, max_treedepth = 13, refresh = 200)

# Two supplied-SDB scenarios (the only thing the scenario changes is s(t)).
SCENARIOS <- c("conflict", "conflict_plusplus")
VARIANTS  <- expand.grid(model = c("rw", "spline"),
                         ipc_clock = c("sdb", "reach"),
                         scenario  = SCENARIOS, stringsAsFactors = FALSE)

# ---- 1. Inputs -------------------------------------------------------------
drc_prep <- readRDS(file.path(DIR_PROCESSED, "DRC_QCurve", "DRC_QCurve_PreppedData.rds"))
anchors  <- drc_prep$anchors
cd_prep  <- readRDS(file.path(DIR_PROCESSED, "DRC_QCurve", "DRC_CommunityDeaths_Prepped.rds"))
cd_obs   <- cd_prep$obs

# Pick the SUPPLIED SDB curve for a scenario (this is the ONLY scenario lever).
qseries_for <- function(scenario) {
  switch(scenario,
    conflict          = drc_prep$conflict_qseries,
    conflict_plusplus = drc_prep$conflict_plusplus_qseries,
    stop("unknown scenario: ", scenario))
}

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

# Supplied SDB ACTUAL success on the grid for a given scenario (success, not the
# normalised q_value):  s_ref = max success => Sclk = s/s_ref = q_value (Model B),
# and the burial term uses the actual success s.
sdb_on_grid <- function(scenario) {
  qs <- qseries_for(scenario)
  list(s_grid = clip01(make_interp(qs$relative_day, qs$success_smoothed)(grid_day)),
       s_ref  = max(qs$success_smoothed, na.rm = TRUE))
}

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
# s_grid / s_ref are scenario-specific and so are injected per-fit (section 6).
base_data <- list(
  G = G, Mc = nrow(cd_obs),
  cd_grid = nearest(cd_obs$relative_day), cd_n = cd_obs$n_comm, cd_N = cd_obs$N_deaths,
  c_hi = cd_prep$c_hi, c_lo = cd_prep$c_lo,
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
fit_variant <- function(model, ipc_clock, scenario) {
  message("\n==== Model C [", model, " | IPC on ", ipc_clock, " | ", scenario, "] ====")
  sdb <- sdb_on_grid(scenario)
  sd  <- base_data
  sd$s_grid     <- sdb$s_grid
  sd$s_ref      <- sdb$s_ref
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

  list(model = model, ipc_clock = ipc_clock, scenario = scenario,
       curve_type = sprintf("%s_IPC-%s", model, ipc_clock),
       curves = bind_rows(th, ufc), reach = reach, fit = fit)
}

# ---- 7. Fit all eight variants ---------------------------------------------
fits <- lapply(seq_len(nrow(VARIANTS)), function(i)
  fit_variant(VARIANTS$model[i], VARIANTS$ipc_clock[i], VARIANTS$scenario[i]))
names(fits) <- sprintf("%s_IPC-%s_%s", VARIANTS$model, VARIANTS$ipc_clock, VARIANTS$scenario)

saveRDS(fits, file.path(DIR_PROCESSED, "DRC_QCurve", "DRC_QCurve_ModelC_fits.rds"))
message("Saved DRC_QCurve_ModelC_fits.rds")

# ---- 8. Overlay the variants' parameter curves vs the DATA (display only) ---
# Colour = model x IPC-clock; linetype = scenario. Ribbons are dropped here
# (eight curves) to keep the panels legible; overlaid on top are the literature
# anchors the endpoints were fit to (orange) and the SDB-derived community
# unsafe-funeral proxy (grey) -- the same data overlays used in 02's plots.
curve_df <- bind_rows(lapply(names(fits), function(nm) {
  f <- fits[[nm]]
  mutate(f$curves, variant = nm, curve_type = f$curve_type, scenario = f$scenario)
})) %>%
  mutate(panel    = factor(PANEL_LOOKUP[parameter], levels = unname(PANEL_LOOKUP)),
         scenario = factor(scenario, levels = SCENARIOS))

# Literature anchors the endpoints were fit to (shared across variants).
anchor_plot_df <- anchors %>%
  filter(parameter %in% ENDPOINT_PARAMS) %>%
  transmute(relative_day, value_used,
            panel = factor(PANEL_LOOKUP[parameter], levels = unname(PANEL_LOOKUP)))

# SDB community unsafe-funeral data points (1 - success) behind the UFC curve.
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
    labs(title = "DRC Model C: parameter curves vs data, eight variants",
         subtitle = paste("colour = model x IPC-clock; linetype = scenario;",
                          "orange = literature anchors; grey = SDB community proxy"),
         x = "Relative outbreak day", y = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position = "top", strip.text = element_text(face = "bold"))
)

# ---- 9. Reach posterior vs the community-death data (the key check) ---------
# Reach is estimated from the community-death binomial ONLY, so it is identical
# across the IPC-clock toggle and the scenario; only RW vs spline differ. We
# therefore draw one representative fit per smoother (conflict / IPC-sdb) and
# overlay the data with point AREA proportional to the denominator (N_deaths):
# big points are well-supported, tiny points are near-noise. This is the lens
# for "is the early shoulder a real peak the curve should chase, or sampling
# noise it is right to smooth through?".
reach_df <- bind_rows(lapply(names(fits), function(nm) {
  f <- fits[[nm]]
  if (f$scenario != "conflict" || f$ipc_clock != "sdb") return(NULL)
  filter(f$reach, metric == "comm_death_p") %>% mutate(model = f$model)
}))

print(
  ggplot(reach_df, aes(relative_day, mean, colour = model, fill = model)) +
    geom_ribbon(aes(ymin = q5, ymax = q95), colour = NA, alpha = 0.15) +
    geom_line(linewidth = 0.9) +
    geom_point(data = cd_obs, aes(relative_day, n_comm / N_deaths, size = N_deaths),
               inherit.aes = FALSE, shape = 21, fill = "grey25",
               colour = "white", stroke = 0.3, alpha = 0.85) +
    scale_size_area(max_size = 7, name = "deaths (N)") +
    labs(title = "Latent community-death proportion vs data (posterior-predictive check)",
         subtitle = "reach is scenario- / clock-invariant: only RW vs spline differ; point area = denominator",
         x = "Relative outbreak day", y = "Community-death proportion") +
    theme_bw(base_size = 11) + theme(legend.position = "top")
)
