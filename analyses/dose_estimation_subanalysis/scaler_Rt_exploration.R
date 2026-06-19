# ============================================================================
# 005_scalar_Rt_exploration.R
# ----------------------------------------------------------------------------
# Computes analytic instantaneous Rt for all 4 NPI scenarios from 01, using
# the PI-specified efficacy values directly (no scalar approach). Quick sanity
# check before running the full FIBER simulations in 004.
#
# Inputs  : outputs/dose_q_curve_extrapolation_scenarios.rds  (from 01)
# Outputs : outputs/005_scalar_Rt_exploration.png
#           outputs/005_scalar_Rt_exploration.csv
# ============================================================================

suppressPackageStartupMessages({
  library(here)
  library(fiber)
  library(dplyr)
  library(ggplot2)
})

source(here("analyses", "dose_estimation_subanalysis", "helpers.R"))
source(here("functions", "setup_model_parameters.R"))
source(here("functions", "calculate_model_approx_r0.R"))
source(here("functions", "calculate_model_approx_rt.R"))

# ----------------------------------------------------------------------------
# 1. Configuration
# ----------------------------------------------------------------------------
R0_FIXED <- 1.60

# Efficacy parameters from PI.
SCALAR_OVERRIDES <- list(
  etu_efficacy                         = 0.84,
  general_hospital_quarantine_efficacy = 0.68,
  ppe_efficacy                         = 0.84,
  safe_funeral_efficacy                = 0.88,
  seeding_cases    = 3L,
  check_final_size = 40000L
)

FUNERAL_FRAC        <- 0.25
EPIDEMIC_START_DATE <- as.Date("2026-02-27")
MATRIX_HORIZON      <- 730L
PHEIC_DATES         <- as.Date(c("2026-05-14", "2026-05-18"))

NPI_SPECS <- list(
  prob_hosp                = list(q0 = 0.00, q1 = 0.80),
  delay_hosp               = list(q0 = 6.00, q1 = 1.50),
  prop_etu                 = list(q0 = 0.00, q1 = 0.90),
  safe_funeral_prop        = list(q0 = 0.00, q1 = 0.90),
  unsafe_funeral_prop_hosp = list(q0 = 0.05, q1 = 0.01),
  ppe_coverage             = list(q0 = 0.00, q1 = 0.90)
)
UNSAFE_FUNERAL_ETU <- 0.0

RT_TIMES <- 0:365L
RT_MC_N  <- 50000L
RT_SEED  <- 1L

scenario_colours <- c(
  linear_to_90 = "#1b9e77",
  logistic     = "#7570b3",
  flat         = "#d95f02",
  conflict     = "#e7298a"
)
scenario_labels <- c(
  linear_to_90 = "1. Linear to 90%",
  logistic     = "2. Logistic projection",
  flat         = "3. Flat at last value",
  conflict     = "4. Conflict episode"
)

# ----------------------------------------------------------------------------
# 2. Load Q-curve scenarios
# ----------------------------------------------------------------------------
scen_path <- file.path(DIR_OUT, "dose_q_curve_extrapolation_scenarios.rds")
if (!file.exists(scen_path))
  stop("Run 01_fit_dose_q_curve.R first to generate ", scen_path, call. = FALSE)

all_scen      <- readRDS(scen_path)
all_scen$date <- as.Date(all_scen$date)
scen_names    <- c("linear_to_90", "logistic", "flat", "conflict")

epi_start    <- as.Date(EPIDEMIC_START_DATE)
q_first_date <- min(all_scen$date)
offset_days  <- as.integer(q_first_date - epi_start)
day_to_date  <- function(d) epi_start + d

pheic_vlines <- lapply(PHEIC_DATES, function(d)
  geom_vline(xintercept = d, linetype = "dashed", colour = "grey40"))

# ----------------------------------------------------------------------------
# 3. Build NPI scenario matrices
# ----------------------------------------------------------------------------
lin         <- function(spec, q) spec$q0 + (spec$q1 - spec$q0) * q
matrix_days <- 0:MATRIX_HORIZON

build_matrix_df <- function(sn) {
  sel    <- all_scen[as.character(all_scen$scenario) == sn, , drop = FALSE]
  q_sim  <- as.integer(sel$date - epi_start)
  q_vals <- approx(q_sim, sel$q, xout = matrix_days, rule = 2)$y
  q_vals[matrix_days < offset_days] <- 0
  q_vals <- clip01(q_vals)
  safe_fp   <- clip01(lin(NPI_SPECS$safe_funeral_prop,        q_vals))
  unsafe_fh <- clip01(lin(NPI_SPECS$unsafe_funeral_prop_hosp, q_vals))
  data.frame(
    scenario                 = sn, scenario_label = sn,
    relative_day             = matrix_days,
    prob_hosp                = clip01(lin(NPI_SPECS$prob_hosp,    q_vals)),
    delay_hosp               = pmax(lin(NPI_SPECS$delay_hosp,     q_vals), 0.01),
    prob_unsafe_funeral_comm = clip01(1 - safe_fp),
    prob_unsafe_funeral_hosp = unsafe_fh,
    prob_unsafe_funeral_etu  = clip01(rep(UNSAFE_FUNERAL_ETU, length(matrix_days))),
    prop_etu                 = clip01(lin(NPI_SPECS$prop_etu,     q_vals)),
    ipc_helper               = clip01(lin(NPI_SPECS$ppe_coverage, q_vals)),
    q_value                  = q_vals,
    stringsAsFactors = FALSE
  )
}

mat_dir <- file.path(DIR_OUT, "005_scenario_matrices")
dir.create(mat_dir, recursive = TRUE, showWarnings = FALSE)

matrix_csvs <- setNames(vapply(scen_names, function(sn) {
  df  <- build_matrix_df(sn)
  pth <- file.path(mat_dir, sprintf("matrix_%s.csv", sn))
  write.csv(df, pth, row.names = FALSE)
  pth
}, character(1)), scen_names)

# ----------------------------------------------------------------------------
# 4. Compute Rt for all 4 scenarios
# ----------------------------------------------------------------------------
ref_sm <- read_scenario_matrix(matrix_csvs[[scen_names[1]]])
mp_ref <- make_model_parameters(scen_names[1], ref_sm, overrides = SCALAR_OVERRIDES)
inv    <- compute_R0_invariants(mp_ref$args, n = 50000L, seed = 42L)
D      <- D_from_invariants(inv, mp_ref$args$etu_efficacy,
                            mp_ref$args$general_hospital_quarantine_efficacy)
F_fun  <- F_from_invariants(inv, mp_ref$args$safe_funeral_efficacy)
means  <- solve_offspring_means(R0_FIXED, FUNERAL_FRAC, D, F_fun)
message(sprintf("R0 = %.2f: mn_genPop = %.4f, mn_funeral = %.4f",
                R0_FIXED, means$mn_genPop, means$mn_funeral))

message("Computing Rt for each scenario...")
rt_all <- do.call(rbind, lapply(scen_names, function(sn) {
  sm <- read_scenario_matrix(matrix_csvs[[sn]])
  mp <- make_model_parameters(sn, sm, overrides = SCALAR_OVERRIDES)
  a  <- mp$args
  a$mn_offspring_genPop  <- means$mn_genPop
  a$mn_offspring_funeral <- means$mn_funeral
  a$seed                 <- NULL
  rt <- Rt_curve_single_type(a, times = RT_TIMES, n = RT_MC_N, seed = RT_SEED)
  data.frame(scenario = sn, day = rt$time,
             R_inst = rt$R_inst, R_case = rt$R_case,
             stringsAsFactors = FALSE)
}))

rt_all$date     <- day_to_date(rt_all$day)
rt_all$scenario <- factor(rt_all$scenario, levels = scen_names)

write.csv(rt_all, file.path(DIR_OUT, "005_scalar_Rt_exploration.csv"), row.names = FALSE)

# ----------------------------------------------------------------------------
# 5. Plot
# ----------------------------------------------------------------------------
p <- ggplot(rt_all, aes(date, R_inst, colour = scenario, group = scenario)) +
  geom_hline(yintercept = 1, colour = "black", linewidth = 0.5, linetype = "dashed") +
  pheic_vlines +
  geom_line(linewidth = 0.85) +
  scale_colour_manual(values = scenario_colours, labels = scenario_labels,
                      name = "NPI scenario") +
  scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") +
  labs(
    title    = sprintf("Instantaneous Rt by NPI scenario (R0 = %.2f)", R0_FIXED),
    subtitle = sprintf("ETU=0.84, GenHosp=0.68, PPE=0.84, SafeFuneral=0.88\nDashed black = Rt 1; dashed grey = PHEIC dates (%s, %s)",
                       format(PHEIC_DATES[1], "%d %b"), format(PHEIC_DATES[2], "%d %b")),
    x = "Date", y = expression(R[t] ~ "(instantaneous)")
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text.x   = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

ggsave(file.path(DIR_OUT, "005_scalar_Rt_exploration.png"), p,
       width = 10, height = 6, dpi = 150)
print(p)
message("005_scalar_Rt_exploration.R complete. Output in outputs/.")
