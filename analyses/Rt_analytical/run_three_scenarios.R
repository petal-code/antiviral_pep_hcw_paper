## =============================================================================
## run_three_scenarios.R
##
## Test-drive the analytical single-type R(t) curves
## (functions/calculate_model_approx_rt.R) on three scenarios, using the STANDARD
## (default) model parameters -- i.e. DEFAULT_SCALAR_INPUTS via make_base_args(),
## which already carry the offspring means (mn_offspring_genPop = 1.25,
## mn_offspring_funeral = 0.25) and the standard efficacies (etu_efficacy = 0.90,
## general_hospital_quarantine_efficacy = 0.30, safe_funeral_efficacy = 0.80).
## No ABC posterior / no R0 inversion -- just the standard parameter set on each
## scenario's time-varying response curves.
##
## Scenarios:
##   * West Africa             -> Worst_WestAfrica
##   * DRC Conflict (smoothed) -> Middle_DRC_ConflictSmoothed
##   * DRC Conflict ++         -> Middle_DRC_ConflictSmoothed_PlusPlus
##
## Run from the repo root:
##   Rscript analyses/Rt_analytical/run_three_scenarios.R
##
## Uses the canonical make_model_parameters() path when fiber is installed; if
## fiber is absent it falls back to a base-R curve loader that mirrors
## build_time_varying_args(), so the script runs either way.
## =============================================================================

## ---- config -----------------------------------------------------------------
SCENARIO_CSV <- "data-processed/final_six_scenario_values_original_approach.csv"
SCENARIOS <- c(
  "West Africa"             = "Worst_WestAfrica",
  "DRC Conflict (smoothed)" = "Middle_DRC_ConflictSmoothed",
  "DRC Conflict ++"         = "Middle_DRC_ConflictSmoothed_PlusPlus"
)
N_MC    <- 50000L
SEED    <- 1L
OUT_DIR <- "outputs/Rt_analytical"

## ---- locate repo root + source ----------------------------------------------
root <- if (file.exists("functions/calculate_model_approx_rt.R")) "." else
        stop("Run from the repo root (functions/ not found).")
source(file.path(root, "functions", "setup_model_parameters.R"))
source(file.path(root, "functions", "calculate_model_approx_rt.R"))
dir.create(file.path(root, OUT_DIR), showWarnings = FALSE, recursive = TRUE)

have_fiber <- requireNamespace("fiber", quietly = TRUE)
if (have_fiber) {
  library(fiber)
  scenario_matrix <- read_scenario_matrix(file.path(root, SCENARIO_CSV))
}

## ---- base-R curve loader (mirrors build_time_varying_args) -------------------
## Only used when fiber is not installed. Builds approxfun curves (rule = 2) with
## the breakpoints exposed via attr(, "times"), and applies the same ETU-weighted
## hospital unsafe-funeral blend as build_time_varying_args().
clip01 <- function(x) pmin(pmax(x, 0), 1)
load_curves_baseR <- function(scenario_id) {
  x <- utils::read.csv(file.path(root, SCENARIO_CSV), check.names = FALSE)
  x <- x[x$scenario == scenario_id, ]
  x <- x[order(x$relative_day), ]
  tt <- x$relative_day
  mk <- function(v) { f <- stats::approxfun(tt, v, method = "linear", rule = 2)
                      attr(f, "times") <- tt; f }
  p_uf_hosp <- clip01((1 - x$prop_etu) * x$prob_unsafe_funeral_hosp +
                      x$prop_etu * x$prob_unsafe_funeral_etu)
  list(
    prob_hospitalised_genPop     = mk(clip01(x$prob_hosp)),
    hospitalisation_delay_factor = mk(pmax(x$delay_hosp, 0.01)),
    p_unsafe_funeral_comm_genPop = mk(clip01(x$prob_unsafe_funeral_comm)),
    p_unsafe_funeral_hosp_genPop = mk(p_uf_hosp),
    prop_etu                     = mk(clip01(x$prop_etu))
  )
}

## Standard-parameter args for one scenario.
build_args <- function(scenario_id) {
  if (have_fiber) {
    make_model_parameters(scenario_id, scenario_matrix)$args     # standard defaults
  } else {
    c(make_base_args(), load_curves_baseR(scenario_id))          # standard defaults
  }
}

## first crossing of R = 1 (linear interp), or NA if it never drops below 1.
cross1 <- function(time, R) {
  i <- which(R < 1)[1]
  if (is.na(i) || i == 1) return(NA_real_)
  stats::approx(c(R[i - 1], R[i]), c(time[i - 1], time[i]), xout = 1)$y
}

## ---- compute curves ---------------------------------------------------------
curves  <- list()
summary <- data.frame()
for (k in seq_along(SCENARIOS)) {
  nm  <- names(SCENARIOS)[k]; id <- SCENARIOS[[k]]
  rt  <- Rt_curve_single_type(build_args(id), n = N_MC, seed = SEED)
  curves[[nm]] <- rt
  summary <- rbind(summary, data.frame(
    scenario      = nm,
    R0_inst_t0    = round(rt$R_inst[1], 3),
    R_case_t0     = round(rt$R_case[1], 3),
    min_R_inst    = round(min(rt$R_inst), 3),
    cross1_inst_d = round(cross1(rt$time, rt$R_inst), 1),
    cross1_case_d = round(cross1(rt$time, rt$R_case), 1),
    max_day       = max(rt$time)
  ))
}
cat("\nStandard-parameter single-type R(t) summary:\n")
print(summary, row.names = FALSE)

## ---- figure: one panel per scenario, R_inst vs R_case -----------------------
cols <- c(inst = "#0072B2", case = "#D55E00")
fig  <- file.path(root, OUT_DIR, "Rt_three_scenarios_standard_params.png")
grDevices::png(fig, width = 1500, height = 1200, res = 165)
op <- par(mfrow = c(3, 1), mar = c(4.0, 4.4, 2.6, 1), oma = c(0, 0, 2, 0))
ylim <- range(1, unlist(lapply(curves, function(d) c(d$R_inst, d$R_case))))
for (nm in names(curves)) {
  d <- curves[[nm]]
  plot(d$time, d$R_inst, type = "l", lwd = 2.4, col = cols["inst"], ylim = ylim,
       xlab = "Day since outbreak start", ylab = "Reproduction number",
       main = sprintf("%s  (%s)", nm, SCENARIOS[[nm]]))
  lines(d$time, d$R_case, lwd = 2.4, col = cols["case"])
  abline(h = 1, col = "grey50", lty = 2)
  if (nm == names(curves)[1])
    legend("topright", bty = "n", lwd = 2.4, col = cols,
           legend = c("R_inst (instantaneous)", "R_case (cohort)"))
}
mtext("Analytical single-type R(t), standard model parameters",
      outer = TRUE, cex = 1.05, font = 2)
par(op); grDevices::dev.off()
cat(sprintf("\nWrote figure: %s\n", fig))

## ---- also save the tidy curves as CSV for downstream plotting ---------------
long <- do.call(rbind, lapply(names(curves), function(nm) {
  d <- curves[[nm]]
  data.frame(scenario = nm, scenario_id = SCENARIOS[[nm]], time = d$time,
             R_inst = d$R_inst, R_case = d$R_case)
}))
csv_out <- file.path(root, OUT_DIR, "Rt_three_scenarios_standard_params.csv")
utils::write.csv(long, csv_out, row.names = FALSE)
cat(sprintf("Wrote curves CSV: %s\n\nDone.\n", csv_out))
