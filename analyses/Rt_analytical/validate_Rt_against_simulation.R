## =============================================================================
## validate_Rt_against_simulation.R
##
## Validation / demonstration harness for the analytical time-resolved
## reproduction-number curves (functions/calculate_model_approx_rt.R).
##
## It does three things, the last two of which run WITHOUT fiber:
##
##   PART A  (needs fiber)  -- the KEY external check. Run many
##     branching_process_main() replicates for a scenario (OBV gate OFF),
##     compute the empirical reproduction numbers with
##     fiber::compute_reproduction_number(out$tdf, type = "both",
##     bin_width = 7), and overlay:
##         analytical R_inst  vs  empirical R_instantaneous
##         analytical R_case  vs  empirical R_case
##     Expect tracking up to (i) the single-type approximation bias and (ii)
##     stochastic noise; trust only early / well-populated weekly bins.
##
##   PART B  (base R only)  -- GI-convolution cross-check. R_case(t) should be
##     approximately  sum_tau w(tau) * R_inst(t + tau), with w the genPop
##     generation-time pmf. It is NOT exact (fiber's hospitalisation gate is
##     infector-anchored, a mixed clock, not a single-clock renewal), so the
##     residual is reported as a diagnostic, not asserted.
##
##   FIGURE  (base R only)  -- save an example R_inst / R_case curve for the
##     scenario with the R = 1 crossing marked.
##
## Run from the repo root:
##   Rscript analyses/Rt_analytical/validate_Rt_against_simulation.R
## Without fiber installed, PART A is skipped automatically.
## =============================================================================

## ---- config -----------------------------------------------------------------
SCENARIO_CSV <- "data-processed/final_six_scenario_values_original_approach.csv"
SCENARIO_ID  <- "Middle_DRC_ConflictSmoothed_PlusPlus"
TARGET_R0    <- 1.8
PROP_FUNERAL <- 0.30
ETU_EFF      <- 0.90
GEN_HOSP_EFF <- 0.30
SAFE_FUN_EFF <- 0.80
N_MC         <- 50000L     # synthetic parents for the analytical MC
SEED         <- 1L
N_REPS       <- 100L       # fiber replicates (PART A)
BIN_WIDTH    <- 7L         # weekly, to match compute_reproduction_number default
OUT_DIR      <- "outputs/Rt_analytical"

## ---- locate repo root + source ----------------------------------------------
root <- if (file.exists("functions/calculate_model_approx_rt.R")) "." else
        stop("Run from the repo root (functions/ not found).")
source(file.path(root, "functions", "setup_model_parameters.R"))
source(file.path(root, "functions", "calculate_model_approx_r0.R"))
source(file.path(root, "functions", "calculate_model_approx_rt.R"))
dir.create(file.path(root, OUT_DIR), showWarnings = FALSE, recursive = TRUE)

have_fiber <- requireNamespace("fiber", quietly = TRUE)

## ---- build a concrete parameter set -----------------------------------------
## Time-varying curves: use the real fiber path when available, else a base-R
## loader that mirrors build_time_varying_args() (approxfun + attr "times").
clip01 <- function(x) pmin(pmax(x, 0), 1)

load_curves_baseR <- function(csv, scenario_id, root = ".") {
  x <- utils::read.csv(file.path(root, csv), check.names = FALSE)
  x <- x[x$scenario == scenario_id, ]
  x <- x[order(x$relative_day), ]
  tt <- x$relative_day
  mk <- function(v) { f <- stats::approxfun(tt, v, method = "linear", rule = 2)
                      attr(f, "times") <- tt; f }
  p_uf_hosp <- clip01((1 - x$prop_etu) * x$prob_unsafe_funeral_hosp +
                      x$prop_etu * x$prob_unsafe_funeral_etu)
  list(
    prob_hospitalised_genPop     = mk(clip01(x$prob_hosp)),
    prob_hospitalised_hcw        = mk(clip01(x$prob_hosp)),
    hospitalisation_delay_factor = mk(pmax(x$delay_hosp, 0.01)),
    p_unsafe_funeral_comm_genPop = mk(clip01(x$prob_unsafe_funeral_comm)),
    p_unsafe_funeral_hosp_genPop = mk(p_uf_hosp),
    p_unsafe_funeral_comm_hcw    = mk(clip01(x$prob_unsafe_funeral_comm)),
    p_unsafe_funeral_hosp_hcw    = mk(p_uf_hosp),
    prop_etu                     = mk(clip01(x$prop_etu)),
    ppe_coverage_hcw             = mk(clip01(x$ipc_helper))
  )
}

overrides <- list(etu_efficacy = ETU_EFF,
                  general_hospital_quarantine_efficacy = GEN_HOSP_EFF,
                  safe_funeral_efficacy = SAFE_FUN_EFF)

if (have_fiber) {
  library(fiber)
  sm <- read_scenario_matrix(file.path(root, SCENARIO_CSV))
  mp <- make_model_parameters(SCENARIO_ID, sm, overrides = overrides)
  args <- mp$args
} else {
  args <- c(make_base_args(overrides = overrides),
            load_curves_baseR(SCENARIO_CSV, SCENARIO_ID, root))
}

## offspring means from the t=0 R0 inversion (same machinery the ABC fit uses)
means <- solve_offspring_means_for_R0(TARGET_R0, args, PROP_FUNERAL,
                                      n = N_MC, seed = SEED)
args$mn_offspring_genPop  <- means$mn_offspring_genPop_required
args$mn_offspring_funeral <- means$mn_offspring_funeral_required

## ---- analytical curves on a daily grid --------------------------------------
rt <- Rt_curve_single_type(args, n = N_MC, seed = SEED)   # daily grid by default
cat(sprintf("Analytical curves: %d daily grid points, t in [%g, %g]\n",
            nrow(rt), min(rt$time), max(rt$time)))
cat(sprintf("R_inst(0) = %.3f  (target R0 = %.2f);  R_case(0) = %.3f\n",
            rt$R_inst[1], TARGET_R0, rt$R_case[1]))

cross1 <- function(time, R) { i <- which(R < 1)[1]
  if (is.na(i) || i == 1) return(NA_real_)
  approx(c(R[i - 1], R[i]), c(time[i - 1], time[i]), xout = 1)$y }
cat(sprintf("R=1 crossing:  R_inst at t = %.1f d,  R_case at t = %.1f d\n",
            cross1(rt$time, rt$R_inst), cross1(rt$time, rt$R_case)))


## =============================================================================
## PART B -- GI-convolution cross-check (base R)
## =============================================================================
## genPop generation-time pmf w(tau) on daily bins, from Gamma(shape, rate).
gi_pmf <- function(shape, rate, max_day) {
  edges <- 0:max_day
  p <- diff(pgamma(edges, shape = shape, rate = rate))
  p / sum(p)
}
w <- gi_pmf(args$Tg_shape_genPop, args$Tg_rate_genPop,
            max_day = ceiling(qgamma(0.999, args$Tg_shape_genPop, args$Tg_rate_genPop)))
## R_case_pred(t) = sum_tau w(tau) * R_inst(t + tau), with R_inst clamped past grid.
Rinst_at <- function(tq) {
  approx(rt$time, rt$R_inst, xout = tq, rule = 2)$y
}
taus <- seq_along(w) - 0.5
Rcase_pred <- vapply(rt$time, function(t0) sum(w * Rinst_at(t0 + taus)), numeric(1))
resid <- rt$R_case - Rcase_pred
cat(sprintf("\nPART B  GI-convolution residual (R_case - conv(R_inst)):\n"))
cat(sprintf("  mean = %.4f, max|.| = %.4f  (diagnostic only; not asserted)\n",
            mean(resid), max(abs(resid))))


## =============================================================================
## FIGURE -- example R(t) curves
## =============================================================================
fig_path <- file.path(root, OUT_DIR, sprintf("Rt_example_%s.png", SCENARIO_ID))
grDevices::png(fig_path, width = 1400, height = 950, res = 170)
op <- par(mar = c(4.2, 4.2, 2.4, 1))
ylim <- range(0.8, rt$R_inst, rt$R_case)
plot(rt$time, rt$R_inst, type = "l", lwd = 2.4, col = "#0072B2", ylim = ylim,
     xlab = "Day since outbreak start", ylab = "Reproduction number",
     main = sprintf("Analytical single-type R(t) -- %s", SCENARIO_ID))
lines(rt$time, rt$R_case, lwd = 2.4, col = "#D55E00")
lines(rt$time, Rcase_pred, lwd = 1.2, lty = 3, col = "#D55E00")
abline(h = 1, col = "grey50", lty = 2)
legend("topright", bty = "n", lwd = c(2.4, 2.4, 1.2), lty = c(1, 1, 3),
       col = c("#0072B2", "#D55E00", "#D55E00"),
       legend = c("R_inst (instantaneous)", "R_case (cohort)",
                  "conv(R_inst) [GI check]"))
par(op); grDevices::dev.off()
cat(sprintf("\nWrote figure: %s\n", fig_path))


## =============================================================================
## PART A -- simulator overlay (needs fiber)
## =============================================================================
if (!have_fiber) {
  cat("\nPART A skipped: fiber not installed.\n")
} else {
  cat(sprintf("\nPART A: running %d fiber replicates (OBV off)...\n", N_REPS))
  args_sim <- args
  args_sim$obv_pep_enabled <- FALSE
  ## accumulate empirical R by weekly bin across replicates.
  emp <- list()
  for (r in seq_len(N_REPS)) {
    a <- args_sim; a$seed <- 1000L + r
    out <- do.call(fiber::branching_process_main, a)
    rn  <- fiber::compute_reproduction_number(out$tdf, type = "both",
                                              bin_width = BIN_WIDTH)
    emp[[r]] <- rn
  }
  ## bind and average per time bin. compute_reproduction_number is expected to
  ## return a per-bin data.frame with a time column plus R_instantaneous / R_case;
  ## detect the time column robustly.
  all_rn  <- do.call(rbind, emp)
  tcol    <- intersect(c("time", "bin", "day", "t", "time_bin"), names(all_rn))[1]
  if (is.na(tcol)) stop("Could not find a time column in compute_reproduction_number() output: ",
                        paste(names(all_rn), collapse = ", "))
  agg <- aggregate(cbind(R_instantaneous, R_case) ~ get(tcol), data = all_rn,
                   FUN = function(z) mean(z[is.finite(z)]))
  names(agg)[1] <- "time"

  ## analytical curves binned to the same weekly grid (bin centre lookup).
  Rinst_b <- approx(rt$time, rt$R_inst, xout = agg$time, rule = 2)$y
  Rcase_b <- approx(rt$time, rt$R_case, xout = agg$time, rule = 2)$y

  ov_path <- file.path(root, OUT_DIR, sprintf("Rt_overlay_%s.png", SCENARIO_ID))
  grDevices::png(ov_path, width = 1500, height = 1000, res = 170)
  op <- par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.6, 1))
  plot(agg$time, agg$R_instantaneous, pch = 16, col = "#0072B255",
       xlab = "Day", ylab = "R_instantaneous", main = "Instantaneous")
  lines(rt$time, rt$R_inst, lwd = 2.4, col = "#0072B2"); abline(h = 1, lty = 2, col = "grey")
  plot(agg$time, agg$R_case, pch = 16, col = "#D55E0055",
       xlab = "Day", ylab = "R_case", main = "Cohort / case")
  lines(rt$time, rt$R_case, lwd = 2.4, col = "#D55E00"); abline(h = 1, lty = 2, col = "grey")
  par(op); grDevices::dev.off()
  cat(sprintf("Wrote overlay: %s\n", ov_path))
}

cat("\nDone.\n")
