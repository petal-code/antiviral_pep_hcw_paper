# ============================================================================
# ODV NHP delay-to-initiation efficacy -- plotting step
# ============================================================================
#
# Purpose
# -------
# Reads the saved fit object(s) from the ODV delay-efficacy analysis and writes
# two figures:
#   1. odv_delay_efficacy_stan_curve     -- the full Bayesian (Stan) posterior
#      median delay-efficacy curve with its 95% credible ribbon, and the
#      empirical hazard-scale points (with bootstrap intervals) overlaid.
#   2. odv_delay_efficacy_method_comparison -- the Bayesian curve overlaid on the
#      profiled-likelihood curve from 01_fit_odv_delay_efficacy.R, to compare the
#      two uncertainty treatments. Produced only if the script-01 fit is present.
#
# This is a display/plotting step: it reads the .rds written by the fitting
# scripts and writes finished figures. It does not refit anything.
#
# Inputs
# ------
#   data-processed/odv_nhp_delay/odv_ebov_rhesus_delay_efficacy_fit_stan.rds  (from 02; required)
#   data-processed/odv_nhp_delay/odv_ebov_rhesus_delay_efficacy_fit.rds       (from 01; optional)
#
# Output
# ------
#   outputs/odv_nhp_delay_efficacy/odv_delay_efficacy_stan_curve.{png,pdf}
#   outputs/odv_nhp_delay_efficacy/odv_delay_efficacy_method_comparison.{png,pdf}
#
# Run from anywhere in the repository (paths use here::here()):
#   Rscript analyses/odv_nhp_delay_efficacy/03_plot_odv_delay_efficacy.R
#
# ============================================================================


# ------------------------------------------------------------
# 0) Load packages
# ------------------------------------------------------------
suppressPackageStartupMessages({
  library(ggplot2)
})


# ------------------------------------------------------------
# 1) Paths
# ------------------------------------------------------------
proc_dir      <- here::here("data-processed", "odv_nhp_delay")
stan_path     <- file.path(proc_dir, "odv_ebov_rhesus_delay_efficacy_fit_stan.rds")
profiled_path <- file.path(proc_dir, "odv_ebov_rhesus_delay_efficacy_fit.rds")

out_dir <- here::here("outputs", "odv_nhp_delay_efficacy")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(stan_path)) {
  stop(
    "Could not find the Stan fit: ", stan_path, "\n",
    "Run analyses/odv_nhp_delay_efficacy/02_fit_odv_delay_efficacy_stan.R first."
  )
}


# ------------------------------------------------------------
# 2) Read fit, set shared plotting elements
# ------------------------------------------------------------
stan_fit   <- readRDS(stan_path)
curve_stan <- stan_fit$fitted_curve      # dpc, efficacy, efficacy_lo, efficacy_hi
emp        <- stan_fit$empirical_points  # dpc, efficacy_hazard_scale, efficacy_lo, efficacy_hi, ...
fs         <- stan_fit$fit_summary

# Largest observed initiation day: efficacy beyond this is extrapolation driven
# by the curve's parametric form and the d_zero constraint, not by data.
max_obs_dpc <- max(emp$dpc)
dpc_max     <- max(curve_stan$dpc)

col_stan <- "#1f77b4"  # Bayesian (Stan)
col_prof <- "#d62728"  # profiled + Laplace (script 01)
col_emp  <- "grey20"   # empirical points

pct_lab <- function(x) paste0(x * 100, "%")

theme_odv <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold"),
    plot.caption     = element_text(size = 8, colour = "grey35"),
    legend.position  = "top"
  )

# Save a figure as both PNG (raster, for slides) and PDF (vector, for the paper).
save_fig <- function(base, plot, width = 7, height = 5) {
  ggsave(file.path(out_dir, paste0(base, ".png")), plot,
         width = width, height = height, dpi = 400, units = "in")
  ggsave(file.path(out_dir, paste0(base, ".pdf")), plot,
         width = width, height = height, units = "in")
}


# ------------------------------------------------------------
# 3) Figure 1: Bayesian posterior curve + empirical points
# ------------------------------------------------------------
cap1 <- sprintf(
  paste0("Posterior median and 95%% credible interval. ",
         "E0 = %.2f [%.2f, %.2f]; d50 = %.1f [%.1f, %.1f]; k %s = %.2f. ",
         "Empirical points: 1 - KM hazard ratio at day %d (bars = bootstrap 95%%)."),
  fs$E0_median, fs$E0_lo, fs$E0_hi,
  fs$d50_median, fs$d50_lo, fs$d50_hi,
  if (isTRUE(fs$k_estimated)) "(fitted)" else "(fixed)", fs$k_median,
  stan_fit$metadata$settings$t_censor_plot
)

p1 <- ggplot() +
  # Shade the extrapolation region beyond the last observed treatment arm.
  annotate("rect", xmin = max_obs_dpc, xmax = dpc_max, ymin = 0, ymax = 1,
           fill = "grey85", alpha = 0.35) +
  geom_ribbon(data = curve_stan, aes(dpc, ymin = efficacy_lo, ymax = efficacy_hi),
              fill = col_stan, alpha = 0.20) +
  geom_line(data = curve_stan, aes(dpc, efficacy), colour = col_stan, linewidth = 1) +
  geom_errorbar(data = emp, aes(x = dpc, ymin = efficacy_lo, ymax = efficacy_hi),
                width = 0.15, colour = col_emp, alpha = 0.8, na.rm = TRUE) +
  geom_point(data = emp, aes(x = dpc, y = efficacy_hazard_scale),
             colour = col_emp, size = 2.4) +
  geom_vline(xintercept = max_obs_dpc, linetype = "dashed", colour = "grey50") +
  scale_y_continuous(limits = c(0, 1), labels = pct_lab,
                     expand = expansion(mult = c(0, 0.02))) +
  scale_x_continuous(breaks = seq(0, dpc_max, by = 3)) +
  labs(
    title    = "ODV delay-to-initiation efficacy (full Bayesian fit)",
    subtitle = "Hazard-scale efficacy vs treatment initiation day; grey band = extrapolation beyond observed arms",
    x        = "ODV initiation day (days post-challenge)",
    y        = "Efficacy (hazard scale)",
    caption  = cap1
  ) +
  theme_odv

save_fig("odv_delay_efficacy_stan_curve", p1)
message("Saved: ", file.path(out_dir, "odv_delay_efficacy_stan_curve.{png,pdf}"))


# ------------------------------------------------------------
# 4) Figure 2: Bayesian vs profiled-likelihood overlay (optional)
# ------------------------------------------------------------
# Only drawn if the script-01 profiled fit is present. The two curves should
# track each other in location; the Bayesian ribbon is expected to be wider
# because it propagates baseline-hazard uncertainty and posterior skew that the
# profiled 2-parameter Laplace step does not.
if (file.exists(profiled_path)) {

  prof_fit   <- readRDS(profiled_path)
  curve_prof <- prof_fit$fitted_curve

  keep <- c("dpc", "efficacy", "efficacy_lo", "efficacy_hi")
  curve_both <- rbind(
    transform(curve_stan[keep], method = "Bayesian (Stan)"),
    transform(curve_prof[keep], method = "Profiled + Laplace (01)")
  )
  curve_both$method <- factor(
    curve_both$method,
    levels = c("Bayesian (Stan)", "Profiled + Laplace (01)")
  )
  method_cols <- c("Bayesian (Stan)" = col_stan, "Profiled + Laplace (01)" = col_prof)

  p2 <- ggplot() +
    geom_ribbon(data = curve_both,
                aes(dpc, ymin = efficacy_lo, ymax = efficacy_hi, fill = method),
                alpha = 0.18, colour = NA, na.rm = TRUE) +
    geom_line(data = curve_both, aes(dpc, efficacy, colour = method), linewidth = 1) +
    geom_errorbar(data = emp, aes(x = dpc, ymin = efficacy_lo, ymax = efficacy_hi),
                  width = 0.15, colour = col_emp, alpha = 0.8, na.rm = TRUE) +
    geom_point(data = emp, aes(x = dpc, y = efficacy_hazard_scale),
               colour = col_emp, size = 2.4) +
    scale_colour_manual(values = method_cols, name = NULL) +
    scale_fill_manual(values = method_cols, name = NULL) +
    scale_y_continuous(limits = c(0, 1), labels = pct_lab,
                       expand = expansion(mult = c(0, 0.02))) +
    scale_x_continuous(breaks = seq(0, max(curve_both$dpc), by = 3)) +
    labs(
      title    = "ODV delay-efficacy curve: Bayesian vs profiled fit",
      subtitle = "Median/point estimate with 95% credible/approximate intervals; points = empirical (KM hazard ratio)",
      x        = "ODV initiation day (days post-challenge)",
      y        = "Efficacy (hazard scale)"
    ) +
    theme_odv

  save_fig("odv_delay_efficacy_method_comparison", p2)
  message("Saved: ", file.path(out_dir, "odv_delay_efficacy_method_comparison.{png,pdf}"))

} else {
  message(
    "Profiled fit not found (", profiled_path, "); skipping the method-comparison overlay.\n",
    "Run analyses/odv_nhp_delay_efficacy/01_fit_odv_delay_efficacy.R to enable it."
  )
}

message("Done. Figures written to: ", out_dir)
