# ============================================================================
# 06_DRC_ConflictSDB_CommunityBlend.R   (EXPLORATORY -- tune before adopting)
# ----------------------------------------------------------------------------
# Reshape the CONFLICT-period SDB reversion using the community-death proportion
# as a data-driven template, instead of the flat "set success -> 0 over 200-300"
# box used by the conflict++ scenario.
#
# Motivation (CW):
#   * the ++ box collapse (success -> 0 across days 200-300) is too sharp -- a
#     step down and back up;
#   * the RAW SDB curve barely reverts, because it is built from deaths that
#     generated alerts, and during the conflict many deaths never generated one
#     (so the denominator misses them and success looks artificially fine).
#   The community-death proportion DOES show a real peak in that window, so we
#   use its SHAPE to drive a smooth, deep, data-shaped reversion, blended back to
#   the raw SDB curve away from the conflict.
#
# Sign convention (matches the pipeline): success_smoothed = safe-burial success;
# unsafe = 1 - success. "Reversion" deepens => success DOWN. At the community-
# death PEAK the response is worst => success -> (1 - DEPTH) (0 at full depth).
#
# Method (daily timeline):
#   c(t)   community-death proportion curve (lightly smoothed for shape)
#   c~(t)  registered so its conflict-era peak lands on REGISTER_TARGET (onset/middle/none)
#   r(t)   = c~(t) / max(c~)                       reversion intensity in [0,1]
#   s_cd   = 1 - DEPTH * r(t)                      community-implied success
#   w(t)   = 1 on [200,300], cosine taper to 0 over TAPER days outside
#   s(t)   = (1 - w) * s_raw + w * s_cd            blended success
#   then q_value = s / max(s), unsafe = 1 - s      (as finalise_q_series)
#
# Output: a drop-in conflict_cdblend qseries (same columns as the other
#   scenarios) saved for inspection; swap it into 02 / 04 / 05 via qseries_for().
#
# Prereqs: 00_DataPreparation_and_Cleaning.R, 00b_DataPreparation_CommunityDeaths.R
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(ggplot2)
})
source(here::here("analyses", "01_latent_response_parameter_estimation", "helpers.R"))
set.seed(123)

# ---- knobs (EDIT) ----------------------------------------------------------
CONFLICT_WINDOW  <- c(200L, 300L)                # matches PLUSPLUS_WINDOW_DAY in 00 (SDB clock)
CONFLICT_MID     <- as.integer(mean(CONFLICT_WINDOW))  # 250
CD_TO_SDB_OFFSET <- 5L      # SDB day0 (2018-08-06, first burial) is 5 days after cd day0 (2018-08-01)
CONFLICT_SEARCH  <- c(100L, 350L)  # search this SDB-day range for the CONFLICT-ERA cd peak
REGISTER_TARGET  <- "onset" # pin the cd peak to: "onset" (window start), "middle", or "none"
DEPTH            <- 1.0     # reversion depth at the peak: success -> 1 - DEPTH (1 => 0)
TAPER_DAYS       <- 75L     # cosine taper half-width outside the window (community -> raw)
CD_SMOOTH_K      <- 3L      # light rolling mean on the community series (shape, not noise)

# ---- inputs ----------------------------------------------------------------
drc_prep <- readRDS(file.path(DIR_PROCESSED, "DRC_QCurve", "DRC_QCurve_PreppedData.rds"))
cd_prep  <- readRDS(file.path(DIR_PROCESSED, "DRC_QCurve", "DRC_CommunityDeaths_Prepped.rds"))
raw      <- drc_prep$conflict_qseries            # the curve we reshape (conflict scenario)
cd_obs   <- cd_prep$obs

horizon <- max(raw$relative_day)
day     <- 0:horizon

# ---- 1. Raw SDB success on the daily grid ----------------------------------
s_raw <- clip01(make_interp(raw$relative_day, raw$success_smoothed)(day))

# ---- 2. Community-death shape on the SDB clock, registered ------------------
# Put the community series on the SDB clock (subtract the day0 offset), then pin
# its CONFLICT-ERA peak to the registration target. The peak is found on the RAW
# proportion INSIDE CONFLICT_SEARCH, so the pre-response Aug-2018 community-death
# high (also ~70%, but not a conflict signal) cannot hijack the registration.
cd_pts <- cd_obs %>% arrange(relative_day) %>%
  mutate(relday_sdb = relative_day - CD_TO_SDB_OFFSET, p = n_comm / N_deaths)
cd_pts$p_s <- if (CD_SMOOTH_K > 1L) {
  sm <- rolling_mean_centered(cd_pts$p, k = CD_SMOOTH_K); sm[!is.finite(sm)] <- cd_pts$p[!is.finite(sm)]; sm
} else cd_pts$p

in_search <- cd_pts$relday_sdb >= CONFLICT_SEARCH[1] & cd_pts$relday_sdb <= CONFLICT_SEARCH[2]
if (!any(in_search)) stop("No community-death points inside CONFLICT_SEARCH.")
t_peak     <- cd_pts$relday_sdb[in_search][which.max(cd_pts$p[in_search])]  # raw argmax in-window
c_peak     <- max(cd_pts$p_s[in_search])         # conflict-era (smoothed) peak -> r = 1 here
target_day <- switch(REGISTER_TARGET, none = t_peak, onset = CONFLICT_WINDOW[1],
                     middle = CONFLICT_MID, stop("REGISTER_TARGET must be none/onset/middle"))
delta      <- as.integer(target_day - t_peak)
message(sprintf("Conflict-era cd peak at SDB day %d; target '%s' day %d; shift %+d.",
                t_peak, REGISTER_TARGET, target_day, delta))

# shift the support by delta so the peak lands on the target, evaluate on the grid
c_reg <- clip01(make_interp(cd_pts$relday_sdb + delta, cd_pts$p_s)(day))
r     <- clip01(c_reg / c_peak)                  # reversion intensity in [0,1]
s_cd  <- clip01(1 - DEPTH * r)                    # community-implied success

# ---- 3. Conflict-window blend weight (1 in core, cosine taper outside) ------
conflict_weight <- function(day, win, taper) {
  lo <- win[1]; hi <- win[2]; w <- numeric(length(day))
  w[day >= lo & day <= hi] <- 1
  lt <- day >= (lo - taper) & day < lo            # left shoulder 0 -> 1
  w[lt] <- 0.5 * (1 - cos(pi * (day[lt] - (lo - taper)) / taper))
  rt <- day > hi & day <= (hi + taper)            # right shoulder 1 -> 0
  w[rt] <- 0.5 * (1 + cos(pi * (day[rt] - hi) / taper))
  w
}
w <- conflict_weight(day, CONFLICT_WINDOW, TAPER_DAYS)

# ---- 4. Blend and recompute the qseries ------------------------------------
s_blend <- clip01((1 - w) * s_raw + w * s_cd)

# Sample the blended success back onto the original support days -> drop-in qseries.
s_at_support <- clip01(make_interp(day, s_blend)(raw$relative_day))
scale_max    <- max(s_at_support, na.rm = TRUE)
conflict_cdblend_qseries <- raw %>%
  mutate(success_smoothed          = s_at_support,
         q_value                   = clip01(s_at_support / scale_max),
         unsafe_funeral_comm_proxy = clip01(1 - s_at_support))

saveRDS(list(qseries = conflict_cdblend_qseries,
             daily = tibble(day, s_raw, s_cd, s_blend, w, r, c_reg),
             knobs = list(window = CONFLICT_WINDOW, mid = CONFLICT_MID,
                          register = REGISTER_TARGET, offset = CD_TO_SDB_OFFSET,
                          search = CONFLICT_SEARCH, depth = DEPTH, taper = TAPER_DAYS,
                          cd_smooth_k = CD_SMOOTH_K, cd_peak_day = t_peak, c_peak = c_peak,
                          delta = delta)),
        file.path(DIR_PROCESSED, "DRC_QCurve", "DRC_ConflictSDB_CommunityBlend.rds"))
message("Saved DRC_ConflictSDB_CommunityBlend.rds")

# ---- 5. Diagnostic plot ----------------------------------------------------
# Success curves on one axis: raw SDB, community-implied (1 - DEPTH*r), and the
# blend; the blend weight w; the conflict window shaded; and the community-death
# data points (registered onto the grid, mapped to their implied success
# 1 - DEPTH*p/c_peak) sized by denominator.
curves <- tibble(day, s_raw, s_cd, s_blend, w) %>%
  pivot_longer(c(s_raw, s_cd, s_blend, w), names_to = "series", values_to = "value") %>%
  mutate(series = recode(series,
                         s_raw   = "raw SDB success",
                         s_cd    = "community-implied (1 - DEPTH*r)",
                         s_blend = "blended success",
                         w       = "blend weight w"))

cd_implied <- cd_pts %>%
  transmute(day = relday_sdb + delta,
            value = clip01(1 - DEPTH * (p_s / c_peak)),
            N_deaths = N_deaths)

print(
  ggplot(curves, aes(day, value, colour = series, linetype = series, linewidth = series)) +
    annotate("rect", xmin = CONFLICT_WINDOW[1], xmax = CONFLICT_WINDOW[2],
             ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.5) +
    geom_vline(xintercept = CONFLICT_MID, colour = "grey60", linetype = "dotted") +
    geom_line() +
    geom_point(data = cd_implied, aes(day, value, size = N_deaths),
               inherit.aes = FALSE, shape = 21, fill = "#ff7f0e",
               colour = "white", stroke = 0.3, alpha = 0.85) +
    scale_size_area(max_size = 6, name = "deaths (N)") +
    scale_colour_manual(values = c("raw SDB success" = "grey55",
                                   "community-implied (1 - DEPTH*r)" = "#ff7f0e",
                                   "blended success" = "#1f77b4",
                                   "blend weight w" = "black")) +
    scale_linetype_manual(values = c("raw SDB success" = "solid",
                                     "community-implied (1 - DEPTH*r)" = "dashed",
                                     "blended success" = "solid",
                                     "blend weight w" = "dotted")) +
    scale_linewidth_manual(values = c("raw SDB success" = 0.8,
                                      "community-implied (1 - DEPTH*r)" = 0.8,
                                      "blended success" = 1.3,
                                      "blend weight w" = 0.6)) +
    labs(title = "Conflict SDB reversion reshaped by the community-death proportion",
         subtitle = sprintf("window %d-%d; register=%s (cd peak SDB day %d, shift %+d); DEPTH=%.2f, taper=%d",
                            CONFLICT_WINDOW[1], CONFLICT_WINDOW[2],
                            REGISTER_TARGET, t_peak, delta, DEPTH, TAPER_DAYS),
         x = "Relative outbreak day", y = "Safe-burial success",
         colour = NULL, linetype = NULL, linewidth = NULL) +
    coord_cartesian(ylim = c(0, 1)) +
    theme_bw(base_size = 11) + theme(legend.position = "top")
)
