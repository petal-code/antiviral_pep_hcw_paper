# ============================================================================
# 07_DRC_ConflictSDB_CommunityBlend_Sweep.R   (EXPLORATORY -- knob sweep)
# ----------------------------------------------------------------------------
# Sweeps the two main knobs of the community-death conflict-blend (see
# 06_DRC_ConflictSDB_CommunityBlend.R) and plots the resulting blended SDB
# success curves so DEPTH x TAPER_DAYS can be compared at a glance.
#
#   DEPTH       reversion depth at the peak: success -> 1 - DEPTH (1 => full
#               collapse to 0, as in conflict++)
#   TAPER_DAYS  cosine taper half-width outside the [200,300] window over which
#               the blend hands back from the community curve to the raw SDB curve
#
# Everything else (community shape, registration, window, smoothing) is held at
# the 06 defaults. Reversion intensity r(t) = c~(t)/max(c~) is the SAME across the
# whole sweep; DEPTH only sets the gain on it and TAPER only sets the wings.
#
# Prereqs: 00_DataPreparation_and_Cleaning.R, 00b_DataPreparation_CommunityDeaths.R
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(ggplot2)
})
source(here::here("analyses", "01_latent_response_parameter_estimation", "helpers.R"))
set.seed(123)

# ---- swept grids + held-fixed knobs (EDIT) ---------------------------------
DEPTH_GRID      <- c(0.7, 0.8, 0.9, 1.0)
TAPER_GRID      <- c(10L, 25L, 40L, 50L, 100L)

CONFLICT_WINDOW  <- c(200L, 300L)
CONFLICT_MID     <- as.integer(mean(CONFLICT_WINDOW))   # 250
CD_TO_SDB_OFFSET <- 5L              # SDB day0 (2018-08-06) is 5 days after cd day0 (2018-08-01)
CONFLICT_SEARCH  <- c(100L, 350L)   # search this SDB-day range for the CONFLICT-ERA cd peak
REGISTER_TARGET  <- "onset"         # pin the cd peak to: "onset", "middle", or "none"
CD_SMOOTH_K      <- 3L

# ---- inputs ----------------------------------------------------------------
drc_prep <- readRDS(file.path(DIR_PROCESSED, "DRC_QCurve", "DRC_QCurve_PreppedData.rds"))
cd_prep  <- readRDS(file.path(DIR_PROCESSED, "DRC_QCurve", "DRC_CommunityDeaths_Prepped.rds"))
raw      <- drc_prep$conflict_qseries
cd_obs   <- cd_prep$obs

horizon <- max(raw$relative_day)
day     <- 0:horizon

# ---- fixed pieces (independent of DEPTH / TAPER) ---------------------------
s_raw <- clip01(make_interp(raw$relative_day, raw$success_smoothed)(day))

# Community series on the SDB clock; CONFLICT-ERA peak found on the RAW proportion
# inside CONFLICT_SEARCH (so the pre-response Aug-2018 high cannot hijack it), then
# registered to the target. r(t) = c~(t) / conflict-era peak; same for all DEPTH.
cd_pts <- cd_obs %>% arrange(relative_day) %>%
  mutate(relday_sdb = relative_day - CD_TO_SDB_OFFSET, p = n_comm / N_deaths)
cd_pts$p_s <- if (CD_SMOOTH_K > 1L) {
  sm <- rolling_mean_centered(cd_pts$p, k = CD_SMOOTH_K); sm[!is.finite(sm)] <- cd_pts$p[!is.finite(sm)]; sm
} else cd_pts$p

in_search  <- cd_pts$relday_sdb >= CONFLICT_SEARCH[1] & cd_pts$relday_sdb <= CONFLICT_SEARCH[2]
if (!any(in_search)) stop("No community-death points inside CONFLICT_SEARCH.")
t_peak     <- cd_pts$relday_sdb[in_search][which.max(cd_pts$p[in_search])]
c_peak     <- max(cd_pts$p_s[in_search])
target_day <- switch(REGISTER_TARGET, none = t_peak, onset = CONFLICT_WINDOW[1],
                     middle = CONFLICT_MID, stop("REGISTER_TARGET must be none/onset/middle"))
delta      <- as.integer(target_day - t_peak)
c_reg      <- clip01(make_interp(cd_pts$relday_sdb + delta, cd_pts$p_s)(day))
r          <- clip01(c_reg / c_peak)             # reversion intensity, same for all DEPTH
message(sprintf("Conflict-era cd peak at SDB day %d; target '%s' day %d; shift %+d.",
                t_peak, REGISTER_TARGET, target_day, delta))

conflict_weight <- function(day, win, taper) {
  lo <- win[1]; hi <- win[2]; w <- numeric(length(day))
  w[day >= lo & day <= hi] <- 1
  lt <- day >= (lo - taper) & day < lo
  w[lt] <- 0.5 * (1 - cos(pi * (day[lt] - (lo - taper)) / taper))
  rt <- day > hi & day <= (hi + taper)
  w[rt] <- 0.5 * (1 + cos(pi * (day[rt] - hi) / taper))
  w
}

depth_lab <- function(d) factor(sprintf("DEPTH = %.1f", d),
                                levels = sprintf("DEPTH = %.1f", DEPTH_GRID))

# ---- sweep DEPTH x TAPER ---------------------------------------------------
sweep <- expand.grid(depth = DEPTH_GRID, taper = TAPER_GRID)
blend_df <- bind_rows(lapply(seq_len(nrow(sweep)), function(i) {
  d <- sweep$depth[i]; tp <- sweep$taper[i]
  w <- conflict_weight(day, CONFLICT_WINDOW, tp)
  tibble(day,
         s_blend     = clip01((1 - w) * s_raw + w * (1 - d * r)),
         depth       = d,
         taper       = factor(tp, levels = TAPER_GRID),
         depth_label = depth_lab(d))
}))

saveRDS(list(sweep = blend_df, day = day, s_raw = s_raw, r = r,
             knobs = list(window = CONFLICT_WINDOW, mid = CONFLICT_MID,
                          register = REGISTER_TARGET, offset = CD_TO_SDB_OFFSET,
                          search = CONFLICT_SEARCH, cd_smooth_k = CD_SMOOTH_K,
                          cd_peak_day = t_peak, c_peak = c_peak, delta = delta,
                          depth_grid = DEPTH_GRID, taper_grid = TAPER_GRID)),
        file.path(DIR_PROCESSED, "DRC_QCurve", "DRC_ConflictSDB_CommunityBlend_Sweep.rds"))
message("Saved DRC_ConflictSDB_CommunityBlend_Sweep.rds")

# ---- reference layers (drawn in every facet) -------------------------------
raw_ref  <- tibble(day, s_raw)                            # no depth_label -> all panels
# community-implied success 1 - DEPTH*r (the w=1 / fully-community-driven limit).
s_cd_ref <- bind_rows(lapply(DEPTH_GRID, function(d)
  tibble(day, s_cd = clip01(1 - d * r), depth_label = depth_lab(d))))

# ---- plot: facet by DEPTH, colour by TAPER ---------------------------------
print(
  ggplot(blend_df, aes(day, s_blend, colour = taper, group = taper)) +
    annotate("rect", xmin = CONFLICT_WINDOW[1], xmax = CONFLICT_WINDOW[2],
             ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.5) +
    geom_vline(xintercept = CONFLICT_MID, colour = "grey60", linetype = "dotted") +
    geom_line(data = raw_ref, aes(day, s_raw), inherit.aes = FALSE,
              colour = "grey55", linewidth = 0.7) +
    geom_line(data = s_cd_ref, aes(day, s_cd), inherit.aes = FALSE,
              colour = "grey20", linetype = "dashed", linewidth = 0.5) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~ depth_label, ncol = 2) +
    scale_colour_viridis_d(name = "TAPER_DAYS", end = 0.92) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = "Conflict-SDB community blend: DEPTH x TAPER_DAYS sweep",
         subtitle = sprintf(paste("blended safe-burial success; grey = raw SDB,",
                                   "dashed = community-implied (w=1) limit; window %d-%d, register=%s (shift %+d)"),
                            CONFLICT_WINDOW[1], CONFLICT_WINDOW[2], REGISTER_TARGET, delta),
         x = "Relative outbreak day", y = "Safe-burial success") +
    theme_bw(base_size = 11) +
    theme(legend.position = "top", strip.text = element_text(face = "bold"))
)
