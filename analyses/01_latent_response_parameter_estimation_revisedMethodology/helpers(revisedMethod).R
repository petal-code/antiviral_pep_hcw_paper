# ============================================================================
# helpers(revisedMethod).R
# ----------------------------------------------------------------------------
# Small shared utilities used across the REVISED-METHODOLOGY latent-response
# pipeline. Keeping them here (rather than copy-pasting into every script) means
# there is exactly ONE definition of each, so they cannot silently drift apart.
#
#   source("helpers(revisedMethod).R")
#
# This is the sibling of helpers.R in the original-methodology folder. The two
# are deliberately kept as close as possible: the generic glue (paths, parameter
# vocabulary, interpolation, the SDB rolling mean) is IDENTICAL. The revised
# methodology adds exactly two pieces of machinery, defined at the bottom:
#   * q_norm()        - the finite-window normalised logistic (Q runs 0 -> 1),
#                       the R twin of the function inside the Stan model; and
#   * lock_endpoints()- the endpoint-locking rule that IS the revised
#                       methodology (see the long comment on that function).
#
# WHAT "REVISED METHODOLOGY" MEANS (one paragraph)
#   In the ORIGINAL methodology each parameter's two magnitude endpoints are
#   ESTIMATED (Model A estimates the curve shape AND the endpoints; Model B fixes
#   the shared Q and estimates the endpoints). In the REVISED methodology the
#   endpoints are instead LOCKED to early/late literature-window extrema and only
#   the curve SHAPE is estimated (West Africa, "Model C"); for DRC the shape is
#   not estimated at all - the empirical conflict Q is mapped directly onto the
#   locked endpoints. A second consequence (handled in 02/03) is that the IPC/PPE
#   index for the DRC conflict scenarios is q-scaled between its locked endpoints
#   rather than read from a fitted latent_IPC curve.
# ============================================================================

# ---- Project paths ---------------------------------------------------------
# Every script resolves paths from the repository root with here::here(), which
# locates obv_hcw_paper.Rproj / the .git directory. This means the scripts run
# correctly no matter which working directory they are launched from. Raw inputs
# and processed outputs live at the TOP LEVEL of the repository (shared across
# analyses, matching the rest of obv_hcw_paper); the Stan model lives in this
# analysis folder alongside the scripts.
#
# NOTE the only difference from the original-methodology helpers.R: ANALYSIS_DIR
# (and therefore DIR_STAN) point at the *_revisedMethodology folder.
ANALYSIS_DIR  <- here::here("analyses", "01_latent_response_parameter_estimation_revisedMethodology")
DIR_RAW       <- here::here("data-raw")
DIR_PROCESSED <- here::here("data-processed")
DIR_STAN      <- file.path(ANALYSIS_DIR, "stan-models")

# ---- The shared model "vocabulary" -----------------------------------------

# The six latent response parameters, in a fixed canonical order. Every script
# uses this ordering so that parameter id 1..6 means the same thing everywhere.
# These are the INTERNAL names (the published column names - prob_hosp, prop_etu,
# ipc_helper, ... - are applied only at the very end, in 03's assemble_scenario,
# exactly as in the original-methodology pipeline).
PARAM_LEVELS <- c(
  "delay_hosp",            # mean delay (days) from symptom onset to hospitalisation
  "p_hosp",                # probability an infected person is hospitalised
  "p_ETU",                 # proportion of hospitalised cases managed in an ETU/ETC
  "latent_IPC",            # latent infection-prevention-and-control / PPE index
  "p_unsafe_funeral_comm", # probability of an unsafe funeral after a community death
  "p_unsafe_funeral_hosp"  # probability of an unsafe funeral after a hospital death
)

# Human-readable panel titles (used when plotting), keyed by parameter.
PANEL_LOOKUP <- c(
  delay_hosp            = "Delay to hospitalisation",
  p_hosp                = "Probability of hospitalisation",
  p_ETU                 = "Proportion in ETU / ETC",
  latent_IPC            = "Latent IPC / PPE index",
  p_unsafe_funeral_comm = "Unsafe funeral after community death",
  p_unsafe_funeral_hosp = "Unsafe funeral after hospital death"
)

# The final scenario matrices are all reported on a common 0..730 day horizon
# (731 daily rows), with tau = relative_day / HORIZON_DAYS.
HORIZON_DAYS <- 730L

# ---- Tiny numeric helpers ---------------------------------------------------

# Clamp a value (or vector) into the closed unit interval [0, 1].
clip01 <- function(x) pmin(1, pmax(0, x))

# Min-max rescale to [0, 1]. Used only for the diagnostic q_value column.
rescale_01 <- function(x) {
  r <- range(x, na.rm = TRUE)
  if (!all(is.finite(r))) stop("rescale_01(): non-finite range.")
  if (diff(r) <= 0)       stop("rescale_01(): zero-width range.")
  (x - r[1]) / diff(r)
}

# ---- Locating raw input files ----------------------------------------------

# Resolve a raw input workbook by name inside the top-level data-raw/. Accepts
# several candidate filenames (the SDB workbook has been distributed under a
# couple of slightly different names) and returns the first that exists.
resolve_input_file <- function(candidates, description = "input file",
                               data_raw_dir = DIR_RAW) {
  for (nm in candidates) {
    p <- file.path(data_raw_dir, nm)
    if (file.exists(p)) return(p)
  }
  stop(
    "Could not find ", description, " in '", data_raw_dir, "'. Looked for: ",
    paste(candidates, collapse = ", ")
  )
}

# ---- Interpolation ----------------------------------------------------------

# Build a linear interpolation function from (x, y) support points. Extrapolates
# by holding the nearest endpoint flat (rule = 2), which is what we want for the
# response curves (before the first / after the last support point the value
# simply holds). Duplicated x values are averaged.
make_interp <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  stats::approxfun(x[ok], y[ok], method = "linear", rule = 2, ties = mean)
}

# Centred rolling mean over an ordered vector (used to smooth the weekly SDB
# reconstruction). For an even window of width k (default 4) the window is
# necessarily slightly asymmetric: floor((k-1)/2) points before, the current
# point, and the remainder after. Positions with no finite values stay NA and
# are filled by the caller.
rolling_mean_centered <- function(y, k = 4) {
  n <- length(y)
  if (n == 0) return(numeric(0))
  if (k < 1) stop("rolling window k must be >= 1.")
  out <- rep(NA_real_, n)
  left_n  <- floor((k - 1) / 2)
  right_n <- (k - 1) - left_n
  for (i in seq_len(n)) {
    idx <- seq.int(max(1, i - left_n), min(n, i + right_n))
    yy  <- y[idx]
    ok  <- is.finite(yy)
    if (any(ok)) out[i] <- mean(yy[ok])
  }
  out
}

# ---- The normalised logistic response-quality curve Q(tau) -----------------

# The raw (un-normalised) logistic, and its finite-window normalisation so that
# Q(0) = 0 and Q(1) = 1 exactly over the observed window, regardless of (t50, k).
# This is the EXACT R twin of the q_raw()/q_norm() functions inside the Stan
# model (modelC_endpointConstrained_estimateShape). Having an R copy lets the
# deterministic DRC mapping (02) and any diagnostics reproduce the same Q the
# sampler uses, with no risk of the two definitions drifting apart.
q_raw_logistic <- function(tau, t50, k) stats::plogis(k * (tau - t50))
q_norm <- function(tau, t50, k) {
  q0 <- q_raw_logistic(0, t50, k)
  q1 <- q_raw_logistic(1, t50, k)
  (q_raw_logistic(tau, t50, k) - q0) / (q1 - q0)
}

# ---- THE endpoint-locking rule (this IS the revised methodology) -----------
#
# In the revised methodology a parameter's two magnitude endpoints are not
# estimated; they are LOCKED to literature-supported extrema, and only the curve
# shape (or, for DRC, nothing at all) is left to the model. lock_endpoints()
# implements that rule on a cleaned anchor table (the SAME anchors prepared in
# 00 and used by the original-methodology fits), so both methodologies start
# from identical evidence.
#
# THE RULE, per parameter:
#   * INCREASING parameter (direction "up"): start at the EARLIEST-window MINIMUM
#     observed value, end at the LATE-window MAXIMUM.
#   * DECREASING parameter (direction "down"): start at the early-window MAXIMUM,
#     end at the late-window MINIMUM.
#   * TERMINAL-ZERO override (unsafe-funeral parameters only): if an explicit zero
#     is observed anywhere, the END endpoint is forced to zero and Q is made to
#     reach 1 at that zero-anchor day (the curve then plateaus at zero). This
#     stops a curve ending above zero just because the zero anchor falls before
#     the scenario end.
#   * FALLBACK: if a parameter has no anchor inside the early (or late) window,
#     the corresponding endpoint falls back to the curated workbook range
#     (lower_bound / upper_bound), mapped by direction. Such rows are flagged in
#     `start_source` / `end_source` so they can be audited.
#
# Inputs
#   anchors                  cleaned anchor table (one row per literature point):
#                            needs columns parameter, relative_day, value_used,
#                            direction ("up"/"down"), lower_bound, upper_bound,
#                            and (optionally) fit_role.
#   scenario_duration_days   the scenario horizon in days; ordinary endpoints are
#                            reached at this day (terminal-zero endpoints earlier).
#   early_window_day         anchors on/before this day form the "early" window.
#   late_start_day           anchors on/after this day form the "late" window.
#   late_start_day_by_param  optional named override of late_start_day for chosen
#                            parameters (e.g. latent_IPC uses a later evidence
#                            window than the other parameters).
#   zero_endpoint_params     parameters eligible for the terminal-zero override.
#   zero_tol                 |value| <= zero_tol counts as an observed zero.
#
# Returns one row per parameter with columns:
#   parameter, direction, theta_start, theta_end, endpoint_day_for_tau,
#   start_source, end_source
# where endpoint_day_for_tau is the day at which Q reaches 1 for that parameter
# (scenario end, or the terminal-zero day) and is used to set tau = day / that.
lock_endpoints <- function(anchors,
                           scenario_duration_days,
                           early_window_day        = 50,
                           late_start_day          = 325,
                           late_start_day_by_param = c(latent_IPC = 275),
                           zero_endpoint_params    = c("p_unsafe_funeral_comm",
                                                       "p_unsafe_funeral_hosp"),
                           zero_tol                = 1e-12) {

  # Pick the start/end extreme value within a window, given the direction. NULL
  # if the window is empty.
  choose_extreme <- function(window, direction, which_end) {
    if (nrow(window) == 0) return(NULL)
    if (direction == "up") {
      idx <- if (which_end == "start") which.min(window$value_used) else which.max(window$value_used)
    } else {                         # "down"
      idx <- if (which_end == "start") which.max(window$value_used) else which.min(window$value_used)
    }
    window[idx[1], , drop = FALSE]
  }

  anchors %>%
    dplyr::group_by(parameter) %>%
    dplyr::group_modify(function(g, key) {
      p   <- key$parameter[[1]]
      dir <- g$direction[[1]]                       # "up" or "down"
      lo_bound <- g$lower_bound[[1]]                # workbook range (fallbacks)
      hi_bound <- g$upper_bound[[1]]

      # Anchors eligible to SET an endpoint: finite, and not reference/plot only.
      eligible <- g %>%
        dplyr::filter(is.finite(value_used), is.finite(relative_day))
      if ("fit_role" %in% names(eligible)) {
        eligible <- eligible %>%
          dplyr::filter(is.na(fit_role) |
                          !tolower(trimws(fit_role)) %in% c("reference_only", "plot_only"))
      }

      # Early window for the start endpoint.
      early <- eligible %>% dplyr::filter(relative_day <= early_window_day)
      start_row <- choose_extreme(early, dir, "start")

      # Late window for the end endpoint (parameter-specific cutoff if supplied).
      late_cut <- if (p %in% names(late_start_day_by_param)) {
        unname(late_start_day_by_param[[p]])
      } else {
        late_start_day
      }
      late <- eligible %>% dplyr::filter(relative_day >= late_cut)

      # Terminal-zero override for unsafe-funeral parameters.
      zero_rows <- if (p %in% zero_endpoint_params) {
        eligible %>% dplyr::filter(abs(value_used) <= zero_tol)
      } else {
        eligible[0, , drop = FALSE]
      }

      if (nrow(zero_rows) > 0) {
        end_day  <- max(zero_rows$relative_day, na.rm = TRUE)   # latest observed zero
        end_val  <- 0
        end_src  <- "terminal_zero_anchor"
        end_day_for_tau <- max(end_day, 1)
      } else {
        end_row <- choose_extreme(late, dir, "end")
        if (!is.null(end_row)) {
          end_val <- end_row$value_used[[1]]
          end_src <- if (p %in% names(late_start_day_by_param)) {
            paste0("late_window_extreme_anchor_fromDay", late_cut)
          } else {
            "late_window_extreme_anchor"
          }
        } else {
          end_val <- if (dir == "up") hi_bound else lo_bound  # workbook fallback
          end_src <- "summary_range_fallback_no_late_anchor"
        }
        end_day_for_tau <- as.numeric(scenario_duration_days)
      }

      # Start endpoint (or workbook fallback).
      if (!is.null(start_row)) {
        start_val <- start_row$value_used[[1]]
        start_src <- "early_window_extreme_anchor"
      } else {
        start_val <- if (dir == "up") lo_bound else hi_bound   # workbook fallback
        start_src <- "summary_range_fallback_no_early_anchor"
      }

      tibble::tibble(
        direction            = dir,
        theta_start          = start_val,
        theta_end            = end_val,
        endpoint_day_for_tau = end_day_for_tau,
        start_source         = start_src,
        end_source           = end_src
      )
    }) %>%
    dplyr::ungroup()
}
