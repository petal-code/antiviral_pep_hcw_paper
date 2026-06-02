# abc_calibration_functions_common.R
# -----------------------------------------------------------------------------
# Generic ABC-SMC helpers shared by BOTH fitting approaches (HCW-risk and
# NPI-efficacy). Nothing here depends on which parameters are fitted, so this
# file is sourced by abc_calibration_functions_hcwRisk.R and
# abc_calibration_functions_npi.R (and by their run scripts and PSOCK workers).
#
# Contents:
#   abc_summarise()          : per-replicate summary (n_cases/deaths/hcw/duration)
#   run_one_abc_replicate()  : do.call(branching_process_main) + abc_summarise()
#   save_abc_config()        : serialise worker config; advertise via FIBER_ABC_CONFIG
#   make_abc_output_dir()    : fresh timestamped per-run output directory
#   with_abc_output_dir()    : run an expression with cwd set to that directory
#   abc_progress()           : print progress from output_step*/tolerance_step* files
#   abc_compare_steps()      : weighted summary across completed steps
#   reconstruct_abc_result() : rebuild an ABC_sequential()-style result from disk
#
# The disk-inspection trio take `param_names` (default = the HCW-risk triple);
# callers with a different parameterisation pass their own, e.g.
# c("R0", "prop_funeral", "npi_scaler"). Requires fiber to be loaded.
# -----------------------------------------------------------------------------

abc_summarise <- function(out) {
  tdf <- out$tdf
  tdf <- tdf[!is.na(tdf$time_infection_absolute), , drop = FALSE]

  if (nrow(tdf) == 0L) {
    return(c(n_cases = 0, n_deaths = 0, n_hcw_deaths = 0, duration = 0))
  }

  deaths <- !is.na(tdf$outcome) & tdf$outcome
  hcw    <- !is.na(tdf$class) & tdf$class == "HCW"

  n_cases      <- nrow(tdf)
  n_deaths     <- sum(deaths)
  n_hcw_deaths <- sum(deaths & hcw)

  if (n_deaths == 0L) {
    duration <- 0
  } else {
    # Duration is measured from the first death to the last event, so it tracks
    # the part of the outbreak that would be observable in surveillance data.
    t_first_death <- min(tdf$time_outcome_absolute[deaths], na.rm = TRUE)
    t_last_event  <- max(tdf$time_outcome_absolute,         na.rm = TRUE)
    duration      <- max(t_last_event - t_first_death, 0)
  }

  c(n_cases      = n_cases,
    n_deaths     = n_deaths,
    n_hcw_deaths = n_hcw_deaths,
    duration     = duration)
}


run_one_abc_replicate <- function(args) {
  out <- do.call(branching_process_main, args)
  abc_summarise(out)
}


save_abc_config <- function(config, file = tempfile(fileext = ".rds")) {
  required <- c("setup_path", "functions_path", "r0_path",
                "scenario_csv", "scenario_id")
  missing_keys <- setdiff(required, names(config))
  if (length(missing_keys) > 0L) {
    stop("save_abc_config(): config is missing required key(s): ",
         paste(missing_keys, collapse = ", "), call. = FALSE)
  }
  saveRDS(config, file = file)
  Sys.setenv(FIBER_ABC_CONFIG = file)
  invisible(file)
}


make_abc_output_dir <- function(base_dir,
                                scenario_id,
                                label = NULL,
                                suffix = NULL,
                                subdir = "abc_outputs",
                                timestamp = TRUE) {
  if (missing(base_dir) || is.null(base_dir) || !nzchar(base_dir)) {
    stop("`base_dir` is required.", call. = FALSE)
  }
  if (missing(scenario_id) || is.null(scenario_id) || !nzchar(scenario_id)) {
    stop("`scenario_id` is required.", call. = FALSE)
  }

  parts <- scenario_id
  if (isTRUE(timestamp)) {
    parts <- paste(parts, format(Sys.time(), "%Y%m%d_%H%M%S"), sep = "_")
  }
  if (!is.null(label) && nzchar(label)) {
    parts <- paste(parts, label, sep = "_")
  }
  # `suffix` is appended last so any run-identifying tag (e.g. the
  # NBREPS_X_NBSIMUL_Y settings tag) lands at the very END of the directory
  # name, after the timestamp and label.
  if (!is.null(suffix) && nzchar(suffix)) {
    parts <- paste(parts, suffix, sep = "_")
  }

  out <- file.path(base_dir, subdir, parts)
  dedup_suffix <- 0L
  candidate <- out
  while (dir.exists(candidate)) {
    dedup_suffix <- dedup_suffix + 1L
    candidate <- paste0(out, "_", sprintf("%02d", dedup_suffix))
  }
  out <- candidate

  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  out
}


with_abc_output_dir <- function(output_dir, expr) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  old <- setwd(output_dir)
  on.exit(setwd(old), add = TRUE)
  force(expr)
}


abc_progress <- function(dir = getwd(),
                         tolerance_target = 1.0,
                         param_names = c("R0", "prop_funeral", "hcw_risk_scalar")) {

  step_num     <- function(f) as.integer(sub(".*_step([0-9]+)$", "\\1", f))
  sort_by_step <- function(f) f[order(step_num(f))]

  out_files <- sort_by_step(list.files(dir, pattern = "^output_step[0-9]+$",       full.names = TRUE))
  tol_files <- sort_by_step(list.files(dir, pattern = "^tolerance_step[0-9]+$",    full.names = TRUE))
  sim_files <- sort_by_step(list.files(dir, pattern = "^n_simul_tot_step[0-9]+$",  full.names = TRUE))

  if (length(out_files) == 0L) {
    cat("No steps completed yet.\n")
    return(invisible(NULL))
  }

  cat(sprintf("Steps completed: %d\n", length(out_files)))

  if (length(sim_files) > 0L) {
    cum_sims <- as.numeric(readLines(sim_files[length(sim_files)]))
    cat(sprintf("Total simulations so far: %s\n", format(cum_sims, big.mark = ",")))
  }

  if (length(tol_files) > 0L) {
    tols <- vapply(tol_files, function(f) as.numeric(readLines(f)), numeric(1))
    cat("\nTolerance trajectory:\n")
    for (i in seq_along(tols)) {
      cat(sprintf("  After step %d: %.3f\n", step_num(tol_files)[i], tols[i]))
    }
    cat(sprintf("\nTarget tolerance: %.3f\n", tolerance_target))

    current <- tols[length(tols)]
    if (current <= tolerance_target) {
      cat("Reached target -- algorithm should terminate after this step.\n")
    } else if (length(tols) >= 2L) {
      ratio <- tols[length(tols)] / tols[length(tols) - 1L]
      if (ratio < 1) {
        remaining <- ceiling(log(tolerance_target / current) / log(ratio))
        cat(sprintf(
          "Projected remaining steps: ~%d (assuming current shrinkage rate continues)\n",
          remaining
        ))
      } else {
        cat("Tolerance is not decreasing -- the algorithm may be stuck.\n")
      }
    }
  }

  last_out <- read.table(out_files[length(out_files)], header = FALSE)
  colnames(last_out) <- c("weight", param_names,
                          "takeoff", "n_deaths", "n_hcw_deaths", "duration")
  cat(sprintf("\nLatest particle cloud (step %d):\n",
              step_num(out_files)[length(out_files)]))
  for (nm in param_names) {
    cat(sprintf("  %-14s median = %.3f, 95%% CI = [%.3f, %.3f]\n",
                paste0(nm, ":"),
                median(last_out[[nm]]),
                quantile(last_out[[nm]], 0.025),
                quantile(last_out[[nm]], 0.975)))
  }

  invisible(last_out)
}


abc_compare_steps <- function(dir = getwd(),
                              param_names = c("R0", "prop_funeral", "hcw_risk_scalar")) {

  step_num     <- function(f) as.integer(sub(".*_step([0-9]+)$", "\\1", f))
  sort_by_step <- function(f) f[order(step_num(f))]

  out_files <- sort_by_step(list.files(dir, pattern = "^output_step[0-9]+$",      full.names = TRUE))
  tol_files <- sort_by_step(list.files(dir, pattern = "^tolerance_step[0-9]+$",   full.names = TRUE))
  sim_files <- sort_by_step(list.files(dir, pattern = "^n_simul_tot_step[0-9]+$", full.names = TRUE))

  if (length(out_files) == 0L) {
    cat("No steps completed yet.\n")
    return(invisible(NULL))
  }

  wq <- function(x, w, p) {
    ord <- order(x); xs <- x[ord]; ws <- w[ord] / sum(w)
    approx(cumsum(ws), xs, p, rule = 2, ties = "ordered")$y
  }

  tol_lookup <- if (length(tol_files) > 0L) {
    setNames(vapply(tol_files, function(f) as.numeric(readLines(f)), numeric(1)),
             as.character(step_num(tol_files)))
  } else c()
  cum_lookup <- if (length(sim_files) > 0L) {
    setNames(vapply(sim_files, function(f) as.numeric(readLines(f)), numeric(1)),
             as.character(step_num(sim_files)))
  } else c()

  rows <- list()
  prev_cum <- 0
  for (f in out_files) {
    s  <- step_num(f)
    df <- read.table(f, header = FALSE)
    colnames(df) <- c("weight", param_names,
                      "takeoff", "n_deaths", "n_hcw_deaths", "duration")
    w <- df$weight / sum(df$weight)

    cum_s     <- if (as.character(s) %in% names(cum_lookup)) cum_lookup[[as.character(s)]] else NA_real_
    sims_step <- if (!is.na(cum_s)) cum_s - prev_cum else NA_real_
    if (!is.na(cum_s)) prev_cum <- cum_s

    row <- data.frame(
      step           = s,
      tol            = if (as.character(s) %in% names(tol_lookup)) {
        round(tol_lookup[[as.character(s)]], 3)
      } else NA_real_,
      sims_this_step = sims_step,
      cum_sims       = cum_s,
      ESS            = round(1 / sum(w^2), 1),
      mean_takeoff   = round(sum(w * df$takeoff), 3),
      mean_deaths    = round(sum(w * df$n_deaths)),
      mean_hcw       = round(sum(w * df$n_hcw_deaths)),
      mean_duration  = round(sum(w * df$duration)),
      stringsAsFactors = FALSE
    )

    for (nm in param_names) {
      row[[paste0(nm, "_med")]] <- round(wq(df[[nm]], w, 0.500), 3)
      row[[paste0(nm, "_lo")]]  <- round(wq(df[[nm]], w, 0.025), 3)
      row[[paste0(nm, "_hi")]]  <- round(wq(df[[nm]], w, 0.975), 3)
    }

    rows[[length(rows) + 1L]] <- row
  }

  do.call(rbind, rows)
}


reconstruct_abc_result <- function(dir = getwd(),
                                   step = NULL,
                                   param_names = c("R0", "prop_funeral", "hcw_risk_scalar"),
                                   stat_names  = c("takeoff", "n_deaths", "n_hcw_deaths", "duration")) {

  step_num <- function(f) as.integer(sub(".*_step([0-9]+)$", "\\1", f))

  out_files <- list.files(dir, pattern = "^output_step[0-9]+$", full.names = TRUE)
  out_files <- out_files[order(step_num(out_files))]
  if (length(out_files) == 0L) stop("No output_step files found in ", dir)

  if (is.null(step)) {
    target_file <- out_files[length(out_files)]
    step <- step_num(target_file)
  } else {
    target_file <- file.path(dir, paste0("output_step", step))
    if (!file.exists(target_file)) {
      stop("output_step", step, " not found in ", dir, call. = FALSE)
    }
  }

  message("Reconstructing result from step ", step)

  df <- read.table(target_file, header = FALSE)
  colnames(df) <- c("weight", param_names, stat_names)

  tol_file <- file.path(dir, paste0("tolerance_step", step))
  epsilon  <- if (file.exists(tol_file)) as.numeric(readLines(tol_file)) else NA_real_

  sim_file <- file.path(dir, paste0("n_simul_tot_step", step))
  nsim     <- if (file.exists(sim_file)) as.numeric(readLines(sim_file)) else NA_real_

  # Normalisation SDs are computed from step 1 by ABC_sequential. Recompute the
  # same way so the reconstructed object behaves like a live one.
  step1_file <- file.path(dir, "output_step1")
  if (file.exists(step1_file)) {
    step1 <- read.table(step1_file, header = FALSE)
    stat_cols <- seq(2L + length(param_names), 1L + length(param_names) + length(stat_names))
    stats_norm <- apply(step1[, stat_cols], 2, sd)
  } else {
    stats_norm <- apply(df[, stat_names], 2, sd)
  }

  list(
    param                   = as.matrix(df[, param_names]),
    stats                   = as.matrix(df[, stat_names]),
    weights                 = df$weight / sum(df$weight),
    stats_normalization     = as.numeric(stats_norm),
    epsilon                 = epsilon,
    nsim                    = nsim,
    computime               = NA_real_,
    step_reconstructed_from = step
  )
}


