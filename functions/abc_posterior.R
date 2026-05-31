# abc_posterior.R
# -----------------------------------------------------------------------------
# Read and use a finished ABC-SMC fit. These helpers turn the on-disk output of
# an ABC_sequential() run (see abc_calibration_functions.R) into a usable
# posterior sample for any downstream analysis.
#
# Contents
#   find_latest_abc_run_dir()  : newest timestamped ABC run dir for a scenario.
#   read_abc_posterior_step()  : read an output_step<k> particle cloud.
#   downsample_posterior()     : weighted resample of the posterior particles.
#   derive_model_parameters()  : 3 fitted params -> 4 fiber model parameters,
#                                as a tidy one-row-per-set record. Delegates the
#                                mapping to map_abc_params_to_model() so the
#                                record can never drift from what is simulated.
#
# Requires abc_calibration_functions.R to be sourced first (for
# map_abc_params_to_model(), used by derive_model_parameters()).
# -----------------------------------------------------------------------------


# Find the most recent timestamped ABC run directory for a scenario. ABC runs
# are written to <abc_outputs>/<scenario_id>_YYYYMMDD_HHMMSS[...]; the directory
# names sort chronologically, so the last one is the newest.
find_latest_abc_run_dir <- function(abc_outputs_dir, scenario_id) {
  if (!dir.exists(abc_outputs_dir)) {
    stop("ABC outputs directory not found: ", abc_outputs_dir, call. = FALSE)
  }
  dirs <- list.dirs(abc_outputs_dir, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[grepl(scenario_id, basename(dirs), fixed = TRUE)]
  if (length(dirs) == 0L) {
    stop("No ABC run directories for scenario '", scenario_id, "' under ",
         abc_outputs_dir, call. = FALSE)
  }
  dirs[order(basename(dirs))][length(dirs)]
}

# Read one output_step<k> file written by EasyABC::ABC_sequential. If `step` is
# NULL the latest completed step is used. Returns a data.frame with the standard
# column names plus a "step" attribute.
read_abc_posterior_step <- function(run_dir,
                                    step = NULL,
                                    param_names = c("R0", "prop_funeral", "hcw_risk_scalar"),
                                    stat_names  = c("takeoff", "n_deaths", "n_hcw_deaths", "duration")) {
  step_of <- function(f) as.integer(sub(".*_step([0-9]+)$", "\\1", f))
  files   <- list.files(run_dir, pattern = "^output_step[0-9]+$", full.names = TRUE)
  if (length(files) == 0L) stop("No output_step files found in ", run_dir, call. = FALSE)
  files <- files[order(step_of(files))]

  target <- if (is.null(step)) {
    files[length(files)]
  } else {
    f <- file.path(run_dir, paste0("output_step", step))
    if (!file.exists(f)) stop("output_step", step, " not found in ", run_dir, call. = FALSE)
    f
  }

  df <- utils::read.table(target, header = FALSE)
  colnames(df) <- c("weight", param_names, stat_names)
  attr(df, "step")      <- step_of(target)
  attr(df, "step_file") <- target
  df
}

# Weighted resample of an ABC particle cloud (posterior). Particles are sampled
# WITH replacement with probability proportional to their ABC weights, which is
# the standard way to turn a weighted ABC-SMC population into an unweighted
# posterior sample. Returns the 3 fitted parameters for the drawn particles.
downsample_posterior <- function(posterior,
                                 n_sets,
                                 seed = 1,
                                 param_names = c("R0", "prop_funeral", "hcw_risk_scalar")) {
  if (!"weight" %in% names(posterior)) stop("`posterior` must contain a 'weight' column.", call. = FALSE)
  w <- posterior$weight / sum(posterior$weight)
  set.seed(seed)
  idx <- sample(seq_len(nrow(posterior)), size = n_sets, replace = TRUE, prob = w)
  out <- posterior[idx, param_names, drop = FALSE]
  out <- cbind(set_id = seq_len(n_sets), particle = idx, out)
  rownames(out) <- NULL
  out
}

# Convert the 3 fitted ABC parameters (R0, prop_funeral, hcw_risk_scalar) into
# the 4 fiber model parameters, as a tidy one-row-per-set data.frame for a
# transparent, saved record. The mapping itself is delegated to
# map_abc_params_to_model() (in abc_calibration_functions.R) -- the SAME function
# build_abc_model_args() uses to build the args fed to the simulator, so this
# record can never drift from what is actually simulated.
#
# D and F_fun are the scenario-level direct / funeral R0 multipliers returned by
# solve_offspring_means_for_R0(); pass the values you computed once for the
# scenario.
derive_model_parameters <- function(theta, D, F_fun, hcw_base_prob = 0.25) {
  model_pars <- map_abc_params_to_model(
    R0 = theta$R0, prop_funeral = theta$prop_funeral,
    hcw_risk_scalar = theta$hcw_risk_scalar,
    D = D, F_fun = F_fun, hcw_base_prob = hcw_base_prob
  )
  data.frame(
    set_id                        = theta$set_id,
    particle                      = theta$particle,
    R0                            = theta$R0,
    prop_funeral                  = theta$prop_funeral,
    hcw_risk_scalar               = theta$hcw_risk_scalar,
    mn_offspring_genPop           = model_pars$mn_offspring_genPop,
    mn_offspring_funeral          = model_pars$mn_offspring_funeral,
    prob_hcw_cond_genPop_hospital = model_pars$prob_hcw_cond_genPop_hospital,
    prob_hcw_cond_hcw_hospital    = model_pars$prob_hcw_cond_hcw_hospital,
    stringsAsFactors              = FALSE
  )
}
