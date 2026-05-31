# setup_model_parameters_new_approach.R
# -----------------------------------------------------------------------------
# NEW-APPROACH build of the parameter set required by fiber's
# branching_process_main(), targeting the *revamped* NPI parameterisation
# (coverage curves + fixed conditional-efficacy scalars + OBV PEP gate).
#
# What changed vs the HCW-risk version (analyses/02_ABC_model_fits_HCWrisk/
# helper_functions/setup_model_parameters.R):
#
#   * Hospital quarantine is now a prop_etu(t)-weighted mixture of TWO fixed
#     conditional-efficacy scalars:
#         hospital_quarantine_efficacy(t)
#             = prop_etu(t) * etu_efficacy
#             + (1 - prop_etu(t)) * general_hospital_quarantine_efficacy
#     (The old `etu_efficacy_baseline` + anchored-floor `ipc_helper` shape is
#     gone.)
#
#   * PPE is now coverage x efficacy: the time-varying lever is
#     `ppe_coverage_hcw(t)` (sourced from the scenario `ipc_helper` column) and
#     the conditional efficacy is the fixed scalar `ppe_efficacy`.
#
#   * The conditional efficacies are scalars in DEFAULT_SCALAR_INPUTS:
#         etu_efficacy                          (ABC-FITTED via the NPI scaler)
#         ppe_efficacy                          (ABC-FITTED via the NPI scaler)
#         general_hospital_quarantine_efficacy  (FIXED)
#         safe_funeral_efficacy                 (FIXED)
#     The two fitted ones are overridden per-particle by the NPI scaler in
#     abc_calibration_functions_new_approach.R; the two fixed ones are held at
#     their DEFAULT_SCALAR_INPUTS values.
#
#   * The HCW-risk scalar machinery is gone. prob_hcw_cond_*_hospital are held
#     fixed at their honest base (0.25); HCW infection volume is now governed
#     mechanistically by ppe_efficacy / etu_efficacy (see the design notes in
#     the run scripts).
#
# >>> PLACEHOLDERS <<<  The four efficacy scalars below are PLACEHOLDER values.
#     UPDATE WITH REAL NUMBERS before any production run. They are flagged again
#     at the point of definition.
#
# Requires fiber (the revamped NPI branch) to be loaded:
#   library(fiber)
# providing branching_process_main(), make_time_varying(), rtrunc_gamma(),
# resolve_time_varying(), offspring_function_*, summarise_output(), etc.
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Default scalar parameters
# -----------------------------------------------------------------------------
# Central / fixed values. Any of these can be overridden via the `overrides`
# list in make_base_args() / make_model_parameters().

DEFAULT_SCALAR_INPUTS <- list(
  # --- NPI conditional efficacies (scalars). ---------------------------------
  # >>> PLACEHOLDER VALUES — UPDATE WITH REAL NUMBERS. <<<
  # etu_efficacy and ppe_efficacy are the ABC-FITTED efficacies; the values
  # here are only used when the NPI scaler is NOT in play (e.g. a plain
  # make_model_parameters() call or the R0-invariant setup). During ABC they
  # are overwritten per-particle by the scaler.
  etu_efficacy                          = 0.70,   # PLACEHOLDER (fitted)
  ppe_efficacy                          = 0.55,   # PLACEHOLDER (fitted)
  # Held FIXED throughout the fit:
  general_hospital_quarantine_efficacy  = 0.30,   # PLACEHOLDER (fixed)
  safe_funeral_efficacy                 = 0.95,   # PLACEHOLDER (fixed)

  # OBV PEP gate off during calibration.
  obv_pep_enabled = FALSE,

  # Transmission means and dispersion.
  mn_offspring_genPop = 1.25,
  overdisp_offspring_genPop = 0.18,
  genpop_generation_mean = 15.4,
  genpop_generation_shape = 2.5,

  mn_offspring_hcw = 0.20,
  overdisp_offspring_hcw = 0.18,
  hcw_generation_mean = 15.4,
  hcw_generation_shape = 2.5,

  mn_offspring_funeral = 0.25,
  overdisp_offspring_funeral = 0.30,
  funeral_generation_shape = 20,
  funeral_generation_rate = 10,

  # Natural history.
  incubation_mean = 8.5,
  incubation_sd = 4.5,
  raw_onset_to_hospitalisation_mean = 1.0,
  raw_onset_to_hospitalisation_sd = 0.35,
  onset_to_death_mean = 9.3,
  onset_to_death_sd = 3.0,
  onset_to_recovery_mean = 13.0,
  onset_to_recovery_sd = 4.0,
  hospitalisation_to_death_mean = 4.5,
  hospitalisation_to_death_sd = 2.0,
  hospitalisation_to_recovery_mean = 8.0,
  hospitalisation_to_recovery_sd = 2.5,

  # Disease severity.
  prob_symptomatic = 1.0,
  prob_death_comm = 0.70,
  prob_death_hosp = 0.50,

  # Conditional class / setting assignment. HCW-given-hospital held at the
  # honest 0.25 base (no more fitted HCW-risk scalar).
  prob_hcw_cond_genPop_comm = 0.005,
  prob_hcw_cond_genPop_hospital = 0.25,
  prob_hcw_cond_hcw_comm = 0.02,
  prob_hcw_cond_hcw_hospital = 0.25,
  prob_hospital_cond_hcw_preAdm = 0.50,

  # Funeral assignment.
  prob_hcw_cond_funeral_hcw = 0.02,
  prob_hcw_cond_funeral_genPop = 0.005,

  # Population / simulation control.
  population = 1000000,
  hcw_per_capita = 0.005,
  seeding_cases = 5,
  check_final_size = 30000,
  initial_immune = 0,
  susceptible_deplete = FALSE
)


# -----------------------------------------------------------------------------
# Distribution helpers
# -----------------------------------------------------------------------------

gamma_shape_rate_from_mean_sd <- function(mean, sd) {
  if (mean <= 0 || sd <= 0) {
    stop("Gamma mean and sd must both be positive.", call. = FALSE)
  }
  list(shape = (mean / sd)^2, rate = mean / (sd^2))
}

make_gamma_sampler <- function(mean, sd) {
  pars <- gamma_shape_rate_from_mean_sd(mean = mean, sd = sd)
  function(n) {
    if (n <= 0L) return(numeric(0))
    rgamma(n = n, shape = pars$shape, rate = pars$rate)
  }
}


# -----------------------------------------------------------------------------
# Scenario matrix helpers
# -----------------------------------------------------------------------------

clip01 <- function(x) pmin(pmax(x, 0), 1)

read_scenario_matrix <- function(matrix_path) {
  if (missing(matrix_path) || is.null(matrix_path)) {
    stop("`matrix_path` is required: pass the path to a scenario CSV ",
         "(e.g. final_four_scenario_values.csv).", call. = FALSE)
  }
  if (!file.exists(matrix_path)) {
    stop("Cannot find scenario matrix CSV: ", matrix_path, call. = FALSE)
  }

  x <- read.csv(matrix_path, stringsAsFactors = FALSE, check.names = FALSE)

  numeric_cols <- setdiff(names(x), c("scenario", "scenario_label"))
  for (nm in numeric_cols) {
    x[[nm]] <- as.numeric(x[[nm]])
  }

  # `ipc_helper` is retained as the source column for ppe_coverage_hcw(t).
  required_cols <- c(
    "scenario", "scenario_label", "relative_day", "prob_hosp", "delay_hosp",
    "prob_unsafe_funeral_comm", "prob_unsafe_funeral_hosp",
    "prob_unsafe_funeral_etu", "prop_etu", "ipc_helper"
  )
  missing_cols <- setdiff(required_cols, names(x))
  if (length(missing_cols) > 0L) {
    stop("Scenario matrix is missing required column(s): ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  x
}

make_curve <- function(times, values, method = "linear") {
  make_time_varying(times = times, values = values, method = method)
}


# -----------------------------------------------------------------------------
# Sanity check: verify the loaded fiber version supports the revamped NPI args.
# -----------------------------------------------------------------------------
# NOTE: the exact branching_process_main() argument names below are inferred
# from the revamped offspring_function_genPop() / offspring_function_funeral()
# signatures. If fiber routes them under different top-level names, update this
# vector (and make_base_args() / build_time_varying_args()) accordingly; this
# check is here precisely so a mismatch fails fast with a clear message.

check_model_function_version <- function() {
  if (!exists("branching_process_main", mode = "function")) {
    stop(
      "branching_process_main() is not on the search path. ",
      "Install the revamped fiber and call library(fiber) first.",
      call. = FALSE
    )
  }
  required_branching_args <- c(
    "prop_etu", "etu_efficacy", "general_hospital_quarantine_efficacy",
    "ppe_coverage_hcw", "ppe_efficacy", "safe_funeral_efficacy"
  )
  missing_branching_args <- setdiff(
    required_branching_args, names(formals(branching_process_main))
  )
  if (length(missing_branching_args) > 0L) {
    stop(
      "The loaded branching_process_main() is missing required argument(s): ",
      paste(missing_branching_args, collapse = ", "), ". ",
      "This setup targets the revamped NPI parameterisation (coverage + fixed ",
      "efficacy scalars). Update fiber, or rename these args to match.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}


# -----------------------------------------------------------------------------
# Scalar / distribution arguments
# -----------------------------------------------------------------------------

make_base_args <- function(overrides = list()) {
  if (!is.list(overrides)) {
    stop("`overrides` must be a (possibly empty) named list.", call. = FALSE)
  }
  if (length(overrides) > 0L && is.null(names(overrides))) {
    stop("`overrides` must be a NAMED list.", call. = FALSE)
  }

  unknown <- setdiff(names(overrides), names(DEFAULT_SCALAR_INPUTS))
  if (length(unknown) > 0L) {
    warning(
      "Override(s) not in DEFAULT_SCALAR_INPUTS will still be applied to ",
      "the returned args list: ", paste(unknown, collapse = ", "), ".",
      call. = FALSE
    )
  }

  scalar_inputs <- utils::modifyList(DEFAULT_SCALAR_INPUTS, overrides)

  args <- list(
    # Transmission.
    mn_offspring_genPop = scalar_inputs$mn_offspring_genPop,
    overdisp_offspring_genPop = scalar_inputs$overdisp_offspring_genPop,
    Tg_shape_genPop = scalar_inputs$genpop_generation_shape,
    Tg_rate_genPop = scalar_inputs$genpop_generation_shape /
      scalar_inputs$genpop_generation_mean,

    mn_offspring_hcw = scalar_inputs$mn_offspring_hcw,
    overdisp_offspring_hcw = scalar_inputs$overdisp_offspring_hcw,
    Tg_shape_hcw = scalar_inputs$hcw_generation_shape,
    Tg_rate_hcw = scalar_inputs$hcw_generation_shape /
      scalar_inputs$hcw_generation_mean,

    mn_offspring_funeral = scalar_inputs$mn_offspring_funeral,
    overdisp_offspring_funeral = scalar_inputs$overdisp_offspring_funeral,
    Tg_shape_funeral = scalar_inputs$funeral_generation_shape,
    Tg_rate_funeral = scalar_inputs$funeral_generation_rate,

    # Natural history (gamma samplers).
    incubation_period = make_gamma_sampler(
      mean = scalar_inputs$incubation_mean,
      sd = scalar_inputs$incubation_sd
    ),
    onset_to_hospitalisation = make_gamma_sampler(
      mean = scalar_inputs$raw_onset_to_hospitalisation_mean,
      sd = scalar_inputs$raw_onset_to_hospitalisation_sd
    ),
    onset_to_death = make_gamma_sampler(
      mean = scalar_inputs$onset_to_death_mean,
      sd = scalar_inputs$onset_to_death_sd
    ),
    onset_to_recovery = make_gamma_sampler(
      mean = scalar_inputs$onset_to_recovery_mean,
      sd = scalar_inputs$onset_to_recovery_sd
    ),
    hospitalisation_to_death = make_gamma_sampler(
      mean = scalar_inputs$hospitalisation_to_death_mean,
      sd = scalar_inputs$hospitalisation_to_death_sd
    ),
    hospitalisation_to_recovery = make_gamma_sampler(
      mean = scalar_inputs$hospitalisation_to_recovery_mean,
      sd = scalar_inputs$hospitalisation_to_recovery_sd
    ),

    # Disease severity and healthcare seeking.
    prob_symptomatic = scalar_inputs$prob_symptomatic,
    prob_death_comm = scalar_inputs$prob_death_comm,
    prob_death_hosp = scalar_inputs$prob_death_hosp,

    # Conditional class probabilities (fixed; no HCW-risk scalar).
    prob_hcw_cond_genPop_comm = scalar_inputs$prob_hcw_cond_genPop_comm,
    prob_hcw_cond_genPop_hospital = scalar_inputs$prob_hcw_cond_genPop_hospital,
    prob_hcw_cond_hcw_comm = scalar_inputs$prob_hcw_cond_hcw_comm,
    prob_hcw_cond_hcw_hospital = scalar_inputs$prob_hcw_cond_hcw_hospital,
    prob_hospital_cond_hcw_preAdm = scalar_inputs$prob_hospital_cond_hcw_preAdm,

    # --- Revamped NPI conditional-efficacy scalars. ---
    # PPE and ETU are overwritten per-particle by the NPI scaler during ABC;
    # general-hospital quarantine and safe-funeral efficacies are fixed.
    etu_efficacy = scalar_inputs$etu_efficacy,
    ppe_efficacy = scalar_inputs$ppe_efficacy,
    general_hospital_quarantine_efficacy =
      scalar_inputs$general_hospital_quarantine_efficacy,
    safe_funeral_efficacy = scalar_inputs$safe_funeral_efficacy,

    # OBV PEP gate (off during calibration).
    obv_pep_enabled = scalar_inputs$obv_pep_enabled,

    # Funeral assignment.
    prob_hcw_cond_funeral_hcw = scalar_inputs$prob_hcw_cond_funeral_hcw,
    prob_hcw_cond_funeral_genPop = scalar_inputs$prob_hcw_cond_funeral_genPop,

    # Simulation controls.
    population = scalar_inputs$population,
    hcw_per_capita = scalar_inputs$hcw_per_capita,
    check_final_size = scalar_inputs$check_final_size,
    initial_immune = scalar_inputs$initial_immune,
    seeding_cases = scalar_inputs$seeding_cases,
    susceptible_deplete = scalar_inputs$susceptible_deplete
  )

  # Attach any unknown overrides so they still reach the model args list.
  for (nm in unknown) {
    args[[nm]] <- overrides[[nm]]
  }

  args
}


# -----------------------------------------------------------------------------
# Time-varying arguments
# -----------------------------------------------------------------------------
# Turns the scenario rows for `scenario_id` into the named list of time-varying
# curves consumed by the revamped branching_process_main(). The conditional
# efficacies are NO LONGER time-varying here — they are scalars in base_args.
# The time variation now lives entirely in the COVERAGE curves prop_etu(t) and
# ppe_coverage_hcw(t).

build_time_varying_args <- function(
    scenario_id,
    matrix,
    curve_method = "linear"
) {
  if (missing(matrix) || is.null(matrix)) {
    stop("`matrix` is required: pass a data frame from read_scenario_matrix().",
         call. = FALSE)
  }

  scenario_matrix <- matrix[matrix$scenario == scenario_id, ]
  if (nrow(scenario_matrix) == 0L) {
    stop("No rows found for scenario_id = ", scenario_id, ". ",
         "Available: ", paste(unique(matrix$scenario), collapse = ", "), ".",
         call. = FALSE)
  }

  scenario_matrix <- scenario_matrix[order(scenario_matrix$relative_day), ]
  times <- scenario_matrix$relative_day

  # Hospital deaths can occur in ordinary hospital care or ETUs. The matrix has
  # an ETU-specific unsafe-funeral probability; weight by p_ETU.
  p_unsafe_funeral_hosp_values <- clip01(
    (1 - scenario_matrix$prop_etu) * scenario_matrix$prob_unsafe_funeral_hosp +
      scenario_matrix$prop_etu * scenario_matrix$prob_unsafe_funeral_etu
  )

  list(
    scenario_label = unique(scenario_matrix$scenario_label)[1L],
    scenario_matrix = scenario_matrix,

    prob_hospitalised_genPop = make_curve(
      times, clip01(scenario_matrix$prob_hosp), curve_method
    ),
    prob_hospitalised_hcw = make_curve(
      times, clip01(scenario_matrix$prob_hosp), curve_method
    ),

    hospitalisation_delay_factor = make_curve(
      times, pmax(scenario_matrix$delay_hosp, 0.01), curve_method
    ),

    p_unsafe_funeral_comm_genPop = make_curve(
      times, clip01(scenario_matrix$prob_unsafe_funeral_comm), curve_method
    ),
    p_unsafe_funeral_comm_hcw = make_curve(
      times, clip01(scenario_matrix$prob_unsafe_funeral_comm), curve_method
    ),
    p_unsafe_funeral_hosp_genPop = make_curve(
      times, p_unsafe_funeral_hosp_values, curve_method
    ),
    p_unsafe_funeral_hosp_hcw = make_curve(
      times, p_unsafe_funeral_hosp_values, curve_method
    ),

    # Coverage levers (time-varying). prop_etu drives the hospital-quarantine
    # mixture; ppe_coverage_hcw (sourced from the scenario `ipc_helper` column)
    # is the probability an HCW recipient has PPE at time t.
    prop_etu = make_curve(
      times, clip01(scenario_matrix$prop_etu), curve_method
    ),
    ppe_coverage_hcw = make_curve(
      times, clip01(scenario_matrix$ipc_helper), curve_method
    )
  )
}


# -----------------------------------------------------------------------------
# High-level wrapper
# -----------------------------------------------------------------------------
# One call produces the ready-to-use model args from a scenario + overrides.
#
# Override routing for entries in `overrides`:
#   * names matching a DEFAULT_SCALAR_INPUTS entry -> applied via make_base_args()
#   * names matching a time-varying-args entry     -> overwrite that curve
#   * everything else                              -> attached (with a warning)

make_model_parameters <- function(
    scenario_id,
    scenario_matrix,
    overrides = list(),
    curve_method = "linear"
) {
  if (!is.list(overrides)) {
    stop("`overrides` must be a (possibly empty) named list.", call. = FALSE)
  }
  if (length(overrides) > 0L && is.null(names(overrides))) {
    stop("`overrides` must be a NAMED list.", call. = FALSE)
  }

  scalar_names <- names(DEFAULT_SCALAR_INPUTS)
  scalar_overrides <- overrides[intersect(names(overrides), scalar_names)]

  base_args <- make_base_args(overrides = scalar_overrides)

  tv_args_full <- build_time_varying_args(
    scenario_id = scenario_id,
    matrix = scenario_matrix,
    curve_method = curve_method
  )
  scenario_label <- tv_args_full$scenario_label
  scenario_matrix_used <- tv_args_full$scenario_matrix
  tv_args <- tv_args_full[setdiff(
    names(tv_args_full), c("scenario_label", "scenario_matrix")
  )]

  tv_names <- names(tv_args)
  tv_overrides <- overrides[intersect(names(overrides), tv_names)]
  if (length(tv_overrides) > 0L) {
    tv_args <- utils::modifyList(tv_args, tv_overrides)
  }

  args <- c(base_args, tv_args)

  remaining_overrides <- overrides[setdiff(
    names(overrides),
    c(scalar_names, tv_names)
  )]
  if (length(remaining_overrides) > 0L) {
    warning(
      "Override(s) not in DEFAULT_SCALAR_INPUTS or the time-varying args ",
      "will still be applied (likely a typo?): ",
      paste(names(remaining_overrides), collapse = ", "), ".",
      call. = FALSE
    )
    args <- utils::modifyList(args, remaining_overrides)
  }

  list(
    args = args,
    base_args = base_args,
    tv_args = tv_args,
    scenario_label = scenario_label,
    scenario_matrix = scenario_matrix_used,
    overrides_applied = overrides
  )
}
