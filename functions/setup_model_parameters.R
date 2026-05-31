# setup_model_parameters.R
# -----------------------------------------------------------------------------
# Build the full parameter set required by fiber's branching_process_main().
#
# Layers:
#   1. DEFAULT_SCALAR_INPUTS                 : literature-informed scalar
#                                              defaults (transmission means,
#                                              natural-history distributions,
#                                              setting/class probabilities,
#                                              simulation controls).
#   2. make_base_args(overrides = list(...)) : merges `overrides` into those
#                                              defaults and assembles the
#                                              scalar / distribution side of
#                                              the model argument list.
#   3. read_scenario_matrix()                : loads the dense time-varying
#                                              matrix (e.g. one row per time
#                                              step per scenario) from a CSV.
#   4. build_time_varying_args()             : turns the scenario rows into
#                                              the named list of time-varying
#                                              curves consumed by the model.
#   5. make_model_parameters()               : convenience wrapper that ties
#                                              (2)-(4) together using a single
#                                              `overrides` list, and returns
#                                              a ready-to-use combined args
#                                              list (`$args`) alongside the
#                                              base / time-varying pieces.
#
# Why overrides as a list. The set of scalar parameters in the model is
# expected to grow and change over time. Threading every overridable knob
# through named function arguments here would mean editing this file every
# time a new parameter is introduced. A single `overrides` list lets new
# entries appear in DEFAULT_SCALAR_INPUTS and immediately be overridable
# with no signature changes downstream.
#
# Requires fiber to be loaded:
#   library(fiber)
# This provides branching_process_main(), make_time_varying(), rtrunc_gamma(),
# prob_hosp_given_symptoms(), prob_death_given_symptoms(),
# offspring_function_*, complete_offspring_info(), summarise_output(),
# resolve_time_varying(), and hcw_loss_function().
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Default scalar parameters
# -----------------------------------------------------------------------------
# Central values from filovirus_model_parameter_table_4_scenarios_with_etu_baseline.xlsx.
# The parameter table itself remains the audit/reference document. Any of these
# can be overridden via the `overrides` list in make_base_args() /
# make_model_parameters().

DEFAULT_SCALAR_INPUTS <- list(
  # Hospital quarantine efficacy (fixed scalars). The offspring functions /
  # R0 approximation derive hospital_quarantine_efficacy(t) by mixing these two
  # by the (time-varying) prop_etu(t):
  #   hq_eff(t) = prop_etu(t) * etu_efficacy
  #             + (1 - prop_etu(t)) * general_hospital_quarantine_efficacy
  # See .hospital_quarantine_efficacy_t0() in calculate_model_approx_r0.R.
  etu_efficacy = 0.90,
  general_hospital_quarantine_efficacy = 0.30,
  # Per-PPE efficacy at preventing hospital transmission to HCWs (fixed scalar).
  # In fiber, PPE thinning is ppe_coverage_hcw(t) * ppe_efficacy; the time-varying
  # coverage lever ppe_coverage_hcw comes from the scenario ipc_helper curve (see
  # build_time_varying_args()). Does NOT enter the R0 approximation, but does
  # affect simulated HCW deaths.
  ppe_efficacy = 0.70,

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

  # Conditional class / setting assignment.
  # The two prob_hcw_cond_*_hospital values were previously asymmetric
  # (0.12 and 0.20). Those defaults were not literature-derived and have
  # been replaced with a symmetric 0.25 because (i) we now fit these via
  # ABC, with hcw_risk_scalar multiplying the symmetric base (see
  # hcw_base_prob in abc_calibration_functions_hcwRisk.R), and (ii) the previous
  # asymmetry was arbitrary.
  prob_hcw_cond_genPop_comm = 0.005,
  prob_hcw_cond_genPop_hospital = 0.25,   # was 0.12; symmetric ABC base, see note above
  prob_hcw_cond_hcw_comm = 0.02,
  prob_hcw_cond_hcw_hospital = 0.25,      # was 0.20; symmetric ABC base, see note above
  prob_hospital_cond_hcw_preAdm = 0.50,

  # Funeral control / assignment.
  safe_funeral_efficacy = 0.90,
  prob_hcw_cond_funeral_hcw = 0.02,
  prob_hcw_cond_funeral_genPop = 0.005,

  # Population / simulation control.
  population = 1000000,
  hcw_per_capita = 0.005,
  seeding_cases = 5,
  check_final_size = 30000,
  initial_immune = 0,
  susceptible_deplete = FALSE,

  # OBV post-exposure-prophylaxis gate. Off by default / during calibration;
  # switched on by the OBV impact analyses.
  obv_pep_enabled = FALSE
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

  required_cols <- c(
    "scenario", "scenario_label", "relative_day", "prob_hosp", "delay_hosp",
    "prob_unsafe_funeral_comm", "prob_unsafe_funeral_hosp",
    "prob_unsafe_funeral_etu", "prop_etu",
    "ipc_helper"   # drives the time-varying PPE coverage lever (ppe_coverage_hcw)
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
# Sanity check: verify the loaded fiber version supports the time-varying
# hospital_quarantine_efficacy(t) derivation used downstream.
# -----------------------------------------------------------------------------

check_model_function_version <- function() {
  if (!exists("branching_process_main", mode = "function")) {
    stop(
      "branching_process_main() is not on the search path. ",
      "Install fiber (devtools::install_github(\"petal-code/fiber\")) and ",
      "call library(fiber) before running this script.",
      call. = FALSE
    )
  }
  required_branching_args <- c("prop_etu", "etu_efficacy",
                               "general_hospital_quarantine_efficacy",
                               "ppe_coverage_hcw", "ppe_efficacy",
                               "safe_funeral_efficacy")
  missing_branching_args <- setdiff(
    required_branching_args, names(formals(branching_process_main))
  )
  if (length(missing_branching_args) > 0L) {
    stop(
      "The loaded branching_process_main() is missing required argument(s): ",
      paste(missing_branching_args, collapse = ", "), ". ",
      "Update fiber to the current time-varying branch.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}


# -----------------------------------------------------------------------------
# Scalar / distribution arguments
# -----------------------------------------------------------------------------
# `overrides` is a named list. Anything present in DEFAULT_SCALAR_INPUTS is
# replaced before the args list is assembled. Names that are NOT in
# DEFAULT_SCALAR_INPUTS still get attached to the returned list (so that, for
# example, a custom override of a derived field can pass through), but a
# warning is emitted to flag possible typos.

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
    # Raw hospitalisation delay has mean 1; the scenario delay_hosp curve
    # rescales this to the desired calendar-day mean.
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

    # Conditional class probabilities.
    prob_hcw_cond_genPop_comm = scalar_inputs$prob_hcw_cond_genPop_comm,
    prob_hcw_cond_genPop_hospital = scalar_inputs$prob_hcw_cond_genPop_hospital,
    prob_hcw_cond_hcw_comm = scalar_inputs$prob_hcw_cond_hcw_comm,
    prob_hcw_cond_hcw_hospital = scalar_inputs$prob_hcw_cond_hcw_hospital,
    prob_hospital_cond_hcw_preAdm = scalar_inputs$prob_hospital_cond_hcw_preAdm,

    # Hospital quarantine + PPE efficacies (fixed scalars). PPE coverage is a
    # separate time-varying lever (ppe_coverage_hcw) added in build_time_varying_args().
    etu_efficacy = scalar_inputs$etu_efficacy,
    general_hospital_quarantine_efficacy = scalar_inputs$general_hospital_quarantine_efficacy,
    ppe_efficacy = scalar_inputs$ppe_efficacy,

    # Safe funerals.
    safe_funeral_efficacy = scalar_inputs$safe_funeral_efficacy,
    prob_hcw_cond_funeral_hcw = scalar_inputs$prob_hcw_cond_funeral_hcw,
    prob_hcw_cond_funeral_genPop = scalar_inputs$prob_hcw_cond_funeral_genPop,

    # Simulation controls.
    population = scalar_inputs$population,
    hcw_per_capita = scalar_inputs$hcw_per_capita,
    check_final_size = scalar_inputs$check_final_size,
    initial_immune = scalar_inputs$initial_immune,
    seeding_cases = scalar_inputs$seeding_cases,
    susceptible_deplete = scalar_inputs$susceptible_deplete,
    obv_pep_enabled = scalar_inputs$obv_pep_enabled
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
# Reads the scenario rows for `scenario_id` and turns them into the named list
# of time-varying curves consumed by branching_process_main(). The returned
# list also carries `scenario_label` and `scenario_matrix` for traceability;
# those two are removed before the list is passed to the model.

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

  # Hospital deaths can occur in ordinary hospital care or ETUs. The matrix
  # has an ETU-specific unsafe-funeral probability; weight by p_ETU.
  p_unsafe_funeral_hosp_values <- clip01(
    (1 - scenario_matrix$prop_etu) * scenario_matrix$prob_unsafe_funeral_hosp +
      scenario_matrix$prop_etu * scenario_matrix$prob_unsafe_funeral_etu
  )

  # NOTE: hospital_quarantine_efficacy(t) is derived internally inside the
  # offspring functions by mixing the scalar etu_efficacy and
  # general_hospital_quarantine_efficacy (both emitted by make_base_args()) by
  # the time-varying prop_etu(t) curve below; we therefore pass prop_etu forward
  # rather than pre-computing a single hq curve here. PPE protection is split into
  # the scalar efficacy ppe_efficacy (in make_base_args()) and the time-varying
  # coverage lever ppe_coverage_hcw(t) built below from the ipc_helper curve;
  # fiber applies the two as ppe_coverage_hcw(t) * ppe_efficacy.

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

    prop_etu = make_curve(
      times, clip01(scenario_matrix$prop_etu), curve_method
    ),

    # PPE coverage lever: probability a relevant HCW is wearing PPE, ramping
    # with IPC/response maturity (the scenario ipc_helper curve). Paired with the
    # scalar ppe_efficacy from make_base_args() as ppe_coverage_hcw(t)*ppe_efficacy.
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
#   * names that match a DEFAULT_SCALAR_INPUTS entry  -> applied via
#     make_base_args() (these include the efficacy scalars etu_efficacy,
#     general_hospital_quarantine_efficacy and ppe_efficacy).
#   * names that match a time-varying-args entry      -> overwrite the
#     corresponding curve / scalar after build_time_varying_args().
#   * everything else                                 -> attached as extra
#     entries on the final combined args list (with a warning, since these
#     usually indicate a typo).

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
