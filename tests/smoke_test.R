#!/usr/bin/env Rscript
# =============================================================================
# smoke_test.R
# -----------------------------------------------------------------------------
# A self-contained smoke test for the refactors in this branch. It exercises
# everything that was changed without needing a full ABC run, and degrades
# gracefully when the fiber package (the model engine) is not installed.
#
# What it covers
#   1. All functions/*.R files parse and the public functions exist.
#   2. make_base_args(): the new fixed efficacy scalars are emitted
#      (etu_efficacy 0.90, general_hospital_quarantine_efficacy 0.30,
#      ppe_efficacy 0.70) and the dropped legacy ones are gone.
#   3. read_scenario_matrix(): reads the real CSV, requires ipc_helper (it now
#      drives the PPE coverage lever) plus the other scenario columns, and
#      errors when any required column is missing.
#   4. The ABC parameter -> model mapping, and the "no drift" invariant that
#      map_abc_params_to_model() / build_abc_model_args() / derive_model_
#      parameters() all agree.
#   5. The reformulated hospital_quarantine_efficacy(t) mixing and the R0
#      solver round-trip (solve -> plug back -> recover target R0). Pure base R.
#   6. io_helpers: find_latest_file() / list_files_matching() across name- vs
#      mtime-sorting, substring vs regex, case, no-match, dir exclusion.
#   7. find_latest_abc_run_dir() + the on-disk posterior readers, against the
#      real abc_outputs/ if present (skipped on a fresh clone, where it's
#      gitignored).
#   8. Output-location helpers: make_abc_output_dir() timestamp format,
#      with_abc_output_dir() cwd restoration, and the new %Y%m%d_%H%M%S result
#      filename sorting chronologically via find_latest_file().
#   9. simulation_helpers: bin_counts(), q_summary(), %||%.
#  10. (fiber only) The live model interface: check_model_function_version(),
#      ppe_efficacy/ppe_coverage_hcw in the signature, and one tiny end-to-end run
#      make_model_parameters() -> solve -> build_abc_model_args() ->
#      branching_process_main() -> abc_summarise()/simulate_one().
#
# Run it:
#   Rscript tests/smoke_test.R
# Exit status is 0 if every executed check passed, 1 if any failed (skips are
# not failures). Sourcing it interactively prints the same report but does not
# quit the session.
# =============================================================================


# ----------------------------------------------------------------------------
# Locate the repo root (independent of the working directory it's launched from)
# ----------------------------------------------------------------------------
find_repo_root <- function() {
  marker <- "obv_hcw_paper.Rproj"
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
  start <- if (length(file_arg) == 1L && nzchar(file_arg)) {
    dirname(normalizePath(file_arg, mustWork = FALSE))
  } else {
    getwd()
  }
  d <- normalizePath(start, mustWork = FALSE)
  repeat {
    if (file.exists(file.path(d, marker))) return(d)
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  if (requireNamespace("here", quietly = TRUE)) return(here::here())
  stop("Could not locate repo root (looked for ", marker, ").", call. = FALSE)
}

REPO <- find_repo_root()
FN   <- file.path(REPO, "functions")


# ----------------------------------------------------------------------------
# Minimal test harness
# ----------------------------------------------------------------------------
.smoke <- new.env(parent = emptyenv())
.smoke$pass <- 0L
.smoke$fail <- 0L
.smoke$skip <- 0L
.smoke$failed <- character(0)

section <- function(title) cat(sprintf("\n──── %s ────\n", title))

pass <- function(desc) { .smoke$pass <- .smoke$pass + 1L; cat(sprintf("  [PASS] %s\n", desc)); invisible(TRUE) }
fail <- function(desc, msg = "") {
  .smoke$fail <- .smoke$fail + 1L
  .smoke$failed <- c(.smoke$failed, desc)
  cat(sprintf("  [FAIL] %s%s\n", desc, if (nzchar(msg)) sprintf("\n         └─ %s", msg) else ""))
  invisible(FALSE)
}
skip <- function(desc, msg = "") {
  .smoke$skip <- .smoke$skip + 1L
  cat(sprintf("  [SKIP] %s%s\n", desc, if (nzchar(msg)) sprintf("  (%s)", msg) else ""))
  invisible(TRUE)
}

# PASS iff `expr` evaluates to TRUE with no error; FAIL on error or non-TRUE.
test <- function(desc, expr) {
  v <- tryCatch(expr, error = function(e) e)
  if (inherits(v, "condition")) return(fail(desc, conditionMessage(v)))
  if (isTRUE(v)) return(pass(desc))
  fail(desc, sprintf("expected TRUE, got: %s",
                     paste(utils::head(as.character(v), 3L), collapse = ", ")))
}

# PASS iff evaluating `expr` raises an error (confirms a guard fires).
test_error <- function(desc, expr) {
  raised <- tryCatch({ expr; FALSE }, error = function(e) TRUE)
  if (isTRUE(raised)) pass(desc) else fail(desc, "expected an error, none raised")
}

# PASS iff evaluating `expr` raises a warning.
test_warns <- function(desc, expr) {
  got <- FALSE
  withCallingHandlers(
    tryCatch(expr, error = function(e) NULL),
    warning = function(w) { got <<- TRUE; invokeRestart("muffleWarning") }
  )
  if (isTRUE(got)) pass(desc) else fail(desc, "expected a warning, none raised")
}

approx_eq <- function(a, b, tol = 1e-8) {
  is.numeric(a) && is.numeric(b) && length(a) == length(b) &&
    all(is.finite(a)) && all(is.finite(b)) && all(abs(a - b) <= tol)
}

# Reused objects (pre-declared so the compute-and-assert idiom can <<- into them).
# NB: avoid names that exist on the search path (e.g. stats::ar) — <<- would
# walk past globalenv and hit the locked package binding. Hence `r0args`.
ba <- sm <- r0args <- sol <- r0 <- rd <- ps <- mpf <- solf <- NULL


# ----------------------------------------------------------------------------
cat("=============================================================\n")
cat(" obv_hcw_paper :: smoke test\n")
cat("=============================================================\n")
cat(sprintf("  R           : %s\n", getRversion()))
cat(sprintf("  repo root   : %s\n", REPO))
has_fiber <- requireNamespace("fiber", quietly = TRUE)
cat(sprintf("  fiber pkg   : %s\n", if (has_fiber) "available (live model checks ON)"
                                    else "NOT installed (live model checks will SKIP)"))


# ----------------------------------------------------------------------------
section("1. Source functions/ + public functions exist")
# ----------------------------------------------------------------------------
# Order matters: abc_calibration needs setup + r0; abc_posterior needs the
# mapping from abc_calibration. io_helpers is dependency-free.
src_files <- c("io_helpers.R", "setup_model_parameters.R",
               "calculate_model_approx_r0.R", "abc_calibration_functions.R",
               "abc_posterior.R", "simulation_helpers.R")
src_ok <- TRUE
for (f in src_files) {
  ok <- test(sprintf("source functions/%s", f),
             { source(file.path(FN, f)); TRUE })
  src_ok <- src_ok && isTRUE(ok)
}

expected_fns <- c(
  "make_base_args", "make_model_parameters", "read_scenario_matrix",
  "build_time_varying_args", "make_curve", "check_model_function_version",
  "map_abc_params_to_model", "build_abc_model_args", "make_abc_output_dir",
  "with_abc_output_dir", "reconstruct_abc_result", "abc_summarise",
  "find_latest_abc_run_dir", "read_abc_posterior_step", "downsample_posterior",
  "derive_model_parameters", "find_latest_file", "list_files_matching",
  "simulate_one", "bin_counts", "q_summary",
  "R0_single_type_from_args", "solve_offspring_means_for_R0",
  ".hospital_quarantine_efficacy_t0"
)
test("all expected public functions are defined",
     all(vapply(expected_fns, exists, logical(1), mode = "function")))

if (!src_ok) {
  cat("\n  Core files failed to source — cannot run the rest. See failures above.\n")
  section("Summary")
  cat(sprintf("\n  PASS: %d   FAIL: %d   SKIP: %d\n", .smoke$pass, .smoke$fail, .smoke$skip))
  if (!interactive()) quit(status = 1L, save = "no")
}


# ----------------------------------------------------------------------------
section("2. make_base_args(): new efficacy scalars wired, legacy ones dropped")
# ----------------------------------------------------------------------------
test("make_base_args() builds without error", { ba <<- make_base_args(); is.list(ba) })
test("emits etu_efficacy = 0.90",                          approx_eq(ba$etu_efficacy, 0.90))
test("emits general_hospital_quarantine_efficacy = 0.30",  approx_eq(ba$general_hospital_quarantine_efficacy, 0.30))
test("emits ppe_efficacy = 0.70",                          approx_eq(ba$ppe_efficacy, 0.70))
test("does NOT emit legacy ppe_efficacy_hcw",              is.null(ba[["ppe_efficacy_hcw"]]))
test("does NOT emit legacy etu_efficacy_baseline",         is.null(ba$etu_efficacy_baseline))
test("does NOT emit legacy ipc_helper (it's a scenario column, not a base arg)",
     is.null(ba$ipc_helper))
test("scalar override of etu_efficacy is applied",
     approx_eq(suppressWarnings(make_base_args(list(etu_efficacy = 0.5)))$etu_efficacy, 0.5))
test_warns("unknown override name warns", make_base_args(list(not_a_real_param = 1)))


# ----------------------------------------------------------------------------
section("3. read_scenario_matrix(): required columns incl. ipc_helper")
# ----------------------------------------------------------------------------
csv <- file.path(REPO, "data-processed", "final_four_scenario_values.csv")
test("scenario CSV exists", file.exists(csv))
test("read_scenario_matrix() reads the real CSV",
     { sm <<- read_scenario_matrix(csv); is.data.frame(sm) && nrow(sm) > 0L })
req_cols <- c("scenario", "scenario_label", "relative_day", "prob_hosp",
              "delay_hosp", "prob_unsafe_funeral_comm", "prob_unsafe_funeral_hosp",
              "prob_unsafe_funeral_etu", "prop_etu", "ipc_helper")
test("matrix has all required columns", all(req_cols %in% names(sm)))
test("run-script scenarios present in matrix$scenario",
     all(c("Worst_WestAfrica", "Middle_DRC_ConflictSmoothed") %in% sm$scenario))

raw <- read.csv(csv, check.names = FALSE, stringsAsFactors = FALSE)
f_no_ipc <- tempfile(fileext = ".csv")
write.csv(raw[setdiff(names(raw), "ipc_helper")], f_no_ipc, row.names = FALSE)
test_error("errors WITHOUT ipc_helper (now required: drives ppe_coverage_hcw)",
           read_scenario_matrix(f_no_ipc))

f_no_req <- tempfile(fileext = ".csv")
write.csv(raw[setdiff(names(raw), "prop_etu")], f_no_req, row.names = FALSE)
test_error("still errors when a genuinely required column (prop_etu) is missing",
           read_scenario_matrix(f_no_req))


# ----------------------------------------------------------------------------
section("4. ABC param mapping + the single-source-of-truth invariant")
# ----------------------------------------------------------------------------
# Hand-computed expectations for R0=2, prop_funeral=0.25, D=0.8, F=0.5,
# hcw_base_prob=0.25, hcw_risk_scalar=2.
m <- map_abc_params_to_model(R0 = 2, prop_funeral = 0.25, hcw_risk_scalar = 2,
                             D = 0.8, F_fun = 0.5, hcw_base_prob = 0.25)
test("map: mn_offspring_genPop  = (1-pf)*R0/D",  approx_eq(m$mn_offspring_genPop,  0.75 * 2 / 0.8))
test("map: mn_offspring_funeral =     pf *R0/F",  approx_eq(m$mn_offspring_funeral, 0.25 * 2 / 0.5))
test("map: prob_hcw_*_hospital = base*scalar",    approx_eq(m$prob_hcw_cond_genPop_hospital, 0.5))
test("map: both hcw-hospital probs equal",        approx_eq(m$prob_hcw_cond_genPop_hospital, m$prob_hcw_cond_hcw_hospital))
test("map: hcw prob capped at 1.0",
     approx_eq(map_abc_params_to_model(2, 0.25, 10, 0.8, 0.5, 0.25)$prob_hcw_cond_genPop_hospital, 1.0))

# Invariant: the simulator path (build_abc_model_args) and the saved-record
# path (derive_model_parameters) must derive identical numbers from the mapping.
bargs <- build_abc_model_args(R0 = 2, prop_funeral = 0.25, hcw_risk_scalar = 2,
                              base = list(dummy = 1, seed = 999), tv = list(dummy2 = 2),
                              D = 0.8, F_fun = 0.5, seeding_cases = 7, hcw_base_prob = 0.25)
test("build_abc_model_args splices the same genPop mean",  approx_eq(bargs$mn_offspring_genPop,  m$mn_offspring_genPop))
test("build_abc_model_args splices the same funeral mean", approx_eq(bargs$mn_offspring_funeral, m$mn_offspring_funeral))
test("build_abc_model_args sets seeding_cases + clears seed",
     isTRUE(bargs$seeding_cases == 7) && is.null(bargs[["seed"]]))

theta <- list(set_id = 1L, particle = 1L, R0 = 2, prop_funeral = 0.25, hcw_risk_scalar = 2)
dp <- derive_model_parameters(theta, D = 0.8, F_fun = 0.5, hcw_base_prob = 0.25)
test("derive_model_parameters matches mapping (genPop)",  approx_eq(dp$mn_offspring_genPop,        m$mn_offspring_genPop))
test("derive_model_parameters matches mapping (funeral)", approx_eq(dp$mn_offspring_funeral,       m$mn_offspring_funeral))
test("derive_model_parameters matches mapping (hcw)",     approx_eq(dp$prob_hcw_cond_hcw_hospital, m$prob_hcw_cond_hcw_hospital))


# ----------------------------------------------------------------------------
section("5. Reformulated hospital_quarantine_efficacy(t) + R0 round-trip")
# ----------------------------------------------------------------------------
# Build a fiber-free args list: scalar defaults + gamma samplers from
# make_base_args(), plus the time-varying inputs as plain scalars (the R0 code
# resolves scalars or function(t) via .at_t0()).
test("assemble fiber-free R0 args from make_base_args()",
     {
       a <- make_base_args()
       a$prob_hospitalised_genPop     <- 0.5
       a$hospitalisation_delay_factor <- 1.0
       a$p_unsafe_funeral_comm_genPop <- 0.3
       a$p_unsafe_funeral_hosp_genPop <- 0.2
       a$prop_etu                     <- 0.4
       r0args <<- a
       is.list(r0args)
     })
# Mix = prop_etu*etu_efficacy + (1-prop_etu)*general = .4*.9 + .6*.3 = 0.54
test("hospital_quarantine_efficacy(0) = mix of the new triplet",
     approx_eq(.hospital_quarantine_efficacy_t0(r0args), 0.54))
test_error("errors if etu_efficacy missing (new triplet required)",
           { bad <- r0args; bad$etu_efficacy <- NULL; .hospital_quarantine_efficacy_t0(bad) })
test_error("errors if general_hospital_quarantine_efficacy missing",
           { bad <- r0args; bad$general_hospital_quarantine_efficacy <- NULL; .hospital_quarantine_efficacy_t0(bad) })

test("solve_offspring_means_for_R0() runs",
     { sol <<- solve_offspring_means_for_R0(R0 = 1.5, args = r0args,
                                            proportion_transmission_from_funerals = 0.3,
                                            n = 20000, seed = 123); is.list(sol) })
test("solver D multiplier in (0, 1]",
     is.finite(sol$D_direct_multiplier) && sol$D_direct_multiplier > 0 && sol$D_direct_multiplier <= 1)
test("solver F multiplier > 0",
     is.finite(sol$F_funeral_multiplier) && sol$F_funeral_multiplier > 0)
test("solver reports hq(0) = 0.54", approx_eq(sol$hospital_quarantine_efficacy_t0, 0.54))
# Round-trip: plug solved means back in (same seed + n => identical MC draws).
test("R0 round-trip recovers the target R0",
     {
       a2 <- r0args
       a2$mn_offspring_genPop  <- sol$mn_offspring_genPop_required
       a2$mn_offspring_funeral <- sol$mn_offspring_funeral_required
       r0 <<- R0_single_type_from_args(a2, n = 20000, seed = 123)
       approx_eq(r0$R0, 1.5, tol = 1e-6)
     })
test("R0 output hq(0) matches the mix", approx_eq(r0$hospital_quarantine_efficacy_t0, 0.54))


# ----------------------------------------------------------------------------
section("6. io_helpers: find_latest_file() / list_files_matching()")
# ----------------------------------------------------------------------------
iod <- file.path(tempdir(), sprintf("smoke_io_%d", as.integer(Sys.time())))
dir.create(iod, recursive = TRUE, showWarnings = FALSE)
io_names <- c(
  "fiber_ABC_SMC_Worst_WestAfrica_HCWrisk_20260101_090000.rds",
  "fiber_ABC_SMC_Worst_WestAfrica_HCWrisk_20260203_120000.rds",
  "fiber_ABC_SMC_Worst_WestAfrica_HCWrisk_20260203_130000.rds",
  "fiber_ABC_SMC_Middle_DRC_ConflictSmoothed_HCWrisk_20260204_000000.rds"
)
for (n in io_names) writeLines("x", file.path(iod, n))
dir.create(file.path(iod, "a_subdir"), showWarnings = FALSE)

test("by name: picks the newest WestAfrica file",
     grepl("20260203_130000", basename(find_latest_file(iod, "Worst_WestAfrica"))))
test("list_files_matching: 3 WestAfrica hits, newest first",
     { L <- list_files_matching(iod, "Worst_WestAfrica"); length(L) == 3L && grepl("20260203_130000", basename(L[1])) })
test("substring match selects the DRC file",
     grepl("Middle_DRC", basename(find_latest_file(iod, "DRC"))))
test("regex match (fixed = FALSE)",
     grepl("20260204_000000", basename(find_latest_file(iod, "[0-9]{8}_000000", fixed = FALSE))))
test("ignore_case match",
     grepl("Worst_WestAfrica", basename(find_latest_file(iod, "westafrica", ignore_case = TRUE))))
test("sub-directories are excluded from results",
     !any(grepl("a_subdir", basename(list_files_matching(iod, NULL)))))
test_error("errors on no match when error = TRUE (default)",
           find_latest_file(iod, "NoSuchPattern"))
test("returns NA on no match when error = FALSE",
     is.na(find_latest_file(iod, "NoSuchPattern", error = FALSE)))
test("list_files_matching returns length 0 on no match (no error)",
     length(list_files_matching(iod, "NoSuchPattern")) == 0L)
# Make the oldest-named file the newest by mtime, to prove name vs mtime diverge.
Sys.setFileTime(file.path(iod, io_names[1]), Sys.time() + 1000)
test("by = 'mtime' picks the newest mtime (diverges from name order)",
     grepl("20260101_090000", basename(find_latest_file(iod, "Worst_WestAfrica", by = "mtime"))))
test("by = 'name' still picks the newest name",
     grepl("20260203_130000", basename(find_latest_file(iod, "Worst_WestAfrica", by = "name"))))


# ----------------------------------------------------------------------------
section("7. find_latest_abc_run_dir() + on-disk posterior readers")
# ----------------------------------------------------------------------------
abc_dir <- file.path(REPO, "analyses", "02_ABC_model_fits_HCWrisk", "abc_outputs")
run_dirs <- if (dir.exists(abc_dir)) list.dirs(abc_dir, recursive = FALSE) else character(0)
scen_present <- c("Worst_WestAfrica", "Middle_DRC_ConflictSmoothed")
scen_present <- scen_present[vapply(scen_present,
                                    function(s) any(grepl(s, basename(run_dirs), fixed = TRUE)),
                                    logical(1))]
if (length(run_dirs) == 0L || length(scen_present) == 0L) {
  skip("find_latest_abc_run_dir() + posterior readers",
       "abc_outputs/ not present — it's gitignored; produced by an ABC run")
} else {
  scen <- scen_present[1]
  test(sprintf("find_latest_abc_run_dir() finds a %s run", scen),
       { rd <<- find_latest_abc_run_dir(abc_dir, scen); dir.exists(rd) })
  test("read_abc_posterior_step() reads latest step (weight+3 params+4 stats)",
       { ps <<- read_abc_posterior_step(rd); is.data.frame(ps) && ncol(ps) == 8L && "R0" %in% names(ps) })
  test("downsample_posterior() draws the requested n with the param columns",
       { ds <- downsample_posterior(ps, n_sets = 50L);
         nrow(ds) == 50L && all(c("set_id", "R0", "prop_funeral", "hcw_risk_scalar") %in% names(ds)) })
  test("reconstruct_abc_result() rebuilds a result object",
       { rr <- reconstruct_abc_result(rd); is.list(rr) && all(c("param", "stats", "weights") %in% names(rr)) })
}


# ----------------------------------------------------------------------------
section("8. Output-location helpers + new timestamped filename")
# ----------------------------------------------------------------------------
od_base <- file.path(tempdir(), sprintf("smoke_odir_%d", as.integer(Sys.time())))
test("make_abc_output_dir() creates a %Y%m%d_%H%M%S-stamped dir",
     { od <- make_abc_output_dir(od_base, "Worst_WestAfrica");
       dir.exists(od) && grepl("^Worst_WestAfrica_[0-9]{8}_[0-9]{6}$", basename(od)) })
cwd0 <- getwd()
invisible(tryCatch(with_abc_output_dir(file.path(od_base, "x"), stop("boom")),
                   error = function(e) NULL))
test("with_abc_output_dir() restores cwd even after an error",
     identical(normalizePath(getwd()), normalizePath(cwd0)))
setwd(cwd0)  # defensive: ensure we're back regardless of the result above

fmt_dir <- file.path(tempdir(), sprintf("smoke_fmt_%d", as.integer(Sys.time())))
dir.create(fmt_dir, showWarnings = FALSE)
for (st in c("20260101_090000", "20260101_093000", "20260102_010000")) {
  writeLines("x", file.path(fmt_dir, sprintf("fiber_ABC_SMC_Worst_WestAfrica_HCWrisk_%s.rds", st)))
}
test("new full-timestamp filenames sort chronologically via find_latest_file()",
     grepl("20260102_010000", basename(find_latest_file(fmt_dir, "Worst_WestAfrica"))))


# ----------------------------------------------------------------------------
section("9. simulation_helpers: bin_counts() / q_summary() / %||%")
# ----------------------------------------------------------------------------
test("bin_counts() bins day-indices into fixed-width bins",
     identical(bin_counts(c(0, 1, 7, 8, 13), bin_width = 7, n_bins = 2), c(2L, 3L)))
test("bin_counts() empty input -> zeros",
     identical(bin_counts(integer(0), 7, 3), integer(3)))
test("bin_counts() clips overflow into the last bin",
     identical(bin_counts(100, bin_width = 7, n_bins = 2), c(0L, 1L)))
qs <- q_summary(c(1, 2, 3, 4, 5), probs = c(0.5))
test("q_summary() mean is correct",   approx_eq(unname(qs["mean"]), 3))
test("q_summary() median is correct", approx_eq(unname(qs["q50"]), 3))
test("q_summary() drops non-finite values",
     approx_eq(unname(q_summary(c(1, 2, 3, NA, Inf), probs = c(0.5))["mean"]), 2))
test("%||% returns lhs when non-null", identical(("a" %||% "b"), "a"))
test("%||% returns rhs when null",     identical((NULL %||% "b"), "b"))


# ----------------------------------------------------------------------------
section("10. Live model interface + end-to-end (requires fiber)")
# ----------------------------------------------------------------------------
if (!has_fiber) {
  skip("check_model_function_version()",          "fiber not installed")
  skip("ppe_efficacy + ppe_coverage_hcw in branching_process_main() signature", "fiber not installed")
  skip("make_model_parameters() builds tv curves", "fiber not installed")
  skip("end-to-end branching_process_main() run",  "fiber not installed")
  skip("simulate_one() one replicate",             "fiber not installed")
  cat("\n  NOTE: this is the one part static checks cannot cover. To exercise it,\n")
  cat("        install fiber and re-run:\n")
  cat("          devtools::install_github(\"petal-code/fiber\")\n")
  cat("          Rscript tests/smoke_test.R\n")
} else {
  suppressWarnings(suppressMessages(library(fiber)))
  test("branching_process_main() is on the search path after library(fiber)",
       exists("branching_process_main", mode = "function"))
  test("check_model_function_version() passes (prop_etu/etu_efficacy/general_* in formals)",
       { check_model_function_version(); TRUE })
  test("ppe_efficacy + ppe_coverage_hcw accepted by branching_process_main()",
       all(c("ppe_efficacy", "ppe_coverage_hcw") %in% names(formals(branching_process_main))))

  test("make_model_parameters() builds time-varying curves for Worst_WestAfrica",
       { mpf <<- make_model_parameters(scenario_id = "Worst_WestAfrica",
                                       scenario_matrix = sm,
                                       overrides = list(check_final_size = 2000L,
                                                        seeding_cases = 5L));
         is.list(mpf) && !is.null(mpf$args$prop_etu) })
  test("solve_offspring_means_for_R0() runs on the real built args",
       { solf <<- solve_offspring_means_for_R0(R0 = 1.5, args = mpf$args,
                                               proportion_transmission_from_funerals = 0.3,
                                               n = 5000, seed = 1);
         is.finite(solf$D_direct_multiplier) && solf$D_direct_multiplier > 0 })
  build_one <- function() build_abc_model_args(
    R0 = 1.5, prop_funeral = 0.3, hcw_risk_scalar = 1.5,
    base = mpf$base_args, tv = mpf$tv_args,
    D = solf$D_direct_multiplier, F_fun = solf$F_funeral_multiplier,
    seeding_cases = 5
  )
  test("build_abc_model_args() + branching_process_main(): one run returns a tdf",
       { set.seed(1); out1 <- do.call(branching_process_main, build_one());
         is.list(out1) && !is.null(out1$tdf) && is.data.frame(out1$tdf) })
  test("abc_summarise() summarises a run",
       { set.seed(2); s1 <- abc_summarise(do.call(branching_process_main, build_one()));
         all(c("n_cases", "n_deaths", "n_hcw_deaths", "duration") %in% names(s1)) })
  test("simulate_one() runs one base-arm replicate",
       { r <- simulate_one(list(set_id = 1L, rep_id = 1L, arm = "base", seed = 1L),
                           list(build_one()), obv_cfg = list());
         all(c("n_cases", "n_deaths", "n_hcw_deaths", "death_days") %in% names(r)) })
}


# ----------------------------------------------------------------------------
section("Summary")
# ----------------------------------------------------------------------------
cat(sprintf("\n  PASS: %d    FAIL: %d    SKIP: %d\n", .smoke$pass, .smoke$fail, .smoke$skip))
if (.smoke$fail > 0L) {
  cat("\n  Failing checks:\n")
  for (d in .smoke$failed) cat(sprintf("    - %s\n", d))
  cat("\n  RESULT: ✗ SMOKE TEST FAILED\n\n")
  if (!interactive()) quit(status = 1L, save = "no")
} else {
  cat(sprintf("\n  RESULT: ✓ all executed checks passed%s\n\n",
              if (.smoke$skip > 0L) "  (some checks skipped — see [SKIP] lines)" else ""))
  if (!interactive()) quit(status = 0L, save = "no")
}
