## =============================================================================
## test-calculate_model_approx_rt.R
##
## Tests for the analytical time-resolved reproduction-number curves
## (functions/calculate_model_approx_rt.R).
##
## Runs EITHER under testthat (e.g. when dropped into an R package) OR standalone
## with plain `Rscript tests/testthat/test-calculate_model_approx_rt.R` (a tiny
## base-R shim below stands in for testthat when the package is not installed).
## No `fiber` is required: time-varying curves are mimicked with approxfun +
## attr(, "times"), exactly the interface make_time_varying() exposes.
## =============================================================================

## ---- testthat, or a minimal shim so the file runs without it -----------------
if (requireNamespace("testthat", quietly = TRUE)) {
  library(testthat)
} else {
  .fail <- function(msg) stop(msg, call. = FALSE)
  test_that    <- function(desc, code) { cat("TEST:", desc, "... "); force(code); cat("PASS\n") }
  expect_true  <- function(x, info = NULL) if (!isTRUE(x)) .fail(info %||% "expected TRUE")
  expect_equal <- function(object, expected, tolerance = 1e-8, info = NULL)
    if (!all(abs(object - expected) <= tolerance)) .fail(info %||% "values not within tolerance")
  expect_lt    <- function(object, expected) if (!(object < expected)) .fail("expected object < expected")
  expect_gt    <- function(object, expected) if (!(object > expected)) .fail("expected object > expected")
}
`%||%` <- function(a, b) if (is.null(a)) b else a

## ---- locate + source the functions under test (standalone use) --------------
if (!exists("Rt_curve_single_type")) {
  root <- if (file.exists("functions/calculate_model_approx_rt.R")) "." else
          if (file.exists(file.path("..", "..", "functions",
                                    "calculate_model_approx_rt.R"))) file.path("..", "..") else
          stop("Cannot locate functions/ from getwd() = ", getwd())
  source(file.path(root, "functions", "setup_model_parameters.R"))
  source(file.path(root, "functions", "calculate_model_approx_r0.R"))
  source(file.path(root, "functions", "calculate_model_approx_rt.R"))
}

## ---- helpers ----------------------------------------------------------------
## Mimic fiber::make_time_varying(): a linear approxfun that clamps outside its
## range (rule = 2) and exposes its breakpoints via attr(, "times").
mk_curve <- function(times, values) {
  f <- stats::approxfun(times, values, method = "linear", rule = 2)
  attr(f, "times") <- times
  f
}

## Build a model-args list (scalars from make_base_args() + chosen curves), with
## the two offspring means the R0 inversion would normally set.
build_args <- function(curves) {
  args <- c(make_base_args(), curves)
  args$mn_offspring_genPop  <- 1.25
  args$mn_offspring_funeral <- 0.25
  args
}

FLAT_SCALARS <- list(
  prob_hospitalised_genPop     = 0.45,
  hospitalisation_delay_factor = 3.0,
  prop_etu                     = 0.40,
  p_unsafe_funeral_comm_genPop = 0.50,
  p_unsafe_funeral_hosp_genPop = 0.08
)


## =============================================================================
## Test 1 -- degeneracy: with flat (scalar) curves, R_inst(t) = R_case(t) = const
## =============================================================================
test_that("flat curves: R_inst == R_case, constant across the grid, == R0", {
  args <- build_args(FLAT_SCALARS)
  grid <- c(0, 10, 50, 100)
  rt <- Rt_curve_single_type(args, times = grid, n = 50000, seed = 7)

  ## R_inst and R_case coincide at every t (same variates, clock is irrelevant
  ## when every curve is a constant) -- exact up to floating point.
  expect_equal(rt$R_inst, rt$R_case, tolerance = 1e-9,
               info = "R_inst and R_case must coincide under flat curves")

  ## both curves are flat across the grid.
  expect_equal(rt$R_inst, rep(rt$R_inst[1], length(grid)), tolerance = 1e-9,
               info = "R_inst must be constant under flat curves")
  expect_equal(rt$R_case, rep(rt$R_case[1], length(grid)), tolerance = 1e-9,
               info = "R_case must be constant under flat curves")

  ## and they match the independent R0 approximation (MC tolerance).
  r0 <- R0_single_type_from_args(args, n = 100000, seed = 11)$R0
  expect_true(abs(rt$R_inst[1] - r0) < 0.01,
              info = sprintf("flat R_inst (%.4f) ~ R0 (%.4f)", rt$R_inst[1], r0))
})


## =============================================================================
## Test 2 -- R_inst(0) reduces to R0_single_type_from_args()$R0 (MC tolerance)
## =============================================================================
test_that("R_inst(0) matches the t=0 R0 approximation under moving curves", {
  tt <- c(0, 50, 100, 200)
  args <- build_args(list(
    prob_hospitalised_genPop     = mk_curve(tt, c(0.20, 0.40, 0.60, 0.70)),
    hospitalisation_delay_factor = mk_curve(tt, c(6.0, 4.0, 3.0, 2.0)),
    prop_etu                     = mk_curve(tt, c(0.00, 0.20, 0.50, 0.75)),
    p_unsafe_funeral_comm_genPop = mk_curve(tt, c(0.90, 0.60, 0.40, 0.30)),
    p_unsafe_funeral_hosp_genPop = mk_curve(tt, c(0.20, 0.10, 0.05, 0.02))
  ))
  r0  <- R0_single_type_from_args(args, n = 100000, seed = 42)$R0
  ri0 <- Rt_curve_single_type(args, times = 0, n = 100000, seed = 42)$R_inst
  ## Bit-exactness is NOT expected (CRN changes the draw order); tolerance is.
  expect_true(abs(ri0 - r0) < 0.01,
              info = sprintf("R_inst(0)=%.5f vs R0=%.5f (absdiff=%.5f)",
                             ri0, r0, abs(ri0 - r0)))
})


## =============================================================================
## Test 3 -- moving curves: R_case(0) < R_inst(0) during the improving ramp
## =============================================================================
test_that("fast intervention ramp over a GI makes R_case(0) < R_inst(0)", {
  ## conditions improve sharply over the first ~generation interval, then flat.
  tt <- c(0, 15, 400)
  args <- build_args(list(
    prob_hospitalised_genPop     = mk_curve(tt, c(0.20, 0.80, 0.80)),
    hospitalisation_delay_factor = mk_curve(tt, c(7.0, 2.0, 2.0)),
    prop_etu                     = mk_curve(tt, c(0.00, 0.90, 0.90)),
    p_unsafe_funeral_comm_genPop = mk_curve(tt, c(0.90, 0.20, 0.20)),
    p_unsafe_funeral_hosp_genPop = mk_curve(tt, c(0.30, 0.05, 0.05))
  ))
  rt0 <- Rt_curve_single_type(args, times = 0, n = 100000, seed = 3)
  ## forward-looking cohort already benefits from the imminent improvement.
  expect_lt(rt0$R_case, rt0$R_inst)
  expect_gt(rt0$R_inst - rt0$R_case, 0.01)

  ## and once everything is flat (well past the ramp), the two coincide again.
  rt_late <- Rt_curve_single_type(args, times = 380, n = 100000, seed = 3)
  expect_true(abs(rt_late$R_case - rt_late$R_inst) < 1e-6,
              info = "curves flat past the ramp -> R_case == R_inst")
})

if (!requireNamespace("testthat", quietly = TRUE)) cat("\nAll standalone tests passed.\n")
