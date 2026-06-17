# dose_estimation_subanalysis

Fits a logistic "Q curve" to a short series of dose-coverage observations,
projects it forward up to a year, turns it into time-varying NPI inputs, and
runs the `fiber` branching-process model across a grid of baseline R0 values.

## Scripts (run in order)

1. **`01_fit_dose_q_curve.R`** — fits a single logistic Q curve to the five
   `(date, coverage %)` observations using the same family of model as the
   previous Q-curve analyses (`stan-models/logistic_qcurve_single.stan`):

   ```
   q(day) = L + (U - L) * inv_logit( k * (day - t50) )
   ```

   * Percentages → proportions; calendar dates → **relative days** from a
     configurable `START_DATE` (day 0).
   * **Front padding:** set `START_DATE` earlier than 18 May 2026 to pad the
     front with daily zeros up to 18 May (encodes a longer 0% period).
   * **Very informative endpoint priors** pin the minimum at 0% (`L ≈ 0`) and
     the maximum at 100% (`U ≈ 1`); the data inform only the shape (`t50`, `k`).
   * The logistic is in raw day units (not renormalised to the window), so it
     **extrapolates** cleanly: the curve is projected forward `FORWARD_DAYS`
     (default 365) past the last observation at the estimated growth rate.
   * **Output:** `outputs/dose_q_curve.rds` is *literally the Q curve* (one row
     per day: mean + 90% band, calendar date, observed-vs-extrapolated flag).
     Also writes `dose_q_curve_fit.rds`, `dose_q_curve.csv`, and a diagnostic PNG.

2. **`02_npi_inputs_and_fiber_runs.R`** — uses the Q curve to drive the model.

   * Define a **min/max value for each NPI parameter** in `NPI_SPECS`
     (`q0` = value at Q = 0 / worst response, `q1` = value at Q = 1 / best
     response). Each parameter is mapped along the Q curve as
     `value(t) = q0 + (q1 - q0) * Q(t)`. Direction is set purely by `q0` vs `q1`
     (e.g. `delay_hosp` improves *downward*). `safe_funeral_prop` is converted
     to the model's unsafe-funeral probability as `1 - safe_funeral_prop`.
   * The resulting time-varying inputs are saved to
     `outputs/dose_npi_scenario_matrix.csv` / `.rds` (+ a tidy long CSV).
   * Those inputs are used as **fixed** model inputs while baseline R0 is swept
     over `R0_GRID` (1.30 → 1.60). For each R0, `N_STOCH` stochastic replicates
     are run in parallel (`future`), seeded with `SEEDING_CASES` (= 5) infections,
     **re-running any replicate that fails to reach `TAKEOFF_N` (= 100) infections**
     (advancing the seed, up to `MAX_RETRIES`).
   * For each R0, summarised **from the individual raw runs**:
       * median cumulative cases at each timepoint in `TIMEPOINTS`
         (day 10, 20, 30, …) → `outputs/dose_r0_grid_cumulative_cases.csv`;
       * median time to reach each case amount in `AMOUNTS`
         → `outputs/dose_r0_grid_time_to_amounts.csv`.
     Per-replicate metrics and a bundled results object are also saved.

## Configuration

Both scripts expose a clearly-marked configuration block near the top
(observations/start date/priors in script 1; NPI min–max values, R0 grid,
replicate count, takeoff threshold, timepoints and case amounts in script 2).

## Dependencies

`cmdstanr` (+ a working CmdStan toolchain) for script 1; `fiber`, `future`,
`future.apply` for script 2; plus `dplyr`, `tidyr`, `ggplot2`, `here`. Run from
anywhere — paths resolve from the repo root via `here::here()`.

## Note on `outputs/`

`outputs/` is generated at runtime and is not tracked in git (consistent with
the repo-wide `.gitignore`); the scripts create it on first run.
