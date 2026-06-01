# obv_hcw_paper

Code and analysis for [PAPER TITLE].

## Repo structure

- `functions/` — shared model code used across analyses. Parameter setup
  (`setup_model_parameters.R`) and the approximate-R0 solver
  (`calculate_model_approx_r0.R`) are shared by every analysis and target the
  current fiber NPI interface. ABC calibration is split so two fitting schemes
  can be run side by side without duplicated plumbing:
    - `abc_calibration_functions_common.R` — generic ABC helpers (summaries,
      output dirs, disk readers) used by both schemes.
    - `abc_calibration_functions_hcwRisk.R` — the HCW-risk scheme (fits a
      `hcw_risk_scalar`).
    - `abc_calibration_functions_npi.R` — the NPI-efficacy scheme (fits a single
      `npi_scaler` over PPE/ETU conditional efficacies).
  Plus generic helpers like `io_helpers.R` (`find_latest_file()`). Analyses
  source these files from here rather than keeping their own copies.
- `data-raw/` — raw, read-only input data.
- `data-processed/` — processed data, regenerable from `data-raw/`.
- `analyses/` — analysis scripts. Each subfolder is one data-processing step or
  figure; the `NN_*` folders are numbered in run order.
- `outputs/` — curated, manuscript-ready artefacts (figures, tables, the single
  final fit per scenario), **tracked** in git and organised one subfolder per
  analysis, e.g. `outputs/02_ABC_model_fits_HCWrisk/`.
- `manuscript/` — paper source.

## Running the analysis

1. Place raw data in `data-raw/`.
2. Run the scripts in each `analyses/` subfolder (numbered folders in order).
3. Outputs land in `outputs/<analysis_name>/` and processed data in
   `data-processed/`.

## Conventions within an `analyses/` subfolder

- Shared model code lives in the repo-level `functions/` folder and is sourced
  from there — not duplicated per analysis.
- `<name>_helper_functions.R` (or `helper_functions_<fig>.R`) — helpers
  specific to that one analysis.
- `01_*.R`, `02_*.R`, … — scripts numbered in run order; the heavy compute step
  saves an intermediate `.rds` that the later plotting step reads, and writes
  finished figures/tables to `outputs/<analysis_name>/`.
- Resolve every path with `here::here()` so scripts run the same regardless of
  the working directory they are launched from.
