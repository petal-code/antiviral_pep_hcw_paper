# obv_hcw_paper

Code and analysis for [PAPER TITLE].

## Repo structure

- `functions/` — shared model code used across analyses (parameter setup, the
  approximate-R0 solver, and the ABC parameter→model-args mapping). Analyses
  source these files from here rather than keeping their own copies.
- `data-raw/` — raw, read-only input data.
- `data-processed/` — processed data, regenerable from `data-raw/`.
- `analyses/` — analysis scripts. Each subfolder is one data-processing step or
  figure; the `NN_*` folders are numbered in run order.
- `outputs/` — figures, tables, manuscript-ready artefacts.
- `manuscript/` — paper source.

## Running the analysis

1. Place raw data in `data-raw/`.
2. Run the scripts in each `analyses/` subfolder (numbered folders in order).
3. Outputs land in `outputs/` and `data-processed/`.

## Conventions within an `analyses/` subfolder

- Shared model code lives in the repo-level `functions/` folder and is sourced
  from there — not duplicated per analysis.
- `<name>_helper_functions.R` (or `helper_functions_<fig>.R`) — helpers
  specific to that one analysis.
- `01_*.R`, `02_*.R`, … — scripts numbered in run order; the heavy compute step
  saves an intermediate `.rds` that the later plotting step reads.
