# obv_hcw_paper

Code and analysis for [PAPER TITLE].

## Repo structure

- `functions/` — shared model code used across analyses (parameter setup, the
  approximate-R0 solver, and the ABC parameter→model-args mapping), plus generic
  helpers like `io_helpers.R` (e.g. `find_latest_file()` to pull the most recent
  timestamped artefact out of a folder). Analyses source these files from here
  rather than keeping their own copies.
- `data-raw/` — raw, read-only input data.
- `data-processed/` — processed data, regenerable from `data-raw/`.
- `analyses/` — analysis scripts. Each subfolder is one data-processing step or
  figure; the `NN_*` folders are numbered in run order.
- `outputs/` — curated, manuscript-ready artefacts (figures, tables, the single
  final fit per scenario), **tracked** in git and organised one subfolder per
  analysis, e.g. `outputs/02_ABC_model_fits_HCWrisk/`.
- `manuscript/` — paper source.

Heavy, regenerable **intermediates** (per-step ABC particle clouds, the
compute→plot handoff `.rds`) are *not* tracked: they live next to the analysis
that makes them — in `abc_outputs/` or a `_intermediate/` folder — both
gitignored. Only the curated finals in `outputs/` are committed.

## Running the analysis

1. Place raw data in `data-raw/`.
2. Run the scripts in each `analyses/` subfolder (numbered folders in order).
3. Curated outputs land in `outputs/<analysis_name>/` (and processed data in
   `data-processed/`); heavy regenerable intermediates stay next to each
   analysis in a gitignored folder.

## Conventions within an `analyses/` subfolder

- Shared model code lives in the repo-level `functions/` folder and is sourced
  from there — not duplicated per analysis.
- `<name>_helper_functions.R` (or `helper_functions_<fig>.R`) — helpers
  specific to that one analysis.
- `01_*.R`, `02_*.R`, … — scripts numbered in run order; the heavy compute step
  saves an intermediate `.rds` (to a gitignored `_intermediate/` folder) that
  the later plotting step reads, and writes finished figures/tables to
  `outputs/<analysis_name>/`.
- Resolve every path with `here::here()` so scripts run the same regardless of
  the working directory they are launched from.
