# obv_hcw_paper

Code and analysis for [PAPER TITLE].

## Repo structure

- `R/` — reusable functions used across analyses
- `data/` — raw, read-only input data
- `data-derived/` — processed data, regenerable from `data/`
- `analyses/` — analysis scripts, numbered in run order. Each subfolder is one data-processing step or figure.
- `outputs/` — figures, tables, manuscript-ready artefacts
- `manuscript/` — paper source

## Running the analysis

1. Place raw data in `data/`
2. Run scripts in `analyses/` in numerical order
3. Outputs land in `outputs/` and `data-derived/`

## Conventions within each `analyses/XX_*/` folder

- `helper_functions.R` — analysis-specific helpers
- `01_compute.R` — analysis, saves intermediate `.rds`
- `02_plot.R` — reads intermediate, produces figure
