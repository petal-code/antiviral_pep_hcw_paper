# Latent response-parameter estimation (West Africa + DRC, REVISED methodology)

Clean reimplementation of the filovirus Q-curve fitting for the **revised
methodology** scenarios. It is the sibling of
`analyses/01_latent_response_parameter_estimation/` (the original methodology):
the pipeline shape, the data preparation, and the combine/consolidation code are
deliberately kept as close as possible, so the two folders can be compared
script-for-script. Every file here has `(revisedMethod)` appended to its name.

The six response parameters (fixed canonical order; see `helpers(revisedMethod).R`):
`delay_hosp`, `p_hosp`, `p_ETU`, `latent_IPC`, `p_unsafe_funeral_comm`,
`p_unsafe_funeral_hosp` (internal names; the published column names `prob_hosp`,
`prop_etu`, `ipc_helper`, … are applied only in `03`).

## What "revised methodology" means (the two differences from the original)

1. **Endpoint-constrained fitting (Model C).** In the original methodology each
   parameter's two magnitude endpoints are *estimated* (Model A estimates the
   shape **and** the endpoints; Model B fixes the shared `Q` and estimates the
   endpoints). In the revised methodology the endpoints are instead **locked** to
   early/late literature-window extrema (`lock_endpoints()` in `helpers`), and:
   - **West Africa** estimates only the curve **shape** (`t50`, `k`, partially
     pooled), with the locked endpoints supplied as data — this is **Model C**
     (`stan-models/modelC_endpointConstrained_estimateShape(revisedMethod).stan`).
     There are **no "tweak" priors** (nothing to tweak), hence no
     with/without-tweaks comparison and no `notweaks` scenario.
   - **DRC** does not fit anything: the fixed empirical conflict `Q` is mapped
     **deterministically** onto the locked endpoints,
     `θ(t) = start + (end − start)·Q(t)` (a smooth fit would erase the
     conflict-interrupted shape, the same reason Model B never smoothed `Q`).
2. **q-scaled IPC for DRC.** Under the revised methodology `ipc_helper` for the
   DRC conflict scenarios is the q-scaled endpoint mapping
   `ipc_helper(t) = ipc_low + (ipc_high − ipc_low)·q_value(t)` with
   `0.071 → 0.746` — exactly `latent_IPC`'s DRC summary range in the
   parameter-table workbook. `02` sets `latent_IPC`'s DRC endpoints to that full
   range, so mapping the shared `Q` reproduces the rule (and the `++` collapse
   drives it to `ipc_low`). (Under the original methodology `ipc_helper` is the
   Model B fitted `latent_IPC` instead.)

## Pipeline (run in order)

All paths are resolved from the repository root with `here::here()`, so the
scripts run regardless of the working directory. Source in order:

```r
d <- here::here("analyses", "01_latent_response_parameter_estimation_revisedMethodology")
source(file.path(d, "00_DataPreparation_and_Cleaning(revisedMethod).R")) # data-raw -> data-processed
source(file.path(d, "01_WestAfrica_QCurve_Fitting(revisedMethod).R"))    # Model C -> WestAfrica_QCurve_Fit(revisedMethod).rds
source(file.path(d, "02_DRC_QCurve_Fitting(revisedMethod).R"))           # endpoint mapping -> DRC_QCurve_*_Fit(revisedMethod).rds
source(file.path(d, "03_Combine_QCurves(revisedMethod).R"))              # -> combined_revised_methodology_outputs.csv
```

`03` produces **three** scenarios on a common 0–730 day grid: `worst_west_africa`,
`drc_conflict`, and `drc_conflict_plusplus` (the latter is the same mapping with
the success→0 collapse over days 200–300 already baked into the `Q` series by
`00`).

## Folders

Inputs and outputs live at the **repository top level** (shared across analyses),
not inside this folder:

- `obv_hcw_paper/data-raw/` — the source workbooks. The revised methodology reads
  its anchors and literature ranges from the **parameter-table workbook**
  (`filovirus_model_parameter_table_4_scenarios_with_etu_baseline.xlsx`, sheets
  "Worst West Africa" / "DRC conflict-smoothed") — the same workbook the
  github_upload revised scripts use — plus the shared DRC SDB line-list. (The
  original methodology reads a *different* anchor workbook, so the two
  methodologies' locked endpoints differ.)
- `obv_hcw_paper/data-processed/` — cleaned inputs, fitted `.rds`, and the final
  CSV. The revised pipeline writes into its **own** subfolders
  (`WestAfrica_QCurve_revisedMethod/`, `DRC_QCurve_revisedMethod/`) and
  `combined_revised_methodology_outputs.csv`, so it never clobbers the
  original-methodology outputs.

Inside this analysis folder:

- `stan-models/` — the single Stan file (Model C).
- `helpers(revisedMethod).R` — the shared utilities and `here::here()` path
  constants, plus the two revised-methodology additions: `q_norm()` (the
  normalised logistic, the R twin of the Stan function) and `lock_endpoints()`
  (the endpoint-locking rule that defines the methodology).

## Notes

- Requires `cmdstanr` + a working CmdStan toolchain (for `01` only — `02` is pure
  R), and `readxl`, `dplyr`, `tidyr`, `readr`, `tibble`, `ggplot2`.
- `00`'s SDB reconstruction (PART 2) is identical to the original-methodology
  `00`; its anchor reader (PART 1) differs because the revised methodology reads
  the parameter-table workbook's layout (summary ranges + description-encoded
  anchors) rather than the original's flat anchor sheet.
- The DRC no-conflict scenario is prepared by `00` but, as in the original
  methodology, is **not** built into the combined output.
