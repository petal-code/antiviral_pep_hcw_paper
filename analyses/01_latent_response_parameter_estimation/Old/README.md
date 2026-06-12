# Latent response-parameter estimation (West Africa + DRC, original methodology)

Clean reimplementation of the filovirus Q-curve fitting for the **original
methodology** scenarios. It estimates how a set of latent outbreak-response
parameters evolve over time and writes a combined scenario matrix.

The six response parameters (fixed canonical order; see `helpers.R`):
`delay_hosp`, `p_hosp`, `p_ETU`, `latent_IPC`, `p_unsafe_funeral_comm`,
`p_unsafe_funeral_hosp`.

## Two Bayesian models (`stan-models/`)

- **Model A** (`modelA_partialpool_estimateQ_*.stan`) — estimates BOTH the
  response-quality curve shape `Q_j(tau)` (a finite-window logistic, partially
  pooled across parameters on `t50` and `k`) AND the per-parameter magnitude
  endpoints. Used for **West Africa**. Provided in two forms, **with** and
  **without** the targeted "tweak" priors, so their effect can be compared.
- **Model B** (`modelB_fixedQ_boundsOnly.stan`) — the shared `Q` is supplied as
  fixed DATA (the empirical SDB success curve); only the per-parameter endpoints
  are estimated (no shape estimation, no pooling). Used for **DRC conflict** and
  **DRC conflict++**.

## Pipeline (run in order)

All paths are resolved from the repository root with `here::here()` (it locates
`obv_hcw_paper.Rproj` / `.git`), so the scripts run regardless of the working
directory. Open `obv_hcw_paper.Rproj` (or set the working directory anywhere in
the repo) and source in order:

```r
d <- here::here("analyses", "01_latent_response_parameter_estimation")
source(file.path(d, "00_DataPreparation_and_Cleaning.R"))      # data-raw -> data-processed
source(file.path(d, "01_WestAfrica_QCurve_Fitting_Original.R")) # Model A (+tweaks) -> wa_fit.rds
source(file.path(d, "02_DRC_QCurve_Fitting_Original.R"))        # Model B -> drc_conflict(_plusplus)_fit.rds
source(file.path(d, "03_Combine_QCurves.R"))                    # -> combined_original_methodology_outputs.csv
```

`03` produces **five** scenarios on a common 0–730 day grid: `worst_west_africa`,
`drc_conflict`, `drc_conflict_plusplus`, `worst_west_africa_conflict` (built by
modulating the West Africa curves with the DRC conflict curve and stretching the
timeline), and `worst_west_africa_conflict_plusplus` (a forced response collapse
over days 200–300).

## Diagnostic / comparison scripts (not part of the main output)

- `west_africa_checking.R` — fits Model A with and without the tweaks and
  overlays the curves, so "what the tweaks do" is visible per parameter.
- `DRC_no_conflict_checking.R` — compares Model B against Model A (using all
  data, with the dense SDB series at full weight vs. down-weighted), to see how
  much the SDB series dominates the partially-pooled fit. The DRC no-conflict
  scenario is intentionally **not** in the main pipeline pending this decision.

## Folders

Inputs and outputs live at the **repository top level** (shared across analyses,
matching the rest of `obv_hcw_paper`), not inside this analysis folder:

- `obv_hcw_paper/data-raw/` — the two source workbooks (curve anchors; DRC SDB
  line-list).
- `obv_hcw_paper/data-processed/` — cleaned inputs, fitted `.rds`, and the final
  CSV produced by this pipeline (alongside the other analyses' processed data).

Inside this analysis folder:

- `stan-models/` — the three Stan files.
- `helpers.R` — small shared utilities and the `here::here()` path constants
  (`DIR_RAW`, `DIR_PROCESSED`, `DIR_STAN`); one definition each, no duplication.

## Notes

- Requires `cmdstanr` + a working CmdStan toolchain, and `readxl`, `dplyr`,
  `tidyr`, `readr`, `tibble`, `ggplot2`.
- Under the **original** methodology `ipc_helper` is the fitted `latent_IPC`
  parameter; there is no separate q-scaling patch (that is a revised-methodology
  step). The conflict++ collapse is built into the SDB Q series / the combine
  step, not applied as a post-hoc patch.
- The DRC no-conflict scenario's only role in the main pipeline is its horizon
  (used to set the West-Africa-with-conflict time stretch).
