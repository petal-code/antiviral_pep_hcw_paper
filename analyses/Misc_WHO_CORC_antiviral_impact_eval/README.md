# 04_antiviral_impact

Estimate the impact of antiviral post-exposure prophylaxis (PEP) on a
DRC-like Ebola outbreak, using the posterior from the ABC-SMC calibration in
`analyses/02_model_fits`.

## What it does

1. **Posterior** — reads the 3 fitted parameters (`R0`, `prop_funeral`,
   `hcw_risk_scalar`) from the final completed ABC step (step 7) of the latest
   `Middle_DRC_ConflictSmoothed` run under
   `analyses/02_model_fits/abc_outputs/`.
2. **Parameter conversion** — converts the 3 fitted parameters into the 4 fiber
   model parameters (`mn_offspring_genPop`, `mn_offspring_funeral`,
   `prob_hcw_cond_genPop_hospital`, `prob_hcw_cond_hcw_hospital`) using the same
   mapping as the calibration (`build_abc_model_args()`).
3. **Downsample** — weighted resample to **100 parameter sets**, with **5
   stochastic replicates** each.
4. **Simulate** — runs fiber **with** and **without** antiviral (80% efficacy,
   100% coverage, 100% adherence; modelled as PEP for HCWs exposed in hospital),
   parallelised with `future` (`multisession`, Windows-compatible; defaults to
   10 workers).
5. **Outputs** — total deaths and total HCW deaths per arm; deaths averted
   (total and HCW) by antiviral; epidemic-curve trajectories.

## Files

- `WHO_CORC_helper_functions.R` — analysis-specific helpers (posterior IO, parameter
  conversion, the per-replicate simulator run on workers, binning/summaries).
- `01_WHO_CORC_run_simulations.R` — runs the simulations and saves
  `WHO_CORC_prelim_antiviral_simulation_results.rds` plus summary CSVs in `outputs/`.
- `02_WHO_CORC_plot.R` — reads the intermediate and writes the figures to `outputs/`.

## Running

```r
# From the repo root or this folder:
source("analyses/Misc_WHO_CORC_antiviral_impact_eval/01_WHO_CORC_run_simulations.R")  # heavy compute
source("analyses/Misc_WHO_CORC_antiviral_impact_eval/02_WHO_CORC_plot.R")             # figures
```

Key knobs live in the CONFIG block at the top of `01_run_simulations.R`
(`N_SETS`, `N_REPS`, `N_WORKERS`, `OBV_CFG`, `BIN_WIDTH_DAYS`, seeds).

## Figures (in `outputs/`)

- `obeldesivir_deaths_over_time_individual.png` — every parameter set's mean
  deaths/week trajectory, with vs without obeldesivir.
- `obeldesivir_deaths_over_time_median_iqr.png` — median + 25–75% band.
- `obeldesivir_hcw_deaths_over_time_individual.png` — as above, HCW deaths.
- `obeldesivir_hcw_deaths_over_time_median_iqr.png` — median + 25–75% band.
- `obeldesivir_epidemic_curves_combined.png` — 2×2 combined panel.
- `obeldesivir_pct_hcw_deaths_averted_bar.png` — % of HCW deaths averted (with %
  of all deaths for context).

## Note on the paired comparison

Obeldesivir runs consume extra random draws, so a with/without pair seeded
identically diverges once the first obeldesivir draw fires (see the
reproducibility caveat in `branching_process_main()`). We reuse the seed across
arms for each `(set, rep)` pair to share early dynamics (variance reduction) and
summarise over many replicates and parameter sets rather than relying on
exact trajectory pairing.
