# Filovirus Q-curve fitting analysis

This folder contains the Q-curve fitting scripts, raw/clean input data, and final reference outputs for the 14 filovirus response scenarios.

## Current pipeline

```r
source('R/pipeline/00_run_all_14_generate_csvs_and_print_native_plots.R')
source('R/pipeline/08_patch_revised_drc_ipc_qscaled_final.R')
source('R/pipeline/12_patch_WA_conflict_plusplus_direction_final.R')
source('R/pipeline/13_print_final_14_WAplusplusFixed_plots_v2.R')
```

## Candidate latest West Africa conflict++ timing update

```r
source('R/pipeline/14_retime_WA_conflict_plusplus_to_WA_conflict_window.R')
```

## Notes

- Scenario-level Q-curve scripts are in `R/scenarios/`.
- Input workbooks and DRC SDB data are in `data/inputs/`.
- Final reference matrices are in `outputs/final/`.
- Revised DRC conflict IPC/PPE uses: `ipc_helper = 0.071 + (0.746 - 0.071) * q_value`.
- West Africa conflict++ is hypothetical and applies a DRC++-style severe response collapse to the West Africa conflict timing.
