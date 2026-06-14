# 04_SI_sensitivity_DRC — DRC-like sensitivity analyses (SI)

Two supplementary sensitivity analyses that stress-test the DRC-like ("Middle",
PlusPlus) conclusions against the possibility that vaccination in the observed
North Kivu / Ituri outbreak made the calibrated archetype conservative.

Both reuse the **exact** posterior, scenario inputs and model machinery that
generated the main figures (see `analyses/03_figure_template/01_analysis_figure2.R`);
the only change is a per-particle multiplicative scaling of one calibrated
quantity, injected at the `build_abc_model_args_decoupled()` step.

## The two analyses

1. **Vaccine-free stress test (transmissibility ↑).** Baseline transmissibility
   is scaled up by **+10 % / +20 % / +30 %** via `R0`. Because
   `solve_offspring_means()` gives `mn_offspring_genPop = (1-prop_funeral)·R0/D`
   and `mn_offspring_funeral = prop_funeral·R0/F` with `D`, `F` independent of
   `R0`, multiplying `R0` by *f* multiplies **both offspring means** by exactly
   *f* — i.e. scaling `R0` ≡ scaling the offspring means. PEP scenarios are held
   fixed.
   *Reports:* no-PEP baseline HCW deaths, HCW deaths averted, % reduction.

2. **HCW-exposure upshift (`hcw_risk_scalar` ↑).** The HCW exposure scalar is
   scaled up by **+10 % / +25 % / +50 %**. It enters as
   `prob_hcw_cond_*_hospital = pmin(hcw_base_prob · hcw_risk_scalar · f, 1)`.
   *Reports:* baseline HCW deaths, HCW deaths averted, HCW-days lost averted,
   % reduction.

The as-fitted DRC archetype (×1.00) is simulated once and used as the shared
reference level in both analyses.

## PEP scenarios (arms)

The full-coverage antiviral-efficacy sweep from Figure 2 (100 % coverage,
50–90 % efficacy; `full_obv80` is the headline arm). The no-PEP baseline is the
matched counterfactual reconstructed as `tdf + prevented_completed` (same
convention as the figures). Extend `ARMS` / `COVERAGE_FNS` in
`01_run_SI_sensitivity_DRC.R` to add the Figure-3 ramp-coverage scenarios.

## Running

```r
# 1. Heavy compute (needs the `fiber` package). ~7 scalings × 5 arms ×
#    200 particles × 10 reps = 70,000 branching-process simulations.
source("analyses/04_SI_sensitivity_DRC/01_run_SI_sensitivity_DRC.R")

# 2. Tables + SI figures (no `fiber` needed).
source("analyses/04_SI_sensitivity_DRC/02_summarise_SI_sensitivity_DRC.R")
```

For a fast plumbing check, set `N_PARTICLES <- 20L; N_REPS <- 2L` near the top of
script 1 (commented override provided).

## Outputs

| File | Tracked | Contents |
|---|---|---|
| `outputs/04_SI_sensitivity_DRC/SI_sensitivity_DRC_run_summary.csv` | no (gitignored) | run-level: one row per scaling × arm × particle × rep |
| `output_figgen/SI_sensitivity_DRC_particle_df.csv` | yes | per-particle means (the analysable unit) |
| `output_figgen/SI_sensitivity_DRC_hcw_saturation.csv` | yes | fraction of particles whose scaled `prob_hcw` hits the `pmin(·,1)` cap |
| `output_figgen/SI_sensitivity_DRC_summary_long.csv` | yes | tidy median + 95 % CrI per metric/level/efficacy |
| `output_figgen/SI_sensitivity_DRC_summary_table.csv` | yes | formatted SI table (median [95 % CrI]) |
| `figures/figure_S_DRC_transmissibility_stress_test.{pdf,png}` | yes | SI figure 1 |
| `figures/figure_S_DRC_hcw_exposure_upshift.{pdf,png}` | yes | SI figure 2 |

Posterior uncertainty is summarised as **median + 95 % credible interval**
(2.5/97.5 % quantiles across particles), matching the main figures.

## Notes / caveats

- **`hcw_risk_scalar` cap.** With the fitted DRC posterior
  (`hcw_risk_scalar ∈ [1.25, 2.93]`) and `hcw_base_prob = 0.25`, no particles
  saturate the cap at ×1.10 or ×1.25, and only ~5 % do at ×1.50. The saturation
  CSV records this so the +50 % level can be interpreted honestly.
- **`check_final_size` is held at the main-analysis value (10 000)** so only
  transmissibility/exposure changes between cells. If the stressed epidemics
  saturate this ceiling, raise `DRC$check_final_size` (and note it).
- **Scenario CSV.** Uses `final_six_scenario_values_original_approach.csv` with
  `id = "Middle_DRC_ConflictSmoothed_PlusPlus"` — the same inputs the posterior
  was fitted against (the *revised*-methodology CSV uses different, lowercase
  scenario IDs and is **not** interchangeable with this posterior).

## What to look for

If absolute HCW deaths averted rise with the stress level while the
**proportional** reductions and the efficacy gradient stay stable, the main
conclusions are robust — and likely conservative — for vaccine-unavailable,
DRC-like outbreaks, and do not hinge on the exact fitted HCW-exposure level.
