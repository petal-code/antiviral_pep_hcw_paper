# 05_SI_leaky_onward — "leaky obeldesivir" sensitivity (reviewer response)

Addresses the reviewer comment that assuming obeldesivir (OBV) **fully** stops a
treated health worker's (HCW) onward transmission is optimistic. We ask what a
partially-effective ("leaky") drug does to OBV's estimated **death** impact — the
paper's main outcome, focused on **HCW deaths** — for the DRC and West Africa
archetypes, using the same posterior draws and the same 80%-coverage /
80%-efficacy programme as the main figures (the `cov80_obv80` arm).

## What it does

For each posterior draw (200, weighted) × replicate (10), per archetype, it runs
the WITH-OBV base model once, takes its `prevented_completed` set of averted
infections, and forward-simulates them with `fiber::estimate_leaky_onward()`
(requires the `claude/antiviral-efficacy-obeldesivir-e4tewv` fiber branch, which
`01_analysis_leaky_onward.R` installs):

- **Analysis 1 — the current reporting is conservative.** The paper credits OBV
  only with the *directly* averted (index) infections/deaths. Forward-simulating
  the averted lineages with no drug (the no-OBV counterfactual) shows the true
  averted burden is index **+ downstream**, i.e. larger. So the current number is
  a lower bound.
- **Analysis 2 — deaths are robust to leakiness.** Re-simulating the averted
  lineages under a leaky drug (residual transmissibility `r = 0 … 1`), while
  keeping OBV's protection against **death** on everyone it treats, shows how
  many deaths still occur. Deaths averted at leakiness `r` = no-OBV deaths −
  leaked deaths. Because the drug keeps treating the downstream HCWs it reaches
  (same coverage), the HCW-death benefit holds up across the whole range and
  never falls below what the paper reports.

Death protection on the treated is applied as a **counting** rule on the
simulated tree (a case counts as a death only if it died *and* was not
treated-effectively). See `leaky_onward_helpers.R` for the full rationale; the
method's validity rests on fiber's current no-depletion regime.

## Files

| file | role |
|---|---|
| `leaky_onward_helpers.R` | pure, unit-tested death-accounting functions |
| `01_analysis_leaky_onward.R` | installs fiber branch, runs base sims + forward sims, saves `_intermediate/leaky_onward_per_run.rds` |
| `02_plot_leaky_onward.R` | aggregates (median over reps, then across draws) and writes tables + figures |

## Run

```r
# heavy: 200 x 10 x 2 base runs + 12 forward sims each. Smoke-test first:
#   set QUICK_TEST <- TRUE in 01_analysis_leaky_onward.R
Rscript analyses/05_SI_leaky_onward/01_analysis_leaky_onward.R
Rscript analyses/05_SI_leaky_onward/02_plot_leaky_onward.R
```

Outputs: `figures/fig_leaky_onward_hcw_deaths.*` (headline), `…_accruing.*`
(companion), and `output_figgen/leaky_onward_summary.csv` +
`…_analysis1_table.csv`. The forward sims grow near `r = 1`; they are capped at
each archetype's `check_final_size` and warn rather than truncating silently.
