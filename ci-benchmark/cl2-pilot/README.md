# CL2 Pilot (Low n)

Small exploratory study for whether a CL2-style cluster correction is likely
to improve CI calibration for `pffr` at low numbers of curves.

## What it runs

- DGP: gaussian, `ff + linear`
- Sample sizes: `n = 20, 40`
- Correlation: `iid`, `ar1(0.9)`
- Methods:
  - `pffr` (default)
  - `pffr_hc` (observation-level HC)
  - `pffr_sandwich` (current cluster CR1-style)
  - `pffr_cl2` (experimental CL2-style adjustment in this pilot)

## Script

```bash
Rscript ci-benchmark/cl2-pilot/sim-study-cl2-pilot.R [mode]
```

Modes:

- `smoke`: 2 reps per DGP
- `pilot`: 8 reps per DGP
- `full`: 20 reps per DGP
- or pass a positive integer

## Non-Gaussian script

```bash
Rscript ci-benchmark/cl2-pilot/sim-study-cl2-pilot-nongaussian.R [mode]
```

Modes:

- `smoke`: 1 rep per DGP
- `pilot`: 4 reps per DGP
- `full`: 10 reps per DGP
- or pass a positive integer

## Outputs

- `ci-benchmark/cl2-pilot/results/`: per-cell incremental RDS files
- `ci-benchmark/cl2-pilot/results/results_combined.rds`
- `ci-benchmark/cl2-pilot/summary.rds`
- `ci-benchmark/cl2-pilot/summary_by_cell.csv`
- `ci-benchmark/cl2-pilot/summary_cluster_vs_cl2.csv`
- `ci-benchmark/cl2-pilot/summary_overall.csv`
- `ci-benchmark/cl2-pilot/pilot_report.md`
- `ci-benchmark/cl2-pilot/nongaussian/`: analogous outputs for Poisson/Binomial pilot

## Caveat

The pilot `pffr_cl2` implementation is still experimental. It supports
standard single-predictor GLM families (e.g., gaussian, poisson, binomial,
Gamma) and `gaulss` via a CL2-style leverage adjustment, but does not
currently handle other families with custom `family$sandwich` logic.
