# CL2 Pilot Handover Notes

Date: 2026-02-17

## Scope completed

1. Added a shared experimental CL2 implementation:
   - `ci-benchmark/cl2-pilot/cl2-utils.R`
   - Supports:
     - standard GLM families with one linear predictor (`gaussian`, `poisson`, `binomial`, `Gamma`, etc.)
     - `gaulss` via a pseudo-observation expansion of the two LP score components
   - Still not implemented for other custom `family$sandwich` families.

2. Refactored Gaussian pilot to use shared CL2 code:
   - `ci-benchmark/cl2-pilot/sim-study-cl2-pilot.R`

3. Added non-Gaussian pilot (Poisson/Binomial, low n):
   - `ci-benchmark/cl2-pilot/sim-study-cl2-pilot-nongaussian.R`

4. Added gaulss quick-test runner:
   - `ci-benchmark/cl2-pilot/quick-test-gaulss.R`

5. Updated documentation:
   - `ci-benchmark/cl2-pilot/README.md`

## Commands run

1. Gaussian pilot (already existing results reused):
   - `Rscript ci-benchmark/cl2-pilot/sim-study-cl2-pilot.R pilot`

2. Non-Gaussian pilot:
   - `Rscript ci-benchmark/cl2-pilot/sim-study-cl2-pilot-nongaussian.R pilot`

3. gaulss quick test:
   - `Rscript ci-benchmark/cl2-pilot/quick-test-gaulss.R pilot`

## Result summary

### Gaussian pilot (`n={20,40}`, `corr={iid,ar1}`)

From `ci-benchmark/cl2-pilot/summary_cluster_vs_cl2.csv`:

- `n=20, ar1`: cluster `0.665` -> cl2 `0.858` (`+0.193`), width ratio `1.61`
- `n=40, ar1`: cluster `0.797` -> cl2 `0.865` (`+0.069`), width ratio `1.19`
- IID gains are smaller (`+0.042` to `+0.056`), width ratio `1.11` to `1.19`

### Non-Gaussian pilot (Poisson/Binomial, non-identity links)

From `ci-benchmark/cl2-pilot/nongaussian/summary_cluster_vs_cl2.csv`:

- Binomial, AR1: cluster `0.807` -> cl2 `0.890` (`+0.083`), width ratio `1.26`
- Binomial, IID: cluster `0.843` -> cl2 `0.870` (`+0.026`), width ratio `1.06`
- Poisson, AR1: cluster `0.706` -> cl2 `0.810` (`+0.104`), width ratio `1.31`
- Poisson, IID: cluster `0.817` -> cl2 `0.857` (`+0.040`), width ratio `1.11`

### gaulss quick test

From `ci-benchmark/cl2-pilot/gaulss-quick/summary_compare.csv`:

- AR1 aggregate: cluster `0.785` -> cl2 `0.811` (`+0.026`), width ratio `1.06`
- IID aggregate: cluster `0.833` -> cl2 `0.843` (`+0.010`), width ratio `1.02`

Status:

- gaulss path runs end-to-end and returns valid metrics (`72/72` rows converged).
- CL2 effect for gaulss is directionally similar but smaller in this tiny quick test.

## Output locations

1. Gaussian:
   - `ci-benchmark/cl2-pilot/pilot_report.md`
   - `ci-benchmark/cl2-pilot/summary_cluster_vs_cl2.csv`

2. Non-Gaussian:
   - `ci-benchmark/cl2-pilot/nongaussian/pilot_report.md`
   - `ci-benchmark/cl2-pilot/nongaussian/summary_cluster_vs_cl2.csv`

3. gaulss quick:
   - `ci-benchmark/cl2-pilot/gaulss-quick/results.rds`
   - `ci-benchmark/cl2-pilot/gaulss-quick/summary_compare.csv`

## Important implementation notes

1. CL2 numerical safeguards:
   - leverage eigenvalue cap: `0.999`
   - inverse-square-root eigen floor: `1e-8`

2. gaulss handling:
   - uses mu and tau score components (same formulas as existing cluster gaulss score code in core package),
   - converts each observation into 2 pseudo-rows (one per LP component),
   - cluster IDs are duplicated accordingly.

3. Current scope is benchmark-only code under `ci-benchmark/cl2-pilot/`; core package API is unchanged.

## Recommended next steps

1. If promoting to package code, add `sandwich = "cl2"` to core pffr API and route through `R/pffr-core.R`.
2. Add unit tests in `tests/testthat/test-pffr.R`:
   - CL2 runs for gaussian/poisson/binomial/gaulss,
   - CL2 covariance is PSD/symmetric and differs from CR1 in low-n cases.
3. Run a larger non-Gaussian + gaulss study (`>=20` reps/cell) before defaulting behavior.
4. Benchmark runtime overhead vs cluster CR1 (CL2 is currently slower due per-cluster eigendecompositions).
