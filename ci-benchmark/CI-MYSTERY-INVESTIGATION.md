# pffr CI Mystery Investigation

## Scope
This note investigates why `ci-benchmark/analysis.html` suggests very poor CI coverage for `pffr`/`mgcv` intervals, including whether the issue is:

1. benchmark contamination / analysis artifacts,
2. data regime difficulty (small n, strong dependence, high wiggliness),
3. fundamental CI construction problems.

The goal was to separate true statistical limits from pipeline issues.

## Executive Summary
Main conclusion: there is **no single fundamental bug in mgcv/pffr pointwise CIs** that explains all bad coverage.

What is happening is a combination of:

1. **Result-set contamination and analysis selection issues**.
2. **Expected undercoverage under strong within-curve correlation when using non-cluster variance estimators**.
3. **Finite-sample + smoothing-bias limitations in hard regimes even for robust methods**.

In easy IID Gaussian regimes, default `pffr` coverage is close to nominal.

## What Was Found

### 1) `analysis.html` is currently analyzing pilot files by default
`ci-benchmark/analysis.R` selected `ci-benchmark/pilot-results` first if present.

Current `analysis.html` confirms this:
- loaded `64` files
- only `7` unique DGPs
- replications range `4-10`

Those pilot sets are stress-oriented and not representative of broad benign behavior.

### 2) `ci-benchmark/results` contains mixed experiments under reused `dgp_id`
`ci-benchmark/results` has stale files from multiple designs mixed together.

Detected:
- `10` mixed `dgp_id`s with multiple parameter definitions:
  - `1, 12, 13, 24, 25, 36, 37, 48, 49, 61`
- legacy method present in same pool: `pffr_gls`

This makes pooled summaries unreliable if untreated.

### 3) Clean easy-regime results are near nominal for default `pffr`
From `ci-benchmark/easy-results/results_combined.rds` (IID, homoskedastic, Gaussian):

- `pffr` smooth coverage: `0.900`, `z_sd=1.011`
- `pffr` E(Y) coverage: `0.897`, `z_sd=1.005`
- `pffr` ff coverage: `0.878`, `z_sd=1.054`

This is not consistent with a universally broken CI implementation.

### 4) In clean pilot/expanded/targeted sets, failures track correlation structure
By `corr_type` (means over terms):

- `pilot-results`:
  - `pffr` AR1: `cov=0.771`, `z_sd=1.382`
  - `pffr` IID: `cov=0.899`, `z_sd=0.932`
- `pilot-results/expanded`:
  - `pffr` AR1: `0.716`, `1.479`
  - `pffr` IID: `0.919`, `0.911`
- `pilot-results/targeted`:
  - `pffr` AR1: `0.700`, `1.584`
  - `pffr` IID: `0.897`, `0.982`

Pattern: when dependence is strong and positive, naive/HC SEs understate uncertainty.

### 5) Under strong AR1 (`rho=0.9`), default and HC collapse as expected
From `ci-benchmark/results` subset (`ar1`, `rho=0.9`, homoskedastic, Gaussian):

- `pffr` ff: `cov=0.498`, `z_sd=2.469`
- `pffr_hc` ff: `cov=0.487`, `z_sd=2.569`
- `pffr_sandwich` ff: `cov=0.811`, `z_sd=1.244`
- `pffr_ar` ff: `cov=0.838`, `z_sd=1.126`

So robust/AR methods improve substantially, but still miss nominal in difficult finite-sample settings.

### 6) `seWithMean` is not the driver here
Targeted checks showed `seWithMean=TRUE` and `FALSE` are identical for key non-intercept terms in this benchmark setup (smooth/ff/linear/concurrent), with only slight intercept differences.

### 7) In IID Gaussian settings, model complexity is not the main failure mode
From exact-parameter IID subsets (`corr_type=iid`, `hetero=none`, `error_dist=gaussian`) in `ci-benchmark/results`:

- `pffr` mean coverage over terms stays close to nominal across `n`, `snr`, and `wiggliness`:
  - `n=50, snr=25, wig=0.01`: `cov=0.909`, `z_sd=0.974`
  - `n=50, snr=25, wig=5`: `cov=0.905`, `z_sd=0.995`
  - `n=200, snr=25, wig=0.01`: `cov=0.911`, `z_sd=0.939`
  - `n=200, snr=25, wig=5`: `cov=0.901`, `z_sd=0.979`

Only the smooth term at high wiggliness shows mild undercoverage (`~0.89`) with `z_sd > 1`, consistent with finite-sample smoothing bias, not a global CI breakdown.

### 8) Sample size matters under strong correlation for robust methods
For the exact strong-correlation subset (`ar1`, `rho=0.9`, `snr=25`, `wiggliness=5`, homoskedastic Gaussian):

- `n=50`:
  - `pffr`: `cov=0.567`, `z_sd=2.02`
  - `pffr_hc`: `cov=0.552`, `z_sd=2.16`
  - `pffr_sandwich`: `cov=0.739`, `z_sd=1.45`
  - `pffr_ar`: `cov=0.822`, `z_sd=1.13`
- `n=200`:
  - `pffr`: `cov=0.659`, `z_sd=1.65`
  - `pffr_hc`: `cov=0.651`, `z_sd=1.72`
  - `pffr_sandwich`: `cov=0.864`, `z_sd=1.06`
  - `pffr_ar`: `cov=0.872`, `z_sd=1.02`

So under strong dependence, low `n` plus complex smooth structure still hurts coverage even after robust correction, but increasing `n` clearly improves calibration for `pffr_sandwich` / `pffr_ar`.

## Theoretical Interpretation

### Why default `pffr` / mgcv pointwise CIs can fail in these simulations
- `pffr` is fit in long format (`n*T` rows).
- Without cluster-aware covariance, uncertainty calculations effectively treat rows as independent.
- Under strong within-curve positive dependence, effective sample size is much closer to `n` than `n*T`.
- Result: underestimated SEs (`z_sd >> 1`) and undercoverage.

This is expected behavior under misspecified independence assumptions, not a coding bug in CI extraction.

### Why sandwich is better but still not perfect
- Cluster-robust sandwich fixes the main dependence misspecification.
- But with modest cluster counts and complex smooths, CR1-style finite-sample behavior can still be anti-conservative.
- Remaining smoothing bias is not corrected by covariance-only adjustments.

## What Was Actually Broken

### Broken / misleading benchmark mechanics
1. Automatic use of pilot data in `analysis.R`.
2. Reuse of existing result files by filename only (no DGP-signature validation), allowing cross-design contamination.
3. Mixed legacy/new methods in pooled `results` (`pffr_gls` included).

### Not fundamentally broken
- mgcv/pffr pointwise CI machinery under IID benign settings.

## Code Changes Made

### `ci-benchmark/confint-benchmark.R`
`run_benchmark()` now validates existing result files against current DGP signatures before skipping.

- Mismatched stale files are **recomputed** (not skipped).
- Prevents silent contamination when DGP definitions change but filenames are reused.

### `ci-benchmark/analysis.R`
- Results directory resolution now prefers production results by default.
- Added env override: `REFUND_BENCHMARK_RESULTS_DIR`.
- Added fail-fast guard for mixed `dgp_id` definitions (override with `REFUND_ALLOW_MIXED_DGP=1`).

## Recommended Next Steps

1. Re-run benchmark in a fresh output directory with one design only.
2. Rebuild `analysis.html` from that directory explicitly (`REFUND_BENCHMARK_RESULTS_DIR=...`).
3. Report IID-benign and AR1-high-correlation regimes separately (do not pool into one headline number).
4. If nominality under dependence is required at small `n`, consider stronger small-sample cluster corrections (CR2/CR3) and/or bootstrap calibration.
