# What seems to be going wrong (analysis notes)

This note summarizes why the current CI benchmark can behave very differently
from intuitive expectations (e.g., sandwich not helping; GLS failing), even
after recent DGP changes (balanced term magnitudes, correlation recalibration).

## 1) Sandwich is targeting the wrong failure mode

`pffr_sandwich` uses `mgcv::vcov.gam(..., sandwich = TRUE)` (via `coef.pffr`),
which is essentially a heteroskedasticity-robust covariance for *independent*
observations.

In this benchmark, the hard cases are largely *within-curve dependence across t*
(AR(1), Gaussian, Fourier). A plain Huber–White style sandwich does not make
inference valid under strong serial correlation unless it is explicitly
clustered / dependence-aware. So it is not surprising when sandwich correction
often provides little or no improvement under correlated errors.

## 2) “Half-grid correlation ≈ 0.4” implies *extremely* high adjacent correlation

Calibrating correlation so that `Corr(t, t + 0.5 * range(t_grid)) ≈ 0.4` makes
adjacent-time correlation very close to 1 on typical grids:

- AR(1): `rho_adj` can be ~0.96+ (e.g., `nygrid ≈ 45`)
- Gaussian: with `phi` chosen to hit 0.4 at half-range, `Corr(dt)` is often
  ~0.99+

This creates covariance matrices that are very ill-conditioned / near-singular.
Downstream inference methods that require inversion or square roots of Σ become
numerically fragile.

## 3) `pffr_gls` is feasible GLS with a noisy, high-dimensional Σ estimate

`pffr_gls` relies on `hatSigma`, which is currently estimated from pilot-fit
residuals as an unstructured `nygrid × nygrid` covariance matrix.

This is difficult even when the mean is correctly specified. Here it is harder
because:

- the pilot mean fit is penalized/smoothed and can be biased;
- bias leaks into residuals and contaminates Σ estimation;
- estimating an unstructured covariance in moderate sample sizes is unstable;
- strong correlation makes Σ close to singular.

Once Σ is unstable, GLS can easily degrade performance relative to plain `pffr`
or appear to “fail” (extreme standard errors, erratic coverage, warnings).

## 4) The GLS implementation explicitly *modifies* Σ when it is ill-conditioned

`compute_sqrt_sigma_inv()` checks Σ’s condition number and, if it exceeds a
cutoff, projects Σ into the positive-definite cone via `Matrix::nearPD(...)`.

That changes the intended covariance model. When this happens frequently
(which is likely under strong correlation), the method being evaluated is no
longer “GLS under the DGP covariance”; it is “GLS under a regularized surrogate
covariance”. This can strongly affect both point estimates and standard errors.

## 5) Some undercoverage is bias/identifiability, not SE choice

If coverage problems are driven by smoothing bias, term confounding, or weak
identifiability (e.g., low effective rank for `X1` relative to `k`), then
changing the SE estimator (sandwich vs non-sandwich) will not fix the problem.

In that regime, the benchmark may mainly be measuring bias/regularization
effects, not the robustness properties of the covariance estimator.

## 6) The current coverage metric can mask localized improvements

The benchmark computes `coverage` per replicate as a *fraction of grid points*
covered, then averages those proportions.

If an estimator helps mostly in specific regions (e.g., high-variance portions
of a heteroskedastic pattern), that improvement can be diluted by many grid
points where coverage is already near nominal.

## Bottom line

The observed behavior is consistent with the experiment being dominated by:

1. strong within-curve correlation (by construction),
2. unstable feasible GLS (full Σ estimation + inversion/regularization),
3. bias/identifiability and penalty effects,

rather than the marginal impact of a heteroskedasticity-robust covariance
estimator like the default mgcv sandwich.

