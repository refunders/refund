# gaulss CL2 Method and Package Integration

## 1) Scope

This note documents the `sandwich = "cl2"` implementation integrated into
`pffr`, with emphasis on `family = mgcv::gaulss()` (Gaussian location-scale).
It also summarizes behavior for non-Gaussian single-linear-predictor families
(e.g., Poisson, Binomial).

## 2) Target Estimator

For cluster `g` (curve-level clustering), CL2 uses

\[
\widehat V_{CL2} = c_{HC1}\, V_p \left( \sum_{g=1}^G U_g U_g^\top \right) V_p + B_2,
\]

with

- \(c_{HC1} = G/(G-1)\),
- \(B_2 = 0\) for frequentist covariance (`freq = TRUE`),
- \(B_2 = V_p - V_e\) for Bayesian/posterior covariance (`freq = FALSE`),
- \(U_g = X_{w,g}^\top A_g z_g\),
- \(A_g = (I - H_{gg})^{-1/2}\),
- \(H_{gg} = X_{w,g} V_p X_{w,g}^\top\).

`H_{gg}` is symmetrized and its eigenvalues are capped at `0.999` for numeric
stability; inverse-square-root eigenvalues are floored at `1e-8`.

## 3) Working Representation (Standard Families)

For standard one-LP families, per-row score contribution is factorized as

\[
s_i = x_i \cdot \frac{\mu'(\eta_i)(y_i-\mu_i) w_i}{\phi\,V(\mu_i)}
    = (x_i\,\mu'(\eta_i)\sqrt{w_i/(\phi V(\mu_i))}) \cdot ((y_i-\mu_i)\sqrt{w_i/(\phi V(\mu_i))}).
\]

This gives:

- `Xw[i, ] = X[i, ] * mu_eta[i] * sqrt_common[i]`
- `z[i] = (y[i] - mu[i]) * sqrt_common[i]`
- `sqrt_common[i] = sqrt(prior_weights[i] / (sig2 * variance(mu[i])))`

Non-finite rows are zeroed.

## 4) gaulss Working Representation

`gaulss` has two LPs (location and scale), with fitted values arranged as:

- `mu` block (first `n_obs`),
- `tau = 1/sigma` block (second `n_obs`).

Using the same score components as existing `compute_gaulss_scores()`:

- \(w_1 = (\tau^2 r)\,\mu'_{1}(\eta_1)\), with \(r = y-\mu\)
- \(w_2 = (1/\tau - \tau r^2)\,\mu'_{2}(\eta_2)\)

(optionally multiplied by prior weights).

Because these weights can be negative, we use pseudo-observation expansion:

1. Build block 1 rows from LP1 coefficients only, with scaling `sqrt(abs(w1))`.
2. Build block 2 rows from LP2 coefficients only, with scaling `sqrt(abs(w2))`.
3. Set pseudo-response vector pieces to `sign(wk) * sqrt(abs(wk))`.
4. Stack both blocks:
   - `Xw = rbind(Xw_block1, Xw_block2)`
   - `z = c(z_block1, z_block2)`
   - `cluster_id = rep(cluster_id_original, times = 2)`.

Then `crossprod(Xw, z)` recovers the original summed score contribution.

## 5) Package Integration Points

Implemented in package source:

- `R/pffr.R`
  - `sandwich` options now include `"cl2"`.
- `R/pffr-core.R`
  - Added `sym_inv_sqrt()`,
    `build_cl2_working_standard()`,
    `build_cl2_working_gaulss()`,
    `gam_sandwich_cluster_cl2()`.
  - `apply_sandwich_correction()` now routes `"cl2"` to CL2 covariance.
- `R/pffr-methods.R`
  - `coef.pffr(..., sandwich = "cl2")` supported.
  - Recomputes CL2 on-the-fly if model was fit with a different sandwich type.

## 6) Fallback Policy (Custom `family$sandwich`)

For families with custom `family$sandwich` and not `gaulss`, CL2 is currently
not implemented in cluster form. The implementation warns and falls back to:

`mgcv::vcov.gam(..., sandwich = TRUE, freq = ...)`

This matches the agreed behavior for unsupported custom families.

## 7) Validation Added in `testthat`

New/extended tests in `tests/testthat/test-pffr.R` cover:

- Gaussian CL2 integration path (`pffr(..., sandwich="cl2")` + `coef(..., sandwich="cl2")`).
- Non-Gaussian one-LP families: Poisson and Binomial CL2 covariance.
- gaulss CL2 covariance and end-to-end fit path.
- Recompute logic in `coef.pffr` when requested sandwich type differs from fit.
- Fallback-to-HC behavior for unsupported custom `family$sandwich`.

Heavy tests are guarded with `skip_on_cran()`.

## 8) Known Limits

- CL2 is not yet implemented for arbitrary custom multi-LP families beyond
  `gaulss`.
- Runtime overhead is higher than CR1 due to per-cluster eigen decompositions
  of `H_{gg}`/`I-H_{gg}` blocks.
- Very high-leverage clusters require eigenvalue capping; this is expected in
  very small-sample/high-edf settings.

## 9) Practical Interpretation

- `cluster` (CR1) remains a robust default.
- `cl2` is intended for small-to-moderate sample sizes where leverage
  correction improves finite-sample behavior (at the cost of wider intervals).
- For gaulss, pseudo-observation expansion yields a practical extension that
  preserves score structure and integrates cleanly with existing `pffr` code.
