# Exact-CL2 implementation notes

Date: 2026-07-21

## Change made

`gam_sandwich_cluster_cl2()` now implements the exact Bell--McCaffrey
leverage block

```
B_g = I - 2 H_gg + (H_t^2)_gg
```

where `(H_t^2)_gg` is evaluated as
`(Xw_g Vp) (Xw' Xw) (Xw_g Vp)'`. Its eigendecomposition floors eigenvalues at
`(1 - leverage_cap)^2`; the historical shortcut continues to cap the
eigenvalues of `H_gg` and uses `(I - H_gg)^(-1/2)`.

`pffr()` and `coef.pffr()` have a new `cl2_adjustment` argument:
`"auto"` (default), `"exact"`, or `"shortcut"`. Fit metadata records the
resolved adjustment plus `n_adjusted`, `min_block_eig`, and
`max_block_kappa` for exact CL2. Legacy shortcut diagnostics remain available.

## Default rule

Within the existing CL2 path only, `"auto"` selects exact CL2 when:

* `G <= 100`; and
* `G * max(D_g) * p^2 <= 5e9` (raised from 5e8, see Addendum below).

The first bound confines exact CL2 to the small-to-moderate-cluster regime
where its finite-sample correction is relevant. The second is a conservative
dense-matrix proxy for the per-block multiplication
`(Xw_g Vp) (Xw'Xw)`. Outside either bound the package visibly records and uses
the shortcut. The package-level `sandwich = "auto"` policy is unchanged.

## Tests run

All R invocations used `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1` and
`nice -n 15`.

* `Rscript hardening-exactcl2.R` — passed; wrote
  `hardening-exactcl2-results.csv`.
* Targeted test files run with `NOT_CRAN=true` — passed:
  `test-pffr-exactcl2.R`, `test-pffr-s2-autopolicy.R`,
  `test-pffr-sandwich-storage.R`, `test-pffr-sandwich-refit.R`,
  `test-pffr-satterthwaite.R`, and `test-pffr.R`.
* `devtools::load_all()` — passed.
* `devtools::document()` — completed and regenerated the affected Rd files.
  It printed pre-existing roxygen warnings about unrelated tags/method exports
  in `pffr-core.R`, `poridge.R`, `pffr-methods.R`, and plot methods.

## Hardening output summary

The 16-cell table covers families Gaussian, Poisson, binomial, and Gamma;
small/moderate response bases `k = 4, 7`; and unbalanced and
influential-cluster designs. Each has five unbalanced clusters with curve
counts `1, 1, 2, 3, 5`; the influential design gives one curve an extreme
covariate value while retaining an ordinary-range outcome.

Actual results:

* All 16 exact covariances and coefficient SE vectors were finite and strictly
  positive.
* The median exact/shortcut SE ratio ranged from 1.00 to 1.02.
* Three influential Poisson/Gamma blocks hit the shipped exact floor
  (`n_adjusted = 1`); all other cells reported zero adjusted blocks.
* The smallest observed pre-floor block eigenvalue was `3.20e-07`.
* Exact versus the `cl2_exact_nocap` analogue (floor `tol = 1e-8`) differed by
  at most `8.76e-06`; the nonzero differences occur in the floored cells.
* The largest reported block condition proxy was `2.98e+06`.

## Open questions

The pointwise Satterthwaite degrees-of-freedom routine remains the established
working-iid approximation based on the shortcut adjustment; this change
updates the CL2 covariance block only. Exact Bell--McCaffrey df remains a
separate methodological task, as already documented in the package.

The requested local commits could not be created in this execution environment:
Git returned `Unable to create .../.git/worktrees/refund-wt-s1/index.lock:
Read-only file system`. The branch `exact-cl2-default` exists and all changes
remain unstaged in this worktree; no unrelated `tests/testthat/Rplots.pdf`
change was touched.

## Addendum 2026-07-21 (PI-side review): cost cap raised 5e8 -> 5e9

Timing at the shipped cap showed the exact block's marginal cost over the
shortcut is ~0.03 s (both paths pay the same per-cluster eigendecomposition;
the marginal cost is only the extra multiplications), and ~0.7 s total at
10x the cap. 5e8 could therefore deny a sub-second correctness upgrade in
the saturated small-G regime where Study EX measured the largest gain
(+1.3pp at G=20). Cap raised to 5e9; G <= 100 unchanged (EX-validated
range; gain < 0.15pp beyond).
