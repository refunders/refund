# refund 0.1-38

## Breaking changes

* **The default sandwich option for `pffr()` is now `sandwich = "cluster"`
  (previously `"none"`).** Existing code calling `pffr()` without explicit
  `sandwich=` will now get cluster-robust standard errors and confidence
  intervals. Set `sandwich = "none"` to restore the previous behavior.
  A one-time informational `message()` is printed when `sandwich` is not
  explicitly supplied.
* **`pffrGLS()` / `pffr_gls()` are deprecated and now error.** Their
  GLS-based covariance correction produced poorly calibrated inference.
  Use `pffr()` with `sandwich = "cluster"` (default) or `sandwich = "cl2"`
  instead.

## Function renames (old names deprecated)

* `pffrSim()` → `pffr_simulate()` (old name warns via `.Deprecated()`)
* `coefboot.pffr()` → `pffr_coefboot()`
* `qq.pffr()` → `pffr_qq()`
* `pffr.check()` → `pffr_check()`

## New features

* `pffr()` modularization: internal refactor into prepare → fit → postprocess
  pipeline (`pffr_prepare()`, `pffr_build_label_map()`, etc.) for better
  maintainability.
* `pffr_simulate()` (formerly `pffrSim()`) now supports a formula-based
  interface for specifying simulation models (e.g.,
  `pffr_simulate(Y ~ ff(X1) + xlin, effects = list(X1 = "cosine"))`).
  The new interface provides:
  - Customizable effect functions via preset libraries or user-defined functions
  - Access to true coefficient functions via the `truth` attribute
  - Support for non-Gaussian responses via the `family` argument
* The `scenario` argument in `pffr_simulate()` is deprecated. Use the formula
  interface instead. Legacy code using `scenario` will continue to work but
  will emit a deprecation warning.
* `pffr(..., sandwich = "cl2")` and `coef.pffr(..., sandwich = "cl2")` now
  support leverage-adjusted cluster-robust covariance (Bell-McCaffrey style
  CL2), including for `family = mgcv::gaulss()`.
* `coef.pffr()` now supports confidence intervals via `ci = "pointwise"` or
  `ci = "simultaneous"` in addition to standard errors. Simultaneous intervals
  are computed with a coefficient-level Gaussian simulation and max-|t|
  calibration over each smooth term's evaluation grid.
* `coef.pffr(ci = "simultaneous")` now uses a finite-sample
  `t_{G-1}` multiplier reference by default via `ci_ref = "t"`, which widens
  simultaneous bands at small numbers of independent curves or clusters. Set
  `ci_ref = "normal"` for the previous Gaussian-multiplier behavior.
* `pffr(..., sandwich = "cl2")` and `coef.pffr(..., sandwich = "cl2")` now
  report CL2 leverage diagnostics: the returned covariance carries
  `max_leverage` and `n_capped_clusters` attributes, and a warning is emitted
  when one or more clusters hit the leverage cap.
* `ff(..., check.ident = TRUE)` (the default) now also warns when the
  effective rank of the functional covariate's covariance is below
  `1.5 * k_s`, where `k_s` is the marginal basis dimension along `s`. The
  previous check only fired at the much weaker `rank < k_s`, so designs in
  which a sizeable part of the coefficient surface lies outside the span of
  the observed curves -- a bias floor that affects point estimates and
  interval coverage alike, without any fitting failure -- passed silently.
  The hard `rank < k_s` case is now reported as part of the same warning.
  See the "Effective rank and weak identifiability" section of `?ff` and the
  `ff-identifiability` vignette.
* The functional covariate simulated by `pffr_simulate()`'s legacy
  `scenario =` path is now drawn from a 12- rather than 7-dimensional spline
  basis, **doubling its effective rank from 5 to 10**. At the old rank the
  package's own simulated covariate was weakly identified against `ff()`'s
  default `k = 5` (which wants at least 7.5) and tripped the new warning
  above. The formula interface was already well clear of the threshold
  (effective rank 13-22, depending on `nxgrid`) and is unchanged. Simulated
  data from the `scenario =` path therefore differ from earlier versions.
* AR(1) support improvements: `pffr()` now automatically switches to
  `algorithm = "bam"` and `method = "fREML"` when `rho` is supplied, and
  sets `discrete = TRUE` for non-Gaussian families.
* New `vcov.pffr()` method. Plain `vcov()` returns the fit's stored
  covariance (robust, for sandwich-corrected fits); `vcov(..., sandwich =
  TRUE)` computes mgcv's HC sandwich from the model-based matrices.

## Bug fixes

* **Fixed double application of sandwich corrections on recomputation.**
  Since `pffr()` defaults to `sandwich = "cluster"`, fitted objects carry
  robust matrices in `Vp`/`Vc`/`Ve` — which the sandwich estimators also use
  as their (penalized) bread. Recomputing via `coef(fit, sandwich = ...)` or
  `coef(fit, cluster = ...)` therefore applied the correction on top of
  itself (SEs inflated ~1.5-2x on a small test example), and
  `coef(fit, sandwich = "none")` silently returned robust instead of
  model-based SEs. `pffr()` now stashes the model-based matrices in
  `$pffr$model_cov` and all recomputation paths restore them first.
  Objects fitted with earlier development versions lack the stash and now
  warn; refit to get correct recomputed results.
* The CR1 cluster sandwich (`sandwich = "cluster"`) now includes prior
  weights in the standard-GLM score, matching the CL2 and gaulss paths;
  CR1 SEs for weighted fits were wrong before.
* `pffr_coefboot(type = "norm")` returned all-NA intervals because
  `boot::boot.ci()` names its result element `$normal`, not `$norm`; the
  element names are now mapped correctly. `type = "stud"` errors up front
  (the bootstrap statistic provides no replicate variances) instead of
  silently yielding all-NA intervals.
* `coef.pffr(cluster = ...)` with a single cluster now errors cleanly
  instead of dividing by zero in the CR1 small-sample factor.

# refund 0.1-37

* Fixed invalid email address for Yakuan Chen

# refund 0.1-36

* Fix threshold for small sigma2 in mfpca.face function to ensure accurate score estimation for level1 scores
* Added parameters npc2 and pve2 to mfpca.face to allow for user to specify separate npc or pve for level 2 decomposition


# refund 0.1-35

* One line fix to pfr that allows pfr to be called from within another function
* Pull request on URL updates in documentation

# refund 0.1-34

* Added pve to what is returned by mfpca.face
* Change names in list of what is returned by mfpca.face to be the same as mfpca.sc
* Updated email address for Erjia Cui
* Changed Erjia from contributor to author


# refund 0.1-33

* Added pve to what is returned by fpca.sc and fpca.face

# refund 0.1-32

* Simon Wood removed "pers" argument from mgcv::plot.random.effect. Removed reference to pers from plot_pfr.gam to avoid documentation and code being inconsistent.

# refund 0.1-31

* added COVID19 datasets
* minor bug fix for mfpca.face
* updated email address for package maintainer


# refund 0.1-30

* removed import lattice::qq from pffr-methods.R

# refund 0.1-28	

* added periodic spline option for `fpca.face`
* changed if(class(object) != "string") to if(!inherits(object, "string")) in `ccb.fpc.R` and `fosr.perm.test.R` files to fix Note.


# refund 0.1-27	

* bug fix for `pfr` models without intercept (thx, @ZheyuanLi)

# refund 0.1-26	

* New function, `mfpca.face`, which is a faster version of `mfpca.sc`


# refund 0.1-25	

* Minor bug fixes in `fpca.face` to patch error when npc = 1


# refund 0.1-24	

* Minor bug fixes in `pfr` to patch error in R version 4.1
* Updated documentation URLs

# refund 0.1-23

* Minor bug fixes in `fgam` examples for upcoming R release

# refund 0.1-22

* Fixes bugs 
  * Commented out option `useSymm = TRUE` in tests for `fpca.sc`

# refund 0.1-21

* Fixes a bug in `fpc()` due to release of R 4.0.0 that changes the following:

```
 R> class(matrix(1 : 4, 2, 2))
 [1] "matrix" "array" 

(and no longer just "matrix" as before), and that conditions of length
greater than one in 'if' and 'while' statements executing in the package
being checked give an error.
```

* change email address for maintainer

# refund 0.1-20

* fix minor bug in rlrt.pfr.R

* change maintainer from Rayman Huang to Julia Wrobel 

# refund 0.1-19


* updates for compatibility with mgcv 1.8-23  (#69 etc.)

* fixed fpcr for scalar covariates (#76)

* now re-exports cmdscale_lanczos

# refund 0.1-15

* homogenized inputs/outputs to most fpca.XXX functions

* fix documentation for fpca.face

* added pco ridge regression see ?poridge

# refund 0.1-14

* add `fpca.lfda()` function

* add functions and data set from refund.shiny package

# refund 0.1-13

* add `mfpca.sc()` function.

* add example for DTI data.

* export `pfr_old()`.

* fix documentation and add warning message for `rlrt.pfr()`.

# refund 0.1-12

* new `pfr()` function. The new `pfr()` function merged the old pfr and fgam functionality.

* Switch to roxygen2 documentation

* allow an "argvals" argument to all relevant functions for consistent strucuture throughout the package

* refactor input and output argument for `fosr()`.

* add vignettes.

* bug fix in `fpca.face()`.

* export `expand.call()`.

* add cross-sectional FoSR using GLS, variational Bayes and Gibbs sampler
