# refund (development version)

* **`pffr()` inference now implements one procedure.** By default
  (`sandwich = TRUE`) standard errors and intervals use the curve-clustered
  CL2 sandwich with the full-block leverage adjustment in its Bayesian form,
  and pointwise intervals from `coef()`, `predict(se.fit = TRUE)` (which now
  returns `crit` and `df`) and the new `pffr_predict_ci()` use Satterthwaite
  critical values with a separate df for every point. `sandwich = FALSE`
  gives the model-based intervals, which are invalid under within-curve
  dependence. `cluster =` clusters several curves per subject; `vcov()`
  returns the fit's covariance. The previous `G <= 150` / `max D_g <= 500`
  limits are gone. pffr warns below 40 curves (the procedure was evaluated
  down to 20), messages once per session that intervals for intercepts and
  coefficient functions of scalar covariates of binary responses can
  undercover, and messages once that intervals around NCV estimates are not
  for inference (use NCV for the shape of coefficient surfaces, REML for
  intervals). Families without a cluster score path, `gamm`/`gamm4` and AR(1)
  fits fall back to model-based intervals with a message.
  Removed (never released): `sandwich = "auto"`, `cl2_adjustment`,
  `dof_correction`, `edf_type`, `crit`, `ci_ref`, `df_gram`, `bias_ref`,
  `se_method`, `pffr_jackknife_se()`, `pffr_dependence_check()`,
  `pffr_upgrade_fit()` and the `refund.pffr.autopolicy` option. Deprecated:
  the character values `"cluster"`, `"hc"`, `"cl2"`, `"none"` of `sandwich`
  (the first two now give CL2), `coef(freq = )`, and `pffr_coefboot()` /
  `coefboot.pffr()` (valid but 1.4-1.6 times wider intervals at a much higher
  cost). `ffpc()` now keeps the FPCs explaining 95% of the covariate's
  variance (`pve = 0.95`, was 0.99, which gave very noisy estimates).
* `pffr(method = "NCV")` now leaves out whole curves (or `cluster` groups)
  by default. `ncv_blocks = "point"` provides pointwise comparisons; explicit
  `nei` takes precedence. Dual mgcv neighbourhood names, a cached behavioural
  check, and model-frame alignment checks prevent silent pointwise fallbacks.
  NCV supports `algorithm = "gam"` (the default even for large data) and
  explicit `algorithm = "bam"`, which sets `discrete = TRUE` and rejects
  `discrete = FALSE`. Both use retained model-frame row indices for `nei`;
  bam requires `a/ma/d/md` names. Discrete NCV inverts a neighbourhood-sized
  matrix per neighbourhood and is usually slower than gam for whole-curve
  blocks; `nei$sample` allows sub-sampling. Nonzero AR1 `rho` is unavailable
  with NCV and errors. With bam, neighbourhoods must have equal size (mgcv's
  discrete NCV corrupts memory otherwise), so missing responses and irregular
  grids need `algorithm = "gam"`.

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
* **The Simpson integration weights used by `ff()` and `sff()` are fixed.**
  The `integration = "simpson"` weights were scaled by
  `(b - a) / (3 * nxgrid)` instead of `(b - a) / (3 * (nxgrid - 1))`, and for
  even `nxgrid` the `[1, 4, 2, ..., 4, 1]` alternation ended in `2` before the
  closing `1`, which composite Simpson does not allow. The weights therefore
  summed to less than the length of the integration domain: a constant on
  `[0, 1]` integrated to 0.956 at `nxgrid = 30`, 0.968 at `nxgrid = 31`, 0.978
  at `nxgrid = 60`, 0.984 at `nxgrid = 61` and 0.989 at `nxgrid = 93` (the DTI
  CCA grid) instead of 1. Estimated `ff()` coefficient surfaces were rescaled
  by the reciprocal of that factor, i.e. inflated by up to ~4.5% on typical
  grids. Simulation studies in which the same weights generated *and* fitted
  the data are unaffected; real-data fits are. The weights now implement
  composite Simpson's rule with `h = (b - a) / (nxgrid - 1)`, using Simpson's
  3/8 rule on the last three intervals when `nxgrid` is even, so that a
  constant integrates to exactly `b - a` and cubics are integrated exactly for
  every `nxgrid >= 3`. The old behaviour is still reachable as
  `integration = "simpson_legacy"` for reproducing results from earlier
  versions; it is deprecated and should not be used for new analyses.

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
* The CL2 sandwich now checks the penalized hat matrix against its own bounds
  (`0 <= h_ii <= 1`, `0 <= eigen(H_gg) <= 1`, and positive semi-definiteness of
  the exact Bell--McCaffrey block) and **warns** when they are violated by more
  than round-off. Such a violation means the model-based bread and the weighted
  design have become numerically inconsistent -- a barely converged or
  extremely ill-conditioned fit -- so the cluster-robust covariance is
  meaningless however plausible it looks. Previously this produced silently
  exploded standard errors (interval widths up to 1e133 were observed on
  degenerate Poisson fits, with nothing to distinguish them from a legitimately
  wide interval). The covariance is still returned, now carrying
  `max_obs_leverage` and `hat_invariant_violation` attributes. Inspect and
  refit the model.
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
## Bug fixes

* **Fixed double application of sandwich corrections on recomputation.**
  Since `pffr()` defaults to `sandwich = "cluster"`, fitted objects carry
  robust matrices in `Vp`/`Vc`/`Ve` — which the sandwich estimators also use
  as their (penalized) bread. Recomputing via `coef(fit, sandwich = ...)` or
  `coef(fit, cluster = ...)` therefore applied the correction on top of
  itself (SEs inflated ~1.5-2x on a small test example), and
  `coef(fit, sandwich = "none")` silently returned robust instead of
  model-based SEs. `pffr()` now keeps `$Vp`/`$Vc`/`$Ve` model-based and
  stores the robust covariance in `$pffr$Vsandwich`. Objects fitted with
  refund 0.1-40 lack the model-based matrices; their stored covariance is used
  with Gaussian critical values and a warning; refit for the current
  intervals.
* `pffr_coefboot(type = "norm")` returned all-NA intervals because
  `boot::boot.ci()` names its result element `$normal`, not `$norm`; the
  element names are now mapped correctly. `type = "stud"` errors up front
  (the bootstrap statistic provides no replicate variances) instead of
  silently yielding all-NA intervals.
* `coef.pffr(cluster = ...)` with a single cluster now errors cleanly
  instead of dividing by zero in the small-sample factor.

# refund 0.1-40

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
* AR(1) support improvements: `pffr()` now automatically switches to
  `algorithm = "bam"` and `method = "fREML"` when `rho` is supplied, and
  sets `discrete = TRUE` for non-Gaussian families.
* Removed dependency on `mgcv::plot.random.effect`

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
