# refund 0.1-38

* `pffrSim()` now supports a formula-based interface for specifying simulation
  models (e.g., `pffrSim(Y ~ ff(X1) + xlin, effects = list(X1 = "cosine"))`).
  The new interface provides:
  - Customizable effect functions via preset libraries or user-defined functions
  - Access to true coefficient functions via the `truth` attribute
  - Support for non-Gaussian responses via the `family` argument
* The `scenario` argument in `pffrSim()` is deprecated. Use the formula
  interface instead. Legacy code using `scenario` will continue to work but
  will emit a deprecation warning.
* Added internal helper functions for formula-driven data simulation:
  `parse_pffr_formula()`, `compute_ff_effect()`, `compute_linear_effect()`,
  `compute_smooth_effect()`, `compute_const_effect()`, `compute_intercept()`,
  `generate_functional_covariate()`, `generate_scalar_covariate()`.
* `pffr(..., sandwich = "cl2")` and `coef.pffr(..., sandwich = "cl2")` now
  support leverage-adjusted cluster-robust covariance (Bell-McCaffrey style
  CL2), including for `family = mgcv::gaulss()`.
* The default sandwich option for `pffr()` and `coef.pffr()` is now
  `sandwich = "cluster"` (previously `"none"`).
* `coef.pffr()` now supports confidence intervals via `ci = "pointwise"` or
  `ci = "simultaneous"` in addition to standard errors. Simultaneous intervals
  are computed with a coefficient-level Gaussian simulation and max-|t|
  calibration over each smooth term's evaluation grid.

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
