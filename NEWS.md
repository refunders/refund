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
