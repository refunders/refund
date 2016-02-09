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