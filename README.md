# `refund`
[![](https://travis-ci.org/refunders/refund.svg?branch=master)](https://travis-ci.org/refunders/refund)

[![](http://cranlogs.r-pkg.org/badges/refund)](https://cran.rstudio.com/web/packages/refund/index.html)

## Methods for regression with functional data

These packages implement various approaches to functional data regression. 

Regression with scalar responses and functional predictors is implemented in functions `pfr`, `peer`, `lpeer`, `fpcr` and `fgam`. For regression with functional responses, see `pffr`, `fosr`, and `fosr2s`.

Regularized covariance and FPC estimation is implemented in functions `fpca.sc`,
`fpca.ssvd`, `fpca.face`, `fpca2s`.


Shiny-based interactive graphics for visualizing results from `fpca` and regression methods in `refund` can be generated using the `plot_shiny()` function in the `refund.shiny` package.


Wavelet-based functional regression methods with scalar responses and functional predictors can be found in the `wcr` and `wnet` functions in the `refund.wave` package.

---------------

### Installation

To install the latest patched version directly from Github, please use `devtools::install_github("refunders/refund")` for `refund` and `devtools::install_github("refunders/refund.shiny")` for `refund.shiny` and `devtools::install_github("refunders/refund.wave")` for `refund.wave`.

To install the developer version with experimental features directly from Github, please use `devtools::install_github("refunders/refund", ref="devel")`.

