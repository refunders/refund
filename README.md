# `refund`
## Methods for regression with functional data

These packages implement various approaches to functional data regression. 

Regression with scalar responses and functional predictors is implemented in functions `pfr`, `peer`, `lpeer`, `fpcr` and `fgam`. For regression with functional responses, see `pffr`, `fosr`, and `fosr2s`.

Regularized covariance and FPC estimation is implemented in functions `fpca.sc`,
`fpca.ssvd`, `fpca.face`, `fpca2s`.

Wavelet-based functional regression methods with scalar responses and functional predictors can be found in the `wcr` and `wnet` functions in the `refund.wave` package.

---------------

### Installation

To install the latest patched version directly from Github, please use `devtools::install_github("refunders/refund/refund")` for `refund` and `devtools::install_github("refunders/refund/refund.wave")` for `refund.wave`.

To install the developer version `refundDevel` with experimental features directly from Github, please use `devtools::install_github("refunders/refund/refund", ref="devel")`.

