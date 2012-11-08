\name{pffrGLS}
\alias{pffrGLS}
\title{Penalized function-on-function regression with non-i.i.d. residuals}
\usage{
  pffrGLS(formula, yind, hatSigma, algorithm = NA,
    method = "REML", tensortype = c("te", "t2"),
    bs.yindex = list(bs = "ps", k = 5, m = c(2, 1)),
    bs.int = list(bs = "ps", k = 20, m = c(2, 1)), ...)
}
\arguments{
  \item{formula}{a formula with special terms as for
  \code{\link[mgcv]{gam}}, with additional special terms
  \code{\link{ff}()} and \code{c()}. See
  \code{\link[refund]{pffr}}.}

  \item{yind}{a vector with length equal to the number of
  columns of the matrix of functional responses giving the
  vector of evaluation points \eqn{(t_1, \dots ,t_{G})}.
  see \code{\link[refund]{pffr}}}

  \item{algorithm}{the name of the function used to
  estimate the model. Defaults to \code{\link[mgcv]{gam}}
  if the matrix of functional responses has less than
  \code{2e5} data points and to \code{\link[mgcv]{bam}} if
  not. "gamm" (see \code{\link[mgcv]{gamm}}) and "gamm4"
  (see \code{\link[gamm4]{gamm4}}) are valid options as
  well.}

  \item{hatSigma}{(an estimate of) the within-observation
  covariance (along the responses' index), evaluated at
  \code{yind}. See Details.}

  \item{method}{See \code{\link[refund]{pffr}}}

  \item{bs.yindex}{See \code{\link[refund]{pffr}}}

  \item{bs.int}{See \code{\link[refund]{pffr}}}

  \item{tensortype}{See \code{\link[refund]{pffr}}}

  \item{...}{additional arguments that are valid for
  \code{\link[mgcv]{gam}} or \code{\link[mgcv]{bam}}. See
  \code{\link[refund]{pffr}}.}
}
\value{
  a fitted \code{pffr}-object, see
  \code{\link[refund]{pffr}}.
}
\description{
  Implements additive regression for functional and scalar
  covariates and functional responses. This function is a
  wrapper for \code{mgcv}'s \code{\link[mgcv]{gam}} and its
  siblings to fit models of the general form \cr
  \eqn{Y_i(t) = \mu(t) + \int X_i(s)\beta(s,t)ds +
  f(z_{1i}, t) + f(z_{2i}) + z_{3i} \beta_3(t) + \dots +
  E_i(t))}\cr with a functional (but not necessarily
  continuous) response \eqn{Y(t)}, (optional) smooth
  intercept \eqn{\mu(t)}, (multiple) functional covariates
  \eqn{X(t)} and scalar covariates \eqn{z_1}, \eqn{z_2},
  etc. The residual functions \eqn{E_i(t) \sim GP(0,
  K(t,t'))} are assumed to be i.i.d. realizations of a
  Gaussian process. An estimate of the covariance operator
  \eqn{K(t,t')} evaluated on \code{yind} has to be supplied
  in the \code{hatSigma}-argument. Note that this has to be
  positive definite. If \code{hatSigma} is close to
  positive \emph{semi-}definite or badly conditioned,
  estimated standard errors become unstable (typically much
  too small). \code{pffrGLS} will try to diagnose this and
  issue a warning. The danger is especially big if the
  number of functional observations is smaller than the
  number of gridpoints (i.e, \code{length(yind)}), since
  the raw covariance estimate will not have full rank.
}
\section{Details}{
  Please see \code{\link[refund]{pffr}} for details on
  model specification and implementation. THIS IS AN
  EXPERIMENTAL VERSION AND NOT WELL TESTED YET -- USE AT
  YOUR OWN RISK.
}
\author{
  Fabian Scheipl
}
\seealso{
  pffr, fpca.sc
}

