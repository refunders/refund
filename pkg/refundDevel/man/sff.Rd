\name{sff}
\alias{sff}
\title{Construct a smooth function-on-function regression term}
\usage{
  sff(X, yind, xind = seq(0, 1, l = ncol(X)),
    basistype = c("te", "t2", "s"),
    integration = c("simpson", "trapezoidal"), L = NULL,
    limits = NULL,
    splinepars = list(bs = "ps", m = c(2, 2, 2)))
}
\arguments{
  \item{X}{an n by \code{ncol(xind)} matrix of function
  evaluations \eqn{X_i(s_{i1}),\dots, X_i(s_{iS})};
  \eqn{i=1,\dots,n}.}

  \item{yind}{matrix (or vector) of indices of evaluations
  of \eqn{Y_i(t)}; i.e. matrix with rows
  \eqn{(t_{i1},\dots,t_{iT})}}

  \item{xind}{matrix (or vector) of indices of evaluations
  of \eqn{X_i(s)}; i.e. matrix with rows
  \eqn{(s_{i1},\dots,s_{iS})}}

  \item{basistype}{defaults to "\code{\link[mgcv]{te}}",
  i.e. a tensor product spline to represent \eqn{f(X_i(s),
  t)}.  Alternatively, use \code{"s"} for bivariate basis
  functions (see \code{\link[mgcv]{s}}) or \code{"t2"} for
  an alternative parameterization of tensor product splines
  (see \code{\link[mgcv]{t2}}).}

  \item{integration}{method used for numerical integration.
  Defaults to \code{"simpson"}'s rule.  Alternatively and
  for non-equidistant grids, \code{"trapezoidal"}.}

  \item{L}{optional: an n by \code{ncol(xind)} giving the
  weights for the numerical integration over \eqn{s}.}

  \item{limits}{(NOT YET IMPLEMENTED)}

  \item{splinepars}{optional arguments supplied to the
  \code{basistype}-term. Defaults to a cubic tensor product
  B-spline with marginal second differences, i.e.
  \code{list(bs="ps", m=c(2,2,2))} See
  \code{\link[mgcv]{te}} or \code{\link[mgcv]{s}} for
  details}
}
\value{
  a list containing \itemize{ \item \code{call} a "call" to
  \code{\link[mgcv]{te}} (or \code{\link[mgcv]{s}},
  \code{\link[mgcv]{t2}}) using the appropriately
  constructed covariate and weight matrices (see
  \code{\link[mgcv]{linear.functional.terms}}) \item
  \code{data} a list containing the necessary covariate and
  weight matrices }
}
\description{
  Defines a term \eqn{\int^{s_{hi, i}}_{s_{lo, i}}
  f(X_i(s), s, t) ds} for inclusion in an
  \code{mgcv::gam}-formula (or \code{bam} or \code{gamm} or
  \code{gamm4:::gamm}) as constructed by
  \code{\link{pffr}}. Defaults to a cubic tensor product
  B-spline with marginal second differences penalties for
  \eqn{f(X_i(s), s, t)} and integration over the entire
  range \eqn{[s_{lo, i}, s_{hi, i}] = [\min(s_i),
  \max(s_i)]}. Can't deal with any missing \eqn{X(s)},
  unequal lengths of \eqn{X_i(s)} not (yet?) possible.
  Unequal ranges for different \eqn{X_i(s)} should work.
  \eqn{X_i(s)} is assumed to be numeric.\cr \code{sff()} IS
  AN EXPERIMENTAL FEATURE AND NOT WELL TESTED YET -- USE AT
  YOUR OWN RISK.
}
\author{
  Fabian Scheipl, based on Sonja Greven's trick for fitting
  functional responses.
}

