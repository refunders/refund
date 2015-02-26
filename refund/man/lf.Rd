\name{lf}
\Rdversion{1.0}
\alias{lf}
\title{
Construct an FLM regression term}
\description{
 Defines a term \eqn{\int_{T}\beta(t)X_i(t)dt} for inclusion in an \code{mgcv::gam-formula} (or \code{bam} or \code{gamm} or \code{gamm4:::gamm}) as constructed by \code{\link{fgam}}, where \eqn{\beta(t)} is an unknown coefficient function and \eqn{X_i(t)} is a functional predictor on the closed interval \eqn{T}. Defaults to a cubic B-spline with second-order difference penalties for estimating \eqn{\beta(t)}.  The functional predictor must be fully observed on a regular grid.}
\usage{
lf(X, xind = seq(0, 1, l = ncol(X)), 
   integration = c("simpson", "trapezoidal", "riemann"), 
   L = NULL, 
   splinepars = list(bs = "ps", k= min(ceiling(n/4), 40), m = c(2, 2)), 
   presmooth = TRUE)
}
\arguments{
  \item{X}{
an \code{N} by \code{J=ncol(xind)} matrix of function evaluations \eqn{X_i(t_{i1}),., X_i(t_{iJ}); i=1,.,N.}}
  \item{xind}{
matrix (or vector) of indices of evaluations of \eqn{X_i(t)}; i.e. a matrix with \emph{i}th row \eqn{(t_{i1},.,t_{iJ})}.}
\item{integration}{method used for numerical integration. Defaults to "\code{simpson}"'s rule for calculating entries in \code{L}. Alternatively and for non-equidistant grids, "\code{trapezoidal}" or "\code{riemann}". "\code{riemann}" integration is always used if \code{limits} is specified.}
  \item{L}{
   optional: an \code{N} by \code{ncol(xind)} matrix giving the weights for the numerical integration over \code{t}.}
     \item{splinepars}{
   optional arguments specifying options for representing and penalizing the functional coefficient \eqn{\beta(t)}. Defaults to a cubic B-spline with second-order difference penalties, i.e. \code{list(bs="ps", m=c(2, 1))} See \code{\link{te}} or \code{\link{s}} in \code{mgcv} for details.}
     \item{presmooth}{
   If true, the functional predictor is pre-smoothed prior to fitting.  See \code{\link{smooth.basisPar}} in package \code{\link{fda}}}
}

\value{
a list with the following entries
\item{call}{a "call" to \code{te} (or \code{s}, \code{t2}) using the appropriately constructed covariate and weight matrices.}
\item{xind}{ the \code{xind} argument supplied to \code{lf}.}
\item{L}{ the  matrix of weights used for the integration.}
\item{xindname}{ the name used for the functional predictor variable in the \code{formula} used by \code{mgcv}.}
\item{tindname}{the name used for \code{xind} variable in the \code{formula} used by \code{mgcv}.}
\item{LXname}{ the name used for the \code{L} variable in the \code{formula} used by \code{mgcv}.}
\item{presmooth}{the \code{presmooth} argument supplied to \code{lf}.}
\item{Xfd}{an \code{fd} object from presmoothing the functional predictors using \code{\link{smooth.basisPar}}.  Only present if \code{presmooth=TRUE}.  See \code{\link{fd}}.}
}
\author{Mathew W. McLean \email{mwm79@cornell.edu} and Fabian Scheipl}
\seealso{
\code{\link{fgam}}, \code{\link{af}}, mgcv's \code{\link{linear.functional.terms}}, \code{\link{fgam}} for examples.}
