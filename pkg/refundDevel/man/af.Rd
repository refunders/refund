\name{af}
\Rdversion{1.0}
\alias{af}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Construct an FGAM regression term}
\description{
   Defines a term \eqn{\int_{T}F(X_i(t),t)dt} for inclusion in an \code{mgcv::gam-formula} (or \code{bam} or \code{gamm} or \code{gamm4:::gamm}) as constructed by \code{\link{fgam}}, where \eqn{F(x,t)}$ is an unknown smooth bivariate function and \eqn{X_i(t)} is a functional predictor on the closed interval \eqn{T}. Defaults to a cubic tensor product B-spline with marginal second-order difference penalties for estimating \eqn{F(x,t)}.  The functional predictor must be fully observed on a regular grid.}
\usage{
af(X, xind = seq(0, 1, l = ncol(X)), basistype = c("te","t2", "s"), 
               integration = c("simpson", "trapezoidal", "riemann"), 
                L = NULL, splinepars = list(bs = "ps", k= rep(min(ceiling(sqrt(nrow(X))),20),2),
                m = list(c(2, 2), c(2, 2))), presmooth = TRUE, Xrange=range(X), Qtransform=FALSE)
}
\arguments{
  \item{X}{
an \code{N} by \code{J=ncol(xind)} matrix of function evaluations \eqn{X_i(t_{i1}),., X_i(t_{iJ}); i=1,.,N.}}
  \item{xind}{
matrix (or vector) of indices of evaluations of \eqn{X_i(t)}; i.e. a matrix with \emph{i}th row \eqn{(t_{i1},.,t_{iJ})}.}
\item{basistype}{
defaults to "\code{te}", i.e. a tensor product spline to represent \eqn{F(x,t)} Alternatively, use "\code{s}" for bivariate basis functions (see mgcv's \code{\link{s}}) or "\code{t2}" for an alternative parameterization of tensor product splines (see mgcv's \code{\link{t2}}).
}
\item{integration}{method used for numerical integration. Defaults to "\code{simpson}"'s rule for calculating entries in \code{L}. Alternatively and for non-equidistant grids, "\code{trapezoidal}" or "\code{riemann}". "\code{riemann}" integration is always used if \code{limits} is specified.}
  \item{L}{
   optional: an \code{N} by \code{ncol(xind)} matrix giving the weights for the numerical integration over \code{t}.}
     \item{splinepars}{
   optional arguments specifying options for representing and penalizing the function \eqn{F(x,t)}. Defaults to a cubic tensor product B-spline with marginal second-order difference penalties, i.e. \code{list(bs="ps", m=list(c(2, 2), c(2, 2))}, see \code{\link{te}} or \code{\link{s}} in \code{mgcv} for details.}
     \item{presmooth}{
   If true, the functional predictor is pre-smoothed prior to fitting.  See \code{\link{smooth.basisPar}} in package \code{\link{fda}}}
     \item{Xrange}{
   Range to use when specifying the marginal basis for the \emph{x}-axis.  It may be desired to increase this slightly over the default of \code{range(X)} if concerned about predicting for future observed curves that take values outside of \code{range(X)}.}
   \item{Qtransform}{Should the functional be transformed using the empirical cdf and applying a quantile transformation on each column of \code{X} prior to fitting?  This ensures \code{Xrange=c(0,1)}.  If \code{Qtransform=TRUE} and \code{presmooth=TRUE}, presmoothing is done prior to transforming the functional predictor.}
}

\value{
a list with the following entries
\item{call}{a "call" to \code{te} (or \code{s}, \code{t2}) using the appropriately constructed covariate and weight matrices.}
\item{xind}{ the \code{xind} argument supplied to \code{af}.}
\item{L}{ the  matrix of weights used for the integration.}
\item{xindname}{ the name used for the functional predictor variable in the \code{formula} used by \code{mgcv}.}
\item{tindname}{the name used for \code{xind} variable in the \code{formula} used by \code{mgcv}.}
\item{Lname}{ the name used for the \code{L} variable in the \code{formula} used by \code{mgcv}.}
\item{presmooth}{the \code{presmooth} argument supplied to \code{af}.}
\item{Qtranform}{the \code{Qtransform} argument supplied to \code{af}.}
\item{Xrange}{the \code{Xrange} argument supplied to \code{af}.}
\item{ecdflist}{ a list containing one empirical cdf function from applying \code{\link{ecdf}} to each (possibly presmoothed) column of \code{X}.  Only present if \code{Qtransform=TRUE}.}
\item{Xfd}{an \code{fd} object from presmoothing the functional predictors using \code{\link{smooth.basisPar}}.  Only present if \code{presmooth=TRUE}.  See \code{\link{fd}}.}
}
\author{Mathew W. McLean \email{mwm79@cornell.edu} and Fabian Scheipl}
\references{
McLean M.W. et al. (2012) \emph{Journal of Computational and Graphical Statistics}, to appear.}
\seealso{
\code{\link{fgam}}, \code{\link{lf}}, mgcv's \code{\link{linear.functional.terms}}, \code{\link{fgam}} for examples.}
\keyword{methods}
\keyword{models}
\keyword{smooth}
