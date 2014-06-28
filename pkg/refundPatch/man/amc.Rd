\name{amc}
\Rdversion{1.1}
\alias{amc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Additive model with constraints}
\description{
An internal function, called by \code{fosr()}, that fits additive models with linear constraints via a call to \code{\link[mgcv]{gam}} or \code{\link[mgcv]{bam}}
in the \pkg{mgcv} package.}
\usage{
amc(y, Xmat, S, gam.method = "REML", C = NULL, lambda = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
response vector.}
  \item{Xmat}{
design matrix.}
  \item{S}{
list of penalty matrices.}
  \item{gam.method}{
smoothing parameter selection method: "REML" for restricted maximum likelihood, "GCV.Cp" for generalized cross-validation.}
  \item{C}{
matrix of linear constraints.  Dimension should be number of constraints times \code{ncol(Xmat)}.}
  \item{lambda}{
smoothing parameter value.  If \code{NULL}, the smoothing parameter(s) will be estimated.}
\item{\dots}{other arguments, passed to \code{\link[mgcv]{gam}} or \code{\link[mgcv]{bam}}.}
}
\details{
The additive model is fitted using \code{\link[mgcv]{gam}}, unless there are more than 10000 responses; in that case \code{\link[mgcv]{bam}} is used.}
\value{
A list with the following elements:
\item{gam}{the \code{gam} object returned by \code{gam} or \code{bam}.}
\item{coefficients}{coefficients with respect to design matrix \code{Xmat}, derived from the \code{gam()} fit.}
\item{Vp, GinvXt}{outputs used by \code{fosr}.}
\item{method}{the \code{gam.method} argument of the call to \code{amc}.}
}

\author{
Philip Reiss \email{phil.reiss@nyumc.org}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{fosr}}
}

