\name{pre-/post-process}
\alias{Data2wd}
\alias{wd2fhat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generic functions for transforming scalar-on-function regression to the wavelet domain and back}
\description{
These two functions enable the user to transform the predictors to the wavelet domain (\code{Data2wd}); apply his/her favorite high-dimensional regression method to the wavelet coefficients (see the example); and inverse-transform the result to obtain a coefficient function estimate in the original domain (\code{wd2fhat}).}
\usage{
Data2wd(y, xfuncs, covt = NULL, min.scale = 0, nfeatures = NULL, filter.number = 10,
 wavelet.family = "DaubLeAsymm")

wd2fhat(est, info)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y, xfuncs, covt,
  	min.scale, nfeatures,
  	filter.number, wavelet.family}{
  		see \code{\link{wnet}}.
}
  \item{est}{
an estimate produced by regressing on the wavelet-transformed predictors.
}
  \item{info}{
the wavelet infomation, output from \code{Data2wd}.
}
}
\value{
\code{Data2wd} returns a list with components:
\item{X}{the design matrix, with columns of covariates and functional predictors converted to wavelet domain.}
\item{info}{the infomation needed to reconstruct coefficient functions}
 
\code{wd2fhat} returns the coefficient function. It could be a curve or a 2D/3D image.
}

\author{
Lan Huo \email{lan.huo@nyumc.org} and Philip Reiss \email{phil.reiss@nyumc.org}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[wavethresh]{wd}}
}
\examples{
# MCP-penalized regression via ncvreg (which can also apply a SCAD penalty)
data(gasoline)
res <- Data2wd(gasoline$octane, xfuncs = gasoline$NIR[,1:256])
obje = ncvreg::cv.ncvreg(res$X, gasoline$octane)
est = obje$fit$beta[,which.min(obje$cve)]
names(est) <- rownames(obje$fit$beta)
beta <- wd2fhat(est, res$info)
}
