\name{wnet}
\alias{wnet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generalized elastic net in the wavelet domain
}
\description{
Performs generalized linear scalar-on-function or scalar-on-image regression in the wavelet domain, by  (naive) elastic net.
}
\usage{
wnet(y, xfuncs, min.scale, alpha, lambda = NULL, standardize = FALSE, covt = NULL, 
     pen.covt = FALSE, filter.number = 10, wavelet.family = "DaubLeAsymm",
     family = "gaussian", nfold = 5, nsplit = 1, compare.fits = FALSE, 
     store.cv = FALSE, seed = NULL, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
scalar outcome vector.}
  \item{xfuncs}{
functional predictors. For 1D predictors, an \eqn{n \times d} matrix of signals, where \eqn{n} is the length of \code{y} and \eqn{d} is the number of sites at which each signal is observed. For 2D predictors, an \eqn{n \times d \times d} array comprising \eqn{n} images of dimension \eqn{d \times d}. For 3D predictors, an \eqn{n \times d \times d \times d} array comprising \eqn{n} images of dimension \eqn{d \times d \times d}. Note that \eqn{d} must be a power of 2.
}
  \item{min.scale}{
either a scalar, or a vector of candidate values to be compared. Used to control the coarseness level of the wavelet decomposition.  The maximum value of \code{min.scale} is \eqn{log_2(d) - 1}.
}
  \item{alpha}{
elastic net mixing parameter, used by \code{\link[glmnet]{glmnet}}: either a scalar, or a vector of values to be compared. \code{alpha=1} gives the lasso penalty, while \code{alpha=0} yields the ridge penalty.
}
  \item{lambda}{
a vector of candidate regularization parameter values. If not supplied, \code{\link[glmnet]{glmnet}} automatically generates a sequence of candidate values.
}
  \item{standardize}{
logical, passed to \code{\link[glmnet]{glmnet}}: should the predictor variables be standardized?  Defaults to \code{FALSE}, but either way, the coefficients are returned on the original scale.
}
  \item{covt}{
covariates, if any: an \eqn{n}-row matrix, or a vector of length \eqn{n}.
}
  \item{pen.covt}{
logical: should the scalar covariates be penalized? If \code{FALSE} (the default), penalization is suppressed by setting the appropriate components of \code{penalty.factor} to 0 in the call to \code{\link[glmnet]{glmnet}}.
}


  \item{filter.number}{
passed to \code{\link[wavethresh]{wd}}, \code{\link[wavethresh]{imwd}}, or \code{\link[wavethresh]{wd3D}}, in the \pkg{wavethresh} package, to select the smoothness of the wavelets; defaults to \code{10}. 
}
  \item{wavelet.family}{
family of wavelets: passed to \code{\link[wavethresh]{wd}}, \code{\link[wavethresh]{imwd}}, or \code{\link[wavethresh]{wd3D}}. Defaults to \code{"DaubLeAsymm"}.
}
  \item{family}{
generalized linear model family. Current version supports \code{"gaussian"} (the default) and \code{"binomial"}.
}

  \item{nfold}{
the number of validation sets ("folds") into which the data are divided.
}
  \item{nsplit}{
the number of independent splits.
}
  \item{compare.fits}{
logical: if \code{TRUE}, the pairwise distances among training-set coefficient function estimates are computed.
}
  \item{store.cv}{
logical: should the output include a CV result table?
}
  \item{seed}{
the seed for random data division. If \code{seed = NULL}, a random seed is used.
}
\item{\dots}{
other arguments passed to \code{\link[glmnet]{glmnet}}.
}
}
\details{
This function supports only the standard discrete wavelet transform (see argument \code{type} in \code{\link[wavethresh]{wd}}) with periodic boundary handling (see argument \code{bc} in \code{\link[wavethresh]{wd}}).
}
\value{
A list with the following components:
\item{fitted.value, residuals}{the fitted values and residuals.}
\item{param.coef}{parametric coefficient estimates, for the scalar covariates.}
\item{fhat}{coefficient function estimate.}
\item{Rsq}{coefficient of determination.}
\item{min.scale,alpha}{the supplied value, or the one chosen by CV, for each of these arguments.}
\item{lambda}{the \code{lambda} chosen by CV (or the specified \code{lambda}, if only a single value was supplied).}
\item{cv.table}{if \code{store.cv = TRUE}, a table giving the CV criterion for each combination of \code{min.scale}, \code{alpha} and \code{lambda}. Otherwise, just the value for the optimized combination of \code{min.scale}, \code{alpha} and \code{lambda}.}
\item{stability}{if \code{compare.fits = TRUE}, a table of stability measures for each combination of \code{min.scale}, \code{alpha} and \code{lambda}. The stability measure is derived by computing the sample variance of the training-set coefficient function estimates at each point, then averaging over all points.}
}
\references{
Zhao, Y., Ogden, R. T., and Reiss, P. T. (2012). Wavelet-based LASSO in functional linear regression.  \emph{Journal of Computational and Graphical Statistics}, to appear.
}
\author{
Lan Huo \email{lan.huo@nyumc.org} and Yihong Zhao
}

\seealso{
\code{\link{wcr}}, \code{\link{fpcr}}
}
\examples{
\dontrun{
### 1D functional predictor example ###

data(gasoline)

# input a single value of each tuning parameters
gas.wnet1 <- wnet(gasoline$octane, xfuncs = gasoline$NIR[,1:256], 
             min.scale = 0, alpha = 1)

# input vectors of candidate tuning parameter values
gas.wnet2 <- wnet(gasoline$octane, xfuncs = gasoline$NIR[,1:256], 
                  min.scale = 0:3, alpha = c(0.9, 1))

### 2D functional predictor example ###

n = 200; d = 64

# Create true coefficient function
ftrue = matrix(0,d,d)
ftrue[40:46,34:38] = 1

# Generate random functional predictors, and scalar responses
ii = array(rnorm(n*d^2), dim=c(n,d,d))
iimat = ii; dim(iimat) = c(n,d^2)
yy = iimat \%*\% as.vector(ftrue) + rnorm(n, sd=.3)

mm.wnet <- wnet(yy, xfuncs = ii, min.scale = 4, alpha = 1)

image(ftrue)
contour(mm.wnet$fhat, add=TRUE)

### 3D functional predictor example ###

n = 200; d = 16

# Create true coefficient function
ftrue = array(0,dim = rep(d, 3))
ftrue[10:16,12:15, 4:8] = 1

# Generate random functional predictors, and scalar responses
ii = array(rnorm(n*d^3), dim=c(n,rep(d,3)))
iimat = ii; dim(iimat) = c(n,d^3)
yy = iimat \%*\% as.vector(ftrue) + rnorm(n, sd=.3)

mmm.wnet <- wnet(yy, xfuncs = ii, min.scale = 2, alpha = 1)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{wavelet}  % __ONLY ONE__ keyword per line
