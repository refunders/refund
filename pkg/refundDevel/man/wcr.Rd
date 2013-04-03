\name{wcr}
\alias{wcr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Principal component regression and partial least squares in the wavelet domain
}
\description{
Performs generalized linear scalar-on-function or scalar-on-image regression in the wavelet domain, by sparse principal component regression  (PCR) and sparse partial least squares (PLS).
}
\usage{
wcr(y, xfuncs, min.scale, nfeatures, ncomp, method = c("pcr", "pls"), 
    mean.signal.term = FALSE, covt = NULL, filter.number = 10, 
    wavelet.family = "DaubLeAsymm", family = "gaussian", 
    cv1 = FALSE, nfold = 5, nsplit = 1, compare.fits = FALSE, 
    store.cv = FALSE, store.glm = TRUE, seed = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
scalar outcome vector.
}
  \item{xfuncs}{
functional predictors. For 1D predictors, an \eqn{n \times d} matrix of signals, where \eqn{n} is the length of \code{y} and \eqn{d} is the number of sites at which each signal is defined. For 2D predictors, an \eqn{n \times d \times d} array comprising \eqn{n} images of dimension \eqn{d \times d}. For 3D predictors, an \eqn{n \times d \times d \times d} array comprising \eqn{n} images of dimension \eqn{d \times d \times d}. Note that \eqn{d} must be a power of \eqn{2}.}
  \item{min.scale}{
either a scalar, or a vector of values to be compared. Used to control the coarseness level of wavelet decomposition.  The maximum value of \code{min.scale} is \eqn{log_2(d) - 1}.
}
  \item{nfeatures}{
number(s) of features, i.e. wavelet coefficients, to retain for prediction of \code{y}: either a scalar, or a vector of values to be compared.
}
  \item{ncomp}{
number(s) of principal components (if \code{method="pcr"}) or PLS components (if \code{method="pcr"}): either a scalar, or a vector of values to be compared.
}
  \item{method}{
either "\code{pcr}" (principal component regression) (the default) or "\code{pls}" (partial least squares).
}
  \item{mean.signal.term}{
logical: should the mean of each subject's signal be included as a covariate? By default, \code{FALSE}.
}
  \item{covt}{
covariates: an \eqn{n}-row matrix, or a vector of length \eqn{n}; defaults to \code{NULL}.
}
  \item{filter.number}{
argument passed to function \code{\link[wavethresh]{wd}}, \code{\link[wavethresh]{imwd}}, or \code{\link[wavethresh]{wd3D}} in the \pkg{wavethresh} package. Used to select the smoothness of wavelet in the decomposition; defaults to \code{10}. 
}
  \item{wavelet.family}{
family of wavelets: passed to functions \code{\link[wavethresh]{wd}}, \code{\link[wavethresh]{imwd}}, or\code{\link[wavethresh]{wd3D}}; defaults to \code{"DaubLeAsymm"}.
}
  \item{family}{
generalized linear model family. Current version supports \code{"gaussian"} (the default) and \code{"binomial"}.
}
  \item{cv1}{
logical: should cross-validation be performed (to estimate prediction error) even if a single value is provided for each of \code{min.scale}, \code{nfeatures} and \code{ncomp}? By default, \code{FALSE}. Note that if multiple candidate values are provided for at least one of these tuning parameters, CV is always performed to select the best model. 
}
  \item{nfold}{
the number of validation sets ("folds") into which the data are divided for a split.
}
  \item{nsplit}{
the number of independent splits.
}
  \item{compare.fits}{
logical: if \code{TRUE}, the pairwise differences among training-set coefficient function estimates are computed.
}
  \item{store.cv}{
logical: should the output include a CV result table?
}
  \item{store.glm}{
logical: should the output include the fitted \code{\link{glm}}? Defaults to \code{TRUE}, but may need to be set to \code{FALSE} for simulations.
}
  \item{seed}{
the seed for random data division. If \code{seed = NULL}, a random seed is used.
}
}
\details{
Essentially, the algorithm works by (1) applying the discrete wavelet transform (DWT) to the predictors; (2) retaining only the \code{nfeatures} wavelet coefficients having the highest variance (for PCR; cf. Johnstone and Lu, 2009) or highest covariance with \code{y} (for PLS); (3) regressing \code{y} on the leading \code{ncomp} PCs or PLS components; and (4) applying the inverse DWT to the result to obtain the coefficient function estimate \code{fhat}.

This function supports only the standard DWT (see argument \code{type} in \code{\link[wavethresh]{wd}}) with periodic boundary handling (see argument \code{bc} in \code{\link[wavethresh]{wd}}).

For 2D predictors, setting \code{min.scale=1} will lead to an error, due to a technical detail regarding \code{\link[wavethresh]{imwd}}. Please contact the author if a workaround is needed.

See the Details for \code{\link{fpcr}} for a note regarding decorrelation.
}
\value{
\item{glm}{if \code{store.glm = TRUE}, the fitted \code{\link{glm}} object.}
\item{fhat}{coefficient function estimate.}
\item{undecor.coef}{coefficients for covariates \emph{without} decorrelation. The model is fitted after decorrelating the functional predictors from any scalar covariates; but for CV, one needs the "undecorrelated" coefficients from the training-set models.}
\item{min.scale, nfeatures, ncomp}{the values of these arguments specified in the call, if a single value was supplied; otherwise, the values chosen by CV.}
\item{cv.table}{a table giving the CV criterion for each combination of \code{min.scale}, \code{nfeatures} and \code{ncomp}, if \code{store.cv = TRUE}; otherwise, the CV criterion only for the optimized combination of these parameters.  Set to \code{NULL} if CV is not performed.}
\item{stability}{if \code{compare.fits = TRUE}, a table of stability measures for each combination of \code{min.scale}, \code{nfeatures} and \code{ncomp}. The stability measure is derived by computing the sample variance of the training-set coefficient function estimates at each point, then averaging over all points.}
%% ...
}
\references{
Johnstone, I. M., and Lu, Y. (2009). On consistency and sparsity for principal components analysis in high dimensions. \emph{Journal of the American Statistical Association}, 104, 682--693.
}
\author{
Lan Huo \email{lan.huo@nyumc.org}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{wnet}}, \code{\link{fpcr}}
}
\examples{
\dontrun{
### 1D functional predictor example ###
data(gasoline)

# input a single value of each tuning parameter
gas.wpcr1 <- wcr(gasoline$octane, xfuncs = gasoline$NIR[,1:256], min.scale = 0,
                 nfeatures = 20, ncomp = 15)
gas.wpls1 <- wcr(gasoline$octane, xfuncs = gasoline$NIR[,1:256], min.scale = 0, 
                 nfeatures = 20, ncomp = 15, method = "pls")

# input vectors of tuning parameter values
gas.wpcr2 <- wcr(gasoline$octane, xfuncs = gasoline$NIR[,1:256], min.scale = 0:3,
                 nfeatures = c(16, 18, 20), ncomp = 10:15)
gas.wpls2 <- wcr(gasoline$octane, xfuncs = gasoline$NIR[,1:256], min.scale = 0:3,
                 nfeatures = c(16, 18, 20), ncomp = 10:15, method = "pls")

### 2D functional predictor example ###

n = 200; d = 64

# Create true coefficient function
ftrue = matrix(0,d,d)
ftrue[40:46,34:38] = 1

# Generate random functional predictors, and scalar responses
ii = array(rnorm(n*d^2), dim=c(n,d,d))
iimat = ii; dim(iimat) = c(n,d^2)
yy = iimat \%*\% as.vector(ftrue) + rnorm(n, sd=.3)

mm.wpls <- wcr(yy, xfuncs = ii, min.scale = 4, nfeatures = 20, ncomp = 6,
                method = "pls")
image(ftrue)
contour(mm.wpls$fhat, add=TRUE)

### 3D functional predictor example ###

n = 200; d = 16

# Create true coefficient function
ftrue = array(0,dim = rep(d, 3))
ftrue[10:16,12:15, 4:8] = 1

# Generate random functional predictors, and scalar responses
ii = array(rnorm(n*d^3), dim=c(n,rep(d,3)))
iimat = ii; dim(iimat) = c(n,d^3)
yy = iimat \%*\% as.vector(ftrue) + rnorm(n, sd=.3)

mmm.wpls <- wcr(yy, xfuncs = ii, min.scale = 2, nfeatures = 20, ncomp = 6,
                method = "pls")

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{wavelet}
