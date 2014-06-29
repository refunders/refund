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
wnet(y, xfuncs, covt = NULL, min.scale = 0, nfeatures = NULL, alpha = 1, 
     lambda = NULL, standardize = FALSE, pen.covt = FALSE, 
     filter.number = 10, wavelet.family = "DaubLeAsymm", family = "gaussian", 
     nfold = 5, nsplit = 1, store.cv = FALSE, store.glmnet = FALSE, 
     seed = NULL, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
scalar outcome vector.}
  \item{xfuncs}{
functional predictors. For 1D predictors, an \eqn{n \times d} matrix of signals, where \eqn{n} is the length of \code{y} and \eqn{d} is the number of sites at which each signal is observed. For 2D predictors, an \eqn{n \times d \times d} array comprising \eqn{n} images of dimension \eqn{d \times d}. For 3D predictors, an \eqn{n \times d \times d \times d} array comprising \eqn{n} images of dimension \eqn{d \times d \times d}. Note that \eqn{d} must be a power of 2.
}
  \item{covt}{
covariates, if any: an \eqn{n}-row matrix, or a vector of length \eqn{n}.
}
  \item{min.scale}{
either a scalar, or a vector of candidate values to be compared. Used to control the coarseness level of the wavelet decomposition.  Possible values are \eqn{0,1,\dots,log_2(d) - 1}.
}
  \item{nfeatures}{
number(s) of features, i.e. wavelet coefficients, to retain for prediction of \code{y}: either a scalar, or a vector of values to be compared.    
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
  \item{pen.covt}{
logical: should the scalar covariates be penalized? If \code{FALSE} (the default), penalization is suppressed by setting the appropriate components of \code{penalty.factor} to 0 in the call to \code{\link[glmnet]{glmnet}}.
}
  \item{filter.number}{
passed to \code{\link[wavethresh]{wd}}, \code{\link[wavethresh]{imwd}}, or \code{\link[wavethresh]{wd3D}}, in the \pkg{wavethresh} package, to select the smoothness of the wavelets. 
}
  \item{wavelet.family}{
family of wavelets: passed to \code{\link[wavethresh]{wd}}, \code{\link[wavethresh]{imwd}}, or \code{\link[wavethresh]{wd3D}}.
}
  \item{family}{
generalized linear model family. Current version supports \code{"gaussian"} (the default) and \code{"binomial"}.
}

  \item{nfold}{
the number of validation sets ("folds") into which the data are divided.
}
  \item{nsplit}{
number of splits into \code{nfold} validation sets; CV is computed by averaging over these splits.
}
%  \item{decorrelate}{
%logical: if \code{TRUE},  \code{xfuncs} is regressed on \code{covt}, and %the residuals  are used instead of the original \code{xfuncs} as the %functional predictors.
%}
  \item{store.cv}{
logical: should the output include a CV result table?
}
  \item{store.glmnet}{
logical: should the output include the fitted \code{\link[glmnet]{glmnet}}? }
  \item{seed}{
the seed for random data division. If \code{seed = NULL}, a random seed is used.
}
\item{\dots}{
other arguments passed to \code{\link[glmnet]{glmnet}}.
}
}
\details{
This function supports only the standard discrete wavelet transform (see argument \code{type} in \code{\link[wavethresh]{wd}}) with periodic boundary handling (see argument \code{bc} in \code{\link[wavethresh]{wd}}).

For 2D predictors, setting \code{min.scale=1} will lead to an error, due to a technical detail regarding \code{\link[wavethresh]{imwd}}. Please contact the authors if a workaround is needed.

}
\value{
An object of class \code{"wnet"}, which is a list with the following components:
\item{glmnet}{if \code{store.glmnet = TRUE}, the object returned by \code{\link[glmnet]{glmnet}}.}
\item{fitted.value}{the fitted values.}
\item{coef.param}{parametric coefficient estimates, for the scalar covariates.}
\item{fhat}{coefficient function estimate.}
\item{Rsq}{coefficient of determination.}
\item{lambda.table}{array giving the candidate lambda values, chosen automatically by \code{\link[glmnet]{glmnet}}, for each combination of the other tuning parameters.}
\item{tuning.params}{a \eqn{2\times 4} table giving the indices and values of \code{min.scale}, \code{nfeatures}, \code{alpha} and \code{lambda} that minimize the CV. For example, if \code{alpha=c(0,0.5,1)} is specified and the CV-minimizing tuning parameter combination takes \code{alpha} to be the 2nd of these values, then the third column of the table is \code{c(2, 0.5)}.}
\item{cv.table}{if \code{store.cv = TRUE}, a table giving the CV criterion for each combination of \code{min.scale}, \code{nfeatures}, \code{alpha} and \code{lambda}. Otherwise, just the minimized CV criterion.}
\item{se.cv}{if \code{store.cv = TRUE}, the standard error of the CV estimate for each combination of \code{min.scale}, \code{nfeatures}, \code{alpha} and \code{lambda}.}
\item{family}{generalized linear model family.}
}
\references{
Zhao, Y., Ogden, R. T., and Reiss, P. T. (2012). Wavelet-based LASSO in functional linear regression.  \emph{Journal of Computational and Graphical Statistics}, 21(3), 600--617.
}
\author{
Lan Huo \email{lan.huo@nyumc.org} and Yihong Zhao
}

\seealso{
\code{\link{wcr}}, \code{\link{wnet.perm}}
}
\examples{
### 1D functional predictor example ###

data(gasoline)

# input a single value of each tuning parameters
gas.wnet1 <- wnet(gasoline$octane, xfuncs = gasoline$NIR[,1:256], 
             nfeatures= 20, min.scale = 0, alpha = 1)
gas.wpcr1 <- wcr(gasoline$octane, xfuncs = gasoline$NIR[,1:256], min.scale = 0,
                 nfeatures = 20, ncomp = 15)
gas.wpls1 <- wcr(gasoline$octane, xfuncs = gasoline$NIR[,1:256], min.scale = 0, 
                 nfeatures = 20, ncomp = 15, method = "pls")
plot(gas.wnet1)
plot(gas.wpcr1)
plot(gas.wpls1)

# input vectors of candidate tuning parameter values
gas.wnet2 <- wnet(gasoline$octane, xfuncs = gasoline$NIR[,1:256], 
                  nfeatures= 20, min.scale = 0:3, alpha = c(0.9, 1))
gas.wpcr2 <- wcr(gasoline$octane, xfuncs = gasoline$NIR[,1:256], min.scale = 0:3,
                 nfeatures = c(16, 18, 20), ncomp = 10:15)
gas.wpls2 <- wcr(gasoline$octane, xfuncs = gasoline$NIR[,1:256], min.scale = 0:3,
                 nfeatures = c(16, 18, 20), ncomp = 10:15, method = "pls")
plot(gas.wnet2)
plot(gas.wpcr2)
plot(gas.wpls2)

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

mm.wpls <- wcr(yy, xfuncs = ii, min.scale = 4, nfeatures = 20, ncomp = 6,
                method = "pls")

plot(mm.wnet)
plot(mm.wpls)

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

mmm.wpls <- wcr(yy, xfuncs = ii, min.scale = 2, nfeatures = 20, ncomp = 6,
                method = "pls")
plot(mmm.wnet)
plot(mmm.wpls)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{wavelet}  % __ONLY ONE__ keyword per line
