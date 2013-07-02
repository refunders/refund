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
    cv1 = FALSE, nfold = 5, nsplit = 1, store.cv = FALSE, 
    store.glm = TRUE, seed = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
scalar outcome vector.
}
  \item{xfuncs}{
functional predictors. For 1D predictors, an \eqn{n \times d} matrix of signals, where \eqn{n} is the length of \code{y} and \eqn{d} is the number of sites at which each signal is defined. For 2D predictors, an \eqn{n \times d \times d} array comprising \eqn{n} images of dimension \eqn{d \times d}. For 3D predictors, an \eqn{n \times d \times d \times d} array comprising \eqn{n} images of dimension \eqn{d \times d \times d}. Note that \eqn{d} must be a power of \eqn{2}.}
  \item{min.scale}{
either a scalar, or a vector of values to be compared. Used to control the coarseness level of wavelet decomposition.  Possible values are \eqn{0,1,\dots,log_2(d) - 1}.
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
logical: should cross-validation be performed (to estimate prediction error) even if a single value is provided for each of \code{min.scale}, \code{nfeatures} and \code{ncomp}? By default, \code{FALSE}. Note that whenever multiple candidate values are provided for one or more of these tuning parameters, CV is performed to select the best model. 
}
  \item{nfold}{
the number of validation sets ("folds") into which the data are divided.
}
  \item{nsplit}{
number of splits into \code{nfold} validation sets; CV is computed by averaging over these splits.
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
Briefly, the algorithm works by (1) applying the discrete wavelet transform (DWT) to the functional/image predictors; (2) retaining only the \code{nfeatures} wavelet coefficients having the highest variance (for PCR; cf. Johnstone and Lu, 2009) or highest covariance with \code{y} (for PLS); (3) regressing \code{y} on the leading \code{ncomp} PCs or PLS components, along with any scalar covariates; and (4) applying the inverse DWT to the result to obtain the coefficient function estimate \code{fhat}.

This function supports only the standard DWT (see argument \code{type} in \code{\link[wavethresh]{wd}}) with periodic boundary handling (see argument \code{bc} in \code{\link[wavethresh]{wd}}).

For 2D predictors, setting \code{min.scale=1} will lead to an error, due to a technical detail regarding \code{\link[wavethresh]{imwd}}. Please contact the author if a workaround is needed.

See the Details for \code{\link{fpcr}} for a note regarding decorrelation.
}
\value{A list that, if \code{store.glm = TRUE}, includes all components of the fitted \code{\link{glm}} object. The following components are included even if \code{store.glm = FALSE}:
\item{fitted}{the fitted values.}
\item{param.coef}{coefficients for covariates with decorrelation. The model is fitted after decorrelating the functional predictors from any scalar covariates; but for CV, one needs the "undecorrelated" coefficients from the training-set models.}
\item{undecor.coef}{coefficients for covariates \emph{without} decorrelation. See \code{param.coef}.}
\item{fhat}{coefficient function estimate.}
\item{tuning.params}{\eqn{2 \times 4} table the indices and values of \code{min.scale}, \code{nfeatures} and \code{ncomp} chosen by CV.}
\item{cv.table}{a table giving the CV criterion for each combination of \code{min.scale}, \code{nfeatures} and \code{ncomp}, if \code{store.cv = TRUE}; otherwise, the CV criterion only for the optimized combination of these parameters.  Set to \code{NULL} if CV is not performed.}
\item{se.cv}{if \code{store.cv = TRUE}, the standard error of the CV estimate for each combination of \code{min.scale}, \code{nfeatures} and \code{ncomp}.}
\item{Rsq}{coefficient of determination.}
\item{family}{generalized linear model family.}
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
% \dontrun{
# See example for wnet

% }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{wavelet}
