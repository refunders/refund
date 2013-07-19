\name{wcr/wnet.perm}
\alias{wnet.perm}
\alias{wcr.perm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Permutation test for wavelet-domain scalar-on-function regression
}
\description{
This function assesses statistical significance of a \code{\link{wcr}} or  \code{\link{wnet}} fit by referring the cross-validation criterion to a permutation distribution.
}
\usage{
wcr.perm(y, xfuncs, min.scale = 0, nfeatures, ncomp, method = c("pcr", "pls"), 
         covt = NULL, nrep = 1, nsplit = 1, nfold = 5, nperm = 20, 
         perm.method = c('responses', 'y.residuals', 'x.residuals'), 
         family = "gaussian", seed.real = NULL, seed.perm = NULL, ...)
wnet.perm(y, xfuncs, min.scale = 0, nfeatures = NULL, alpha = 1, lambda, 
         covt = NULL, nrep = 1, nsplit = 1, nfold = 5, nperm = 20, 
         perm.method = c('responses', 'y.residuals', 'x.residuals'), 
         family = "gaussian", seed.real = NULL, seed.perm = NULL, ...)         
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y, xfuncs, min.scale, nfeatures, method, covt, family, nsplit, nfold}{arguments passed to \code{\link{wcr}} or \code{\link{wnet}}.}
  \item{ncomp}{number of components; passed to \code{\link{wcr}}.}
  \item{alpha, lambda}{tuning parameters, passed to \code{\link{wnet}}.}
  \item{nrep}{number of replicates to average over, when computing the real-data CV criterion.}
  \item{nperm}{
number of permutations. The default is suitable for toy applications only.
}
  \item{perm.method}{
one of
\itemize{
\item {\code{"responses"}: permute the response vector \code{y}.}
\item {\code{"y.residuals"}: permute the residuals upon regressing \code{y ~ covt}.}
\item {\code{"x.residuals"}: permute the residuals upon regressing \code{xfuncs ~ covt}.}
}
See Details.}
  \item{seed.real}{
the seed for random data division for real data. If \code{seed.real = NULL}, a random seed is used.
}
  \item{seed.perm}{
the seed for random data division for permuted data. If \code{seed.perm = NULL}, a random seed is used.
}
  \item{seed.real}{
the seed for random data division for real data. If \code{seed.real = NULL}, a random seed is used.
}
  \item{seed.perm}{
the seed for random data division for permuted data. If \code{seed.perm = NULL}, a random seed is used.
}
  \item{\dots}{other arguments passed to \code{\link{wcr}} or \code{\link{wnet}}.}
}
\details{
	Permutation tests of this type are discussed, in a classification setting, by Ojala and Garriga (2010). Permuting the responses is appropriate when regressing on functions/images only, with no scalar covariates. For linear regression with covariates, it is preferable to first regress on the covariates, and then permute the residuals. For logistic regression this is not feasible; but, following Potter (2005), one can instead permute the residuals from a regression of the functions/images on the covariates.
}
\value{
  \item{cv}{CV value for the real data (averaged over \code{nrep} replications).
}
  \item{cv.perm}{CV values for the permuted data.
}
  \item{pvalue}{p-value for the permutation test.
}
}
\references{
		Ojala, M., and Garriga, G. C. (2010). Permutation tests for studying classifier performance. \emph{Journal of Machine Learning Research}, 11, 1833--1863.
	
	Potter, D. M. (2005). A permutation test for inference in logistic regression with small- and moderate-sized data sets. \emph{Statistics in Medicine}, 24, 693--708.
	}
\author{
Lan Huo \email{lan.huo@nyumc.org}
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{wcr}}, \code{\link{wnet}}
}
\examples{

n = 200; d = 64

# Create true coefficient function
ftrue = matrix(0,d,d)
ftrue[40:46,34:38] = 1

# Generate random functional predictors, and scalar responses
ii = array(rnorm(n*d^2), dim=c(n,d,d))
iimat = ii; dim(iimat) = c(n,d^2)
yy = iimat \%*\% as.vector(ftrue) + rnorm(n, sd=.3)

obj.wnet.perm <- wnet.perm(yy, xfuncs = ii, min.scale = 4, nfeatures = 200, alpha = 1, 
                           nperm = 10)

obj.wcr.perm <- wcr.perm(yy, xfuncs = ii, min.scale = 4, nfeatures = 20, ncomp = 6, 
                         cv1 = TRUE, method = "pls", nperm = 10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
