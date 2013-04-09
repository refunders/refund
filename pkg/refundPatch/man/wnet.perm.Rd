\name{wnet.perm}
\alias{wnet.perm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Permutation test for (naive) elastic net in the wavelet domain
}
\description{
This function assesses statistical significance of a \code{\link{wnet}} fit by referring the cross-validation criterion to a permutation distribution.
}
\usage{
wnet.perm(y, xfuncs, min.scale, nfeatures, alpha, lambda, covt = NULL, nsplit = 10, nperm = 20, 
         perm.method = c('responses', 'y.residuals', 'x.residuals'), 
         family = "gaussian", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y, xfuncs, min.scale, nfeatures, alpha, lambda, covt, family}{arguments passed to \code{\link{wnet}}.}
  \item{nsplit}{
number of splits to estimate the true cv.
}
  \item{nperm}{
number of permutations. The default is meant for toy applications.
}
  \item{perm.method}{
one of
\itemize{
\item {\code{"responses"}: permute the response vector \code{y}.}
\item {\code{"y.residuals"}: permute the residuals upon regressing \code{y ~ covt}.}
\item {\code{"x.residuals"}: permute the residuals upon regressing \code{xfuncs ~ covt}.}
}
See Details.}
  \item{\dots}{other parameters passed to \code{\link{wnet}}.}
}
\details{
see \code{\link{wcr.perm}}
}
\value{
  \item{cv}{ the CV criterion for the real data.
}
  \item{cv.perm}{ the CV values for the permuted data.
}
  \item{pvalue}{ the p-value for the permutation test.
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
\code{\link{wnet}}, \code{\link{wcr.perm}}
}
\examples{
\dontrun{
n = 200; d = 64

# Create true coefficient function
ftrue = matrix(0,d,d)
ftrue[40:46,34:38] = 1

# Generate random functional predictors, and scalar responses
ii = array(rnorm(n*d^2), dim=c(n,d,d))
iimat = ii; dim(iimat) = c(n,d^2)
yy = iimat \%*\% as.vector(ftrue) + rnorm(n, sd=.3)

cv.obj <- wnet.perm(yy, xfuncs = ii, min.scale = 4, nfeatures = 200, alpha = 1, nperm = 10)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
