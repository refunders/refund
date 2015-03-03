\name{lofocv}
\Rdversion{1.1}
\alias{lofocv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Leave-one-function-out cross-validation}
\description{
This internal function, called by \code{fosr()} when \code{method="OLS"},
performs efficient leave-one-function-out cross-validation using Demmler-Reinsch orthogonalization to choose the
smoothing parameter.}
\usage{
lofocv(Y, X, S1, argvals=NULL, lamvec = NULL, constr = NULL, maxlam = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
matrix of responses, e.g. with columns corresponding to basis function coefficients.}
  \item{X}{
model matrix.}
  \item{S1}{
penalty matrix.}
  \item{argvals}{The argument values the functions are evaluated}
  \item{lamvec}{
vector of candidate smoothing parameter values.  If \code{NULL}, smoothing parameter is chosen by \code{\link{optimize}}.}
  \item{constr}{
matrix of linear constraints.}
  \item{maxlam}{
maximum smoothing parameter value to consider (when \code{lamvec=NULL}). 
}
}

\value{
if \code{lamvec=NULL}, a list (returned by \code{optimize}) with elements \code{minimum} and \code{objective} giving, respectively, the chosen smoothing parameter and the associated cross-validation score.  Otherwise a 2-column table with the candidate smoothing parameters in the first column and the corresponding cross-validation scores in the second.
}

\author{
Philip Reiss \email{phil.reiss@nyumc.org}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{fosr}}, \code{\link{pwcv}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
