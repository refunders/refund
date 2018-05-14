\name{pred.sigcom}
\alias{pred.sigcom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Prediction for linear function-on-function regression using signal compression
%%  ~~function to do ... ~~
}
\description{Make predition for functional response from the CV object obtained by \code{\link{cv.sigcom}}.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pred.sigcom(fit.cv, X.test, Z.test = NULL, t.y.test=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit.cv }{ the CV object obtained by \code{\link{cv.sigcom}}.
%%     ~~Describe \code{fit.cv} here~~
}
  \item{X.test }{new observations for the functional predictors. It is a list of length \eqn{p},  the number of functional predcitors. Each element is the observed  matrix from a functional predictor, with rows repsenting observation vectors and columns corresponding to the observation time points.
%%     ~~Describe \code{X.test} here~~
}
  \item{Z.test }{new observations for the scalar predictors. It is a matrix with rows representing observation vectors and columns respresenting scalar variables. Default is NULL, indicating no scalar predictors.
%%     ~~Describe \code{Z.test} here~~
}
  \item{t.y.test}{a vector of observation time points where values of predicted response curves are to be calculated. If \code{t.y.test}=NULL (default), \code{t.y} in \code{\link{cv.sigcom}} will be used.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{A matrix containing the predicted response for the new
  observations. The number of rows is equal to  the sample size of the
  new data set, and the number of columns is equal to the length of
  \code{t.y.test} or \code{t.y} when \code{t.y.test}=NULL. 
}
\references{Ruiyan Luo and Xin Qi, (2017) Function-on-Function Linear Regression by Signal Compression, Journal of the American Statistical Association.  112(518), 690-705. \url{http://www.tandfonline.com/doi/abs/10.1080/01621459.2016.1164053}
%% ~put references to the literature/web site here ~
}
\author{Ruiyan Luo and Xin Qi
%%  ~~who you are~~
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{cv.sigcom}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{#See the examples in cv.sigcom().
}
