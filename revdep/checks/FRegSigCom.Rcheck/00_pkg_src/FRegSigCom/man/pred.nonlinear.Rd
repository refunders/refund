\name{pred.nonlinear}
\alias{pred.nonlinear}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Prediction for nonlinear function-on-function regression
%%  ~~function to do ... ~~
}
\description{Make predition for response based on new observations of predictor curves from the CV object obtained by \code{\link{cv.nonlinear}}.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pred.nonlinear(fit.cv, X.test, t.y.test=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit.cv}{ the CV object obtained by \code{\link{cv.nonlinear}}.
%%     ~~Describe \code{fit.cv} here~~
}
 \item{X.test}{new observations of functional predictors. It is a list of length \eqn{p}, the number of functional predcitors. Each element is the observed  matrix from a functional predictor, with rows repsenting observation vectors and columns corresponding to the observation time points.
%%     ~~Describe \code{X.test} here~~
}
\item{t.y.test}{a vector of observation time points where values of predicted response curves are to be calculated. If \code{t.y.test}=NULL (default), \code{t.y} in \code{\link{cv.nonlinear}} will be used.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{A matrix which contains the predicted values of response
  curves. The number of rows equal to  the sample size of the new data
  set, and the number of columns is equal to the length of
  \code{t.y.test} or \code{t.y} when \code{t.y.test}=NULL. 
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{Ruiyan Luo and Xin Qi
%%  ~~who you are~~
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{cv.nonlinear}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{#See the examples in cv.nonlinear().}
