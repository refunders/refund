\name{pred.ff.interaction}
\alias{pred.ff.interaction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Prediction for a linear FOF regression
  model with two-way interactions
}
\description{Make predition for new observations of functional
  predictors based on the fitted function-on-function interaction model
  in the output of \code{\link{cv.ff.interaction}} with given main and
  two-way interaction effects, or the model in the output of \code{\link{step.ff.interaction}} with selected main and two-way interaction effects.}
\usage{
pred.ff.interaction(fit.obj, X.test) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit.obj}{the output object of the function \code{\link{cv.ff.interaction}}, or \code{\link{step.ff.interaction}}.}
%%     ~~Describe \code{fit.cv} here~~
 \item{X.test}{new observations of functional predictors. It is a list of the same length and the same structure as the input \code{X} for \code{\link{cv.ff.interaction}} or \code{\link{step.ff.interaction}}. }
 }
\value{A matrix containing the predicted values of response curves
  evaluated at \code{t.y} for the new observations. The number of rows
  of the matrix  equals to the sample size of the new data set, and the
  number of columns equals to the length of \code{t.y}, the input in
  \code{\link{cv.ff.interaction}} or \code{\link{step.ff.interaction}}.}

\author{Xin Qi and Ruiyan Luo}

\seealso{\code{\link{cv.ff.interaction}}, \code{\link{step.ff.interaction}}}

\examples{ #See the examples in cv.ff.interaction() and step.ff.interaction().}
