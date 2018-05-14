\name{getcoef.hd}
\alias{getcoef.hd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Get the estimated intercept and coefficient functions  for sparse linear FOF models
%%  ~~function to do ... ~~
}
\description{Calculate the estimates for \eqn{\mu(t), \beta_i(s,t)}{\mu(t), \beta_i(s,t)} based on the object obtained from \code{\link{cv.hd}}.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getcoef.hd(fit.cv)
%\method{coef}{hd}(fit.cv)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit.cv}{the object obtained from \code{\link{cv.hd}}.
%%     ~~Describe \code{fit.cv} here~~
}
}
\value{a list containing
 \item{mu }{the vector of estimated values for \eqn{\mu(t)} at \code{t.y}, the vector of observation time points for the response function.}
\item{beta }{ a list of length \eqn{p}, the number of functional predictors. Its \eqn{i}-th element is a matrix of the estimated values of \eqn{\beta_i(s,t)} at the full grid of observaion time points created by arguments \code{t.x.list[[i]]} and \code{t.y}  for \code{\link{cv.hd}}. }
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{Ruiyan Luo and Xin Qi
%%  ~~who you are~~
}

\seealso{\code{\link{cv.hd}}}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{ #See the examples in cv.hd().}
