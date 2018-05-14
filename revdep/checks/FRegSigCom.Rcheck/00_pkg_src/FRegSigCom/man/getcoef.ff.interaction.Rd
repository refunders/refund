\name{getcoef.ff.interaction}
\alias{getcoef.ff.interaction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Get the estimated coefficient functions for function-on-function interaction model
}
\description{Calculate the estimates for \eqn{\mu(t), \beta_i(s,t), \gamma_{ij}(u,v,t)}{\mu(t), \beta_i(s,t), \gamma_{ij}(u,v,t)} for function-on-function interaction model (see the model in \code{\link{cv.ff.interaction}}) based on the output object of the function \code{\link{cv.ff.interaction}}, or the function \code{\link{step.ff.interaction}}.
}
\usage{
getcoef.ff.interaction(fit.obj)
%\method{coef}{hd}(fit.cv)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit.obj}{the output object of the function \code{\link{cv.ff.interaction}}, or \code{\link{step.ff.interaction}}}.
}
\value{a list providing the given or selected main effects and interactions, together with the corresponding estimated coefficient functions.
 \item{intercept }{the vector of estimated \eqn{\mu(t)} evaluated at the vector \code{t.y} of the observation points for the response function \eqn{y(t)}.}
\item{main_effects}{the index vector of the input main_effects for \code{\link{cv.ff.interaction}} or the selected main effects by \code{\link{step.ff.interaction}}.}
\item{coef_main }{a list of matrices of the estimated values of the coefficient functions of main effect specified by \code{main_effects}. Each matrix gives the estimated values of \eqn{\beta_i(s,t)} at the two-dimensional grid created by the observation point vectors \code{t.x.list[[i]]} and \code{t.y}, where \eqn{i} is an index in  \code{main_effects}.}
\item{inter_effects}{a matrix of two columns showing the input interactions for \code{\link{cv.ff.interaction}} or the selected interactions by \code{\link{step.ff.interaction}}. Each row shows the indices of the pair of functional variables in an interaction or quadratic effect.}
\item{coef_inter}{a list of three-dimensional arrays  of estimated values of the coefficient functions of interaction or quadratic effects specified by \code{inter_effects}. Each array gives the estimated values of \eqn{\gamma_{ij}(u,v,t)} at the three-dimensional grid created by the observation point vectors \code{t.x.list[[i]]}, \code{t.x.list[[j]]}  and \code{t.y}, where the pair \eqn{i, j} is in  \code{inter_effects}.}
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{Xin Qi and Ruiyan Luo
%%  ~~who you are~~
}

\seealso{\link{cv.ff.interaction}, \link{step.ff.interaction}.}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{ #See the examples in cv.ff.interaction() and step.ff.interaction().}
