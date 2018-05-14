\name{getcoef.nonlinear}
\alias{getcoef.nonlinear}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Get the estimated intercept and nonlinear functions
%%  ~~function to do ... ~~
}
\description{Calculate the estimates for \eqn{\mu(t), F_i(x,s,t)}{\mu(t), \alpha(t), F_i(x,s,t)} based on the object obtained from \code{\link{cv.nonlinear}}.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getcoef.nonlinear(fit.cv, n.x.grid = 100)
%\method{coef}{nonlinear}(fit.cv, n.x.grid=100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit.cv }{the object obtained from \code{\link{cv.nonlinear}}.
%%     ~~Describe \code{fit.cv} here~~
}
  \item{n.x.grid }{the number of grid points of \eqn{x}. The estimated \eqn{F_i(x,s,t)} is calculated in a three-dimensional grid of \eqn{(x,s,t)}. The grid points of \eqn{s} and \eqn{t} are the observation points of \eqn{X_i(s)} and \eqn{Y(t)} used in \code{\link{cv.nonlinear}}, respectively. The grid of \eqn{x} includes \code{n.x.grid} equally spaced values between the minimum and maximum of all the discretely observed values of \eqn{X_i(s)}. Default of  n.x.grid is 100.
%%     ~~Describe \code{n.x.grid} here~~
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{a list containing
 \item{mu }{the vector of estimated values of \eqn{\mu(t)} at the observation points of the response function.}
\item{F }{ a list of length \eqn{p},  the number of functional predictors. Its \eqn{i}-th element is a three dimensional array with estimated values of \eqn{F_i(x,s,t)} on the three-dimensional grid \code{X.grid[[i]]*t.x.list[[i]]*t.y} (see below). }
\item{X.grid }{a list of length \eqn{p}. Its \eqn{i}-th element is the vector of grid points for \eqn{x} and includes  \code{n.x.grid} equally spaced values between the minimum and maximum of all the discretely observed values of \eqn{X_i(s)}.  }
\item{t.x.list}{one of the arguments in \code{\link{cv.nonlinear}}, specifying the list of the vectors of obesrvation points for \eqn{X_i(s)}{X_i(s)}, \eqn{1\le i\le p}{1\le i\le p}. }
\item{t.y}{one of the arguments in \code{\link{cv.nonlinear}}, specifying the vector of obesrvation points of the response curve \eqn{Y(t)}{Y(t)}.}

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{Ruiyan Luo and Xin Qi,
%% ~put references to the literature/web site here ~
}
\author{Ruiyan Luo and Xin Qi
%%  ~~who you are~~
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{cv.nonlinear}}.
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{#See the examples in cv.nonlinear().
}
