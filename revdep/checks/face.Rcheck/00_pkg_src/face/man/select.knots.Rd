\name{select.knots}
\alias{select.knots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Knots selection for P-spline smoothing
}
\description{
Construct knots from either quantiles of observed time points or
equally-spaced time points.
}
\usage{
select.knots(t,knots=10,p=3,option="equally-spaced")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{ Observed time points.
}
\item{knots}{ Number of interior knots.
}
\item{p}{Degrees of B-splines to be used.}
\item{option}{Default "equally-spaced": equally-spaced time points in the range of \code{t}; if "quantile", 
then quantiles of \code{t}.}
}

\details{
The number of knots in the output will be \code{knot} plus  2 times \code{p}; and
the B-spline basis matrix constructed from this vector of knots with degrees \code{p}
will be \code{knots} plus \code{p}.
}
\value{
A vector of knots
}

\author{
Luo Xiao <lxiao5@ncsu.edu>
}


\examples{
t <- rnorm(100)
knots <- select.knots(t)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{~Pspline}