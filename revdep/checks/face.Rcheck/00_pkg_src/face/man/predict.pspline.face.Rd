\name{predict.pspline.face}
\alias{predict.pspline.face}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mean prediction from a P-spline smoothing fit
}
\description{
Predict mean values based on a fit from "pspline".
}
\usage{
\method{predict}{pspline.face}(object, argvals.new,\dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
a fitted object from the R function  "pspline".
}
\item{argvals.new}{
a vector of new time points.
}
\item{...}{further arguments passed to or from other methods.}
}

\value{
Predicted means at \code{argvals.new}.
}

\references{
Luo Xiao, Cai Li, William Checkley and Ciprian Crainiceanu, Fast covariance estimation
for sparse functional data, Stat. Comput., \doi{10.1007/s11222-017-9744-8}.
}
\author{
Luo Xiao <lxiao5@ncsu.edu>
}


\examples{
#See the examples for "pspline".
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{~Pspline}
\keyword{~prediction}% __ONLY ONE__ keyword per line
