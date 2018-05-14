\name{pspline}
\alias{pspline}
\title{
Univariate P-spline smoothing
}
\description{
Univariate P-spline smoothing with the smoothing parameter selected by leave-one-subject-out cross validation.
}
\usage{
pspline(data, argvals.new = NULL, knots = 35,
            p = 3, m = 2, lambda = NULL,
            search.length = 100,
            lower = -20, upper = 20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
a data frame with three arguments: 
(1) \code{argvals}: observation times;
(2) \code{subj}: subject indices;
(3) \code{y}: values of observations.
Missing values not allowed.
}
  \item{argvals.new}{
a vector of observations times for prediction; if NULL, then the same as \code{data$argvals}.
}
  \item{knots}{
a vector of interior knots or the number of knots for B-spline basis functions to be used; defaults to 35.
}


  \item{p}{
the degrees of B-splines; defaults to 3.
}
  \item{m}{
 the order of differencing penalty; defaults to 2.
}
  \item{lambda}{
the value of the smoothing parameter; defaults to NULL.
}
  \item{search.length}{
the number of equidistant (log scale) smoothing parameters to search; defaults to 100.
}
 \item{lower, upper}{
bounds for log smoothing parameter; defaults are -20 and 20.
}
}
\details{
The function is an implementation of the P-spline smoothing in Eilers and Marx (1996). P-splines uses B-splines as basis functions and employs a differencing penalty on the coefficients. Leave-one-subject-out cross validation is used for selecting the smoothing parameter and a fast algorithm is implemented.

}
\value{
\item{fitted.values}{Fitted mean values}
\item{B}{B-spline design matrix}
\item{theta}{Estimated coefficients}
\item{s}{Eigenvalues}
\item{knots}{Knots}
\item{p}{The degrees of B-splines}
\item{m}{The order of differencing penalty}
\item{lambda}{The value of the smoothing parameter}
\item{argvals.new}{A vector of observations times}
\item{mu.new}{Fitted mean values at \code{argvals.new}}
}
\references{
Paul Eilers and Brian Marx, Flexible smoothing with B-splines and penalties,
Statist. Sci., 11, 89-121, 1996.

Luo Xiao, Cai Li, William Checkley and Ciprian Crainiceanu, Fast covariance estimation
for sparse functional data, Stat. Comput., \doi{10.1007/s11222-017-9744-8}.
}
\author{
Luo Xiao <lxiao5@ncsu.edu>
}

\seealso{
 \code{\link{gam}} in \code{\link{mgcv}}.
}
\examples{
## cd4 data
require(refund)
data(cd4)
n <- nrow(cd4)
T <- ncol(cd4)

id <- rep(1:n,each=T)
t <- rep(-18:42,times=n)
y <- as.vector(t(cd4))
sel <- which(is.na(y))

## organize data
data <- data.frame(y=log(y[-sel]),
argvals = t[-sel],
subj = id[-sel])
data <- data[data$y>4.5,]

## smooth
fit <- pspline(data)

## plot
plot(data$argvals,fit$mu.new,type="p")

## prediction
pred <- predict(fit,quantile(data$argvals,c(0.2,0.6)))
pred
}

\keyword{~Pspline}
