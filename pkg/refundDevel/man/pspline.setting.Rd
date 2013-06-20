\name{pspline.setting}
\alias{pspline.setting}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Settings for Univariate P-splines
}
\description{
Creates the settings for univariate \emph{P}-splines (Eilers and Marx, 1996).
}
\usage{
pspline.setting(x, knots = 35, p = 3, m = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  vector of data points
}
  \item{knots}{
  vector of knots or  number of equidistant knots; defaults to 35
}
  \item{p}{
 degrees of B-splines; defaults to 3
}
  \item{m}{
 order of differencing penalty; defaults to 2
}
}
\details{
See Xiao et al. (2013) for details on the list of return values.
}
\value{
\item{A}{matrix with orthonormal columns}
\item{s}{vector of eigenvalues}
\item{Sigi.sqrt}{squared root inverse of t(B)B}
\item{U}{matrix of eigenvectors}
}
\references{
Eilers, P.H.C. and Marx, B.D. (1996),
Flexbile smoothing with B-splines and penalties (with Discussion),
\emph{Statistical Science}, 11(2), 89-121.

Xiao, L., Li, Y., and Ruppert, D. (2013),
Fast bivariate \emph{P}-splines: the sandwich smoother,
\emph{Journal of the Royal Statistical Society: Series B}, 75(3), 577-599.
}
\author{
Luo Xiao \email{lxiao@jhsph.edu}
}

\examples{

## Univariate P-spline smoothing
n <- 20
x <- (1:n)/n
y <- sin(2*pi*x) + 0.1*rnorm(20)

## Creates the setting for P-splines
List <- pspline.setting(x,knots=10)
A <- List$A
s <- List$s
Sigi.sqrt <- List$Sigi.sqrt
U <- List$U

## Fix the value of the smoothing parameter
lambda <- 2
## estimates of the coefficient vector
Theta <- Sigi.sqrt\%*\%U\%*\% diag(1/(1+lambda*s))\%*\%t(A)\%*\%y

## Fitted values
yhat <- A\%*\%diag(1/(1+lambda*s))\%*\%t(A)\%*\%y

## Draw the fitted values
plot(x,yhat,type="l",main="P-spline smoothing",ylab="fitted values")
points(x,y)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{~\emph{P}-splines }
%\keyword{~smoothing}
%\keyword{~regression}
