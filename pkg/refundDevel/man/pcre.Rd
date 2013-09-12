\name{pcre}
\alias{pcre}
\title{pffr-constructor for functional principal component-based functional random effects.}
\usage{
  pcre(id, efunctions, evalues, yind, ...)
}
\arguments{
  \item{id}{grouping variable a factor}

  \item{efunctions}{matrix of eigenfunction evaluations on
  gridpoints \code{yind} (<length of \code{yind}> x <no. of
  used eigenfunctions>)}

  \item{evalues}{eigenvalues associated with
  \code{efunctions}}

  \item{yind}{vector of gridpoints on which
  \code{efunctions} are evaluated.}

  \item{...}{not used}
}
\value{
  a list used internally for constructing an appropriate
  call to \code{mgcv::gam}
}
\description{
  pffr-constructor for functional principal component-based
  functional random effects.
}
\section{Details}{
  Fits functional random effects \eqn{B_i(t)} for a
  grouping variable \code{id} using as a basis the
  functions \eqn{\phi_m(t)} in \code{efunctions} with
  variances \eqn{\lambda_m} in \code{evalues}: \eqn{B_i(t)
  \approx \sum_m^M \phi_m(t)\delta_{im}} with independent
  \eqn{\delta_{im} \sim N(0, \sigma^2\lambda_m)}, where
  \eqn{\sigma^2} is (usually) estimated and controls the
  overall contribution of the \eqn{B_i(t)} while the
  relative importance of the \eqn{M} basisfunctions is
  controlled by the supplied variances \code{lambda_m}. Can
  be used to model smooth residuals if \code{id} is simply
  an index of observations.

  \code{efunctions} and \code{evalues} are typically
  eigenfunctions and eigenvalues of an estimated covariance
  operator for the functional process to be modeled, i.e.,
  they are a functional principal components basis.
}
\examples{
\dontrun{
residualfunction <- function(t){
#generate quintic polynomial error functions
    drop(poly(t, 5)\%*\%rnorm(5, sd=sqrt(2:6)))
}
# generate data Y(t) = mu(t) + E(t) + white noise
set.seed(1122)
n <- 50
T <- 30
t <- seq(0,1, l=T)
# E(t): smooth residual functions
E <- t(replicate(n, residualfunction(t)))
int <- matrix(scale(3*dnorm(t, m=.5, sd=.5) - dbeta(t, 5, 2)), byrow=T, n, T)
Y <- int + E + matrix(.2*rnorm(n*T), n, T)
data <- data.frame(Y=I(Y))
# fit model under independence assumption:
summary(m0 <- pffr(Y ~ 1, yind=t, data=data))
# get first 5 eigenfunctions of residual covariance
# (i.e. first 5 functional PCs of empirical residual process)
Ehat <- resid(m0)
fpcE <- fpca.sc(Ehat, npc=5)
efunctions <- fpcE$efunctions
evalues <- fpcE$evalues
data$id <- factor(1:nrow(data))
# refit model with fpc-based residuals
m1 <- pffr(Y ~ 1 + pcre(id=id, efunctions=efunctions, evalues=evalues, yind=t), yind=t, data=data)
t1 <- predict(m1, type="terms")
summary(m1)
#compare squared errors
mean((int-fitted(m0))^2)
mean((int-t1[[1]])^2)
mean((E-t1[[2]])^2)
# compare fitted & true smooth residuals and fitted intercept functions:
layout(t(matrix(1:4,2,2)))
matplot(t(E), lty=1, type="l", ylim=range(E, t1[[2]]))
matplot(t(t1[[2]]), lty=1, type="l", ylim=range(E, t1[[2]]))
plot(m1, select=1, main="m1", ylim=range(Y))
lines(t, int[1,], col=rgb(1,0,0,.5))
plot(m0, select=1, main="m0", ylim=range(Y))
lines(t, int[1,], col=rgb(1,0,0,.5))
}
}
\author{
  Fabian Scheipl
}

