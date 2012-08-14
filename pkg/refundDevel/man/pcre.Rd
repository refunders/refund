\name{pcre}
\alias{pcre}
\title{pffr-constructor for PC-basis functional random effects}
\usage{
pcre(id, efunctions, evalues, yind)
}
\arguments{
  \item{id}{grouping variable: a factor.}

  \item{efunctions}{matrix of eigenfunction evaluations on
  grid points \code{yind} (<length of \code{yind}> x <no. of
  eigenfunctions used>).}

  \item{evalues}{eigenvalues associated with
  \code{efunctions.}}

  \item{yind}{vector of grid points on which responses
  \eqn{Y(t)} are evaluated.}
}
\value{
  a list used internally for constructing an appropriate
  call to \code{mgcv::gam}
}
\description{Fits functional random effects, by restricting them to the span of the leading eigenfunctions of the functional covariance}
\details{\code{pcre}
  fits functional random effects \eqn{B_i(t)} for a
  grouping variable \code{id}, using as a basis the
  eigenfunctions \eqn{\phi_m(t)} in \code{efunctions} with
  eigenvalues \eqn{\lambda_m} in \code{evalues}:
  \eqn{B_i(t) \approx \sum_m^M \phi_m(t)\delta_{im}} with
  independent \eqn{\delta_{im} \sim N(0,
  \sigma^2\lambda_m)}, where \eqn{\sigma^2} is estimated
  and controls the overall contribution of the \eqn{B_i(t)},
  while the relative importance of the \eqn{M}
  eigenfunctions is controlled by the supplied eigenvalues
  \code{evalues}. Can be used to model smooth residual functions if
  \code{id} is simply an index of functional observations. This is an
  experimental feature and not well tested yet, use at your
  own risk.
}

