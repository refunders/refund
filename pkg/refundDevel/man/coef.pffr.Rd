\name{coef.pffr}
\alias{coef.pffr}
\title{Get estimated coefficients from a pffr}
\usage{
\method{coef}{pffr} (object, raw = FALSE, se = TRUE, freq = FALSE, n1 = 100, n2 = 40, n3 = 20, ...)
}
\arguments{
  \item{object}{a fitted \code{pffr}-object}

  \item{raw}{logical, defaults to FALSE. If TRUE, the
  function simply returns \code{object$coefficients}}

  \item{se}{logical, defaults to TRUE. Return estimated
  standard error of the estimates?}

  \item{freq}{logical, defaults to FALSE. If FALSE, use
  posterior variance \code{object$Vp} for variability
  estimates, else use \code{object$Ve}. See
  \code{\link[mgcv]{gamObject}}}

  \item{n1}{default: 100 see below}

  \item{n2}{default: 40 see below}

  \item{n3}{default: 20. These parameters give the number
  of gridpoints for 1-/2-/3-dimensional smooth terms for
  the marginal equidistant grids over the range of the
  covariates at which the estimated effects are evaluated.}

  \item{...}{other arguments, not used.}
}
\value{
  If \code{raw==FALSE}, a list containing \itemize{ \item
  \code{pterms} a matrix containing the parametric /
  non-functional coefficients (and, optionally, their se's)
  \item \code{smterms} a named list with one entry for each
  smooth term in the model. Each entry contains \itemize{
  \item \code{coef} a matrix giving the grid values over
  the covariates, the estimated effect (and, optionally,
  the se's). The first covariate varies the fastest. \item
  \code{x, y, z} the unique gridpoints used to evaluate the
  smooth/coefficient function/coefficient surface \item
  \code{xlim, ylim, zlim} the extent of the x/y/z-axes
  \item \code{xlab, ylab, zlab} the names of the covariates
  for the x/y/z-axes \item \code{dim} the dimensionality of
  the effect \item \code{main} the label of the smooth term
  (a short label, same as the one used in
  \code{summary.pffr}) }}
}
\description{
  Returns estimated coefficient functions/surfaces
  \eqn{\beta(t), \beta(s,t)} and estimated smooth effects
  \eqn{f(z), f(x,z)} or \eqn{f(x, z, t)} and their
  point-wise estimated standard errors. Not implemented for
  smooths in more than 3 dimensions.
}
\author{
  Fabian Scheipl
}

