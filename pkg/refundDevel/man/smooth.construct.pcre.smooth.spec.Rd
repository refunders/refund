\name{smooth.construct.pcre.smooth.spec}
\alias{smooth.construct.pcre.smooth.spec}
\title{mgcv-style constructor for PC-basis functional random effects}
\usage{
  smooth.construct.pcre.smooth.spec(object, data, knots)
}
\arguments{
  \item{object}{a smooth specification object, see
  \code{\link[mgcv]{smooth.construct}}}

  \item{data}{}

  \item{knots}{}
}
\value{
  An object of class \code{"random.effect"}. See
  \code{\link[mgcv]{smooth.construct}} for the elements
  that this object will contain.
}
\description{
  Sets up design matrix for functional random effects based
  on the PC loadings of the covariance operator of the
  random effect process. See
  \code{\link[mgcv]{smooth.construct.re.smooth.spec}} for
  more details.
}

