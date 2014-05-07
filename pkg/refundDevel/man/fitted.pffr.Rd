\name{fitted.pffr}
\alias{fitted.pffr}
\title{Obtain fitted values for a pffr fit}
\usage{
\method{fitted}{pffr}(object, reformat = TRUE, ...)
}
\arguments{
  \item{object}{a fitted \code{pffr}-object}

  \item{reformat}{logical, defaults to TRUE. Should
  residuals be returned in \code{n x yindex} matrix form
  (default for regular grid data) or in the shape of the
  originally supplied \code{ydata} argument or (default for
  sparse/irregular data), or, if \code{FALSE}, simply in
  the long vector shape returned by \code{resid.gam()}?}

  \item{...}{not used}
}
\value{
A matrix (grid data) or \code{data.frame} (sparse/irregular
data) or vector (if \code{!reformat}) of fitted values
}
\description{
Returns the linear predictor for the model data. See
\code{\link{predict.pffr}} for alternative options.
}
\author{
Fabian Scheipl
}

