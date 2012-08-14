\name{fitted.pffr}
\alias{fitted.pffr}
\title{Obtain fitted values for a pffr fit}
\usage{
  \method{fitted}{pffr} (object, reformat = TRUE, ...)
}
\arguments{
  \item{object}{a fitted \code{pffr}-object}

  \item{reformat}{logical, defaults to TRUE. Should
  residuals be returned in \code{n x yindex} matrix form
  (default) or in the long vector shape returned by
  \code{fitted.default()}?}

  \item{...}{not used}
}
\value{
  A matrix or vector of fitted values
}
\description{
  Returns the linear predictor for the model data. See
  \code{\link{predict.pffr}} for alternative options.
}
\author{
  Fabian Scheipl
}

