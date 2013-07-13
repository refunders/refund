\name{predict.wnet}
\alias{predict.wnet}
\title{Prediction method for generalized elastic net in the wavelet domain}
\usage{
  \method{predict}{wnet}(object, newx, newcovt = NULL, ...)
}
\arguments{
  \item{object}{a fitted \code{\link{wnet}} object.}

  \item{newx}{matrix of new values for functional predictors.}

  \item{newcovt}{matrix of new covariate values.}
  
  \item{\dots}{not currently used.}
}
\value{
  A vector of predicted values.
}
\description{
  Produces predicted values for a wnet object, given new functional predictors and scalar covariates.
}
\author{
  Lan Huo
}
\seealso{
  \code{\link{wnet}}
}

