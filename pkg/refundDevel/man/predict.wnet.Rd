\name{predict.wnet}
\alias{predict.wnet}
\title{Prediction method for generalized elastic net in the wavelet domain}
\usage{
  \method{predict}{wnet}(object, newx, newcovt = NULL, ...)
}
\arguments{
  \item{object}{a fitted \code{wnet}-object}

  \item{newx}{matrix of new values for functional predictors.}

  \item{newcovt}{matrix of new covariate values.}
  
  \item{...}{Other argument.}
}
\value{
  The predicted value.
}
\description{
  Prediction function fro wnet object.
}
\author{
  Lan Huo
}
\seealso{
  \code{wnet}
}

