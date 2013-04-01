\name{predict.wnet}
\alias{predict.wnet}
\title{Prediction for generalized elastic net in the wavelet domain}
\usage{
  \method{predict}{wnet}(object, newx, covt = NULL, ...)
}
\arguments{
  \item{object}{a fitted \code{wnet}-object}

  \item{newx}{Matrix of new values for x at which predictions are to be made.}

  \item{covt}{Matrix of new covariates.}
  
  \item{...}{used}
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

