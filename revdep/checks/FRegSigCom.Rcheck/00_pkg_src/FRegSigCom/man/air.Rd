\name{air}
\alias{air}
\docType{data}
\title{Air quality data}
\description{Data collected hourly in 355 days (days with missing values
  removed) in a significantly polluted area within an Italian city.}
\usage{data("air")}
\format{A list of 7 matrices with 355 rows and 24 columns:
 \describe{
  \item{NO2}{Hourly observation of concentration level of NO2 in 355 days}
  \item{CO}{Hourly observation of concentration level of CO in 355 days}
  \item{NMHC}{Hourly observation of concentration level of NMHC in 355 days}
  \item{NOx}{Hourly observation of concentration level of NOx in 355 days}
  \item{C6H6}{Hourly observation of concentration level of C6H6 in 355 days}
  \item{temperature}{Hourly observation of concentration level of temperature in 355 days}
  \item{humidity}{Hourly observation of concentration level of humidity in 355 days}
  }
}
\keyword{datasets}
\references{
  De Vito, S. etal (2008) Sensors and Actuators B: Chemical, 129: 50-757.

 Xin Qi and Ruiyan Luo. (Accepted) Nonlinear function on function additive model with multiple predictor curves. Stistica Sinica.
 }
 \source{ \href{https://archive.ics.uci.edu/ml/datasets/Air+quality}{Data link}}
\seealso{ \code{\link{cv.nonlinear}}}
\examples{
 data(air)
 str(air)
}

