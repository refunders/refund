\name{ocean}
\alias{ocean}
\docType{data}
\title{Hawaii ocean data}
\description{Five variables, Salinity, Potential Density, Temperature,  Oxygen  and Chloropigment,  were measured every two meters between 0 and 200 meters below the sea surface at a site located 100 kilometers north of Oahu, Hawaii in 116 different days.}
\usage{data("ocean")}
\format{A list of 5 matrices:
 \describe{
  \item{Salinity}{a matrix with 116 rows and 101 columns.}
  \item{Potential.density}{a matrix with 116 rows and 101 columns.}
  \item{Temperature}{a matrix with 116 rows and 101 columns.}
  \item{Oxygen}{a matrix with 116 rows and 101 columns.}
  \item{Chloropigment}{a matrix with 116 rows and 101 columns.}
  }
}
\source{ \href{http://hahana.soest.hawaii.edu/hot/hot-dogs/cextraction.html}{Hawaii Ocean Time-series}}
\seealso{ \code{\link{cv.ff.interaction}}}
\examples{
 data(ocean)
 str(ocean)
}

