\name{cor.face}
\alias{cor.face}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extraction of correlation and mean from a \code{face.sparse} object
}
\description{
Extraction of correlation and mean from a \code{face.sparse} object
}
\usage{
cor.face(object,argvals.new,option="raw")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ A \code{face.sparse} object.
}
\item{argvals.new}{Where to evalulate correlation and mean.
}
\item{option}{Defaults to "raw"; if "smooth", then extract correlation from smoothed covariance function.}
}

\value{
\item{argvals.new}{Where to evaluate correlation and mean.}
\item{option}{Defaults to "raw"; if "smooth", then extract correlation from smoothed covariance function.}
\item{Cor}{estimated correlation matrix at \code{argvals.new}}
\item{mu}{estimated group/population mean at \code{argvals.new}}
}

\references{
Luo Xiao, Cai Li, William Checkley and Ciprian Crainiceanu, Fast covariance estimation
for sparse functional data, Stat. Comput., \doi{10.1007/s11222-017-9744-8}.

}
\author{
Luo Xiao <lxiao5@ncsu.edu>
}


\examples{
# See the examples for "face.sparse".
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{~face.sparse}
