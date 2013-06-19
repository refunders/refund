\name{plot.peer}
\alias{plot.peer}
\title{Plotting of estimated regression functions obtained through \code{peer()}}
\description{
Plots the estimate of components of estimated regression function obtained from an \code{peer} object along with the point-wise confidence band.
}
\usage{
\method{plot}{peer}(x, conf=0.95, ylab='Estimated regression function', main=expression(gamma), ...)
}
\arguments{
\item{x}{An peer object}
\item{conf}{Confidene level of point-wise confidence band.}
\item{ylab}{Label on Y-axis.}
\item{main}{Title for the plot.}
\item{...}{Additional arguments passed to the \code{plot} function.}
}

\details{
It plots the components of estimated regression function over sampling points (X-axis). Confidence band is displayed only if the peer object is produced with \code{se=T} in \code{peer()} and is based on pointwise confidence interval without any multiplicity correction. 
}

\author{Madan Gopal Kundu \email{mgkundu@iupui.edu}}

\references{
Kundu, M. G., Harezlak, J., and Randolph, T. W. (2012). Longitudinal functional models with structured penalties. (please contact J. Harezlak at harezlak@iupui.edu)

Randolph, T. W., Harezlak, J, and Feng, Z. (2012). Structured penalties for functional linear models - partially empirical eigenvectors for regression. Electronic Journal of Statistic, 6:323--353.
}

\seealso{
\code{peer}, \code{lpeer}, \code{plot.lpeer}
}

\examples{
data(DTI)
cca = DTI$cca[which(DTI$case == 1),]
DTI = DTI[which(DTI$case == 1),]
fit.cca.peer1 = peer(Y=DTI$pasat, funcs = cca)
plot(fit.cca.peer1)
}
