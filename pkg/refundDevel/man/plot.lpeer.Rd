\name{plot.lpeer}
\alias{plot.lpeer}
\title{Plotting of estimated regression functions obtained through \code{lpeer()}}
\description{
Plots the estimate of components of estimated regression function obtained from an \code{lpeer} object along with the point-wise confidence band.
}
\usage{
\method{plot}{lpeer}(x, conf=0.95, ...)
}
\arguments{
\item{x}{An lpeer object}
\item{conf}{Confidene level of point-wise confidence band.}
\item{...}{Additional arguments passed to the \code{plot} function.}
}

\details{
It plots the components of estimated regression functions over sampling points (X-axis). Confidence bands are displayed only if the lpeer object is produced with \code{se=T} in \code{lpeer()} and are based on pointwise confidence interval without any multiplicity correction.
}

\author{ Madan Gopal Kundu \email{mgkundu@iupui.edu}}

\references{

Kundu, M. G., Harezlak, J., and Randolph, T. W. (2012). Longitudinal functional models with structured penalties. (please contact J. Harezlak at harezlak@iupui.edu)

Randolph, T. W., Harezlak, J, and Feng, Z. (2012). Structured penalties for functional linear models - partially empirical eigenvectors for regression. Electronic Journal of Statistic, 6:323--353.
}

\seealso{
\code{peer}, \code{lpeer}, \code{plot.peer}
}

\examples{
data(DTI)
cca = DTI$cca[which(DTI$case == 1),]
DTI = DTI[which(DTI$case == 1),]
fit.cca.lpeer1 = lpeer(Y=DTI$pasat, t=DTI$visit, subj=DTI$ID, funcs = cca)
plot(fit.cca.lpeer1)
}
