\name{plot.peer}
\alias{plot.peer}
\title{Plotting of estimated regression functions obtained through \code{\link{peer}()}}
\description{
Plots the estimate of components of estimated regression function obtained from a \code{\link{peer}} object along with the point-wise confidence band.
}
\usage{
\method{plot}{peer}(x, conf=0.95, ...)
}
\arguments{
\item{x}{an object of class \code{"\link{peer}"}.}
\item{conf}{confidence level of pointwise confidence band.}
\item{...}{additional arguments passed to the \code{\link{plot}} function.}
}

\details{
Plots the estimated regression function over sampling points (x-axis). The confidence interval produced is pointwise and no multiplicity correction is made. 
}

\author{ Madan Gopal Kundu \email{mgkundu@iupui.edu}}

\references{
Kundu, M. G., Harezlak, J., and Randolph, T. W. (2012). Longitudinal functional models with structured penalties. (please contact J. Harezlak at harezlak@iupui.edu)

Randolph, T. W., Harezlak, J, and Feng, Z. (2012). Structured penalties for functional linear models - partially empirical eigenvectors for regression. \emph{Electronic Journal of Statistics}, 6, 323--353.
}

\seealso{
\code{\link{peer}}, \code{\link{lpeer}}, \code{\link{plot.lpeer}}
}

\examples{
data(DTI)
cca = DTI$cca[which(DTI$case == 1),]
DTI = DTI[which(DTI$case == 1),]
fit.cca.peer1 = peer(Y=DTI$pasat, funcs = cca)
plot(fit.cca.peer1)
}
