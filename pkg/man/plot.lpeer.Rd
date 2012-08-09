\name{plot.lpeer}
\alias{plot.lpeer}
\title{Plotting of estimated regression functions obtained through \code{\link{lpeer}}()}
\description{
Plots the estimate of components of estimated regression function obtained from an \code{\link{lpeer}} object along with the point-wise confidence band.
}
\usage{
\method{plot}{lpeer}(x, conf=0.95, ...)
}
\arguments{
\item{x}{an object of class \code{"\link{lpeer}"}.}
\item{conf}{confidence level of pointwise confidence band.}
\item{...}{additional arguments passed to the \code{\link{plot}} function.}
}

\details{
Plots the components the of estimated regression function over sampling points (x-axis). Confidence interval produced is pointwise and no multiplicity correction is made. 
}

\author{ Madan Gopal Kundu \email{mgkundu@iupui.edu}}

\references{

Kundu, M. G., Harezlak, J., and Randolph, T. W. (2012). Longitudinal functional models with structured penalties. (please contact J. Harezlak at harezlak@iupui.edu)

Randolph, T. W., Harezlak, J, and Feng, Z. (2012). Structured penalties for functional linear models - partially empirical eigenvectors for regression. \emph{Electronic Journal of Statistics}, 6, 323--353.
}

\seealso{
\code{\link{peer}}, \code{\link{lpeer}}, \code{\link{plot.peer}}
}

\examples{
data(DTI)
cca = DTI$cca[which(DTI$case == 1),]
DTI = DTI[which(DTI$case == 1),]
fit.cca.lpeer1 = lpeer(Y=DTI$pasat, t=DTI$visit, subj=DTI$ID, funcs = cca)
plot(fit.cca.lpeer1)
}
