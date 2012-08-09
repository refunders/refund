\name{fosr.perm}  
\Rdversion{1.1}
\alias{fosr.perm}
\alias{fosr.perm.fit}
\alias{fosr.perm.test}
\alias{plot.fosr.perm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Permutation testing for function-on-scalar regression}
\description{
\code{fosr.perm()} is a wrapper function calling \code{fosr.perm.fit()}, which fits models to permuted data, followed by \code{fosr.perm.test()}, which performs the actual simultaneous hypothesis test.  Calling the latter two functions separately may be useful for performing tests at different significance levels.  By default, \code{fosr.perm()} produces a plot using the plot function for class \code{fosr.perm}.
}
\usage{
fosr.perm(Y=NULL, fdobj=NULL, X, con = NULL, X0 = NULL, con0 = NULL,
          argvals = NULL, lambda = NULL, lambda0 = NULL, multi.sp = FALSE, 
          nperm, level = 0.05, plot = TRUE, xlabel = "", title = NULL, 
          prelim = 15, ...)

fosr.perm.fit(Y = NULL, fdobj = NULL, X, con = NULL, X0 = NULL, con0 = NULL,                     
              argvals = NULL, lambda = NULL, lambda0 = NULL, multi.sp = FALSE, 
              nperm, prelim, ...)

fosr.perm.test(x, level=.05)

\method{plot}{fosr.perm}(x, level = .05, xlabel = "", title = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y,fdobj}{
the functional responses, given as either an \eqn{n\times d} matrix \code{Y} or a functional data object (class \code{"\link[fda]{fd}"}) as in the \pkg{fda} package.}
  \item{X}{
the design matrix, whose columns represent scalar predictors.}
  \item{con}{
a row vector or matrix of linear contrasts of the coefficient functions, to be restricted to equal zero.}
  \item{X0}{
design matrix for the null-hypothesis model.  If \code{NULL}, the null hypothesis is the intercept-only model.}
  \item{con0}{
linear constraints for the null-hypothesis model.}
  \item{argvals}{
the \eqn{d} argument values at which the coefficient functions will be evaluated.}
  \item{lambda}{
smoothing parameter value.  If \code{NULL}, the smoothing parameter(s) will be estimated.  See \code{\link{fosr}} for details.}
  \item{lambda0}{
smoothing parameter for null-hypothesis model.}
  \item{multi.sp}{a logical value indicating whether separate smoothing parameters should be estimated for each coefficient function.  Currently must be \code{FALSE} if \code{method = "OLS"}.}
  \item{nperm}{
number of permutations.}
  \item{prelim}{
number of preliminary permutations.  The smoothing parameter in the main permutations will be fixed to the median value from these preliminary permutations.  If \code{prelim=0}, this is not done.}
\item{level}{significance level for the simultaneous test.}\item{plot}{logical value indicating whether to plot the real- and permuted-data pointwise F-type statistics.}\item{xlabel}{x-axis label for plots.}\item{title}{title for plot.}
  \item{x}{object of class \code{fosr.perm}, outputted by \code{fosr.perm}, \code{fosr.perm.fit}, or \code{fosr.perm.test}.}
  \item{\dots}{
for \code{fosr.perm} and \code{fosr.perm.fit}, additional arguments passed to \code{\link{fosr}}.  These arguments may include \code{max.iter}, \code{method}, \code{gam.method}, and \code{scale}.  For \code{plot.fosr.perm}, graphical parameters (see \code{\link{par}}) for the plot.}
}

\value{
\code{fosr.perm} or \code{fosr.perm.test} produces an object of class \code{fosr.perm}, which is a list with the elements below.  \code{fosr.perm.fit} also outputs an object of this class, but without the last five elements.
\item{F}{pointwise F-type statistics at each of the points given by \code{argvals}.}
\item{F.perm}{a matrix, each of whose rows gives the pointwise F-type statistics for a permuted data set.}
\item{argvals}{points at which F-type statistics are computed.}
\item{lambda.real}{smoothing parameter(s) for the real-data fit.}
\item{lambda.prelim}{smoothing parameter(s) for preliminary permuted-data fits.}
\item{lambda.perm}{smoothing parameter(s) for main permuted-data fits.}
\item{lambda0.real, lambda0.prelim, lambda0.perm}{as above, but for null hypothesis models.}
\item{level}{significance level of the test.}
\item{critval}{critical value for the test.}
\item{signif}{vector of logical values indicating whether significance is attained at each of the points \code{argvals}.}
\item{n2s}{subset of \{1, \dots, \code{length(argvals)}\} identifying the points at which the test statistic changes from non-significant to significant.}
\item{s2n}{points at which the test statistic changes from significant to non-significant.}
}
\references{
Reiss, P. T., Huang, L., and Mennes, M. (2010).  Fast function-on-scalar regression with penalized basis expansions.  \emph{International Journal of Biostatistics}, 6(1), article 28.  Available at \url{http://works.bepress.com/phil_reiss/16/}
}
\author{
Philip Reiss \email{phil.reiss@nyumc.org} and Lan Huo}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{fosr}}
}
\examples{
\dontrun{
# Test effect of region on mean temperature in the Canadian weather data  
# The next two lines are taken from the fRegress.CV help file (package fda)
smallbasis  <- create.fourier.basis(c(0, 365), 25)
tempfd <- smooth.basis(day.5,
          CanadianWeather$dailyAv[,,"Temperature.C"], smallbasis)$fd
  
Xreg = cbind(1, model.matrix(~factor(CanadianWeather$region)-1))
conreg = matrix(c(0,1,1,1,1), 1)   # constrain region effects to sum to 0

# This is for illustration only; for a real test, must increase nperm 
# (and probably prelim as well)
regionperm = fosr.perm(fdobj=tempfd, X=Xreg, con=conreg, method="OLS", nperm=10, prelim=3)

# Redo the plot, using axisIntervals() from the fda package
plot(regionperm, axes=FALSE, xlab="")
box()
axis(2)
axisIntervals(1)
}
}
