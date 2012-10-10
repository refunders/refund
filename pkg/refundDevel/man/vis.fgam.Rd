\name{vis.fgam}
\Rdversion{1.0}
\alias{vis.fgam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visualization of FGAM objects
}
\description{
Produces perspective or contour plot views of an estimated surface corresponding to \code{af} terms fit using \code{\link{fgam}} or plots "slices"" of the estimated surface or estimated second derivative surface with one of its arguments fixed and corresponding twice-standard error "Bayesian" confidence bands constructed using the method in Marra and Wood (2012).  See the details.
}
\usage{
vis.fgam(object, af.term, xval = NULL, tval = NULL, deriv2 = FALSE, theta = 50, plot.type = "persp", ticktype = "detailed", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an \code{fgam} object, produced by \code{\link{gam}}().
}
  \item{af.term}{
the name of the functional predictor to be plotted.  Only important if multiple \code{af} terms are fit.  Defaults to the first \code{af} term in \code{object$call}.
}
  \item{xval}{
a number in the range of functional predictor to be plotted.  The surface will be plotted with the first argument of the estimated surface fixed at this value.
}
  \item{tval}{
a number in the domain of the functional predictor to be plotted.  The surface will be plotted with the second argument of the estimated surface fixed at this value. Ignored if \code{xval} is specified.
}
  \item{deriv2}{
If \code{TRUE}, plot the estimated second derivative surface along with Bayesian confidence bands.  Only implemented for the "slices" plot from either \code{xval} or \code{tval} being specified.
}
  \item{theta}{
%%     ~~Describe \code{theta} here~~
}
  \item{plot.type}{
one of "\code{contour}" or "\code{persp}".  Ignored if either \code{xval} or \code{tval} is specified.
}
  \item{ticktype}{
how to draw the tick marks if \code{plot.type="persp"}.  Defaults to \code{"detailed"}.
}
  \item{\dots}{
other options to be passed to \code{\link{persp}}, \code{\link{levelplot}}, or \code{\link{plot}}.
}
}
\details{
The confidence bands used when plotting slices of the estimated surface or second derivative surface are the ones proposed in Marra and Wood (2012).  These are a generalization of the "Bayesian" intevals of Wahba (1983) with an adjustment for the uncertainty about the model intercept. The estimated covariance matrix of the model parameters is obtained from assuming a particular Bayesian model on the parameters.
}
\value{
Simply produces a plot.
}
\references{
McLean M.W. et al. (2012) \emph{Journal of Computational and Graphical Statistics}, to appear.

Marra, G and S.N. Wood (2012) Coverage Properties of Confidence Intervals for Generalized Additive Model Components. \emph{Scandinavian Journal of Statistics}, \bold{39}(1), 53--74.

Wabha, G (1983) \emph{Journal of the Royal Statistical Society, Series B}, \bold{45}, 133--150. }
\author{
Mathew W. McLean \email{mwm79@cornell.edu}
}

\seealso{
\code{\link{vis.gam}}, \code{\link{plot.gam}}, \code{\link{fgam}}, \code{\link{persp}}, \code{\link{levelplot}}
}
\examples{
################# DTI Example #####################
data(DTI)

## only consider first visit and cases ( since no PASAT scores for controls)
y <- DTI$pasat[DTI$visit==1 & DTI$case==1]
X <- DTI$cca[DTI$visit==1 & DTI$case==1,]

## remove samples containing missing data
ind <- rowSums(is.na(X))>0

y <- y[!ind]
X <- X[!ind,]

## fit the fgam using measurements along cca tract as functional predictor with PASAT as response
## using 8 cubic B-splines for each marginal bases with third order marginal difference penalties
## specifying gamma>1 enforces more smoothing when using GCV to choose smoothing parameters
fit <- fgam(y~af(X,splinepars=list(k=c(8,8),m=list(c(2,3),c(2,3)))),gamma=1.2)

## contour plot of the fitted surface
vis.fgam(fit,plot.type='contour')

## similar to Figure 5 from McLean et al.
## Bands seem too conservative in some cases
xval=runif(1,min(fit$fgam$ft[[1]]$Xrange),max(fit$fgam$ft[[1]]$Xrange))
tval=runif(1,min(fit$fgam$ft[[1]]$xind),max(fit$fgam$ft[[1]]$xind))
par(mfrow=c(4,1))
vis.fgam(fit,af.term='X',deriv2=FALSE,xval=xval)
vis.fgam(fit,af.term='X',deriv2=FALSE,tval=tval)
vis.fgam(fit,af.term='X',deriv2=TRUE,xval=xval)
vis.fgam(fit,af.term='X',deriv2=TRUE,tval=tval)

}
\keyword{hplot}
