##' Default plotting for functional principal component regression output
##'
##' Inputs an object created by \code{\link{fpcr}}, and plots the estimated
##' coefficient function.
##'
##'
##' @param x an object of class \code{"\link{fpcr}"}.
##' @param se if \code{TRUE} (the default), upper and lower lines are added at
##' 2 standard errors (in the Bayesian sense; see Wood, 2006) above and below
##' the coefficient function estimate.  If a positive number is supplied, the
##' standard error is instead multiplied by this number.
##' @param col color for the line(s).  This should be either a number, or a
##' vector of length 3 for the coefficient function estimate, lower bound, and
##' upper bound, respectively.
##' @param lty line type(s) for the coefficient function estimate, lower bound,
##' and upper bound.
##' @param xlab,ylab x- and y-axis labels.
##' @param \dots other arguments passed to the underlying plotting function.
##' @return None; only a plot is produced.
##' @author Philip Reiss \email{phil.reiss@@nyumc.org}
##' @seealso \code{\link{fpcr}}, which includes an example.
##' @references Wood, S. N. (2006). \emph{Generalized Additive Models: An
##' Introduction with R}. Boca Raton, FL: Chapman & Hall.
##' @export
##' @importFrom graphics matplot
plot.fpcr = function(x, se=TRUE, col=1, lty=c(1,2,2), xlab="", ylab="Coefficient function", ...) {
    if (se) {
        if (is.numeric(se)) se.mult <- se
        else se.mult <- 2
        se.mult = max(se.mult, 0)
    }
    if (se) matplot(x$argvals, cbind(x$fhat, x$fhat-se.mult*x$se, x$fhat+se.mult*x$se), type="l", lty=lty, col=col, xlab=xlab, ylab=ylab, ...)
    else plot(x$argvals, x$fhat, type="l", col=col, xlab=xlab, ylab=ylab, ...)
}
