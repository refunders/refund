##' Default plotting of function-on-scalar regression objects
##'
##' Plots the coefficient function estimates produced by \code{fosr()}.
##'
##'
##' @param x an object of class \code{"\link{fosr}"}.
##' @param split value, or vector of values, at which to divide the set of
##' coefficient functions into groups, each plotted on a different scale.
##' E.g., if set to 1, the first function is plotted on one scale, and all
##' others on a different (common) scale.  If \code{NULL}, all functions are
##' plotted on the same scale.
##' @param titles character vector of titles for the plots produced, e.g.,
##' names of the corresponding scalar predictors.
##' @param xlabel label for the x-axes of the plots.
##' @param ylabel label for the y-axes of the plots.
##' @param set.mfrow logical value: if \code{TRUE}, the function will try to
##' set an appropriate value of the \code{mfrow} parameter for the plots.
##' Otherwise you may wish to set \code{mfrow} outside the function call.
##' @param \dots graphical parameters (see \code{\link{par}}) for the plot.
##' @author Philip Reiss \email{phil.reiss@@nyumc.org}
##' @seealso \code{\link{fosr}}, which includes examples.
##' @export
plot.fosr <-
function(x, split=NULL, titles=NULL, xlabel="", ylabel="Coefficient function", set.mfrow=TRUE, ...) {
	nplots = ncol(x$fd$coef)
	if (set.mfrow) {
	    nro = floor(sqrt(nplots))
	    nco = ceiling(nplots / nro)
	    par(mfrow=c(nro, nco))
	}
	firsts = c(1, split+1)
	lasts = c(split, nplots)
	ngps = length(firsts)
	rng = matrix(NA, ngps, 2)
	for (i in 1:ngps) {
	    rng[i, ] = c(min(x$est[ , firsts[i]:lasts[i]]-2*x$se[ , firsts[i]:lasts[i]]), max(x$est[ , firsts[i]:lasts[i]]+2*x$se[ , firsts[i]:lasts[i]]))
	    for (k in firsts[i]:lasts[i]) {
	    	plot(x$argvals, x$est[ , k], type='l', ylim=rng[i, ], main=titles[k], xlab=xlabel, ylab=ylabel, ...)
	    	lines(x$argvals, x$est[ , k]-2*x$se[ , k], lty=3, lwd=1.5)
	    	lines(x$argvals, x$est[ , k]+2*x$se[ , k], lty=3, lwd=1.5)
	    	abline(h=0, col='grey')
	    }
	}
}
