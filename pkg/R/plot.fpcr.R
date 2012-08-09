plot.fpcr = function(x, se=TRUE, col=1, lty=c(1,2,2), xlab="", ylab="Coefficient function", ...) {
    if (se) {
        if (is.numeric(se)) se.mult <- se
        else se.mult <- 2
        se.mult = max(se.mult, 0) 
    }
    if (se) matplot(x$argvals, cbind(x$fhat, x$fhat-se.mult*x$se, x$fhat+se.mult*x$se), type="l", lty=lty, col=col, xlab=xlab, ylab=ylab, ...)
    else plot(x$argvals, x$fhat, type="l", col=col, xlab=xlab, ylab=ylab, ...)
}
