plot.fosr <-
function(x, split=NULL, titles=NULL, xlabel="", ylabel="Coefficient function", set.mfrow=TRUE, ...) {
	nplots = nrow(x$B)
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

