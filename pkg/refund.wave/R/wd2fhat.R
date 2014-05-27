wd2fhat <-
function(est, info){
	info$temp$coef <- rep(0, length(info$temp$coef))
	info$temp$coef[as.numeric(names(est)[-(1:(1+info$ncovt))])] <- 
	                                     est[-(1:(1+info$ncovt))]
	fhat <- as.vector(info$rec(info$temp))
	dim(fhat) = rep(info$d, info$dim.sig)
	return(fhat)
}
