plot.wnet <- plot.wcr <- function(x, xlabel = "", ylabel = "Coefficient function", which.dim = 1, slices = NULL, set.mfrow = TRUE, image.axes = FALSE, ...){
	fhat <- x$fhat
	ndim.fhat = length(dim(fhat))
	if (ndim.fhat == 1){
		plot(fhat, type = 'l', xlab = xlabel, ylab = ylabel, ...)
	} else if (ndim.fhat == 2){
		image(fhat, axes = image.axes, ...)
	} else if (ndim.fhat == 3){
		if (which.dim == 2){
			fhat = aperm(fhat, c(2,1,3))
		} else if (which.dim == 3){
			fhat = aperm(fhat, c(3,1,2))
		}
		if (is.null(slices)){
			slices = 1 : dim(fhat)[1]
		}
		
		if (set.mfrow) {
			nplots = length(slices)
			nro = floor(sqrt(nplots))
			nco = ceiling(nplots / nro)
			par(mfrow = c(nro, nco))
		}
		for (i in slices){
			image(fhat[i,,], axes = image.axes, ...)
		}
	}
}