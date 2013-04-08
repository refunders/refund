decomp3d <- function(wd3Dobj, min.scale) {
	if (class(wd3Dobj) != "wd3D") stop("wd3Dobj must be of class 'wd3D'")
	l <- list(coef = as.vector(wd3Dobj$a), rowNum = 2 ^ wd3Dobj$nlevels, 
	          family = wd3Dobj$family, filter.number = wd3Dobj$filter.number)
	class(l) <- "decomp3d"
	return(l)
}

