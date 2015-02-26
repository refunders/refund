reconstr3d <- function(decomp3dobj) {
	ctmp <- class(decomp3dobj)
	if (is.null(ctmp)) stop("decomp3dobj has no class.")
	else if (ctmp != "decomp3d") stop("decomp3dobj is not of class 'decomp3d'.")
	wd3dobj <- wd3D(array(0, dim = rep(decomp3dobj$rowNum, 3)), 
	                filter.number = decomp3dobj$filter.number,
	                family = decomp3dobj$family)
	wd3dobj$a <- decomp3dobj$coef
	dim(wd3dobj$a) <- rep(decomp3dobj$rowNum, 3)
	return(wr3D(wd3dobj))
}

