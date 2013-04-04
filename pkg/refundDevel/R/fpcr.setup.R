fpcr.setup <- function(y, xfuncs = NULL, nbasis = NULL, basismat = NULL, penmat = NULL, argvals = NULL,  
                       covt = NULL, mean.signal.term = FALSE, spline.order = NULL, fdobj = NULL, pen.order = 2) {
    if (!is.null(fdobj)) {
    	# cat("performing scalar on functional regression\n")
    	basis <- fdobj$basis
    	if (!is.null(nbasis)){
    		if (length(nbasis) > 1){
    			warning("nbasis =", nbasis[-1], "will be ignored. Only the first nbasis, nbasis =", 
    			        nbasis[1],  "will be considerred")
    		}
    		if (nbasis[1] != basis$nbasis){
    			warning(paste("'nbasis =", nbasis, "'overridden since the supplied 'fdobj' 
    			              has basis dimension", basis$nbasis))
    		}
    	}
    	if (!is.null(spline.order)) {
    	    if(spline.order != norder(basis)){
    	    	warning(paste("'spline.order =", spline.order, "'overridden since the supplied 
    	    	              'fdobj' has a basis of order", norder(basis)))
    	    }
    	}
    	if (is.null(argvals)){
    		argvals <- seq(basis$rangeval[1], basis$rangeval[2], length.out = 401)
    	}
    	xfuncs <- t(eval.fd(argvals, fdobj))
    }
    
    if (!is.array(xfuncs) || !length(dim(xfuncs)) %in% 2:3) 
    	stop("xfuncs must either be a 3D or 2D array")
    
    dim.sig <- length(dim(xfuncs)) - 1    
    if (dim.sig == 1) {    	
    	if (is.null(fdobj)){
    		if (!is.null(nbasis) && length(nbasis) > 1) {
    		    warning("nbasis = ", nbasis[-1], " will be ignored. Only the first nbasis, nbasis = ", 
    			        nbasis[1], " will be considerred")
    		}
    		if (is.null(argvals)) argvals <- seq(0, 1, , ncol(xfuncs))
    		if (is.null(nbasis)) nbasis <- 40
    		if (is.null(spline.order)) spline.order <- 4
    		basis <- create.bspline.basis(rangeval = c(min(argvals), max(argvals)), 
    		                             nbasis = nbasis[1], norder = spline.order)
    	}
    	if (is.null(basismat)) basismat <- eval.basis(argvals, basis)
    	if (is.null(penmat)) penmat <- getbasispenalty(basis, pen.order)
    } 
    else {
    	d1 <- dim(xfuncs)[2]
    	d2 <- dim(xfuncs)[3]
    	if (is.null(argvals)) argvals <- list(seq(0, 1,, d1), seq(0,1,,d2))
    	if (is.null(nbasis)) {
    		nbasis <- c(15, 15)
    	} else if (length(nbasis) == 1){
    		warning("nbasis = ", nbasis[1], " will be applied to both direction")
    		nbasis <- c(nbasis, nbasis)
    	} else if (length(nbasis) > 2){
    	    warning("only nbasis = ", nbasis[1:2], " will be considerred")
    	    nbasis <- nbasis[1:2]	
        }
    	if (is.null(spline.order)) spline.order <- 4
    	basis1 <- create.bspline.basis(rangeval = c(min(argvals[[1]]), max(argvals[[1]])), 
    	                              nbasis = nbasis[1], norder = spline.order)
    	basis2 <- create.bspline.basis(rangeval = c(min(argvals[[2]]), max(argvals[[2]])), 
    	                              nbasis = nbasis[2], norder = spline.order)
    	if (is.null(basismat)){
    		basismat <- eval.basis(argvals[[2]], basis2) %x% eval.basis(argvals[[1]], basis1)
    	}
    	if (is.null(penmat)) penmat <- osplinepen2d(basis1, basis2)
    	dim(xfuncs) <- c(dim(xfuncs)[1], d1 * d2)
    }
    X0 <- if (mean.signal.term) cbind(1, apply(xfuncs, 1, mean))
          else matrix(1, length(y), 1)
    if (!is.null(covt)) X0 <- cbind(X0, as.matrix(covt))
    sigs.decor <- if (ncol(X0) == 1) scale(xfuncs, center = TRUE, scale = FALSE)
                  else lm(xfuncs ~ X0 - 1)$resid
    SB <- sigs.decor %*% basismat
    return(list(X0 = X0, SB = SB, penmat = penmat, basismat = basismat, xfuncs = xfuncs, 
                nbasis = nbasis, argvals = argvals, dim.sig = dim.sig))
}


