fpcr <- function(y, xfuncs = NULL, fdobj = NULL, ncomp=NULL, pve = 0.99, nbasis = NULL, basismat = NULL, 
                 penmat = NULL, argvals = NULL, covt = NULL, mean.signal.term = FALSE, 
                 spline.order = NULL, family = "gaussian", method = "REML", sp = NULL, 
                 pen.order = 2, cv1 = FALSE, nfold = 5, store.cv = FALSE, store.gam = TRUE, ...){
    require(mgcv)
    require(fda)
    require(MASS)
    # Error handling
    n <- length(y)
    do.cv <- FALSE
    if (!nfold %in% 1:n) 
        stop("Argument 'nfold' is invalid: must be an integer between 1 and ", n, ".") 
    if (!family %in% c("gaussian", "binomial")) 
        stop("Only 'gaussian' and 'binomial' models are implemented in the current version.")
    if (is.null(fdobj)) {
        if (!is.array(xfuncs) || !length(dim(xfuncs)) %in% 2:3)
    	    stop("xfuncs must either be a 2D or 3D array")  
        dim.sig <- length(dim(xfuncs)) - 1
        if (is.null(nbasis)) nbasis <- ifelse(dim.sig == 1, 40, 15)
        if (dim.sig == 1){
    	    if (!is.list(nbasis) && is.vector(nbasis))
    		    nbs <- matrix(nbasis, ncol = 1)
    	    else stop("for 1D predictors, 'nbasis' must be a vector")   
    	    siglength <- ncol(xfuncs)
        } else{
    	    if (is.list(nbasis) && length(nbasis) == 2)
    		    nbs <- cbind(rep(nbasis[[1]], length(nbasis[[2]])), rep(nbasis[[2]], 
    		                 each = length(nbasis[[1]])))
    	    else if (!is.list(nbasis) && is.vector(nbasis))
    		    nbs <- matrix(rep(nbasis,2), ncol=2)
    	    else stop("for 2D predictors, 'nbasis' must either be a vector or a list of length 2")
    	    d1 <- dim(xfuncs)[2L]
            d2 <- dim(xfuncs)[3L]
            siglength <- d1 * d2
        }
        if (cv1 || nrow(nbs) > 2 || (!is.null(ncomp) && length(ncomp) > 1)) do.cv <- TRUE
    } 
    else if (cv1 || (!is.null(ncomp) && length(ncomp) > 1)) do.cv <- TRUE 
    if (!do.cv) {
        store.cv <- FALSE
        nfold <- 1
    }
    groups <- split(sample(1:n), rep(1:nfold, length = n))   
    n.unpen.cols <- 1 + mean.signal.term + ifelse(is.null(covt), 0, ncol(as.matrix(covt)))
    cv.table <- array(0, dim = c(ifelse(is.null(fdobj),nrow(nbs), 1), ifelse(is.null(ncomp), 1, length(ncomp))))
    for (inb in 1 : dim(cv.table)[1]){
    	st <- fpcr.setup(y = y, xfuncs = xfuncs, fdobj = fdobj, nbasis = if (is.null(fdobj)) nbs[inb,] else NULL, 
    	                 basismat = basismat, penmat = penmat, argvals = argvals, covt = covt, 
    	                 mean.signal.term = mean.signal.term, spline.order = spline.order, 
    	                 pen.order = pen.order)
    	argvals <- st$argvals
    	if(is.null(fdobj)) cat("nbasis:", st$nbasis, "\n")
    	for (ifold in 1 : nfold){
    		if (do.cv) {
    		    idxTest <- groups[[ifold]]
    		    idxTrain <- (1:n)[-idxTest]
    		}
    		else idxTest <- idxTrain <- 1:n
    		X0.tr <- st$X0[idxTrain,]
    		SB.tr <- st$SB[idxTrain,]
    		svdSB <- svd(SB.tr)
    		if (is.null(ncomp)) ncomp <- min(which(cumsum(svdSB$d) > pve * sum(svdSB$d)))
    		for (incomp in 1 : length(ncomp)) {
    			V.ncomp <- svdSB$v[,1:ncomp[incomp]]
    			X <- cbind(X0.tr, SB.tr %*% V.ncomp)
    			S <- list(matrix(0, ncol(X), ncol(X)))
    			S[[1]][-(1:n.unpen.cols), -(1:n.unpen.cols)] <- 
    			                                 crossprod(V.ncomp, st$penmat %*% V.ncomp)
    			obje <- gam(y[idxTrain]~X - 1, paraPen = list(X=S), family = get(family), 
    			            method = method, sp = sp, ...)
    			BV <- st$basismat %*% V.ncomp
    			fhat <- BV %*% obje$coef[-(1:n.unpen.cols)]
    			undecor.coef <- obje$coef[1:n.unpen.cols] - 
    			                ginv(X0.tr) %*% st$xfuncs[idxTrain,] %*% fhat
    			yhat <- st$X0[idxTest,] %*% undecor.coef + st$xfuncs[idxTest,] %*% fhat
    			if (family == "gaussian")
    			    cv.table[inb, incomp] <- cv.table[inb, incomp] + 
    			                             mean((yhat - y[idxTest]) ^ 2)
    			else if (family == "binomial"){
                    phat <- exp(yhat) / (1 + exp(yhat))
   					phat <- replace(phat, exp(yhat) == Inf, 1)
   					cv.table[inb, incomp] <- cv.table[inb, incomp] + 
   					                         mean((phat > mean(y[idxTrain])) != y[idxTest]) 
   				}
    		}
    	}                 
    }
    if (do.cv) {
        idxmin <- which(cv.table == min(cv.table[cv.table!=0], na.rm = TRUE), arr.ind = TRUE)
        if (nrow(idxmin) > 1) idxmin <- idxmin[1,]
   	    if (is.list(nbasis)){
    	    dim(cv.table) <- c(length(nbasis[[1]]), length(nbasis[[2]]), length(ncomp))
            dimnames(cv.table) <- list(paste("nbasis1=", nbasis[[1]]), 
                                       paste("nbasis2=", nbasis[[2]]), paste("ncomp", ncomp))
        } 
        else dimnames(cv.table) <- list(paste("nbasis=", nbasis), paste("ncomp", ncomp))   
        if (dim.sig == 1) nbasis <- nbs[idxmin[1],]
        else {
        	nbasis <- list()
        	nbasis[[1]] <- nbs[idxmin[1],1]
        	nbasis[[2]] <- nbs[idxmin[1],2]
        }
        obje <- fpcr(y = y, xfuncs = xfuncs, fdobj = fdobj, ncomp = ncomp[idxmin[2]], 
                     nbasis = nbasis, basismat = basismat, penmat = penmat, 
                     argvals = argvals, covt = covt, mean.signal.term = mean.signal.term, 
                     spline.order = spline.order, family = family, method = method, sp = sp, 
                     pen.order = pen.order, cv1 = FALSE, store.gam = store.gam, ...)
        ncomp <- ncomp[idxmin[2]]
        if (store.cv) obje$cv.table <- cv.table
        else obje$cv.table <- min(cv.table[cv.table!=0])
    } 
    else {
        se <- sqrt(rowSums((BV %*% obje$Vp[-(1:n.unpen.cols), -(1:n.unpen.cols)]) * BV))
        if (!store.gam) {
        	yhat <- obje$fitted.values
        	obje <- list()
        	obje$fitted.values <- yhat
        }
        obje$se <- se
        obje$argvals <- argvals
    	obje$undecor.coef <- as.matrix(undecor.coef, nrow = 1)
    	colnames(obje$undecor.coef) <- if (is.null(dimnames(covt)) || is.null(dimnames(covt)[[2]])) 
                                          paste("X", 0:(n.unpen.cols-1), sep="")
    	                               else c("Intercept", dimnames(covt)[[2]])
        if (st$dim.sig == 2) dim(fhat) <- c(d1, d2)
        obje$fhat <- fhat
        obje$nbasis <- nbasis
        obje$ncomp <- ncomp
    }
    class(obje)="fpcr"
    return(obje)    
}
