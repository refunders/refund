checkError <- function(y, xfuncs, covt, vars){
	if (!is.array(xfuncs) || !length(dim(xfuncs)) %in% 2:4){
		stop("Argument <xfuncs> must be a 2D, 3D or 4D array.")
	}
	if (length(y) != dim(xfuncs)[1]){
		stop("Argument <y> and <xfuncs> are not of the same length.")
	}
	if ((!is.null(covt)) && (length(y) != dim(covt)[1])){
		stop("Argument <y> and <covt> are not of the same length.")
	}
	if (length(unique(dim(xfuncs)[-1])) != 1){
		stop("Argument <xfuncs> is invalid: dimensions are not identical.")
	}
	if (is.na(IsPowerOfTwo(dim(xfuncs)[2]))){
		stop("Argument <xfuncs> is invalid: each dimension must be of power of 2.")
	}
	for (i in 1 : length(vars)){
		varName = names(vars[i])
		varVal = vars[[i]]
        if ((varName == "min.scale") && 
            (max(varVal) > IsPowerOfTwo(dim(xfuncs)[2]) - 1)){
			stop("<min.scale> must in [0, ", IsPowerOfTwo(dim(xfuncs)[2]) - 1, "].")
		} else if ((varName == "nfeatures") && (!is.null(varVal)) && 
		           (max(varVal) > prod(dim(xfuncs)[-1]) + 1)){
		    stop("<nfeatures> must in [1, ", prod(dim(xfuncs)[-1]) + 1, "].")          	
		} else if ((varName == "alpha") && ((min(varVal) < 0) || (max(varVal) > 1))){
			stop("<alpha> must in [0, 1].")
		} 
	}
}

decorrelate <- function(xfuncs, covt){
    covt1 = cbind(1,as.matrix(covt))
    xfmat = matrix(xfuncs, nrow = dim(xfuncs)[1])
    xfuncs = xfuncs - 
    	     array(covt1 %*% solve(crossprod(covt1), crossprod(covt1, xfmat)), dim = dim(xfuncs))
    return(xfuncs)
}

waveletSetup <- function(xfuncs, filter.number, family, min.scale){
	dim.sig <- length(dim(xfuncs)[-1])
	wave.decomp <- switch(dim.sig, wd, imwd, wd3D)
    dec <- switch(dim.sig, decomp, decomp2d, decomp3d)
    rec <- switch(dim.sig, reconstr, reconstr2d, reconstr3d)  
    wdobjs <- apply(xfuncs, 1, wave.decomp, filter.number = filter.number, family = family)
    temp <- dec(wdobjs[[1]]) 
    p <- length(temp$coef)
    coef <- vector("list", length(min.scale))
    for (ims in 1 : length(min.scale)){
    	coef[[ims]] <- t(array(unlist(sapply(wdobjs, dec, min.scale = min.scale[ims])[1,]), dim = c(p, dim(xfuncs)[1])))
    	dimnames(coef[[ims]]) <- list(NULL, 1:p)
    }
	return(list(coef = coef, temp = temp, rec = rec))
}

waveletGetResult <- function(cv.table, family){
	nVars <- length(dim(cv.table)) - 2
	if (family == 'gaussian'){
		cv.splits <- apply(cv.table, c(1, 3:(nVars+2)), mean, na.rm = TRUE)
	} else if (family == 'binomial'){
		cv.splits <- apply(cv.table, c(1, 3:(nVars+2)), median, na.rm = TRUE)
	}
	se.cv <- apply(cv.splits, 2:(nVars+1), sd, na.rm = TRUE) / sqrt(dim(cv.table)[1])
	cv.table <- apply(cv.splits, 2:(nVars+1), mean, na.rm = TRUE)
	idxmin <- which(cv.table == min(cv.table[cv.table != 0], na.rm = TRUE), arr.ind = TRUE)[1,]
	return(list(cv.table = cv.table, se.cv = se.cv, idxmin = idxmin))
}

getRsq <- function(y, yhat, family){
	pp <- mean(y)
	if (family == 'gaussian'){
		Rsq <- 1 - var(y - yhat) / var(y)
	} else if (family == 'binomial'){
		logL0 <- length(y) * (pp * log(pp) + (1 - pp) * log(1 - pp))
		Rsq <- 1 - (sum(log(yhat)[y == 1]) + sum(log(1 - yhat)[y == 0])) / logL0
	}
	return(Rsq)
}
