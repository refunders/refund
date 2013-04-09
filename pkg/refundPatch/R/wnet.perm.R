wnet.perm <- function(y, xfuncs, min.scale, nfeatures, alpha, lambda = NULL, covt = NULL, nsplit = 10, nperm = 20, 
                      perm.method = c("responses", "y.residuals", "x.residuals"), family = "gaussian",...){
    perm.method = match.arg(perm.method)
    if (is.null(covt) && perm.method == "x.residuals"){
    	stop("'x.residuals' method is unavailable when 'covt' is NULL.")
    }
    cv <- mean(replicate(nsplit, expr = {
    	            wnet(y = y, xfuncs = xfuncs, min.scale = min.scale, nfeatures = nfeatures, alpha = alpha, 
    	                 lambda = lambda, covt = covt, family = family, ...)$cv.table
                    }))
    if (perm.method == "y.residuals") {
        obje <- wnet(y = y, xfuncs = xfuncs, min.scale = min.scale, nfeatures = nfeatures, alpha = alpha, 
                     lambda = lambda, covt = covt, ...)
        y.resid <- obje$fitted - y
    }  else if (perm.method == "x.residuals") {
    	X = as.matrix(covt)
    	Y = matrix(xfuncs, nrow = length(y))
        XtX.inv = solve(crossprod(X))
        coef = XtX.inv %*% crossprod(X, Y)
        fitted= X %*% coef
    	x.resid = xfuncs - array(fitted, dim = dim(xfuncs))
    }
    cv.perm <- rep(0, nperm)
    for (i in 1 : nperm){
    	cat("nperm:", i, "\n")
        if (perm.method == "responses"){
        	yperm <- sample(y)
        	xperm <- xfuncs
        } else if (perm.method == "y.residuals"){
            yperm <- obje$fitted + sample(y.resid)
            xperm <- xfuncs
        } else if (perm.method == "x.residuals"){
        	yperm <- y
        	xperm <- xfuncs + x.resid[sample(1:dim(xfuncs)[1]),,]
        }
        cv.perm[i] <- min(wnet(y = yperm, xfuncs = xperm, min.scale = min.scale, nfeatures = nfeatures, alpha = alpha, 
                               lambda = lambda, covt = covt, family = family, ...)$cv.table)                      
    }    
    pvalue <- (1 + sum(cv.perm < cv)) / (1 + nperm)
    list(cv = cv, cv.perm = cv.perm, pvalue = pvalue)                   	
}
