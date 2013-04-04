wcr.perm <- function(y, xfuncs, min.scale, nfeatures, ncomp, method = c("pcr", "pls"), 
                     covt = NULL, nperm = 20, perm.method = "responses", ...){
    if (is.null(covt) && perm.method == "x.residuals"){
    	stop("'x.residuals' method is unavailable when 'covt' is NULL.")
    }
    cv <- wcr(y = y, xfuncs = xfuncs, min.scale = min.scale, nfeatures = nfeatures, ncomp = ncomp, 
              method = method, covt = covt, store.glm = FALSE, ...)$cv
    if (perm.method == "y.residuals") {
        obje <- wcr(y = y, xfuncs = xfuncs, min.scale = min.scale, nfeatures = nfeatures, ncomp = ncomp, 
                    method = method, covt = covt, store.glm = FALSE, ...)
        y.resid <- obje$fitted - y
    } 
    else if (perm.method == "x.residuals") {
    	X = covt
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
        cv.perm[i] <- min(wcr(y = yperm, xfuncs = xperm, min.scale = min.scale, nfeatures = nfeatures, 
                         ncomp = ncomp, method = method, covt = covt, store.glm = FALSE, ...)$cv)
    }
    pvalue <- (1 + sum(cv.perm < cv)) / (1 + nperm)
    list(cv = cv, cv.perm = cv.perm, pvalue = pvalue)
}


