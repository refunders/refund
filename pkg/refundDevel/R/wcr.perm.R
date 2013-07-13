wcr.perm <- function(y, xfuncs, min.scale = 0, nfeatures, ncomp, method = c("pcr", "pls"), 
                      covt = NULL, nrep = 1, nsplit=1, nfold = 5, nperm = 20, 
                      perm.method = c("responses", "y.residuals", "x.residuals"),
                      family = "gaussian", seed.real=NULL, seed.perm=NULL, ...){
    perm.method = match.arg(perm.method)
    if (is.null(covt) && perm.method == "x.residuals"){
    	stop("'x.residuals' method is unavailable when 'covt' is NULL.")
    }
    cat("******* Real-data model *******\n")
    replicate_count <- 1
    res <- apply(replicate(nrep, expr = {
    	            if (nrep != 1) cat("replicate", replicate_count, "\n")
    	            replicate_count <<- replicate_count + 1
    	            obj <- wcr(y = y, xfuncs = xfuncs, min.scale = min.scale, nfeatures = nfeatures, 
    	                       ncomp = ncomp, method = method, covt = covt, family = family,  
                               nsplit = nsplit, nfold = nfold, 
                               seed = if(!is.null(seed.real)) seed.real * replicate_count else NULL, ...)
                    c(obj$cv.table, obj$fhat)
                    }), 1, mean, na.rm=TRUE)
    if (perm.method == "y.residuals") {
        obje <- wcr(y = y, xfuncs = xfuncs, min.scale = min.scale, nfeatures = nfeatures, ncomp = ncomp, 
                    method = method, covt = covt, family = family, nsplit = nsplit, ...)
        y.resid <- obje$fitted - y
    }  else if (perm.method == "x.residuals") {
    	X = as.matrix(covt)
    	Y = matrix(xfuncs, nrow = length(y))
        XtX.inv = solve(crossprod(X))
        coef = XtX.inv %*% crossprod(X, Y)
        fitted= X %*% coef
    	x.resid = xfuncs - array(fitted, dim = dim(xfuncs))
    }
    cat("***** Permuted-data models *****\n")
    cv.perm <- rep(0, nperm)
    if (!is.null(seed.perm)) set.seed(seed.perm)
    for (i in 1 : nperm){
    	cat("perm", i, "\n")
        if (perm.method == "responses"){
        	    yperm <- sample(y)
        	    xperm <- xfuncs
        } else if (perm.method == "y.residuals"){
            yperm <- obje$fitted + sample(y.resid)
            xperm <- xfuncs
        } else if (perm.method == "x.residuals"){
          	yperm <- y
            if (length(dim(xfuncs)) == 3) xperm <- xfuncs + x.resid[sample(1:dim(xfuncs)[1]),,]
            else if (length(dim(xfuncs)) == 4) xperm <- xfuncs + x.resid[sample(1:dim(xfuncs)[1]),,,]       	
        }
        obj <- wcr(y = yperm, xfuncs = xperm, min.scale = min.scale, nfeatures = nfeatures, ncomp = ncomp, 
                   method = method, covt = covt, family = family, nsplit = nsplit, nfold = nfold,
                   seed = if(!is.null(seed.perm)) seed.perm * i else NULL, ...)
        cv.perm[i] <- min(obj$cv.table, na.rm=TRUE)                     
    }    
    pvalue <- (1 + sum(cv.perm < res[1])) / (1 + nperm)    
    list(cv = res[1], cv.perm = cv.perm, pvalue = pvalue)                  
}
