fosr.perm.fit <-
function(Y = NULL, fdobj = NULL, X, con = NULL, X0 = NULL, con0 = NULL, argvals = NULL, lambda = NULL, lambda0 = NULL, 
         multi.sp = FALSE, nperm, prelim, ...) {
    if (is.null(Y) == is.null(fdobj)) stop("Please specify 'Y' or 'fdobj', but not both") 
	if (length(lambda) > 1) stop("'lambda' must be a scalar or NULL")
	if (is.null(X0) & !is.null(con0)) stop("If 'con0' is given, 'X0' must be given too")
	df1 = ncol(X) - 1 
	if (!is.null(con)) df1 = df1 - nrow(con)
	df2 = nrow(X) - df1 - 1
	if (!is.null(X0)) df1 = df1 - ncol(X0) + 1
	if (!is.null(con0)) df1 = df1 + nrow(con0) 
	realfit = fosr(Y = Y, fdobj = fdobj, X = X, con = con, argvals = argvals, lambda = lambda, multi.sp = multi.sp, ...)
	lambda.real = realfit$lambda
	lambda0.real = NULL
	argvals = realfit$argvals
	if (!is.null(X0)) {
		realfit0 = fosr(Y = Y, fdobj = fdobj, X = X0, con = con0, argvals = argvals, lambda = lambda0, multi.sp = multi.sp, ...)
	    lambda0.real = realfit0$lambda
    }
    
    if (realfit$resp.typ == "fd"){
        yvec = eval.fd(argvals, fdobj)
        yhatvec = eval.fd(argvals, realfit$yhat)
    } else {
    	yvec = t(Y)
    	yhatvec = t(realfit$yhat)
    }
    
    if (is.null(X0)) F = (apply(yhatvec, 1, var) / df1) / (apply((yvec - yhatvec)^2, 1, mean) / df2)
    else {
    	yhatvec0 = if (realfit$resp.type == "fd") eval.fd(argvals, realfit0$yhat) else realfit0$yhat
    	F = ((apply((yvec - yhatvec0)^2, 1, mean) - apply((yvec - yhatvec)^2, 1, mean)) / df1) / (apply((yvec - yhatvec)^2, 1, mean) / df2)
    }

    lambda.prelim = lambda0.prelim = NULL
    if (prelim > 0) {
    	if (length(lambda) == 1) stop("If 'lambda' is specified, 'prelim' should be set to 0")
    	lambda.prelim = matrix(NA, prelim, ncol(X) ^ multi.sp)
    	if (!is.null(X0)) lambda0.prelim = matrix(NA, prelim, ncol(X0) ^ multi.sp)

    	begin.prelim = proc.time()
    	for (ee in 1:prelim) {
    		if (realfit$resp.type == "fd"){
    			fdobj.perm = fdobj
    			fdobj.perm$coefs = fdobj.perm$coefs[, sample(ncol(fdobj$coef))]
    			Y.perm = Y
    		} else {
    			Y.perm = Y[sample(nrow(Y)), ]
    			fdobj.perm = fdobj
    		}
    	    lambda.prelim[ee, ] = fosr(Y = Y.perm, fdobj = fdobj.perm , X = X, con = con, argvals=argvals, multi.sp=multi.sp, ...)$lambda 
    	    if (!is.null(X0)) lambda0.prelim[ee, ] = fosr(Y = Y.perm, fdobj = fdobj.perm, X = X0, con = con0, argvals=argvals,
    	                                             multi.sp=multi.sp, ...)$lambda
    	    if (ee==1) {
    	    	elapsed.time <- max(proc.time() - begin.prelim, na.rm = TRUE)
               if (elapsed.time > 10/prelim) cat("\n***** Estimated computing time for preliminary permutations:", round(prelim * elapsed.time), "seconds *****\n")
            }
    	}
    	lambda = apply(lambda.prelim, 2, median)
    	if (!is.null(X0)) lambda0 = apply(lambda0.prelim, 2, median)
        cat("***** Computing time for preliminary permutations:", max(proc.time() - begin.prelim, na.rm = TRUE), "seconds *****\n\n")
    }
    
    Y.perm = Y
    fdobj.perm = fdobj
    F.perm = matrix(NA, nperm, length(argvals))
    lambda.perm = matrix(NA, nperm, ncol(X)^multi.sp)
    lambda0.perm = if (!is.null(X0)) matrix(NA, nperm, ncol(X0)^multi.sp) else NULL

    begin.perm = proc.time()    
	for (i in 1:nperm) {
		if (i/20==floor(i/20)) cat('Permutation', i, '\n')
		if(realfit$resp.type == "fd"){
			fdobj.perm$coefs = fdobj.perm$coefs[, sample(ncol(fdobj$coefs))]
		} else {
			Y.perm = Y[sample(nrow(Y)), ]
		}
	    fit.perm = fosr(Y = Y.perm, fdobj = fdobj.perm, X = X, con = con, argvals=argvals, lambda=lambda, multi.sp=multi.sp, ...)
	    if (!is.null(X0)) fit.perm0 = fosr(Y = Y.perm, fdobj = fdobj.perm, X = X0, con = con0, argvals=argvals, lambda=lambda0,
	                                       multi.sp=multi.sp, ...)
	    lambda.perm[i, ] = fit.perm$lambda
	    if (!is.null(X0)) lambda0.perm[i, ] = fit.perm0$lambda
	    if (realfit$resp.type == "fd"){
	        yvec.perm = eval.fd(argvals, fdobj.perm)
            yhatvec.perm = eval.fd(argvals, fit.perm$yhat)   	
	    } else {
	    	yvec.perm = t(Y.perm)
	    	yhatvec.perm = t(fit.perm$yhat)
	    }
        if (is.null(X0)) F.perm[i, ] = (apply(yhatvec.perm, 1, var) / df1) / (apply((yvec.perm - yhatvec.perm)^2, 1, mean) / df2)
        else {
        	if (realfit$resp.type == "fd"){
        		yhatvec.perm0 = eval.fd(argvals, fit.perm0$yhat)
        	} else {
        		yhatvec.perm0 = fit.perm0$yhat
        	}
        	F.perm[i, ] = ((apply((yvec.perm - yhatvec.perm0)^2, 1, mean) - apply((yvec.perm - yhatvec.perm)^2, 1, mean)) / df1) / (apply((yvec.perm - yhatvec.perm)^2, 1, mean) / df2)
        }
    	if (i==1) {
    	    elapsed.time <- max(proc.time() - begin.perm, na.rm = TRUE)
            if (elapsed.time > 10/nperm) cat("\n***** Estimated computing time for permuted-data models:", round(nperm * elapsed.time), "seconds *****\n")
        }
    }
    cat("***** Computing time for permuted-data models:", max(proc.time() - begin.perm, na.rm = TRUE), "seconds *****\n\n")
    ll = list(F=F, F.perm=F.perm, argvals=argvals, lambda.real=lambda.real, lambda.prelim=lambda.prelim, lambda.perm=lambda.perm, lambda0.real=lambda0.real, lambda0.prelim=lambda0.prelim, lambda0.perm=lambda0.perm)
    class(ll) = "fosr.perm"
    ll
}

