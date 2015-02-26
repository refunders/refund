pwcv <- function(fdobj, Z, L=NULL, lambda, eval.pts=seq(min(fdobj$basis$range),max(fdobj$basis$range), length.out=201), scale=FALSE) {
    Z = scale(Z, center=FALSE, scale=scale)
	bss = fdobj$basis
    q = ncol(Z)

    J = getbasispenalty(bss, 0)
    svdJ = svd(J)
    J12 = svdJ$u %*% diag(sqrt(svdJ$d)) %*% t(svdJ$u)
    
    if (length(lambda) %in% c(1,q)) S = diag(lambda, q) %x% getbasispenalty(bss, 2)
    else stop("lambda must be either a scalar or a vector of length ncol(Z)")

    C = t(fdobj$coefs)
    N = NROW(C); K = NCOL(C)
    coefs.t = as.vector(J12 %*% t(C))

    if (!is.null(L)) {
        constr =  L %x% diag(bss$nbasis)  
	    n.con = dim(constr)[1]
	    Z. = qr.Q(qr(t(constr)), complete=TRUE)[ , -(1:n.con)]
	    X. = (Z %x% J12) %*% Z.
	    S. = crossprod(Z., S %*% Z.)
	}
	else {
		X. = Z %x% J12
		S. = S
	}

    A = X. %*% solve(crossprod(X.)+S., t(X.))
    resmat = t(matrix(coefs.t - A %*% coefs.t, K))

	discreps = matrix(NA, K, N)
	for (i in 1:N) {
	    ith = ((i-1)*K+1):(i*K)
	    discreps[ , i] = solve(diag(K)-A[ith,ith], resmat[i, ]) 
	}
	
	pw.discreps = eval.basis(eval.pts, bss) %*% solve(J12, discreps)
	pw.prederr = rowMeans(pw.discreps^2) 
	pw.prederr
}
