amc <- function(y, Xmat, S, gam.method='REML', C=NULL, lambda=NULL, ...) {
	n.p = length(S)
	if (!is.null(C)) {
		# The following is based on Wood (2006), p. 186
	    n.con = dim(C)[1]
	    Z. = qr.Q(qr(t(C)), complete=TRUE)[ , -(1:n.con)]
	    Xmat. = Xmat %*% Z.
	    S. = vector("list", n.p)
	    for (i in 1:n.p) S.[[i]] = crossprod(Z., S[[i]] %*% Z.)
	}
	else {
		Z. = diag(ncol(Xmat))
		Xmat. = Xmat
		S. = S
	}

    fitter = if (length(y) > 10000) bam else gam
    if (is.null(lambda)) fitobj = fitter(y ~ Xmat.-1, method=gam.method, paraPen=list(Xmat.=S.), ...)
    else fitobj = fitter(y ~ Xmat.-1, paraPen=list(Xmat.=S.), sp=lambda, ...)
    
	lambdavec = if (!is.null(fitobj$full.sp)) fitobj$full.sp else fitobj$sp
	fullpen = 0
	for (i in 1:n.p) fullpen = lambdavec[i] * S.[[i]]
	list(gam = fitobj, 
	     coefficients = Z. %*% fitobj$coef, 
	     Vp = Z. %*% fitobj$Vp %*% t(Z.), 
	     GinvXT = Z. %*% solve(crossprod(Xmat.) + fullpen, t(Xmat.)),
	     method = gam.method)
}

