fosr2s <- function (Y, X, argvals = seq(0,1,,ncol(Y)), nbasis = 15, norder = 4, 
                       pen.order=norder-2, basistype="bspline") 
{
    # Stage 1: raw estimates
    n = dim(X)[1]
    p = dim(X)[2]
    XtX.inv = solve(crossprod(X))
    raw.coef = t(XtX.inv %*% crossprod(X, Y))
    resmat = Y - X %*% t(raw.coef)
    covmat = cov(resmat)  # TODO: add more sophisticated covariance methods?
    sigma2 = apply(resmat, 2, crossprod) / (n-p)
    raw.se = t(sqrt(diag(XtX.inv) %o% sigma2))

    # Stage 2: smooth raw estimates to get coefficient functions
    if (basistype=="bspline") bss = create.bspline.basis(range(argvals), nbasis = nbasis, norder = norder)
    else if (basistype=="fourier") bss = create.fourier.basis(range(argvals), nbasis = nbasis)
    
    Bmat = eval.basis(argvals, bss)
    P = getbasispenalty(bss, pen.order)
    Bt.Sig.B <- crossprod(Bmat, covmat %*% Bmat)
    
    lambda=c()
    coefmat <- matrix(NA, nbasis, p)
    est.func <- se.func <- matrix(NA,length(argvals),p)
    for (j in 1:p) {
    	swmod <- gam(raw.coef[ , j] ~ Bmat-1, paraPen=list(Bmat=list(P)), method="REML")
    	lambda[j] <- swmod$sp
    	coefmat[ , j] <- swmod$coef
    	est.func[ , j] <- fitted(swmod)
    	m1 <- Bmat %*% solve(crossprod(Bmat) + swmod$sp * P)
    	se.func[ , j] <- sqrt(XtX.inv[j,j] * rowSums(m1 * (m1 %*% Bt.Sig.B)))
    }
    
    list(fd=fd(coef=coefmat, basisobj=bss), raw.coef=raw.coef, raw.se=raw.se,
         yhat=tcrossprod(X, est.func),
         est.func=est.func, se.func=se.func, argvals=argvals, lambda=lambda)
}
