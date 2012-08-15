# Incorporates Fabian's pcre idea for random effect functions (method="mix")
fosr <- function (Y=NULL, fdobj=NULL, X, con = NULL, argvals = NULL, 
        method = c("OLS","GLS","mix"),
        gam.method = c("REML", "ML", "GCV.Cp", "GACV.Cp", "P-REML", "P-ML"), 
        cov.method = c("naive", "mod.chol"),
        lambda = NULL, nbasis=15, norder=4, 
        pen.order=2, multi.sp = ifelse(method=="OLS", FALSE, TRUE), pve=.99,
        max.iter = 1, maxlam = NULL, cv1 = FALSE, scale = FALSE)
{
    if (is.null(Y)==is.null(fdobj)) stop("Please specify 'Y' or 'fdobj', but not both") 
    resp.type <- if (is.null(Y)) "fd" else "raw"
    if (is.null(argvals)) 
        argvals <- if (is.null(fdobj)) seq(0,1, length=ncol(Y)) 
                else seq(min(fdobj$basis$range), max(fdobj$basis$range), length=201)                
    method <- match.arg(method)
    cov.method <- match.arg(cov.method)
    gam.method <- match.arg(gam.method)
    
    if (method != "OLS" & (length(lambda) > 1))
        stop("Vector-valued lambda allowed only if method = 'OLS'")
    if (!is.null(lambda) & multi.sp)
        stop("Fixed lambda not implemented with multiple penalties")
    if (method == "OLS" & multi.sp)
        stop("OLS not implemented with multiple penalties")
    
    if (resp.type=="raw") {
        bss = create.bspline.basis(range(argvals), nbasis=nbasis, norder=norder)
        bmat <- Theta <- eval.basis(argvals, bss)
        respmat <- Y 
    }
    else if (resp.type=="fd") {
        if (!is.fd(fdobj)) stop("'fdobj' must be a functional data object")
        bss = fdobj$basis
        nbasis = bss$nbasis
        Theta <- eval.basis(argvals, bss)
        C = t(fdobj$coefs)
        J = getbasispenalty(bss, 0)
        svdJ = svd(J)
        bmat <- J12 <- svdJ$u %*% diag(sqrt(svdJ$d)) %*% t(svdJ$u)
        respmat <- C %*% J12
    }
    
    newfit = U = pca.resid = NULL
    X.sc = scale(X, center = FALSE, scale = scale)
    q = ncol(X)
    
    ncurve <- nrow(respmat)
    if (multi.sp) {
        pen = vector("list", q)
        for (j in 1:q) {
            one1 = matrix(0, q, q)
            one1[j, j] = 1
            pen[[j]] = one1 %x% getbasispenalty(bss, pen.order)
        }
    }
    else pen = list(diag(q) %x% getbasispenalty(bss, pen.order))
    
    constr = if (!is.null(con)) con %x% diag(nbasis) else NULL
    cv = NULL
    if (method == "OLS") {
        if (length(lambda) != 1 | cv1) {
            lofo <- lofocv(respmat, X.sc %x% bmat, S1 = pen[[1]],
                    lamvec = lambda, constr = constr, maxlam = maxlam)
            cv = if (is.null(lambda))
                        lofo$objective
                    else min(lofo[, 2])
            lambda = if (is.null(lambda))
                        lofo$min
                    else lofo[which.min(lofo[, 2]), 1]
        }
    }
    
    firstfit <- amc(as.vector(t(respmat)), X.sc %x% bmat,
            gam.method = gam.method, S = pen, C = constr, lambda = lambda)
    
    B = B.ols = t(matrix(firstfit$coef, ncol = q))
    se = NULL
    if (method != "OLS") {
        iter = 0
        B.old = 3 * B.ols
        newfit = NULL
        if (!is.null(lambda) & max.iter > 0)
            warning("Given lambda used for initial fit only")
        while (any(abs((B - B.old)/B.old) > 0.001) & (iter < max.iter)) {
            iter = iter + 1
            if (max.iter > 1) cat("Refit", iter, "\n")
            oldfit = if (!is.null(newfit)) newfit else firstfit
            B.old = B
            residvec <- as.vector(t(respmat)) - (X.sc %x% bmat) %*% oldfit$coef[1:(q*nbasis)]
            residmat = t(matrix(residvec, ncol = ncurve))
            if (method == "GLS") {
                # Estimate symmetric square root of the precision matrix
                if (cov.method=="mod.chol") {
                    # browser() 
                    p = ncol(residmat)
                    res.cent = scale(residmat, TRUE, FALSE)
                    sqrt.prec.list = list()
                    lwstat = lwpval = c()
                    for (nband in 1:(p-1)) {
                        # cat(nband, if (nband==1) "band...\n" else "bands...\n")
                        TT = diag(p); Ddiag = rep(0,p)
                        Ddiag[1] = var(res.cent[ , 1])
                        for (k in 2:p) {
                            qrResCent <- qr(res.cent[ , max(1,k-nband):(k-1)])
                            TT[k, max(1,k-nband):(k-1)] <- (-qr.coef(qrResCent, res.cent[ , k]))
                            Ddiag[k] <- var(qr.resid(qrResCent, res.cent[ , k]))
                        }
                        prec = scale(t(TT), FALSE, Ddiag) %*% TT
                        svdp = eigen(prec, symmetric=TRUE)
                        sqrt.prec.list[[nband]] = svdp$vectors %*% tcrossprod(diag(sqrt(svdp$values)), svdp$vectors)
                        lwprec = lw.test(residmat %*% sqrt.prec.list[[nband]])	
                        lwstat[nband] = lwprec$stat; lwpval[nband] = lwprec$pvalue  
                        if (lwstat[nband] < -5) break      
                    }
                    
                    nband.best = which.max(lwpval)
                    cat("Using half-bandwidth", nband.best, "for precision matrix of residuals\n")
                    sqrt.prec <- sqrt.prec.list[[nband.best]]
                } else if (cov.method=="naive") {
                    if (nrow(residmat) < ncol(residmat)) stop("Sample covariance matrix of residuals is singular.")
                    svd.cov.mle <- svd(cov(residmat) * (ncurve-1)/ncurve)
                    sqrt.prec <- tcrossprod(scale(svd.cov.mle$u, FALSE, sqrt(svd.cov.mle$d)), svd.cov.mle$u)
                }    
                    # Next line comes from eq. (18) of Reiss et al. (2010)
                    #F start iterations at last solution for quicker convergence
                    #F (maybe? not sure this has much effect, but it surely won't hurt):
                    newfit <- amc(as.vector(tcrossprod(sqrt.prec, respmat)), 
                            X.sc %x% (sqrt.prec %*% bmat), 
                            gam.method = gam.method, S = pen, C = constr, 
                            start = if (is.null(con)) as.vector(t(B)) else NULL)
                    B = t(matrix(newfit$coef, ncol = q))	
                } #end GLS
                else if (method == "mix") {                     
                    if (resp.type=="fd") {
                        resid.fd <- fd(solve(J12, t(residmat)), bss)
                        if (iter==1) {
                            pca.resid <- pca.fd(resid.fd, nharm=min(ncurve-1, nbasis)) 
                            npc <- min(which(cumsum(pca.resid$varprop) > pve))
                        }
                        else pca.resid <- pca.fd(resid.fd, nharm=npc)
                        evalues <- pca.resid$values[1:npc]
                        efuncmat.scaled <- bmat %*% t(t(pca.resid$harmonics$coef[ , 1:npc]) * sqrt(evalues))
                    }
                    else if (resp.type=="raw") {
                        if (iter==1) {                           
                            pca.resid <- fpca.sc(residmat, pve=pve) 
                            npc <- pca.resid$npc
                        }
                        else pca.resid <- fpca.sc(residmat, npc=npc)
                        evalues <- pca.resid$evalues
                        efuncmat.scaled <- t(t(pca.resid$efunctions) * sqrt(evalues))
                    }
                    if (iter==1) cat("Using leading", npc, "PCs of residual functions for random effects\n")
                    
                    npen <- length(pen); pendim <- ncol(pen[[1]])
                    pen.aug = vector("list", npen+1)
                    for (l in 1:npen) {
                        pen.aug[[l]] <- matrix(0, pendim+npc*ncurve, pendim+npc*ncurve)
                        pen.aug[[l]][1:pendim, 1:pendim] <- pen[[l]]
                    }
                    pen.aug[[npen+1]] <- diag(rep(0:1, c(pendim, npc*ncurve))) 
                    constr.aug <- if (is.null(constr)) NULL else cbind(constr, matrix(0, nrow(constr), npc*ncurve))
                    #F
                    startB <- if(iter==1){
                                c(as.vector(t(B)), rep(0, ncurve*npc))
                            } else {
                                newfit$coefficients
                            }
                    #/F
                    #browser()

                    newfit <- amc(as.vector(t(respmat)), 
                            cbind(X.sc %x% bmat,diag(ncurve) %x% efuncmat.scaled),
                            gam.method = gam.method, S = pen.aug, C = constr.aug, 
                            start=if (is.null(constr.aug)) startB else NULL)      
                    
                    
                    vecBt = newfit$coef[1:(q*nbasis)]
                    vecUt = newfit$coef[(q*nbasis+1):(q*nbasis+npc*ncurve)]
                    B = t(matrix(vecBt, ncol=q))
                    U <- t(matrix(vecUt, ncol=ncurve))
                } #end mix
        } #end while 
    } # end !OLS
    
    if (method == "OLS" | max.iter == 0) {
        residvec <- as.vector(t(respmat)) - (X.sc %x% bmat) %*% firstfit$coef
        covmat = ((ncurve - 1)/ncurve) * cov(t(matrix(residvec, ncol = ncurve)))
        var.b = firstfit$GinvXT %*% (diag(ncurve) %x% covmat) %*% t(firstfit$GinvXT)
    }
    else var.b = newfit$Vp
    se.func = matrix(NA, length(argvals), q)
    for (j in 1:q) {
        # Pointwise SE estimate for j-th coefficient function,
        # derived from variance of basis coefficients as given by
        # eq. (23) of Reiss et al. (2010)
        se.func[ , j] = sqrt(rowSums((Theta %*% var.b[(nbasis * (j - 1) + 1):(nbasis * j), 
                                            (nbasis * (j - 1) + 1):(nbasis * j)]) * Theta))
    }
    est.func = eval.fd(argvals, fd(t(B), bss))
    fit <- if (method == "mix" & max.iter > 0) newfit else firstfit
    roughness = diag(B %*% getbasispenalty(bss, pen.order) %*% t(B))
    skale = attr(X.sc, "scaled:scale")
    if (!is.null(skale)) {
        B = t(scale(t(B), center = FALSE, scale = skale))
        est.func = scale(est.func, center = FALSE, scale = skale)
        se.func = scale(se.func, center = FALSE, scale = skale)
        roughness = roughness/skale^2
    }
    yhat = if (resp.type=="raw") X %*% tcrossprod(B, Theta) else fd(t(X %*% B), bss)
    llist = list(B = B, U = U, pca.resid = pca.resid,
            yhat = yhat, resid = if (resp.type == 'raw') Y - yhat else fdobj - yhat,
            est.func = est.func,  se.func = se.func, argvals = argvals, fit = fit, 
            edf = sum(fit$gam$edf), 
            lambda = if (length(fit$gam$sp) > 0) fit$gam$sp else fit$gam$full.sp,
            cv = cv, roughness = roughness, resp.type = resp.type)
    class(llist) = "fosr"
    llist
}
