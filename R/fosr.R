##' Function-on-scalar regression
##'
##' Fit linear regression with functional responses and scalar predictors, with
##' efficient selection of optimal smoothing parameters.
##'
##' The GLS method requires estimating the residual covariance matrix, which
##' has dimension \eqn{d\times d} when the responses are given by \code{Y}, or
##' \eqn{nbasis\times nbasis} when they are given by \code{fdobj}. When
##' \code{cov.method = "naive"}, the ordinary sample covariance is used. But
##' this will be singular, or nonsingular but unstable, in high-dimensional
##' settings, which are typical. \code{cov.method = "mod.chol"} implements the
##' modified Cholesky method of Pourahmadi (1999) for estimation of covariance
##' matrices whose inverse is banded. The number of bands is chosen to maximize
##' the p-value for a sphericity test (Ledoit and Wolf, 2002) applied to the
##' "prewhitened" residuals. Note, however, that the banded inverse covariance
##' assumption is sometimes inappropriate, e.g., for periodic functional
##' responses.
##'
##' There are three types of values for argument \code{lambda}:\enumerate{\item
##' if \code{NULL}, the smoothing parameter is estimated by
##' \code{\link[mgcv]{gam}} (package \pkg{mgcv}) if \code{method = "GLS"}, or
##' by \code{optimize} if \code{method = "OLS"}; \item if a scalar, this value
##' is used as the smoothing parameter (but only for the initial model, if
##' \code{method = "GLS"}); \item if a vector, this is used as a grid of values
##' for optimizing the cross-validation score (provided \code{method = "OLS"};
##' otherwise an error message is issued).}
##'
##' Please note that currently, if \code{multi.sp = TRUE}, then \code{lambda}
##' must be \code{NULL} and \code{method} must be \code{"GLS"}.
##'
##' @param formula Formula for fitting fosr. If used, data argument must not be null.
##' @param Y,fdobj the functional responses, given as either an \eqn{n\times d}
##' matrix \code{Y} or a functional data object (class \code{"\link[fda]{fd}"})
##' as in the \pkg{fda} package.
##' @param data data frame containing the predictors and responses.
##' @param X the model matrix, whose columns represent scalar predictors.
##' Should ordinarily include a column of 1s.
##' @param con a row vector or matrix of linear contrasts of the coefficient
##' functions, to be constrained to equal zero.
##' @param argvals the \eqn{d} argument values at which the coefficient
##' functions will be evaluated.
##' @param method estimation method: \code{"OLS"} for penalized ordinary least
##' squares, \code{"GLS"} for penalized generalized least squares, \code{"mix"}
##' for mixed effect models.
##' @param gam.method smoothing parameter selection method, to be passed to
##' \code{\link[mgcv]{gam}}: \code{"REML"} for restricted maximum likelihood,
##' \code{"GCV.Cp"} for generalized cross-validation.
##' @param cov.method covariance estimation method: the current options are
##' naive or modified Cholesky. See Details.
##' @param lambda smoothing parameter value.  If \code{NULL}, the smoothing
##' parameter(s) will be estimated.  See Details.
##' @param nbasis,norder number of basis functions, and order of splines (the
##' default, 4, gives cubic splines), for the B-spline basis used to represent
##' the coefficient functions. When the functional responses are supplied using
##' \code{fdobj}, these arguments are ignored in favor of the values pertaining
##' to the supplied object.
##' @param pen.order order of derivative penalty.
##' @param multi.sp a logical value indicating whether separate smoothing
##' parameters should be estimated for each coefficient function.  Currently
##' must be \code{FALSE} if \code{method = "OLS"}.
##' @param pve if \code{method = 'mix'}, the percentage of variance explained
##' by the principal components; defaults to 0.99.
##' @param max.iter maximum number of iterations if \code{method = "GLS"}.
##' @param maxlam maximum smoothing parameter value to consider (when
##' \code{lamvec=NULL}; see \code{\link{lofocv}}).
##' @param cv1 logical value indicating whether a cross-validation score should
##' be computed even if a single fixed \code{lambda} is specified (when
##' \code{method = "OLS"}).
##' @param scale logical value or vector determining scaling of the matrix
##' \code{X} (see \code{\link{scale}}, to which the value of this argument is
##' passed).
##' @return An object of class \code{fosr}, which is a list with the following
##' elements: \item{fd}{object of class \code{"\link{fd}"} representing the
##' estimated coefficient functions. Its main components are a basis and a
##' matrix of coefficients with respect to that basis.} \item{pca.resid}{if
##' \code{method = "mix"}, an object representing a functional PCA of the
##' residuals, performed by \code{\link{fpca.sc}} if the responses are in raw
##' form or by \code{\link[fda]{pca.fd}} if in functional-data-object form.}
##' \item{U}{if \code{method = "mix"}, an \eqn{n\times m} matrix of random
##' effects, where \eqn{m} is the number of functional PC's needed to explain
##' proportion \code{pve} of the residual variance. These random effects can be
##' interpreted as shrunken FPC scores.} \item{yhat, resid}{objects of the same
##' form as the functional responses (see arguments \code{Y} and \code{fdobj}),
##' giving the fitted values and residuals.} \item{est.func}{matrix of values
##' of the coefficient function estimates at the points given by
##' \code{argvals}.} \item{se.func}{matrix of values of the standard error
##' estimates for the coefficient functions, at the points given by
##' \code{argvals}.} \item{argvals}{points at which the coefficient functions
##' are evaluated.} \item{fit}{fit object outputted by \code{\link{amc}}.}
##' \item{edf}{effective degrees of freedom of the fit.}
##' \item{lambda}{smoothing parameter, or vector of smoothing parameters.}
##' \item{cv}{cross-validated integrated squared error if \code{method="OLS"},
##' otherwise \code{NULL}.} \item{roughness}{value of the roughness penalty.}
##' \item{resp.type}{\code{"raw"} or \code{"fd"}, indicating whether the
##' responses were supplied in raw or functional-data-object form.}
##' @author Philip Reiss \email{phil.reiss@@nyumc.org}, Lan Huo, and Fabian
##' Scheipl
##' @seealso \code{\link{plot.fosr}}
##' @references Ledoit, O., and Wolf, M. (2002). Some hypothesis tests for the
##' covariance matrix when the dimension is large compared to the sample size.
##' \emph{Annals of Statistics}, 30(4), 1081--1102.
##'
##' Pourahmadi, M. (1999). Joint mean-covariance models with applications to
##' longitudinal data: unconstrained parameterisation. \emph{Biometrika},
##' 86(3), 677--690.
##'
##' Ramsay, J. O., and Silverman, B. W. (2005).  \emph{Functional Data
##' Analysis}, 2nd ed., Chapter 13.  New York: Springer.
##'
##' Reiss, P. T., Huang, L., and Mennes, M. (2010).  Fast function-on-scalar
##' regression with penalized basis expansions.  \emph{International Journal of
##' Biostatistics}, 6(1), article 28.  Available at
##' \url{http://works.bepress.com/phil_reiss/16/}
##' @examples
##' \dontrun{
##' require(fda)
##' # The first two lines, adapted from help(fRegress) in package fda,
##' # set up a functional data object representing daily average
##' # temperatures at 35 sites in Canada
##' daybasis25 <- create.fourier.basis(rangeval=c(0, 365), nbasis=25,
##'                   axes=list('axesIntervals'))
##' Temp.fd <- with(CanadianWeather, smooth.basisPar(day.5,
##'                 dailyAv[,,'Temperature.C'], daybasis25)$fd)
##'
##' modmat = cbind(1, model.matrix(~ factor(CanadianWeather$region) - 1))
##' constraints = matrix(c(0,1,1,1,1), 1)
##'
##' # Penalized OLS with smoothing parameter chosen by grid search
##' olsmod = fosr(fdobj = Temp.fd, X = modmat, con = constraints, method="OLS", lambda=100*10:30)
##' plot(olsmod, 1)
##' 
##' # Test use formula to fit fosr
##' set.seed(2121)
##' data1 <- pffrSim(scenario="ff", n=40)
##' formod = fosr(Y~xlin+xsmoo, data=data1)
##' plot(formod, 1)
##' 
##' # Penalized GLS
##' glsmod = fosr(fdobj = Temp.fd, X = modmat, con = constraints, method="GLS")
##' plot(glsmod, 1)
##' }
##' @importFrom fda create.bspline.basis eval.basis getbasispenalty fd pca.fd is.fd
##' @export

fosr <- function (formula=NULL, Y=NULL, fdobj=NULL, data=NULL, X, con = NULL, argvals = NULL, 
        method = c("OLS","GLS","mix"),
        gam.method = c("REML", "ML", "GCV.Cp", "GACV.Cp", "P-REML", "P-ML"), 
        cov.method = c("naive", "mod.chol"),
        lambda = NULL, nbasis=15, norder=4, 
        pen.order=2, multi.sp = ifelse(method=="OLS", FALSE, TRUE), pve=.99,
        max.iter = 1, maxlam = NULL, cv1 = FALSE, scale = FALSE)
{
    # parse formula
    if (!is.null(formula)) {
      if (is.null(data)) stop("Please specify the data.")
      tf <- terms.formula(formula)
      trmstrings <- attr(tf, "term.labels")
      terms <- sapply(trmstrings, function(trm) as.call(parse(text=trm))[[1]], simplify=FALSE)
      responsename <- as.character(attr(tf,"variables")[2][[1]])
      Y = data[,responsename]
      X = model.matrix(formula, data=data)
    }
  
  
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
        Bmat <- Theta <- eval.basis(argvals, bss)
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
        Bmat <- J12 <- svdJ$u %*% diag(sqrt(svdJ$d)) %*% t(svdJ$u)
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
            lofo <- lofocv(respmat, X.sc %x% Bmat, S1 = pen[[1]], argvals = argvals,
                    lamvec = lambda, constr = constr, maxlam = maxlam)
            cv = if (is.null(lambda))
                        lofo$objective
                    else min(lofo[, 2])
            lambda = if (is.null(lambda))
                        lofo$min
                    else lofo[which.min(lofo[, 2]), 1]
        }
    }

    firstfit <- amc(as.vector(t(respmat)), X.sc %x% Bmat,
            gam.method = gam.method, S = pen, C = constr, lambda = lambda)
    coefmat = coefmat.ols = t(matrix(firstfit$coef, ncol = q))
    se = NULL
    if (method != "OLS") {
        iter = 0
        coefmat.old = 3 * coefmat.ols
        newfit = NULL
        if (!is.null(lambda) & max.iter > 0)
            warning("Given lambda used for initial fit only")
        while (any(abs((coefmat - coefmat.old)/coefmat.old) > 0.001) & (iter < max.iter)) {
            iter = iter + 1
            if (max.iter > 1) cat("Refit", iter, "\n")
            oldfit = if (!is.null(newfit)) newfit else firstfit
            coefmat.old = coefmat
            residvec <- as.vector(t(respmat)) - (X.sc %x% Bmat) %*% oldfit$coef[1:(q*nbasis)]
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
                        #svdp = eigen(prec, symmetric=TRUE)
                        #sqrt.prec.list[[nband]] = svdp$vectors %*% tcrossprod(diag(sqrt(svdp$values)), svdp$vectors)
                        #F avoid unnecessary additional eigen-decomp
                        sqrt.prec.list[[nband]] = scale(t(TT), FALSE, sqrt(Ddiag))
                        lwprec = lw.test(residmat %*% sqrt.prec.list[[nband]])
                        lwstat[nband] = lwprec$stat; lwpval[nband] = lwprec$pvalue
                        if (lwstat[nband] < -5) break
                        #F insert abort to avoid running into numerical problems for thick bands with
                        #F  poor fit (i.e.: worse than 1-band and the previous iteration).
                        if(nband > 5){
                          if (lwstat[nband] > lwstat[1] && lwstat[nband] > lwstat[nband-1]) break
                        }
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
                            X.sc %x% (sqrt.prec %*% Bmat),
                            gam.method = gam.method, S = pen, C = constr, 
                            start = if (is.null(con)) as.vector(t(coefmat)) else NULL)
                    coefmat = t(matrix(newfit$coef, ncol = q))
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
                        efuncmat.scaled <- Bmat %*% t(t(pca.resid$harmonics$coef[ , 1:npc]) * sqrt(evalues))
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
                                c(as.vector(t(coefmat)), rep(0, ncurve*npc))
                            } else {
                                newfit$coefficients
                            }
                    #/F
                    #browser()
                    newfit <- amc(as.vector(t(respmat)),
                            cbind(X.sc %x% Bmat,diag(ncurve) %x% efuncmat.scaled),
                            gam.method = gam.method, S = pen.aug, C = constr.aug, 
                            start=if (is.null(constr.aug)) startB else NULL)      
                    

                    vecBt = newfit$coef[1:(q*nbasis)]
                    vecUt = newfit$coef[(q*nbasis+1):(q*nbasis+npc*ncurve)]
                    coefmat = t(matrix(vecBt, ncol=q))
                    U <- t(matrix(vecUt, ncol=ncurve))
                } #end mix
        } #end while 
    } # end !OLS
    if (method == "OLS" | max.iter == 0) {
        residvec <- as.vector(t(respmat)) - (X.sc %x% Bmat) %*% firstfit$coef
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
    fd = fd(t(coefmat), bss)
    est.func = eval.fd(argvals, fd)
    fit <- if (method == "mix" & max.iter > 0) newfit else firstfit
    roughness = diag(coefmat %*% getbasispenalty(bss, pen.order) %*% t(coefmat))
    skale = attr(X.sc, "scaled:scale")
    if (!is.null(skale)) {
        coefmat = t(scale(t(coefmat), center = FALSE, scale = skale))
        est.func = scale(est.func, center = FALSE, scale = skale)
        se.func = scale(se.func, center = FALSE, scale = skale)
        roughness = roughness/skale^2
    }
    yhat = if (resp.type=="raw") X %*% tcrossprod(coefmat, Theta) else fd(t(X %*% coefmat), bss)
    llist = list(fd = fd, pca.resid = pca.resid, U = U,
            yhat = yhat, resid = if (resp.type == 'raw') Y - yhat else fdobj - yhat,
            est.func = est.func,  se.func = se.func, argvals = argvals, fit = fit,
            edf = sum(fit$gam$edf),
            lambda = if (length(fit$gam$sp) > 0) fit$gam$sp else fit$gam$full.sp,
            cv = cv, roughness = roughness, resp.type = resp.type)
    class(llist) = "fosr"
    llist
}
