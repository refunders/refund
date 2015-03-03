##' Functional principal component regression
##'
##' Implements functional principal component regression (Reiss and Ogden,
##' 2007, 2010) for generalized linear models with scalar responses and
##' functional predictors.
##'
##' One-dimensional functional predictors can be given either in functional
##' data object form, using argument \code{fdobj} (see the \pkg{fda} package of
##' Ramsay, Hooker and Graves, 2009, and Method 1 in the example below), or
##' explicitly, using \code{xfuncs} (see Method 2 in the example).  In the
##' latter case, arguments \code{basismat} and \code{penmat} can also be used
##' to specify the basis and/or penalty matrices (see Method 3).
##'
##' For two-dimensional predictors, functional data object form is not
##' supported.  Instead of radial B-splines as in Reiss and Ogden (2010), this
##' implementation employs tensor product cubic B-splines, sometimes known as
##' bivariate O-splines (Ormerod, 2008).
##'
##' For purposes of interpreting the fitted coefficients, please note that the
##' functional predictors are decorrelated from the scalar predictors before
##' fitting the model (when there are no scalar predictors other than the
##' intercept, this just means the columns of the functional predictor matrix
##' are de-meaned); see section 3.2 of Reiss (2006) for details.
##'
##' @param y scalar outcome vector.
##' @param xfuncs for 1D predictors, an \eqn{n \times d} matrix of
##' signals/functional predictors, where \eqn{n} is the length of \code{y} and
##' \eqn{d} is the number of sites at which each signal is defined.  For 2D
##' predictors, an \eqn{n \times d1 \times d2} array representing \eqn{n}
##' images of dimension \eqn{d1 \times d2}.
##' @param fdobj functional data object (class "\code{\link[fda]{fd}}") giving
##' the functional predictors. Allowed only for 1D functional predictors.
##' @param ncomp number of principal components. If \code{NULL}, this is chosen
##' by \code{pve}.
##' @param pve proportion of variance explained: used to choose the number of
##' principal components.
##' @param nbasis number(s) of B-spline basis functions: either a scalar, or a
##' vector of values to be compared.  For 2D predictors, tensor product
##' B-splines are used, with the given basis dimension(s) in each direction;
##' alternatively, \code{nbasis} can be given in the form \code{list(v1,v2)},
##' in which case cross-validation will be performed for each combination of
##' the first-dimension basis sizes in \code{v1} and the second-dimension basis
##' sizes in \code{v2}. Ignored if \code{fdobj} is supplied. If \code{fdobj} is
##' \emph{not} supplied, this defaults to 40 (i.e., 40 B-spline basis
##' functions) for 1D predictors, and 15 (i.e., tensor product B-splines with
##' 15 functions per dimension) for 2D predictors.
##' @param basismat a \eqn{d \times K} matrix of values of \eqn{K} basis
##' functions at the \eqn{d} sites.
##' @param penmat a \eqn{K \times K} matrix defining a penalty on the basis
##' coefficients.
##' @param argvals points at which the functional predictors and the
##' coefficient function are evaluated.  By default, if 1D functional
##' predictors are given by the \eqn{n \times d} matrix \code{xfuncs},
##' \code{argvals} is set to \eqn{d} equally spaced points from 0 to 1; if they
##' are given by \code{fdobj}, \code{argvals} is set to 401 equally spaced
##' points spanning the domain of the given functions. For 2D (image)
##' predictors supplied as an \eqn{n \times d1 \times d2} array, \code{argvals}
##' defaults to a list of (1) \eqn{d1} equally spaced points from 0 to 1; (2)
##' \eqn{d2} equally spaced points from 0 to 1.
##' @param covt covariates: an \eqn{n}-row matrix, or a vector of length
##' \eqn{n}.
##' @param mean.signal.term logical: should the mean of each subject's signal
##' be included as a covariate?
##' @param spline.order order of B-splines used, if \code{fdobj} is not
##' supplied; defaults to \code{4}, i.e., cubic B-splines.
##' @param family generalized linear model family. Current version supports
##' \code{"gaussian"} (the default) and \code{"binomial"}.
##' @param method smoothing parameter selection method, passed to function
##' \code{\link[mgcv]{gam}}; see the \code{\link[mgcv]{gam}} documentation for
##' details.
##' @param sp a fixed smoothing parameter; if \code{NULL}, an optimal value is
##' chosen (see \code{method}).
##' @param pen.order order of derivative penalty applied when estimating the
##' coefficient function; defaults to \code{2}.
##' @param cv1 logical: should cross-validation be performed to select the best
##' model if only one set of tuning parameters provided? By default,
##' \code{FALSE}. Note that, if there are multiple sets of tuning parameters
##' provided, cv is always performed.
##' @param nfold the number of validation sets ("folds") into which the data
##' are divided; by default, 5.
##' @param store.cv logical: should a CV result table be in the output? By
##' default, \code{FALSE}.
##' @param store.gam logical: should the \code{\link[mgcv]{gam}} object be
##' included in the output? Defaults to \code{TRUE}.
##' @param \dots other arguments passed to function \code{\link[mgcv]{gam}}.
##' @return A list with components \item{gamObject}{if \code{store.gam = TRUE},
##' an object of class \code{gam} (see \code{\link[mgcv]{gamObject}} in the
##' \pkg{mgcv} package documentation).} \item{fhat}{coefficient function
##' estimate.} \item{se}{pointwise Bayesian standard error.}
##' \item{undecor.coef}{undecorrelated coefficient for covariates.}
##' \item{argvals}{the supplied value of \code{argvals}.} \item{cv.table}{a
##' table giving the CV criterion for each combination of \code{nbasis} and
##' \code{ncomp}, if \code{store.cv = TRUE}; otherwise, the CV criterion only
##' for the optimized combination of these parameters.  Set to \code{NULL} if
##' CV is not performed.} \item{nbasis, ncomp}{when CV is performed, the values
##' of \code{nbasis} and \code{comp} that minimize the CV criterion.}
##' @author Philip Reiss \email{phil.reiss@@nyumc.org}, Lan Huo
##' \email{lan.huo@@nyumc.org}, and Lei Huang \email{huangracer@@gmail.com}
##' @references Ormerod, J. T. (2008).  On semiparametric regression and data
##' mining.  Ph.D. dissertation, School of Mathematics and Statistics,
##' University of New South Wales.
##'
##' Ramsay, J. O., Hooker, G., and Graves, S. (2009). \emph{Functional Data
##' Analysis with R and MATLAB}. New York: Springer.
##'
##' Reiss, P. T. (2006).  Regression with signals and images as predictors.
##' Ph.D. dissertation, Department of Biostatistics, Columbia University.
##' Available at \url{http://works.bepress.com/phil_reiss/11/}.
##'
##' Reiss, P. T., and Ogden, R. T. (2007).  Functional principal component
##' regression and functional partial least squares.  \emph{Journal of the
##' American Statistical Association}, 102, 984--996.
##'
##' Reiss, P. T., and Ogden, R. T. (2010).  Functional generalized linear
##' models with images as predictors.  \emph{Biometrics}, 66, 61--69.
##'
##' Wood, S. N. (2006). \emph{Generalized Additive Models: An Introduction with
##' R}. Boca Raton, FL: Chapman & Hall.
##' @examples
##'
##' require(fda)
##' ### 1D functional predictor example ###
##'
##' ######### Octane data example #########
##' data(gasoline)
##'
##' # Create the requisite functional data objects
##' bbasis = create.bspline.basis(c(900, 1700), 40)
##' wavelengths = 2*450:850
##' nir = matrix(NA, 401, 60)
##' for (i in 1:60) nir[, i] = gasoline$NIR[i, ]
##' # Why not just take transpose of gasoline$NIR above?
##' # Because for some reason it leads to an error in the following statement
##' gas.fd = smooth.basisPar(wavelengths, nir, bbasis)$fd
##'
##' # Method 1: Call fpcr with fdobj argument
##' gasmod1 = fpcr(gasoline$octane, fdobj = gas.fd, ncomp = 30)
##' plot(gasmod1, xlab="Wavelength")
##' \dontrun{
##' # Method 2: Call fpcr with explicit signal matrix
##' gasmod2 = fpcr(gasoline$octane, xfuncs = gasoline$NIR, ncomp = 30)
##' # Method 3: Call fpcr with explicit signal, basis, and penalty matrices
##' gasmod3 = fpcr(gasoline$octane, xfuncs = gasoline$NIR,
##'                basismat = eval.basis(wavelengths, bbasis),
##'                penmat = getbasispenalty(bbasis), ncomp = 30)
##'
##' # Check that all 3 calls yield essentially identical estimates
##' all.equal(gasmod1$fhat, gasmod2$fhat, gasmod3$fhat)
##' # But note that, in general, you'd have to specify argvals in Method 1
##' # to get the same coefficient function values as with Methods 2 & 3.
##' }
##'
##' ### 2D functional predictor example ###
##'
##' n = 200; d = 70
##'
##' # Create true coefficient function
##' ftrue = matrix(0,d,d)
##' ftrue[40:46,34:38] = 1
##'
##' # Generate random functional predictors, and scalar responses
##' ii = array(rnorm(n*d^2), dim=c(n,d,d))
##' iimat = ii; dim(iimat) = c(n,d^2)
##' yy = iimat %*% as.vector(ftrue) + rnorm(n, sd=.3)
##'
##' mm = fpcr(yy, ii, ncomp=40)
##'
##' image(ftrue)
##' contour(mm$fhat, add=TRUE)
##'
##' \dontrun{
##' ### Cross-validation ###
##' cv.gas = fpcr(gasoline$octane, xfuncs = gasoline$NIR,
##'                  nbasis=seq(20,40,5), ncomp = seq(10,20,5), store.cv = TRUE)
##' image(seq(20,40,5), seq(10,20,5), cv.gas$cv.table, xlab="Basis size",
##'       ylab="Number of PCs", xaxp=c(20,40,4), yaxp=c(10,20,2))
##'
##' }
##' @export
##' @importFrom mgcv gam
##' @importFrom MASS ginv
fpcr <- function(y, xfuncs = NULL, fdobj = NULL, ncomp=NULL, pve = 0.99, nbasis = NULL, basismat = NULL, 
                 penmat = NULL, argvals = NULL, covt = NULL, mean.signal.term = FALSE, 
                 spline.order = NULL, family = "gaussian", method = "REML", sp = NULL, 
                 pen.order = 2, cv1 = FALSE, nfold = 5, store.cv = FALSE, store.gam = TRUE, ...){
    # Error handling
    n <- length(y)
    do.cv <- FALSE
    if (!nfold %in% 1:n) 
        stop("Argument 'nfold' is invalid: must be an integer between 1 and ", n, ".") 
    if (!family %in% c("gaussian", "binomial")) 
        stop("Only 'gaussian' and 'binomial' models are implemented in the current version.")
    if (is.null(fdobj)) {
        if (!is.array(xfuncs) || !length(dim(xfuncs)) %in% 2:3)
    	    stop("xfuncs must either be a 2D or 3D array")  
        dim.sig <- length(dim(xfuncs)) - 1
        if (is.null(nbasis)) nbasis <- ifelse(dim.sig == 1, 40, 15)
        if (dim.sig == 1){
    	    if (!is.list(nbasis) && is.vector(nbasis))
    		    nbs <- matrix(nbasis, ncol = 1)
    	    else stop("for 1D predictors, 'nbasis' must be a vector")   
    	    siglength <- ncol(xfuncs)
        } else{
    	    if (is.list(nbasis) && length(nbasis) == 2)
    		    nbs <- cbind(rep(nbasis[[1]], length(nbasis[[2]])), rep(nbasis[[2]], 
    		                 each = length(nbasis[[1]])))
    	    else if (!is.list(nbasis) && is.vector(nbasis))
    		    nbs <- matrix(rep(nbasis,2), ncol=2)
    	    else stop("for 2D predictors, 'nbasis' must either be a vector or a list of length 2")
    	    d1 <- dim(xfuncs)[2L]
            d2 <- dim(xfuncs)[3L]
            siglength <- d1 * d2
        }
        if (cv1 || nrow(nbs) > 2 || (!is.null(ncomp) && length(ncomp) > 1)) do.cv <- TRUE
    } 
    else if (cv1 || (!is.null(ncomp) && length(ncomp) > 1)) do.cv <- TRUE 
    if (!do.cv) {
        store.cv <- FALSE
        nfold <- 1
    }
    groups <- split(sample(1:n), rep(1:nfold, length = n))   
    n.unpen.cols <- 1 + mean.signal.term + ifelse(is.null(covt), 0, ncol(as.matrix(covt)))
    cv.table <- array(0, dim = c(ifelse(is.null(fdobj),nrow(nbs), 1), ifelse(is.null(ncomp), 1, length(ncomp))))
    for (inb in 1 : dim(cv.table)[1]){
    	st <- fpcr.setup(y = y, xfuncs = xfuncs, fdobj = fdobj, nbasis = if (is.null(fdobj)) nbs[inb,] else NULL, 
    	                 basismat = basismat, penmat = penmat, argvals = argvals, covt = covt, 
    	                 mean.signal.term = mean.signal.term, spline.order = spline.order, 
    	                 pen.order = pen.order)
    	argvals <- st$argvals
    	if(is.null(fdobj)) cat("nbasis:", st$nbasis, "\n")
    	for (ifold in 1 : nfold){
    		if (do.cv) {
    		    idxTest <- groups[[ifold]]
    		    idxTrain <- (1:n)[-idxTest]
    		}
    		else idxTest <- idxTrain <- 1:n
    		X0.tr <- st$X0[idxTrain,]
    		SB.tr <- st$SB[idxTrain,]
    		svdSB <- svd(SB.tr)
    		if (is.null(ncomp)) ncomp <- min(which(cumsum(svdSB$d) > pve * sum(svdSB$d)))
    		for (incomp in 1 : length(ncomp)) {
    			V.ncomp <- svdSB$v[,1:ncomp[incomp]]
    			X <- cbind(X0.tr, SB.tr %*% V.ncomp)
    			S <- list(matrix(0, ncol(X), ncol(X)))
    			S[[1]][-(1:n.unpen.cols), -(1:n.unpen.cols)] <- 
    			                                 crossprod(V.ncomp, st$penmat %*% V.ncomp)
    			obje <- gam(y[idxTrain]~X - 1, paraPen = list(X=S), family = get(family), 
    			            method = method, sp = sp, ...)
    			BV <- st$basismat %*% V.ncomp
    			fhat <- BV %*% obje$coef[-(1:n.unpen.cols)]
    			undecor.coef <- obje$coef[1:n.unpen.cols] - 
    			                ginv(X0.tr) %*% st$xfuncs[idxTrain,] %*% fhat
    			yhat <- st$X0[idxTest,] %*% undecor.coef + st$xfuncs[idxTest,] %*% fhat
    			if (family == "gaussian")
    			    cv.table[inb, incomp] <- cv.table[inb, incomp] + 
    			                             mean((yhat - y[idxTest]) ^ 2)
    			else if (family == "binomial"){
                    phat <- exp(yhat) / (1 + exp(yhat))
   					phat <- replace(phat, exp(yhat) == Inf, 1)
   					cv.table[inb, incomp] <- cv.table[inb, incomp] + 
   					                         mean((phat > mean(y[idxTrain])) != y[idxTest]) 
   				}
    		}
    	}                 
    }
    if (do.cv) {
        idxmin <- which(cv.table == min(cv.table[cv.table!=0], na.rm = TRUE), arr.ind = TRUE)
        if (nrow(idxmin) > 1) idxmin <- idxmin[1,]
   	    if (is.list(nbasis)){
    	    dim(cv.table) <- c(length(nbasis[[1]]), length(nbasis[[2]]), length(ncomp))
            dimnames(cv.table) <- list(paste("nbasis1=", nbasis[[1]]), 
                                       paste("nbasis2=", nbasis[[2]]), paste("ncomp", ncomp))
        } 
        else dimnames(cv.table) <- list(paste("nbasis=", nbasis), paste("ncomp", ncomp))   
        if (dim.sig == 1) nbasis <- nbs[idxmin[1],]
        else {
        	nbasis <- list()
        	nbasis[[1]] <- nbs[idxmin[1],1]
        	nbasis[[2]] <- nbs[idxmin[1],2]
        }
        obje <- fpcr(y = y, xfuncs = xfuncs, fdobj = fdobj, ncomp = ncomp[idxmin[2]], 
                     nbasis = nbasis, basismat = basismat, penmat = penmat, 
                     argvals = argvals, covt = covt, mean.signal.term = mean.signal.term, 
                     spline.order = spline.order, family = family, method = method, sp = sp, 
                     pen.order = pen.order, cv1 = FALSE, store.gam = store.gam, ...)
        ncomp <- ncomp[idxmin[2]]
        if (store.cv) obje$cv.table <- cv.table
        else obje$cv.table <- min(cv.table[cv.table!=0])
    } 
    else {
        se <- sqrt(rowSums((BV %*% obje$Vp[-(1:n.unpen.cols), -(1:n.unpen.cols)]) * BV))
        if (!store.gam) {
        	yhat <- obje$fitted.values
        	obje <- list()
        	obje$fitted.values <- yhat
        }
        obje$se <- se
        obje$argvals <- argvals
    	obje$undecor.coef <- as.matrix(undecor.coef, nrow = 1)
    	colnames(obje$undecor.coef) <- if (is.null(dimnames(covt)) || is.null(dimnames(covt)[[2]])) 
                                          paste("X", 0:(n.unpen.cols-1), sep="")
    	                               else c("Intercept", dimnames(covt)[[2]])
        if (st$dim.sig == 2) dim(fhat) <- c(d1, d2)
        obje$fhat <- fhat
        obje$nbasis <- nbasis
        obje$ncomp <- ncomp
    }
    class(obje)="fpcr"
    return(obje)    
}
