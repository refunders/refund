##' Pointwise cross-validation for function-on-scalar regression
##'
##' Estimates prediction error for a function-on-scalar regression model by
##' leave-one-function-out cross-validation (CV), at each of a specified set of
##' points.
##'
##' Integrating the pointwise CV estimate over the function domain yields the
##' \emph{cross-validated integrated squared error}, the standard overall model
##' fit score returned by \code{\link{lofocv}}.
##'
##' It may be desirable to derive the value of \code{lambda} from an
##' appropriate call to \code{\link{fosr}}, as in the example below.
##'
##' @param fdobj a functional data object (class \code{fd}) giving the
##' functional responses.
##' @param Z the model matrix, whose columns represent scalar predictors.
##' @param L a row vector or matrix of linear contrasts of the coefficient
##' functions, to be restricted to equal zero.
##' @param lambda smoothing parameter: either a nonnegative scalar or a vector,
##' of length \code{ncol(Z)}, of nonnegative values.
##' @param eval.pts argument values at which the CV score is to be evaluated.
##' @param scale logical value or vector determining scaling of the matrix
##' \code{Z} (see \code{\link{scale}}, to which the value of this argument is
##' passed).
##' @return A vector of the same length as \code{eval.pts} giving the CV
##' scores.
##' @author Philip Reiss \email{phil.reiss@@nyumc.org}
##' @seealso \code{\link{fosr}}, \code{\link{lofocv}}
##' @references Reiss, P. T., Huang, L., and Mennes, M. (2010).  Fast
##' function-on-scalar regression with penalized basis expansions.
##' \emph{International Journal of Biostatistics}, 6(1), article 28.  Available
##' at \url{https://pubmed.ncbi.nlm.nih.gov/21969982/}
##' @export
##' @importFrom fda getbasispenalty eval.basis
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
