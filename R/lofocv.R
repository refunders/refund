##' Leave-one-function-out cross-validation
##'
##' This internal function, called by \code{fosr()} when \code{method="OLS"},
##' performs efficient leave-one-function-out cross-validation using
##' Demmler-Reinsch orthogonalization to choose the smoothing parameter.
##'
##'
##' @param Y matrix of responses, e.g. with columns corresponding to basis
##' function coefficients.
##' @param X model matrix.
##' @param S1 penalty matrix.
##' @param lamvec vector of candidate smoothing parameter values.  If
##' \code{NULL}, smoothing parameter is chosen by \code{\link{optimize}}.
##' @param constr matrix of linear constraints.
##' @param maxlam maximum smoothing parameter value to consider (when
##' \code{lamvec=NULL}).
##' @return if \code{lamvec=NULL}, a list (returned by \code{optimize}) with
##' elements \code{minimum} and \code{objective} giving, respectively, the
##' chosen smoothing parameter and the associated cross-validation score.
##' Otherwise a 2-column table with the candidate smoothing parameters in the
##' first column and the corresponding cross-validation scores in the second.
##' @author Philip Reiss \email{phil.reiss@@nyumc.org}
##' @seealso \code{\link{fosr}}, \code{\link{pwcv}}
##' @keywords internal
lofocv <-
function(Y, X, S1, lamvec=NULL, constr=NULL, maxlam=NULL) {
	nn = nrow(X)
	N = NROW(Y); K = NCOL(Y)
	if (N*K!=nn) stop('Number of elements of Y must equal number of rows of X')
	y = as.vector(t(Y))

	if (!is.null(constr)) {
	    # The following is based on Wood (2006), p. 186
	    n.con = dim(constr)[1]
	    Z. = qr.Q(qr(t(constr)), complete=TRUE)[ , -(1:n.con)]
	    X. = X %*% Z.
	    S1. = crossprod(Z., S1 %*% Z.)
	}
	else {
		X. = X
		S1. = S1
	}

	qrX = qr(X.)
	Rinv = solve(qr.R(qrX))
	svd211 = svd(crossprod(Rinv, S1. %*% Rinv))  # see p. 211 of Wood
	QU = qr.Q(qrX) %*% svd211$u

	cvfcn = function(lam) {
		A = tcrossprod(scale(QU, center=FALSE, scale=1+lam*svd211$d), QU)
		resmat = t(matrix(y - A %*% y, K))

		MSEp = 0
		for (i in 1:N) {
			ith = ((i-1)*K+1):(i*K)
			MSEp = MSEp + crossprod(solve(diag(K)-A[ith,ith], resmat[i, ])) / N
		}

	    MSEp
	}

	if (is.null(lamvec)) {  # minimize LOFO-CV criterion
		if (is.null(maxlam)) {  # use GCV-minimizing lambda
		    model.gcv = gam(y~X.-1, paraPen=list(X.=list(S1.)), method="GCV.Cp")
	        maxlam = model.gcv$sp
	    }

        cat("Finding optimal lambda by optimize()...\n")
	    opt = optimize(cvfcn, c(0, maxlam), tol=.01)
	    if (round(opt$minimum)==maxlam) warning("maxlam may be set too low")
	    return(opt)
	}

	else {  # calculate LOFO-CV for given values
		cvvals = c()
        cat("Calculating CV for candidate smoothing parameter values...\n")
		for (i in 1:length(lamvec)) cvvals[i] = cvfcn(lamvec[i])
		cvtable = cbind(lamvec, cvvals)
		dimnames(cvtable)[[2]] = c('lambda', 'LOFO-CV')
		print(cvtable)
		return(cvtable)
		if (which.min(cvvals)==1) warning("CV minimized at lowest lambda considered")
		if (which.min(cvvals)==length(lamvec)) warning("CV minimized at highest lambda considered")
	}
}

