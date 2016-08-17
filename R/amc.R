##' Additive model with constraints
##'
##' An internal function, called by \code{fosr()}, that fits additive models
##' with linear constraints via a call to \code{\link[mgcv]{gam}} or
##' \code{\link[mgcv]{bam}} in the \pkg{mgcv} package.
##'
##' The additive model is fitted using \code{\link[mgcv]{gam}}, unless there
##' are more than 10000 responses; in that case \code{\link[mgcv]{bam}} is
##' used.
##'
##' @param y response vector.
##' @param Xmat design matrix.
##' @param S list of penalty matrices.
##' @param gam.method smoothing parameter selection method: "REML" for
##' restricted maximum likelihood, "GCV.Cp" for generalized cross-validation.
##' @param C matrix of linear constraints.  Dimension should be number of
##' constraints times \code{ncol(Xmat)}.
##' @param lambda smoothing parameter value.  If \code{NULL}, the smoothing
##' parameter(s) will be estimated.
##' @param \dots other arguments, passed to \code{\link[mgcv]{gam}} or
##' \code{\link[mgcv]{bam}}.
##' @return A list with the following elements: \item{gam}{the \code{gam}
##' object returned by \code{gam} or \code{bam}.}
##' \item{coefficients}{coefficients with respect to design matrix \code{Xmat},
##' derived from the \code{gam()} fit.} \item{Vp, GinvXt}{outputs used by
##' \code{fosr}.} \item{method}{the \code{gam.method} argument of the call to
##' \code{amc}.}
##' @author Philip Reiss \email{phil.reiss@@nyumc.org}
##' @seealso \code{\link{fosr}}
##' @keywords internal
##' @importFrom mgcv bam gam
amc <- function(y, Xmat, S, gam.method='REML', C=NULL, lambda=NULL, ...) {
  n.p = length(S)
  stopifnot( is.null(lambda) | length(lambda)==n.p )
  if (!is.null(C)) {
    # The following is based on Wood (2006), p. 186
      n.con = dim(C)[1]
      Z. = qr.Q(qr(t(C)), complete=TRUE)[ , -(1:n.con)]
      Xmat. = Xmat %*% Z.
      S. = vector("list", n.p)
      #F: for (i in 1:n.p) S.[[i]] = crossprod(Z., S[[i]] %*% Z.)
      #F: seems like mgcv >= 1.8.13
      for (i in 1:n.p) {
        if(is.null(lambda)) {
          S.[[i]] = list(crossprod(Z., S[[i]] %*% Z.))
        } else {
          S.[[i]] = list(crossprod(Z., S[[i]] %*% Z.), sp = lambda[i])
        }
      }
  }
  else {
    Z. = diag(ncol(Xmat))
    Xmat. = Xmat
    if(is.null(lambda)) {
      S. = list(list(S[[1]]))
    } else {
      S. = list(list(S[[1]], sp=lambda[1]))
  }
  }

    fitter = if (length(y) > 10000) bam else gam
    #if (is.null(lambda)) fitobj = fitter(y ~ Xmat.-1, method=gam.method, paraPen=list(Xmat.=S.), ...)
    #else
    #fitobj = fitter(y ~ Xmat.-1, paraPen=list(Xmat.=S.), sp=lambda, ...)
    fitobj = fitter(y ~ Xmat.-1, paraPen=list(Xmat.=S.[[1]]), ...)

  lambdavec = if (!is.null(fitobj$full.sp)) fitobj$full.sp else fitobj$sp
  fullpen = 0
  for (i in 1:n.p) fullpen = lambdavec[i] * S.[[i]][[1]]
  list(gam = fitobj,
       coefficients = Z. %*% fitobj$coef,
       Vp = Z. %*% fitobj$Vp %*% t(Z.),
       GinvXT = Z. %*% solve(crossprod(Xmat.) + fullpen, t(Xmat.)),
       method = gam.method)
}

