#' Construct an FLM regression term
#'
#' Defines a term \eqn{\int_{T}\beta(t)X_i(t)dt} for inclusion in an \code{\link[mgcv]{gam}}-formula
#' (or \code{\link{bam}} or \code{\link{gamm}} or \code{\link[gamm4]{gamm4}}) as constructed by
#' \code{\link{fgam}}, where \eqn{\beta(t)} is an unknown coefficient function and \eqn{X_i(t)}
#' is a functional predictor on the closed interval \eqn{T}. Defaults to a cubic B-spline with
#' second-order difference penalties for estimating \eqn{\beta(t)}.  The functional predictor must
#' be fully observed on a regular grid.
#' @param X an \code{N} by \code{J=ncol(xind)} matrix of function evaluations
#' \eqn{X_i(t_{i1}),., X_i(t_{iJ}); i=1,.,N.}
#' @param xind matrix (or vector) of indices of evaluations of \eqn{X_i(t)}; i.e. a matrix with
#' \emph{i}th row \eqn{(t_{i1},.,t_{iJ})}
#' @param integration method used for numerical integration. Defaults to \code{"simpson"}'s rule
#' for calculating entries in \code{L}. Alternatively and for non-equidistant grids,
#' \dQuote{\code{trapezoidal}} or \code{"riemann"}. \code{"riemann"} integration is always used if
#' \code{L} is specified
#' @param L an optional \code{N} by \code{ncol(xind)} matrix giving the weights for the numerical
#' integration over \code{t}
#' @param splinepars optional arguments specifying options for representing and penalizing the
#' functional coefficient \eqn{\beta(t)}. Defaults to a cubic B-spline with second-order difference
#' penalties, i.e. \code{list(bs="ps", m=c(2, 1))} See \code{\link{te}} or \code{\link{s}} for details
#' @param presmooth logical; if true, the functional predictor is pre-smoothed prior to fitting.  See
#' \code{\link{smooth.basisPar}}
#' @return a list with the following entries
#' \enumerate{
#' \item \code{call} - a \code{call} to \code{te} (or \code{s}, \code{t2}) using the appropriately
#' constructed covariate and weight matrices
#' \item \code{xind} - the \code{xind} argument supplied to \code{lf}
#' \item \code{L} - the  matrix of weights used for the integration
#' \item{xindname} - the name used for the functional predictor variable in the \code{formula}
#' used by \code{mgcv}
#' \item \code{tindname} - the name used for \code{xind} variable in the \code{formula} used by \code{mgcv}
#' \item \code{LXname} - the name used for the \code{L} variable in the \code{formula} used by \code{mgcv}
#' \item \code{presmooth} - the \code{presmooth} argument supplied to \code{lf}
#' \item \code{Xfd} - an \code{fd} object from presmoothing the functional predictors using
#' \code{\link{smooth.basisPar}}.  Only present if \code{presmooth=TRUE}.  See \code{\link{fd}}
#' }
#' @author Mathew W. McLean \email{mathew.w.mclean@@gmail.com} and Fabian Scheipl
#' @seealso \code{\link{fgam}}, \code{\link{af}}, mgcv's \code{\link{linear.functional.terms}},
#' \code{\link{fgam}} for examples
#' @importFrom fda create.bspline.basis smooth.basisPar eval.fd
#' @importFrom utils getFromNamesapce modifyList
#' @export
lf <- function(X, xind = seq(0, 1, l = ncol(X)),
               integration = c("simpson", "trapezoidal", "riemann"),
               L = NULL, splinepars = list(bs = "ps", k= min(ceiling(n/4),40),
                                  m = c(2, 2)), presmooth = TRUE) {

  n=nrow(X)
  nt=ncol(X)
  integration <- match.arg(integration)
  if(is.null(splinepars$bs)) splinepars$bs <- 'ps'
  if(is.null(splinepars$k)) splinepars$k <- min(ceiling(n/4),40)
  if(is.null(splinepars$m)) splinepars$m = c(2, 2)

  tindname <- paste(deparse(substitute(X)), ".tmat", sep = "")
  LXname <- paste("L.", deparse(substitute(X)), sep = "")
  basistype = "s"

  if (is.null(dim(xind))) {
    xind <- t(xind)
    stopifnot(ncol(xind) == nt)
    if (nrow(xind) == 1) {
      xind <- matrix(as.vector(xind), nrow = n, ncol = nt,
                     byrow = T)
    }
    stopifnot(nrow(xind) == n)
  }

  Xfd=NULL
  if(presmooth){
    bbt=create.bspline.basis(rangeval=range(xind),nbasis=ceiling(nt/4),
                             norder=splinepars$m[1]+2, breaks=NULL)

    # pre-smooth functional predictor
    temp <- smooth.basisPar(t(xind),t(X),bbt,int2Lfd(splinepars$m[2]))
    Xfd <- temp$fd
    Xfd$y2cMap <-temp$y2cMap
    X <- t(sapply(1:n,function(i){eval.fd(xind[i,],Xfd[i])}))
  }

  if (!is.null(L)) {
    stopifnot(nrow(L) == n, ncol(L) == nt)
  }else {
    L <- switch(integration, simpson = {
      ((xind[, nt] - xind[, 1])/nt)/3 * matrix(c(1,rep(c(4, 2), length = nt - 2), 1), nrow = n,
                                               ncol = nt, byrow = T)
    }, trapezoidal = {
      diffs <- t(apply(xind, 1, diff))
      0.5 * cbind(diffs[, 1], t(apply(diffs, 1, filter,filter = c(1, 1)))[, -(nt - 1)],
                  diffs[,(nt - 1)])
    }, riemann = {
      diffs <- t(apply(xind, 1, diff))
      cbind(rep(mean(diffs), n), diffs)
    })
  }
  LX <- L*X
  data <- list(xind, LX)
  names(data) <- c(tindname, LXname)
  splinefun <- as.symbol(basistype)
  frmls <- formals(getFromNamespace(deparse(splinefun), ns = "mgcv"))
  frmls <- modifyList(frmls[names(frmls) %in% names(splinepars)],
                      splinepars)
  call <- as.call(c(list(splinefun, x = as.symbol(substitute(tindname)),
                         by = as.symbol(substitute(LXname))),frmls))
  res <-list(call = call, data = data, xind = xind[1,], L = L, tindname=tindname,
             LXname=LXname,presmooth=presmooth)
  if(presmooth) res$Xfd <- Xfd
  return(res)
}

