#' Construct an FGAM regression term
#'
#' Defines a term \eqn{\int_{T}F(X_i(t),t)dt} for inclusion in an \code{mgcv::gam}-formula (or
#' \code{\link{bam}} or \code{\link{gamm}} or \code{gamm4:::gamm}) as constructed by
#' \code{\link{fgam}}, where \eqn{F(x,t)}$ is an unknown smooth bivariate function and \eqn{X_i(t)}
#' is a functional predictor on the closed interval \eqn{T}. Defaults to a cubic tensor product
#' B-spline with marginal second-order difference penalties for estimating \eqn{F(x,t)}.  The
#' functional predictor must be fully observed on a regular grid
#' @param X an \code{N} by \code{J=ncol(xind)} matrix of function evaluations
#' \eqn{X_i(t_{i1}),., X_i(t_{iJ}); i=1,.,N.}
#' @param xind matrix (or vector) of indices of evaluations of \eqn{X_i(t)}; i.e. a matrix with
#' \emph{i}th row \eqn{(t_{i1},.,t_{iJ})}
#' @param basistype defaults to \code{"te"}, i.e. a tensor product spline to represent \eqn{F(x,t)} Alternatively,
#' use \code{"s"} for bivariate basis functions (see \code{\link{s}}) or \code{"t2"} for an alternative
#' parameterization of tensor product splines (see \code{\link{t2}})
#' @param integration method used for numerical integration. Defaults to \code{"simpson"}'s rule for
#' calculating entries in \code{L}. Alternatively and for non-equidistant grids, \code{"trapezoidal"}
#' or \code{"riemann"}. \code{"riemann"} integration is always used if \code{L} is specified
#' @param splinepars optional arguments specifying options for representing and penalizing the
#' function \eqn{F(x,t)}. Defaults to a cubic tensor product B-spline with marginal second-order
#' difference penalties, i.e. \code{list(bs="ps", m=list(c(2, 2), c(2, 2))}, see \code{\link{te}} or
#' \code{\link{s}} for details
#' @param presmooth logical; if true, the functional predictor is pre-smoothed prior to fitting; see
#' \code{\link{smooth.basisPar}}
#' @param Xrange numeric; range to use when specifying the marginal basis for the \emph{x}-axis.  It may
#' be desired to increase this slightly over the default of \code{range(X)} if concerned about predicting
#' for future observed curves that take values outside of \code{range(X)}
#' @param Qtransform logical; should the functional be transformed using the empirical cdf and
#' applying a quantile transformation on each column of \code{X} prior to fitting?  This ensures
#' \code{Xrange=c(0,1)}.  If \code{Qtransform=TRUE} and \code{presmooth=TRUE}, presmoothing is done prior
#' to transforming the functional predictor
#' @param L optional weight matrix for the linear functional
#' @return A list with the following entries:
#' \enumerate{
#' \item \code{call} - a \code{"call"} to \code{te} (or \code{s}, \code{t2}) using the appropriately
#' constructed covariate and weight matrices.
#' \item \code{xind} - the \code{xind} argument supplied to \code{af}
#' \item \code{L}{ the  matrix of weights used for the integration
#' \item \code{xindname}{ the name used for the functional predictor variable in the \code{formula} used by \code{mgcv}.}
#' \item \code{tindname} - the name used for \code{xind} variable in the \code{formula} used by \code{mgcv}
#' \item \code{Lname} - the name used for the \code{L} variable in the \code{formula} used by \code{mgcv}
#' \item \code{presmooth} - the \code{presmooth} argument supplied to \code{af}
#' \item \code{Qtranform} - the \code{Qtransform} argument supplied to \code{af}
#' \item \code{Xrange} - the \code{Xrange} argument supplied to \code{af}
#' \item \code{ecdflist} - a list containing one empirical cdf function from applying \code{\link{ecdf}}
#' to each (possibly presmoothed) column of \code{X}.  Only present if \code{Qtransform=TRUE}
#' \item \code{Xfd} - an \code{fd} object from presmoothing the functional predictors using
#' \code{\link{smooth.basisPar}}.  Only present if \code{presmooth=TRUE}.  See \code{\link{fd}}.}
#' }
#' @author Mathew W. McLean \email{mathew.w.mclean@@gmail.com} and Fabian Scheipl
#' @references McLean, M. W., Hooker, G., Staicu, A.-M., Scheipl, F., and Ruppert, D. (2014). Functional
#' generalized additive models. \emph{Journal of Computational and Graphical Statistics}, \bold{23 (1)},
#' pp. 249-269.  Available at \url{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3982924}.
#' @seealso \code{\link{fgam}}, \code{\link{lf}}, mgcv's \code{\link{linear.functional.terms}},
#' \code{\link{fgam}} for examples
#' @importFrom stats ecdf
#' @importFrom fda int2Lfd smooth.basisPar eval.fd create.bspline.basis
#' @importFrom utils modifyList getFromNamespace
af <- function(X, xind = seq(0, 1, l = ncol(X)), basistype = c("te","t2", "s"),
               integration = c("simpson", "trapezoidal", "riemann"),
               L = NULL, 
               splinepars = list(bs = "ps",
                                 k= c(min(ceiling(nrow(X)/5),20),min(ceiling(ncol(X)/5),20)),
                                 m = list(c(2, 2), c(2, 2))), presmooth = TRUE,Xrange=range(X),Qtransform=FALSE) {
  
  n=nrow(X)
  nt=ncol(X)
  basistype <- match.arg(basistype)
  integration <- match.arg(integration)
  if(is.null(splinepars$bs)) splinepars$bs <- 'ps'
  if(is.null(splinepars$k)) splinepars$k <- c(min(ceiling(nrow(X)/5),20),
                                                min(ceiling(ncol(X)/5),20))
  if(is.null(splinepars$m)) splinepars$m = list(c(2, 2), c(2, 2))

  xindname <- paste(deparse(substitute(X)), ".omat", sep = "")
  tindname <- paste(deparse(substitute(X)), ".tmat", sep = "")
  Lname <- paste("L.", deparse(substitute(X)), sep = "")

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
                             norder=splinepars$m[[2]][1]+2, breaks=NULL)

    # pre-smooth functional predictor
    temp <- smooth.basisPar(t(xind),t(X),bbt,int2Lfd(splinepars$m[[2]][2]))
    Xfd <- temp$fd
    Xfd$y2cMap <-temp$y2cMap
    X <- t(sapply(1:n,function(i){eval.fd(xind[i,],Xfd[i])}))

    # need to check that smoothing didn't change range of data
    if(!Qtransform){
      if(max(X)>Xrange[2]){
        Xrange[2] <- max(X)
      }
      if(min(X)<Xrange[1]){
        Xrange[1] <- min(X)
      }
    }
  }

  ecdf=NULL
  if(Qtransform){
    Xrange <- c(0,1)
    X <- apply(X, 2, function(x){ (rank(x)-1)/(length(x)-1)} )
    # need to keep ecdf's for prediction later
    ecdflist <- apply(X, 2, ecdf)
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

  data <- list(X, xind, L)
  names(data) <- c(xindname, tindname, Lname)
  splinefun <- as.symbol(basistype)
  frmls <- formals(getFromNamespace(deparse(splinefun), ns = "mgcv"))
  frmls <- modifyList(frmls[names(frmls) %in% names(splinepars)],
                      splinepars)
  call <- as.call(c(list(splinefun, x = as.symbol(substitute(xindname)),
                         z = as.symbol(substitute(tindname)), by = as.symbol(substitute(Lname))),
                    frmls))
  res <-list(call = call, data = data, xind = xind[1,], L = L, xindname = xindname, tindname=tindname,
             Lname=Lname,Qtransform=Qtransform,presmooth=presmooth,Xrange=Xrange)
  if(Qtransform) res$ecdflist <- ecdflist
  if(presmooth) res$Xfd <- Xfd
  return(res)
}


