#' Construct a smooth function-on-function regression term
#'
#' Defines a term \eqn{\int^{s_{hi, i}}_{s_{lo, i}} f(X_i(s), s, t) ds} for
#' inclusion in an \code{mgcv::gam}-formula (or \code{bam} or \code{gamm} or
#' \code{gamm4:::gamm}) as constructed by \code{\link{pffr}}. Defaults to a
#' cubic tensor product B-spline with marginal second differences penalties for
#' \eqn{f(X_i(s), s, t)} and integration over the entire range \eqn{[s_{lo, i},
#' s_{hi, i}] = [\min(s_i), \max(s_i)]}. Can't deal with any missing \eqn{X(s)},
#' unequal lengths of \eqn{X_i(s)} not (yet?) possible. Unequal ranges for
#' different \eqn{X_i(s)} should work. \eqn{X_i(s)} is assumed to be numeric.\cr
#' \code{sff()} IS AN EXPERIMENTAL FEATURE AND NOT WELL TESTED YET -- USE AT
#' YOUR OWN RISK.
#'
#' @param X an n by \code{ncol(xind)} matrix of function evaluations
#'   \eqn{X_i(s_{i1}),\dots, X_i(s_{iS})}; \eqn{i=1,\dots,n}.
#' @param yind \emph{DEPRECATED} matrix (or vector) of indices of evaluations of
#'   \eqn{Y_i(t)}; i.e. matrix with rows \eqn{(t_{i1},\dots,t_{iT})}; no longer
#'   used.
#' @param xind vector of indices of evaluations of \eqn{X_i(s)},
#'   i.e, \eqn{(s_{1},\dots,s_{S})}
#' @param basistype defaults to "\code{\link[mgcv]{te}}", i.e. a tensor product
#'   spline to represent \eqn{f(X_i(s), t)}. Alternatively, use \code{"s"} for
#'   bivariate basis functions (see \code{\link[mgcv]{s}}) or \code{"t2"} for an
#'   alternative parameterization of tensor product splines (see
#'   \code{\link[mgcv]{t2}}).
#' @param integration method used for numerical integration. Defaults to
#'   \code{"simpson"}'s rule. Alternatively and for non-equidistant grids,
#'   \code{"trapezoidal"}.
#' @param L optional: an n by \code{ncol(xind)} giving the weights for the
#'   numerical integration over \eqn{s}.
#' @param limits defaults to NULL for integration across the entire range of
#'   \eqn{X(s)}, otherwise specifies the integration limits \eqn{s_{hi, i},
#'   s_{lo, i}}: either one of \code{"s<t"} or \code{"s<=t"} for \eqn{(s_{hi,
#'   i}, s_{lo, i}) = (0, t)} or a function that takes \code{s} as the first and
#'   \code{t} as the second argument and returns TRUE for combinations of values
#'   \code{(s,t)} if \code{s} falls into the integration range for the given
#'   \code{t}. This is an experimental feature and not well tested yet; use at
#'   your own risk.
#' @param splinepars optional arguments supplied to the \code{basistype}-term.
#'   Defaults to a cubic tensor product B-spline with marginal second
#'   differences, i.e. \code{list(bs="ps", m=list(c(2,2),c(2,2),c(2,2))}. See
#'   \code{\link[mgcv]{te}} or \code{\link[mgcv]{s}} for details
#' @param constant_over_yind fit a term \eqn{\int^{s_{hi, i}}_{s_{lo, i}} f(X_i(s), s) ds}
#' that is constant over t? Defaults to `FALSE`. Integration limits can still depend on the
#' index of the response and induce estimated effects that vary over t.
#' @param use_lag if `limits` defines a window of effectiveness with constant width h around t,
#'  it can make sense to fit effects like  \eqn{\int^{h}_{0} f(X_i(t-h), h) dh}, i.e.,
#'  only the time distance ("lag") between exposure and response affects their partial associations, not the
#'  timescales themselves.
#' @return a list containing \itemize{ \item \code{call} a "call" to
#'   \code{\link[mgcv]{te}} (or \code{\link[mgcv]{s}}, \code{\link[mgcv]{t2}})
#'   using the appropriately constructed covariate and weight matrices (see
#'   \code{\link[mgcv]{linear.functional.terms}}) \item \code{data} a list
#'   containing the necessary covariate and weight matrices }
#'
#' @author Fabian Scheipl, based on Sonja Greven's trick for fitting functional
#'   responses.
#' @export
#' @importFrom utils getFromNamespace modifyList
# FIXME: figure out weights for simpson's rule on non-equidistant grids
# TODO: allow X to be of class fd (?)
# TODO: by variables (?)
sff <- function(X,
                yind,
                xind=seq(0, 1, l=ncol(X)),
                basistype= c("te", "t2", "s"),
                integration=c("simpson", "trapezoidal"),
                L=NULL,
                limits=NULL,
                constant_over_yind = FALSE,
                use_lag = FALSE,
                splinepars=list(bs="ps", m=list(c(2,2),c(2,2),c(2,2)))
){
  n <- nrow(X)
  nxgrid <- ncol(X)
  stopifnot(all(!is.na(X)))

  # check & format index for X
  if(is.null(dim(xind))){
    xind <- t(as.matrix(xind))
  }
  stopifnot(ncol(xind) == nxgrid)
  if(nrow(xind)== 1){
    xind <- matrix(as.vector(xind), nrow=n, ncol=nxgrid, byrow=T)
  } else {
    stop("<xind> has to be supplied as a vector or matrix with a single row.")
  }
  stopifnot(nrow(xind) == n)
  stopifnot(all.equal(order(xind[1,]), 1:nxgrid))

  basistype <- match.arg(basistype)
  integration <- match.arg(integration)

  # scale xind to [0, 1] and check for reasonably equidistant gridpoints
  xind.sc <- xind - min(xind)
  xind.sc <- xind.sc/max(xind.sc)
  diffXind <- t(round(apply(xind.sc, 1, diff), 3))
  if(is.null(L) & any(apply(diffXind, 1, function(x) length(unique(x))) != 1) && # gridpoints for any  X_i(s) not equidistant?
       integration=="simpson"){
    warning("Non-equidistant grid detected for ", deparse(substitute(X)), ".\n Changing to trapezoidal rule for integration.")
    integration <- "trapezoidal"

  }
  # FIXME: figure out weights for simpson's rule on non-equidistant grids instead of all this...

  #make weight matrix for by-term
  if(!is.null(L)){
    stopifnot(nrow(L) == n, ncol(L) == nxgrid)
    #TODO: check whether supplied L is compatibel with limits argument
  } else {
    L <- switch(integration,
                "simpson" = {
                  # int^b_a f(t) dt = (b-a)/gridlength/3 * [f(a) + 4*f(t_1) + 2*f(t_2) + 4*f(t_3) + 2*f(t_3) +...+ f(b)]
                  ((xind[,nxgrid]-xind[,1])/nxgrid)/3 *
                    matrix(c(1, rep(c(4, 2), length=nxgrid-2), 1), nrow=n, ncol=nxgrid, byrow=T)
                },
                "trapezoidal" = {
                  # int^b_a f(t) dt = .5* sum_i (t_i - t_{i-1}) f(t_i) + f(t_{i-1}) =
                  #   (t_2 - t_1)/2 * f(a=t_1) + sum^{nx-1}_{i=2} ((t_i - t_i-1)/2 + (t_i+1 - t_i)/2) * f(t_i) + ... +
                  #           + (t_nx - t_{nx-1})/2 * f(b=t_n)
                  diffs <- t(apply(xind, 1, diff))
                  .5 * cbind(diffs[,1], t(apply(diffs, 1, filter,
                                                filter=c(1,1)))[,-(nxgrid-1)], diffs[,(nxgrid-1)])
                })
  }
  if(!is.null(limits)){
    if(!is.function(limits)){
      if(!(limits %in% c("s<t","s<=t"))){
        stop("supplied <limits> argument unknown")
      }
      if(limits=="s<t"){
        limits <- function(s, t){
          s < t
        }
      } else {
        if(limits=="s<=t"){
          limits <- function(s, t){
            (s < t) | (s == t)
          }
        }
      }
    }
  }
  #assign unique names based on the given args
  xname <- paste(deparse(substitute(X)), ".mat", sep="")
  xindname <- paste(deparse(substitute(X)), ".smat", sep="")
  yindname <- paste(deparse(substitute(X)), ".tmat", sep="")
  LXname <- paste("L.", deparse(substitute(X)), sep="")

  # make call
  splinefun <- as.symbol(basistype) # if(basistype=="te") quote(te) else quote(s)
  frmls <- formals(getFromNamespace(deparse(splinefun), ns="mgcv"))
  frmls <- modifyList(frmls[names(frmls) %in% names(splinepars)], splinepars)
  call <- as.call(c(
    list(splinefun,
      x = as.symbol(substitute(yindname)),
      y = as.symbol(substitute(xindname)),
      z = as.symbol(substitute(xname)),
      by =as.symbol(substitute(LXname))),
    frmls))
  if(constant_over_yind) {
    call <- as.call(c(
      list(splinefun,
        y = as.symbol(substitute(xindname)),
        z = as.symbol(substitute(xname)),
        by =as.symbol(substitute(LXname))),
      frmls))
  }


  return(list(call=call, xind=xind[1,], L=L, X=X,
              xname=xname, xindname=xindname, yindname=yindname, LXname=LXname,
              limits = limits, use_lag = use_lag))
}#end sff()


