#' Construct an FLM regression term
#'
#' Defines a term \eqn{\int_{T}\beta(t)X_i(t)dt} for inclusion in an \code{mgcv::gam}-formula (or
#' \code{\link{bam}} or \code{\link{gamm}} or \code{gamm4:::gamm}) as constructed by
#' \code{\link{pfr}}, where \eqn{\beta(t)} is an unknown coefficient
#' function and \eqn{X_i(t)} is a functional predictor on the closed interval
#' \eqn{T}. See
#' \code{\link{smooth.terms}} for a list of basis and penalty options; the
#' default is thin-plate regression splines, as this is the default option
#' for \code{\link[mgcv]{s}}.
#' 
#' @param X functional predictors, expressed as an \code{N} by \code{J} matrix,
#'   where \code{N} is the number of columns and \code{J} is the number of
#'   evaluation points. May include missing/sparse functions, which are
#'   indicated by \code{NA} values.
#' @param argvals indices of evaluation of \code{X}, i.e. \eqn{(t_{i1},.,t_{iJ})} for
#'   subject \eqn{i}. May be entered as either a length-\code{J} vector, or as
#'   an \code{N} by \code{J} matrix. Indices may be unequally spaced. Entering
#'   as a matrix allows for different observations times for each subject.
#' @param xind same as argvals. It will not be supported in the next version of refund.
#' @param integration method used for numerical integration. Defaults to \code{"simpson"}'s rule
#'   for calculating entries in \code{L}. Alternatively and for non-equidistant grids,
#'   \code{"trapezoidal"} or \code{"riemann"}.
#' @param L an optional \code{N} by \code{ncol(argvals)} matrix giving the weights for the numerical
#'   integration over \code{t}. If present, overrides \code{integration}.
#' @param presmooth string indicating the method to be used for preprocessing functional predictor prior 
#'   to fitting. Options are \code{fpca.sc}, \code{fpca.face}, \code{fpca.ssvd}, \code{fpca.bspline}, and 
#'   \code{fpca.interpolate}. Defaults to \code{NULL} indicateing no preprocessing. See
#'   \code{\link{create.prep.func}}.
#' @param presmooth.opts list including options passed to preprocessing method
#'   \code{\link{create.prep.func}}.
#' @param ... optional arguments for basis and penalization to be passed to
#'   \code{mgcv::s}. These could include, for example,
#'   \code{"bs"}, \code{"k"}, \code{"m"}, etc. See \code{\link[mgcv]{s}} for details.
#' 
#' @return a list with the following entries
#'   \item{\code{call}}{a \code{call} to \code{te} (or \code{s}, \code{t2}) using the appropriately
#'     constructed covariate and weight matrices}
#'   \item{\code{argvals}}{the \code{argvals} argument supplied to \code{lf}}
#'   \item{\code{L}}{the  matrix of weights used for the integration}
#'   \item{\code{xindname}}{the name used for the functional predictor variable in the \code{formula}
#'     used by \code{mgcv}}
#'   \item{\code{tindname}}{the name used for \code{argvals} variable in the \code{formula} used by \code{mgcv}}
#'   \item{\code{LXname}}{the name used for the \code{L} variable in the \code{formula} used by \code{mgcv}}
#'   \item{\code{presmooth}}{the \code{presmooth} argument supplied to \code{lf}}
#'   \item{\code{prep.func}}{a function that preprocesses data based on the preprocessing method specified in \code{presmooth}. See
#'     \code{\link{create.prep.func}}}
#' @author Mathew W. McLean \email{mathew.w.mclean@@gmail.com}, Fabian Scheipl,
#'   and Jonathan Gellar
#' 
#' @references
#' Goldsmith, J., Bobb, J., Crainiceanu, C., Caffo, B., and Reich, D. (2011).
#' Penalized functional regression. \emph{Journal of Computational and Graphical
#' Statistics}, 20(4), 830-851.
#' 
#' Goldsmith, J., Crainiceanu, C., Caffo, B., and Reich, D. (2012). Longitudinal
#' penalized functional regression for cognitive outcomes on neuronal tract
#' measurements. \emph{Journal of the Royal Statistical Society: Series C},
#' 61(3), 453-469.
#' 
#' @seealso \code{\link{pfr}}, \code{\link{af}}, mgcv's \code{\link{smooth.terms}}
#'  and \code{\link{linear.functional.terms}}, \code{\link{pfr}} for additonal examples
#' @importFrom fda create.bspline.basis smooth.basisPar eval.fd
#' @importFrom utils getFromNamesapce modifyList
lf <- function(X, argvals = seq(0, 1, l = ncol(X)), xind = NULL,
               integration = c("simpson", "trapezoidal", "riemann"),
               L = NULL, presmooth = NULL, presmooth.opts = NULL, ...) {
  
  # Catch if lf_old syntax is used
  dots <- list(...)
  dots.unmatched <- names(dots)[!(names(dots) %in% names(formals(s)))]
  if (any(dots.unmatched %in% names(formals(lf_old))) | is.logical(presmooth)) {
    warning(paste0("The interface for lf() has changed, see ?lf for details. ",
                   "This interface will not be supported in the next ",
                   "refund release."))
    # Call lf_old()
    call <- sys.call()
    call[[1]] <- as.symbol("lf_old")
    ret <- eval(call, envir=parent.frame())
    return(ret)
  }
  
  if (!is.null(xind)) {
    cat("Argument xind is placed by argvals. xind will not be supported in the next
        version of refund.")
    argvals = xind
  }
  xind = argvals
  
  n=nrow(X)
  nt=ncol(X)
  integration <- match.arg(integration)

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

  if(!is.null(presmooth)){
    # create preprocessing function
    prep.func = create.prep.func(X = X, argvals = xind[1,], method = presmooth, options = presmooth.opts)
    
    # preprocess data
    X <- prep.func(newX = X)$processed
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
  call <- as.call(c(list(splinefun, x = as.symbol(substitute(tindname)),
                         by = as.symbol(substitute(LXname))), dots))
  res <-list(call = call, data = data, xind = xind[1,], L = L, tindname=tindname,
             LXname=LXname,presmooth=presmooth)
  if(!is.null(presmooth)) {res$prep.func <- prep.func} 
  return(res)
}
