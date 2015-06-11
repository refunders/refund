#' Construct an FGAM regression term
#'
#' Defines a term \eqn{\int_{T}F(X_i(t),t)dt} for inclusion in an \code{mgcv::gam}-formula (or
#' \code{\link{bam}} or \code{\link{gamm}} or \code{gamm4:::gamm}) as constructed by
#' \code{\link{pfr}}, where \eqn{F(x,t)} is an unknown smooth bivariate function and \eqn{X_i(t)}
#' is a functional predictor on the closed interval \eqn{T}. See \code{\link{smooth.terms}}
#' for a list of bivariate basis and penalty options; the default is a tensor
#' product basis with marginal cubic regression splines for estimating \eqn{F(x,t)}.
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
#' @param basistype defaults to \code{"te"}, i.e. a tensor product spline to represent \eqn{F(x,t)} Alternatively,
#'   use \code{"s"} for bivariate basis functions (see \code{\link{s}}) or \code{"t2"} for an alternative
#'   parameterization of tensor product splines (see \code{\link{t2}})
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
#' @param Xrange numeric; range to use when specifying the marginal basis for the \emph{x}-axis.  It may
#'   be desired to increase this slightly over the default of \code{range(X)} if concerned about predicting
#'   for future observed curves that take values outside of \code{range(X)}
#' @param Qtransform logical; should the functional be transformed using the empirical cdf and
#'   applying a quantile transformation on each column of \code{X} prior to fitting?  This ensures
#'   \code{Xrange=c(0,1)}.  If \code{Qtransform=TRUE} and \code{presmooth=TRUE}, presmoothing is done prior
#'   to transforming the functional predictor
#' @param ... optional arguments for basis and penalization to be passed to the
#'   function indicated by \code{basistype}. These could include, for example,
#'   \code{"bs"}, \code{"k"}, \code{"m"}, etc. See \code{\link{te}} or
#'   \code{\link{s}} for details.
#' 
#' @return A list with the following entries:
#'   \item{\code{call}}{a \code{"call"} to \code{te} (or \code{s}, \code{t2}) using the appropriately
#'     constructed covariate and weight matrices.}
#'   \item{\code{argvals}}{the \code{argvals} argument supplied to \code{af}}
#'   \item{\code{L}}{the  matrix of weights used for the integration}
#'   \item{\code{xindname}}{the name used for the functional predictor variable in the \code{formula} used by \code{mgcv}}
#'   \item{\code{tindname}}{the name used for \code{argvals} variable in the \code{formula} used by \code{mgcv}}
#'   \item{\code{Lname}}{the name used for the \code{L} variable in the \code{formula} used by \code{mgcv}}
#'   \item{\code{presmooth}}{the \code{presmooth} argument supplied to \code{af}}
#'   \item{\code{Qtranform}}{the \code{Qtransform} argument supplied to \code{af}}
#'   \item{\code{Xrange}}{the \code{Xrange} argument supplied to \code{af}}
#'   \item{\code{ecdflist}}{a list containing one empirical cdf function from applying \code{\link{ecdf}}
#'     to each (possibly presmoothed) column of \code{X}.  Only present if \code{Qtransform=TRUE}}
#'   \item{\code{prep.func}}{a function that preprocesses data based on the preprocessing method specified in \code{presmooth}. See
#'     \code{\link{create.prep.func}}}
#' 
#' @examples
#' 
#' data(DTI)
#' ## only consider first visit and cases (no PASAT scores for controls)
#' DTI1 <- DTI[DTI$visit==1 & DTI$case==1,]
#' DTI2 <- DTI1[complete.cases(DTI1),]
#'
#' ## fit FGAM using FA measurements along corpus callosum
#' ## as functional predictor with PASAT as response
#' ## using 8 cubic B-splines for marginal bases with third
#' ## order marginal difference penalties
#' ## specifying gamma > 1 enforces more smoothing when using
#' ## GCV to choose smoothing parameters
#' fit1 <- pfr(pasat ~ af(cca, k=c(8,8), m=list(c(2,3), c(2,3))),
#'             method="GCV.Cp", gamma=1.2, data=DTI1)
#' vis.pfr(fit1)
#'
#'
#' ## fgam term for the cca measurements plus an flm term for the rcst measurements
#' ## leave out 10 samples for prediction
#' test <- sample(nrow(DTI2), 10)
#' fit2 <- pfr(pasat ~ af(cca, k=c(7,7), m=list(c(2,3), c(2,3))) +
#'                     lf(rcst, k=7, m=c(2,2), bs="ps"),
#'             method="GCV.Cp", gamma=1.2, data=DTI2[-test,])
#' par(mfrow=c(1,2))
#' plot(fit2, scheme=2, rug=FALSE)
#' vis.pfr(fit2, af.term="cca", xval=.6)
#' pred <- predict(fit2, newdata = DTI2[test,], type='response', PredOutOfRange = TRUE)
#' sqrt(mean((DTI2$pasat[test] - pred)^2))
#' 
#' ## Try to predict the binary response disease status (case or control)
#' ##   using the quantile transformed measurements from the rcst tract
#' ##   with a smooth component for a scalar covariate that is pure noise
#' DTI3 <- DTI[DTI$visit==1,]
#' DTI3 <- DTI3[complete.cases(DTI3$rcst),]
#' z1 <- rnorm(nrow(DTI3))
#' fit3 <- pfr(case ~ af(rcst, k=c(7,7), m = list(c(2, 1), c(2, 1)), Qtransform=TRUE) +
#'                     s(z1, k = 10), family="binomial", select=TRUE, data=DTI3)
#' par(mfrow=c(1,2))
#' plot(fit3, scheme=2, rug=FALSE)
#' abline(h=0, col="green")
#' vis.pfr(fit3, af.term="rcst", plot.type="contour")
#' 
#' @author Mathew W. McLean \email{mathew.w.mclean@@gmail.com}, Fabian Scheipl,
#'   and Jonathan Gellar
#' @references McLean, M. W., Hooker, G., Staicu, A.-M., Scheipl, F., and Ruppert, D. (2014). Functional
#' generalized additive models. \emph{Journal of Computational and Graphical Statistics}, \bold{23 (1)},
#' pp. 249-269.  Available at \url{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3982924}.
#' @seealso \code{\link{pfr}}, \code{\link{lf}}, mgcv's \code{\link{linear.functional.terms}},
#' \code{\link{pfr}} for examples
#' @importFrom stats ecdf
#' @importFrom fda int2Lfd smooth.basisPar eval.fd create.bspline.basis
#' @importFrom utils modifyList getFromNamespace

af <- function(X, argvals = seq(0, 1, l = ncol(X)), xind = NULL,
               basistype = c("te","t2", "s"),
               integration = c("simpson", "trapezoidal", "riemann"),
               L = NULL, presmooth = NULL, presmooth.opts = NULL, 
               Xrange=range(X), Qtransform=FALSE, ...) {
  
  # Catch if af_old syntax is used
  dots <- list(...)
  dots.unmatched <- names(dots)[!(names(dots) %in%
                                    names(formals(eval(basistype))))]
  if (any(dots.unmatched %in% names(formals(af_old))) |
        is.logical(presmooth)) {
    warning(paste0("The interface for af() has changed, see ?af for details. ",
                   "This interface will not be supported in the next ",
                   "refund release."))
    # Call af_old()
    call <- sys.call()
    call[[1]] <- as.symbol("af_old")
    ret <- eval(call, envir=parent.frame())
    return(ret)
  }
  
  if (!is.null(xind)) {
    argvals = xind
    cat("Warnings: xind argument is renamed as argvals and will not be supported
        in the next version of refund.")
  }
  
  xind = argvals
  n=nrow(X)
  nt=ncol(X)
  basistype <- match.arg(basistype)
  integration <- match.arg(integration)

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

  if(!is.null(presmooth)){
    # create and executepreprocessing function
    prep.func = create.prep.func(X = X, argvals = xind[1,], method = presmooth,
                                 options = presmooth.opts)
    X <- prep.func(newX = X)

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
  } else {
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
  call <- as.call(c(list(splinefun, z = as.symbol(substitute(tindname)),
                         x = as.symbol(substitute(xindname)),
                         by = as.symbol(substitute(Lname))), dots))
  res <-list(call = call, data = data, xind = xind[1,], L = L, xindname = xindname, tindname=tindname,
             Lname=Lname,Qtransform=Qtransform,presmooth=presmooth,Xrange=Xrange)
  if(Qtransform) res$ecdflist <- ecdflist
  if(!is.null(presmooth)) {res$prep.func <- prep.func} 
  return(res)
}
