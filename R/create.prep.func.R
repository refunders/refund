#' Construct a function for preprocessing functional predictors
#'
#' Prior to using functions \code{X} as predictors in a scalar-on-function regression, it is often
#' necessary to presmooth curves to remove measurement error or interpolate to a common grid. This
#' function creates a function to do this preprocessing depending on the method specified.
#' 
#' @param X an \code{N} by \code{J=ncol(argvals)} matrix of function evaluations
#'   \eqn{X_i(t_{i1}),., X_i(t_{iJ}); i=1,.,N.} For FPCA-based processing methods, these functions are
#'   used to define the eigen decomposition used to preprocess current and future data (for example, in
#'   \code{\link{predict.pfr}})
#' @param argvals matrix (or vector) of indices of evaluations of \eqn{X_i(t)}; i.e. a matrix with
#'   \emph{i}th row \eqn{(t_{i1},.,t_{iJ})}
#' @param method character string indicating the preprocessing method. Options
#'   are \code{"fpca.sc"}, \code{"fpca.face"}, \code{"fpca.ssvd"}, \code{"bspline"},
#'   and \code{"interpolate"}. The first three use the corresponding existing function;
#'   \code{"bspline"} uses an (unpenalized) cubic bspline smoother with \code{nbasis} basis 
#'   functions; \code{"interpolate"} uses linear interpolation.
#' @param options list of options passed to the preprocessing method; as an example, options for \code{fpca.sc}
#'   include \code{pve}, \code{nbasis}, and \code{npc}.
#' 
#' @return a function that returns the preprocessed functional predictors, with arguments
#'   \item{newX}{The functional predictors to process}
#'   \item{argvals.}{Indices of evaluation of \code{newX}}
#'   \item{options.}{Any options needed to preprocess the predictor functions}
#' 
#' @author Jeff Goldsmith \email{ajg2202@@cumc.columbia.edu}
#' @seealso \code{\link{pfr}}, \code{\link{fpca.sc}}, \code{\link{fpca.face}}, \code{\link{fpca.ssvd}}

create.prep.func = function(X, argvals=seq(0, 1, length = ncol(X)),
                            method = c("fpca.sc", "fpca.face", "fpca.ssvd",
                                       "bspline", "interpolate"),
                            options=NULL) {
  method <- match.arg(method)
  
  if (method == "fpca.sc"){
    function(newX = NULL, argvals. = argvals, options. = options) {
      args. = c(as.list(options.), list(Y=X, argvals=argvals., Y.pred=newX))
      eval(do.call(fpca.sc, args = args.))$Yhat
    }
  } else if (method == "fpca.face"){
    function(newX = NULL, argvals. = argvals, options. = options){
      args. = c(as.list(options.), list(Y=X, argvals=argvals., Y.pred=newX))
      eval(do.call(fpca.face, args = args.))$Yhat
    }
  } else if (method == "fpca.ssvd"){
    warning("Preprocssing method `fpca.ssvd` has not been implemented for prediction. Argument argvals not used.")
    function(newX = NULL, argvals. = argvals, options. = options){
      #; args.$argvals = argvals.
      #; args.$Y.pred = newX
      args. = c(as.list(options.), list(Y=X))
      eval(do.call(fpca.ssvd, args = args.))$Yhat
    }
  } else if (method == "bspline"){
    function(newX = NULL, argvals. = argvals, options. = options){
      nbasis = ifelse(is.null(options.$nbasis), 10, options.$nbasis)
      Bspline = bs(x = argvals., degree=3, df=nbasis)
      processed = matrix(NA, nrow = nrow(newX), ncol = ncol(newX))
      for(i in 1:nrow(newX)){
        processed[i,] = Bspline %*% coef(lm(newX[i, ]~ 0+Bspline ))
      }
      processed
    }
  } else if (method == "interpolate") {
    function(newX = NULL, argvals. = argvals){
      processed = matrix(NA, nrow = nrow(newX), ncol = ncol(newX))
      for(i in 1:nrow(newX)){
        obs.pts = which(!is.na(newX[i,]))
        x = argvals.[obs.pts]; y = newX[i,obs.pts]
        processed[i,] = approx(x = x, y = y, xout = argvals., rule = 2)$y
      }
      processed
    }
  } else {
    stop("Unrecognized preprocessing method")
  }
}