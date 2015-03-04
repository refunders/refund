#' Construct a function for preprocessing functional predictors
#'
#' Prior to using functions \code{X} as predictors in a scalar-on-function regression, it is often
#' necessary to presmooth curves to remove measurement error or interpolate to a common grid. This
#' function creates a function to do this preprocessing depending on the method specified.
#' @param X an \code{N} by \code{J=ncol(argvals)} matrix of function evaluations
#' \eqn{X_i(t_{i1}),., X_i(t_{iJ}); i=1,.,N.} For FPCA-based processing methods, these functions are
#' used to define the eigen decomposition used to preprocess current and future data (for example, in
#' \code{\link{predict.pfr}})
#' @param argvals matrix (or vector) of indices of evaluations of \eqn{X_i(t)}; i.e. a matrix with
#' \emph{i}th row \eqn{(t_{i1},.,t_{iJ})}
#' @param method approach used to preprocess curves. Options are \code{fpca.sc}, \code{fpca.face}, 
#' \code{fpca.ssvd}, \code{fpca.bspline}, and \code{fpca.interpolate}. The first three are use existing
#' functions; \code{fpca.bspline} uses an (unpenalized) cubic bspline smoother with \code{nbasis} basis 
#' funcitons; \code{fpca.interpolate} uses linear interpolation.
#' @param options list of options passed to the preprocessing method; as an example, options for \code{fpca.sc}
#' include \code{pve}, \code{nbasis}, and \code{npc}.
#' @return a list with the following entries
#' \enumerate{
#' \item \code{prep.func} - a function that will preprocess functional predictors using on the method and options 
#' specified and, where appropriate, the original data. This function can be used for current and future 
#' functional predictors.
#' }
#' @author Jeff Goldsmith \email{ajg2202@@cumc.columbia.edu}
#' @seealso \code{\link{pfr}}, \code{\link{fpca.sc}}, \code{\link{fpca.face}}, \code{\link{fpca.ssvd}}
create.prep.func = function(X, argvals=NULL, method = "fpca.sc", options=NULL){
  
  if(is.null(argvals)){
    argvals = seq(0, 1, length = dim(X)[2])
  }
  
  prep.func <- if (method == "fpca.sc"){
    function(newX = NULL, argvals. = argvals, options. = options) {
      args. = as.list(options.)
      args.$Y = X; args.$argvals = argvals.; args.$Y.pred = newX
      call = do.call(fpca.sc, args = args.)
      processed = eval(call)$Yhat
      ret = list(processed = processed)
      return(ret)
    }
  } else if (method == "fpca.face"){
    function(newX = NULL, argvals. = argvals, options. = options){
      args. = as.list(options.)
      args.$Y = X; args.$argvals = argvals.; args.$Y.pred = newX
      call = do.call(fpca.face, args = args.)
      processed = eval(call)$Yhat
      ret = list(processed = processed)
      return(ret)
    }
  } else if (method == "fpca.ssvd"){
    warning("Preprocssing method `fpca.ssvd` has not been implemented for prediction. Argument argvals not used.")
    function(newX = NULL, argvals. = argvals, options. = options){
      args. = as.list(options.)
      args.$Y = X; args.$argvals = argvals.; args.$Y.pred = newX
      call = do.call(fpca.ssvd, args = args.)
      processed = eval(call)$Yhat
      ret = list(processed = processed)
      return(ret)
    }
  } else if (method == "bspline"){
    function(newX = NULL, argvals. = argvals, options. = options){
      nbasis = ifelse(is.null(options.$nbasis), 10, options.$nbasis)
      Bspline = bs(x = argvals., degree=3, df=nbasis)
      processed = matrix(NA, nrow = nrow(newX), ncol = ncol(newX))
      for(i in 1:nrow(newX)){
        processed[i,] = Bspline %*% coef(lm(newX[i, ]~ 0+Bspline ))
      }
      ret = list(processed = processed)
      return(ret)
    }
  } else if (method == "interpolate") {
    function(newX = NULL, argvals. = argvals){
      processed = matrix(NA, nrow = nrow(newX), ncol = ncol(newX))
      for(i in 1:nrow(newX)){
        obs.pts = which(!is.na(newX[i,]))
        x = argvals.[obs.pts]; y = newX[i,obs.pts]
        processed[i,] = approx(x = x, y = y, xout = argvals., rule = 2)$y
      }
      ret = list(processed = processed)
      return(ret)
    }
  }
  return(prep.func)
}