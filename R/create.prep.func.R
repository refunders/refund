create.prep.func = function(X, argvals=NULL, method = "fpca.sc", options=NULL, nbasis=10){
  
  ## creates a function to preprocess predictor functions. this function is
  ## later used to actually preprocess data prior to analysis in scalar-on-function
  ## regression; creating a preprocessor function facilitates prediction.
  ##
  ## To Do:
  ## - unify inputs with other functions
  ## - add options
  ##
  ## Assumes:
  ## - argvals is the grid on which X (and newX) are observed
  
  if(is.null(argvals)){
    argvals = seq(0, 1, length = dim(X)[2])
  }
  
  if(method == "fpca.sc"){
    prep.func = function(newX = NULL, argvals. = argvals){
      processed = fpca.sc(Y=X, Y.pred=newX, argvals = argvals.)$Yhat
      ret = list(processed = processed)
      return(ret)
    }
  } else if(method == "fpca.face"){
    prep.func = function(newX = NULL, argvals. = argvals){
      processed = fpca.face(Y=X, Y.pred=newX, argvals = argvals.)$Yhat
      ret = list(processed = processed)
      return(ret)
    }
  } else if(method == "fpca.ssvd"){
    warning("Preprocssing method `fpca.ssvd` has not been implemented for prediction. Argument argvals not used.")
    prep.func = function(newX = NULL, argvals. = argvals){
      processed = fpca.ssvd(Y=newX)$Yhat
      ret = list(processed = processed)
      return(ret)
    }
  } else if(method == "bspline"){
    prep.func = function(newX = NULL, argvals. = argvals, nbasis.=nbasis){
      Bspline = bs(x = argvals., degree=3, df=nbasis.)
      processed = matrix(NA, nrow = nrow(newX), ncol = ncol(newX))
      for(i in 1:nrow(newX)){
        processed[i,] = Bspline %*% coef(lm(newX[i, ]~ 0+Bspline ))
      }
      ret = list(processed = processed)
      return(ret)
    }
  } else if(method == "interpolate") {
    prep.func = function(newX = NULL, argvals. = argvals){
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