create.prep.func = function(X, argvals=NULL, method = "fpca.sc", options=NULL){
  
  ## creates a function to preprocess predictor functions. this function is
  ## later used to actually preprocess data prior to analysis in scalar-on-function
  ## regression; creating a preprocessor function facilitates prediction.
  ##
  ## To Do:
  ## - unify inputs with other functions
  ##
  
  if(is.null(argvals)){
    argvals = seq(0, 1, length = dim(X)[2])
  }
  
  if(method == "fpca.sc"){
    prep.func = function(newX = NULL, argvals. = argvals, options. = options){
      args. = as.list(options.)
      args.$Y = X; args.$argvals = argvals.; args.$Y.pred = newX
      call = do.call(fpca.sc, args = args.)
      processed = eval(call)$Yhat
      ret = list(processed = processed)
      return(ret)
    }
  } else if(method == "fpca.face"){
    prep.func = function(newX = NULL, argvals. = argvals, options. = options){
      args. = as.list(options.)
      args.$Y = X; args.$argvals = argvals.; args.$Y.pred = newX
      call = do.call(fpca.face, args = args.)
      processed = eval(call)$Yhat
      ret = list(processed = processed)
      return(ret)
    }
  } else if(method == "fpca.ssvd"){
    warning("Preprocssing method `fpca.ssvd` has not been implemented for prediction. Argument argvals not used.")
    prep.func = function(newX = NULL, argvals. = argvals, options. = options){
      args. = as.list(options.)
      args.$Y = X; args.$argvals = argvals.; args.$Y.pred = newX
      call = do.call(fpca.ssvd, args = args.)
      processed = eval(call)$Yhat
      ret = list(processed = processed)
      return(ret)
    }
  } else if(method == "bspline"){
    prep.func = function(newX = NULL, argvals. = argvals, options. = options){
      nbasis = ifelse(is.null(options.$nbasis), 10, options.$nbasis)
      Bspline = bs(x = argvals., degree=3, df=nbasis)
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