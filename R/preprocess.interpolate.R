preprocess.interpolate <- function (funcs, argvals=NULL){
  
  ## preprocesses predictor functions using splines, exports processed data for additional analysis.
  ## uses linear interpolation to preprocess data. if data are interpolated beyond the grid on which a
  ## curve is observed (argvals is wider than the observation for a subject), the closes value is carried
  ## to the end of the domain.
  ##
  ## handles functions in the form of a matrix
  ##
  ## To Do:
  ## - unify inputs with other functions
  
  if(is.null(argvals)){
    argvals = seq(0, 1, length = dim(funcs)[2])
  }
  
  ## use approx function to interpolate curves
  funcs.processed = matrix(NA, nrow = nrow(funcs), ncol = ncol(funcs))
  for(i in 1:nrow(funcs)){
    obs.pts = which(!is.na(funcs[i,]))
    x = argvals[obs.pts]; y = funcs[i,obs.pts]
    funcs.processed[i,] = approx(x = x, y = y, xout = argvals, rule = 2)$y
  }
  
  ret <- list(funcs.processed = funcs.processed)
  return(ret)
  
}