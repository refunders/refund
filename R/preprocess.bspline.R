preprocess.bspline <- function (funcs, argvals=NULL, nbasis=10){
  
  ## preprocesses predictor functions using splines, exports processed data for additional analysis.
  ## uses unpenalized splines to smooth each curve independently.
  ## handles functions in the form of a matrix
  ##
  ## To Do:
  ## - unify inputs with other functions
  
  if(is.null(argvals)){
    argvals = seq(0, 1, length = dim(funcs)[2])
  }
  
  ## set up design matrix, outcome vector for regression
  Bspline = bs(x = argvals, degree=3, df=nbasis)
  
  ## estimate using OLS
  funcs.processed = matrix(NA, nrow = nrow(funcs), ncol = ncol(funcs))
  for(i in 1:nrow(funcs)){
    funcs.processed[i,] = fitted(lm(funcs[u, ]~ 0+Bspline ))
  }
  
  ret <- list(funcs.processed = funcs.processed)
  return(ret)
  
}