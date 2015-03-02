preprocess.fpca <- function (funcs, argvals=NULL, kz = NULL, nbasis=10, smooth.option="fpca.sc", pve=0.99){
  
  ## preprocesses predictor functions using fpca, exports processed data for additional analysis.
  ## handles functions in the form of a matrix
  ##
  ## To Do:
  ## - unify inputs with other functions
  
  if(is.null(argvals)){
    argvals = seq(0, 1, length = dim(funcs)[2])
  }
  
  ## obtain FPCA expansion each predictor
  if(smooth.option=="fpca.sc"){
    funcs.processed = fpca.sc(Y=funcs, argvals=argvals, pve=pve, nbasis=nbasis, npc=kz)$Yhat 
  }

  if (smooth.option=="fpca.face"){ 
    funcs.processed = fpca.face(Y=funcs, argvals=argvals, knots=nbasis, pve=pve, npc=kz)$Yhat
  }

  if (smooth.option=="fpca.ssvd"){ 
    funcs.processed = fpca.ssvd(Y=funcs)$Yhat
  }
  
  ret <- list(funcs.processed = funcs.processed)
  return(ret)

}

