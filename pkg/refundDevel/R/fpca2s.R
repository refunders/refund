fpca2s <-
function(Y,npc=NA,center = TRUE, argvals = NULL,smooth=TRUE){
  
  ## data: X, I by J data matrix 
  ## argvals: vector of J
  X <- Y
  data_dim <- dim(X)
  I <- data_dim[1] 
  J <- data_dim[2]  
  
  if(is.na(npc)){
    npc <- getNPC.DonohoGavish(X)
  }
  
  if(is.null(argvals)) argvals <- seq(0, 1, length=J)
  
  meanX <- rep(0,J)
  if(center) {
    meanX <- apply(X,2,function(x) mean(x,na.rm=TRUE))
    meanX <- smooth.spline(argvals,meanX,all.knots =TRUE)$y
    X <- t(t(X)-meanX)
  }
  ### SVD decomposition
  if(J>I){
    VV <- X%*%t(X)
    Eigen <- eigen(VV)
    D <- Eigen$values
    sel <- (D>0)
    V <- Eigen$vectors[,sel==1]
    D <- D[sel==1]
    D <- sqrt(D)
    U <- t(X)%*%V%*%diag(1/D)
  }
  
  if(J<=I){
    UU <- t(X)%*%X
    Eigen <- eigen(UU)
    D <- Eigen$values
    U <- Eigen$vectors[,D>0]
    D <- D[D>0]
    D <- sqrt(D)
  }
  
  lambda <- D^2/(I-1)/J
    
  if(!is.numeric(npc)) stop("Invalid <npc>.")
  if(npc<1 | npc>min(I,J)) stop("Invalid <npc>.")
  #### end: borrowed from Fabian's code
  message("Extracted ", npc, " smooth components.")
  
  if(smooth==TRUE){
    #### smoothing
    for(j in 1:npc){ 
      #temp = pspline(U[,j],argvals,knots=knots,p=p,m=m)$fitted.values
      temp = smooth.spline(argvals,U[,j],all.knots =TRUE)$y
      U[,j] = temp/sqrt(sum(temp^2))
    }
  }
  eigenvectors=U[,1:npc]
  eigenvalues = lambda[1:npc]
  scores = X%*%U[,1:npc]/sqrt(J)
  
  Yhat = t(eigenvectors%*%t(scores)*sqrt(J) + meanX)
  return(list(Yhat = Yhat,scores = scores, mu = meanX,eigenvectors=eigenvectors,eigenvalues=eigenvalues, npc=npc))
}
