fpca2s <-
function(X,npc=NA,center = TRUE, t = NULL,smooth=TRUE){
  
  if(is.na(npc)){
    npc <- getNPC.DonohoGavish(X)
  }
  
  ## data: X, I by J data matrix 
  ## t: vector of J
  data_dim <- dim(X)
  J <- data_dim[1] 
  I <- data_dim[2]  
  
  if(is.null(t)) t <- seq(0, 1, length=J)
  if(center) X <- apply(X,2,function(x){x-mean(x,na.rm=TRUE)})
  #### start: borrowed from Fabian's code
  n <- I
  m <- J
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
  if(npc<1 | npc>min(m,n)) stop("Invalid <npc>.")
  #### end: borrowed from Fabian's code
  message("Extracted ", npc, " smooth components.")
  
  if(smooth==TRUE){
  #### smoothing
  for(j in 1:npc){ 
    #temp = pspline(U[,j],t,knots=knots,p=p,m=m)$fitted.values
    temp = smooth.spline(t,U[,j],all.knots =TRUE)$y
    U[,j] = temp/sqrt(sum(temp^2))
  }
}
return(list(npc=npc,eigenvectors=U,eigenvalues=lambda,scores = X%*%U[,1:npc]/sqrt(J)))
}
