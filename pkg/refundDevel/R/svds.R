svds <-
function(Y,t=NULL,smooth=TRUE, N=NULL){
  

  ## data: Y, I by J data matrix
  ## t: vector of J
  
  Y <- t(Y) ## becomes J by I
  data_dim <- dim(Y)
  J <- data_dim[1] 
  I <- data_dim[2]  
   
  ### SVD decomposition
  if(J>I){
    VV <- t(Y)%*%Y
    Eigen <- eigen(VV)
    D <- Eigen$values
    sel <- (D>0)
    V <- Eigen$vectors[,sel==1]
    D <- D[sel==1]
    D <- sqrt(D)
    U <- Y%*%V%*%diag(1/D)
  }
  
  if(J<=I){
   UU <- Y%*%t(Y)
   Eigen <- eigen(UU)
   D <- Eigen$values
   U <- Eigen$vectors[,D>0]
   D <- D[D>0]
   D <- sqrt(D)
  }

  lambda <- D^2/I/J
  
  if(smooth==TRUE){
  #### smoothing
  if(is.null(N)) N <- dim(U)[2]
  for(j in 1:N){ 
    #temp = pspline(U[,j],t,knots=knots,p=p,m=m)$fitted.values
    temp = smooth.spline(t,U[,j],all.knots =TRUE)$y
    U[,j] = temp/sqrt(sum(temp^2))
  }
}
return(list(N=N,eigenvectors=U,eigenvalues=lambda,scores = t(Y)%*%U/sqrt(J)))
}
