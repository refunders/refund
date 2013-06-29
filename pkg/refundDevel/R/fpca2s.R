fpca2s <-
function(X,npc=NA,center = TRUE, t = NULL,smooth=TRUE){
  
 
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
  
  if(is.na(npc)){
    # use Donoho, Gavish (2013) for estimating suitable number of sv's to extract: 
    beta <- n/m
    
    if(beta > 1 | beta < 1e-3){
      warning("Approximation for \\beta(\\omega) may be invalid.")
    }
    #            ## approx for omega.beta below eq. (25):
    #            betaplus  <- 1 + sqrt(beta)^2 
    #            betaminus <- 1 - sqrt(beta)^2
    #            marcenkopastur <- function(x){
    #                abs(integrate(function(t){
    #                          sqrt((betaplus - t) * (t - betaminus))/(2*pi*t)  
    #                        }, lower=betaminus, upper=x, subdivisions=1e5, 
    #                        stop.on.error = FALSE)$value - 0.5)
    #            }
    #            mu.beta <- optimize(marcenkopastur, 
    #                    interval=c(betaminus, betaplus))$minimum
    #            # eq. (10)
    #            lambda.beta <- sqrt(2*(beta+1) + (8*beta)/(beta + 1 + sqrt(beta^2 + 14*beta +1)))
    #            omega.beta <- lambda.beta/sqrt(mu.beta)
    
    omega.beta <- .56*beta^3 - 0.95*beta^2 + 1.82*beta + 1.43
    y <- sqrt((I-1)*J)*sqrt(lambda)
    rankY <- min(which(cumsum(y[y>0])/sum(y[y>0]) > .995))
    y.med <- median(y)
    
    npc <- min(max(1, sum(y > omega.beta*y.med)),  rankY)

  }
  
  
  if(!is.numeric(npc)) stop("Invalid <npc>.")
  if(npc<1 | npc>min(m,n)) stop("Invalid <npc>.")
  #### end: borrowed from Fabian's code
 print(npc)
  
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
