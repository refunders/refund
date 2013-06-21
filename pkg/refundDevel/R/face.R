face <-
function(Y,t=NULL,knots=35,p=3,m=2,lambda=NULL, percentage = 0.95, 
         score.method = "int", search.grid=TRUE,search.length=100,
         method="L-BFGS-B", lower=-20,upper=20, control=NULL){
  
  ## data: Y, I by J data matrix
  ## t:  vector of J
  ## knots: to specify either the number of knots or the vectors of knots;
  ##        defaults to 35;
  ## p: the degree of B-splines;
  ## m: the order of difference penalty
  ## lambda: user-selected smoothing parameter
  ## method: see R function "optim" 
  ## lower, upper, control: see R function "optim"
  #require(splines)  
  #source("pspline.setting.R") 
  Y <- t(Y) ## becomes J by I
  data_dim <- dim(Y)
  J <- data_dim[1] 
  I <- data_dim[2]  
  
  if(is.null(t))  t <- (1:J)/J-1/2/J ## if NULL, assume equally spaced
  
  p.p <- p
  m.p <- m
  if(length(knots)==1){
    K.p <- knots
    knots <- seq(-p.p,K.p+p.p,length=K.p+1+2*p.p)/K.p
    knots <- knots*(max(t)-min(t)) + min(t)
  }
  if(length(knots)>1) K.p <- length(knots)-2*p.p-1
  c.p <- K.p + p.p
  ######### precalculation for smoothing #############
  List <- pspline.setting(t,knots,p.p,m.p)
  AS <- List$A
  s <- List$s
  Sigi.sqrt <- List$Sigi.sqrt
  U <- List$U
  
  MM <- function(A,s,option=1){
    if(option==2)
      return(A*(s%*%t(rep(1,dim(A)[2]))))
    if(option==1)
      return(A*(rep(1,dim(A)[1])%*%t(s)))
  }
  
  ######## precalculation for missing data ########
  imputation <- FALSE
  Niter.miss <- 1
  
  Index.miss <- is.na(Y)
  if(sum(Index.miss)>0){
  Y[Index.miss] <- 0
  Y0 <- matrix(NA,c.p,I)

  
  imputation <- TRUE
  Niter.miss <- 100
  }
  convergence.vector <- rep(0,Niter.miss);
  iter.miss <- 1
  
  while(iter.miss <= Niter.miss&&convergence.vector[iter.miss]==0){
  ###################################################
  ######## Transform the Data           #############
  ###################################################
  Ytilde <- as.matrix(t(AS)%*%Y)
  C_diag <- rowSums(Ytilde^2)
  if(iter.miss==1) Y0 = Ytilde
  ###################################################
  ########  Select Smoothing Parameters #############
  ###################################################
  Y_square <- sum(Y^2)
  Ytilde_square <- sum(Ytilde^2)
  
  face_gcv <- function(x){
    lambda <- exp(x)
    lambda_s <- (lambda*s)^2/(1 + lambda*s)^2
    gcv <- sum(C_diag*lambda_s) - Ytilde_square + Y_square
    trace <- sum(1/(1+lambda*s))
    gcv <- gcv/(1-trace/J)^2
    return(gcv)
  }
  
  if(is.null(lambda)){
    
    if(!search.grid){
      fit <- optim(0,face_gcv,method=method,lower=lower,upper=upper,control=control)
      if(fit$convergence>0) {
        expression <- paste("Smoothing failed! The code is:",fit$convergence)
        print(expression)
      }
      
      lambda <- exp(fit$par)
    }
    
    if(search.grid){
      Lambda <- seq(lower,upper,length=search.length)
      Length <- length(Lambda)
      Gcv <- rep(0,Length)
      for(i in 1:Length) 
        Gcv[i] <- face_gcv(Lambda[i])
      i0 <- which.min(Gcv)
      lambda <- exp(Lambda[i0])
    }
  }

    YS <- MM(Ytilde,1/(1+lambda*s),2)  
  
  ###################################################
  ####  Eigendecomposition of Smoothed Data #########
  ###################################################
  if(c.p <= I){
    temp <- YS%*%t(YS)/I
    Eigen <- eigen(temp)
    A <- Eigen$vectors
    Sigma <- Eigen$values/J
  }
  
  if(c.p > I){
    SVD <- svd(YS)
    A <- SVD$u
    Sigma <- (SVD$d)^2/J/I
  }
  
  if(iter.miss>1&&iter.miss< Niter.miss) {
    
    diff <- norm(YS-YS.temp,"F")/norm(YS,"F");
    if(diff <= 10^(-2)) 
      convergence.vector[iter.miss+1] <- 1
  }
  YS.temp <- YS
  iter.miss <- iter.miss + 1
  N <- min(I,c.p)
  d <- Sigma[1:N]
  d <- d[d>0]
  per <- cumsum(d)/sum(d)
  N <- 1
  while(per[N]<percentage) N <- N + 1
  #########################################
  #######     Principal  Scores   #########
  ########   data imputation      #########
  #########################################
  if(imputation==T){
    
    A.N <- A[,1:N]
    d <- Sigma[1:N]
    if(score.method=="int") sigmahat2 <- 0
    if(score.method=="blup") sigmahat2 <- sum(Y[!Index.miss]^2)/(J*I - sum(Index.miss)) -sum(Sigma[Sigma>0])
    sigmahat2 <- max(0,sigmahat2)
   
    
    Xi <- t(A.N)%*%Ytilde
    Xi <- as.matrix(AS%*%((A.N%*%diag(d/(d+sigmahat2/J)))%*%Xi))
    Y <- Y*(1-Index.miss) + Xi*Index.miss
 
    if(sum(is.na(Y))>0) print("error")
  }
  #if(iter.miss%%10==0) print(iter.miss)
   
  }## end of while loop
  #if(sum(Index.miss)>0) cat("The number of iterations is:",iter.miss,"\n")
  
  ### now calculate scores
  Ytilde <- as.matrix(t(AS)%*%Y)
  if(score.method=="int") Xi <- t(Ytilde)%*%A[,1:N]/sqrt(J)
  if(score.method=="blup"){Xi <- t(Ytilde)%*%A[,1:N]/sqrt(J)
                          Xi <- Xi%*%diag(Sigma[1:N]/(Sigma[1:N] + sigmahat2/J))
  }
    
  results <- list(N = N, eigenvectors=as.matrix(AS%*%A[,1:N]), eigenvalues = Sigma[1:N],scores = Xi, lambda=lambda)
  return(results)      	                	
}
