fpca.face <-
function(Y=NULL,ydata=NULL,Y.pred = NULL,argvals=NULL,pve = 0.99, npc  = NULL,
         var = FALSE, simul = FALSE, sim.alpha = 0.95,
         center=TRUE,knots=35,p=3,m=2,lambda=NULL,alpha = 1,
         search.grid=TRUE,search.length=100,
         method="L-BFGS-B", lower=-20,upper=20, control=NULL){
  
  ## data: Y, I by J data matrix, functions on rows
  ## argvals:  vector of J
  ## knots: to specify either the number of knots or the vectors of knots;
  ##        defaults to 35;
  ## p: the degree of B-splines;
  ## m: the order of difference penalty
  ## lambda: user-selected smoothing parameter
  ## method: see R function "optim" 
  ## lower, upper, control: see R function "optim"
  #require(Matrix)
  #source("pspline.setting.R")
  stopifnot(!is.null(Y))
  stopifnot(is.matrix(Y))
  data_dim <- dim(Y)
  I <- data_dim[1] ## number of subjects
  J <- data_dim[2] ## number of obs per function
  
  if(is.null(argvals))  argvals <- (1:J)/J-1/2/J ## if NULL, assume equally spaced
  
  meanX <- rep(0,J)
  if(center) {##center the functions
    meanX <- apply(Y,2,function(x) mean(x,na.rm=TRUE))
    meanX <- smooth.spline(argvals,meanX,all.knots =TRUE)$y
    Y <- t(t(Y)- meanX)
  }
  
  ## specify the B-spline basis: knots
  p.p <- p
  m.p <- m
  if(length(knots)==1){
    if(knots+p.p>=J) cat("Too many knots!\n")
    stopifnot(knots+p.p<J)
    
    K.p <- knots
    knots <- seq(-p.p,K.p+p.p,length=K.p+1+2*p.p)/K.p
    knots <- knots*(max(argvals)-min(argvals)) + min(argvals)
  }
  if(length(knots)>1) K.p <- length(knots)-2*p.p-1
  if(K.p>=J) cat("Too many knots!\n")
  stopifnot(K.p <J)
  c.p <- K.p + p.p
  ######### precalculation for smoothing #############
  List <- pspline.setting(argvals,knots,p.p,m.p)
  B <- List$B
  Bt <- Matrix(t(as.matrix(B)))
  s <- List$s
  Sigi.sqrt <- List$Sigi.sqrt
  U <- List$U
  A0 <- Sigi.sqrt%*%U
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
    num.miss <- apply(Y,2,function(x){sum(is.na(x))})
    for(i in 1:I){
      if(num.miss[i]>0){
        y <- Y[i,]
        seq <- (1:J)[!is.na(y)]
        seq2 <-(1:J)[is.na(y)]
        t1 <- argvals[seq]
        t2 <- argvals[seq2]
        fit <- smooth.spline(t1,y[seq])
        temp <- predict(fit,t2,all.knots=TRUE)$y
        if(max(t2)>max(t1)) temp[t2>max(t1)] <- mean(y[seq])#y[seq[length(seq)]]
        if(min(t2)<min(t1)) temp[t2<min(t1)] <- mean(y[seq])#y[seq[1]]
        Y[i,seq2] <- temp
      }
    }
    
  Y0 <- matrix(NA,c.p,I)
  imputation <- TRUE
  Niter.miss <- 100
  }
  convergence.vector <- rep(0,Niter.miss);
  iter.miss <- 1
  totalmiss <- mean(Index.miss)

  while(iter.miss <= Niter.miss&&convergence.vector[iter.miss]==0){
  ###################################################
  ######## Transform the Data           #############
  ###################################################
  Ytilde <- as.matrix(t(A0)%*%as.matrix(Bt%*%t(Y)))
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
    gcv <- gcv/(1-alpha*trace/J/(1-totalmiss))^2
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
    Eigen <- eigen(temp,symmetric=TRUE)
    A <- Eigen$vectors
    Sigma <- Eigen$values/J
   }

  if(c.p > I){
    temp <- t(YS)%*%YS/I
    Eigen <- eigen(temp,symmetric=TRUE)
    Sigma <- Eigen$values/J
    N <- sum(Sigma>0.0000001)
    A <- YS%*%(Eigen$vectors[,1:N]%*%diag(1/sqrt(Eigen$values[1:N])))/sqrt(I)
    }
  if(iter.miss>1&&iter.miss< Niter.miss) {
    diff <- norm(YS-YS.temp,"F")/norm(YS,"F");
    if(diff <= 0.02) 
      convergence.vector[iter.miss+1] <- 1
  }
  YS.temp <- YS
  iter.miss <- iter.miss + 1
  N <- min(I,c.p)
  d <- Sigma[1:N]
  d <- d[d>0]
  N <- length(d)
  per <- cumsum(d)/sum(d)
  N <- 1
  while(per[N]<pve) N <- N + 1;

  #print(c(iter.miss,convergence.vector[iter.miss+1],lambda,diff))
  #########################################
  #######     Principal  Scores   #########
  ########   data imputation      #########
  #########################################
 
  if(imputation==T){
    A.N <- A[,1:N]
    d <- Sigma[1:N]
    sigmahat2  <-  max(mean(Y[!Index.miss]^2) -sum(Sigma),0)
    
    Xi <- t(A.N)%*%Ytilde
    Xi <- t(as.matrix(B%*%(A0%*%((A.N%*%diag(d/(d+sigmahat2/J)))%*%Xi))))
    Y <- Y*(1-Index.miss) + Xi*Index.miss
    if(sum(is.na(Y))>0) print("error")
 
  }
  #if(iter.miss%%10==0) print(iter.miss)
   
  }## end of while loop         
  ### now calculate scores
  if(is.null(Y.pred)) Y.pred = Y
  else {Y.pred = t(t(as.matrix(Y.pred))-meanX)}


  Ytilde <- as.matrix(t(A0)%*%(Bt%*%t(Y.pred)))
  sigmahat2 <- max(mean(Y[!Index.miss]^2) -sum(Sigma),0)
  Xi <- t(Ytilde)%*%(A[,1:N]/sqrt(J))
  Xi <- MM(Xi,Sigma[1:N]/(Sigma[1:N] + sigmahat2/J))

  eigenvectors = as.matrix(B%*%(A0%*%A[,1:N]))
  eigenvalues = Sigma[1:N] #- sigmahat2/J
  
  Yhat <- t(A[,1:N])%*%Ytilde
  Yhat <- as.matrix(B%*%(A0%*%A[,1:N]%*%diag(eigenvalues/(eigenvalues+sigmahat2/J))%*%Yhat))
  Yhat <- t(Yhat + meanX)
  
  if(!is.null(npc)) {
    if(npc>N){
      cat("Warnig! The number of PCs is ", N,"s maller than ", npc,"\n");
      cat("Will use",N,"PCs\n");
      npc = N
    }
  }
  if(is.null(npc)) npc = N
  
  scores <- sqrt(J)*Xi[,1:npc]
  mu <- meanX
  efunctions <- eigenvectors[,1:npc]
  evalues <- J*eigenvalues[1:npc]
  
  
  ret.objects <- c("Yhat", "scores", "mu", "efunctions", "evalues", "npc")
  if(var){
    sigma2 = sigmahat2
    VarMats = vector("list",I)
    diag.var = matrix(NA, nrow=I,ncol=J)
    crit.val = rep(0,I)
    for(i.subj in 1:I){
      temp = sigma2*eigenvectors%*%solve(t(eigenvectors)%*%eigenvectors + sigma2*diag(eigenvalues))%*%t(eigenvectors)
      VarMats[[i.subj]] = temp
      diag.var[i.subj,] = diag(temp)
      if (simul & sigma2 != 0) {
        norm.samp = mvrnorm(2500, mu = rep(0, J), Sigma = VarMats[[i.subj]])/matrix(sqrt(diag(VarMats[[i.subj]])), 
                                                  nrow = 2500, ncol = J, byrow = TRUE)
        crit.val[i.subj] = quantile(apply(abs(norm.samp), 1, max), sim.alpha)
      }
    }
    ret.objects = c(ret.objects,"sigma2","diag.var","VarMats")
    if (simul) {
      #require(MASS)
      ret.objects = c(ret.objects, "crit.val")
    }
  }
  ret = lapply(1:length(ret.objects), function(u) get(ret.objects[u]))
  names(ret) = ret.objects
  return(ret)      	                	
}
