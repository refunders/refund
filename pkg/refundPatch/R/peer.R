peer<- function(Y, funcs, pentype='Ridge', L.user=NULL, Q=NULL, 
                a=10^3, se=FALSE, ...)
{
  require(nlme)
  require(magic)
  lobj=TRUE
  
  W<- as.matrix(funcs)
  norm<- sqrt(diag(W%*%t(W)))
  W<- W/norm
  K<- ncol(W) 
  
  Y<- as.matrix(Y)
  
  #Check 1:Making sure Y has only 1 column
  if(dim(Y)[2]>1) return(cat("Error: No. of column for Y cannot be greater than 1. \nThe peer() will not proceed further.\n"))
  
  Yl<- dim(Y)[1]
  
  #Check 3: Check the dimension of Y, id, t, W and X
  chk.eq<- ifelse(Yl==nrow(W), 0 ,1)
  if(chk.eq==1) return(cat("Error: Length of Y and number of rows of funcs are not equal.\n The peer() will not proceed further.\n"))
  
  #Removal of missing and infinite observations
  tdata<- data.frame(Y, W)
  tdata<- tdata[which(apply(is.na(tdata), 1, sum)==0),]
  tdata<- tdata[which(apply(is.infinite(as.matrix(tdata)), 1, sum)==0),]
  tdata<- tdata[!is.na(tdata$Y),]
  Y<- tdata$Y
  W.e<- dim(W)[2]+1; W<- as.matrix(tdata[,c(2:W.e)])
  
  N<- length(Y)
  
  #Check 4.1: Compatibility of Q and W matrix
  if(toupper(pentype)=='DECOMP' | toupper(pentype)=='DECOMPOSITION'){
    if(ncol(Q)!=ncol(W)) return(cat('Error: number of columns of Q need to be equal with number of columns of funcs.\nThe peer() will not proceed further.\n'))
    
    #Check 4.2: Removing rows containing missing and infinite values
    Q<- Q[which(apply(is.na(Q), 1, sum)==0),]
    Q<- Q[which(apply(is.infinite(Q), 1, sum)==0),]
    Q<- matrix(Q, ncol=K)
    
    #Check 4.3: Singularity of Q matrix
    Q.eig<- abs(eigen(Q %*% t(Q))$values)
    if(any(Q.eig<1e-12)) return(cat('Error: Q matrix is singular or near singular.\nThe peer() will not proceed further.\n'))
  }
  
  #Check 5.1: Dimension of L matrix
  if(toupper(pentype)=='USER'){
    if(ncol(L)!=ncol(W)) return(cat('Error: number of columns of L need to be equal with number of columns of funcs.\nThe peer() will not proceed further.\n'))
    
    #Check 5.2: Removing rows containing missing and infinite values
    L<- L[which(apply(is.na(L), 1, sum)==0),]
    L<- L[which(apply(is.infinite(L), 1, sum)==0),]
    L<- matrix(L, ncol=K)
    
    #Check 5.3: Singularity of L'L matrix
    LL<- t(L)%*%L
    LL.eig<- abs(eigen(LL %*% t(LL))$values)
    if(any(LL.eig<1e-12)) return(cat("Error: L'L matrix is singular or near singular.\nThe peer() will not proceed further.\n"))
  }
  
  #Generate L matrix for D2 penalty
  if(toupper(pentype)=='D2'){
    Left<- cbind(diag(rep(1,K-2)),rep(0,K-2),rep(0,K-2))
    Middle<- cbind(rep(0,K-2),diag(rep(-2,K-2)),rep(0,K-2))
    Right<- cbind(rep(0,K-2),rep(0,K-2),diag(rep(1,K-2)))
    D.2<- rbind(Left+Middle+Right, c(rep(0, K-2), 1, -2), c(rep(0, K-2), 0, 1))
  }
  
  #Generate W* matrix
  if(toupper(pentype)=='DECOMP' | toupper(pentype)=='DECOMPOSITION'){
    P_Q <- t(Q) %*% solve(Q %*% t(Q)) %*% Q
    L_PEER<- a*(diag(K)- P_Q) + 1*P_Q
  } else
    if(toupper(pentype)=='RIDGE'){
      L_PEER<- diag(K)
    } else
      if(toupper(pentype)=='D2'){
        L_PEER<- D.2
      } else
        if(toupper(pentype)=='USER'){
          L_PEER<- L.user
        } 
  
  v<- diag(K)
  if(K>N) v<-  svd((data.matrix(W))%*% solve(L_PEER))$v
  W1_PEER<- (data.matrix(W))%*% solve(L_PEER) %*% v
  
  #Fitting the model
  id.bd1<- factor(rep(1, length(Y)))
  out_PEER<- lme(fixed=Y~1,
                 random=list(id.bd1=pdIdent(~W1_PEER-1)),
                 ...
  ) 
  
  
  cat('The fit is successful.\n')
  
  #Estimate
  Gamma.PEER.hat<-matrix(out_PEER$coeff$random$id.bd1, ncol=1)
  
  GammaHat <- solve(L_PEER) %*% v %*%Gamma.PEER.hat
  
  colnames(GammaHat)<- c('Gamma')
  
  fitted.vals<- summary(out_PEER)$tTable
  logLik<- summary(out_PEER)$logLik
  AIC<- summary(out_PEER)$AIC
  BIC<- summary(out_PEER)$BIC
  
  tVarCorr<- nlme:::VarCorr(out_PEER, rdig=4)[,2]
  r<- ncol(v)
  lambda<- 1/ as.numeric(unique(tVarCorr[1:r]))
  
  print(lambda)
  
  sigma<- out_PEER$sigma
  
  ###---- Standard Error
  sigma.e<- sigma
  
  tsigma<- as.numeric(unique(tVarCorr[1:K]))
  LL.inv<- lambda^(-2)*solve(t(L_PEER)%*%L_PEER) 
  
  if(!se){
    status<- 0
    ret <- list(out_PEER, fitted.vals, GammaHat,  lobj,
                AIC, BIC, logLik,
                lambda, N, K, sigma.e, status)
    names(ret)<- c("fit", "fitted.vals", "GammaHat", "peerobj",
                   "AIC", "BIC", "logLik", 
                   "lambda", "N", "K", "sigma", "status")
    class(ret) <- "peer"
    return(ret)
  }
  
  V.1<- W%*%LL.inv%*%t(W)+sigma.e^2*diag(rep(1, N))
  V<- sigma.e^2*diag(rep(1, N))
  
  X<- as.matrix(rep(1,N))
  X.Vinv.X.inv<- solve(t(X)%*%solve(V.1)%*%X)
  X.Vinv.Y<- t(X)%*%solve(V.1)%*%Y
  
  p1<- LL.inv%*%t(W)%*%solve(V.1)
  p2<- V.1 - X%*%X.Vinv.X.inv%*%t(X)
  p3<- solve(V.1)
  p<- p1%*%p2%*%p3
  
  se.Gamma<- sqrt(diag(p%*%V%*%t(p)))
  Gamma<- cbind(GammaHat, se.Gamma)
  colnames(Gamma)<- c('Estimate', 'SE')
  
  status<- 1
  ret <- list(out_PEER, fitted.vals, Gamma,  lobj,
              GammaHat, se.Gamma, AIC, BIC, logLik,
              lambda, N, K, sigma.e, status)
  names(ret)<- c("fit", "fitted.vals", "Gamma", "peerobj",
                 "GammaHat", "se.Gamma", "AIC", "BIC", "logLik", 
                 "lambda", "N", "K", "sigma", "status")
  class(ret) <- "peer"
  ret
}







