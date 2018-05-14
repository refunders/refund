weight.construct <- function(fit,data,include.diag=TRUE){

 
  r <- fit$r 
  t <- data$argvals
  subj <- data$subj
  
  subj_unique <- unique(subj)
  n <- length(subj_unique)
  W <- list(length=n)
  
  Theta <- fit$Theta
  sigma2 <- fit$sigma2
  knots <- fit$knots
  p <- fit$p
  
  B <- spline.des(knots=knots, x=t, ord = p+1,
                  outer.ok = TRUE,sparse=TRUE)$design
  
  for(i in 1:n){

  r1 <- r[subj==subj_unique[i]]
  t1 <- t[subj==subj_unique[i]] # no sort
  m1 <- length(t1) # m_i
  B1 <- matrix(B[subj==subj_unique[i],],nrow=sum(subj==subj_unique[i]))
  Chati <- as.matrix(tcrossprod(B1%*%Matrix(Theta),B1))
  if(m1>1){
       ##still input r w/o any deletion, CovZ.new deletes diagonal entires
    re <- covZ(t1, r1, Chati, sigma2) # <------ 
    if(include.diag){
    sel <- c()
    for(j in 1:m1) sel <- c(sel,(j-1)*m1 + j:m1)
    }
    if(!include.diag){
      sel <- c()
      for(j in 1:(m1-1)) sel <- c(sel,(j-1)*m1 + (j+1):m1)
    }
  
    CovZZ.tri <- re$CovZZ[sel,sel]
    if(length(sel)>1){
    V1 <- forceSymmetric(CovZZ.tri)
    V1 <- 0.95*V1 + 0.05*diag(diag(V1))
    # avoid inv 
    eSig <- eigen(V1)
    eV = eSig$vectors
    eE = eSig$values
    eE[eE<0] <- 0
    eE = eE + 0.000001*max(eE)
    W[[i]] = eV %*% diag(1/eE) %*% t(eV)
    }
    if(length(sel)==1){
      W[[i]] = 1/CovZZ.tri
    }
    }## for if(m1>1) 
  if(m1==1){
    if(include.diag){
      W[[i]]=matrix(1/(Chati + sigma2)^2/2,1,1)
    }
    if(!include.diag){W[[i]]=NULL}
  }
}##for i

#scale W 
temp <- sapply(W,max)
W <- lapply(W,function(x) x/mean(temp,na.rm=TRUE))

return(W)
}
