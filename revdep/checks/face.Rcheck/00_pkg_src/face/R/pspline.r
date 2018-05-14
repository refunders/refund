pspline <- function(data,argvals.new= NULL,knots=35,
                        # knots.option="equally-spaced",
                        p=3,m=2,lambda=NULL,
                        search.length = 100,
                        lower=-20,upper=20){
  
  
  #require(splines)  
  #source("pspline.setting.weight.R")
  
  ## read in data
  check.data(data)
  Y <- data$y
  t <- data$argvals
  subj <- data$subj
  tnew <- argvals.new
  if(is.null(tnew)){tnew <- t}
  ###
  J <- length(Y)
  if(is.null(t))  t <- (1:J)/J-1/2/J ## if NULL, assume equally spaced
  subj_unique <- unique(subj)
  n <- length(subj_unique)
  WY <- Y
  N <- rep(NA,n)
  W <- 0*Y
  for(i in 1:n){
    seq <- (1:J)[subj==subj_unique[i]]
    N[i] <- length(seq)
    WY[seq] <- Y[seq]/N[i]
    W[seq] <- 1/N[i]
  }
  
  p.p <- p
  m.p <- m

  ######### precalculation for smoothing ############
  knots.option="equally-spaced"
  knots <- construct.knots(t,knots,knots.option,p)
  
  
  List <- pspline.setting(t,knots=knots,p.p,m.p,weight=W,knots.option=knots.option)
  AS <- as.matrix(List$A)
  s <- List$s
  Sigi.sqrt <- List$Sigi.sqrt
  U <- List$U
  B <- List$B
  Bnew <- spline.des(knots=knots, x=tnew, ord = p.p+1,outer.ok = TRUE,sparse=TRUE)$design
    
  AStY <- as.vector(t(AS)%*%Y)
  AStWY <- as.vector(t(AS)%*%WY)
  AStAS <- as.matrix(t(AS)%*%AS)
  AStAS_eigen <- eigen(AStAS)
  AStAS_eigen$values[AStAS_eigen$values<=0] <- 0
  AStAS_half <- AStAS_eigen$vectors%*%diag(sqrt(AStAS_eigen$values))%*%t(AStAS_eigen$vectors)

  ##calculate Term1
  AStAS_N <- matrix(NA,n,length(s)^2)
  ASitYi <- matrix(NA,n,length(s))
  Term1 <- matrix(0,length(s),length(s))
  Term0 <- rep(0,length(s))
  for(i in 1:n){
    seq <- (1:J)[subj==subj_unique[i]]
    ASi <- matrix(AS[seq,],N[i])
    ASitASi <- t(ASi)%*%ASi
    AStAS_N[i,] <- as.vector(ASitASi)
    temp <- as.vector(t(ASi)%*%Y[seq]/N[i])
    ASitYi[i,] <- temp
    Term1 <- Term1 + diag(temp)%*%as.matrix(ASitASi%*%diag(AStWY))
    Term0 <- Term0 + temp^2*N[i] 
  }
  

  pspline_gcv <- function(x){
    lambda <- exp(x)
    lambda_s <- 1/(1 + lambda*s)
    gcv <- sum((AStAS_half%*%(lambda_s*AStWY))^2)- 2*sum(lambda_s*(AStY*AStWY))
    gcv <- gcv + 2*sum(lambda_s*Term0)
    gcv <- gcv - 4*sum(lambda_s*(Term1%*%lambda_s))
    for(i in 1:n){
      gcv <- gcv + 2*sum((sqrt(lambda_s)*(matrix(AStAS_N[i,],length(s),length(s))%*%(lambda_s*AStWY)))^2)/N[i]
    }
    return(gcv) 
  }
  
  if(is.null(lambda)){
      Lambda <- seq(lower,upper,length=search.length)
      Length <- length(Lambda)
      Gcv <- rep(0,Length)
      for(i in 1:Length) 
        Gcv[i] <- pspline_gcv(Lambda[i])
      i0 <- which.min(Gcv)
      lambda <- exp(Lambda[i0])
  }
  
  theta <- (Sigi.sqrt %*% U)%*%(1/(1+lambda*s)*AStWY)

  res <- list(fitted.values = as.vector(B%*%theta), B = B,theta=theta,s = s,
              knots=knots,p=p,m=m,
              lambda=lambda,argvals.new = tnew, mu.new = as.vector(Bnew%*%theta))
  class(res) <- "pspline.face"
  return(res)
  
}




