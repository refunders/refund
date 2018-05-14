raw.construct <- function(data,include.diag=TRUE){
  
  y <- data$y
  t <- data$argvals
  subj <- data$subj
  
  subj_unique <- unique(subj)
  n <- length(subj_unique)
  C <- c()
  st <- matrix(NA,ncol=2,nrow=0)
  N <- c()
  N2 <- c()
  n0 <- 0
  W <- list(length=n)
  for(i in 1:n){
    
    r1 <- y[subj==subj_unique[i]]
    t1 <- t[subj==subj_unique[i]]
    m1 <- length(t1)
    n0 <- n0 + 1
    if(m1>1){
      if(include.diag) {
        N2 <-c(N2,m1*(m1+1)/2)  # <------
        sel = 1:N2[n0]
      }
      
      if(!include.diag) {
        N2 <-c(N2,m1*(m1-1)/2)  # <------
        sel = setdiff(1:(m1*(m1+1)/2), c(1,1 + cumsum(m1:1)[1:(m1-1)]))
      }
      
      st <- rbind(st,cbind(vech(kronecker(t1,t(rep(1,m1)))),
                           vech(kronecker(rep(1,m1),t(t1))))[sel,])
      C <- c(C,vech(kronecker(r1,t(r1)))[sel]) 
      

      N <- c(N,m1)
      # N2 <-c(N2,m1^2)

      W[[i]] <- sparseMatrix(1:N2[n0],1:N2[n0],x=rep(1,N2[n0]))# <----
      #if(include.diag) diag(W[[i]])[c(1,1 + cumsum(m1:1)[1:(m1-1)])] <- 1/2
    }## for if(m1>1)
    if(m1==1){
      if(include.diag){
      N2 <- c(N2,1)
      st <- rbind(st,c(t1,t1))
      C <- c(C,r1^2)
      N <- c(N,1)
      W[[i]] <- matrix(1,1,1)
      }
      if(!include.diag){
        N2 <- c(N2,0)
        N <- c(N,1)
        W[[i]] <- NULL
      }
    }
  }##for i
  
  res <- list("C" = C,
              "st" = st,
              "N" = N,
              "N2" = N2,
              "W" = W,
              "n0" = n0)
  return(res)
}

