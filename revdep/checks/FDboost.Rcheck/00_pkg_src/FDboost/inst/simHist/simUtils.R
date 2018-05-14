###############################################################################
# Utility functions for functional regression with pffr and FDboost 
# code based on code of Fabian Scheipl: utility functions for lfpr simulations 
# Author: Sarah Brockhaus
###############################################################################


#################################
# Function of Fabian "superUtils.R"
expandList <- function(...) {
  ## expand.grid for lists
  dots <- list(...)
  #how many settings per entry
  dims <- lapply(dots, length)
  #make all combinations of settings
  sets <- do.call(expand.grid, lapply(dims, function(x) seq(1:x)))
  ret <- apply(sets, 1, function(x) {
    l <- list()
    for (i in 1:length(x)) l[[i]] <- dots[[i]][[as.numeric(x[i])]]
    names(l) <- colnames(sets)
    return(l)
  })
  for (c in 1:ncol(sets)) sets[, c] <- factor(sets[, c], labels = paste(dots[[colnames(sets)[c]]]))
  attr(ret, "settings") <- sets
  return(ret)
}


makeSettings <- function(dgpsettings, algorithms){
  # generates a full factorial design of all settings 
  #
  # dgpsettings:  a named list of parameters for the data generating process
  # algorithms: a named list of parameters for the algorithms/estimators
  #
  # sets up the simulation study so that algorithms are compared on the same data
  # for each combination of dgpsettings by re-using the same seed for all
  # algorithms i.e., replication <x> of a "datasetting" will produce the same data
  # regardless of the "algorithm".
  
  datasettings <- do.call(expand.grid, c(dgpsettings, stringsAsFactors =
                                           FALSE))
  datasettings$seed <- as.integer(sample(1L:5e7L, nrow(datasettings)))
  datasettings$datasetting <- 1:nrow(datasettings)
  
  algsettings <- do.call(expand.grid, c(algorithms, stringsAsFactors =
                                          FALSE))
  algsettings$algorithm <- 1:nrow(algsettings)
  
  settings <- cbind(datasettings[rep(1:nrow(datasettings),
                                     each  = nrow(algsettings)), ],
                    algsettings[rep(1:nrow(algsettings),
                                    times = nrow(datasettings)), ],
                    set=1:(nrow(datasettings)*nrow(algsettings)), 
                    combination=rep(1:(nrow(algsettings)*nrow(datasettings)/max(datasettings$rep)), 
                                    times= max(datasettings$rep)) )
  
  settingsList <- lapply(1:nrow(settings), function(j) sapply(1:ncol(settings), 
                                                              function(i) settings[j,][i]))
  
  attr(settingsList, "settings") <- settings
  
  return(settingsList)
}

# generate colors
alpha <- function(x, alpha=25){
  tmp <- sapply(x, col2rgb)
  tmp <- rbind(tmp, rep(alpha, length(x)))/255
  return(apply(tmp, 2, function(x) do.call(rgb, as.list(x))))
} 

#################################
# Do the iterations over oneRep()
doSim <- function(settings, cores=40){
    library(plyr)
    
    # only works on Linux -> with try() no error on windows
    try(library(doMC))
    try(registerDoMC(cores=cores))
    
    split <- sample(rep(1:cores, length=length(settings)))
    
    settingsSplit <- alply(1:cores, 1, function(i) {
                settings[split==i]   
            })
    
    ret <- ldply(settingsSplit, function(settings) {
                return(ldply(settings, oneRep, .parallel=FALSE))
             }, .parallel=TRUE)
    
    return(ret)    
}
# 
# # Do the iterations over oneRep2()
# doSim2 <- function(settings, cores=40){
#   library(plyr)
#   
#   # only works on Linux -> with try() no error on windows
#   try(library(doMC))
#   try(registerDoMC(cores=cores))
#   
#   split <- sample(rep(1:cores, length=length(settings)))
#   
#   settingsSplit <- alply(1:cores, 1, function(i) {
#     settings[split==i]   
#   })
#   
#   ret <- ldply(settingsSplit, function(settings) {
#     return(ldply(settings, oneRep2, .parallel=FALSE))
#   }, .parallel=TRUE)
#   
#   return(ret)    
# }
# 

# Do the iterations over oneRepPffr()
doSimPffr <- function(settings, cores=10){
  library(plyr)
  
  # only works on Linux -> with try() no error on windows
  try(library(doMC))
  try(registerDoMC(cores=cores))
  
  split <- sample(rep(1:cores, length=length(settings)))
  
  settingsSplit <- alply(1:cores, 1, function(i) {
    settings[split==i]   
  })
  
  ret <- ldply(settingsSplit, function(settings) {
    return(ldply(settings, oneRepPffr, .parallel=FALSE))
  }, .parallel=TRUE)
  
  return(ret)    
}

# Do the iterations over oneRepFDboost() 
# slightly modified doSafeSim() 
doSimFDboost <- function(settings){  #, savefile
  library(plyr)
  
  ret <- data.frame()
  failed <- list()
  
  for(s in 1:length(settings)){
    res <- try(do.call(oneRepFDboost, settings[s]), silent = TRUE)
    if(any(class(res)=="try-error")){
      cat("\n some of ", s, "failed:\n")
      print(do.call(rbind, settings[s]))
      failed <- c(failed, settings[s])
    }else{     
      ret <- rbind(ret, res)
      save(ret, file="savefile.Rdata")
      save(ret, file=glue("cpy", "savefile.Rdata"))
      cat(format(Sys.time(), "%b %d %X"), ": ", s, " of ", length(settings), "\n")
    }
  }
  #ret <- list(ret=ret, failed=failed)
  attr(ret, "failed") <- failed
  #save(ret, file=savefile)
  #save(ret, file=glue("cpy", savefile))
  return(ret)    
}

## Try function ldply()
#test.list <- list(set1=1:4, set2=5:8, set3=9:12)
#test.fun <- function(l) data.frame(l, x=rep(99,4))
#test.fun(test.list[[1]])
#ldply(test.list, test.fun


glue <- function(..., collapse = NULL) {
  paste(..., sep = "", collapse)
}
# 
# doSafeSim <- function(settings, savefile){
#   library(plyr)
#   
#   ret <- data.frame()
#   failed <- list()
#   
#   for(s in 1:length(settings)){
#     res <- try(do.call(oneRep, settings[s]), silent = TRUE)
#     if(any(class(res)=="try-error")){
#       cat("\n some of ", s, "failed:\n")
#       print(do.call(rbind, settings[s]))
#       failed <- c(failed, settings[s])
#     }else{     
#       ret <- rbind(ret, res)
#       save(ret, file=savefile)
#       save(ret, file=glue("cpy", savefile))
#       cat(format(Sys.time(), "%b %d %X"), ": ", s, " of ", length(settings), "\n")
#     }
#   }
#   ret <- list(ret=ret, failed=failed)
#   save(ret, file=savefile)
#   save(ret, file=glue("cpy", savefile))
#   return(ret)    
# }



#################################
# Funcitons for generating variables and coefficients

# function for generating a functional covariable as (spline-basis)%*%(random coefficient)
rf <- function(x=seq(0,1,length=100), k=15) {
  
  # parallel increasing lines with random slope 
  if(k==0) ret <- rnorm(1, mean=-2, sd=1) + 0.3*x
  
  # lines with random slopes
  if(k==1) ret <- rnorm(1, mean=0, sd=0.1)*x
  
  # lines with random intercepts and random slopes 
  if(k==2) ret <- rnorm(1, mean=-2, sd=1) + rnorm(1, mean=0, sd=0.1)*x
  
  # data simulated by splines with random coefficients
  if(k>2) ret <- drop(bs(x, k, int=TRUE) %*% runif(k, -3, 3))
  
  return(ret)
}



# function for generating a functional covariable with special properties
rf2 <- function(x=seq(0,1,length=100), k=15, type="bsplines") {
  
  if(type=="lines"){
    ret <- rf(x=x, k=k)
  }
  
  if(type=="bsplines"){
    ret <- rf(x=x, k=k)
  }
  
  if(type == "local"){
    ## locale but with different random bases
    temp <- c( min(x)-0.05*(range(x)[2] - range(x)[1]), max(x)+0.05*(range(x)[2] - range(x)[1]) )
    ret <- bs(x, 5*k, intercept=TRUE, Boundary.knots=temp)[,sample(2:(5*k-1), k)] %*% c(runif(k, -3, 3))
    # ret<- cbind(ret, bs(x, 5*k, intercept=TRUE)[,sample(2:(5*k-1), k)] %*% c(runif(k, -3, 3)))
    # funplot(x, t(ret))
  }
  
  ### all information is at the same locations
  #   if(type == "locale0"){
  #     ## local and all functions have the same local bases
  #     ret <- bs(x, 5*k, intercept=TRUE)[,round(seq(5, 5*k-5, l=k))] %*% c(runif(k, -3, 3))
  #     # ret <- cbind(ret, bs(x, 5*k, intercept=TRUE)[,round(seq(5, 5*k-5, l=k))] %*% c(runif(k-2, -3, 3), runif(2, -2, 2)))
  #   }
   
  ### information in the beginning
  if(type == "start"){  
    ## information in the beginning
    # ret <- cbind(bs(x, 10*k, intercept=TRUE)[,sample(2:(3*k), k-1)], 1) %*% c(runif(k-1, -3, 2), runif(1, -1, 2))
    ret <- bs(x, 10*k, intercept=TRUE)[,sample(2:(3*k), k)] %*% c(runif(k, -3, 3))
  }
  
  ### some information in the end
  if(type == "end"){
    ## only information in the end
    #ret <- cbind(bs(x, 10*k, intercept=TRUE)[,sample((7*k):(10*k-1), k-1)], 1) %*% c(runif(k-1, -3, 2), runif(1, -1, 2))
    ret <- bs(x, 10*k, intercept=TRUE)[ ,sample((7*k):(10*k-1), k)] %*% c(runif(k, -3, 3))
  }
  
  ## use a fourier basis with linear decreasing eigenvalues
  if(type == "fourierLin"){
      EFourier <- eval.basis(x, create.fourier.basis(rangeval=c(0, 1), 
                                                     nbasis=ifelse(k%%2, k, k+1)))
      loadings <- replicate(1, rnorm(k, sd=sqrt(((k+1)-(1:k))/k)) )
      ret <- EFourier[,1:k] %*% loadings
  }  
  
  ## use a fourier basis with exponentailly decreasing eigenvalues
  if(type == "fourierExp"){
    EFourier <- eval.basis(x, create.fourier.basis(rangeval=c(0, 1), 
                                                   nbasis=ifelse(k%%2, k, k+1)))
    loadings <- replicate(1, rnorm(k, sd=sqrt(exp(-(0:(k-1))/2)) ))
    ret <- EFourier[,1:k] %*% loadings
  } 
  
  ## use a fourier basis with constant eigenvalues
  if(type == "fourier"){
    #ret <- cos(runif(1, 1, 2)*pi*x)*runif(1, 1, 2) #*seq(1, runif(1,2,5), l=length(x))
    #ret <- cbind(ret, cos(runif(1, 1, 2)*pi*x)*runif(1, 1, 2))
    #funplot(x, t(ret))
    EFourier <- eval.basis(x, create.fourier.basis(rangeval=c(0, 1), 
                                                   nbasis=ifelse(k%%2, k, k+1)))
    loadings <- replicate(1, rnorm(k, sd=1))
    ret <- EFourier[,1:k] %*% loadings
  } 
  
  return(drop(ret))
}

types <- c("bsplines","locale","locale0","start","start0","end","end0","fourier")


### plot examples for all datasettings
if(FALSE){
  library(FDboost)
  library(splines)
  s <- seq(0, 1, l=100)*15 + 1
  generateX <- function(k, type){
    centerX <- TRUE
    X1 <- t(replicate(10, rf2(x=s, k=k, type=type)))
    if(centerX) X1 <- sweep(X1, 2, apply(X1, 2, mean))
    funplot(s, X1, type="l", rug = FALSE, 
            lwd=2, col=1, lty=1, xlab="s", ylab="", 
            main=paste(type, "-", k))
    return(X1)
  }
  pdf("dataSim.pdf")
  par(mar=c(3.5, 3, 2.5, 1), cex=2, mgp = c(2, 1, 0), cex.main=1.5)
  set.seed(124)
  X1 <- generateX(k=5, type="local")
  X1 <- generateX(k=10, type="local")
  X1 <- generateX(k=5, type="bsplines")
  X1 <- generateX(k=10, type="bsplines")
  X1 <- generateX(k=5, type="end")
  X1 <- generateX(k=10, type="end")
  X1 <- generateX(k=0, type="lines")
  X1 <- generateX(k=1, type="lines")
  X1 <- generateX(k=2, type="lines")
  frame()
  dev.off()
}

#  
# # rf1 <- function(x=seq(0,1,length=100)) { 
# #     m <- ceiling(runif(1)*5) ## number of components 
# #     f <- x*0 + rnorm(1); 
# #     mu <- runif(m,min(x),max(x))
# #     sig <- (runif(m)+.5)*(max(x)-min(x))/10 
# #     for (i in 1:m) f <-  f + dnorm(x,mu[i],sig[i]) 
# #     f 
# # } 
# # # matlplot(seq(0,1,length=100), (replicate(50, rf1())))
# # 
# # rf2 <- function(x=seq(0,1,length=100)) { 
# #     m <- sample(2:4, 1) ## number of components 
# #     f <- x*0; 
# #     for (i in 1:m) f <- f + rnorm(1, sd=1/i) * sin(i/2*pi*x) + rnorm(1, sd=1/i) * cos(i/2*pi*x) 
# #     f 
# # } 
# # # matlplot(seq(0,1,length=100), (replicate(50, rf2())))
# 
# rf3 <- function(x=seq(0,1,length=100)) {
#     rnorm(1)/2 + rnorm(1)*(x-.5) +  rnorm(1)*(x-.2)^2 + rnorm(1)*(x-.8)^3
# } 
# # matlplot(seq(0,1,length=100), (replicate(50, rf3())))
# 
# rf4 <- function(x=seq(0,1,length=100)) {
#     drop(bs(x, 5, int=TRUE) %*% ((-2:2)/2+rnorm(5)))
# } 
# #matlplot(seq(0,1,length=100), (replicate(50, rf4())))
#  
# # rf5 <- function(x=seq(0,1,length=100)) {
# #     m <- ceiling(runif(1)*5)+1
# #     w <- rnorm(m, sd=1/(1:m))
# #     drop(poly(x, m+1)[,-1]%*%w)
# # } 
# # matlplot(seq(0,1,length=100), scale(replicate(50, rf5())))

# # standardize matrix, so that colSums are zero
# zeroConstraint <- function(x){
#   stopifnot(is.matrix(x))
#   t(t(x) - colMeans(x))
# }

# # generate zero centered scalar covariates
# cenScalarCof <- function(n){
#   z <- runif(n) - 0.5
#   z <- z - mean(z)   
#   z
# }

## built a regular grid over the range of a scalar covariable
#zgrid <- function (z, l=40) seq(min(z), max(z), l=l)

# # generate a smooth global intercept
# # global intercept
# intf <- function(t){
#   0.5*cos(3*pi*t^2) 
# }

# generate smooth intercept
intf1 <- function(t){ 
 # 1 + log(t+0.5)
  1 + 2*sqrt(t)
}


# ## beta(s,t)
# (from Sonja's example)
test1 <- function(s, t){
  
  stopifnot(length(s)==length(t))
  ret <- 1/2* sin(s*2*t*pi)*log(1+s+t) + s*t*2 + exp(s)*t^2
  ret[s>t] <- 0 

  return(ret)
}
# persp(seq(0,1, l=20), seq(0,1, l=20), outer(seq(0,1, l=20),seq(0,1, l=20), test1)) 
# image(seq(0,1, l=31), seq(0,1, l=20), outer(seq(0,1, l=31), seq(0,1, l=20), test1)) 


# (from ?gam example)   
test2 <- function(s, t, ss=0.3, st=0.4){ 
  
  stopifnot(length(s)==length(t))
#   ret <- 4*((pi^ss*st)*(1.2*exp(-(s-0.2)^2/ss^2-(t-0.3)^2/st^2) +
#                 0.8*exp(-(s-0.7)^2/ss^2-(t-0.8)^2/st^2)))
  ret <- 1.5*sin(pi*t+0.3) * sin(pi*s)
  ret[s>t] <- 0 
  
  return(ret)
}

# persp(seq(0,1, l=20), seq(0,1, l=20), outer(seq(0,1, l=20),seq(0,1, l=20), test2))


# # Harezlak etAl. 2007
# test3 <- function(s, t, a1=3, a2=1, a3=1, a4=1){ 
#   stopifnot(length(s)==length(t))
#   
#   ret <- a4*sin(a1*t - a2*s) + a3 # Harezlak etAl. 2007
#   
#   ret[s>t] <- 0 
#   
#   return(ret)
# }
# 
# # persp(seq(0,1, l=20), seq(0,1, l=20), outer(seq(0,1, l=20),seq(0,1, l=20), test3))
# # persp(seq(0,1, l=20), seq(0,1, l=20), outer(seq(0,1, l=20),seq(0,1, l=20), test3, 6, 3, 0, 10), ticktype="detailed")
# 
# # g2zt <- function(z, t){ 2*(-t^2-0.1)*sin(pi*z+0.5) } 
# # #g2zt <- function(z, t){ 2*(-t^2-0.1)*cos(pi*z + pi/4) }


# similar to Harezlak etAl. 2007 with reparametrization to 1, ..., 16
test3 <- function(s, t, a1=3, a2=1, a3=1, a4=1){ 
  stopifnot(length(s)==length(t))
  
  s <- s/15-1
  t <- t/15-1
  
  ret <- a4*sin(a1*t + a2*s) + a3 + t + s 
  
  ret[s>t] <- 0 
  
  return(ret)
}

### function to set lower triangular to 0 or NA
lowerTo <- function(x, repl=0){
  stopifnot(ncol(x)==nrow(x))
  #x*outer(1:ncol(x), 1:nrow(x), "<=") # gives the same if repl=0
  x[ outer(1:ncol(x), 1:nrow(x), "<=")==FALSE] <- repl
  x
}

#persp(1:16, 1:16, lowerTo(outer(1:16, 1:16, test3, 2, 2, 1, 3), NA), ticktype="detailed", zlab="", theta=30, phi=30)
#persp(1:16, 1:16, lowerTo(outer(1:16, 1:16, test3, 1, 1, 1, 3), NA), ticktype="detailed", zlab="", theta=30, phi=30)
#persp(1:16, 1:16, lowerTo(outer(1:16, 1:16, test3, 0, 3, 1, 3), NA), ticktype="detailed", zlab="", theta=30, phi=30)
#persp(1:16, 1:16, lowerTo(outer(1:16, 1:16, test3, 3, 0, 1, 3), NA), ticktype="detailed", zlab="", theta=30, phi=30)



## function written by Fabian Scheipl for functional effect
## <SB> changed coefficient surface to historical effect
randomcoef <- function(s, t, coef=NULL, seed=NULL, df=5, pen=c(1,1), lambda=c(1,1)){
  if(!is.null(seed)) set.seed(seed)
  require(splines)
  Bs <- bs(s, df=df, intercept = TRUE)
  Bt <- bs(t, df=df, intercept = TRUE)
  
  # Recursion for difference operator matrix
  makeDiffOp <- function(degree, dim){
    if(degree==0){
      return(diag(dim))  
    } else {
      return(diff(makeDiffOp(degree-1, dim)))
    }    
  }
  Pt <- lambda[1] * kronecker(crossprod(makeDiffOp(pen[1], df)), diag(df))
  Ps <- lambda[2] * kronecker(diag(df), crossprod(makeDiffOp(pen[2], df)))            
  P <- .1*diag(df^2) + Pt + Ps
  
  if(is.null(coef)){
    coef <- matrix(solve(P, rnorm(df^2)), df, df) 
  }
  
  
  ret <- Bs%*%coef%*%t(Bt)
  
  rownames(ret) <- round(s, 2)
  colnames(ret) <- round(t, 2)
  
  for(i in 1:length(s)){
    for(j in 1:length(t)){
      if(s[i] > t[j]){
        ret[i, j] <- 0
      } 
    }
  }
  
  attr(ret, "coef") <- coef
    
  return(ret)
  
}


## coefficient functions with FIRST order differences
pen1coef4 <- function(s, t, coef=NULL) randomcoef(s,t, coef=coef, df=4, lambda=c(1,1), pen=c(1,1))
#pen.1coef4 <- function(s,t, coef=NULL) randomcoef(s,t, coef=coef, df=4, lambda=c(.1,.1))
#pen1coef6 <- function(s, t, coef=NULL) randomcoef(s, t, coef=coef, df=6)
#pen.1coef6 <- function(s,t, coef=NULL) randomcoef(s,t, coef=coef, df=6, lambda=c(.1,.1))


## coefficient functions with SECOND order differences
#pen1coef4s <- function(s, t, coef=NULL) randomcoef(s,t, coef=coef, df=4, pen=c(2,2))
pen2coef4 <- function(s,t, coef=NULL) randomcoef(s,t, coef=coef, df=4, lambda=c(1,1), pen=c(2,2))

#pen1coef2 <- function(s, t, coef=NULL) randomcoef(s,t, coef=coef, df=2)

if(FALSE){
  s <- seq(0, 1, l=25)
  t <- seq(0, 1, l=25)
  test <- pen1coef6(s, t, coef=NULL)
  #persp(s, t, test, ticktype="detailed")
  #test2 <- pen1coef6(sgrid, tgrid, coef=attr(test, "coef"))
  #persp(sgrid, tgrid, test2, ticktype="detailed")

  par(mfrow=c(1,2))
  persp(s, t, lowerTo(pen1coef4(s,t), NA), zlab="", theta=30, phi=30, ticktype="detailed", main="pen1coef4")
  persp(s, t, lowerTo(pen2coef4(s,t), NA), zlab="", theta=30, phi=30, ticktype="detailed", main="pen.1coef4")
  
  
}



###############################
dlv1 <- function(A, B, tol=1e-10){
  ## A, B orthnormal!!
  #Rolf Larsson, Mattias Villani (2001) 
  #"A distance measure between cointegration spaces"
  if(NCOL(A)==0 | NCOL(B)==0){
    return(1.0)  
  } 
  
  if(NROW(A) != NROW(B) | NCOL(A) > NROW(A) | NCOL(B) > NROW(B)){
    return(NA)
  }
  
  if(NCOL(B)<=NCOL(A)){
    Aorth <- MASS::Null(A)
    dist <- sum(diag(t(B) %*% Aorth %*% t(Aorth) %*% B)) / 
      min(NCOL(B), NROW(B) - NCOL(B))
  } else {
    Borth <- MASS::Null(B)
    dist <- sum(diag(t(A) %*% Borth %*% t(Borth) %*% A)) / 
      min(NCOL(A), NROW(A)-NCOL(A))
  } 
  return(dist)
}

dlv2 <- function(A, B, tol=1e-10){
  ## A, B orthnormal!!
  
  #Rolf Larsson, Mattias Villani (2001) 
  #"A distance measure between cointegration spaces"
  
  
  if(NCOL(A)==0 | NCOL(B)==0){
    return(1.0)  
  } 
  
  if(NROW(A) != NROW(B) | NCOL(A) > NROW(A) | NCOL(B) > NROW(B)){
    return(NA)
  }
  
  trace <- if(NCOL(B)<=NCOL(A)){
    sum(diag(t(B) %*% A %*% t(A) %*% B))
  } else {
    sum(diag(t(A) %*% B %*% t(B) %*% A)) 
  }
  
  dist <-  (min(NCOL(B), NCOL(A)) - trace) 
  
  dist / min(NCOL(A), NCOL(B), 
             NROW(A) - NCOL(A), NROW(B) - NCOL(B))
}

## measure degree of overlap between the spans of X and Y using A=svd(X)$u, B=svd(Y)$u
## code written by Fabian Scheipl
trace_lv <- function(A, B, tol=1e-10){
  ## A, B orthnormal!!
  
  #Rolf Larsson, Mattias Villani (2001)
  #"A distance measure between cointegration spaces"
  
  if(NCOL(A)==0 | NCOL(B)==0){
    return(0)
  }
  
  if(NROW(A) != NROW(B) | NCOL(A) > NROW(A) | NCOL(B) > NROW(B)){
    return(NA)
  }
  
  trace <- if(NCOL(B)<=NCOL(A)){
    sum(diag(t(B) %*% A %*% t(A) %*% B))
  } else {
    sum(diag(t(A) %*% B %*% t(B) %*% A))
  }
  trace
}



#####################
## function to compute identifiability checks, part of FDboost 0.0-13
check_ident <- function(X1, L, Bs, K, xname, penalty, 
                        cumOverlap=FALSE, 
                        limits=NULL, yind=NULL, 
                        t_unique=NULL, 
                        id=NULL, 
                        X1des=NULL, ind0=NULL, xind=NULL, 
                        giveWarnings = TRUE){
  
  ## center X1 per column
  X1 <- scale(X1, scale=FALSE)
  
  #print("check.ident")
  ## check whether (number of basis functions in Bs) < (number of relevant eigenfunctions of X1)
  evls <- svd(X1, nu=0, nv=0)$d^2 # eigenvalues of centered fun. cov.
  evls[evls<0] <- 0
  maxK <- max(1, min(which((cumsum(evls)/sum(evls)) >= .995)))
  bsdim <- ncol(Bs) # number of basis functions in Bs
  #if(maxK < bsdim){
  #  warning("<k> (" , bsdim , ") larger than effective rank of <", xname, "> (", maxK, "). ", 
  #          "Effect identifiable only through penalty.")
  #}
  ## automatically use less basis-functions in case of problems?
  ## you would have to change args$knots accordingly
  
  ### compute condition number of Ds^t Ds
  Ds <- (X1 * L) %*% Bs
  DstDs <- crossprod(Ds)
  e_DstDs <- try(eigen(DstDs))
  e_DstDs$values <- pmax(0, e_DstDs$values) # set negative eigenvalues to 0
  logCondDs <- log10(e_DstDs$values[1]) - log10(tail(e_DstDs$values, 1))
  if(giveWarnings & logCondDs > 6 & is.null(limits)){
    warning("condition number for <", xname, "> greater than 10^6. ", 
            "Effect identifiable only through penalty.")
  }
  
  ### compute condition number of Ds^t Ds for subsections of Ds accoring to limits
  logCondDs_hist <- NULL
  
  # look at condition number of Ds for all values of yind for historical effect
  # use X1des, as this is the marginal design matrix using the limits
  if(!is.null(limits)){ 
    ind0Bs <- ((!ind0)*1) %*% Bs # matrix to check for 0 columns
    ## implementation is suitable for common grid of t, maybe with some missings
    ## common grid is assumed if Y(t) is observed at least in 80% for each point 
    if( all(table(yind)/max(id)>0.8) ){
      if(is.null(t_unique)) t_unique <- sort(unique(yind))
      logCondDs_hist <- rep(NA, length=length(t_unique))
      for(k in 1:length(t_unique)){
        Ds_t <- X1des[yind==t_unique[k], ] # get rows of Ds corresponding to yind
        ind0Bs_t <- ind0Bs[yind==t_unique[k], ] # get rows of ind0Bs corresponding to yind
        # only keep columns that are not completely 0, otherwise matrix is always rank deficient
        # idea: only this part is used to model y(t) at this point
        # also delete if not perfectly but almost zero, for all spline bases
        Ds_t <- Ds_t[ , apply(ind0Bs_t, 2, function(x) !all(abs(x)<10^-1) ), drop=FALSE ] 
        if(dim(Ds_t)[2]!=0){ # for matrix with 0 columns does not make sense
          DstDs_t <- crossprod(Ds_t)
          e_DstDs_t <- try(eigen(DstDs_t))
          e_DstDs_t$values <- pmax(0, e_DstDs_t$values) # set negative eigenvalues to 0
          logCondDs_t <- log10(e_DstDs_t$values[1]) - log10(tail(e_DstDs_t$values, 1))
          logCondDs_hist[k] <- logCondDs_t
        }
        ## matplot(xind, Bs, type="l", lwd=2, ylim=c(-2,2)); rug(xind); rug(yind, col=2, lwd=2)
        ## matplot(knots[1:ncol(Ds_t)], t(Ds_t), type="l", lwd=1, add=TRUE)
        ## lines(t_unique, logCondDs_hist-6, col=2, lwd=4)
      }
      names(logCondDs_hist) <- round(t_unique,2)
      
      ### implementation for seriously irregular observation points t
    }else{
      # use the mean number grid points, in the case of irregular t
      #t_unique <- seq(min(yind), max(yind), length=round(mean(table(id))))
      ### use quantiles of yind, as only at places with observations effect can be identifiable
      ### using quntiles prevents Ds_t from beeing completely empty
      if(is.null(t_unique))  t_unique <- quantile(yind, probs=seq(0,1,length=round(mean(table(id)))) )
      names(t_unique) <- NULL
      logCondDs_hist <- rep(NA, length=length(t_unique)-1)
      for(k in 1:(length(t_unique)-1)){
        # get rows of Ds corresponding to t_unique[k] <= yind < t_unique[k+1]
        Ds_t <- X1des[(t_unique[k] <= yind) & (yind < t_unique[k+1]), ] 
        ind0Bs_t <- ind0Bs[(t_unique[k] <= yind) & (yind < t_unique[k+1]), ]
        # for the last interval: include upper limit
        if(k==length(t_unique)-1){
          Ds_t <- X1des[(t_unique[k] <= yind) & (yind <= t_unique[k+1]), ]
          ind0Bs_t <- ind0Bs[(t_unique[k] <= yind) & (yind <= t_unique[k+1]), ]
        } 
        # only keep columns that are not completely 0, otherwise matrix is always rank deficient
        # idea: only this part is used to model y(t) at this point
        # also delete if not perfectly but almost zero, for all spline bases
        Ds_t <- Ds_t[ , apply(ind0Bs_t, 2, function(x) !all(abs(x)<10^-1) ), drop=FALSE] 
        if(dim(Ds_t)[2]!=0){ # for matrix with 0 columns does not make sense
          DstDs_t <- crossprod(Ds_t)
          e_DstDs_t <- try(eigen(DstDs_t))
          e_DstDs_t$values <- pmax(0, e_DstDs_t$values) # set negative eigenvalues to 0
          logCondDs_t <- log10(e_DstDs_t$values[1]) - log10(tail(e_DstDs_t$values, 1))
          logCondDs_hist[k] <- logCondDs_t
        }
      }
      names(logCondDs_hist) <- round(t_unique[-length(t_unique)],2)
    }
    if(giveWarnings & any(logCondDs_hist > 6)){
      # get the last entry of t, for which the condition number is >10^6
      temp <- names(which.max(which(logCondDs_hist > 6)))
      warning("condition number for <", xname, "> considering limits of historical effect ", 
              "greater than 10^6, for some time-points up to ", temp, ". ",
              "Effect in this region identifiable only through penalty.")
    }
  } ## end of computation of logCondDs_hist for historical effects
  
  
  ## measure degree of overlap between the spans of ker(t(X1)) and W%*%Bs%*%ker(K)
  ## overlap after Larsson and Villani 2001, Scheipl and Greven, 2014
  
  tryNA <- function(expr){
    ret <- try(expr, silent = TRUE)
    if(any(class(ret)=="try-error")) return(NA)
    return(ret)
  }
  tryNull <- function(expr){
    ret <- try(expr, silent = TRUE)
    if(any(class(ret)=="try-error")) return(matrix(NA, 0, 0))
    return(ret)
  }
  
  ### get special measures for kernel overlap of WB_s(P_s) with subset of Xobs
  ### overlap measure of Larsson and Villani 2001
  ### as proposed by Scheipl and Greven 2015
  getOverlap <- function(subset, X1, L, Bs, K){
    # In the case that all observations are 0, kernel is everything -> kernel overlap
    if(all(X1[ , subset]==0)){
      return(5)
    }
    KeXsub <- tryNull(Null(t(X1[ , subset])))
    if(ncol(KeXsub)==0){ # no null space
      return(0)
    }
    KePen2sub <- tryNull(diag(L[1,subset]) %*% Bs[subset,] %*% Null(K))
    overlapSub <- tryNA(trace_lv(svd(KeXsub)$u, svd(KePen2sub)$u))
    return(overlapSub)
  }
  
  cumOverlapKe <- NULL
  overlapKe <- NULL
  overlapKeComplete <- NULL
  
  #   ## cumulative overlap for historical model in the special case of s<t
  #   if(cumOverlap){   
  #     restm <- ncol(X1) %% 10 # rest of modulo calculation 
  #     ntemp <- (ncol(X1)-restm)/10 # group-size without rest
  #     ## subset with 1/10, 2/10, ..., 10/10 of the observation points
  #     if(restm > ntemp){ # case that rest is bigger than group size
  #       subs <- c(list(1:restm), lapply(1:8, function(i) 1:(restm+i*ntemp)), list(1:ncol(X1)))
  #     }else{
  #       subs <- c(lapply(1:9, function(i) 1:(restm+i*ntemp)), list(1:ncol(X1)))
  #     }
  #     cumOverlapKe <- sapply(subs, getOverlap, X1=X1, L=L, Bs=Bs, K=K)
  #     overlapKe <- max(cumOverlapKe, na.rm = TRUE) #cumOverlapKe[[length(cumOverlapKe)]]
  #     
  #   }else{ # overlap between whole matrix X and penalty
  #     overlapKe <- getOverlap(subset=1:ncol(X1), X1=X1, L=L, Bs=Bs, K=K)
  #   } 
  #   print("overlapKe")
  #   print(overlapKe)
  #   plot( seq(min(t_unique), max(t_unique), l=10), cumOverlapKe, ylim=c(0,1))
  
  
  ## sequential overlap for historical model with general integraion limits
  if(!is.null(limits)){  
    
    subs <- list()
    for(k in 1:length(t_unique)){
      subs[[k]] <- which(limits(s=xind, t=t_unique[k]))
    }
    cumOverlapKe <- sapply(subs, getOverlap, X1=X1, L=L, Bs=Bs, K=K)
    overlapKe <- max(cumOverlapKe, na.rm = TRUE) #cumOverlapKe[[length(cumOverlapKe)]]
    
  }else{ # overlap between whole matrix X and penalty
    overlapKe <- getOverlap(subset=1:ncol(X1), X1=X1, L=L, Bs=Bs, K=K)
  }
  
  overlapKeComplete  <- getOverlap(subset=1:ncol(X1), X1=X1, L=L, Bs=Bs, K=K)

  if(giveWarnings & overlapKe >= 1){
    warning("Kernel overlap for <", xname, "> and the specified basis and penalty detected. ",
            "Changing basis for X-direction to <penalty='pss'> to make model identifiable through penalty. ", 
            "Coefficient surface estimate will be inherently unreliable.") 
    penalty <- "pss"
  }
  
  return(list(logCondDs=logCondDs, logCondDs_hist=logCondDs_hist,  
              overlapKe=overlapKe, cumOverlapKe=cumOverlapKe, 
              overlapKeComplete=overlapKeComplete,
              maxK=maxK, penalty=penalty))
}




## get diagnostic measures of the FDboost model
## code by Fabian Scheipl, from identifiability paper
getDiags <- function(m, data, bl=2, cut=.995){
  
  
  if(is.null(m)){
    diags <- vector(6, mode="list")
    names(diags) <- c("logCondDs", "logCondDs_hist", "overlapKe", 
                      "cumOverlapKe", "maxK", "penalty" )
    return( diags )
  }
  
  
  tryNA <- function(expr){
    ret <- try(expr, silent = TRUE)
    if(any(class(ret)=="try-error")) return(NA)
    return(ret)
  }
  tryNull <- function(expr){
    ret <- try(expr, silent = TRUE)
    if(any(class(ret)=="try-error")) return(matrix(NA, 0, 0))
    return(ret)
  }
  
  #   #get Ds=XWBs & Ps
  #   Bs <- (m$smooth[[1]]$margin[[1]]$X[seq(1, 
  #                                          ncol(data$Y)*ncol(data$X1), 
  #                                          by=ncol(data$Y)), ])
  #   Ds <- (data$X1 * m$pffr$ff[[1]]$L) %*% Bs
  #   
  #   DstDs <- crossprod(Ds)
  #   Ps <- m$sig2 * m$sp[1] *  m$smooth[[1]]$margin[[1]]$S[[1]]
  
  limitsDefault <- function(s, t) {
    (s < t) | (s == t)
  }
  
  
  if(any(class(m)=="FDboost")){
    #   Ds <- X1des
    #   Ps <- K1
    #   Bs <- Bs
    Bs <- get("args", environment(m$baselearner[[bl]]$dpp))$Bs
    L <- get("args", environment(m$baselearner[[bl]]$dpp))$L
    Xobs <- scale(m$baselearner[[bl]]$get_data()[[1]], scale=FALSE)
    ### use Ds like for functional model (NOT historical model)
    Ds <- (Xobs * L) %*% Bs
    #Ds <- get("args", environment(m$baselearner[[bl]]$dpp))$X1des
    Ps <- get("args", environment(m$baselearner[[bl]]$dpp))$K1
    
    D <- extract(m, "design", which=bl)[[1]]
    ### multiply the penalty matrix with the corresponding lambda??
    ### should be irrelevant for boosting, as there is only one lambda in both directions
    P <- extract(m, "penalty", which=bl)[[1]]
    # P <- extract(m, "lambda")[[bl]]
    
    
    ## get information necessary for ident_check()
    xind <- get("args", environment(m$baselearner[[bl]]$dpp))$s
    yind <- m$yind
    id <- m$id
    
    ind0 <- !t(outer( xind, yind, limitsDefault) )
    
    X1 <- m$baselearner[[bl]]$get_data()[[1]]
    X1des <- X1[id, ] 
    X1des[ind0] <- 0    
    #X1des <- Matrix(X1des, sparse=TRUE) # convert into sparse matrix
    # X1des <- X1des * Lnew
    X1des <- X1des %*% Bs
    
    
  }else{
    #     if(m$long){
    #       Bs <- unique(m$smooth[[bl]]$margin[[1]]$X)
    #     }else{
    #       Bs <- (m$smooth[[bl]]$margin[[1]]$X[seq(1, 
    #                                               ncol(data$Y)*ncol(data$X1), 
    #                                               by=ncol(data$Y)), ])
    #     }

    Bs <- unique(m$smooth[[bl]]$margin[[1]]$X)
    L <- m$pffr$ff[[bl-1]]$L
    Xobs <- scale(m$pffr$ff[[bl-1]]$LX / L, scale=FALSE) 
    Ds <- (Xobs * L) %*% Bs
    #Ps <- m$sig2 * m$sp[3] *  m$smooth[[bl]]$margin[[1]]$S[[1]]
    Ps <- m$smooth[[bl]]$margin[[1]]$S[[1]]
    
    ### get desig and penalty matrix of historical effect
    D <- predict(m, type="lpmatrix", reformat = FALSE)[ ,(m$smooth[[bl]]$first.para):(m$smooth[[bl]]$last.para)]
    
    if(bl==2){
      where.sp <- c(2,3)
    }else{
      where.sp <- c(4,5)
    }
    
    Pt <- m$smooth[[bl]]$margin[[2]]$S[[1]] 
    #P <- (m$sig2*m$sp[where.sp[1]]) * kronecker(Ps , diag(ncol(Pt))) +
    #  (m$sig2*m$sp[where.sp[2]]) * kronecker(diag(ncol(Ps)), Pt)
    
    P <- kronecker(Ps , diag(ncol(Pt))) + kronecker(diag(ncol(Ps)), Pt)
    
    ## get information necessary for ident_check()
    xind <- data$s
    yind <- data$tlong
    id <- data$id
    
    ind0 <- !t(outer( xind, yind, limitsDefault) )
    
    X1 <- data[[paste0("X", bl-1)]]
    X1des <- X1[id, ] 
    X1des[ind0] <- 0    
    #X1des <- Matrix(X1des, sparse=TRUE) # convert into sparse matrix
    # X1des <- X1des * Lnew
    X1des <- X1des %*% Bs
    
  }
  
  
  ########### use check_ident() from package FDboost
  #browser()
  ident_check <- check_ident(X1=Xobs, L=L, Bs=Bs, K=Ps, xname="test", penalty="ps", 
                                cumOverlap = FALSE, limits = limitsDefault, 
                        yind = yind, t_unique=sort(unique(yind)), 
                        id = id, X1des = X1des, ind0 = ind0, xind = xind, 
                        giveWarnings=FALSE)
  
  ident_check
  
}





