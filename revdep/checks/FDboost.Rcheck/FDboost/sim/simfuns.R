###############################################################################
# changed code of Fabian Scheipl lfpr3.R for data generation and fit of pffr
# Author: Sarah Brockhaus
###############################################################################
   

# M=50; ni=1; Gy=30; Gx=20; snrEps=5; snrE=0; snrB=1; balanced=TRUE; nuisance=0
# scenario="1"
makeData <- function(M=50, # number of subjects
        ni=1, # mean number of observations per subject
        Gy=30, # number of grid points for t 
        Gx=30, # number of grid points for s
        snrEps=5, # signal-to-noise ratio
        snrE=0, # smooth errors per curve E_ij? - currently not used
        snrB=1, # relative importance of random effects
        scenario="1",
        balanced=TRUE, # in the case of reapeated measures - currently always TRUE
        nuisance=0, # number of nuisance variables
        ...){
  
  # set df in base learners
  df <- 2.25
  
  # function for generating a functional variable as (spline-basis)%*%(random coefficient)
  rf <- function(x=seq(0,1,length=100), k=15) {
    drop(ns(x, k, int=TRUE) %*% runif(k, -3, 3))
  }
  
  # standardize matrix, so that colSums are zero
  zeroConstraint <- function(x){
    stopifnot(is.matrix(x))
    t(t(x) - colMeans(x))
  }
  
  # generate zero centered scalar covariates
  cenScalarCof <- function(n){
    z <- runif(n) - 0.5
    z <- z - mean(z)   
    z
  }
  
  # generate a smooth global intercept
  intf <- function(t){
    cos(3*pi*t^2) + 2
  }
  
  ## smooth function of scalar covariate and time
  g2zt <- function(z, t){ 2*(-t^2-0.1)*sin(pi*z+0.5) }
  
  # build a regular grid over the range of a scalar covariable
  zgrid <- function (z, l=40) seq(min(z), max(z), l=l)
  
  # functions for scalar effects
  beta2 <- function(z,t) 3*sin(-2*t) + 1
  beta3 <- function(z,t) -4*10^(z^2) + 5
  
  # generate bivariate coefficient surface
  test1 <- function(s, t, ss=0.3, st=0.4)
  { 
    dnorm(s, mean = 0.2, sd = 0.3)*dnorm(t, mean = 0.2, sd = 0.3) + 
      dnorm(s, mean = 0.6, sd = 0.3)*dnorm(t, mean = 0.8, sd = 0.25) - 0.5
  }
  #persp(outer(seq(0,1, l=20), seq(0,1, l=20), test1 ), ticktype="detailed", theta=30, phi=30 )
 
  # generate bivariate coefficient surface
  test2 <- function(s, t){
    1.5*sin(pi*t+0.3) * sin(pi*s)
  }
  #persp(outer(seq(0,1, l=20), seq(0,1, l=20), test2 ), ticktype="detailed", theta=30, phi=30 )
  
  
  #### function that generates nuisance variables, uses the variables within the functions
  # dgp1, ..., dgp5
  generateNuisance <- function(M, ni, Gy, Gx, snrEps, snrE, snrB, balanced, 
                               nuisance, scenario){
    
    if(nuisance==0) return(NULL)
    
    dataN <- data.frame(dummy=rep(NA, M*ni))
    
    # Define number of scalar and functional nuisance variables
    if(scenario %in% c("1")){
      nrX <- ceiling(nuisance/2)
      nrz <- nuisance - nrX      
    }    
    if(scenario %in% c("2")){
      nrX <- 0
      nrz <- nuisance     
    }    
    if(scenario %in% c("3")){
      nrX <- nuisance
      nrz <- 0    
    }
    
    pffrN <- c()
    boostN <- c()   
    namesVariables <- c()
    
    # Generate scalar and functional nuisance variables
    if(nrz>0){
      for(i in 1:nrz){
        dataN <- data.frame(dataN, cenScalarCof(M*ni))
      }
      pffrN <- paste("s(", paste("zn", 1:nrz, sep=""), ", k=5)", sep="")
      boostN <- paste("bbsc(", paste("zn", 1:nrz, sep=""), ", df=",df, ")", sep="")
      namesVariables <- c(namesVariables, paste("zn", 1:nrz, sep=""))
      
      # use factor nuisance variables in scenario 2
      if(scenario=="2"){        
        for(i in 1:ceiling(nrz/4)){
          dataN[,i+1] <- sample(c(-1,1), M*ni, replace=TRUE) # +1 because of dummy variable in first column
        }
        pffrN[1:ceiling(nrz/4)] <- paste( paste("zn", 1:ceiling(nrz/4), sep=""), sep="")
        boostN[1:ceiling(nrz/4)] <- paste("bols(", paste("zn", 1:ceiling(nrz/4), sep=""), 
                                          ", df=", df ,", intercept=FALSE)", sep="") 
        namesVariables[1:ceiling(nrz/4)] <- paste("zn", 1:ceiling(nrz/4), "_bin", sep="")
      }
    }    
    if(nrX>0){
      s <- seq(0, 1, l=Gx)
      # Mix nuisance variables of different complexity in X
      for(i in 1:nrX){
        dataN <- data.frame(dataN, I(zeroConstraint( t(replicate(M*ni, rf(s, k=5))) )) )
        #if(i%%2==0) dataN <- data.frame(dataN, I(zeroConstraint( t(replicate(M*ni, rf(s, k=5))) )) )
      }
      pffrN <- c(pffrN, paste("ff(", paste("Xn", 1:nrX, sep=""), ", yind=t)", sep=""))
      boostN <- c(boostN, paste("bsignal(", paste("Xn", 1:nrX, sep=""), ", s, df=", df ,")", sep="") )
      namesVariables <- c(namesVariables, paste("Xn", 1:nrX, sep=""))
    }  
    dataN <- dataN[,-1]  
    attr(dataN, "pffrN") <- pffrN
    attr(dataN, "boostN") <- boostN
    attr(dataN, "namesVariables") <- namesVariables
    names(dataN) <- c(  if(nrz>0) paste("zn", 1:nrz, sep=""), 
                        if(nrX>0) paste("Xn", 1:nrX, sep=""))    
    return(dataN)    
  }

  # call the function generateNuisance() and add the nuisance variables to the dataset  
  addNuisance <- function(data, M, ni, Gy, Gx, snrEps, snrE, snrB, balanced, 
                          nuisance, scenario){
    dataN <- generateNuisance(M, ni, Gy, Gx, snrEps, snrE, snrB, 
                              balanced, nuisance, scenario)
    data <- c(data, dataN)   
    attr(data, "pffrN") <- attr(dataN, "pffrN")
    attr(data, "boostN") <-  attr(dataN, "boostN")
    attr(data, "namesVariables") <-  attr(dataN, "namesVariables")
    return(data)
  }
  
  # Scenario with one scalar and one functional covariable
  dgp1 <- function(M, ni, Gy, Gx, snrEps, snrE, snrB, balanced, nuisance){
    id <- if(balanced){
      gl(M, ni)  # generate factor levels  
    }  else {
      factor(c(1:M, sample(1:M, M*ni-M, repl=TRUE, prob=sqrt(1:M))))
    }
    t <- seq(0, 1, l=Gy)
    s <- seq(0, 1, l=Gx)
    tgrid <- sgrid <- seq(0,1, l=40)
    tgrid100 <- seq(0,1, l=100)
    
    # Intercept
    int <- matrix(intf(t), nrow=ni*M, ncol=Gy, byrow=TRUE)
    
    # Scalar covariate
    z1 <- cenScalarCof(M*ni)      
    fz1 <- zeroConstraint(outer(z1, t, g2zt))  # g(z1,t)             
    
    # Functional covariate
    X1 <- zeroConstraint( t(replicate(M*ni, rf(s, k=10))) )
    L <- integrationWeights(X1=X1, xind=s)  # Riemann integration weights      
    beta1.st <- outer(s, t, test1)
    X1f <- (L*X1)%*%beta1.st
    
    # Compute response
    Ytrue <- int + fz1 + X1f            
    Y <- Ytrue + sd(as.vector(Ytrue))/snrEps * matrix(scale(rnorm(M*ni*Gy)), 
                                                      nrow=M*ni, ncol=Gy)      
    true_g0 <- intf(tgrid100)
    true_fz1 <- zeroConstraint(outer(zgrid(z1), tgrid, g2zt))
    true_bst1 <- outer(sgrid, tgrid, test1)
    
    data <- list(id=id, Y=Y, z1=z1, X1=X1, s=s, t=t, 
                 Ytrue=I(Ytrue), int=int, fz1=fz1, X1f=X1f,
                 true_g0=true_g0, true_fz1=true_fz1, true_bst1=true_bst1)      
    data <- addNuisance(data, M, ni, Gy, Gx, snrEps, snrE, snrB, balanced, nuisance, scenario)
    return(data)
  }  
  
  # Scenario with three different effects of scalar covariables 
  dgp2 <- function(M, ni, Gy, snrEps, snrE, snrB, balanced, nuisance){
    id <- if(balanced){
      gl(M, ni)  # generate factor levels
    }else{
      factor(c(1:M, sample(1:M, M*ni-M, repl=TRUE, prob=sqrt(1:M))))
    }
    t <- seq(0, 1, l=Gy)
    tgrid <- sgrid <- seq(0,1, l=40)
    tgrid100 <- seq(0,1, l=100)
    
    # Intercept
    int <- matrix(intf(t), nrow=ni*M, ncol=Gy, byrow=TRUE)
    
    z1 <- cenScalarCof(M*ni)
    z2 <- cenScalarCof(M*ni)
    z3 <- cenScalarCof(M*ni)
    z4 <- sample(c(-1,1), M*ni, replace=TRUE)
    
    # first two effects are centered around 0 by construction
    fz1 <-  zeroConstraint(outer(z1, t, g2zt))  # g(z1,t)                   
    fz2 <-  zeroConstraint(outer(z2, t, beta3))  # g(z2)     
    fz3 <-  outer(z3, t, function(z,t) beta2(z,t)*z) # beta(t)*z3
    fz4 <-  outer(z4, t, function(z,t) beta2(z,t)*z) # beta(t)*z4
    
    true_g0 <- intf(t=tgrid100)      
    true_fz1 <- zeroConstraint(outer(zgrid(z1), tgrid, g2zt))
    true_fz2 <- beta3(z=zgrid(z3, l=100), t=NA) 
    true_fz2 <- true_fz2 - mean(true_fz2)
    true_fz3 <-beta2(z=NA, t=tgrid100)
    true_fz4 <- beta2(z=NA, t=tgrid100) 
    
    Ytrue <- int + fz1 + fz2 + fz3 + fz4
    
    Y <- Ytrue + sd(as.vector(Ytrue))/snrEps * matrix(scale(rnorm(M*ni*Gy)), nrow=M*ni, ncol=Gy)
    data <- list(id=id, Y=Y, t=t, z1=z1, z2=z2, z3=z3, z4=z4,
                 Ytrue=Ytrue, int=int,
                 fz1=fz1, fz2=fz2, fz3=fz3, fz4=fz4,
                 true_g0=true_g0, true_fz1=true_fz1, true_fz2=true_fz2, 
                 true_fz3=true_fz3, true_fz4=true_fz4)
    data <- addNuisance(data, M, ni, Gy, Gx, snrEps, snrE, snrB, balanced, nuisance, scenario)    
    return(data)      
  }
  
  # Scenario with TWO functional covariates
  dgp3 <- function(M, ni, Gy, Gx, snrEps, snrE, snrB, balanced, nuisance){
    id <- if(balanced){
      gl(M, ni)  # generate factor levels
    }  else {
      factor(c(1:M, sample(1:M, M*ni-M, repl=TRUE, prob=sqrt(1:M))))
    }
    t <- seq(0, 1, l=Gy)
    s <- seq(0, 1, l=Gx)
    tgrid <- sgrid <- seq(0,1, l=40)
    tgrid100 <- seq(0,1, l=100)
    
    # Intercept
    int <- matrix(intf(t), nrow=ni*M, ncol=Gy, byrow=TRUE)
    
    # Functional covariates
    X1 <- zeroConstraint( t(replicate(M*ni, rf(s, k=5))) )
    L <- integrationWeights(X1=X1, xind=s)  # Riemann integration weights      
    betast_1 <- outer(s, t, test1)
    X1f <- (L*X1)%*%betast_1
    
    X2 <- zeroConstraint( t(replicate(M*ni, rf(s, k=5))) )
    L <- integrationWeights(X1=X2, xind=s)  # Riemann integration weights      
    #betast_2 <- outer(s, t, betast2)
    betast_2 <- outer(s, t, test2)
    X2f <- (L*X2)%*%betast_2
        
    # Compute response
    Ytrue <- int + X1f + X2f     
    Y <- Ytrue + sd(as.vector(Ytrue))/snrEps * matrix(scale(rnorm(M*ni*Gy)), 
                                                      nrow=M*ni, ncol=Gy)      
    true_g0 <- intf(tgrid100)
    true_bst1 <- outer(sgrid, tgrid, test1)
    true_bst2 <- outer(sgrid, tgrid, test2)
    
    data <- list(id=id, Y=Y, X1=X1, X2=X2, 
                 s=s, t=t,
                 Ytrue=I(Ytrue), int=int, X1f=X1f, X2f=X2f, 
                 true_g0=true_g0, true_bst1=true_bst1, true_bst2=true_bst2)  
    data <- addNuisance(data, M, ni, Gy, Gx, snrEps, snrE, snrB, balanced, nuisance, scenario)    
    return(data)
  }
  
  
  #### generate data
  data <- switch(scenario,
                 "1" = dgp1(M, ni, Gy, Gx, snrEps, snrE, snrB, balanced, nuisance),
                 "2" = dgp2(M, ni, Gy, snrEps, snrE, snrB, balanced, nuisance),
                 "3" = dgp3(M, ni, Gy, Gx, snrEps, snrE, snrB, balanced, nuisance))
   
  # Formula of function pffr()
  formulaPffr <- switch(scenario,
                        "1" = "s(z1, k=5) + ff(X1, yind=t)",
                        "2" = "s(z1, k=5) + c(s(z2)) + z3 + z4",
                        "3" = "ff(X1, yind=t) + ff(X2, yind=t)")
  
  # effects of nuisance variables  
  formulaPffr <- as.formula(paste("Y ~ ", formulaPffr, 
                                  if(!is.null(attr(data, "pffrN"))){
                                    paste(c(" ", attr(data, "pffrN")), collapse= " + ")
                                  } ))
 
  # Formula of function FDboost()  
  formulaFDboost <- switch(scenario,
                           "1" = paste("1 + bbsc(z1, df=", df,") + bsignal(X1, s, df=", df,")", sep=""),
                           "2" = paste("1 + bbsc(z1, df=", df,") + c(bbs(z2, knots=10, df=", df,")) + bols(z3, df=", df,", intercept=FALSE) + bols(z4, df=", df,", intercept=FALSE)", sep=""),
                           "3" = paste("1 + bsignal(X1, s, df=", df,") + bsignal(X2, s, df=", df,")" , sep=""))
  # effects of nuisance variables    
  formulaFDboost <- as.formula(paste("Y ~ ", formulaFDboost, 
                                     if(!is.null(attr(data, "boostN"))){
                                    paste(c(" ", attr(data, "boostN")), collapse= " + ")
                                  } )) 
  timeformula <- formula(paste("~ bbs(t, knots=10, df=", df, ")", sep=""))
  
  var1 <- c("Int", "z1", "X1") 
  namesVariables <- switch(scenario,
                             "1" = var1,
                             "2" = c("Int", paste("z", 1:3, sep=""), "z4_bin"),
                             "3" = c("Int", paste("X", 1:2, sep="")),                      
                             "4" = c(var1, "id"),
                             "5" = c(var1, "id", "id_u"))
  namesVariables <- c(namesVariables, attr(data, "namesVariables"))
      
  return(structure(data,
                   sigmaEps = sd(as.vector(data$Ytrue))/snrEps,
                   formulaPffr=formulaPffr,
                   formulaFDboost=formulaFDboost,
                   timeformula=timeformula,
                   namesVariables=namesVariables,
                   call=match.call()))
}

if(FALSE){
    data1 <- makeData(scenario="1", seed=45)
    data2 <- makeData(scenario="2", seed=14)  # , M=100, Gx=30, Gy=30
    data3 <- makeData(scenario="3", seed=33)
    
    data1 <- makeData(scenario="1", nuisance=2, seed=12)
    data2 <- makeData(scenario="2", nuisance=2, seed=34)
    data3 <- makeData(scenario="3", nuisance=2, seed=37)
}


### Fit model using pffr() of package refund using the data and information of data from makeData()
fitModelPffr <- function(data, 
        bs.yindex = list(bs = "ps", k = 5, m = c(2, 1)), # default for expasion along t
        bs.int = list(bs = "ps", k = 20, m = c(2, 1)), ...){ # default for expansion of intercept
      
  formula <- as.formula(attr(data, "formulaPffr"))  
  t <- data$t
  
  time <- system.time(
    mem <- mem_change(m <- try(pffr(eval(formula), yind=t, data=data, 
                   bs.yindex=bs.yindex, bs.int=bs.int)))
  )[3]

  if(any(class(m) != "try-error")){
    m$runTime <- time 
    m$memory <- mem 
    return(m)
  }else return(NULL)
}

if(FALSE){
    m1 <- fitModelPffr(data1)
    m2 <- fitModelPffr(data2)
    m3 <- fitModelPffr(data3)
}


### Fit model using FDboost() that is based on package mboost
fitModelMboost <- function(data, 
                           control=boost_control(mstop=100, nu=0.1), # settings of mboost
                           grid=seq(10, 200, by=10),
                           m_max=1000, 
                           optimizeMstop=TRUE,
                           useArray=TRUE,
                           ...){ 
  
  #scenario <- attr(data, "call")$scenario    
  formula <- attr(data, "formulaFDboost")
  timeformula <- attr(data, "timeformula")
  
  ## do not use the array structure -> fit the model in long format
  if(!useArray){
    data$Y <- as.vector(t(data$Y)) 
    data$myid <- rep(1:nrow(data$X1), each=length(data$t))
    data$t <- rep(data$t, times=nrow(data$X1))
  }else{
    data$myid <- NULL
  }
  
  time <- system.time({ 
    
    if(useArray){
      mem <- mem_change(m <- try(FDboost(eval(formula), data=data, control=control,
                                         timeformula = timeformula))) 
    }else{ ## data in long format
      mem <- mem_change(m <- try(FDboost(eval(formula), data=data, control=control,
                                         timeformula = timeformula, id=~myid))) 
    }
 
    # plot(m$risk())
    if(any(class(m)=="FDboost") & optimizeMstop){
      
      # set up two different splittings of the data
      if(attr(data, "call")$scenario!=5){ id <- data$id        
      }else id <- NULL
      
      set.seed(attr(data, "call")$seed)
      folds1 <- cvMa(ydim=m$ydim, type="bootstrap", B=10) 
      
      set.seed(attr(data, "call")$seed + 100)
      folds2 <- cvMa(ydim=m$ydim, type="bootstrap", B=10)
            
      # Increase maximal number of boosting iterations
      increaseGrid <- TRUE
      while(increaseGrid){
        
        # cvm <- try(suppressMessages(cvrisk(m, papply = lapply, grid=grid)))
        if(Sys.info()["sysname"]=="Linux"){ # use 10 cores on Linux
          cvm <- try(cvrisk(m, folds = folds1, grid=grid, mc.cores=10), silent=TRUE)        
        }else cvm <- try(cvrisk(m, folds = folds1, grid=grid))          
        
        # Try for a second time if cvrisk() stops with error
        if(class(cvm)!="cvrisk"){
          if(Sys.info()["sysname"]=="Linux"){
            cvm <- try(cvrisk(m, folds = folds2, grid=grid, mc.cores=10), silent=TRUE)        
          }else cvm <- try(cvrisk(m, folds = folds2, grid=grid)) 
          cat(paste("2nd calculation of cvrisk, class(cvm) =", class(cvm), "\n", sep=" "))
        }
        
        print(max(grid))
        
        if(class(cvm) =="cvrisk"){
          # stop while-loop if optimal mstop is not the largest possible value
          # and if mstop is grater or equal to m_max 
          increaseGrid <- mstop(cvm)==max(grid) & !mstop(cvm) >= m_max
        }else increaseGrid <- FALSE # stop while-loop if cross validation was not possible         
        
        # Continue to search th optimal mstop in the next 200 iterations
        grid <- seq(max(grid)+10, max(grid)+200, by=10)
      }
      # plot(cvm)      
    }else cvm <- NULL 
  })[3]
  
  if(!optimizeMstop){
    m$runTime <- time
    m$memory <- mem
    return(m)
  } 
  
  if(any(class(m)=="FDboost") & class(cvm)=="cvrisk"){
    #if(mstop(cvm)==max(grid)) warning("mstop is maximal number in grid.")
    #if(mstop(cvm)==min(grid)) warning("mstop is minimal number in grid.")
    m <- m[mstop(cvm)]
    m$runTime <- time 
    m$memory <- mem 
    return(m)
  }else return(NULL)    
}

if(FALSE){
  m1m <- fitModelMboost(data1, m_max=500) # 300
  m2m <- fitModelMboost(data2, m_max=500)
  m3m <- fitModelMboost(data3, m_max=500)
}



### calculate some errors/ measures for goodness of fit of the model 
# from fitModelMboost() and fitModelPffr()
getErrors <- function(data, m=NULL, plotModel=FALSE){
  
  scenario <- attr(data, "call")$scenario
  
  # In case that model was not fitted: m=NULL
  errorE <- c(msey=NA, mseg0=NA, 
              msefx1=NA, msefx2=NA,
              msefz1=NA, msefz2=NA, msefz3=NA, msefz4=NA) 
  
  relerrorE <- errorE
  names(relerrorE) <- paste("rel", names(errorE), sep="")

  if(is.null(m)) return(as.list(c(errorE, relerrorE, relmseyT=NA, 
                                  funRsqrt=NA, mstop=NA, time.elapsed=NA)) )
    
  classM <- class(m)[1]
  
  # Save number of iterations
  mstop <- switch(classM,
                 "pffr" = NA,
                 "FDboost" = mstop(m),
                 "NULL" = NA)  
  
  
  # function to predict each component of linear predictor separately for FDboost() and pffr()
  # intercept =TRUE indicates that the first variable is an intercept 
  # to which the offset should be added 
  
  predictComponents <- function(object, intercept=TRUE){
    fit <- NULL
    if(any(class(object)=="FDboost")){
      fit <- predict(object=object, which=1:length(object$baselearner))
      if(intercept) fit[[1]] <- fit[[1]] + object$offset  # add offset
    }else{
      fit <- predict(object=object, type="terms")
      if(intercept) fit[[1]] <- fit[[1]] + coef(object, se=FALSE, seWithMean=FALSE)$pterms[1]  # add global constan intercept
    }
    return(fit)    
  }
  
  # Perdiciton for pffr and boostFD
  fit <- switch(classM,
                "pffr" = predictComponents(m),
                "FDboost" = predictComponents(m),
                "NULL" = NULL) 
  
  yhat <- switch(classM,
                 "pffr" = fitted(m),
                 "FDboost" = predict(m),
                 "NULL" = NULL)  

  ############### Calculate errors of estimated coefficients directly on coefficients
  if(any(class(m)=="pffr")){    
    cm <- coef(m, se=FALSE, seWithMean=FALSE)

    # functional intercept g0
    est_g0 <- drop(cm$smterms[[1]]$value) + cm$pterms[1]
    true_g0 <- data$true_g0
    #plot(true_g0, col=2, ylim=range(true_g0, est_g0)); points(est_g0)
    #mean((est_g0-true_g0)^2)
        
    if(scenario %in% c("1", "2")){
      true_fz1 <- data$true_fz1
      est_fz1 <- matrix(drop(cm$smterms[[2]]$value), ncol=40)
    }else est_fz1 <- true_fz1 <- NA
    
    if(scenario %in% c("2")){
      true_fz2 <- data$true_fz2
      true_fz3 <- data$true_fz3
      true_fz4 <- data$true_fz4
      est_fz2 <- drop(cm$smterms[[3]]$value)
      est_fz3 <- drop(cm$smterms[[4]]$value)
      est_fz4 <- drop(cm$smterms[[5]]$value)           
    }else{
      est_fz2 <- true_fz2 <- NA
      est_fz3 <- true_fz3 <- NA
      est_fz4 <- true_fz4 <- NA
    } 
    
    true_bst1 <- est_bst1 <- true_bst2 <- est_bst2 <- NA        
    if(scenario %in% c("1", "3")){
      whereX1 <- if(scenario %in% c("1")) 3 else 2
      true_bst1 <- data$true_bst1
      est_bst1 <- matrix(drop(cm$smterms[[whereX1]]$value), ncol=40)      
      if(scenario %in% c("3")){
        true_bst2 <- data$true_bst2
        est_bst2 <- matrix(drop(cm$smterms[[3]]$value), ncol=40)       
      }     
    }
  }

  ################ Calculate errors of estimated coefficients directly on coefficients
  if(any(class(m)=="FDboost")){  
    cm <- coef(m, n1=100)
    
    # functional intercept g0
    est_g0 <- cm$offset$value + cm$smterms[[1]]$value
    true_g0 <- data$true_g0
    
    if(scenario %in% c("1", "2")){
      true_fz1 <- data$true_fz1
      est_fz1 <- matrix(drop(cm$smterms[[2]]$value), ncol=40)
    }else est_fz1 <- true_fz1 <- NA
    
    if(scenario %in% c("2")){
      true_fz2 <- data$true_fz2
      true_fz3 <- data$true_fz3
      true_fz4 <- data$true_fz4
      est_fz2 <- cm$smterms[[3]]$value
      est_fz3 <- cm$smterms[[4]]$value
      est_fz3 <- drop(predict(m, which=4, newdata=list(z3=1, t=seq(0, 1, l=100 ))))
      est_fz4 <- drop(predict(m, which=5, newdata=list(z4=1, t=seq(0, 1, l=100 ))))
    }else{
      est_fz2 <- true_fz2 <- NA
      est_fz3 <- true_fz3 <- NA
      est_fz4 <- true_fz4 <- NA
    } 
    
    true_bst1 <- est_bst1 <- true_bst2 <- est_bst2 <- NA        
    if(scenario %in% c("1", "3", "4", "5")){
      whereX1 <- if(scenario %in% c("1")) 3 else 2
      true_bst1 <- data$true_bst1
      est_bst1 <- cm$smterms[[whereX1]]$value 
      if(scenario %in% c("3")){
        true_bst2 <- data$true_bst2
        est_bst2 <- cm$smterms[[3]]$value
        true_bst3 <- data$true_bst3     
      }     
    }
  }
            
  # Calculate MSE
  calcError <- function(x, xhat){
    if((length(x)==1 & is.na(x[1]))|(length(xhat)==1 & is.na(xhat[1]))){
      return(NA)
    }else{
      mean((x - xhat)^2)
    }
  }
  
  errorE <- c(msey = calcError(data$Ytrue, yhat),
                mseg0 = calcError(true_g0, est_g0),              
                msefx1 = calcError(true_bst1, est_bst1),
                msefx2 = calcError(true_bst2, est_bst2),
                msefz1 = calcError(true_fz1, est_fz1),
                msefz2 = calcError(true_fz2, est_fz2),
                msefz3 = calcError(true_fz3, est_fz3),
                msefz4 = calcError(true_fz4, est_fz4)) 
  
  # Calculate relative error for y, scaled with var_y(t)
  # irMSE with standardized with local variance at each point t
  relcalcErrorT <- function(x, xhat){  
    if( (length(x)==1 & is.na(x[1])) |  (length(xhat)==1 & is.na(xhat[1]))) return(NA)    

    # Calculation like functional R^2 - standardize with mu(t) (Ramsay, Silverman Chap 16)
      mut <- matrix(colMeans(x), nrow=nrow(x), ncol=ncol(x), byrow=TRUE)
    
      if (sum((x - mut)^2)==0){
        warning("Error is scaled by sigmaEps")
        stand <- attr(data, "sigmaEps")
      } else{
        stand <- colMeans((x-mut)^2)
        # If you would have to divide with 0 replace 0 with the minimum
        stand[stand==0] <- min(stand[stand!=0])
      }
      # Standardize with "variability" at each time-point
      err <- mean( colMeans((x-xhat)^2)/ stand)       
  return(err)
  }
  
  # calculate irMSE that is the MSE standardized by the global variance
  relcalcError <- function(x, xhat){
    if((length(x)==1 & is.na(x[1]))|(length(xhat)==1 & is.na(xhat[1]))) return(NA)    
    stopifnot(dim(x)==dim(xhat))  
    
    # Calculation like functional R^2 - standardize with global mu
    mu <- mean(x)
    stand <- mean((x-mu)^2)
    if (stand==0){
      warning("Error is scaled by sigmaEps")
      stand <- attr(data, "sigmaEps")
    }
    # Standardize with global "variability" 
    err <- mean( (x-xhat)^2) / stand   
    return(err)
  }
  
  relerrorE <- c(msey = relcalcError(data$Ytrue, yhat),
              mseg0 = relcalcError(true_g0, est_g0),              
              msefx1 = relcalcError(true_bst1, est_bst1),
              msefx2 = relcalcError(true_bst2, est_bst2),
              msefz1 = relcalcError(true_fz1, est_fz1),
              msefz2 = relcalcError(true_fz2, est_fz2),
              msefz3 = relcalcError(true_fz3, est_fz3),
              msefz4 = relcalcError(true_fz4, est_fz4))
  
  names(relerrorE) <- paste("rel", names(errorE), sep="")
    

 ### calculate errors on scale of response?
  
  # functional R^2
  funRsqrt <- function(x, xhat){
    stopifnot(dim(x)==dim(xhat))
    mut <- matrix(colMeans(x), nrow=nrow(x), ncol=ncol(x), byrow=TRUE)
    if (sum((x - mut)^2)==0) return(NaN)    
#    1 - ( sum((x - xhat)^2) / sum((x - mut)^2) )  # over all
    1 - mean( colSums((x - xhat)^2) / colSums((x - mut)^2) ) # over t
#    1 - mean( rowSums((x - xhat)^2) / rowSums((x - mut)^2) ) # over subjects
  }
  
  ret <- as.list(c(errorE, relerrorE, relmseyT = relcalcErrorT(data$Ytrue, yhat), 
                   funRsqrt=funRsqrt(x=data$Ytrue, xhat=yhat), 
                   mstop=mstop, time=m$runTime))
  
  # Save information for plotting
  if(plotModel==TRUE){
    est <- list(Ytrue=data$Ytrue, yhat=yhat, true_g0=true_g0, est_g0=est_g0,
                true_bst1=true_bst1, est_bst1=est_bst1, 
                true_bst2=true_bst2, est_bst2=est_bst2,
                true_fz1=true_fz1, est_fz1=est_fz1,
                true_fz2=true_fz2, est_fz2=est_fz2,
                true_fz3=true_fz3, est_fz3=est_fz3,
                true_fz4=true_fz4, est_fz4=est_fz4
                       ) 
    attr(ret, "est") <- est
  }
          
  return(ret)
}

if(FALSE){
    eNULL <- getErrors(NULL, NULL)
    length(eNULL)
    e1 <- getErrors(data1, m1)
    e2 <- getErrors(data2, m2) #, plotModel=TRUE
    e3 <- getErrors(data3, m3)
    length(e3)
    #cbind(eNULL, e1)
}

#true_fz2, est_fz2

if(FALSE){
  e1m <- getErrors(data1, m1m)
  e2m <- getErrors(data2, m2m)  #, plotModel=TRUE
  e3m <- getErrors(data3, m3m)
  cbind(e3, e3m)
}

# matplot that defaults to plotting lines
matlplot <- function(x,y, lty=1, col=1, lwd=.1, ...){
  if (missing(x)) {
    if (missing(y)) 
      stop("must specify at least one of 'x' and 'y'")
    else x <- 1L:NROW(y)
  }
  else if (missing(y)) {
    y <- x
    x <- 1L:NROW(y)
  }
  matplot(x,y, type="l", lty=lty, col=col, lwd=lwd, ...)
}


# function to plot true values, estimates of FDboost and estimates of pffr
# depends on results of function getErrors()
# use only a subset of 10 observations for plotting them as example
plotModel <- function(err, errB=NULL, data, theseSettings, subset=sample(1:theseSettings$M, size=10)){
  
  # build a regular grid over the range of a scalar covariable
  zgrid <- function (z, l=40) seq(min(z), max(z), l=l)
  
  
  est <- attr(err, "est")
  estB <- attr(errB, "est")
  
  model <- "pffr"
  modelB <- "mboost"
  
  # set errors to NULL if no model was fitted
  if(sum(is.na(err))==length(err)) err <- NULL
  if(sum(is.na(errB))==length(errB)) errB <- NULL
  
  # Plot nothing if none of the models could be fitted
  if(is.null(err) & is.null(errB)){
    return(NULL)
  }
        
  settingString <- paste(names(theseSettings), unlist(theseSettings), sep=":")
  settingString <- paste(settingString[!grepl("seed:", settingString) &
    !grepl("set:", settingString)], collapse=" ")
  #opar <- par()
  #on.exit(try(par(opar), silent=TRUE))
  
  par(xpd=TRUE, mar=par()$mar/2)
  #clrs <- alpha( rainbow(length(subset)), 200 )
  clrs <- rainbow(length(subset))
  
  #  plotIntercept <- function(m=m, mb=mb, fit=fit, fitB=fitB, errors=errors, errorsB=errorsB, data=data){
  plotIntercept <- function(){
    range <- range(est$est_g0, estB$est_g0, data$true_g0)
    plot(data$true_g0, main=bquote(g[0](t)), xlab="", ylab="", xaxt="n", ylim=range, type="l")
    if(!is.null(err)){
      plot(est$est_g0, xlab="", ylab="", xaxt="n", ylim=range, type="l", 
           main=bquote(paste("FAMM: ", hat(g[0])(t), ": reliMSE", phantom(x)%~~% .(round(err$relmseg0, 4)))))      
    }
    if(!is.null(errB)){
      plot(estB$est_g0, xlab="", ylab="", xaxt="n", ylim=range, type="l",
           main=bquote(paste("FLAM: ", hat(g[0])(t), ": reliMSE", phantom(x)%~~% .(round(errB$relmseg0, 4)))))
    }   
  }
  
  twoModels <- !( is.null(errB)|is.null(err) )

  if(theseSettings$scenario %in% c("1", "3"))
    layout(matrix(1:(10+5*twoModels), ncol=2+1*twoModels, byrow=TRUE))
    
  if(theseSettings$scenario %in% c("2"))
    layout(matrix(1:(14+7*twoModels), ncol=2+1*twoModels, byrow=TRUE))
  
  plotIntercept()
  
  # If pffr was not fitted but FDboost was, then change err and errB
  if(is.null(err) & !is.null(errB)){
    err <- errB
    est <- estB
    errB <- NULL
    estB <- NULL
    model <- modelB
  }
    
  # Plot g_1(z_1,t)
  if(theseSettings$scenario %in% c("1","2")){
    tgrid <- seq(0,1, l=40)
    # Only plot part of the points -> grid is visible
    k=(2*1:20)
    #k=1:40
    range <- range(est$est_fz1, estB$est_fz1, data$true_fz1)
    persp(zgrid(data$z1)[k], tgrid[k], data$true_fz1[k,k], main=expression(g[1](z[1],t)), 
          xlab="z", ylab="t", zlab="", zlim=range,  
          theta=30, phi=30, ticktype="detailed")
    persp(zgrid(data$z1)[k], tgrid[k], est$est_fz1[k,k], 
          xlab="z", ylab="t", zlab="", zlim=range, theta=30, phi=30, ticktype="detailed",
          main=bquote(paste(hat(g[1])(z[1],t), ": reliMSE", phantom(x)%~~% .(round(err$relmsefz1, 4)))))
    if(!is.null(errB)){
      persp(zgrid(data$z1)[k], tgrid[k], estB$est_fz1[k,k], 
            xlab="z", ylab="t", zlab="", zlim=range, theta=30, phi=30, ticktype="detailed",
            main=bquote(paste(hat(g[1])(z[1],t), ": reliMSE", phantom(x)%~~% .(round(errB$relmsefz1, 4)))))      
    }
  }
    
  # Plot \beta_1(s,t)
  if(theseSettings$scenario %in% c("1","3")){    
    tgrid <- seq(0,1, l=40)
    sgrid <- seq(0,1, l=40)    
    # Only plot part of the points -> grid is visible
    k=(2*1:20)
    #k=1:40  
    range <- range(est$est_bst1, estB$est_bst1, est$true_bst1)
    persp(tgrid[k], sgrid[k], z=est$true_bst1[k,k], theta = 30, phi = 30, 
          ticktype="detailed", zlim=range, xlab="s", ylab="t", zlab="",
          main=bquote(paste(beta[1](s,t))))    
    persp(tgrid[k], sgrid[k], est$est_bst1[k,k], theta = 30, phi = 30, 
          ticktype="detailed", zlim=range, xlab="s", ylab="t", zlab="",
          main=bquote(paste(hat(beta[1])(s,t), ": reliMSE",  phantom(x)%~~% .(round(err$relmsefx1, 4)))))    
    if(!is.null(errB)){
      persp(tgrid[k], sgrid[k], estB$est_bst1[k,k], theta=30, phi=30,    
            ticktype="detailed", zlim=range, xlab="s", ylab="t", zlab="",
            main=bquote(paste(hat(beta[1])(s,t), ": reliMSE",  phantom(x)%~~% .(round(errB$relmsefx1, 4)))))      
    }    
  }
    
  # Plot g(z)
  if(theseSettings$scenario %in% c("2")){
    range <- range(est$est_fz2, estB$est_fz2, est$true_fz2)
    plot(est$true_fz2~zgrid(data$z2, l=100), main=expression(g[2](z[2])), xlab="z", ylab="", ylim=range, type="l", xaxt="n")
    plot(est$est_fz2~zgrid(data$z2, l=100), xlab="z", ylab="", ylim=range, type="l", xaxt="n",
         main=bquote(paste(hat(g[2])(z[2]), ": reliMSE", phantom(x)%~~% .(round(err$relmsefz2, 4)))))
    if(!is.null(errB)){
      plot(estB$est_fz2~zgrid(data$z2, l=100), xlab="z", ylab="", ylim=range, type="l", xaxt="n",
           main=bquote(paste(hat(g[2])(z[2]), ": reliMSE", phantom(x)%~~% .(round(errB$relmsefz2, 4)))))      
    }
  }
  
  # Plot \beta(t)
  if(theseSettings$scenario %in% c("2")){
    tgrid100 <- seq(0,1, l=100)
    range <- range(est$est_fz3, estB$est_fz3, est$true_fz3)
    plot(est$true_fz3~tgrid100, main=expression(beta[3](t)), xlab="t", ylab="", ylim=range, type="l", xaxt="n")
    plot(est$est_fz3~tgrid100, xlab="t", ylab="", ylim=range, type="l", xaxt="n",
          main=bquote(paste(hat(beta[3])(t), ": reliMSE", phantom(x)%~~% .(round(err$relmsefz3, 4)))))
    if(!is.null(errB)){
      plot(estB$est_fz3~tgrid100, xlab="t", ylab="", ylim=range, type="l", xaxt="n",
           main=bquote(paste(hat(beta[3])(t), ": reliMSE", phantom(x)%~~% .(round(errB$relmsefz3, 4)))))      
    }

    range <- range(est$est_fz4, estB$est_fz4, est$true_fz4)
    plot(est$true_fz4~tgrid100, main=expression(beta[4](t)), xlab="t", ylab="", ylim=range, type="l", xaxt="n")
    plot(est$est_fz4~tgrid100, xlab="t", ylab="", ylim=range, type="l", xaxt="n",
         main=bquote(paste(hat(beta[4])(t), ": reliMSE", phantom(x)%~~% .(round(err$relmsefz4, 4)))))
    if(!is.null(errB)){
      plot(estB$est_fz4~tgrid100, xlab="t", ylab="", ylim=range, type="l", xaxt="n",
           main=bquote(paste(hat(beta[4])(t), ": reliMSE", phantom(x)%~~% .(round(errB$relmsefz4, 4)))))      
    }
  }

  # Plot \beta_2(s,t)
  if(theseSettings$scenario %in% c("3")){    
    tgrid <- seq(0,1, l=40)
    sgrid <- seq(0,1, l=40)    
    # Only plot part of the points -> grid is visible
    k=(2*1:20)
    #k=1:40 
    
    range <- range(est$est_bst2, estB$est_bst2, est$true_bst2)
    persp(tgrid[k], sgrid[k], z=est$true_bst2[k,k], theta = 30, phi = 30, 
          ticktype="detailed", zlim=range, xlab="s", ylab="t", zlab="",
          main=bquote(paste(beta[2](s,t))))    
    persp(tgrid[k], sgrid[k], est$est_bst2[k,k], theta = 30, phi = 30, 
          ticktype="detailed", zlim=range, xlab="s", ylab="t", zlab="",
          main=bquote(paste(hat(beta[2])(s,t), ": reliMSE",  phantom(x)%~~% .(round(err$relmsefx2, 4)))))    
    if(!is.null(errB)){
      persp(tgrid[k], sgrid[k], estB$est_bst2[k,k], theta=30, phi=30,    
            ticktype="detailed", zlim=range, xlab="s", ylab="t", zlab="",
            main=bquote(paste(hat(beta[2])(s,t), ": reliMSE",  phantom(x)%~~% .(round(errB$relmsefx2, 4)))))      
    }
  }
    
  # Plot Ytrue, reliMSE(Yfitted), Y and Ytrue-Yfitted 
  ylim <- range(est$Ytrue, est$Yhat, data$Y)
  matlplot(t(data$Ytrue[subset,]), col=clrs, xlab="", ylab="", xaxt="n", 
           main=expression(E(Y[i](t))), ylim=ylim, lty=1:length(subset))
  matlplot(t(est$yhat[subset,]), col=clrs, xlab="", ylab="", xaxt="n", ylim=ylim, lty=1:length(subset),
           main=bquote(paste(hat(Y[i])(t), ": reliMSE",  phantom(x)%~~%.(round(err$relmsey, 4))))) 
  if(!is.null(errB)) matlplot(t(estB$yhat[subset,]), col=clrs, xlab="", ylab="", 
                              xaxt="n", ylim=ylim, lty=1:length(subset),
                            main=bquote(paste(hat(Y[i])(t), ": reliMSE",  phantom(x)%~~%.(round(errB$relmsey, 4)))))

  matlplot(t(data$Y[subset,]), lwd=.5, lty=1:length(subset), col=clrs, xlab="", ylab="", xaxt="n", 
           main=expression(Y[i](t)))
  
  if(!is.null(errB)){
    range <- range(data$Ytrue - est$yhat, data$Ytrue - estB$yhat)
  }else{
    range <- range(data$Ytrue - est$yhat)
  }                      
  matlplot(t((data$Ytrue - est$yhat)[subset,]),  lwd=.5, lty=1:length(subset), col=clrs, xlab="", 
           ylab="", xaxt="n", main=expression(E(Y[i](t)) - hat(Y[i])(t)), ylim=range)
  segments(x0=1, y0=0, x1 = ncol(data$Y)+1, col="grey")
  if(!is.null(errB)) matlplot(t((data$Ytrue - estB$yhat)[subset,]), lwd=.5, lty=1:length(subset), col=clrs, xlab="", 
                            ylab="", xaxt="n", main=expression(E(Y[i](t)) - hat(Y[i])(t)), ylim=range)
  segments(x0=1, y0=0, x1 = ncol(data$Y)+1, col="grey")

  par(xpd=FALSE, mar=par()$mar*2)

}


#M=10; ni=1; Gy=30; Gx=20; snrEps=1; snrE=0; snrB=1; balanced=TRUE
if(FALSE){
  str(data1 <- makeData(scenario=1, seed=123))
  summary(m1 <- fitModelPffr(data1), freq=FALSE)
  #plot(m1, pers=TRUE, pages=1)
  summary(m1m <- fitModelMboost(data1, m_max=300), freq=FALSE)
  e1 <- getErrors(data1, m1, plotModel=TRUE)
  e1m <- getErrors(data1, m1m, plotModel=TRUE)  
  plotModel(err=e1, errB=e1m, data=data1, theseSettings=list(M=50, Gy=30, Gx=30, snrEps=1, scenario="1"))    
  
  str(data2 <- makeData(scenario=2, seed=123))
  summary(m2 <- fitModelPffr(data2), freq=FALSE)
  #plot(m2, pers=TRUE, pages=1)
  summary(m2m <- fitModelMboost(data2, m_max=300), freq=FALSE)
  e2 <- getErrors(data2, m2, plotModel=TRUE)
  e2m <- getErrors(data2, m2m, plotModel=TRUE)  
  plotModel(err=e2, errB=e2m, data=data2, theseSettings=list(M=50, Gy=30, Gx=30, snrEps=1, scenario="2"))    
  
  str(data3 <- makeData(scenario=3, seed=123))
  summary(m3 <- fitModelPffr(data3), freq=FALSE)
  #plot(m3, pers=TRUE, pages=1)
  summary(m3m <- fitModelMboost(data3, m_max = 250), freq=FALSE)
  e3 <- getErrors(data3, m3, plotModel=TRUE)
  e3m <- getErrors(data3, m3m, plotModel=TRUE)  
  plotModel(err=e3, errB=e3m, data=data3, theseSettings=list(M=50, Gy=30, Gx=30, snrEps=1, scenario="3"))    
     
}


##############################################################################################


### one replication of simulation to plot the models of boosting and pffr
# theseSettings=settings[[2]]
# theseSettings=settingsSplit[[9]][[1]]
oneRep <- function(theseSettings){
  #browser()
  #print(data.frame(theseSettings))
    
  # Generate data
  set.seed(theseSettings$seed)
  data <- do.call(makeData, theseSettings)
  args <- theseSettings
  args$data <- data
  
  #data5 <- makeData(scenario="5")
  
  print(unlist(theseSettings)[c(1,3,5,8,10,13)])
    
  # Fit models
  #modPffr <- NULL
  modPffr <- suppressMessages(do.call(fitModelPffr, args))
  modMboost <- suppressMessages(do.call(fitModelMboost, args))
  
  if(is.null(modPffr)) print(paste("modPffr, set " , theseSettings$set, ", is NULL", sep=""))
  if(is.null(modMboost)) print(paste("modMboost, set " , theseSettings$set, ", is NULL", sep=""))
  
  err <- getErrors(data, modPffr, plotModel=TRUE)
  errB <- getErrors(data, modMboost, plotModel=TRUE)
  
  # Save models of first rep
  if(theseSettings$rep==1){     
    nm <- paste(names(theseSettings)[c(1,3,4,5,7,8,10)], 
                theseSettings[c(1,3,4,5,7,8,10)], sep="",collapse="_")    
    pdf(paste(pathModels, nm, ".pdf", sep=""), width=8, height=12)
    try(plotModel(err=err, errB=errB, data=data, theseSettings=theseSettings))
    dev.off()   
  }
    
  # Calculate errors for both models 
  resMboost <- try(c(model=2, getErrors(data, modMboost)))
  resPffr <- try(c(model=1, getErrors(data, modPffr)))
  #cbind(resMboost, resPffr)
  
  rm(modMboost); rm(modPffr); rm(data); rm(args)
    
  res <- rbind(do.call(data.frame, c(theseSettings, resPffr)), 
        do.call(data.frame, c(theseSettings,resMboost)))
  
  return(res)
}

#  test <- oneRep(theseSettings)


### one replication of simulation: only fit the models with FDboost
# theseSettings=settings[[4]]
# theseSettings=settingsSplit[[9]][[1]]
oneRepFDboost <- function(theseSettings){
  #browser()
  #print(data.frame(theseSettings))
  
  # Generate data
  set.seed(theseSettings$seed)
  data <- do.call(makeData, theseSettings)
  args <- theseSettings
  args$data <- data
  
  print(unlist(theseSettings)[c(1,3,5,8,10,13)])
  
  # Fit models
  modPffr <- NULL
  ###  modPffr <- suppressMessages(do.call(fitModelPffr, args))
  
  args$control <- boost_control(mstop=2000, nu=0.1)
  args$grid <- seq(10, 2000, by = 10)
  args$m_max <- 2000
  
  modMboost <- suppressMessages(do.call(fitModelMboost, args)) # args$m_max=100
    
  if(is.null(modMboost)) print(paste("modMboost, set " , theseSettings$set, ", is NULL", sep=""))
    
  # Calculate errors for both models 
  resMboost <- try(c(model=2, getErrors(data, modMboost)))
  resPffr <- try(c(model=1, getErrors(data, modPffr)))
  #cbind(resMboost, resPffr)
  
  rm(modMboost); rm(modPffr); rm(data); rm(args)
    
  res <- rbind(do.call(data.frame, c(theseSettings, resMboost)))
  
  return(res)
}

#  test <- oneRep(theseSettings)


#s# one replication of simulation for pffr
# theseSettings=settings[[1]]
# theseSettings=settingsSplit[[9]][[1]]
oneRepPffr <- function(theseSettings){
  
  # Generate data
  set.seed(theseSettings$seed)
  data <- do.call(makeData, theseSettings)
  args <- theseSettings
  args$data <- data
    
  print(unlist(theseSettings)[c(1,3,5,8,10,13)])
    
  # Fit models
  modPffr <- suppressMessages(try(do.call(fitModelPffr, args)))
  
  # Should not be necessary!
  if(any(class(modPffr)=="try-error")){
    warning(paste("fitModelPffr, set " , theseSettings$set, ", gives error!", sep=""))
    modPffr <- NULL
  }
  modMboost <- NULL
  
  if(is.null(modPffr)) print(paste("modPffr, set " , theseSettings$set, ", is NULL", sep=""))
    
  # Calculate errors for both models 
  resMboost <- try(c(model=2, getErrors(data, modMboost)))
  resPffr <- try(c(model=1, getErrors(data, modPffr)))
  #cbind(resMboost, resPffr)
  
  #print(paste("after getErrors"))
  
  rm(modMboost); rm(modPffr); rm(data); rm(args)

  res <- rbind(do.call(data.frame, c(theseSettings, resPffr)))
    
  return(res)
}



### one replication of simulation: fit model with FDboost without optimizing mstop!
### to compare time with and without array model 
# theseSettings=settings[[4]]
# theseSettings=settingsSplit[[9]][[1]]
oneRepTime <- function(theseSettings){
  #browser()
  #print(data.frame(theseSettings))
  
  # Generate data
  set.seed(theseSettings$seed)
  data <- do.call(makeData, theseSettings)
  args <- theseSettings
  args$data <- data
  
  print(unlist(theseSettings)[c(1,3,5,8,10,13)])
  
  # Fit models
  modPffr <- NULL
  ###  modPffr <- suppressMessages(do.call(fitModelPffr, args))
  
  # do 1000 boosting iterations for model with array structure
  args$control <- boost_control(mstop = 1000, nu = 0.1)
  # no further search for the optimal mstop
  args$optimizeMstop <- FALSE
  modMboost <- suppressMessages(do.call(fitModelMboost, args)) 
  
  # do 1000 boosting iterations for model with array structure
  args$useArray <- FALSE
  modMboost2 <- suppressMessages(do.call(fitModelMboost, args))  
  
  if(is.null(modMboost)) print(paste("modMboost, set " , theseSettings$set, ", is NULL", sep=""))
  
  # Calculate errors for both models 
  resMboost <- try(c(model=2, time=modMboost$runTime, memory=modMboost$memory))
  resMboost2 <- try(c(model=3, time=modMboost2$runTime, memory=modMboost2$memory))
  # cbind(resMboost, resMboost2)
  
  rm(modMboost); rm(modMboost2); rm(data); rm(args)
  
  res <- rbind(do.call(data.frame, c(theseSettings, resMboost)), 
               do.call(data.frame, c(theseSettings, resMboost2)))
  
  return(res)
}

###########################################################################################

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

# make a list containing the settings
makeSettings <- function(...){
  settings <- expandList(...)
  seeds <- sample(1e5:3e6, length(settings))
  for(i in 1:length(settings)){
    settings[[i]]$seed <- seeds[i]
    settings[[i]]$set <- i
  }
  return(settings)
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
#   library(plyr)
#   
#   # only works on Linux -> with try() no error on windows
#   try(library(doMC))
#   try(registerDoMC(cores=cores))
  
  split <- sample(rep(1:cores, length=length(settings)))
  
  settingsSplit <- alply(1:cores, 1, function(i) {
    settings[split==i]   
  })
  
  ret <- ldply(settingsSplit, function(settings) {
    return(ldply(settings, oneRep, .parallel=FALSE))
  }, .parallel=TRUE)
  
  return(ret)    
}


# Do the iterations over oneRepPffr()
doSimPffr <- function(settings, cores=40){
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

# Do the iterations over fun(), use e.g. oneRepFDboost() or oneRepTime()
# slightly modified doSafeSim() 
doSimFDboost <- function(settings, fun){  #, savefile
  library(plyr)
  
  ret <- data.frame()
  failed <- list()
  
  for(s in 1:length(settings)){
    res <- try(do.call(fun, settings[s]), silent = TRUE)
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



