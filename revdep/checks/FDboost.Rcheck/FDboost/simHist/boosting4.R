###############################################################################
# changed code of Fabian Scheipl lfpr3.R for data generation and fit of pffr
# Author: Sarah Brockhaus
###############################################################################


# source(paste0(getwd(),"/Code/simUtils.R")) 
source("simUtils.R")

### Generate data of different scenarios
### save dataset, formula and some extra information

# M=12; ni=1; p=0.5; Gy=21; Gx=31; snrEps=5; snrE=0; snrB=1; balanced=TRUE;
# scenario="2"; centerX=TRUE; k=5; type="bsplines"; addNoise=FALSE; nuisance=0
# a="pen1coef4"; regularS=TRUE; regularT=TRUE; seed=123

makeData <- function(M=25, # number of subjects
        ni=1, # mean number of observations per subject
        p=0.5, # percentage of observations per curve
        Gy=30, # maximal number of grid points for t 
        Gx=55, # number of grid points for s
        snrEps=2, # signal-to-noise ratio
        snrE=0, # smooth errors per curve E_ij? - currently not used
        snrB=1, # relative importance of random effects
        balanced=TRUE, # not used
        scenario="2",
        centerX=TRUE, # center the functional predictor in each s?
        k=5, # number of spline bases to simulate the functional variables, c.f. lookatX1.R
        type="bsplines", # type of splines, e.g. "local", "start", "end", ...
        addNoise=FALSE, # add a small amount of noise to the covariates?
        nuisance=0, # number of nuisance variables
        a="pen1coef4", # function for coefficient surface 
        regularS=TRUE, regularT=TRUE, # are s and t regular?
        seed=123,
        ...){
  
  # set df in base learners
  df <- 2.24  ## round(sqrt(5), 2)
  
  dgp1 <- function(){
    stop("No scenario 1 implemented!")
  }
  
  # Scenario with TWO functional covariates
  dgp2 <- function(){
    
    idvar <- if(balanced){
      gl(M, ni)  # generate factor levels
    }  else {
      factor(c(1:M, sample(1:M, M*ni-M, repl=TRUE, prob=sqrt(1:M))))
    }
    
    # generate time variables with domain 1, ..., 16
    if(regularS){
      s <- seq(0, 1, l=Gx)*15 + 1
    }else{
      s <- (1:Gx-1)^2/(Gx-1)^2*15 + 1
    }
    
    if(regularT){
      t <- seq(0, 1, l=Gy)*15 + 1
    }else{
      t <- (1:Gy-1)^2/(Gy-1)^2*15 + 1
    }
    
    tgrid <- sgrid <- seq(1, 16, l=40)
    
    # Intercept
    int <- matrix(intf1(t), nrow=ni*M, ncol=Gy, byrow=TRUE)
    
    # Functional covariates
    #X1 <- t(replicate(M*ni, (3*s-2))) + rnorm(M*ni) + rnorm(M*ni*Gx, mean=0, sd=0.01)
    #set.seed(seed + 23) 
    X1 <- t(replicate(M*ni, rf2(s, k=k, type=type)))
    if(centerX) X1 <- sweep(X1, 2, apply(X1, 2, mean))
    L <- integrationWeightsLeft(X1=X1, xind=s)  # Riemann integration weights      
    #betast_1 <- outer(s, t, test3, a1=a1, a2=a2, a3=a3, a4=a4)
    betast_1 <- get(a)(s, t, coef=NULL)
    X1f <- (L*X1)%*%betast_1
    if(addNoise){
      X1 <- X1 + sd(as.vector(X1))/20*matrix(scale(rnorm(M*ni*Gx)), nrow=M*ni, ncol=Gx) 
    }
    
    #X2 <- t(replicate(M*ni, (rnorm(1, sd=2)*s-2))) + rnorm(M*ni, sd=0.1) + rnorm(M*ni*Gx, mean=0, sd=0.01)
    #set.seed(seed + 2233) 
    X2 <- t(replicate(M*ni, rf2(s, k=k, type=type)))
    if(centerX) X2 <- sweep(X2, 2, apply(X2, 2, mean))
    L <- integrationWeightsLeft(X1=X2, xind=s)  # Riemann integration weights      
    #betast_2 <- outer(s, t, test3, a1=a1, a2=a2, a3=a3, a4=a4)
    betast_2 <- get(a)(s, t, coef=NULL)
    X2f <- (L*X2)%*%betast_2
    if(addNoise){
      X2 <- X2 + sd(as.vector(X2))/20*matrix(scale(rnorm(M*ni*Gx)), nrow=M*ni, ncol=Gx) 
    }
    
    # Compute response
    Ytrue <- int + X1f + X2f     
    Y <- Ytrue + sd(as.vector(Ytrue))/snrEps * matrix(scale(rnorm(M*ni*Gy)), 
                                                      nrow=M*ni, ncol=Gy)      
    true_g0 <- intf1(tgrid)
    true_bst1 <- get(a)(sgrid, tgrid, coef=attr(betast_1, "coef"))
    true_bst2 <- get(a)(sgrid, tgrid, coef=attr(betast_2, "coef"))
    #true_bst1 <- outer(sgrid, tgrid, test3, a1=a1, a2=a2, a3=a3, a4=a4)
    #true_bst2 <- outer(sgrid, tgrid, test3, a1=a1, a2=a2, a3=a3, a4=a4)
    
    ## response in long format
    temp <- data.frame(Ylong=as.vector(t(Y)), Ytruelong=as.vector(t(Ytrue)),
                       X1flong=as.vector(t(X1f)), X2flong=as.vector(t(X2f)),
                       tlong=rep(t, times=M*ni), id=rep(1:(M*ni), each=Gy))
    
    # delete part of observations to obtain irregular observations in long format
    if(p<1){
      ########## delete proportion p of the observations
      #with(temp, funplot(tvec, Yvec, id))
      temp <- temp[ sort( sample(1:(M*ni*Gy), round(p*M*ni*Gy)) ), ]
      #with(temp, funplot(tvec, Yvec, id))
      #table(temp$id)
    }
    

    data <- list(Ylong=temp$Ylong, Ytruelong=temp$Ytruelong, tlong=temp$tlong, id=temp$id, 
                 idvar=idvar, 
                 Y=Y, X1=X1, X2=X2, 
                 s=s, t=t,
                 Ytrue=I(Ytrue), int=int, X1f=X1f, X2f=X2f, 
                 true_g0=true_g0, true_bst1=true_bst1, true_bst2=true_bst2,
                 X1flong=temp$X1flong, X2flong=temp$X2flong) 
    
    ## add nuisance variables
    if(nuisance>0){
      nuVars <- vector("list", nuisance)
      for(i in 1:nuisance){
        nuVars[[i]] <- t(replicate(M*ni, rf2(s, k=k, type=type)))
        if(centerX) nuVars[[i]] <- sweep(nuVars[[i]], 2, apply(nuVars[[i]], 2, mean))
        if(addNoise){
          nuVars[[i]] <- nuVars[[i]] + sd(as.vector(nuVars[[i]]))/
            20*matrix(scale(rnorm(M*ni*Gx)), nrow=M*ni, ncol=Gx) 
        }
      }
      names(nuVars) <- paste0("X", 1:nuisance+2)
      
      data <- c(data, nuVars)
    }
    
    return(data)
  }
    
  #### generate data
  data <- switch(scenario,
                 "1" = dgp1(),
                 "2" = dgp2())
  
  formulaPffr <- switch(scenario, 
                        "1" = "dummy",
                        "2" = paste0("ff(X", 1:(nuisance+2),", xind=s, limits=\"s<=t\", integration=\"riemann\", check.ident=FALSE, splinepars=list(bs = penaltyS, k=9, m = list(c(2, diffPen), c(2, diffPen))))", 
                                     collapse=" + "))
  
  # add response and do formula  
  formulaPffr <- as.formula(paste("Y ~ ", formulaPffr))  
  
 
  ### Formula for call to FDboost()  
  # arginS <- paste("\"smooth\"")
  formulaFDboost <- switch(scenario,
                           "1" = "dummy",
                           "2" = paste0("bhist(X", 1:(nuisance+2),", s, t, df=", df*df, ", knots=5, inS=arginS, penalty=penaltyS, differences = diffPen, check.ident=FALSE)", 
                                        collapse=" + "))
  formulaFDboost <- as.formula(paste("Y ~ 1 + ", formulaFDboost)) 
  
  
  timeformula <- formula(paste("~ bbs(t, knots=5, df=", df, ")", sep=""))
  
  ### Formula for call to FDboost() in long format  
  formulaFDboostLong <- switch(scenario,
                           "1" = "dummy",
                           "2" = paste0("bhist(X", 1:(nuisance+2),", s, tlong, df=", df*df, ", knots=5, inS=arginS, penalty=penaltyS, differences = diffPen, check.ident=FALSE)", 
                                        collapse=" + "))
  
  formulaFDboostLong <- as.formula(paste("Ylong ~ 1 + ", formulaFDboostLong)) 
  
  
  timeformulaLong <- formula(paste("~ bbs(tlong, knots=5, df=", df, ", differences = 1)", sep=""))
  
  namesVariables <- switch(scenario,
                             "1" = NULL,
                             "2" = c("Int", paste("X",  1:(nuisance+2), sep="")))
                           
  namesVariables <- c(namesVariables, attr(data, "namesVariables"))
      
  return(structure(data,
                   sigmaEps = sd(as.vector(data$Ytrue))/snrEps,
                   formulaPffr=formulaPffr,
                   formulaFDboost=formulaFDboost,
                   timeformula=timeformula,
                   formulaFDboostLong=formulaFDboostLong,
                   timeformulaLong=timeformulaLong,
                   namesVariables=namesVariables,
                   call=match.call()))
}

if(FALSE){
    data2 <- makeData(scenario="2", seed=14)  # , M=100, Gx=30, Gy=30
}


### Fit model using pffr() of package refund using the data and information of data from makeData()
fitModelPffr <- function(data, 
        ...){ 
      
  formula <- as.formula(attr(data, "formulaPffr"))  
  t <- data$t
  
  time <- system.time(
     m <- try(pffr(eval(formula), yind=t, data=data))
  )[3]
  
  rm(data)

  if(any(class(m) != "try-error")){
    m$runTime <- time 
    m$long <- FALSE
    return(m)
  }else return(NULL)
}

if(FALSE){
    m2 <- fitModelPffr(data2)
}


### Fit irregular model using pffr() 
fitModelPffrLong <- function(data, 
                         ...){ 
    
  formula <- as.formula(attr(data, "formulaPffr"))
  
  # save necessary data
  ydata <- data.frame(.obs=data$id, .index=data$tlong, .value=data$Ylong)
  
  data1 <-  data.frame(lapply(data[c(attr(data, "namesVariables")[-1])], I))
  tlong <<- data$tlong
  s <<- data$s
  
  time <- system.time(
    m <- try(pffr(eval(formula), yind=tlong, data=data1, ydata=ydata))  
    )[3]
  
  rm(data1)
  
  if(any(class(m) != "try-error")){
    m$runTime <- time 
    m$long <- TRUE
    return(m)
  }else return(NULL)
}

if(FALSE){
  m2 <- fitModelPffrLong(data2)
}


### Fit model using FDboost() that is based on package mboost
fitModelMboost <- function(data, 
                           control=boost_control(mstop=100, nu=0.1), # settings of mboost
                           grid=seq(10, 100, by=10),
                           m_max=2000, 
                           nuisance=0, 
                           ...){ 
     
  formula <- attr(data, "formulaFDboost")
  timeformula <- attr(data, "timeformula")
    
  time <- system.time({
    m <- try(FDboost(eval(formula), timeformula = timeformula, data=data, 
                     control=control, numInt="Riemann", 
                     offset=NULL, offset_control=o_control(k_min=10)))
    
    
    ##################################################################
    # do stability selection
    if(doStabsel && nuisance>0 & length(m$baselearner) > 5){
      
      mAll <- m[1] # save model with all baselearners
      
      # fix the cutoff at 0.9 and the PFER at 1
      print(stabsel_parameters(p=length(m$baselearner), PFER=0.1*length(m$baselearner), cutoff=0.9, 
                               sampling.type="SS"))
      
      folds0 <- cvLong(id=m$id, weights=model.weights(m), B=50, type="subsampling")
      
      m[200]
      #table(selected(m))
      
      if(Sys.info()["sysname"]=="Linux"){        
        stab1 <- try(stabsel(m, cutoff=0.9, PFER=0.1*length(m$baselearner), 
                             folds=folds0, sampling.type="SS", mc.cores=10)) 
      }else{ # Windows
        stab1 <- try(stabsel(m, cutoff=0.9, PFER=0.1*length(m$baselearner), 
                             folds=folds0, sampling.type="SS")) 
      }
      
      #if(grepl("stabsel", as.character(warnings()))) print("stab1: INCREASE mstop!!!")
      
      print(stab1$selected)
      
      ## Function to obtain effects of original formula
      shortnames2 <- function(x){
        if(substr(x,1,1)=="\"") x <- substr(x, 2, nchar(x)-1)
        
        sign <- "%O%"
        if( !grepl(sign, x) ) sign <- "%X%"
        
        xpart <- unlist(strsplit(x, sign)) 
        for(i in 1:length(xpart)){
          xpart[i] <- gsub("\\\"", "'", xpart[i], fixed=TRUE)
          commaSep <- unlist(strsplit(xpart[i], ","))
          xpart[i] <- paste(paste(commaSep[!grepl("index", commaSep)], collapse=","), sep="")
        }
        ret <- xpart
        if(length(xpart)>1) ret <- xpart[1]
        ret
      }
      
      # create a new formula containing only selected effects
      newForm <- formula(paste("Y ~ 1 + ", paste(lapply(names(m$baselearner)[stab1$selected], shortnames2), 
                                                 collapse=" + ")))
      print(newForm)
      
      # fit the model with the effects selected by stabsel
      m <- try(FDboost(newForm, timeformula=timeformula, 
                       data=data, control=control, numInt="Riemann", 
                       offset=NULL, offset_control=o_control(k_min=10)))    
            
      #summary(m)
      attr(m, "mAll") <- mAll
    }
    
    
    ################################################################## 
    # search for optimal stopping iteration
    if(any(class(m)=="FDboost")){
      
      # set up two different splittings of the data  
      set.seed(attr(data, "call")$seed)
      folds1 <- cvMa(ydim=m$ydim, type="bootstrap", B=10) 
      
      set.seed(attr(data, "call")$seed + 100)
      folds2 <- cvMa(ydim=m$ydim, type="bootstrap", B=10)
                  
      # cvm <- try(suppressMessages(cvrisk(m, papply = lapply, grid=grid)))
      if(Sys.info()["sysname"]=="Linux"){ # use 10 cores on Linux
        suppressWarnings(cvm <- try(cvrisk(m, folds = folds1, grid=grid, mc.cores=10), silent=TRUE))        
      }else cvm <- try(cvrisk(m, folds = folds1, grid=grid))          
      
      # Try for a second time if cvrisk() stops with error
      if(class(cvm)!="cvrisk"){
        if(Sys.info()["sysname"]=="Linux"){
          suppressWarnings(cvm <- try(cvrisk(m, folds = folds2, grid=grid, mc.cores=10), silent=TRUE))        
        }else cvm <- try(cvrisk(m, folds = folds2, grid=grid)) 
        cat(paste("2nd calculation of cvrisk, class(cvm) =", class(cvm), "\n", sep=" "))
      }
      # print(max(grid))
      # plot(cvm)    
    }else cvm <- NULL 
  })[3]
  
  if(any(class(m)=="FDboost") & class(cvm)=="cvrisk"){
    print(mstop(cvm))
    m <- m[mstop(cvm)]
    m$runTime <- time 
    m$long <- FALSE
    return(m)
  }else return(NULL)    
}

if(FALSE){
  m2m <- fitModelMboost(data2, m_max=500)
}


### Fit model to data in long format
fitModelMboostLong <- function(data, 
                           control=boost_control(mstop=100, nu=0.1), # settings of mboost
                           grid=seq(10, 100, by=10),
                           m_max=2500, 
                           nuisance=0,
                           ...){ 
  
  #scenario <- attr(data, "call")$scenario    
  formula <- attr(data, "formulaFDboostLong")
  timeformula <- attr(data, "timeformulaLong")
   
  #      m <- FDboost(Ylong ~ 1 + bhist(X1, s, tlong, df = 5.0625) + bhist(X2, s, tlong, df = 5.0625),
  #                   timeformula=~bbs(tlong, knots = 10, df = 2.25), 
  #                   id=data$id, data=data)
  
  time <- system.time({
    m <- try(FDboost(eval(formula), timeformula = timeformula, id=~id, 
                     data=data, control=control, numInt="Riemann", 
                     offset=NULL, offset_control=o_control(k_min=10)))
    
    
    ##################################################################
    # do stability selection
    if(doStabsel && nuisance>0 & length(m$baselearner) > 5){
      
      mAll <- m[1] # save model with all baselearners
      
      # fix the cutoff at 0.9 and the PFER at 1
      print(stabsel_parameters(p=length(m$baselearner), PFER=0.1*length(m$baselearner), cutoff=0.9, 
                               sampling.type="SS"))
      
      folds0 <- cvLong(id=m$id, weights=model.weights(m), B=50, type="subsampling")
      
      m[200]
      #table(selected(m))
      
      if(Sys.info()["sysname"]=="Linux"){        
        stab1 <- try(stabsel(m, cutoff=0.9, PFER=0.1*length(m$baselearner), 
                             folds=folds0, sampling.type="SS", mc.cores=10)) 
      }else{ # Windows
        stab1 <- try(stabsel(m, cutoff=0.9, PFER=0.1*length(m$baselearner), 
                             folds=folds0, sampling.type="SS")) 
      }
      
      #if(grepl("stabsel", as.character(warnings()))) print("stab1: INCREASE mstop!!!")
      print(stab1$selected)
      
      ## Function to obtain effects of original formula
      shortnames2 <- function(x){
        if(substr(x,1,1)=="\"") x <- substr(x, 2, nchar(x)-1)
        
        sign <- "%O%"
        if( !grepl(sign, x) ) sign <- "%X%"
        
        xpart <- unlist(strsplit(x, sign)) 
        for(i in 1:length(xpart)){
          xpart[i] <- gsub("\\\"", "'", xpart[i], fixed=TRUE)
          commaSep <- unlist(strsplit(xpart[i], ","))
          xpart[i] <- paste(paste(commaSep[!grepl("index", commaSep)], collapse=","), sep="")
        }
        ret <- xpart
        if(length(xpart)>1) ret <- xpart[1]
        ret
      }
      
      # create a new formula containing only selected effects
      newForm <- formula(paste("Ylong ~ 1 + ", 
                               paste(lapply(names(m$baselearner)[stab1$selected], shortnames2), 
                                                 collapse=" + ")))
      
      print(newForm)
      
      # fit the model with the effects selected by stabsel
      m <- try(FDboost(newForm, timeformula=timeformula, id=~id, 
                       data=data, control=control, numInt="Riemann", 
                       offset=NULL, offset_control=o_control(k_min=10)))
      #summary(m)
      attr(m, "mAll") <- mAll
    }
    
    
    ################################################################## 
    # search for optimal stopping iteration
    if(any(class(m)=="FDboost")){
      
      # set up two different splittings of the data
      set.seed(attr(data, "call")$seed) 
      folds1 <- cvLong(id=m$id, weights=model.weights(m), type="bootstrap", B=10)
      
      set.seed(attr(data, "call")$seed + 100)
      folds2 <- cvLong(id=m$id, weights=model.weights(m), type="bootstrap", B=10)
      
      # cvm <- try(suppressMessages(cvrisk(m, papply = lapply, grid=grid)))
      if(Sys.info()["sysname"]=="Linux"){ # use 10 cores on Linux
        # results are not the same as in validateFDboost the offset is refitted
        # in cvrisk the smooth offset of the model with all data is used
        suppressWarnings(cvm <- try(cvrisk(m, folds = folds1, grid=grid, mc.cores=10), silent=TRUE)) 
        #cvm2 <- try(validateFDboost(m, folds = folds1, grid=grid, 
        #                            getCoefCV=FALSE, mc.cores=10), silent=TRUE)          
      }else suppressWarnings(cvm <- try(cvrisk(m, folds = folds1, grid=grid)))        
      
      # Try for a second time if cvrisk() stops with error
      if(class(cvm)!="cvrisk"){
        if(Sys.info()["sysname"]=="Linux"){
          cvm <- try(cvrisk(m, folds = folds2, grid=grid, mc.cores=10), silent=TRUE)        
        }else cvm <- try(cvrisk(m, folds = folds2, grid=grid)) 
        cat(paste("2nd calculation of cvrisk, class(cvm) =", class(cvm), "\n", sep=" "))
      }
      
      #print(max(grid))
      # plot(cvm)      
    }else cvm <- NULL 
  })[3]
  
  if(any(class(m)=="FDboost") & class(cvm)=="cvrisk"){
    print(mstop(cvm))
    m <- m[mstop(cvm)]
    m$runTime <- time 
    m$long <- TRUE
    return(m)
  }else return(NULL)    
}

if(FALSE){
  m2m <- fitModelMboostLong(data2, m_max=200)
}



### calculate some errors/ measures for goodness of fit of the model 
getErrors <- function(data, m=NULL, plotModel=FALSE){
  
  scenario <- attr(data, "call")$scenario
  
  # In case that model was not fitted: m=NULL
  errorE <- c(msey=NA, mseg0=NA, 
              msefx1=NA, msefx2=NA, 
              msefx1f=NA, msefx2f=NA, 
              msefz1=NA
              ) 
  relerrorE <- errorE
  names(relerrorE) <- paste("rel", names(errorE), sep="")
  relerrorE[["relmsefx1b"]] <- NA
  relerrorE[["relmsefx2b"]] <- NA
  
  ret <- as.list(c(errorE, relerrorE, #relmseyT=NA, 
                   funRsqrt=NA, mstop=NA, time.elapsed=NA, long=NA))

  if(is.null(m)) return(ret)
    
  classM <- class(m)[1]
  
  # Save number of iterations
  mstop <- switch(classM,
                 "pffr" = NA,
                 "FDboost" = mstop(m),
                 "FDboostLong" = mstop(m),
                 "NULL" = NA)  
  
  
  # function to predict each component of linear predictor separately for FDboost() and pffr()
  # intercept =TRUE indicates that the first variable is an intercept 
  # to which the offset should be added 
  
  predictComponents <- function(object, intercept=TRUE){
    fit <- NULL
    if(any(class(object)=="FDboost")){
      fit <- predict(object=object, which=1:length(object$baselearner))
      if(any(class(object)=="FDboostLong")){ 
        if(intercept) fit[,1] <- fit[,1] + object$offset  # add offset 
      }else{
        if(intercept) fit[[1]] <- fit[[1]] + object$offset  # add offset 
      }
    }else{
      fit <- predict(object=object, type="terms")
      if(intercept) fit[[1]] <- fit[[1]] + coef(object, se=FALSE, seWithMean=FALSE)$pterms[1]  # add global constan intercept
    }
    return(fit)    
  }
  
  # Perdiciton for pffr or FDboost, each effect separately
  fit <- switch(classM,
                "pffr" = predictComponents(m),
                "FDboost" = predictComponents(m),
                "FDboostLong" = predictComponents(m),
                "NULL" = NULL) 
  
  # prediction of all effects together
  yhat <- switch(classM,
                 "pffr" = if(m$long) fitted(m)$.value else fitted(m),
                 "FDboost" = predict(m),
                 "FDboostLong" = predict(m),
                 "NULL" = NULL)  

  ############### Calculate errors of estimated coefficients directly on coefficients
  if(any(class(m)=="pffr")){    
    cm <- coef(m, se=FALSE, seWithMean=FALSE, n1=40)

    # functional intercept g0
    if(scenario %in% c("2")){
      est_g0 <- drop(cm$smterms[[1]]$value) + cm$pterms[1]
      true_g0 <- data$true_g0
    } else{
      est_g0 <- true_g0 <- NA
    }
    
    #plot(true_g0, col=2, ylim=range(true_g0, est_g0)); points(est_g0)
    #mean((est_g0-true_g0)^2)
        
    if(scenario %in% c("1")){
      true_fz1 <- data$true_fz1
      est_fz1 <- t(cbind(cm$smterms[[1]]$value + cm$pterms[1], 
                         cm$smterms[[2]]$value + cm$pterms[2], 
                         cm$smterms[[3]]$value + cm$pterms[3]))
    }else est_fz1 <- true_fz1 <- NA
    #     funplot(seq(0,1,l=40), (true_fz1), ylim=range(true_fz1, est_fz1))
    #     funplot(seq(0,1,l=40), (est_fz1), ylim=range(true_fz1, est_fz1), lwd=1.5, add=TRUE) 
    
    ### funciton to set lower triangular to 0 or NA
    lowerTo <- function(x, repl=0){
      stopifnot(ncol(x)==nrow(x))
      #x*outer(1:ncol(x), 1:nrow(x), "<=") # gives the same if repl=0
      x[ outer(1:ncol(x), 1:nrow(x), "<=")==FALSE] <- repl
      x
    }
    
    true_bst1 <- est_bst1 <- true_bst2 <- est_bst2 <- NA #<- true_bst3 <- est_bst3 <- true_bst4 <- est_bst4 <- NA        
    if(scenario %in% c("1", "2")){
      whereX1 <- if(scenario %in% c("1")) 4 else 2 # 1 intercept, 2-4 z1
      true_bst1 <- data$true_bst1
      est_bst1 <- lowerTo(matrix(drop(cm$smterms[[whereX1]]$value), ncol=40))   
      if(scenario %in% c("2")){
        true_bst2 <- data$true_bst2
        est_bst2 <- lowerTo(matrix(drop(cm$smterms[[3]]$value), ncol=40))      
      }     
    }
    #     par(mfrow=c(1,2))
    #     persp(seq(0,1,l=40), seq(0,1,l=40), true_bst1, zlim=range(true_bst1, est_bst1), theta=30, ticktype="detailed")
    #     persp(seq(0,1,l=40), seq(0,1,l=40), est_bst1, zlim=range(true_bst1, est_bst1), theta=30, ticktype="detailed")        
    #     persp(seq(0,1,l=40), seq(0,1,l=40), true_bst2, zlim=range(true_bst2, est_bst2), theta=30, ticktype="detailed")
    #     persp(seq(0,1,l=40), seq(0,1,l=40), est_bst2, zlim=range(true_bst2, est_bst2), theta=30, ticktype="detailed")        
    
    ## get estimated effects
    if(scenario %in% c("2")){ 
      if(get("p", environment(attr(data, "timeformula"))) < 1){ # irregualr response
        true_x1f <- data$X1flong
        true_x2f <- data$X2flong
        est_x1f <- fit[[2]]$.value
        est_x2f <- fit[[3]]$.value
      }else{ # regualr response
        true_x1f <- data$X1f
        true_x2f <- data$X2f
        est_x1f <- fit[[2]]
        est_x2f <- fit[[3]]
      }
    }else{
      true_x1f <- true_x2f <- est_x1f <- est_x2f <- NA
    }
    
  }

  ################ Calculate errors of estimated coefficients directly on coefficients
  if(any(class(m)=="FDboost")){
    
    ## compute coefficients
    cm <- coef(m, which=1:3)  
    
    if(scenario %in% c("2")){
      # functional intercept g0
      est_g0 <- drop(cm$smterms[[1]]$value) + cm$offset$value
      true_g0 <- data$true_g0
      # tgrid <- seq(0,1,l=40)
      # plot(true_g0 ~ tgrid, col=2, ylim=range(true_g0, est_g0)); points(est_g0~tgrid)
      # mean((est_g0-true_g0)^2)
    }else{
      est_g0 <- true_g0 <- NA
    }

    est_fz1 <- true_fz1 <- NA
    
    true_bst1 <- est_bst1 <- true_bst2 <- est_bst2 <- NA #<- true_bst3 <- est_bst3 <- true_bst4 <- est_bst4 <- NA        
    if(scenario == "2"){
      whereX1 <- if(scenario %in% c("1")) 3 else 2
      true_bst1 <- data$true_bst1
      est_bst1 <- matrix(cm$smterms[[whereX1]]$value, ncol=40, nrow=40)
      if(scenario %in% c("2")){
        true_bst2 <- data$true_bst2
        est_bst2 <- matrix(cm$smterms[[3]]$value, ncol=40, nrow=40)
        true_bst3 <- data$true_bst3     
      }     
    }
    # par(mfrow=c(1,2))
    # persp(seq(0,1,l=40), seq(0,1,l=40), true_bst1, zlim=range(true_bst1, est_bst1), theta=30, ticktype="detailed")
    # persp(seq(0,1,l=40), seq(0,1,l=40), est_bst1, zlim=range(true_bst1, est_bst1), theta=30, ticktype="detailed")        
    # persp(seq(0,1,l=40), seq(0,1,l=40), true_bst2, zlim=range(true_bst2, est_bst2), theta=30, ticktype="detailed")
    # persp(seq(0,1,l=40), seq(0,1,l=40), est_bst2, zlim=range(true_bst2, est_bst2), theta=30, ticktype="detailed")            
    
    ## get estimated effects
    if(scenario %in% c("2")){ 
      if(get("p", environment(attr(data, "timeformula"))) < 1){ # irregular response
        true_x1f <- data$X1flong
        true_x2f <- data$X2flong
        est_x1f <- fit[,2]
        est_x2f <- fit[,3]
      }else{ # regular response
        true_x1f <- data$X1f
        true_x2f <- data$X2f
        est_x1f <- fit[[2]]
        est_x2f <- fit[[3]]
      }
    }else{
      true_x1f <- true_x2f <- est_x1f <- est_x2f <- NA
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
  
  ## only look at part of coefficient surface that was really fitted:
  limits <- function(s, t) {
    (s < t) 
  }
  sgrid <- cm$smterms[[2]]$x
  tgrid <- cm$smterms[[2]]$y
  ind0 <- !t(outer( sgrid, tgrid, limits) )
  
  errorE <- c(msey = if(m$long) calcError(data$Ytruelong, yhat) else calcError(data$Ytrue, yhat),
                mseg0 = calcError(true_g0, est_g0),              
                msefx1 = calcError(true_bst1[ind0], est_bst1[ind0]),
                msefx2 = calcError(true_bst2[ind0], est_bst2[ind0]),
                msefx1f= calcError(true_x1f, est_x1f),
                msefx2f= calcError(true_x2f, est_x2f),
                msefz1 = calcError(true_fz1, est_fz1)) 
  
  ### look at the error for each place of beta
  diffSurface <- true_bst1 - est_bst1
    
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

  relerrorE <- c(msey = if(m$long) relcalcError(data$Ytruelong, yhat)  else relcalcError(data$Ytrue, yhat),
              mseg0 = relcalcError(true_g0, est_g0),              
              msefx1 = relcalcError(true_bst1[ind0], est_bst1[ind0]),
              msefx2 = relcalcError(true_bst2[ind0], est_bst2[ind0]),
              msefx1f= relcalcError(true_x1f, est_x1f),
              msefx2f= relcalcError(true_x2f, est_x2f),
              msefz1 = relcalcError(true_fz1, est_fz1))
  
  names(relerrorE) <- paste("rel", names(errorE), sep="")
  
  ### compute relativeMSE as in identifiability paper:
  relerrorE[["relmsefx1b"]] <- mean((est_bst1[ind0] - true_bst1[ind0])^2)/mean(true_bst1[ind0]^2)
  relerrorE[["relmsefx2b"]] <- mean((est_bst2[ind0] - true_bst2[ind0])^2)/mean(true_bst2[ind0]^2)
  
  
  # return list with errors and relative errors
  ret <- as.list(c(errorE, relerrorE, 
                   mstop=mstop, time=m$runTime, diffSurface = c(diffSurface)))
  
  # Save information for plotting
  if(plotModel==TRUE){
    est <- list(Ytrue=data$Ytrue, Ytruelong=data$Ytruelong, id=data$id, 
                yhat=yhat, 
                true_g0=true_g0, est_g0=est_g0,
                true_bst1=true_bst1, est_bst1=est_bst1, 
                true_bst2=true_bst2, est_bst2=est_bst2,
                true_fz1=true_fz1, est_fz1=est_fz1) 
    attr(ret, "est") <- est
  }
  
  ret$long <- if(m$long) TRUE else FALSE
          
  return(ret)
}

if(FALSE){
    eNULL <- getErrors(NULL, NULL)
    length(eNULL)
    e2 <- getErrors(data2, m2) #, plotModel=TRUE
    length(e2)
    #cbind(eNULL, e2)
}

#true_fz2, est_fz2

if(FALSE){
  e2m <- getErrors(data2, m2m)  #, plotModel=TRUE
  cbind(e2, e2m)
}


# function to plot true values, estimates of FDboost and estimates of pffr
# depends on results of function getErrors()
# use only a subset of 10 observations for plotting them as example
plotModel <- function(err, errB=NULL, data, theseSettings, 
                      subset=sample(1:theseSettings$M, size=10)){
  
  ### funciton to set lower triangular to 0 or NA
  lowerTo <- function(x, repl=0){
    stopifnot(ncol(x)==nrow(x))
    #x*outer(1:ncol(x), 1:nrow(x), "<=") # gives the same if repl=0
    x[ outer(1:ncol(x), 1:nrow(x), "<=")==FALSE] <- repl
    x
  }
  
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
  settingString <- paste(settingString[names(theseSettings) %in% c("a1","a2","a3","M","k")], collapse=" ")
  #opar <- par()
  #on.exit(try(par(opar), silent=TRUE))
  
  par(xpd=TRUE, mar=par()$mar/2)
  
  #clrs <- try(alpha(rainbow(theseSettings$M), 4 * 255/ (length(unique(data$id))/5) ), silent=TRUE )
  
  #clrs <- alpha(rainbow(theseSettings$M), 4*255/(min(theseSettings$M, 40)/5))
  clrs <- alpha( rainbow(length(subset)), 200 )
  if(class(clrs)=="try-error")  clrs <- rainbow(theseSettings$M)
  
  #  plotIntercept <- function(m=m, mb=mb, fit=fit, fitB=fitB, errors=errors, errorsB=errorsB, data=data){
  plotIntercept <- function(){
    range <- range(est$est_g0, estB$est_g0, data$true_g0)
    plot(data$true_g0, main=bquote(beta[0](t)), xlab="", ylab="", xaxt="n", ylim=range, type="l")
    if(!is.null(err)){
      plot(est$est_g0, xlab="", ylab="", xaxt="n", ylim=range, type="l", 
           main=bquote(paste("FAMM: ", hat(beta)[0](t), ": reliMSE", phantom(x)%~~% .(round(err$relmseg0, 4)))))      
    }
    if(!is.null(errB)){
      plot(estB$est_g0, xlab="", ylab="", xaxt="n", ylim=range, type="l",
           main=bquote(paste("FDboost: ", hat(beta)[0](t), ": reliMSE", phantom(x)%~~% .(round(errB$relmseg0, 4)))))
    }   
  }
  
  twoModels <- !( is.null(errB)|is.null(err) )

  layout(matrix(1:(10+5*twoModels), ncol=2+1*twoModels, byrow=TRUE))
  
  # plot the intercept
  if(theseSettings$scenario %in% c("2")) plotIntercept()
  
  # If pffr was not fitted but FDboost was, then change err and errB
  if(is.null(err) & !is.null(errB)){
    err <- errB
    est <- estB
    errB <- NULL
    estB <- NULL
    model <- modelB
  }
    
  # Plot \beta_1(s,t)
  if(theseSettings$scenario %in% c("1","2")){    
    tgrid <- seq(0,1, l=40)
    sgrid <- seq(0,1, l=40)    
    # Only plot part of the points -> grid is visible
    k=(2*1:20)
    #k=1:40  
    range <- range(lowerTo(est$est_bst1, NA), lowerTo(estB$est_bst1, NA), lowerTo(est$true_bst1, NA), na.rm=TRUE)
    persp(tgrid[k], sgrid[k], z=lowerTo(est$true_bst1[k,k], NA), theta = 30, phi = 30, 
          ticktype="detailed", zlim=range, xlab="s", ylab="t", zlab="",
          main=bquote(paste(beta[1](s,t))))    
    persp(tgrid[k], sgrid[k], lowerTo(est$est_bst1[k,k], NA), theta = 30, phi = 30, 
          ticktype="detailed", zlim=range, xlab="s", ylab="t", zlab="",
          main=bquote(paste(hat(beta)[1](s,t), ": reliMSE",  phantom(x)%~~% .(round(err$relmsefx1b, 4)))))    
    if(!is.null(errB)){
      persp(tgrid[k], sgrid[k], lowerTo(estB$est_bst1[k,k], NA), theta=30, phi=30,    
            ticktype="detailed", zlim=range, xlab="s", ylab="t", zlab="",
            main=bquote(paste(hat(beta)[1](s,t), ": reliMSE",  phantom(x)%~~% .(round(errB$relmsefx1b, 4)))))      
    }    
  }
  
  # Plot \beta_2(s,t), \beta_3(s,t), \beta_4(s,t)
  if(theseSettings$scenario %in% c("2")){    
    tgrid <- seq(0,1, l=40)
    sgrid <- seq(0,1, l=40)    
    # Only plot part of the points -> grid is visible
    k=(2*1:20)
    #k=1:40 
    range <- range(lowerTo(est$est_bst2, NA), lowerTo(estB$est_bst2, NA), lowerTo(est$true_bst2, NA), na.rm=TRUE)
    persp(tgrid[k], sgrid[k], z=lowerTo(est$true_bst2[k,k], NA), theta = 30, phi = 30, 
          ticktype="detailed", zlim=range, xlab="s", ylab="t", zlab="",
          main=bquote(paste(beta[2](s,t))))    
    persp(tgrid[k], sgrid[k], lowerTo(est$est_bst2[k,k],NA), theta = 30, phi = 30, 
          ticktype="detailed", zlim=range, xlab="s", ylab="t", zlab="",
          main=bquote(paste(hat(beta)[2](s,t), ": reliMSE",  phantom(x)%~~% .(round(err$relmsefx2b, 4)))))    
    if(!is.null(errB)){
      persp(tgrid[k], sgrid[k], lowerTo(estB$est_bst2[k,k], NA), theta=30, phi=30,    
            ticktype="detailed", zlim=range, xlab="s", ylab="t", zlab="",
            main=bquote(paste(hat(beta)[2](s,t), ": reliMSE",  phantom(x)%~~% .(round(errB$relmsefx2b, 4)))))      
    } 
  }

  # Plot Ytrue, reliMSE(Yfitted), Y and Ytrue-Yfitted 
  ylim <- range(est$Ytrue, est$Yhat, data$Y)
  subsetID <- data$id %in% subset
  if(err$long && errB$long){ ### plot the irregualr data
    funplot(data$tlong[subsetID], data$Ytruelong[subsetID], id=data$id[subsetID], 
            col=clrs, xlab="", ylab="", xaxt="n",
            pch="", main=expression(E(Y[i](t))), ylim=ylim)
  }else{
    funplot(data$t, data$Ytrue[subset,], col=clrs, xlab="", ylab="", xaxt="n", pch="", 
            main=expression(E(Y[i](t))), ylim=ylim)
  }
  
  if(err$long){
    plot(data$tlong[subsetID], est$yhat[subsetID], ylim=ylim, col="white", xlab="", ylab="", xaxt="n",
         main=bquote(paste(hat(Y)[i](t), ": reliMSE",  phantom(x)%~~%.(round(err$relmsey, 4)))))
    funplot(data$tlong[subsetID], est$yhat[subsetID], id=data$id[subsetID], 
            col=clrs, xlab="", ylab="", xaxt="n",
            pch="", main=expression( hat(Y)[i](t) ), ylim=ylim, add=TRUE)
  }else{
    matlplot(data$t, t(est$yhat[subset,]), col=clrs, xlab="", ylab="", xaxt="n", ylim=ylim,
             main=bquote(paste(hat(Y)[i](t), ": reliMSE",  phantom(x)%~~%.(round(err$relmsey, 4)))))
    rug(data$t, 0.01)
  }
  
  if(!is.null(errB)){
    if(errB$long){
      plot(data$tlong[subsetID], estB$yhat[subsetID], ylim=ylim, col="white", xlab="", ylab="", xaxt="n",
           main=bquote(paste(hat(Y)[i](t), ": reliMSE",  phantom(x)%~~%.(round(errB$relmsey, 4)))))
      funplot(data$tlong[subsetID], estB$yhat[subsetID], id=data$id[subsetID], 
              col=clrs, xlab="", ylab="", xaxt="n",
              pch="", main=expression( hat(Y)[i](t) ), ylim=ylim, add=TRUE)
    }else{
      matlplot(data$t, t(estB$yhat[subset,]), col=clrs, xlab="", ylab="", xaxt="n", ylim=ylim,
               main=bquote(paste(hat(Y)[i](t), ": reliMSE",  phantom(x)%~~%.(round(errB$relmsey, 4)))))
      rug(data$t, 0.01)
    }
  } 

  # plot data with error
  if(err$long){
    funplot(data$tlong[subsetID], data$Ylong[subsetID], id=data$id[subsetID], 
            col=clrs, xlab="", ylab="", xaxt="n",
            pch="", main=expression(Y[i](t)), ylim=ylim)
  }else{
    matlplot(data$t, t(data$Y[subset,]), lwd=.5, lty=1, col=clrs, xlab="", ylab="", xaxt="n", 
             main=expression(Y[i](t)))
    rug(data$t, 0.01)
  }

  ### plot residuals
  
  if(err$long){
    ylim <- range(data$Ytruelong - est$yhat, data$Ytruelong - estB$yhat)
    funplot(data$tlong[subsetID], (data$Ytruelong - est$yhat)[subsetID], id=data$id[subsetID], 
            col=clrs, xlab="", ylab="", xaxt="n",
            pch="", main=expression(E(Y[i](t)) - hat(Y)[i](t)), ylim=ylim)
    
    if(!is.null(errB)){
      funplot(data$tlong[subsetID], (data$Ytruelong - estB$yhat)[subsetID], id=data$id[subsetID], 
              col=clrs, xlab="", ylab="", xaxt="n",
              pch="", main=expression(E(Y[i](t)) - hat(Y)[i](t)), ylim=ylim)
    }
  }else{
    range <- range(data$Ytrue - est$yhat, data$Ytruelong - estB$yhat)
    matlplot(data$t, t((data$Ytrue - est$yhat)[subset,]),  lwd=.5, lty=1, col=clrs, xlab="", 
             ylab="", xaxt="n", main=expression(E(Y[i](t)) - hat(Y)[i](t)), ylim=range)
    
    if(!is.null(errB)) matlplot(t((data$Ytrue - estB$yhat)[subset,]),  lwd=.5, lty=1, col=clrs, xlab="", 
                                ylab="", xaxt="n", main=expression(E(Y[i](t)) - hat(Y)[i](t)), ylim=range)
  }
  
  ##title(sub=settingString, cex.sub=1.5, line=-2, outer=TRUE)
  par(xpd=FALSE, mar=par()$mar*2)
  
}    


### CEHCK
#M=10; ni=1; Gy=30; Gx=20; snrEps=1; snrE=0; snrB=1; balanced=TRUE
if(FALSE){
  
 ### models on all data - regular grid
  str(data2 <- makeData(scenario=2, seed=123)) #, nuisance=4
  summary(m2 <- fitModelPffr(data2), freq=FALSE)
  #plot(m2, pers=TRUE, pages=1)
  summary(m2m <- fitModelMboost(data2, m_max=300), freq=FALSE) #, nuisance=4
  e2 <- getErrors(data2, m2, plotModel=TRUE)
  e2m <- getErrors(data2, m2m, plotModel=TRUE)  
  # cbind(e2, e2m)
  plotModel(err=e2, errB=e2m, data=data2, theseSettings=list(M=25, Gy=30, Gx=35, snrEps=2, scenario="2"))    
  
  
  ### models on irregular data
  str(data2 <- makeData(scenario=2, seed=123, nuisance=4)) #
  summary(m2l <- fitModelPffrLong(data2), freq=FALSE)
  #plot(m2l, pers=TRUE, pages=1)
  summary(m2ml <- fitModelMboostLong(data2, m_max=300, nuisance=4), freq=FALSE)
  e2l <- getErrors(data2, m2l, plotModel=TRUE)
  e2ml <- getErrors(data2, m2ml, plotModel=TRUE)
  # cbind(e2l, e2ml)
  plotModel(err=e2l, errB=e2ml, data=data2, theseSettings=list(M=25, Gy=30, Gx=35, snrEps=2, scenario="2"))      
}


# save proportions of selected variables
selVar <- function(m, nuisance=0){
  
  allVariables <- c("ONEx, t", "ONEx, tlong",
                    #paste("z", 1, sep=""), 
                    paste("X", 1:(nuisance+2), sep=""))
  
  retlist <- vector("list",  length(allVariables))
  names(retlist) <- allVariables
  
  if(is.null(m)) return(retlist)
  
  xs <- selected(m)
  nm <- variable.names(m)
  selprob <- tabulate(xs, nbins = length(nm)) / length(xs)
  #names(selprob) <- names(nm) 
  
  sel <- data.frame(vars=nm, selprob)
  d1 <- data.frame(vars=allVariables)
  
  ret <- merge(d1, sel, by="vars", all.x=TRUE)
  
  # transform dataframe to a list
  retlist <- list()
  for(i in 1:nrow(ret)){
    retlist[[i]] <- ret[i,2]
  }
  names(retlist) <- ret[,1]
  
  return(retlist)
}


### one replication of simulation to plot the models of boosting and pffr
# theseSettings=settings[[2]]
# theseSettings=settingsSplit[[9]][[1]]
oneRep <- function(theseSettings){
  #browser()
  #print(data.frame(theseSettings))
    
  # Generate data
  arginS <- theseSettings$inS
  arginS <<- theseSettings$inS
  penaltyS <- theseSettings$penaltyS
  penaltyS <<- theseSettings$penaltyS
  diffPen <- theseSettings$diffPen
  diffPen <<- theseSettings$diffPen
  # theseSettings$type <- "locale"
  # theseSettings$type <- "locale0"
  # theseSettings$type <- "end"
  # theseSettings$type <- "end0"
  # theseSettings$type <- "start"
  # theseSettings$type <- "start0"
  # theseSettings$type <- "fourier"
  set.seed(theseSettings$seed)
  data <- do.call(makeData, theseSettings)
  # funplot(data$s, data$X1)
  args <- theseSettings
  args$data <- data
  
  # look for optimal mstop up to 1000 or 2000 depending on inS
  args$control <- boost_control(mstop=100, nu=0.1)
  if(arginS=="smooth"){
    args$grid <- seq(10, 1000, by = 10)
  }else{
    args$grid <- seq(10, 2000, by = 10)
  } 
  args$m_max <- 2000

  print(unlist(theseSettings)[c(1,3,4,6,9,11,19)])
    
  # Fit models

  if(args$p<1){
    #print("long format")
    modMboost <- suppressMessages(do.call(fitModelMboostLong, args)) # args$grid <- seq(10, 100, by = 10)
    modPffr <- suppressMessages(try(do.call(fitModelPffrLong, args)))
  }else{
    #print("wide format")
    modMboost <- suppressMessages(do.call(fitModelMboost, args)) # args$grid <- seq(10, 100, by = 10)
    modPffr <- suppressMessages(try(do.call(fitModelPffr, args)))
  }
  
  if(is.null(modPffr)) print(paste("modPffr, set " , theseSettings$set, ", is NULL", sep=""))
  if(is.null(modMboost)) print(paste("modMboost, set " , theseSettings$set, ", is NULL", sep=""))
  
  err <- getErrors(data, modPffr, plotModel=TRUE)
  errB <- getErrors(data, modMboost, plotModel=TRUE)
  
  # Save models of first rep
  if(theseSettings$rep==1){ 
    temp <- c("type", "a", "scenario", "k", "centerX", "penaltyS", "diffPen")
    temp2 <- theseSettings[temp]
    temp[temp %in% "scenario"] <- "sc"
    temp[temp %in% "type"] <- ""
    temp[temp %in% "penaltyS"] <- ""
    nm <- paste(temp, temp2, sep="", collapse="_")  
    
    pdf(paste(pathModels, nm, ".pdf", sep=""), width=7, height=12)
    try(plotModel(err=err, errB=errB, data=data, theseSettings=theseSettings))
    
    par(mfrow=c(4,2))
    funplot(data$s, data$X1, xlab="s", ylab="", main="X1")
    funplot(data$s, data$X2, xlab="s", ylab="", main="X2")
    
    predPffr2 <- predict(modPffr, type="terms")[[2]][[".value"]]
    predPffr3 <- predict(modPffr, type="terms")[[3]][[".value"]]
    
    predMboost2 <- predict(modMboost, which=2)#[,1]
    predMboost3 <- predict(modMboost, which=3)#[,1]
    
    if(length(predMboost2)==1) predMboost2 <- rep(0, length(data$tlong))
    if(length(predMboost3)==1) predMboost3 <- rep(0, length(data$tlong))
    
    range <- range(data$X1f, data$X2f, 
                   predMboost2, predMboost3, 
                   predPffr2, predPffr3)
    
    funplot(data$t, data$X1f, xlab="t", ylab="", main="effect of X1", ylim=range)
    funplot(data$t, data$X2f, xlab="t", ylab="", main="effect of X2", ylim=range)
    
    funplot(data$tlong, predMboost2, data$id, xlab="t", ylab="", main="FDboost: effect of X1", ylim=range)
    funplot(data$tlong, predMboost3, data$id, xlab="t", ylab="", main="FDboost: effect of X2", ylim=range)
    funplot(data$tlong, predPffr2, data$id, xlab="t", ylab="", main="pffr: effect of X1", ylim=range)
    funplot(data$tlong, predPffr3, data$id, xlab="t", ylab="", main="pffr: effect of X2", ylim=range)
    
    dev.off()   
  }
    
  # Calculate errors for both models 
  resMboost <- try(c(model=2, getErrors(data, modMboost)))
  resPffr <- try(c(model=1, getErrors(data, modPffr)))
  #cbind(resMboost, resPffr)
  
  rm(modMboost); rm(modPffr); rm(data); rm(args)
    
  res <- rbind(do.call(data.frame, c(theseSettings, resPffr)), 
        do.call(data.frame, c(theseSettings, resMboost)))
  
  return(res)
}


### one replication of simulation: only fit the models of mboost!!
# theseSettings=settings[[1]]
# theseSettings=settingsSplit[[9]][[1]]
oneRepFDboost <- function(theseSettings){
  #browser()
  #print(data.frame(theseSettings))
  
  # Generate data
  arginS <- theseSettings$inS
  arginS <<- theseSettings$inS
  penaltyS <- theseSettings$penaltyS
  penaltyS <<- theseSettings$penaltyS
  diffPen <- theseSettings$diffPen
  diffPen <<- theseSettings$diffPen
  
  set.seed(theseSettings$seed)
  data <- do.call(makeData, theseSettings)
  args <- theseSettings
  args$data <- data
  
  # look for optimal mstop up to 1000 or 2000 depending on inS
  args$control <- boost_control(mstop=100, nu=0.1)
  if(arginS=="smooth"){
    args$grid <- seq(10, 1000, by = 10)
  }else{
    args$grid <- seq(10, 2000, by = 10)
  } 
  args$m_max <- 2000
  
  # Fit models
  modPffr <- NULL
  ###  modPffr <- suppressMessages(do.call(fitModelPffr, args))
  
  if(FALSE){    
    formula <- attr(data, "formulaFDboostLong")
    timeformula <- attr(data, "timeformulaLong")
    
    if(FALSE){
      
      ## lead effect
      mylimits <- function(s, t) {
        (s < t) & (s < t-3)
      }
      
      ## band effect
      mylimits <- function(s, t) {
        (s < t + 5) & (s > t -5)
      }
      
      formula <- Ylong ~ 1 + bhist(X1, s, tlong, df = 5.0176, knots = 5, inS = arginS, 
                                   penalty = "ps", differences = 2, check.ident = TRUE) + 
        bhist(X2, s, tlong, df = 5.0176, knots = 5, inS = arginS, 
              penalty = penaltyS, differences = diffPen, check.ident = FALSE)
      
    }
    
    m <- FDboost(eval(formula), timeformula = timeformula, id=~id, 
                     data=data, control=boost_control(mstop=100, nu=0.1), 
                     numInt="Riemann", offset=NULL, 
                     offset_control=o_control(k_min=10)) 
    # test <- getDiags(m, data)
    
    str(extract(m))
    desMatrix <- extract(m)[[2]]
    dim(desMatrix)
    colSums(desMatrix) # there are columns that are completely zero!
    m$coef()[[2]]
    round(matrix(m$coef()[[2]], ncol=sqrt(ncol(desMatrix))), 2)
    plot(m, which=2, pers=TRUE)
    plot(m, which=2)
    persp(matrix(m$coef()[[2]], ncol=sqrt(ncol(desMatrix))), theta=30, phi=30)
    image(matrix(m$coef()[[2]], ncol=sqrt(ncol(desMatrix))))
  }
  
  if(args$p<1){
    #print("long format")
    modMboost <- suppressMessages(do.call(fitModelMboostLong, args)) # args$m_max=200
  }else{
    #print("wide format")
    modMboost <- suppressMessages(do.call(fitModelMboost, args)) # args$m_max=100
  }
  #browser()
   
  # Calculate errors for both models 
  resMboost <- try(c(model=2, getErrors(data, modMboost)))
  resPffr <- try(c(model=1, getErrors(data, modPffr)))
  # cbind(resMboost, resPffr)
  
  # compute the measures of identifiability
  identMboostX1 <- getDiags(modMboost, data, bl=2)
  identMboostX2 <- getDiags(modMboost, data, bl=3)
  names(identMboostX2) <- paste0(names(identMboostX2), ".2")
  
  #### selected variables
  selectedVariables <- selVar(m=modMboost, nuisance=8)
  
  
  rm(modMboost); rm(modPffr); rm(data); rm(args)
    
  res <- rbind(do.call(data.frame, c(theseSettings, resMboost, 
                                     unlist(identMboostX1), unlist(identMboostX2), 
                                     selectedVariables)))
  
  return(res)
}

#  test <- oneRep(theseSettings)


#s# one replication of simulation
# theseSettings=settings[[1]]
# theseSettings=settingsSplit[[9]][[1]]
oneRepPffr <- function(theseSettings){
  
  # Generate data
  arginS <- theseSettings$inS
  arginS <<- theseSettings$inS
  penaltyS <- theseSettings$penaltyS
  penaltyS <<- theseSettings$penaltyS
  diffPen <- theseSettings$diffPen
  diffPen <<- theseSettings$diffPen
  
  set.seed(theseSettings$seed)
  data <- do.call(makeData, theseSettings)
  #print(data$Y[5,5])
  #print(data$Ytrue[5,5])
  #print(data$X1[5,5])
  #print(data$X2[5,5])
  args <- theseSettings
  args$data <- data
  
  print(unlist(theseSettings)[c(1,3,4,6,9,11,19)])
    
  # Fit models
  #modPffr <- NULL
  
  if(args$p<1){
    #print("long format")
    modPffr <- suppressMessages(try(do.call(fitModelPffrLong, args)))
  }else{
    #print("wide format")
    modPffr <- suppressMessages(try(do.call(fitModelPffr, args)))
  }
  
  # Should not be necessary!
  if(any(class(modPffr)=="try-error")){
    warning(paste("fitModelPffr, set " , theseSettings$set, ", gives error!", sep=""))
    modPffr <- NULL
  }
    
  ###  modMboost <- suppressMessages(do.call(fitModelMboost, args))
  modMboost <- NULL
  
  if(is.null(modPffr)){
    print(paste("modPffr, set " , theseSettings$set, ", is NULL", sep=""))
  } 
    
  # Calculate errors for both models 
  resMboost <- try(c(model=2, getErrors(data, modMboost)))
  resPffr <- try(c(model=1, getErrors(data, modPffr)))
  #cbind(resMboost, resPffr)
  
  # compute the measures of identifiability
  identPffrX1 <- getDiags(modPffr, data, bl=2)
  identPffrX2 <- getDiags(modPffr, data, bl=3)
  names(identPffrX2) <- paste0(names(identPffrX2), ".2")
  
  rm(modMboost); rm(modPffr); rm(data); rm(args)

  res <- rbind(do.call(data.frame, c(theseSettings, resPffr, 
                                     unlist(identPffrX1), unlist(identPffrX2))))
    
  return(res)
}


