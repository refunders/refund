
#' Cross-Validation and Bootstrapping over Curves
#' 
#' Cross-validation and bootstrapping over curves to compute the empirical risk for 
#' hyper-parameter selection.    
#' 
#' @param object fitted FDboost-object
#' @param folds a weight matrix with number of rows equal to the number of observed trajectories.  
#' @param grid the grid over which the optimal number of boosting iterations (mstop) is searched.  
#' @param showProgress logical, defaults to \code{TRUE}.
#' @param compress logical, defaults to \code{FALSE}. Only used to force a meaningful
#' behaviour of \code{applyFolds} with hmatrix objects when using nested resampling.
#' @param papply (parallel) apply function, defaults to \code{\link[parallel]{mclapply}}, 
#' see \code{\link[mboost]{cvrisk}} for details.  
#' @param fun if \code{fun} is \code{NULL}, the out-of-bag risk is returned. 
#' \code{fun}, as a function of \code{object}, 
#' may extract any other characteristic of the cross-validated models. These are returned as is.
#' @param riskFun only exists in \code{applyFolds}; allows to compute other risk functions than the risk 
#' of the family that was specified in object. 
#' Must be specified as function of arguments \code{(y, f, w = 1)}, where \code{y} is the 
#' observed response, \code{f} is the prediciton from the model and \code{w} is the weight. 
#' The risk function must return a scalar numeric value for vector valued imput.  
#' @param numInt only exists in \code{applyFolds}; the scheme for numerical integration, 
#' see \code{numInt} in \code{\link{FDboost}}. 
#' @param corrected see \code{\link[mboost]{cvrisk}}. 
#' @param mc.preschedule Defaults to \code{FALSE}. Preschedule tasks if are parallelized using \code{mclapply}?  
#' For details see \code{\link[parallel]{mclapply}}. 
#' @param ... further arguments passed to \code{\link[parallel]{mclapply}} 
#' 
#' @param id the id-vector as integers 1, 2, ... specifying which observations belong to the same curve, 
#' deprecated in \code{cvMa()}. 
#' @param weights a numeric vector of (integration) weights, defaults to 1.
#' @param type character argument for specifying the cross-validation 
#' method. Currently (stratified) bootstrap, k-fold cross-validation, subsampling and 
#' leaving-one-curve-out cross validation (i.e. jack knife on curves) are implemented.
#' @param B number of folds, per default 25 for \code{bootstrap} and
#' \code{subsampling} and 10 for \code{kfold}.
#' @param prob percentage of observations to be included in the learning samples 
#' for subsampling.
#' @param strata a factor of the same length as \code{weights} for stratification.
#' 
#' @param ydim dimensions of response-matrix
#' 
#' @details The number of boosting iterations is an important hyper-parameter of boosting.   
#' It be chosen using the functions \code{applyFolds} or \code{cvrisk.FDboost}. Those functions 
#' they compute honest, i.e., out-of-bag, estimates of the empirical risk for different 
#' numbers of boosting iterations. 
#' The weights (zero weights correspond to test cases) are defined via the folds matrix, 
#' see \code{\link[mboost]{cvrisk}} in package mboost. 
#' 
#' In case of functional response, we recommend to use \code{applyFolds}. 
#' It recomputes the model in each fold using \code{FDboost}. Thus, all parameters are recomputed, 
#' including the smooth offset (if present) and the identifiability constraints (if present, only 
#' relevant for \code{bolsc}, \code{brandomc} and \code{bbsc}).  
#' Note, that the function \code{applyFolds} expects folds that give weights
#' per curve without considering integration weights.  
#' 
#' The function \code{cvrisk.FDboost} is a wrapper for \code{\link[mboost]{cvrisk}} in package mboost. 
#' It overrides the default for the folds, so that the folds are sampled on the level of curves 
#' (not on the level of single observations, which does not make sense for functional response).  
#' Note that the smooth offset and the computation of the identifiability constraints
#' are not part of the refitting if \code{cvrisk} is used. 
#' Per default the integration weights of the model fit are used to compute the prediction errors 
#' (as the integration weights are part of the default folds). 
#' Note that in \code{cvrisk} the weights are rescaled to sum up to one. 
#' 
#' The functions \code{cvMa} and \code{cvLong} can be used to build an appropriate 
#' weight matrix for functional response to be used with \code{cvrisk} as sampling 
#' is done on the level of curves. The probability for each 
#' curve to enter a fold is equal over all curves.     
#' The function \code{cvMa} takes the dimensions of the response matrix as input argument and thus
#' can only be used for regularly observed response. 
#' The function \code{cvLong} takes the id variable and the weights as arguments and thus can be used
#' for responses in long format that are potentially observed irregularly. 
#'  
#' If \code{strata} is defined 
#' sampling is performed in each stratum separately thus preserving 
#' the distribution of the \code{strata} variable in each fold. 
#' 
#' @note Use argument \code{mc.cores = 1L} to set the numbers of cores that is used in 
#' parallel computation. On Windows only 1 core is possible, \code{mc.cores = 1}, which is the default.
#' 
#' @seealso \code{\link[mboost]{cvrisk}} to perform cross-validation with scalar response.
#' 
#' @return \code{cvMa} and \code{cvLong} return a matrix of sampling weights to be used in \code{cvrisk}. 
#' 
#' The functions \code{applyFolds} and \code{cvrisk.FDboost} return a \code{cvrisk}-object, 
#' which is a matrix of the computed out-of-bag risk. The matrix has the folds in rows and the 
#' number of boosting iteratins in columns. Furhtermore, the matrix has attributes including: 
#' \item{risk}{name of the applied risk function}
#' \item{call}{model call of the model object}
#' \item{mstop}{gird of stopping iterations that is used}
#' \item{type}{name for the type of folds}
#' 
#' @examples
#' Ytest <- matrix(rnorm(15), ncol = 3) # 5 trajectories, each with 3 observations 
#' Ylong <- as.vector(Ytest)
#' ## 4-folds for bootstrap for the response in long format without integration weights
#' cvMa(ydim = c(5,3), type = "bootstrap", B = 4)  
#' cvLong(id = rep(1:5, times = 3), type = "bootstrap", B = 4)
#' 
#' if(require(fda)){
#'  ## load the data
#'  data("CanadianWeather", package = "fda")
#'  
#'  ## use data on a daily basis 
#'  canada <- with(CanadianWeather, 
#'                 list(temp = t(dailyAv[ , , "Temperature.C"]),
#'                      l10precip = t(dailyAv[ , , "log10precip"]),
#'                      l10precip_mean = log(colMeans(dailyAv[ , , "Precipitation.mm"]), base = 10),
#'                      lat = coordinates[ , "N.latitude"],
#'                      lon = coordinates[ , "W.longitude"],
#'                      region = factor(region),
#'                      place = factor(place),
#'                      day = 1:365,  ## corresponds to t: evaluation points of the fun. response 
#'                      day_s = 1:365))  ## corresponds to s: evaluation points of the fun. covariate
#'  
#' ## center temperature curves per day 
#' canada$tempRaw <- canada$temp
#' canada$temp <- scale(canada$temp, scale = FALSE) 
#' rownames(canada$temp) <- NULL ## delete row-names 
#'   
#' ## fit the model  
#' mod <- FDboost(l10precip ~ 1 + bolsc(region, df = 4) + 
#'                  bsignal(temp, s = day_s, cyclic = TRUE, boundary.knots = c(0.5, 365.5)), 
#'                timeformula = ~ bbs(day, cyclic = TRUE, boundary.knots = c(0.5, 365.5)), 
#'                data = canada)
#' mod <- mod[75]
#' 
#' \dontrun{
#'   #### create folds for 3-fold bootstrap: one weight for each curve
#'   set.seed(123)
#'   folds_bs <- cv(weights = rep(1, mod$ydim[1]), type = "bootstrap", B = 3)
#' 
#'   ## compute out-of-bag risk on the 3 folds for 1 to 75 boosting iterations  
#'   cvr <- applyFolds(mod, folds = folds_bs, grid = 1:75)
#' 
#'   ## weights per observation point  
#'   folds_bs_long <- folds_bs[rep(1:nrow(folds_bs), times = mod$ydim[2]), ]
#'   attr(folds_bs_long, "type") <- "3-fold bootstrap"
#'   ## compute out-of-bag risk on the 3 folds for 1 to 75 boosting iterations  
#'   cvr3 <- cvrisk(mod, folds = folds_bs_long, grid = 1:75)
#' }
#' 
#' \dontrun{
#'   ## plot the out-of-bag risk
#'   par(mfrow = c(1,3))
#'   plot(cvr); legend("topright", lty=2, paste(mstop(cvr)))
#'   plot(cvr3); legend("topright", lty=2, paste(mstop(cvr3)))
#' }
#' 
#'}
#' 
#' @aliases cvMa cvLong cvrisk.FDboost
#' 
#' @export
## computes the empirical out-of-bag risk for each fold  
applyFolds <- function(object, folds = cv(rep(1, length(unique(object$id))), type = "bootstrap"),
                       grid = 1:mstop(object), fun = NULL, 
                       riskFun = NULL, numInt = object$numInt, 
                       papply = mclapply, 
                       mc.preschedule = FALSE, showProgress = TRUE, 
                       compress = FALSE,
                       ...) {
  
  if (is.null(folds)) {
    stop("Specify folds.")
  } 
  ## check that folds are given on the level of curves
  if(length(unique(object$id)) != nrow(folds)){
    stop("The folds-matrix must have one row per observed trajectory.")
  }
  
  if(any(class(object) == "FDboostLong")){ # irregular response
    nObs <- length(unique(object$id)) # number of curves
    Gy <- NULL # number of time-points per curve
  }else{ # regular response / scalar response 
    nObs <- object$ydim[1] # number of curves
    Gy <- object$ydim[2] # number of time-points per curve
    if(class(object)[1] == "FDboostScalar"){
      nObs <- length(object$response)
      Gy <- 1
    }
  }
  
  sample_weights <- rep(1, length(unique(object$id))) # length N
  # if(any(sample_weights == 0)) warning("zero weights") # fullfilled per construction
  
  # save integration weights of original model
  if(is.null(numInt)){
    numInt <- "equal"
    warning("'numInt' is NULL. It is set to 'equal' which means that all integration weights are set to 1.")
  }
  if(!numInt %in% c("equal", "Riemann")) 
    warning("argument 'numInt' is ignored as it is none of 'equal' and 'Riemann'.")
  
  if(is.numeric(numInt)){  # use the integration scheme specified in applyFolds
    if(length(numInt) != length(object$yind)) stop("Length of integration weights and time vector are not equal.")
    integration_weights <- numInt
  }else{
    
    if(numInt == "Riemann"){ # use the integration scheme specified in applyFolds
      if(!any(class(object) == "FDboostLong")){
        integration_weights <- as.vector(integrationWeights(X1 = matrix(object$response, 
                                                                        ncol = object$ydim[2]), object$yind))
      }else{
        integration_weights <- integrationWeights(X1 = object$response, object$yind, object$id)
      }
    }else{ ## numInt == "equal"
      integration_weights <- rep(1, length(object$response))
      ## correct integration weights for matrix valued response like possibly in Binomial()
      if( class(object)[1] == "FDboostScalar") integration_weights <- rep(1, NROW(object$response))
    }
  }
  
  ### get yind in long format
  yindLong <- object$yind
  if(!any(class(object) == "FDboostLong")){
    yindLong <- rep(object$yind, each = nObs)
  }
  ### compute ("length of each trajectory")^-1 in the response
  ### more precisely ("sum of integration weights")^-1 is used
  if(numInt == "equal"){
    lengthTi1 <- rep(1, l = length(unique(object$id)))
  }else{
    if(length(object$yind) > 1){
      # lengthTi1 <- 1/tapply(yindLong[!is.na(response)], object$id[!is.na(response)], function(x) max(x) - min(x))
      lengthTi1 <- 1/tapply(integration_weights, object$id, function(x) sum(x))
      if(any(is.infinite(lengthTi1))) lengthTi1[is.infinite(lengthTi1)] <- max(lengthTi1[!is.infinite(lengthTi1)])
    }else{
      lengthTi1 <- rep(1, l = length(object$response))
    }
  }
  
  # Function to suppress the warning of missings in the response
  h <- function(w){
    if( any( grepl( "response contains missing values;", w) ) ) 
      invokeRestart( "muffleWarning" )  
  }
  
  ### start preparing the data 
  dathelp <- object$data
  nameyind <- attr(object$yind, "nameyind")
  dathelp[[nameyind]] <- object$yind
  
  ## try to set up data using $get_data() 
  ## problem with index for bl containing index, and you do not get s for bsignal/bhist
  if(FALSE){
    dathelp2 <- list()
    for(j in 1:length(object$baselearner)){
      dat_bl_j <- object$baselearner[[j]]$get_data() ## object$baselearner[[j]]$model.frame()
      # if the variable is already present, do not add it again
      dathelp2 <- c(dathelp2, dat_bl_j[!names(dat_bl_j) %in% names(dathelp2)])
    }
  }
  
  if(!any(class(object) == "FDboostLong") & !any(class(object) == "FDboostScalar")){
    dathelp[[object$yname]] <- matrix(object$response, ncol=object$ydim[2])
    dathelp$integration_weights <- matrix(integration_weights, ncol=object$ydim[2])
    dathelp$object_id <- object$id
  }else{
    dathelp[[object$yname]] <- object$response
    dathelp$integration_weights <- integration_weights
  }
  
  ## get the names of all variables x_i, i = 1, ... , N
  names_variables <- unlist(lapply(object$baselearner, function(x) x$get_names() ))
  ## check for index
  has_index <- sapply(object$baselearner, function(x) !is.null(x$get_index()))
  if( any( has_index )){
    index_names <- sapply(lapply(object$baselearner[has_index], function(x) x$get_call()), 
                                 function(l) gsub("index[[:space:]*]=[[:space:]*]|\\,","",
                                                  regmatches(l, regexpr('index[[:space:]*]=.*\\,', l)))
    )
  } else index_names <- NULL
  
  names(names_variables) <- NULL
  names_variables <- names_variables[names_variables != nameyind]
  names_variables <- names_variables[names_variables != "ONEx"]
  names_variables <- names_variables[names_variables != "ONEtime"]
  if(!any(class(object) == "FDboostLong")) names_variables <- c(object$yname, "integration_weights", names_variables)
  
  length_variables <- if("FDboostScalar" %in% class(object)) 
    lapply(dathelp[names_variables], length) else
      lapply(dathelp[names_variables], NROW)
  names_variables_long <- names_variables[ length_variables == length(object$id) ]
  nothmatrix <- ! sapply(dathelp[names_variables_long], function(x) any(class(x) == "hmatrix" ))
  names_variables_long <- names_variables_long[ nothmatrix ]
  if(identical(names_variables_long, character(0))) names_variables_long <- NULL
  
  names_variables <- names_variables[! names_variables %in% names_variables_long ]
  if(identical(names_variables, character(0))) names_variables <- NULL
  
  ## check if there is a baselearner without brackets
  # the probelm with such base-learners is that their data is not contained in object$data
  # using object$baselearner[[j]]$get_data() is difficult as this can be blow up by index for %X%
  singleBls <- gsub("\\s", "", unlist(lapply(strsplit(
    strsplit(object$formulaFDboost, "~")[[1]][2], # split formula
                                      "\\+")[[1]], # split additive terms
    function(y) strsplit(y, split = "%.{1,3}%")) # split single baselearners
  )) 
  
  singleBls <- singleBls[singleBls != "1"]
  
  if(any(!grepl("\\(", singleBls))) 
    stop(paste0("applyFolds can not deal with the following base-learner(s) without brackets: ", 
                paste(singleBls[!grepl("\\(", singleBls)], collapse = ", ")))
  
  
  ## check if data includes all variables
  if(any(whMiss <- ! c(names_variables,
             object$yname, 
             nameyind, 
             "integration_weights", 
             names_variables_long) %in% names(dathelp))){
    
    # for each missing variable get the first baselearner, which contains the variable
    blWithMissVars <- lapply(names_variables[whMiss], function(w) 
      unlist(lapply(1:length(object$baselearner), function(i) if(
        any( grepl(w, object$baselearner[[i]]$get_names() ) )) return(i))
      )[1])
    
    stop(paste0("base-learner(s) ", paste(unlist(list(1,2)), collapse = ", "), 
                " contain(s) variables, which are not part of the data object."))
    
  }
  
  ## fitfct <- object$update
  fitfct <- function(weights, oobweights){
    
    ## get data according to weights
    if(any(class(object) == "FDboostLong")){
      dat_weights <- reweightData(data = dathelp, vars = names_variables, 
                                  longvars = c(object$yname, nameyind, "integration_weights", names_variables_long),  
                                  weights = weights, idvars = c(attr(object$id, "nameid"), index_names),
                                  compress = compress)
    }else if(class(object)[1] == "FDboostScalar"){
      dat_weights <- reweightData(data = dathelp, 
                                  vars = c(names_variables, names_variables_long),
                                  weights = weights)
    }else{
      dat_weights <- reweightData(data = dathelp, vars = names_variables, 
                                  longvars = names_variables_long, 
                                  weights = weights, idvars = c("object_id", index_names))
    }

    # check for factors
    
    isFac <- sapply(dathelp, is.factor)
    
    if(any(isFac)){
      namesFac <- names(isFac)[isFac]
      
      for(i in 1:length(namesFac)){
        
      if(length(levels(droplevels(dathelp[[namesFac[i]]]))) != 
         length(levels(droplevels(dat_weights[[namesFac[i]]]))))
        stop(paste0("The factor variable '", namesFac[i], "' has unobserved levels in the training data. ",
                    "Make sure that training data in each fold contains all factor levels."))
        
      }
    }
    
    call <- object$callEval
    call$data <- dat_weights
    if(! is.null(call$weights)) warning("Argument weights of original model is not considered.")
    call$weights <- NULL 
    
    ## fit the model for dat_weights  
    mod <- withCallingHandlers(suppressMessages(eval(call)), warning = h) # suppress the warning of missing responses    
    mod <- mod[max(grid)]
    
    mod
    
  }
  
  ## create data frame for model fit and use the weights vector for the CV
  # oobrisk <- matrix(0, nrow = ncol(folds), ncol = length(grid))
  
  if (!is.null(fun))
    stopifnot(is.function(fun))
  
  fam_name <- object$family@name
  
  call <- deparse(object$call)
  
  
  if (is.null(fun)) {
    dummyfct <- function(weights, oobweights) {

      mod <- fitfct(weights = weights, oobweights = oobweights)
      mod <- mod[max(grid)]
      # mod$risk()[grid]
      
      # get risk function of the family
      if(is.null(riskFun)){
        myfamily <- get("family", environment(mod$update))
        riskfct <- myfamily@risk
      }else{
        stopifnot(is.function(riskFun))
        riskfct <- riskFun 
      }

      ## get data according to oobweights
      if(any(class(object) == "FDboostLong")){
        dathelp$lengthTi1 <- c(lengthTi1)
        dat_oobweights <- reweightData(data = dathelp, vars = c(names_variables, "lengthTi1"),  
                                       longvars = c(object$yname, nameyind, 
                                                    "integration_weights", names_variables_long),  
                                       weights = oobweights, 
                                       idvars = c(attr(object$id, "nameid"), index_names),
                                       compress = compress)
        ## funplot(dat_oobweights[[nameyind]], dat_oobweights[[object$yname]], 
        ##         id = dat_oobweights[[attr(object$id, "nameid")]])
        
        for(v in names_variables){  ## blow up covariates by id so that data can be used with predict()
          if(!is.null(dim(dat_oobweights[[v]]))){
            dat_oobweights[[v]] <- dat_oobweights[[v]][dat_oobweights[[attr(object$id, "nameid")]], ]
          }else{
            dat_oobweights[[v]] <- dat_oobweights[[v]][dat_oobweights[[attr(object$id, "nameid")]]]
          }
        }
        response_oobweights <- c(dat_oobweights[[object$yname]])
        
      }else{ ## scalar or regular response
        
        if(class(object)[1] == "FDboostScalar"){
          dat_oobweights <- reweightData(data = dathelp, 
                                         vars = c(names_variables, names_variables_long),
                                         weights = oobweights)
          response_oobweights <- dat_oobweights[[object$yname]]
        }else{
          dat_oobweights <- reweightData(data = dathelp, vars = names_variables, 
                                         longvars = names_variables_long,
                                         weights = oobweights, idvars = c("object_id", index_names))
          response_oobweights <- c(dat_oobweights[[object$yname]])
        }
        
      }
      
      if(is.character(response_oobweights)) response_oobweights <- factor(response_oobweights)
      ## this check is important for Binomial() as it recodes factor to -1, 1
      response_oobweights <- myfamily@check_y(response_oobweights)
      
      # Function to suppress the warning of extrapolation in bbs / bbsc 
      h2 <- function(w){
        if( any( grepl( "Linear extrapolation used.", w) ) ) 
          invokeRestart( "muffleWarning" )  
      }
      
      if(any(class(object) == "FDboostLong")){
        
        if(numInt == "equal"){ 
          oobwstand <- dat_oobweights$integration_weights * (1/sum(dat_oobweights$integration_weights))
        }else{
          # compute integration weights for standardizing risk 
          oobwstand <- dat_oobweights$lengthTi1[dat_oobweights[[attr(object$id, "nameid")]]] * 
            dat_oobweights$integration_weights * (1/sum(oobweights))
        }
        
        # compute risk with integration weights like in FDboost::validateFDboost
        risk <- sapply(grid, function(g){riskfct( response_oobweights, 
                                                  withCallingHandlers(predict(mod[g], newdata = dat_oobweights, toFDboost = FALSE), warning = h2), 
                                                  w = oobwstand )}) ## oobwstand[oobweights[object$id] != 0 ]
        
      }else{
        
        if(numInt == "equal"){ # oobweights for i = 1, ..., N
          oobwstand <- oobweights[object$id]*(1/sum(oobweights[object$id]))
        }else{
          # compute integration weights for standardizing risk 
          oobwstand <- lengthTi1[object$id]*oobweights[object$id]*integration_weights*(1/sum(oobweights))
        }
        
        # compute risk with integration weights like in FDboost::validateFDboost
        risk <- sapply(grid, function(g){riskfct( response_oobweights, 
                                                  withCallingHandlers(predict(mod[g], newdata = dat_oobweights, toFDboost = FALSE), warning = h2), 
                                                  w = oobwstand[oobweights != 0 ])})
        
      }

      if(showProgress) cat(".")
      
      risk
    }
    
  } else { ## !is.null(fun)
    
    if(!is.null(riskFun)) warning("riskFun is ignored as fun is specified.")
    
    dummyfct <- function(weights, oobweights) {
      mod <- fitfct(weights = weights, oobweights = oobweights)
      mod[max(grid)]
      ## make sure dispatch works correctly
      class(mod) <- class(object)
      
      fun(mod) # Provide an extra argument for dat_oobweights?
    }
  }
  
  ## use case weights as out-of-bag weights (but set inbag to 0)
  OOBweights <- matrix(rep(sample_weights, ncol(folds)), ncol = ncol(folds))
  OOBweights[folds > 0] <- 0
  
  if (all.equal(papply, mclapply) == TRUE) {
    oobrisk <- papply(1:ncol(folds),
                      function(i) try(dummyfct(weights = folds[, i],
                                               oobweights = OOBweights[, i]),
                                      silent = TRUE),
                      mc.preschedule = mc.preschedule,
                      ...)
  } else {
    oobrisk <- papply(1:ncol(folds),
                      function(i) try(dummyfct(weights = folds[, i],
                                               oobweights = OOBweights[, i]),
                                      silent = TRUE),
                      ...)
  }
  
  ## if any errors occured remove results and issue a warning
  if (any(idx <- sapply(oobrisk, is.character))) {
    if(sum(idx) == length(idx)){
      stop("All folds encountered an error.\n",
           "Original error message(s):\n",
           sapply(oobrisk[idx], function(x) x))
    } 
    warning(sum(idx), " fold(s) encountered an error. ",
            "Results are based on ", ncol(folds) - sum(idx),
            " folds only.\n",
            "Original error message(s):\n",
            sapply(oobrisk[idx], function(x) x))
    oobrisk[idx] <- NULL
  }
  
  if (!is.null(fun))
    return(oobrisk)
  
  oobrisk <- t(as.data.frame(oobrisk))
  ## oobrisk <- oobrisk  / colSums(OOBweights[object$id, ]) # is done in dummyfct() 
  colnames(oobrisk) <- grid
  rownames(oobrisk) <- 1:nrow(oobrisk)
  attr(oobrisk, "risk") <- fam_name
  attr(oobrisk, "call") <- call
  attr(oobrisk, "mstop") <- grid
  attr(oobrisk, "type") <- ifelse(!is.null(attr(folds, "type")),
                                  attr(folds, "type"), "user-defined")
  
  class(oobrisk) <- c("cvrisk", "applyFolds")
  
  oobrisk
}


#' Cross-Validation and Bootstrapping over Curves
#' 
#' DEPRECATED! 
#' The function \code{validateFDboost()} is deprecated,  
#' use \code{\link{applyFolds}} and \code{\link{bootstrapCI}} instead. 
#' 
#' @param object fitted FDboost-object
#' @param response optional, specify a response vector for the computation of the prediction errors.  
#' Defaults to \code{NULL} which means that the response of the fitted model is used.
#' @param folds a weight matrix with number of rows equal to the number of observed trajectories.  
#' @param grid the grid over which the optimal number of boosting iterations (mstop) is searched.  
#' @param getCoefCV logical, defaults to \code{TRUE}. Should the coefficients and predictions
#' be computed for all the models on the sampled data?
#' @param riskopt how is the optimal stopping iteration determined. Defaults to the mean, 
#' but median is possible as well. 
#' @param mrdDelete Delete values that are \code{mrdDelete} percent smaller than the mean
#' of the response. Defaults to 0 which means that only response values being 0 
#' are not used in the calculation of the MRD (= mean relative deviation).  
#' @param refitSmoothOffset logical, should the offset be refitted in each learning sample? 
#' Defaults to \code{TRUE}. In \code{\link[mboost]{cvrisk}} the offset of the original model fit in  
#' \code{object} is used in all folds.
#' @param showProgress logical, defaults to \code{TRUE}.
#' @param fun if \code{fun} is \code{NULL}, the out-of-bag risk is returned. 
#' \code{fun}, as a function of \code{object}, 
#' may extract any other characteristic of the cross-validated models. These are returned as is.
#' 
#' @param ... further arguments passed to \code{\link[parallel]{mclapply}} 
#' 
#' @details The number of boosting iterations is an important hyper-parameter of boosting  
#' and can be chosen using the function \code{validateFDboost} as they compute
#' honest, i.e., out-of-bag, estimates of the empirical risk for different numbers of boosting iterations. 
#' 
#' The function \code{validateFDboost} is especially suited to models with functional response. 
#' Using the option \code{refitSmoothOffset} the offset is refitted on each fold. 
#' Note, that the function \code{validateFDboost} expects folds that give weights
#' per curve without considering integration weights. The integration weights of 
#' \code{object} are used to compute the empirical risk as integral. The argument \code{response} 
#' can be useful in simulation studies where the true value of the response is known but for 
#' the model fit the response is used with noise. 
#' 
#' @return The function \code{validateFDboost} returns a \code{validateFDboost}-object, 
#' which is a named list containing: 
#' \item{response}{the response}
#' \item{yind}{the observation points of the response}
#' \item{id}{the id variable of the response}
#' \item{folds}{folds that were used}
#' \item{grid}{grid of possible numbers of boosting iterations}
#' \item{coefCV}{if \code{getCoefCV} is \code{TRUE} the estimated coefficient functions in the folds}
#' \item{predCV}{if \code{getCoefCV} is \code{TRUE} the out-of-bag predicted values of the response}
#' \item{oobpreds}{if the type of folds is curves the out-of-bag predictions for each trajectory}
#' \item{oobrisk}{the out-of-bag risk}
#' \item{oobriskMean}{the out-of-bag risk at the minimal mean risk}
#' \item{oobmse}{the out-of-bag mean squared error (MSE)}
#' \item{oobrelMSE}{the out-of-bag relative mean squared error (relMSE)}
#' \item{oobmrd}{the out-of-bag mean relative deviation (MRD)}
#' \item{oobrisk0}{the out-of-bag risk without consideration of integration weights}
#' \item{oobmse0}{the out-of-bag mean squared error (MSE) without consideration of integration weights}
#' \item{oobmrd0}{the out-of-bag mean relative deviation (MRD) without consideration of integration weights}
#' \item{format}{one of "FDboostLong" or "FDboost" depending on the class of the object}
#' \item{fun_ret}{list of what fun returns if fun was specified}
#' 
#' @examples
#' \dontrun{
#' if(require(fda)){
#'  ## load the data
#'  data("CanadianWeather", package = "fda")
#'  
#'  ## use data on a daily basis 
#'  canada <- with(CanadianWeather, 
#'                 list(temp = t(dailyAv[ , , "Temperature.C"]),
#'                      l10precip = t(dailyAv[ , , "log10precip"]),
#'                      l10precip_mean = log(colMeans(dailyAv[ , , "Precipitation.mm"]), base = 10),
#'                      lat = coordinates[ , "N.latitude"],
#'                      lon = coordinates[ , "W.longitude"],
#'                      region = factor(region),
#'                      place = factor(place),
#'                      day = 1:365,  ## corresponds to t: evaluation points of the fun. response 
#'                      day_s = 1:365))  ## corresponds to s: evaluation points of the fun. covariate
#'  
#' ## center temperature curves per day 
#' canada$tempRaw <- canada$temp
#' canada$temp <- scale(canada$temp, scale = FALSE) 
#' rownames(canada$temp) <- NULL ## delete row-names 
#'   
#' ## fit the model  
#' mod <- FDboost(l10precip ~ 1 + bolsc(region, df = 4) + 
#'                  bsignal(temp, s = day_s, cyclic = TRUE, boundary.knots = c(0.5, 365.5)), 
#'                timeformula = ~ bbs(day, cyclic = TRUE, boundary.knots = c(0.5, 365.5)), 
#'                data = canada)
#' mod <- mod[75]
#' 
#'   #### create folds for 3-fold bootstrap: one weight for each curve
#'   set.seed(123)
#'   folds_bs <- cv(weights = rep(1, mod$ydim[1]), type = "bootstrap", B = 3)
#' 
#'   ## compute out-of-bag risk on the 3 folds for 1 to 75 boosting iterations  
#'   cvr <- applyFolds(mod, folds = folds_bs, grid = 1:75)
#' 
#'   ## compute out-of-bag risk and coefficient estimates on folds  
#'   cvr2 <- validateFDboost(mod, folds = folds_bs, grid = 1:75)
#' 
#'   ## weights per observation point  
#'   folds_bs_long <- folds_bs[rep(1:nrow(folds_bs), times = mod$ydim[2]), ]
#'   attr(folds_bs_long, "type") <- "3-fold bootstrap"
#'   ## compute out-of-bag risk on the 3 folds for 1 to 75 boosting iterations  
#'   cvr3 <- cvrisk(mod, folds = folds_bs_long, grid = 1:75)
#' 
#'   ## plot the out-of-bag risk
#'   par(mfrow = c(1,3))
#'   plot(cvr); legend("topright", lty=2, paste(mstop(cvr)))
#'   plot(cvr2)
#'   plot(cvr3); legend("topright", lty=2, paste(mstop(cvr3)))
#' 
#'   ## plot the estimated coefficients per fold
#'   ## more meaningful for higher number of folds, e.g., B = 100 
#'   par(mfrow = c(2,2))
#'   plotPredCoef(cvr2, terms = FALSE, which = 2)
#'   plotPredCoef(cvr2, terms = FALSE, which = 3)
#'   
#'   ## compute out-of-bag risk and predictions for leaving-one-curve-out cross-validation
#'   cvr_jackknife <- validateFDboost(mod, folds = cvLong(unique(mod$id), 
#'                                    type = "curves"), grid = 1:75)
#'   plot(cvr_jackknife)
#'   ## plot oob predictions per fold for 3rd effect 
#'   plotPredCoef(cvr_jackknife, which = 3) 
#'   ## plot coefficients per fold for 2nd effect
#'   plotPredCoef(cvr_jackknife, which = 2, terms = FALSE)
#' 
#'}
#'}
#' 
#' @export
validateFDboost <- function(object, response = NULL,  
                            #folds=cvMa(ydim=object$ydim, weights=model.weights(object), type="bootstrap"),
                            folds = cv(rep(1, length(unique(object$id))), type = "bootstrap"),
                            grid = 1:mstop(object), 
                            fun = NULL, 
                            getCoefCV = TRUE, riskopt = c("mean","median"), 
                            mrdDelete = 0, refitSmoothOffset = TRUE, 
                            showProgress = TRUE, ...){
  
  .Deprecated(new = "applyFolds", 
              msg = "'validateFDboost' is deprecated. Use 'applyFolds' and 'bootstrapCI' instead.")
  
  names_bl <- names(object$baselearner)
  if(any(grepl("brandomc", names_bl))) message("For brandomc, the transformation matrix Z is fixed over all folds.")
  if(any(grepl("bolsc", names_bl))) message("For bolsc, the transformation matrix Z is fixed over all folds.")
  if(any(grepl("bbsc", names_bl))) message("For bbsc, the transformation matrix Z is fixed over all folds.")
  
  type <- attr(folds, "type")
  if(is.null(type)) type <- "unknown"
  call <- match.call()
  riskopt <- match.arg(riskopt)
  
  ## check that folds are given on the level of curves
  if(length(unique(object$id)) != nrow(folds)){
    stop("The folds-matrix must have one row per observed trajectory.")
  }
  
  if(any(class(object) == "FDboostLong")){ # irregular response
    nObs <- length(unique(object$id)) # number of curves
    Gy <- NULL # number of time-points per curve
  }else{ # regular response / scalar response 
    nObs <- object$ydim[1] # number of curves
    Gy <- object$ydim[2] # number of time-points per curve
    if(class(object)[1] == "FDboostScalar"){
      nObs <- length(object$response)
      Gy <- 1
    }
  }
  
  myfamily <- get("family", environment(object$update))
  
  if(is.null(response)) response <- object$response # response as vector!
  # for Binomial() transform factor to -1/1 coding
  response <- myfamily@check_y(response)
  
  id <- object$id
  
  # save integration weights of original model
  # intWeights <- model.weights(object) 
  # weights are rescaled in mboost, see mboost:::rescale_weights  
  if(!is.null(object$callEval$numInt) && object$callEval$numInt == "Riemann"){
    if(!any(class(object) == "FDboostLong")){
      intWeights <- as.vector(integrationWeights(X1 = matrix(object$response, 
                              ncol = object$ydim[2]), object$yind))
    }else{
      intWeights <- integrationWeights(X1=object$response, object$yind, id)
    }
  }else{
    intWeights <- model.weights(object) 
  }
  
  # out-of-bag-weights: i.e. the left out curve/ the left out observations
  OOBweights <- matrix(1, ncol = ncol(folds), nrow = nrow(folds))
  OOBweights[folds > 0] <- 0 
  
  # Function to suppress the warning of missings in the response
  h <- function(w){
    if( any( grepl( "response contains missing values;", w) ) ) 
      invokeRestart( "muffleWarning" )  
  }
  
  ### get yind in long format
  yindLong <- object$yind
  if(!any(class(object) == "FDboostLong")){
    yindLong <- rep(object$yind, each = nObs)
  }
  ### compute ("length of each trajectory")^-1 in the response
  ### more precisely ("sum of integration weights")^-1 is used
  if(length(object$yind) > 1){
    # lengthTi1 <- 1/tapply(yindLong[!is.na(response)], id[!is.na(response)], function(x) max(x) - min(x))
    lengthTi1 <- 1/tapply(intWeights, id, function(x) sum(x))
    if(any(is.infinite(lengthTi1))) lengthTi1[is.infinite(lengthTi1)] <- max(lengthTi1[!is.infinite(lengthTi1)])
  }else{
    lengthTi1 <- rep(1, l = length(response))
  }
  
  ###### Function to fit the model
  # function working with FDboost, thus the smooth offset is recalculated in each model
  dummyfct <- function(weights, oobweights) {
    
    # create data frame for model fit and use the weights vector for the CV
    dathelp <- object$data
    nameyind <- attr(object$yind, "nameyind")
    dathelp[[nameyind]] <- object$yind
    
    if(!any(class(object) == "FDboostLong") & !any(class(object) == "FDboostScalar")){
      dathelp[[object$yname]] <- matrix(object$response, ncol = Gy)
    }else{
      dathelp[[object$yname]] <- object$response
    }

    call <- object$callEval
    call$data <- dathelp
    
    # use weights of training data expanded by id to suitable length
    call$weights <- weights[id]
    
    # use call$numInt of original model fit, as weights contains only resampling weights
 
    # Using the offset of object with the following settings 
    # call$control <- boost_control(risk="oobag")
    # call$oobweights <- oobweights[id]
    if(refitSmoothOffset == FALSE && is.null(call$offset) ){
      if(!any(class(object) == "FDboostLong")){
        call$offset <- matrix(object$offset, ncol = Gy)[1, ]
      }else{
        call$offset <- object$offset
      }
    } 
    # the model is the same for  
    # mod <- object$update(weights = weights, oobweights = oobweights) # (cvrisk) 
    # and
    # mod <- withCallingHandlers(suppressMessages(eval(call)), warning = h)
    # and then the risk can be computed by
    # risk <- mod$risk()[grid] 
    
    ## compute the model by FDboost() - the offset is computed on learning sample
    mod <- withCallingHandlers(suppressMessages(eval(call)), warning = h) # suppress the warning of missing responses    
    mod <- mod[max(grid)]
        
    # compute weights for standardizing risk, mse, ...
    oobwstand <- lengthTi1[id]*oobweights[id]*intWeights*(1/sum(oobweights))
    
    ############# compute risk, mse, relMSE and mrd
    # get risk function of the family
    riskfct <- get("family", environment(mod$update))@risk
    
    ####################
    ### compute risk and mse without integration weights, like in cvrisk
    risk0 <- sapply(grid, function(g){riskfct(response, mod[g]$fitted(), 
                                               w = oobweights[id])}) / sum(oobweights[id])
    
    mse0 <- simplify2array(mclapply(grid, function(g){
      sum( ((response - mod[g]$fitted())^2*oobweights[id]), na.rm = TRUE )
    }, mc.cores=1) ) /sum(oobweights[id])
    ####################
    
    # oobweights using riskfct() like in mboost, but with different weights!
    risk <- sapply(grid, function(g){riskfct( response, mod[g]$fitted(), w=oobwstand)})
    
    ### mse (mean squared error) equals risk in the case of familiy=Gaussian()
    mse <- simplify2array(mclapply(grid, function(g){
      sum( ((response - mod[g]$fitted())^2*oobwstand), na.rm = TRUE )
    }, mc.cores=1) )
    
    #     ### mse2 equals mse in the case of equal grids without missings at the ends
    #     mse2 <- simplify2array(mclapply(grid, function(g){
    #       sum( ((response - mod[g]$fitted())^2*oobweights[id]*intWeights), na.rm=TRUE )
    #     }, mc.cores=1) ) / (sum(oobweights)* (max(mod$yind)-min(mod$yind) ) )
    
    ### compute overall mean of response in learning sample
    meanResp <- sum(response*intWeights*lengthTi1[id]*weights[id], na.rm = TRUE) / sum(weights)
    
    # # compute overall mean of response in whole sample
    # meanResp <- sum(response*intWeights*lengthTi1[id], na.rm=TRUE) / nObs
    
    ### compute relative mse
    relMSE <- simplify2array(mclapply(grid, function(g){
      sum( ((response - mod[g]$fitted())^2*oobwstand), na.rm = TRUE ) /
        sum( ((response - meanResp)^2*oobwstand), na.rm = TRUE )
    }, mc.cores=1) )
    
    ### mean relative deviation
    resp0 <- response
    resp0[abs(resp0) <= mrdDelete | round(resp0, 1) == 0] <- NA
    mrd <- simplify2array(mclapply(grid, function(g){
      sum( abs(resp0 - mod[g]$fitted())/abs(resp0)*oobwstand, na.rm = TRUE )
    }, mc.cores=1) )

    mrd0 <- simplify2array(mclapply(grid, function(g){
              sum( abs(resp0 - mod[g]$fitted())/abs(resp0)*oobweights[id], na.rm = TRUE)
                 }, mc.cores=1) ) / sum(oobweights[id])

    
    rm(resp0, meanResp)
    
    ####### prediction for all observations, not only oob! 
    # the predictions are in a long vector for all model types (regular, irregular, scalar)
    predGrid <- predict(mod, aggregate = "cumsum", toFDboost = FALSE)
    predGrid <- predGrid[ , grid] # save vectors of predictions for grid in matrix

    if(showProgress) cat(".")
    
    ## user-specified function to use on FDboost-object 
    if(! is.null(fun) ){
      fun_ret <- fun(mod)
    }else{
      fun_ret <- NULL
    }

    return(list(risk = risk, predGrid = predGrid, # predOOB = predOOB, respOOB = respOOB, 
                mse = mse, relMSE = relMSE, mrd = mrd, risk0 = risk0, mse0 = mse0, mrd0 = mrd0, 
                mod = mod, fun_ret = fun_ret))  
  } 
  
  ### computation of models on partitions of data
  if(Sys.info()["sysname"]=="Linux"){
    modRisk <- mclapply(1:ncol(folds),
                        function(i) dummyfct(weights = folds[, i],
                                             oobweights = OOBweights[, i]), ...)
  }else{
    modRisk <- mclapply(1:ncol(folds),
                        function(i) dummyfct(weights = folds[, i], 
                                             oobweights = OOBweights[, i]), mc.cores = 1)
  }
  
  # str(modRisk, max.level=2)
  # str(modRisk, max.level=5)
  
  # check whether model fit worked in all iterations
  modFitted <- sapply(modRisk, function(x) class(x) == "list")
  if(any(!modFitted)){
    
    # stop() or warning()?
    if(sum(!modFitted) > sum(modFitted)) warning("More than half of the models could not be fitted.")
    
    warning("Model fit did not work in fold ", paste(which(!modFitted), collapse = ", "))
    modRisk <- modRisk[modFitted]
    OOBweights <- OOBweights[,modFitted]   
    folds <- folds[,modFitted]
  } 
  
  ####### restructure the results
  
  ## get the out-of-bag risk
  oobrisk <- t(sapply(modRisk, function(x) x$risk)) 
  ## get out-of-bag mse
  oobmse <- t(sapply(modRisk, function(x) x$mse))  
  ## get out-of-bag relMSE
  oobrelMSE <- t(sapply(modRisk, function(x) x$relMSE))  
  ## get out-of-bag mrd
  oobmrd <- t(sapply(modRisk, function(x) x$mrd))
  
  colnames(oobrisk) <- colnames(oobmse) <- colnames(oobrelMSE) <- colnames(oobmrd) <- grid
  rownames(oobrisk) <- rownames(oobmse) <- rownames(oobrelMSE) <- rownames(oobmrd) <- which(modFitted)

  ## get out-of-bag risk without integration weights
  oobrisk0 <- t(sapply(modRisk, function(x) x$risk0))
  ## get out-of-bag mse without integration weights
  oobmse0 <- t(sapply(modRisk, function(x) x$mse0))
  ## get out-of-bag mrd without integration weights
  oobmrd0 <- t(sapply(modRisk, function(x) x$mrd0))
  
  colnames(oobrisk0) <- colnames(oobmse0) <- colnames(oobmrd0)  <- grid
  rownames(oobrisk0) <- rownames(oobmse0) <- rownames(oobmrd0) <- which(modFitted)
  
  ############# check for folds with extreme risk-values at the global median
  riskOptimal <- oobrisk[ , which.min(apply(oobrisk, 2, median))]
  bound <- median(riskOptimal) + 1.5*(quantile(riskOptimal, 0.75) - quantile(riskOptimal, 0.25))
  
  # fold equals curve if type="curves"
  if(any(riskOptimal>bound)){
    message("Fold with high values in oobrisk (median is ", round(median(riskOptimal), 2), "):")
    message(paste("In fold ", which(riskOptimal > bound), ": " ,
                  round(riskOptimal[which(riskOptimal > bound)], 2), collapse = ",  ", sep = "" )  )  
  }
  
  ## only makes sense for type="curves" with leaving-out one curve per fold!!
  if(grepl( "curves", type)){
    # predict response for all mstops in grid out of bag
    # predictions for each response are in a vector!
    oobpreds0 <- lapply(modRisk, function(x) x$predGrid)
    oobpreds <- matrix(nrow = nrow(oobpreds0[[1]]), ncol = ncol(oobpreds0[[1]]))
    
    if(any(class(object) == "FDboostLong")){
      for(i in 1:length(oobpreds0)){ # i runs over observed trajectories, i.e. over id 
        oobpreds[id == i, ]  <- oobpreds0[[i]][id == i, ] 
      }
    }else{
      for(j in 1:length(oobpreds0)){
        oobpreds[folds[ , j] == 0] <- oobpreds0[[j]][folds[ , j] == 0] 
      }
    }
    colnames(oobpreds) <- grid
    rm(oobpreds0)
    
  }else{
    oobpreds <- NULL
  }
    
  #   # alternative OOB-prediction: works for general folds not only oob
  #   predOOB <- lapply(modRisk, function(x) x$predOOB)
  #   predOOB <- do.call('rbind', predOOB)
  #   colnames(predOOB) <- grid
  #   respOOB <- lapply(modRisk, function(x) x$respOOB)
  #   respOOB <- do.call('c', respOOB)
  #   indexOOB <- lapply(modRisk, function(x) attr(x$respOOB, "curves"))
  #   if(is.null(object$id)){
  #     indexOOB <- lapply(indexOOB, function(x) rep(x, times=Gy) )
  #     indexOOB <- unlist(indexOOB)
  #   }else{
  #     indexOOB <- names(unlist(indexOOB))[unlist(indexOOB)]
  #   }
  #   attr(respOOB, "index") <- indexOOB
  
  coefCV <- list()
  predCV <- list()

  if(getCoefCV){
    
    if(riskopt == "median"){
      optimalMstop <- grid[which.min(apply(oobrisk, 2, median))]
    }else{
      optimalMstop <- grid[which.min(apply(oobrisk, 2, mean))]
    }  
    
    attr(coefCV, "risk") <- paste("minimize", riskopt, "risk")
    
    ### estimates of coefficients
    timeHelp <- seq(min(modRisk[[1]]$mod$yind), max(modRisk[[1]]$mod$yind), l = 40)
    for(l in 1:length(modRisk[[1]]$mod$baselearner)){
      # estimate the coefficients for the model of the first fold
      my_coef <- coef(modRisk[[1]]$mod[optimalMstop], 
           which = l, n1 = 40, n2 = 20, n3 = 15, n4 = 10)$smterms[[1]]
      if(is.null(my_coef)){
        my_coef <- list(0)
        my_coef$dim <- 0
      }
      coefCV[[l]] <- my_coef
      
      ## no %X% with several levels in the coefficients
      if(is.null(coefCV[[l]]$numberLevels)){
        attr(coefCV[[l]]$value, "offset") <- NULL # as offset is the same within one model
        
        # add estimates for the models of the other folds
        coefCV[[l]]$value <- lapply(1:length(modRisk), function(g){
          ret <- coef(modRisk[[g]]$mod[optimalMstop], 
                      which = l, n1 = 40, n2 = 20, n3 = 15, n4 = 10)$smterms[[1]]$value
          #         if(l==1){
          #           coefCV[[l]]$offset[g,] <- modRisk[[g]]$mod$predictOffset(time=timeHelp)
          #         }
          attr(ret, "offset") <- NULL # as offset is the same within one model
          return(ret)
        })
      }else{
        ## %X% with numberLevels coefficient values in a list
        ## lapply(1:coefCV[[l]]$numberLevels, function(x) coefCV[[l]][[x]]$value)
        for(j in 1:coefCV[[l]]$numberLevels){
          coefCV[[l]][[j]]$value <- lapply(1:length(modRisk), function(g){
            ret <- coef(modRisk[[g]]$mod[optimalMstop], 
                        which = l, n1 = 40, n2 = 20, n3 = 15, n4 = 10)$smterms[[1]][[j]]$value
            attr(ret, "offset") <- NULL # as offset is the same within one model
            return(ret)
          })
        }
      } # end else for numberLevels
      
    }
    
    ## predict offset
    offset <- sapply(1:length(modRisk), function(g){
      # offset is vector of length yind or numeric of length 1 for constant offset
      ret <- modRisk[[g]]$mod$predictOffset(time = timeHelp)
      if( length(ret) == 1 & length(object$yind) > 1 ) ret <- rep(ret, length(timeHelp))
      return(ret)
    })
    
    attr(coefCV, "offset") <- offset
    
    niceNames <- c("offset", lapply(coefCV, function(x) x$main))
    
    ### predictions of terms based on the coefficients for each model
    # only makes sense for type="curves" with leaving-out one curve per fold!!
    if(grepl("curves", type)){
      for(l in 1:(length(modRisk[[1]]$mod$baselearner)+1)){
        predCV[[l]] <- t(sapply(1:length(modRisk), function(g){
          if(l == 1){ # save offset of model
            # offset is vector of length yind or numeric of length 1 for constant offset
            ret <- modRisk[[g]]$mod[optimalMstop]$predictOffset(object$yind) 
            # regular data or scalar response
            if(!any(class(object) == "FDboostLong")){
              if( length(ret) == 1 ) ret <- rep(ret, modRisk[[1]]$mod$ydim[2])
            # irregular data
            }else{ 
              if( length(ret) == 1 ){ ret <- rep(ret, sum(object$id == g)) }else{ ret <- ret[object$id == g] }
            }
          }else{ # other effects
            ret <- predict(modRisk[[g]]$mod[optimalMstop], which = l-1) # model g
            if(!(l-1) %in% selected(modRisk[[g]]$mod[optimalMstop]) ){ # effect was never chosen
              if(!any(class(object) == "FDboostLong")){
                ret <- matrix(0, ncol=modRisk[[1]]$mod$ydim[2], nrow=modRisk[[1]]$mod$ydim[1])
              }else{
                ret <- matrix(0, nrow = length(object$id), ncol=1)
              } 
            }
            if(!any(class(object) == "FDboostLong")){
              ret <- ret[g,] # save g-th row = preds for g-th observations
            }else{
              ret <- ret[object$id == g] # save preds of g-th observations
            } 
          }
          return(ret)
        }))
        names(predCV)[l] <- niceNames[l]
        #       matplot(modRisk[[1]]$mod$yind, t(predCV[[l]]), type="l", 
        #               main=names(predCV)[l], xlab=attr(modRisk[[1]]$mod$yind, "nameyind"), ylab="coef")
      } 
    }
    
  } # end of if(getCoefCV)
  
  ret <- list(response = response, yind = object$yind, id = object$id,
              folds = folds, grid=grid,
              coefCV = coefCV,
              predCV = predCV,
              oobpreds = oobpreds, 
              oobrisk = oobrisk, 
              oobriskMean = colMeans(oobrisk),
              oobmse = oobmse,
              oobrelMSE = oobrelMSE,
              oobmrd = oobmrd,
              oobrisk0 = oobrisk0, 
              oobmse0 = oobmse0,
              oobmrd0 = oobmrd0, 
              format = if(any(class(object) == "FDboostLong")) "FDboostLong" else "FDboost", 
              fun_ret = if(is.null(fun)) NULL else lapply(modRisk, function(x) x$fun_ret) )
  
  rm(modRisk)

  attr(ret, "risk") <- object$family@name
  attr(ret,  "call") <- deparse(object$call)
  
  class(ret) <- "validateFDboost"
  
  return(ret) 
}



#' @rdname plot.validateFDboost
#' @method mstop validateFDboost
#' @export
#' 
# Function to extract the optimal stopping iteration
mstop.validateFDboost <- function(object, riskopt=c("mean", "median"), ...){
  
  dots <- list(...)
  riskopt <- match.arg(riskopt)
  
  if(riskopt=="median"){
    riskMedian <- apply(object$oobrisk, 2, median) 
    mstop <- object$grid[which.min(riskMedian)]
    attr(mstop, "risk") <- "minimize median risk"
  }else{
    riskMean <- colMeans(object$oobrisk)
    mstop <- object$grid[which.min(riskMean)]
    attr(mstop, "risk") <- "minimize mean risk"
  }  
  return(mstop) 
}


#' @rdname plot.validateFDboost
#' @method print validateFDboost
#' @export
#' 
# Function to print an object of class validateFDboost 
print.validateFDboost <- function(x, ...){
  
  cat("\n\t Cross-validated", attr(x, "risk"), "\n\t", attr(x,  "call"), "\n\n")
  print(colMeans(x$oobrisk))
  cat("\n\t Optimal number of boosting iterations:", mstop(x), 
      "\n")
  return(invisible(x))  
}


#' Methods for objects of class validateFDboost
#' 
#' Methods for objects that are fitted to determine the optimal mstop and the 
#' prediction error of a model fitted by FDboost.
#' 
#' @param object object of class \code{validateFDboost} 
#' @param riskopt how the risk is minimized to obtain the optimal stopping iteration; 
#' defaults to the mean, can be changed to the median.
#' 
#' @param x an object of class \code{validateFDboost}. 
#' @param modObject if the original model object of class \code{FDboost} is given 
#' predicted values of the whole model can be compared to the predictions of the cross-validated models
#' @param predictNA should missing values in the response be predicted? Defaults to \code{FALSE}. 
#' @param which In the case of \code{plotPredCoef()} the subset of base-learners to take into account for plotting. 
#' In the case of \code{plot.validateFDboost()} the diagnostic plots that are given 
#' (1: empirical risk per fold as a funciton of the boosting iterations, 
#' 2: empirical risk per fold, 3: MRD per fold, 
#' 4: observed and predicted values, 5: residuals; 
#' 2-5 for the model with the optimal number of boosting iterations). 
#' @param names.arg names of the observed curves
#' @param ask defaults to \code{TRUE}, ask for next plot using \code{par(ask = ask)}  ? 
#' @param pers plot coefficient surfaces as persp-plots? Defaults to \code{TRUE}.
#' @param commonRange, plot predicted coefficients on a common range, defaults to \code{TRUE}.
#' @param showQuantiles plot the 0.05 and the 0.95 Quantile of coefficients in 1-dim effects.
#' @param showNumbers show number of curve in plot of predicted coefficients, defaults to \code{FALSE}
#' @param terms logical, defaults to \code{TRUE}; plot the added terms (default) or the coefficients?
#' @param probs vector of quantiles to be used in the plotting of 2-dimensional coefficients surfaces,
#' defaults to \code{probs = c(0.25, 0.5, 0.75)}
#' @param ylab label for y-axis
#' @param xlab label for x-axis
#' @param ylim values for limits of y-axis
#' @param ... additional arguments passed to callies.
#' 
#' @details The function \code{mstop.validateFDboost} extracts the optimal mstop by minimizing the 
#' mean (or the median) risk. 
#' \code{plot.validateFDboost} plots cross-validated risk, RMSE, MRD, measured and predicted values 
#' and residuals as determined by \code{validateFDboost}. The function \code{plotPredCoef} plots the 
#' coefficients that were estimated in the folds - only possible if the argument getCoefCV is \code{TRUE} in 
#' the call to \code{validateFDboost}. 
#' 
#' @aliases mstop.validateFDboost
#' 
#' @method plot validateFDboost
#' 
#' @export
plot.validateFDboost <- function(x, riskopt=c("mean", "median"), 
                                 ylab = attr(x, "risk"), xlab = "Number of boosting iterations",
                                 ylim = range(x$oobrisk), 
                                 which = 1, 
                                 modObject = NULL, predictNA = FALSE, 
                                 names.arg = NULL, ask=TRUE, ...){
  
  # get the optimal mstop
  riskopt <- match.arg(riskopt)
  mopt <- mstop(x, riskopt=riskopt)
  # get the position in the grid of mopt
  mpos <- which(x$grid == mopt)
  
  if(length(which) > 1) par(ask = ask)
  
  if(1 %in% which){
    # Plot the cross validated risk
    ##ylim <- c(min(x$oobrisk), quantile(x$oobrisk, 0.97))
    ## ylim <- c(min(x$oobrisk), max(x$oobrisk))
    matplot(colnames(x$oobrisk), t(x$oobrisk), type="l", col="lightgrey", lty=1,
            ylim=ylim, 
            xlab=xlab, ylab=ylab,
            main=attr(x$folds, "type"), 
            sub=attr(mopt, "risk"))
    
    if(riskopt == "mean"){
      riskMean <- colMeans(x$oobrisk)
      lines(colnames(x$oobrisk), riskMean, lty=1)
      mOptMean <- x$grid[which.min(riskMean)]
      lines(c(mOptMean, mOptMean), 
            c(min(c(0, ylim[1] * ifelse(ylim[1] < 0, 2, 0.5))), 
              riskMean[paste(mOptMean)]), lty = 2)
      legend("topright", legend=paste(c(mOptMean)),
             lty=c(2), col=c("black"))
      
    }

    if(riskopt == "median"){
      riskMedian <- apply(x$oobrisk, 2, median) 
      lines(colnames(x$oobrisk), riskMedian, lty=1)
      mOptMedian <- x$grid[which.min(riskMedian)]
      lines(c(mOptMedian, mOptMedian), 
            c(min(c(0, ylim[1] * ifelse(ylim[1] < 0, 2, 0.5))), 
              riskMedian[paste(mOptMedian)]), lty = 2)
      
      legend("topright", legend=paste(c(mOptMedian)),
             lty=c(2), col=c("black"))
    }
    
  }
  
  if(any(c(2,3) %in% which)){
    # Plot RMSE and MRD for optimal mstop
    if(is.null(names.arg)){
      names.arg <- seq(along=x$oobrisk[,mpos])
    }
    stopifnot(length(names.arg)==length(x$oobrisk[,mpos]))
  }
  
  # plot risk
  if(2 %in% which){
    barplot(x$oobrisk[,mpos], main="risk", names.arg = names.arg, las=2) # 
    #abline(h=x$rmse[mpos], lty=2)
  }
  
  # plot MRD
  if(3 %in% which){
    barplot(x$oobmrd[,mpos], main="MRD", names.arg = names.arg, las=2)
    #abline(h=x$mrd[mpos], lty=2)
  }
    
  # Plot the predictions for the optimal mstop
  if(4 %in% which | 5 %in% which){
    
    if(!is.null(x$oobpreds)){
      response <- x$response
      pred <- x$oobpreds[,mpos]

      if(!predictNA){
        pred[is.na(response)] <- NA
      }
      
      predMat <- pred
      responseMat <- response
      
      if(x$format == "FDboost") predMat <- matrix(pred, ncol=length(x$yind))
      if(x$format == "FDboost") responseMat <- matrix(response, ncol=length(x$yind))
      
      if(4 %in% which){
        ylim <- range(response, pred, na.rm = TRUE)
        funplot(x$yind, responseMat, id=x$id, lwd = 1, pch = 1, ylim = ylim,  
                ylab = "", xlab = attr(x$yind, "nameyind"), 
                main="Observed and Predicted Values", ...)
        funplot(x$yind, predMat, id=x$id, lwd = 2, pch = 2, add = TRUE, ...)
        posLegend <- "topleft"
        legend(posLegend, legend=c("observed","predicted"), col=1, pch=1:2)  
      }
      
      if(5 %in% which){
        # Plot residuals for the optimal mstop
        funplot(x$yind, responseMat-predMat, id=x$id, ylab = "", xlab = attr(x$yind, "nameyind"), 
                main="Residuals", ...)
        abline(h=0, lty=2, col="grey")
      }
    }
  }
  
  # Plot coefficients
  
  #   # example: plot coeficients of 5th effect for folds 1-4 each for the optimal mstop
  #   plot(modObject, which=5, n1 = 50, n2 = 20, n3 = 15, n4 = 10, levels=seq(-.4, 1.4, by=0.4))
  #   contour(x$coefCV[[1]][[5]]$x,
  #           x$coefCV[[1]][[5]]$y,
  #           t(x$coefCV[[1]][[5]]$value[[mpos]]), lty=2, add=TRUE, col=2, levels=seq(-.4, 1.4, by=0.4))
  #   contour(x$coefCV[[2]][[5]]$x,
  #           x$coefCV[[2]][[5]]$y,
  #           t(x$coefCV[[2]][[5]]$value[[mpos]]), lty=2, add=TRUE, col=2, levels=seq(-.4, 1.4, by=0.4))
  #   
  #   image(x$coefCV[[1]][[5]]$value[[mpos]], type="l", lty=1)
  
  #   matplot(x$coefCV[[1]][[5]]$value[[mpos]], type="l", lty=1)
  #   matplot(x$coefCV[[2]][[5]]$value[[mpos]], type="l", lty=2, add=TRUE)
  #   matplot(x$coefCV[[3]][[5]]$value[[mpos]], type="l", lty=3, add=TRUE)
  #   matplot(x$coefCV[[4]][[5]]$value[[mpos]], type="l", lty=4, add=TRUE)
  
  par(ask = FALSE)  
}


#' @rdname plot.validateFDboost
#' @export
#' 
plotPredCoef <- function(x, which = NULL, pers = TRUE,
                         commonRange = TRUE, showNumbers = FALSE, showQuantiles = TRUE,
                         ask = TRUE, 
                         terms = TRUE,         
                         probs = c(0.25, 0.5, 0.75), # quantiles of variables to use for plotting
                         ylim = NULL, ...){
  
  stopifnot(any(class(x) == "validateFDboost"))

  if(is.null(which)) which <- 1:length(x$coefCV)
  
  if(length(which) > 1) par(ask = ask)
  
  if(terms){
    
    if(all(which == 1:length(x$coefCV))){
      which <- 1:(length(x$coefCV)+1)
    }else{
      which <- which + 1 
    }
    
    if(length(x$predCV) == 0){
      warning("no bootstrapped predcition, set terms = FALSE, to plot bootstrapped coefficients.")
      return(NULL)
    }
    
    if(commonRange & is.null(ylim)){ 
      ylim <- range(x$predCV[which])
    }
    
    ### loop over base-learners
    for(l in which){
      
      if( x$format == "FDboostLong" ){
        
        funplot(x$yind, unlist(x$predCV[[l]]), id=x$id, col="white", 
                main=names(x$predCV)[l], xlab=attr(x$yind, "nameyind"), ylab="coef", ylim=ylim, ...)
        for(i in 1:length(x$predCV[[l]])){
          lines(x$yind[x$id==i], x$predCV[[l]][[i]], lwd=1, col=i)
          if(showNumbers){
            points(x$yind[x$id==i], x$predCV[[l]][[i]], type="p", pch=paste0(i))
          }
        }
        
      }else{
        funplot(x$yind, unlist(x$predCV[[l]]), 
                id=x$id, type="l", 
                main=names(x$predCV)[l], xlab=attr(x$yind, "nameyind"), ylab="coef", ylim=ylim, ...)
        
        if(showNumbers){
          funplot(x$yind, unlist(x$predCV[[l]]), id=x$id, pch=paste0(x$id), add=TRUE, type="p")
        }
      } 
    } # end for-loop 
    
  }else{ # plot coefficients
    
    if(commonRange & is.null(ylim)){
      if(length(x$yind)>1){
        if(!any(sapply(lapply(x$coefCV[which], function(x) x$value), is.null))){
          ylim <- range(lapply(x$coefCV[which], function(x) x$value)) 
        }else{
          ## in case of composed base-learners with %X% the ylim is determinded using the first value
          ylim <- range(lapply(x$coefCV[which], function(x) x[[1]]$value)) 
        }
      }else{
        ylim <- range(lapply(x$coefCV[which], function(x) x$value))
      }
      if(any(is.infinite(ylim))) ylim <- NULL
    }
    
    for(l in which){ # loop over effects
      
      # coef() of a certain term
      temp <- x$coefCV[l][[1]]
      
      plot_bootstrapped_coef(temp = temp, l = l, 
                             offset = attr(x$coefCV, "offset"), yind = x$yind, 
                             pers = pers, 
                             showNumbers = showNumbers, showQuantiles = showQuantiles, 
                             probs = probs, ylim = ylim, ...)
      
      
    } # end loop over effects
  }
  par(ask=FALSE)
}


## helper function to plot bootstrapped coefficients with basis plots 
plot_bootstrapped_coef <- function(temp, l,  
                                   offset, yind, 
                                   pers, showNumbers, showQuantiles, probs, ylim, ...){
  
  ### Get further arguments passed to the different plot-functions
  dots <- list(...)
  
  if(!is.null(ylim) && all(ylim == 0)){
    ylim <- c(-0.1, 0.1)
  }
  
  getArguments <- function(x, dots=dots){
    if(any(names(dots) %in% names(x))){
      dots[names(dots) %in% names(x)]
    }else list()
  }
  
  argsPlot <- getArguments(x = c(formals(graphics::plot.default), par()), dots = dots)
  argsMatplot  <- getArguments(x = c(formals(graphics::matplot), par()), dots = dots)
  argsFunplot  <- getArguments(x = c(formals(funplot), par()), dots = dots)
  
  argsImage <- getArguments(x=c(formals(graphics::plot.default), 
                                formals(graphics::image.default)), dots=dots)
  dotsContour <- dots
  dotsContour$col <- "black"
  argsContour <- getArguments(x=formals(graphics::contour.default), dots=dotsContour)
  
  argsPersp <- getArguments(x=formals(getS3method("persp", "default")), dots = dots)
  
  plotWithArgs <- function(plotFun, args, myargs){        
    args <- c(myargs[!names(myargs) %in% names(args)], args)        
    do.call(plotFun, args)            
  }  
  
  
  ## helper funciton to plot curves 
  plot_curves <- function(x_i, y_i, xlab_i, main_i, ylim_i, ...){
    
    plotWithArgs(matplot, args = argsMatplot, 
                 myargs = list(x = x_i, y = y_i, type = "l", xlab = xlab_i, ylab = "coef", 
                             main = main_i, ylim = ylim_i, lty = 1,  
                             col = rgb(0.6,0.6,0.6, alpha = 0.5)))
    
    if(showNumbers){
      matplot(x_i, y_i, add = TRUE, col = rgb(0.6,0.6,0.6, alpha = 0.5), ...)
    }
    
    if(showQuantiles){ 
      lines(x_i, rowMeans(y_i), col = 1, lwd = 2)
      lines(x_i, apply(y_i, 1, quantile, 0.95, na.rm = TRUE), col = 2, lwd = 2, lty = 2)
      lines(x_i, apply(y_i, 1, quantile, 0.05, na.rm = TRUE), col = 2, lwd = 2, lty = 2)
    }
    
  }
  
  
  if(!is.null(temp$numberLevels)){
    temp <- temp[[1]]
    warning("Of the composed base-learner ", l, " only the first effect is plotted.")
  }
  
  # set the range for each effect individually
  if(FALSE) ylim <- range(temp$value)
  
  ## plot the estimated offset for functional response 
  if(l == 1 && length(yind) > 1){
    myMat <- offset ## attr(x$coefCV, "offset")
    timeHelp <- seq(min(yind), max(yind), length = nrow(myMat))

    plot_curves(x_i = timeHelp, y_i = myMat, xlab_i = temp$xlab, 
                main_i = temp$main, ylim_i = NULL)
  }
  
  
  if(temp$dim == 0){

    # intercept and offset of model with scalar response 
    temp$value[sapply(temp$value, function(x) is.null(x))] <- 0
    y_i <- unlist(temp$value)
    
    offset[sapply(offset, function(x) is.null(x))] <- 0
    offset_i <- offset
    
    x_i <-rep(0, length(y_i))

    plotWithArgs(plot, args = argsPlot, 
                 myargs = list(x = x_i, y = y_i + offset_i, xlab = "", ylab = "coef", 
                               xaxt = "n", main = "offset + intercept",   
                               col = rgb(0.6,0.6,0.6, alpha = 0.5)))
    
    if(showNumbers){
      matplot(x_i, y_i + offset_i, col = rgb(0.6,0.6,0.6, alpha = 0.5), add = TRUE)
    }
    
    if(showQuantiles){ 
      points(0, mean(y_i + offset_i), col = 1, lwd = 2)
      points(0, quantile(y_i + offset_i, 0.95, na.rm = TRUE), col = 2, lwd = 2, lty = 2)
      points(0, quantile(y_i + offset_i, 0.05, na.rm = TRUE), col = 2, lwd = 2, lty = 2)
    }
    
  }
  
  if(temp$dim == 1){
    # impute vector of 0 if effect was never chosen
    temp$value[sapply(temp$value, function(x) length(x)==1)] <- list(rep(0, 40))
    myMat <- sapply(temp$value, function(x) x) # as one matrix
    
    plot_curves(x_i = temp$x, y_i = myMat, xlab_i = temp$xlab, 
                main_i = temp$main, ylim_i = ylim)
    
  }
  
  if(temp$dim == 2){
    
    if(!is.factor(temp$x)){
      quantx <- quantile(temp$x, probs=probs, type=1)
    } else quantx <- temp$x
    
    quanty <- quantile(temp$y, probs=probs, type=1)
    
    # set lower triangular matrix to NA for historic effect
    if(grepl("bhist", temp$main)){
      for(k in 1:length(temp$value)){
        temp$value[[k]][temp$value[[k]]==0] <- NA
      }
    }
    
    if(!is.factor(temp$x)){ # temp$x is metric
      
      # impute matrix of 0 if effect was never chosen
      temp$value[sapply( temp$value, function(x){ is.null(dim(x)) || dim(x)[2]==1 })] <- list(matrix(0, ncol=20, nrow=20))
      
      # plot coefficient surfaces at different pointwise quantiles
      if(pers){ 
        matvec <- sapply(temp$value, c)
        for(k in 1:length(probs)){
          
          tempZ <- matrix(apply(matvec, 1, quantile, probs=probs[k], na.rm=TRUE), ncol=length(temp$x))
          
          plotWithArgs(persp, args=argsPersp, 
                       myargs=list(x=temp$x, y=temp$y, z=tempZ,
                                   ticktype="detailed", theta=30, phi=30,
                                   xlab=paste("\n", temp$xlab), ylab=paste("\n", temp$ylab), 
                                   zlab=paste("\n", "coef"), 
                                   zlim=if(any(is.null(ylim))) range(matvec, na.rm=TRUE) else ylim,  
                                   main=paste(temp$main, " at ", probs[k]*100, "%-quantile", sep=""), 
                                   col=getColPersp(tempZ)))
        
        }
        
      }else{ # do 2-dim plots

        # for(j in 1:length(quanty)){ 
        #   
        #   myCol <- sapply(temp$value, function(x) x[, quanty[j]==temp$y]) # first column
        #   
        #   plot_curves(x_i = temp$y, y_i = myCol, xlab_i = temp$ylab, 
        #               main_i = paste(temp$main, " at ", probs[j]*100, "% of ", temp$xlab, sep = ""), 
        #               ylim_i = ylim)
        #   
        # } # end loop over quanty
        # 
        # for(j in 1:length(quantx)){  
        #   myRow <- sapply(temp$value, function(x) x[quantx[j]==temp$x, ]) # first column
        #   
        #   plot_curves(x_i = temp$x, y_i = myRow, xlab_i = temp$xlab, 
        #               main_i = paste(temp$main, " at ", probs[j]*100, "% of ", temp$ylab, sep = ""), 
        #               ylim_i = ylim)
        #     
        # }
        
        matvec <- sapply(temp$value, c)
        for(k in 1:length(probs)){
          
          tempZ <- matrix(apply(matvec, 1, quantile, probs=probs[k], na.rm=TRUE), ncol=length(temp$x))
          
          plotWithArgs(image, args=argsImage,
                       myargs=list(x=temp$y, y=temp$x, z=t(tempZ), xlab=paste("\n", temp$xlab), 
                                   ylab=paste("\n", temp$ylab), zlim=c(min(matvec, na.rm=TRUE),
                                                                       max(matvec, na.rm=TRUE)),
                                   main=paste(temp$main, " at ", probs[k]*100, "%-quantile", sep=""), 
                                   col = heat.colors(length(temp$x)^2) 
                                   )
          )
          plotWithArgs(contour, args=argsContour,
                       myargs=list(temp$y, temp$x, z=t(tempZ), col = "black", add = TRUE))
          
        }
         
      } # end else
      
    }else{ # temp$x is factor
      
      for(j in 1:length(quantx)){ 
        
        # impute matrix of 0 if effect was never chosen
        temp$value[sapply(temp$value, function(x) is.null(dim(x)))] <- list(matrix(0, ncol=20, nrow=length(quantx)))
        
        if(is.null(temp$z)){
          # myRow <- sapply(temp$value, function(x) x[quantx[j]==temp$x, ]) # first column
          myRow <- t(temp$value[[which(quantx[j]==temp$x)]])
                     
          plot_curves(x_i = temp$y, y_i = myRow, xlab_i = temp$ylab, 
                      main_i = paste(temp$main, " at ", temp$xlab,"=" ,quantx[j], sep = ""), 
                      ylim_i = ylim)
          
        }else{
          quantz <- temp$z
          myRow <- sapply(temp$value, function(x) x[quantx[j]==temp$x & quantz[j]==temp$z, ]) # first column
          
          plot_curves(x_i = temp$y, y_i = myRow, xlab_i = temp$ylab, 
                      main_i = paste(temp$main, " at ", temp$xlab, "=" , quantx[j], ", " , 
                                     temp$zlab, "=", quantz[j], sep = ""), 
                      ylim_i = ylim)
        } 
      }
    }
    
  } # end if(temp$dim == 2)
  
}



#' @rdname applyFolds
#' @export
## wrapper for cvrisk of mboost, specifying folds on the level of curves
cvrisk.FDboost <- function(object, folds = cvLong(id=object$id, weights=model.weights(object)),
                           grid = 1:mstop(object),
                           papply = mclapply,
                           fun = NULL, corrected = TRUE, mc.preschedule = FALSE, ...){
  
  if(!length(unique(object$offset)) == 1) message("The smooth offset is fixed over all folds.")
  
  names_bl <- names(object$baselearner)
  if(any(grepl("brandomc", names_bl))) message("For brandomc, the transformation matrix Z is fixed over all folds.")
  if(any(grepl("bolsc", names_bl))) message("For bolsc, the transformation matrix Z is fixed over all folds.")
  if(any(grepl("bbsc", names_bl))) message("For bbsc, the transformation matrix Z is fixed over all folds.")
  
  class(object) <- "mboost"
  
  ret <- cvrisk(object = object, folds = folds,
                grid = grid,
                papply = papply,
                fun = fun, corrected = corrected, mc.preschedule = mc.preschedule, ...)
  return(ret) 
}

#' @rdname applyFolds
#' @export
# wrapper for function cv() of mboost, additional type "curves"
# create folds for data in long format
cvLong <- function(id, weights = rep(1, l=length(id)), 
                   type = c("bootstrap", "kfold", "subsampling", "curves"), 
                   B = ifelse(type == "kfold", 10, 25), prob = 0.5, strata = NULL){
  
  stopifnot(length(id) == length(weights))
  
  type <- match.arg(type)
  n <- length(weights)
  
  if(type == "curves"){
    if(!is.null(strata)) warning("Argument strata is ignored for type = 'curves'.")
    # set up folds so that always one curve is left out for the estimation
    folds <- - diag(length(unique(id))) + 1
    foldsLong <- folds[id, ] * weights    
    B <- length(unique(id))
  }else{
    # expand folds over the functional measures of the response
    folds <- cv(weights=rep(1, length(unique(id))), type = type, 
                B = B, prob = prob, strata = strata)
    foldsLong <- folds[id, , drop = FALSE] * weights
  }
  attr(foldsLong, "type") <- paste(B, "-fold ", type, sep = "")  
  return(foldsLong)
  
}

#' @rdname applyFolds
#' @export
# wrapper for function cv() of mboost, additional type "curves"
# add option id to sample on the level of id if there are repeated measures
cvMa <- function(ydim, weights = rep(1, l = ydim[1] * ydim[2]), 
                 type = c("bootstrap", "kfold", "subsampling", "curves"), 
                 B = ifelse(type == "kfold", 10, 25), prob = 0.5, strata = NULL, ...){
  
  dots <- list(...)
  
  if(any(names(dots) == "id")) message("argument id in cvMa is deprecated and ignored.")
  
  ncolY <- ydim[2]
  nrowY <- ydim[1]  
  
  type <- match.arg(type)
  n <- length(weights)
  
  if ( (nrowY * ncolY) != n) stop("The arguments weights and ydim do not match.")
  
  ## cvMa is only a wrapper for cvLong
  foldsMa <- cvLong(id = rep(1:nrowY, times = ncolY), weights = weights, 
                    type = type, B=B, prob = 0.5, strata = NULL) 
  return(foldsMa)
}


