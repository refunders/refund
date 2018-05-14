
#' Function to control estimation of smooth offset
#' 
#' @param k_min maximal number of k in s()
#' @param rule which rule to use in approx() of the response before calculating the 
#' global mean, rule=1 means no extrapolation, rule=2 means to extrapolate the 
#' closest non-missing value, see \code{\link[stats]{approx}}
#' @param silent print error messages of model fit?
#' @param cyclic defaults to FALSE, if TRUE cyclic splines are used
#' @param knots arguments knots passed to \code{\link[mgcv]{gam}}
#' 
#' @export 
o_control <- function(k_min=20, rule=2, silent=TRUE, cyclic=FALSE, knots=NULL) { 
  RET <- list(k_min=k_min, rule=rule, silent=silent, cyclic=cyclic, knots=knots)
  class(RET) <- c("offset_control")
  RET
}



#' Function to truncate time in functional data 
#' 
#' @param funVar names of functional variables that should be truncated
#' @param time name of time variable
#' @param newtime new time vector that should be used. Must be part of the old time-line.
#' @param data list containing all the data
#' @note All variables that are not part if \code{funVar}, or \code{time}
#' are simply copied into the new data list
#' @return A list with the data containing all variables of the original dataset
#' with the variables of \code{funVar} truncated according to \code{newtime}.
#' @examples
#' if(require(fda)){
#'   dat <- fda::growth
#'   dat$hgtm <- t(dat$hgtm[,1:10])
#'   dat$hgtf <- t(dat$hgtf[,1:10])
#'   
#'   ## only use time-points 1:16 of variable age
#'   datTr <- truncateTime(funVar=c("hgtm","hgtf"), time="age", newtime=1:16, data=dat)
#'   
#'   \dontrun{
#'   par(mfrow=c(1,2))
#'   with(dat, funplot(age, hgtm, main="Original data"))
#'   with(datTr, funplot(age, hgtm, main="Yearly data"))
#'   par(mfrow=c(1,1))   
#'   }
#' }
#' @export 
truncateTime <- function(funVar, time, newtime, data){
  
  stopifnot(c(funVar, time) %in% names(data))
  stopifnot(newtime %in% data[[time]])
  
  ret <- data
  ret[[time]] <- newtime
  
  for(i in 1:length(funVar)){
    ret[[funVar[i]]] <- ret[[funVar[i]]][ , data[[time]] %in% newtime] 
  }  
  rm(data)
  return(ret)  
}

#' Plot functional data with linear interpolation of missing values 
#' 
#' @param x optional, time-vector for plotting 
#' @param y matrix of functional data with functions in rows and measured times in columns; 
#' or vector or functional observations, in this case id has to be specified 
#' @param id defaults to NULL for y matrix, is id-variables for y in long format
#' @param rug logical. Should rugs be plotted? Defaults to TRUE.
#' @param ... further arguments passed to \code{\link[graphics]{matplot}}.
#' 
#' @details All observations are marked by a small cross (\code{pch=3}).
#' Missing values are imputed by linear interpolation. Parts that are
#' interpolated are plotted by dotted lines, parts with non-missing values as solid lines.
#' @examples 
#' \dontrun{
#' ### examples for regular data in wide format
#' data(viscosity)
#' with(viscosity, funplot(timeAll, visAll, pch=20))
#' if(require(fda)){
#'   with(fda::growth, funplot(age, t(hgtm)))
#' }
#' }
#' @export
funplot <- function(x, y, id=NULL, rug=TRUE, ...){
  
  ### Get further arguments passed to the matplot-functions
  dots <- list(...)
  
  getArguments <- function(x, dots=dots){
    if(any(names(dots) %in% names(x))){
      dots[names(dots) %in% names(x)]
    }else list()
  }
  
  plotWithArgs <- function(plotFun, args, myargs){        
    args <- c(myargs[!names(myargs) %in% names(args)], args)        
    do.call(plotFun, args)            
  }
  
  argsMatplot  <- getArguments(x=c(formals(graphics::matplot), par(), main="", sub=""), dots=dots)
  argsPlot <- getArguments(x=c(formals(graphics::plot.default), par()), dots=dots)
  
  if( !is.null(dim(y)) ){ # is.null(id)
    
    # Deal with missing values: interpolate data
    if (missing(x)) {
      if (missing(y)) 
        stop("must specify at least one of 'x' and 'y'")
      else{
        x <- seq_len(NCOL(y))
        xlabel <- "index"
      } 
    } else xlabel <- deparse(substitute(x))
    
    ylabel <- if (!missing(y)) 
      deparse(substitute(y))
    
    time <- x
    
    stopifnot(length(x) == ncol(y))
    
    # Check whether there are at least two values per row for the interpolation
    atLeast2values <- apply(y, 1, function(x) sum(is.na(x)) < length(x)-1 )
    if(any(!atLeast2values)) warning(sum(!atLeast2values), " rows contain less than 2 non-missing values.")  
    
    # Linear interpolation per row
    yint <- matrix(NA, ncol=ncol(y), nrow=nrow(y))
    yint[atLeast2values, ] <- t(apply(y[atLeast2values, ], 1, function(x) approx(time, x, xout=time)$y))
    
    # Plot the observed points
    plotWithArgs(matplot, args=argsMatplot, 
                 myargs=list(x=time, y=t(y), xlab=xlabel, ylab=ylabel, type="p", pch=3) )
    
    # Plot solid lines for parts of the function without missing values
    plotWithArgs(matplot, args=argsMatplot, 
                 myargs=list(x=time, y=t(y), type="l", lty=1, add=TRUE) )
    
    # Plot dotted lines for parts of the function without missing values
    plotWithArgs(matplot, args=argsMatplot, 
                 myargs=list(x=time, y=t(yint), type="l", lty=3, add=TRUE) )
    
    if(rug) rug(time, 0.01)
    
  }else{
    
    stopifnot(length(x)==length(y) & length(y)==length(id))
    
    idOrig <- id
    for(i in 1:length(unique(idOrig))){
      id[idOrig==unique(idOrig)[i]] <- i
    }
    
    xlabel <- deparse(substitute(x))
    ylabel <- deparse(substitute(y))
    
    # get color specification
    if("col" %in% names(dots)){
      col <- dots$col
      argsPlot$col <- 1
    }else{
      col <- 1
    }
    
    # expand vector col if it contains one color per trajectory
    if(length(col) == length(unique(id)) ){
      col <- rep(col, table(id) )
    }
    
    # there should be no mising values in long format
    temp <- data.frame(id, y, x, col, stringsAsFactors=FALSE) # dim(temp)
    temp <- na.omit(temp)   
    # order values of temp
    temp <- temp[order(temp$id, temp$x),] 
    id <- temp$id
    x <- temp$x
    y <- temp$y
    col <- temp$col
    rm(temp)
    
    # generate vector of colors for each observation
    if(is.null(argsPlot$col)){
      col <- rep( (unique(id)-1), table(id))
      col <- col %% 6 + 1  # only use the colors 1-6
    }
    argsPlot$col <- NULL
    
    # Plot the observed points
    if(!"add" %in% names(dots)){
      if(is.null(argsPlot$ylim)) argsPlot$ylim <- range(y, na.rm=TRUE) 
      plotWithArgs(plot, args=argsPlot, 
                   myargs=list(x=x[id==1], y=y[id==1], xlab=xlabel, ylab=ylabel, type="p", pch=3,
                               ylim=range(y, na.rm=TRUE), xlim=range(x, na.rm=TRUE),  col=col[id==1] ) )
    }
    
    for(i in unique(id)){
      plotWithArgs(points, args=argsPlot, 
                   myargs=list(x=x[id==i], y=y[id==i], xlab=xlabel, ylab=ylabel, type="p", pch=3,
                               col=col[id==i]) )
      plotWithArgs(lines, args=argsPlot, 
                   myargs=list(x=x[id==i], y=y[id==i], xlab=xlabel, ylab=ylabel, col=col[id==i]) )
    }
    
    if(rug) rug(x, 0.01)
    
    
  }
  
}



#####################################################################################

#' @rdname plot.FDboost
#' @export
#' 
### function to plot the observed response and the predicted values of a model
plotPredicted <- function(x, subset=NULL, posLegend="topleft", lwdObs=1, lwdPred=1, ...){
  
  stopifnot("FDboost" %in% class(x))
  
  if(any(class(x) == "FDboostScalar")){
    
    if(is.null(subset)) subset <- 1:length(x$response)
    response <- x$response[subset, drop=FALSE] 
    pred <- fitted(x)[subset, drop=FALSE]
    pred[is.na(response)] <- NA
    
  }else{
    
    if(!any(class(x) == "FDboostLong")){
      if(is.null(subset)) subset <- 1:x$ydim[1]
      response <- matrix(x$response, nrow=x$ydim[1], ncol=x$ydim[2])[subset, , drop=FALSE] 
      pred <- fitted(x)[subset, , drop=FALSE]
      pred[is.na(response)] <- NA
      yind <- x$yind
      id <- NULL
    }else{
      if(is.null(subset)) subset <- unique(x$id)
      response <- x$response[x$id %in% subset] 
      pred <- fitted(x)[x$id %in% subset]
      pred[is.na(response)] <- NA
      yind <- x$yind[x$id %in% subset] 
      id <- x$id[x$id %in% subset]
    }
    
  }
  
  if(is.character(response) | is.factor(x$response)){

    if(length(x$yind) > 1){
      message("For functional response that is not continuous only the predicted values are plotted.")
      funplot(yind, pred, id=id, pch=2, lwd=lwdPred, ... )
    }else{
      plot(response, pred, ylab="predicted", xlab="observed", ...)
    }
    
    
  } else{
    
    ylim <- range(response, pred, na.rm = TRUE)
    
    if(length(x$yind) > 1){
      # Observed values
      funplot(yind, response, id=id, pch=1, ylim=ylim, lty=3, 
              ylab=x$yname, xlab=attr(x$yind, "nameyind"), lwd=lwdObs, ...)
      funplot(yind, pred, id=id, pch=2, lwd=lwdPred, add=TRUE, ...)
      # predicted values
      legend(posLegend, legend=c("observed","predicted"), col=1, pch=1:2)  
    }else{
      plot(response, pred, ylab="predicted", xlab="observed", ...)
      abline(0,1)
    }
    
  }
  
}


#####################################################################################

#' @rdname plot.FDboost
#' @export
#' 
### function to plot the residuals
plotResiduals <- function(x, subset=NULL, posLegend="topleft", ...){
  
  stopifnot("FDboost" %in% class(x))
  
  if(any(class(x) == "FDboostScalar")){
    
    if(is.null(subset)) subset <- 1:length(x$response)
    response <- x$response[subset, drop=FALSE] 
    resid <- x$resid()[subset, drop=FALSE]
    
  }else{
    
    if(!any(class(x) == "FDboostLong")){ ## wide format
      if(is.null(subset)) subset <- 1:x$ydim[1]
      resid <- matrix(x$resid(), nrow = x$ydim[1])[subset, , drop=FALSE]
      yind <- x$yind
      id <- NULL
      
    }else{ ## long format
      if(is.null(subset)) subset <- unique(x$id)
      resid <- x$resid()[x$id %in% subset, drop=FALSE]
      yind <- x$yind[x$id %in% subset]
      id <- x$id[x$id %in% subset] 
    }
    
  }
  
  # Observed - predicted values
  if(length(x$yind) > 1){
    funplot(yind, resid, id=id, ylab=x$yname, xlab=attr(x$yind, "nameyind"), ...) 
  }else{
    plot(response, resid, ylab="residuals", xlab="observed", ...)
  }
  
}


#####################################################################################
### Goodness of fit

# function to get y, yhat and time
getYYhatTime <- function(object, breaks=object$yind){
  
  y <- matrix(object$response, nrow=object$ydim[1], ncol=object$ydim[2]) 
  time <- object$yind
  
  ### use the original time variable
  if(all(breaks==object$yind)){ 
    yhat <- matrix(object$fitted(), nrow=object$ydim[1], ncol=object$ydim[2]) 
  }else{ ### use a time variables according to breaks
    if(length(breaks)==1){ # length of equidistant time-points
      time <- seq( min(object$yind), max(object$yind), l=breaks)      
    }else{ # time-points ot be evaluated
      time <- breaks
    }
    # Interpolate observed values
    yInter <- t(apply(y, 1, function(x) approx(object$yind, x, xout=time)$y))
    # Get dataframe to predict values at time
    newdata <- list()
    for(j in 1:length(object$baselearner)){
      datVarj <- object$baselearner[[j]]$get_data()
      if(grepl("bconcurrent", names(object$baselearner)[j])){
        datVarj <- t(apply(datVarj[[1]], 1, function(x) approx(object$yind, x, xout=time)$y))
        datVarj <- list(datVarj)
      } 
      names(datVarj) <- names(object$baselearner[[j]]$get_data())
      newdata <- c(newdata, datVarj)
    }
    newdata[[attr(object$yind, "nameyind")]] <- time 
    yhatInter <- predict(object, newdata=newdata)
    
    y <- yInter
    yhat <- yhatInter
  }
  
  return(list(y=y, yhat=yhat, time=time))
}



#' Functional R-squared
#' 
#' Calculates the functional R-squared for a fitted FDboost-object
#' 
#' @param object fitted FDboost-object
#' @param overTime per default the functional R-squared is calculated over time
#' if \code{overTime=FALSE}, the R-squared is calculated per curve
#' @param breaks an optional vector or number giving the time-points at which the model is evaluated.
#' Can be specified as number of equidistant time-points or as vector of time-points.
#' Defaults to the index of the response in the model.
#' @param global logical. defaults to \code{FALSE}, 
#' if TRUE the global R-squared like in a normal linear model is calculated
#' @param ... currently not used
#' 
#' @note \code{breaks} cannot be changed in the case the \code{bsignal()} 
#' is used over the same domain
#' as the response! In that case you would have to rename the index of the response or that 
#' of the covariates.
#' 
#' @details \code{breaks} should be set to some grid, if there are many
#' missing values or time-points with very few observations in the dataset.
#' Otherwise at these points of t the variance will be almost 0 
#' (or even 0 if there is only one observation at a time-point),
#' and then the prediction by the local means \eqn{\mu(t)} is locally very good.
#' The observations are interpolated linearly if necessary.
#' 
#' Formula to calculate R-squared over time, \code{overTime=TRUE}: \cr
#' \eqn{R^2(t) = 1 - \sum_{i}( Y_i(t) - \hat{Y}_i(t))^2 /  \sum_{i}( Y_i(t) - \bar{Y}(t) )^2 } 
#' 
#' Formula to calculate R-squared over subjects, \code{overTime=FALSE}: \cr
#' \eqn{R^2_i = 1 - \int (Y_i(t) - \hat{Y}_i(t))^2 dt /  \int (Y_i(t) - \bar{Y}_i )^2 dt }
#' 
#' @references Ramsay, J., Silverman, B. (2006). Functional data analysis. 
#' Wiley Online Library. chapter 16.3
#' 
#' @return Returns a vector with the calculated R-squared and some extra information in attributes.
#' 
#' @export
funRsquared <- function(object, overTime=TRUE, breaks=object$yind, global=FALSE, ...){
  
  if(length(object$yind)<2 | any(class(object)=="FDboostLong")){
    y <- object$response
    yhat <- object$fitted()
    time <- object$yind
    id <- object$id
    if(is.null(id)) id <- 1:length(y)
    if(overTime & !global) {
      overTime <- FALSE
      message("For scalar or irregualr response the functional R-squared cannot be computed over time.")
    }
    if(length(object$yind)<2) global <- TRUE
  }else{
    # Get y, yhat and time of the model fit
    temp <- getYYhatTime(object=object, breaks=breaks)
    y <- temp$y
    yhat <- temp$yhat
    time <- temp$time
    
    stopifnot(dim(y)==dim(yhat))
    stopifnot(dim(y)[2]==length(time))
    
  }
  
  if(global){
    ret <- 1 - ( sum((y-yhat)^2, na.rm=TRUE)  / sum( (y-mean(y, na.rm=TRUE))^2, na.rm=TRUE) )
    attr(ret, "name") <- "global R-squared"
    return(ret)
  }
  
  ### for each time-point t 
  if(overTime){ 
    # Mean function over time (matrix containing the mean in each t in the whole column)
    mut <- matrix(colMeans(y, na.rm=TRUE), nrow=nrow(y), ncol=ncol(y), byrow=TRUE)
    
    # numerator cannot be 0
    num <- colSums((y - mut)^2, na.rm=TRUE)
    num[round(num, 2)==0] <- NA
    ret <- 1 - ( colSums((y - yhat)^2, na.rm=TRUE) / num ) # over t 
    
    # Set values of R-squared to NA if too many values are missing
    ret[apply(y, 2, function(x) sum(is.na(x))>0.5*length(x) )] <- NA
    
    attr(ret, "name") <- "R-squared over time"
    attr(ret, "time") <- time
    attr(ret, "missings") <- apply(y, 2, function(x) sum(is.na(x))/length(x) )
    
  }else{ ### for each subject i
    if(length(object$yind)<2 | any(class(object)=="FDboostLong")){
      # Mean for each subject
      mut <- tapply(y, id, mean, na.rm=TRUE  )[id]
      # numerator cannot be 0
      num <- tapply((y - mut)^2, id, mean, na.rm=TRUE )[id]
      num[round(num, 2)==0] <- NA 
      ret <- 1 - tapply((y - yhat)^2 /num, id, mean, na.rm=TRUE  )
      attr(ret, "name") <- "MSE over subjects"              
    }else{
      # Mean for each subject
      mut <- matrix(rowMeans(y, na.rm=TRUE), nrow=nrow(y), ncol=ncol(y), byrow=TRUE)
      # numerator cannot be 0
      num <- rowSums((y - mut)^2, na.rm=TRUE)
      num[round(num, 2)==0] <- NA    
      ret <- 1 - ( rowSums((y - yhat)^2, na.rm=TRUE) /  num ) # over subjects
      
      attr(ret, "name") <- "R-squared over subjects"
      attr(ret, "missings") <- apply(y, 1, function(x) sum(is.na(x))/length(x)) 
    }
  }
  return(ret)
}


# R2t <- funRsqrt(mod)
# plot(attr(R2t, "time"), R2t, ylim=c(0,1), type="b")
# points(attr(R2t, "time"), attr(R2t, "missings"), col=2, type="b")
# abline(h=c(0, 1))
# 
# R2i <- funRsqrt(mod, overTime=FALSE)
# plot(R2i, type="b", ylim=c(0,1))
# points(attr(R2i, "missings"), col=2, type="b")
# abline(h=c(0, 1))

#' Functional MSE
#' 
#' Calculates the functional MSE for a fitted FDboost-object
#' 
#' @param object fitted FDboost-object
#' @param overTime per default the functional R-squared is calculated over time
#' if \code{overTime=FALSE}, the R-squared is calculated per curve
#' @param breaks an optional vector or number giving the time-points at which the model is evaluated.
#' Can be specified as number of equidistant time-points or as vector of time-points.
#' Defaults to the index of the response in the model.
#' @param global logical. defaults to \code{FALSE}, 
#' if TRUE the global R-squared like in a normal linear model is calculated
#' @param relative logical. defaults to \code{FALSE}. If \code{TRUE} the MSE is standardized
#' by the global variance of the response \cr
#' \eqn{ n^{-1} \int  \sum_i (Y_i(t) - \bar{Y})^2 dt \approx  G^{-1} n^{-1} \sum_g \sum_i (Y_i(t_g) - \bar{Y})^2 } 
#' @param root take the square root of the MSE
#' @param ... currently not used
#' 
#' @note \code{breaks} cannot be changed in the case the \code{bsignal()} 
#' is used over the same domain
#' as the response! In that case you would have to rename the index of the response or that 
#' of the covariates.
#' 
#' @details 
#' Formula to calculate MSE over time, \code{overTime=TRUE}: \cr
#' \eqn{ MSE(t) = n^{-1} \sum_i (Y_i(t) - \hat{Y}_i(t))^2 } 
#' 
#' Formula to calculate MSE over subjects, \code{overTime=FALSE}: \cr
#' \eqn{ MSE_i = \int (Y_i(t) - \hat{Y}_i(t))^2 dt  \approx G^{-1} \sum_g (Y_i(t_g) - \hat{Y}_i(t_g))^2}
#' 
#' @return Returns a vector with the calculated MSE and some extra information in attributes.
#' 
#' @export
funMSE <- function(object, overTime=TRUE, breaks=object$yind, global=FALSE, 
                   relative=FALSE, root=FALSE, ...){
  
  if(length(object$yind)<2 | any(class(object)=="FDboostLong")){
    y <- object$response
    yhat <- object$fitted()
    time <- object$yind
    id <- object$id
    if(is.null(id)) id <- 1:length(y)
    if(overTime & !global) {
      overTime <- FALSE
      message("For scalar or irregualr response the functional MSE cannot be computed over time.")
    }
  }else{
    # Get y, yhat and time of the model fit
    temp <- getYYhatTime(object=object, breaks=breaks)
    y <- temp$y
    yhat <- temp$yhat
    time <- temp$time
  }
  
  if(global){
    ret <- mean((y-yhat)^2, na.rm=TRUE)
    attr(ret, "name") <- "global MSE"
  }else{
    ### for each time-point t 
    if(overTime){     
      ret <- colMeans((y - yhat)^2, na.rm=TRUE)   
      attr(ret, "name") <- "MSE over time"
      attr(ret, "time") <- time
      attr(ret, "missings") <- apply(y, 2, function(x) sum(is.na(x))/length(x))     
    }else{ 
      ### for each subject i
      if(length(object$yind)<2 | any(class(object)=="FDboostLong")){
        ret <- tapply((y - yhat)^2, id, mean, na.rm=TRUE  )
        attr(ret, "name") <- "MSE over subjects"              
      }else{
        ret <- rowMeans((y - yhat)^2, na.rm=TRUE)
        attr(ret, "name") <- "MSE over subjects"      
        attr(ret, "missings") <- apply(y, 1, function(x) sum(is.na(x))/length(x))
      }
    }    
  }
  
  if(relative){
    variY <- mean( (y - mean(y, na.rm=TRUE))^2, na.rm=TRUE ) # global variance of y
    ret <- ret / variY
    attr(ret, "name") <- paste("relative", attr(ret, "name"))
  }
  
  if(root){
    ret <- sqrt(ret)
    attr(ret, "name") <- paste("root", attr(ret, "name"))
  }
  
  return(ret)
}


#' Functional MRD
#' 
#' Calculates the functional MRD for a fitted FDboost-object
#' 
#' @param object fitted FDboost-object with regular response
#' @param overTime per default the functional MRD is calculated over time
#' if \code{overTime=FALSE}, the MRD is calculated per curve
#' @param breaks an optional vector or number giving the time-points at which the model is evaluated.
#' Can be specified as number of equidistant time-points or as vector of time-points.
#' Defaults to the index of the response in the model.
#' @param global logical. defaults to \code{FALSE}, 
#' if TRUE the global MRD like in a normal linear model is calculated
#' @param ... currently not used
#' 
#' @note \code{breaks} cannot be changed in the case the \code{bsignal()} 
#' is used over the same domain
#' as the response! In that case you would have to rename the index of the response or that 
#' of the covariates.
#' 
#' @details 
#' Formula to calculate MRD over time, \code{overTime=TRUE}: \cr
#' \eqn{ MRD(t) = n^{-1} \sum_i |Y_i(t) - \hat{Y}_i(t)| / |Y_i(t)| } 
#' 
#' Formula to calculate MRD over subjects, \code{overTime=FALSE}: \cr
#' \eqn{ MRD_{i} = \int |Y_i(t) - \hat{Y}_i(t)| / |Y_i(t)| dt  \approx G^{-1} \sum_g |Y_i(t_g) - \hat{Y}_i(t_g)| / |Y_i(t)|}
#' 
#' @return Returns a vector with the calculated MRD and some extra information in attributes.
#' 
#' @export
funMRD <- function(object, overTime=TRUE, breaks=object$yind, global=FALSE,  ...){
  
  if(length(object$yind)<2 | any(class(object)=="FDboostLong")){
    y <- object$response
    yhat <- object$fitted()
    time <- object$yind
    id <- object$id
    if(is.null(id)) id <- 1:length(y)
    if(overTime & !global) {
      overTime <- FALSE
      message("For scalar or irregualr response the functional MRD cannot be computed over time.")
    }
  }else{
    # Get y, yhat and time of the model fit
    temp <- getYYhatTime(object=object, breaks=breaks)
    y <- temp$y
    yhat <- temp$yhat
    time <- temp$time
  }
  
  # You cannot use observations that are 0, so set them to NA
  y1 <- y
  y1[ round(y1, 1) == 0 ] <- NA
  
  if(global){
    ret <- mean( abs((y1 - yhat) / y1), na.rm=TRUE  )
    attr(ret, "name") <- "global MRD"
  }else{
    ### for each time-point t 
    if(overTime){     
      ret <- colMeans( abs((y1 - yhat) / y1), na.rm=TRUE  )   
      attr(ret, "name") <- "MRD over time"
      attr(ret, "time") <- time
      attr(ret, "missings") <- apply(y, 2, function(x) sum(is.na(x))/length(x))     
    }else{ 
      ### for each subject i
      if(length(object$yind)<2 | any(class(object)=="FDboostLong")){
        ret <- tapply( abs((y1 - yhat) / y1), id, mean, na.rm=TRUE  )
        attr(ret, "name") <- "MRD over subjects"              
      }else{
        ret <- rowMeans( abs((y1 - yhat) / y1), na.rm=TRUE  )
        attr(ret, "name") <- "MRD over subjects"      
        attr(ret, "missings") <- apply(y, 1, function(x) sum(is.na(x))/length(x))
      }
    }    
  }
  
  return(ret)
}


###############################################################################
## helper functions for identifiability 


## based on code of function ff() in package refund
## see Scheipl and Greven, 2016: Identifiability in penalized function-on-function regression models 
# X1 matrix of functional covariate x(s)
# L matrix of integration weights
# Bs matrix of spline expansion in s
# K penalty matrix
# xname name of functional covariate
# penalty the type of the penalty one of "ps" or "pps"
# cumOverlap DEPRECATED should a cumulative overlap be computed, 
#   which is especially suited for a historical effect with triangular coefficient surface?
# limits the limits function of the historical effect, default to NULL for unconstrained effect
#   if limits is supplied a cumulative overlap and condition number according to limits is computed 
# yind, id, X1des, ind0, xind passed from X_hist() to compute sequential identifiability measures
# giveWarnings should warnings be printed
check_ident <- function(X1, L, Bs, K, xname, penalty, 
                        ## cumOverlap = FALSE, 
                        limits = NULL, yind = NULL, 
                        t_unique = NULL, 
                        id = NULL, 
                        X1des = NULL, ind0 = NULL, xind = NULL, 
                        giveWarnings = TRUE){
  
  ## center X1 per column
  X1 <- scale(X1, scale = FALSE)
  
  #print("check.ident")
  ## check whether (number of basis functions in Bs) < (number of relevant eigenfunctions of X1)
  evls <- svd(X1, nu = 0, nv = 0)$d^2 # eigenvalues of centered fun. cov.
  evls[evls<0] <- 0
  maxK <- max(1, min(which((cumsum(evls)/sum(evls)) >= .995)))
  bsdim <- ncol(Bs) # number of basis functions in Bs
  if(maxK <= 4)
    warning("Very low effective rank of <", xname,
            "> detected. ", maxK,
            " largest eigenvalues of its covariance alone account for >99.5% of ",
            "variability. <bfpc> might be a better choice here.")
  if(maxK < bsdim){
    warning("<k> (",  bsdim,") larger than effective rank of <", xname,
            "> (", maxK, "). ",
            "Model identifiable only through penalty.")
  }
  
  ### compute condition number of Ds^t Ds
  Ds <- (X1 * L) %*% Bs
  ## DstDs <- crossprod(Ds)
  ## e_DstDs <- try(eigen(DstDs))
  ## e_DstDs$values <- pmax(0, e_DstDs$values) # set negative eigenvalues to 0
  ## logCondDs <- log10(e_DstDs$values[1]) - log10(tail(e_DstDs$values, 1))
  evDs <- svd(Ds, nu = 0, nv = 0)$d^2  ## the same as eigenvalues of DstDs
  logCondDs <- log10(max(evDs)) - log10(min(evDs))
  if(giveWarnings & logCondDs > 6 & is.null(limits)){
    warning("Condition number for <", xname, "> greater than 10^6 (logCondDs = ", round(logCondDs, 2),"). ", 
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
    if( length(yind) < nrow(X1des) | all(table(yind) / max(id) > 0.8) ){
      if(is.null(t_unique)) t_unique <- sort(unique(yind))
      logCondDs_hist <- rep(NA, length=length(t_unique))
      for(k in 1:length(t_unique)){
        Ds_t <- X1des[yind==t_unique[k], ] # get rows of Ds corresponding to yind
        ind0Bs_t <- ind0Bs[yind==t_unique[k], ] # get rows of ind0Bs corresponding to yind
        # only keep columns that are not completely 0, otherwise matrix is always rank deficient
        # idea: only this part is used to model y(t) at this point
        # also delete if not perfectly but almost zero, for all spline bases
        Ds_t <- Ds_t[ , apply(ind0Bs_t, 2, function(x) !all(abs(x)<10^-1) ), drop = FALSE ]
        
        if(dim(Ds_t)[2]!=0){ # for matrix with 0 columns does not make sense
          ## DstDs_t <- crossprod(Ds_t)
          ## e_DstDs_t <- try(eigen(DstDs_t))
          ## e_DstDs_t$values <- pmax(0, e_DstDs_t$values) # set negative eigenvalues to 0
          ## logCondDs_t <- log10(e_DstDs_t$values[1]) - log10(tail(e_DstDs_t$values, 1))
          ## logCondDs_hist[k] <- logCondDs_t
          evDs <- svd(Ds_t, nu = 0, nv = 0)$d^2  ## the same as eigenvalues of DstDs
          logCondDs_hist[k] <- log10(max(evDs)) - log10(min(evDs))
        }
        ## matplot(xind, Bs, type="l", lwd=2, ylim=c(-2,2)); rug(xind); rug(yind, col=2, lwd=2)
        ## matplot(knots[1:ncol(Ds_t)], t(Ds_t), type="l", lwd=1, add=TRUE)
        ## lines(t_unique, logCondDs_hist-6, col=2, lwd=4)
      }
      names(logCondDs_hist) <- round(t_unique,2)
      
      ### implementation for seriously irregular observation points t
    }else{
      ## use the mean number of grid points, in the case of irregular t
      # t_unique <- seq(min(yind), max(yind), length=round(mean(table(id))))
      ### use quantiles of yind, as only at places with observations effect can be identifiable
      ### using quantiles prevents Ds_t from beeing completely empty
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
        Ds_t <- Ds_t[ , apply(ind0Bs_t, 2, function(x) !all(abs(x)<10^-1) ), drop = FALSE] 
        
        if(dim(Ds_t)[2]!=0){ # for matrix with 0 columns does not make sense
          ## DstDs_t <- crossprod(Ds_t)
          ## e_DstDs_t <- try(eigen(DstDs_t))
          ## e_DstDs_t$values <- pmax(0, e_DstDs_t$values) # set negative eigenvalues to 0
          ## logCondDs_t <- log10(e_DstDs_t$values[1]) - log10(tail(e_DstDs_t$values, 1))
          ## logCondDs_hist[k] <- logCondDs_t
          evDs <- svd(Ds_t, nu = 0, nv = 0)$d^2  ## the same as eigenvalues of DstDs
          logCondDs_hist[k] <- log10(max(evDs)) - log10(min(evDs))
        }
      }
      names(logCondDs_hist) <- round(t_unique[-length(t_unique)],2)
    }
    if(giveWarnings & any(logCondDs_hist > 6)){
      # get the first and the last entry of t, for which the condition number is >10^6
      tempL <- names(which.min(which(logCondDs_hist > 6)))
      tempU <- names(which.max(which(logCondDs_hist > 6)))
      warning("Condition number for <", xname, "> considering limits of historical effect ", 
              "greater than 10^6, for time-points between ", tempL, " and ", tempU, ". ",
              "Effect in this region identifiable only through penalty.")
    }
  } ## end of computation of logCondDs_hist for historical effects
  
  getOverlap <- function(subset, X1, L, Bs, K){
    # In the case that all observations are 0, kernel is everything -> kernel overlap
    if(all(X1[ , subset]==0)){
      return(5)
    }
    N.X <- Null(t(X1[ , subset, drop=FALSE]))
    if(NCOL(N.X) != 0){
      N.pen <- diag(L[1, subset]) %*% Bs[subset, , drop=FALSE] %*% Null(K)
    }else N.pen <- 0
    
    if (any(c(NCOL(N.X) == 0, NCOL(N.pen) == 0))) {
      nullOverlap <- 0
    }
    else {
      nullOverlap <- trace_lv(svd(N.X)$u, svd(N.pen)$u)
    }
    return(nullOverlap)
  }
  
  cumOverlapKe <- NULL
  overlapKe <- NULL
  overlapKeComplete <- NULL

  ## sequential overlap for historical model with general integration limits
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
  
  # look at overlap with whole functional covariate 
  overlapKeComplete  <- getOverlap(subset=1:ncol(X1), X1=X1, L=L, Bs=Bs, K=K)
  
  if(giveWarnings & overlapKe >= 1){
    warning("Kernel overlap for <", xname, "> and the specified basis and penalty detected. ",
            "Changing basis for x-direction to <penalty='pss'> to make model identifiable through penalty. ", 
            "Coefficient surface estimate will be inherently unreliable. ", 
            "See Scheipl & Greven (2016) for details and alternatives.") 
    penalty <- "pss"
  }
  
  return(list(logCondDs=logCondDs, logCondDs_hist=logCondDs_hist,  
              overlapKe=overlapKe, cumOverlapKe=cumOverlapKe, 
              overlapKeComplete=overlapKeComplete, 
              maxK=maxK, penalty=penalty))
}


## measure degree of overlap between the spans of X and Y using A=svd(X)$u, B=svd(Y)$u
## code written by Fabian Scheipl
trace_lv <- function(A, B, tol=1e-10){
  ## A, B orthnormal!!
  
  # Rolf Larsson, Mattias Villani (2001)
  # "A distance measure between cointegration spaces"
  
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


## use a penalty matrix with full rank, so-called "shrinkage approach" 
## after Marra and Wood (2011) Practical variable selection for generalized additive models. 
## code taken from smooth.construct.pss.smooth.spec() in package refund (written by Fabian Scheipl)
## which is based on mgcv of Simon Wood, e.g. smooth.construct.tp.smooth.spec() 
## function with the following arguments 
# K sqaured differences penalty matrix
# difference degree of difference
# shrink shrinkage parameter 
penalty_pss <- function(K, difference, shrink){
  
  stopifnot(shrink > 0, shrink < 1)
  
  bsdim <- nrow(K) # ncol(Bs) # number of basis functions in Bs
  ## add shrinkage term to penalty: 
  ## Modify the penalty by increasing the penalty 
  ## on the unpenalized space from zero...
  es <- eigen(K, symmetric=TRUE)
  ## now add a penalty on the penalty null space
  es$values[(bsdim-difference+1):bsdim] <- es$values[bsdim-difference]*shrink
  ## ... so penalty on null space is still less than that on range space.
  K <- es$vectors %*% (as.numeric(es$values)*t(es$vectors)) 
  
  return(K)
}


###################### expand.call() is taken from R package refund 0.1-15
# Return call with all possible arguments
#
# Return a call in which all of the arguments which were supplied or have presets are specified by their 
# full names and their supplied or default values.
#
# @param definition a function. See \code{\link[base]{match.call}}.
# @param call an unevaluated call to the function specified by definition. See \code{\link[base]{match.call}}.
# @param expand.dots logical. Should arguments matching ... in the call be included or 
# left as a ... argument? See \code{\link[base]{match.call}}.
# @return An object of mode "\code{\link[base]{call}}".
# @author Fabian Scheipl
# @seealso \code{\link[base]{match.call}}
expand.call <- function(definition=NULL, call=sys.call(sys.parent(1)), expand.dots = TRUE){
  call <- match.call(definition, call, expand.dots)
  #given args:
  ans <- as.list(call)
  
  #possible args:
  frmls <- formals(safeDeparse(ans[[1]]))
  #remove formal args with no presets:
  frmls <- frmls[!sapply(frmls, is.symbol)]
  
  add <- which(!(names(frmls) %in% names(ans)))
  return(as.call(c(ans, frmls[add])))
}

safeDeparse <- function(expr){
  # turn an expression into a _single_ string, regardless of the expression's length
  ret <- paste(deparse(expr), collapse="")
  #rm whitespace
  gsub("[[:space:]][[:space:]]+", " ", ret)
}



#' Function to Reweight Data
#' 
#' @param data a named list or data.frame.
#' @param argvals character (vector); name(s) for entries in data giving 
#' the index for observed grid points; must be supplied if \code{vars} is not supplied.
#' @param vars character (vector); name(s) for entries in data, which
#' are subsetted according to weights or index. Must be supplied if \code{argvals} is not supplied.
#' @param longvars variables in long format, e.g., a response that is observed at curve specific grids. 
#' @param weights vector of weights for observations. Must be supplied if \code{index} is not supplied.
#' @param index vector of indices for observations. Must be supplied if \code{weights} is not supplied.
#' @param idvars character (vector); index, which is needed to expand \code{vars} to be conform
#' with the \code{hmatrix} structure when using \code{bhistx}-base-learners or to be conform with 
#' variables in long format specified in \code{longvars}. 
#' @param compress logical; whether \code{hmatrix} objects are saved in compressed form or not. Default is \code{TRUE}.
#' Should be set to \code{FALSE} when using \code{reweightData} for nested resampling.
#' 
#' @return A list with the reweighted or subsetted data.
#' 
#' @details \code{reweightData} indexes the rows of matrices and / or positions of vectors by using
#' either the \code{index} or the \code{weights}-argument. To prevent the function from indexing
#' the list entry / entries, which serve as time index for observed grid points of each trajectory of
#' functional observations, the \code{argvals} argument (vector of character names for these list entries) 
#' can be supplied. If \code{argvals} is not supplied, \code{vars} must be supplied and it is assumed that 
#' \code{argvals} is equal to \code{names(data)[!names(data) \%in\% vars]}.
#' 
#' When using \code{weights}, a weight vector of length N must be supplied, where N is the number of observations.
#' When using \code{index}, the vector must contain the index of each row as many times as it shall be included in the
#' new data set.
#' 
#' @examples
#' ## load data
#' data("viscosity", package = "FDboost")
#' interval <- "101"
#' end <- which(viscosity$timeAll == as.numeric(interval))
#' viscosity$vis <- log(viscosity$visAll[ , 1:end])
#' viscosity$time <- viscosity$timeAll[1:end]
#' 
#' ## what does data look like
#' str(viscosity)
#' 
#' ## do some reweighting
#' # correct weights
#' str(reweightData(viscosity, vars=c("vis", "T_C", "T_A", "rspeed", "mflow"), 
#'     argvals = "time", weights = c(0, 32, 32, rep(0, 61))))
#' 
#' str(visNew <- reweightData(viscosity, vars=c("vis", "T_C", "T_A", "rspeed", "mflow"), 
#'     argvals = "time", weights = c(0, 32, 32, rep(0, 61))))
#' # check the result
#' # visNew$vis[1:5, 1:5] ## image(visNew$vis)
#' 
#' # incorrect weights
#' str(reweightData(viscosity, vars=c("vis", "T_C", "T_A", "rspeed", "mflow"), 
#'     argvals = "time", weights = sample(1:64, replace = TRUE)), 1)
#' 
#' # supply meaningful index
#' str(visNew <- reweightData(viscosity, vars = c("vis", "T_C", "T_A", "rspeed", "mflow"), 
#'               argvals = "time", index = rep(1:32, each = 2)))
#' # check the result
#' # visNew$vis[1:5, 1:5]
#' 
#' # errors
#' if(FALSE){
#'    reweightData(viscosity, argvals = "")
#'    reweightData(viscosity, argvals = "covThatDoesntExist", index = rep(1,64))
#'    }
#'    
#' @author David Ruegamer, Sarah Brockhaus
#' 
#' @export 
reweightData <- function(data, argvals, vars, 
                         longvars = NULL, weights, index, 
                         idvars = NULL, compress = FALSE)
{
  
  if(missing(argvals) & missing(vars)) 
    stop("Either argvals or vars must be supplied.")
  if(missing(weights) & missing(index)) 
    stop("Either weights or index must be supplied.")

  # get names of data
  nd <- names(data)
  
  # if(missing(idvars)) idvars <- NULL
  
  # drop not used entries if both argvals and vars are given
  if(!missing(argvals) & !missing(vars)){
    
    data[nd[!nd %in% c(argvals, vars, longvars, idvars)]] <- NULL
    nd <- names(data) # reset names
    
  }
  
  # define argvals or vars if missing exclusively
  if(missing(argvals)) argvals <- nd[!nd %in% c(vars, idvars, longvars)]
  if(missing(vars)) vars <- nd[!nd %in% c(argvals, idvars, longvars)]
  
  whichNot <- which(!c(argvals, vars, idvars, longvars) %in% nd)
  
  # check names
  if(length(whichNot) != 0) 
    stop(paste0("Could not find ", 
                paste(c(argvals, vars, idvars, longvars)[whichNot], collapse = ", "),
                " in data."))
  
  # check for hmatrix and delete in argvals or vars if present
  whichHmat <- sapply(data[vars], function(x) "hmatrix" %in% class(x))
  
  # get dimensions of data
  dimd <- lapply(data, dim)
  isVec <- sapply(dimd, is.null)
  
  if(length(vars) == 1 && sum(whichHmat) == 1){
    
    n <- nrow(attr(data[[vars]],"x"))
    if((nalt <- length(unique(getId.hmatrix(data[[vars[whichHmat]]]))))!=n)
      n <- nalt
      #warning("Dimension of hmatrix is not equal to its corresponding attribute.")
    
  }else{
      
      if(length(whichHmat) == 0){ 
        
        if(is.null(vars)){
          
          if(is.null(idvars)) stop("idvars must be given if vars is NULL.")
          n <- n_variables <- length(unique(data[[idvars[1]]]))
          
        }else{
          
          n <- length(data[[vars[[1]]]])
          n_variables <- sapply(data[vars], NROW)
          
        }
        
      }else{
      
          n <- NROW(data[[vars[!whichHmat][1]]])
          n_variables <- sapply(data[vars][!whichHmat], NROW)
      
      }
      
      if(any(n_variables[1] != n_variables)) stop("variables imply different number of rows.")
      
    }
  
  if(!missing(weights) && length(weights) != n) 
    stop("Length of weights and number of observations do not match!")
  # if(!missing(weights) && sum(weights) != n)
  #   warning("The resulting data will have more / less observations than the original data.")
  
  # transform weights in index and vice versa 
  if(missing(index)) index <- rep(1:n, weights) else index <- sort(index)
  ## computation of longvars needs weights 
  if(missing(weights)) weights <- sapply(1:n, function(i) sum(index == i)) 
  
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  # is.wholenumber(x) # is TRUE
  
  if(!missing(weights) && any(!is.wholenumber(weights))) stop("weights can only contain integers.")
  if(any(!is.wholenumber(index))) stop("index can only contain integers.")
  
  ## check that all idvars are equal
  if(length(idvars)>1)
    if(!all(sapply(data[idvars][-1],function(x)all.equal(data[idvars][[1]],x)=="TRUE")))
      stop("All idvars must be identical.")
  idvars_new <- NULL
  
  # get names of hmatrix variables
  nhm <- vars[whichHmat]
  # get all names of data items, which are neither a hmatrix 
  # nor supplied via the idvars or longvars argument
  remV <- !nd %in% c(nhm, idvars, longvars)  
  # remove those
  nd <- nd[remV] 
  # the same for isVec
  isVec <- isVec[remV]
  
  # if there is a hmatrix in the data, subsetting is done by reconstructing the hmatrix appropriately
  if(any(whichHmat)){
    
    # construct a list for new hmatrices
    newHmats <- vector("list", length(nhm))
    
    ## construct the new hmatrices
    for(j in 1:length(nhm)){
      
      ## check that idvars == idvars[[1]] and match id-variables in all hmatrix-objects
      if(!is.null(idvars) && !(all.equal(c(getId(data[[nhm[j]]])), c(data[[idvars[1]]])) == "TRUE")) 
        stop("id variable in hmatrix object must be equal to idvars")
      
      ## subset hmatrix
      newHmats[[j]] <- subset_hmatrix(data[[nhm[j]]], index = index, compress = compress)
      
      if( any(class(data[[nhm[j]]]) == "AsIs") ){
        newHmats[[j]] <- I(newHmats[[j]])
      }
      
    }
    names(newHmats) <- nhm
    
  }else{ ## if there are no hmatrices, set the list and corresponding names to NULL
    
    newHmats <- NULL
    nhm <- NULL
    
  }
  
  temp_long <- NULL 
  ## do the indexing for the variables in long format 
  if(!is.null(longvars)){
    if(any(idvars %in% longvars)) longvars <- longvars[!longvars %in% idvars]
    ## create weights and index in long format
    weights_long <- weights[data[[idvars[1]]]]
    index_long <- rep(1:length(weights_long), weights_long)
    ## indexing variables in long format
    temp_long <- lapply(longvars, function(nameWithoutDim) data[[nameWithoutDim]][index_long])
    
    #### create new id variable 1, 2, 3, ... that can be used for FDboost()
    #### always generate idvars_new, even though it exists already because of a hmatrix object 
    idvars_new_hmatrix <- NULL
    if(!is.null(idvars_new)){
      idvars_new_hmatrix <- idvars_new
    }
    
    temp_idvars <- data[[idvars[1]]][index_long] # compute new id-variable
    
    ## gives equal numbers to repetitions of the same observation
    ## idvars_new <- c(factor(temp_idvars))
   
    ## gives different numbers to repetitions of the same observation
    my_index_long <- index_long 
    my_temp_idvars <- temp_idvars
    i <- 1
    ## add 0.1^1 to duplicates, 0.1^1 + 0.1^2 = 0.11 to triplicates, ...
    while(any(duplicated(my_index_long))){ # loop until no more duplicates in the data  
      my_temp_idvars[duplicated(my_index_long)] <- my_temp_idvars[duplicated(my_index_long)] + 0.1^i
      my_index_long[duplicated(my_index_long)] <- my_index_long[duplicated(my_index_long)] + 0.1^i
      i <- i + 1
    }
    idvars_new <- c(factor(my_temp_idvars))
    ## check whether id variable of hmatrix-object and id variable of long variables are equal
    if(!is.null(idvars_new_hmatrix)){
      if(!all(idvars_new == idvars_new_hmatrix)) 
        warning("id variable generated for long variables and id variable of hmatrix-object do not match. ",
                "Sort the long variables and the hmatrix-object by the time variable.")
    }
  }
  
  ## compute idvars_new in hmatrix-part or in longvars part, but add to data here 
  ## if idvars exist, subset accordingly;
  ## idvars has to be the same for all hmatrix-objects and response! 
  if(!is.null(idvars)){
    
    ## only works for common observation grid of response
    ## idvars_new <- rep(1:length(index), nc) # index = c(1, 1, 2) -> 1, 2, 3
    
    for(ifr in idvars){
      data[[ifr]] <- idvars_new
    } 
    ## data[[idvars]] <- getId(newHmats[[j]])
    argvals <- c(argvals, idvars)
  }
  
  inAVs <- nd %in% argvals
  
  ## recycle data
  data <- c(lapply(nd[!isVec & !inAVs], function(nameWithDim) data[[nameWithDim]][index, , drop=FALSE]),
            lapply(nd[isVec & !inAVs], function(nameWithoutDim) data[[nameWithoutDim]][index]), 
            newHmats, 
            temp_long, 
            data[argvals])
  names(data) <- c(nd[!isVec & !inAVs], nd[isVec & !inAVs], nhm, longvars, argvals)
  
  return(drop(data))
  
}

