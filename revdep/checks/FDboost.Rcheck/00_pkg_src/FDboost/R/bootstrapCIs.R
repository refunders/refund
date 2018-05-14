#' Function to compute bootstrap confidence intervals
#' 
#' The model is fitted on bootstrapped samples of the data to compute bootstrapped 
#' coefficient estimates. To determine the optimal stopping iteration an inner bootstrap 
#' is run within each bootstrap fold. 
#' As estimation by boosting shrinks the coefficient estimates towards zero, 
#' to bootstrap confidence intervals are biased towards zero. 
#'
#' @param object a fitted model object of class \code{FDboost}, 
#' for which the confidence intervals should be computed.
#' @param which a subset of base-learners to take into account for 
#' computing confidence intervals. 
#' @param resampling_fun_outer function for the outer resampling procedure.
#' \code{resampling_fun_outer} must be a function with arguments \code{object}
#' and \code{fun}, where \code{object} corresponds to the fitted 
#' \code{FDboost} object and \code{fun} is passed to the \code{fun}
#' argument of the resampling function (see examples).
#' If \code{NULL}, \code{\link{applyFolds}} is used with 100-fold boostrap.
#' Further arguments to \code{\link{applyFolds}} can be passed via \code{...}.
#' Although the function can be defined very flexible, it is recommended 
#' to use \code{applyFolds} and, in particular, not \code{cvrisk}, 
#' as in this case, weights of the inner and outer 
#' fold will interact, probably causing the inner 
#' resampling to crash. For bootstrapped confidence intervals
#' the outer function should usually be a bootstrap type of resampling.
#' @param resampling_fun_inner function for the inner resampling procudure,
#' which determines the optimal stopping iteration in each fold of the
#' outer resampling procedure. Should be a function with one argument
#' \code{object} for the fitted \code{FDboost} object. 
#' If \code{NULL}, \code{cvrisk} is used with 25-fold bootstrap.
#' @param B_outer Number of resampling folds in the outer loop.
#' Argument is overwritten, when a custom \code{resampling_fun_outer}
#' is supplied.
#' @param B_inner Number of resampling folds in the inner loop.
#' Argument is overwritten, when a custom \code{resampling_fun_inner}
#' is supplied.
#' @param type_inner character argument for specifying the cross-validation method for
#' the inner resampling level. Default is \code{"bootstrap"}. Currently  
#' bootstrap, k-fold cross-validation and subsampling are implemented.
#' @param levels the confidence levels required. If NULL, the 
#' raw results are returned. 
#' @param ... further arguments passed to \code{\link{applyFolds}} if
#' the default for \code{resampling_fun_outer} is used
#'
#' @author David Ruegamer, Sarah Brockhaus
#' 
#' @note Note that parallelization can be achieved by defining
#' the \code{resampling_fun_outer} or \code{_inner} accordingly.
#' See, e.g., \code{\link{cvrisk}} on how to parallelize resampling
#' functions or the examples below. Also note that by defining
#' a custum inner or outer resampling function the respective
#' argument \code{B_inner} or \code{B_outer} is ignored.
#' For models with complex baselearners, e.g., created by combining
#' several baselearners with the Kronecker or row-wise tensor product,
#' it is also recommended to use \code{levels = NULL} in order to
#' let the function return the raw results and then manually compute
#' confidence intervals.
#' If a baselearner is not selected in any fold, the function
#' treats its effect as constantly zero.
#' 
#'  
#' @return A list containing the elements \code{raw_results}, the 
#' \code{quantiles} and \code{mstops}. 
#' In \code{raw_results} and \code{quantiles}, each baselearner
#' selected with \code{which} in turn corresponds to a list
#' element. The quantiles are given as vector, matrix or list of
#' matrices depending on the nature of the effect. In case of functional
#' effects the list element in\code{quantiles} is a \code{length(levels)} times
#' \code{length(effect)} matrix, i.e. the rows correspond to the quantiles.
#' In case of coefficient surfaces, \code{quantiles} comprises a list of matrices,
#' where each list element corresponds to a quantile.
#' 
#' @examples 
#' if(require(refund)){
#' #########
#' # model with linear functional effect, use bsignal()
#' # Y(t) = f(t) + \int X1(s)\beta(s,t)ds + eps
#' set.seed(2121)
#' data1 <- pffrSim(scenario = "ff", n = 40)
#' data1$X1 <- scale(data1$X1, scale = FALSE)
#' dat_list <- as.list(data1)
#' dat_list$t <- attr(data1, "yindex")
#' dat_list$s <- attr(data1, "xindex")
#' 
#' ## model fit by FDboost 
#' m1 <- FDboost(Y ~ 1 + bsignal(x = X1, s = s, knots = 8, df = 3), 
#'               timeformula = ~ bbs(t, knots = 8), data = dat_list)
#'
#'}
#'               
#' \dontrun{             
#' # a short toy example with to few folds  
#' # and up to 200 boosting iterations 
#' bootCIs <- bootstrapCI(m1[200], B_inner = 2, B_outer = 5) 
#' 
#' # look at stopping iterations
#' bootCIs$mstops
#' 
#' # plot bootstrapped coefficient estimates
#' plot(bootCIs, ask = FALSE)
#' }
#' 
#' ## now speed things up by defining the inner resampling
#' ## function with parallelization based on mclapply (does not work on Windows)
#' 
#' my_inner_fun <- function(object){ 
#' cvrisk(object, folds = cvLong(id = object$id, weights = 
#' model.weights(object), 
#' B = 10 # 10-fold for inner resampling
#' ), mc.cores = 10) # use ten cores
#' }
#' 
#' \dontrun{
#' bootCIs <- bootstrapCI(m1, resampling_fun_inner = my_inner_fun)
#' }
#' 
#' ## We can also use the ... argument to parallelize the applyFolds
#' ## function in the outer resampling 
#' 
#' \dontrun{
#' bootCIs <- bootstrapCI(m1, mc.cores = 30)
#' }
#' 
#' ## Now let's parallelize the outer resampling and use 
#' ## crossvalidation instead of bootstrap for the inner resampling
#' 
#' my_inner_fun <- function(object){ 
#' cvrisk(object, folds = cvLong(id = object$id, weights = 
#' model.weights(object), type = "kfold", # use CV
#' B = 10, # 10-fold for inner resampling
#' ),
#' mc.cores = 10) # use ten cores
#' }
#' 
#' # use applyFolds for outer function to avoid messing up weights
#' my_outer_fun <- function(object, fun){
#' applyFolds(object = object,
#' folds = cv(rep(1, length(unique(object$id))), 
#' type = "bootstrap", B = 100), fun = fun,
#' mc.cores = 10) # parallelize on 10 cores
#' }
#' 
#' ######## Example for scalar-on-function-regression with bsignal() 
#' data("fuelSubset", package = "FDboost")
#' 
#' ## center the functional covariates per observed wavelength
#' fuelSubset$UVVIS <- scale(fuelSubset$UVVIS, scale = FALSE)
#' fuelSubset$NIR <- scale(fuelSubset$NIR, scale = FALSE)
#' 
#' ## to make mboost:::df2lambda() happy (all design matrix entries < 10)
#' ## reduce range of argvals to [0,1] to get smaller integration weights
#' fuelSubset$uvvis.lambda <- with(fuelSubset, (uvvis.lambda - min(uvvis.lambda)) /
#' (max(uvvis.lambda) - min(uvvis.lambda) ))
#' fuelSubset$nir.lambda <- with(fuelSubset, (nir.lambda - min(nir.lambda)) /
#' (max(nir.lambda) - min(nir.lambda) ))
#' 
#' ## model fit with scalar response and two functional linear effects 
#' ## include no intercept as all base-learners are centered around 0    
#' 
#' mod2 <- FDboost(heatan ~ bsignal(UVVIS, uvvis.lambda, knots = 40, df = 4, check.ident = FALSE) 
#'                + bsignal(NIR, nir.lambda, knots = 40, df=4, check.ident = FALSE), 
#'                timeformula = NULL, data = fuelSubset) 
#' 
#' 
#' \dontrun{
#' # takes some time, because of defaults: B_outer = 100, B_inner = 25
#' bootCIs <- bootstrapCI(mod2)
#' }
#' 
#' ## run with a larger number of outer bootstrap samples
#' ## and only 10-fold for validation of each outer fold
#' ## WARNING: This may take very long!
#' \dontrun{
#' bootCIs <- bootstrapCI(mod2, B_outer = 1000, B_inner = 10)
#' }
#' 
#' @export
bootstrapCI <- function(object, which = NULL, 
                        resampling_fun_outer = NULL,
                        resampling_fun_inner = NULL,
                        B_outer = 100,
                        B_inner = 25,
                        type_inner = c("bootstrap", "kfold", "subsampling"),
                        levels = c(0.05, 0.95),
                        ...)
{
  
  type_inner <- match.arg(type_inner)
  
  ########## check for scalar response #########
  scalarResp <- "FDboostScalar" %in% class(object)

  ########## define outer resampling function if NULL #########
  if(is.null(resampling_fun_outer)){
    
    resampling_fun_outer <- function(object, fun) applyFolds(object = object,
                                                             folds = cv(rep(1, length(unique(object$id))), type =
                                                                          "bootstrap", B = B_outer), fun = fun,
                                                             compress = FALSE, ...)

  }else{
    
    message("B_outer is ignored as resampling_fun_outer is not NULL.")
    
  }
      
  ########## define inner resampling function if NULL #########
  if(is.null(resampling_fun_inner)){
    
    if(scalarResp){
      
      resampling_fun_inner <- function(object) cvrisk(object, folds = cvLong(id = object$id, weights =
                                                                               model.weights(object), B = B_inner,
                                                                             type = type_inner))
      
    }else{ 
        
      resampling_fun_inner <- function(object) applyFolds(object, folds = cv(rep(1, length(unique(object$id))), 
                                                                             B = B_inner, 
                                                                             type = type_inner))
        
      }
      
    
  }else{
    
    message("B_inner is ignored as resampling_fun_inner is not NULL.")
    
  }
    
  # 'catch' error caused by using the cvrisk function for inner and outer resampling
  if(identical(resampling_fun_outer, cvrisk) & 
     identical(resampling_fun_inner, cvrisk)) 
    stop("Please specify a different outer resampling function.")
  
  ########## get coefficients ##########
  cat("Start computing bootstrap confidence intervals... \n")
  
  results <- resampling_fun_outer(object, 
                                  fun = function(mod)
                                  {
                                    ms <- mstop( resampling_fun_inner(mod) )
                                    coefs <- coef( mod[ms], which = which )
                                    return(list(coefs = coefs, ms = ms))
                                  }
                                  )
  
  cat("\n")
  
  coefs <- lapply(results, "[[", "coefs")
  mstops <- sapply(results, "[[", "ms")
  offsets <- t(sapply(coefs, function(x) x$offset$value))
  
  ########## format coefficients #########
  # number of baselearners
  nrEffects <- max(sapply(1:length(coefs), 
                          function(i) length(coefs[[i]]$smterms)))
  
  isFacSpecEffect <- sapply(1:nrEffects, 
                            function(i) "numberLevels" %in% names(coefs[[1]]$smterms[[i]]))
  # check for time-varying factor effects
  isFacEffect <- sapply(1:nrEffects, function(i) is.factor(coefs[[1]]$smterms[[i]]$x))

  # check for intercept
  withIntercept <- any(names(coefs[[1]]) == "intercept")
  
  if(withIntercept){
    intercept <- sapply(coefs, "[[", "intercept")
    if(!is.matrix(intercept)) intercept <- matrix(intercept, ncol = 1)
  } 
  
  # extract values
  listOfCoefs <- lapply(1:nrEffects, function(i)
  { 
    if(isFacSpecEffect[i]){
      # factor specific effect
      lapply(1:length(coefs), function(j) lapply(1:(coefs[[1]]$smterms[[i]]$numberLevels), 
                                                 function(k) coefs[[j]]$smterms[[i]][[k]]$value))
    }else{
      lapply(1:length(coefs), function(j) coefs[[j]]$smterms[[i]]$value)
    }
  })
  
  nnob <- names(object$baselearner)
  # exclude intercept in names
  names(listOfCoefs) <- nnob[(1+withIntercept):length(nnob)]
  
  # check for effect surfaces
  isSurface <- sapply(1:nrEffects, function(i) !is.null(coefs[[1]]$smterms[[i]]$y) )
  # do not treat time-varying factor effects as surfaces
  isSurface[isFacEffect] <- FALSE
  
  # reduce lists for non surface effects
  listOfCoefs[!isFacSpecEffect & !isSurface & !isFacEffect] <- 
    lapply(listOfCoefs[!isFacSpecEffect & !isSurface & !isFacEffect], 
           function(x) do.call("rbind", x))
  
  listOfCoefs[isFacSpecEffect | isSurface] <- 
    lapply(listOfCoefs[isFacSpecEffect | isSurface],
           function(x) lapply(x, function(y) do.call("rbind", lapply(y, c))))

  # add information about the values of the covariate
  # and change format
  for(i in 1:length(listOfCoefs)){
    
    if(isFacSpecEffect[i]){
      
      atx <- coefs[[1]]$smterms[[i]][[1]]$x
      
    }else{
      
      atx <- coefs[[1]]$smterms[[i]]$x 
      
    }

    aty <- NA
    if(isSurface[i] | isFacEffect[i]) aty <- coefs[[1]]$smterms[[i]]$y 
    if(isFacSpecEffect[i]) aty <- coefs[[1]]$smterms[[i]][[1]]$y 

    # format functional factors
    if(is.list(listOfCoefs[[i]]) & is.factor(atx)){
      
      # combine each factor level
      listOfCoefs[[i]] <- lapply(1:length(levels(droplevels(atx))),
                                 function(faclevnr) t(sapply(listOfCoefs[[i]], function(x) x[faclevnr,])))
      isSurface[i] <- FALSE
      
    }else if(is.list(listOfCoefs[[i]]) & !isFacSpecEffect[i]){ # effect surfaces
      
      listOfCoefs[[i]] <- do.call("rbind", lapply(listOfCoefs[[i]],c))
      
    }
    
    attr(listOfCoefs[[i]], "x") <- atx
    if(!is.na(sum(aty))) attr(listOfCoefs[[i]], "y") <- aty

    # add all plotting infos as attribute
    my_plot_info <- coefs[[1]]$smterms[[i]]
    my_plot_info$value <- NA
    attr(listOfCoefs[[i]], "plot_info") <- my_plot_info 

  }

  # add intercept and offset separately
  if(withIntercept){
    
    listOfCoefs <- c(offsets = list(offsets), 
                     intercept = list(intercept), 
                     listOfCoefs)
  
  }else{
    
    listOfCoefs <- c(offsets = list(offsets), 
                     listOfCoefs)
    
  }  
  
  my_plot_info <- coefs[[1]]$offset
  my_plot_info$value <- NA
  attr(listOfCoefs[[1]], "plot_info") <- my_plot_info
  attr(listOfCoefs[[1]], "x") <- coefs[[1]]$offset$x
  if(withIntercept){
    # for functional response, the intercept is a vector 
    attr(listOfCoefs[[2]], "plot_info") <- list(dim = 1)
    # for scalar response, the intercept is a scalar 
    if(class(object)[1] == "FDboostScalar") attr(listOfCoefs[[2]], "plot_info") <- list(dim = 0)
  } 
  
  # return raw results
  if(is.null(levels)) return(listOfCoefs)

  # define isSurface for quantile calculations
  isSurface <- c(rep(FALSE, 1 + withIntercept), isSurface)
  
  ########## calculate quantiles #########
  # create list for quantiles
  listOfQuantiles <- vector("list", length(listOfCoefs))
  
  # calculate quantiles
  for(i in 1:length(listOfCoefs)){
    
    # for matrix object
    if(is.matrix(listOfCoefs[[i]]) & !is.list(listOfCoefs[[i]])){
      
      listOfQuantiles[[i]] <- apply(listOfCoefs[[i]], 2, quantile, probs = levels)
      attr(listOfQuantiles[[i]], "x") <- attr(listOfCoefs[[i]], "x")
      if(!is.null(attr(listOfCoefs[[i]], "y")))
        attr(listOfQuantiles[[i]], "y") <- attr(listOfCoefs[[i]], "y")

    }else if(is.list(listOfCoefs[[i]])){ # functional factor variables
      
      listOfQuantiles[[i]] <- lapply(listOfCoefs[[i]], function(x) apply(t(x), 1, quantile, probs = levels))
      
    }else{# scalar case
      
      listOfQuantiles[[i]] <- quantile(listOfCoefs[[i]], probs = levels)
      
    }
         
  }
  
  # since coefficient surfaces are saved as vectors, reconstruct quantiles 
  # as coefficient surfaces and return a list of matrices, where each 
  # matrix corresponds to a quantile in levels
  if(sum(isSurface)!=0) listOfQuantiles[which(isSurface)] <- 
    lapply(listOfQuantiles[isSurface], 
           function(x){ 
             
             retL <- lapply(1:nrow(x), function(i) 
               matrix(x[i,], nrow = length(attr(x, "y"))))
             names(retL) <- levels
             return(retL)
             
           })
  
  # name rows of the matrices for non-surface effects
  if(sum(!isSurface)!=0) listOfQuantiles[which(!isSurface)] <- 
    lapply(listOfQuantiles[!isSurface],
           function(x){
             
             if(is.list(x)){
               
               for(j in 1:length(x)){
                 
                 if(!is.null(dim(x[[j]]))){
                   rownames(x[[j]]) <- levels
                 }else{
                   names(x[[j]]) <- levels
                 }
                 
               }
               
             }else{
             
               if(!is.null(dim(x))){
                 rownames(x) <- levels
               }else{
                 names(x) <- levels
               }
             
             }
             
             return(x)
             
           })
  
  # save names of baselearners
  names(listOfQuantiles) <- names(listOfCoefs)


  ########## return results #########
  ret <- list(raw_results = listOfCoefs,
              quantiles = listOfQuantiles,
              mstops = mstops,
              resampling_fun_outer = resampling_fun_outer,
              resampling_fun_inner = resampling_fun_inner,
              B_outer = B_outer,
              B_inner = B_inner,
              which = which,
              levels = levels, 
              yind = object$yind,
              family = object$family@name)
  
  class(ret) <- "bootstrapCI"
  
  return(ret)
  
}



#' Methods for objects of class bootstrapCI
#' 
#' Methods for objects that are fitted to compute bootstrap confidence intervals.
#' 
#' @param x an object of class \code{bootstrapCI}. 
#' @param which base-learners that are plotted 
#' @param pers plot coefficient surfaces as persp-plots? Defaults to \code{TRUE}.
#' @param commonRange, plot predicted coefficients on a common range, defaults to \code{TRUE}.
#' @param showQuantiles plot the 0.05 and the 0.95 Quantile of coefficients in 1-dim effects.
#' @param showNumbers show number of curve in plot of predicted coefficients, defaults to \code{FALSE}
#' @param ask defaults to \code{TRUE}, ask for next plot using \code{par(ask = ask)}? 
#' @param probs vector of quantiles to be used in the plotting of 2-dimensional coefficients surfaces,
#' defaults to \code{probs = c(0.25, 0.5, 0.75)}
#' @param ylim values for limits of y-axis
#' @param ... additional arguments passed to callies.
#' 
#' @details \code{plot.bootstrapCI} plots the bootstrapped coefficients.
#' 
#' @aliases print.bootstrapCI
#' 
#' @method plot bootstrapCI
#' 
#' @export
#' 
plot.bootstrapCI <- function(x, which = NULL, pers = TRUE,
                             commonRange = TRUE, showNumbers = FALSE, showQuantiles = TRUE,
                             ask = TRUE, 
                             probs = c(0.25, 0.5, 0.75), 
                             ylim = NULL, ...)
{
  
  stopifnot(class(x) == "bootstrapCI")
  
  boot_offset <- 0 
  
  if( names(x$raw_results)[1] == "offsets" ){
    
    boot_offset <- t(x$raw_results$offsets)
    x$raw_results$offsets <- NULL
    
    ## for scalar response, keep only one value per offset
    if(length(x$yind) == 1) boot_offset <- boot_offset[1, ]
    
  }
  
  if(is.null(which)) which <- 1:length(x$raw_results)
  
  if(length(which)>1) par(ask=ask)
  
  # find common range for all effects 
  if(commonRange & is.null(ylim)){
    ylim <- range(x$raw_results) 
    if(any(is.infinite(ylim))) ylim <- NULL
  }
  
  for(l in which){ # loop over effects

    ### prepare objects
    # coef() of a certain term
    temp_CI <- x$raw_results[[l]]
    
    temp <- attr(temp_CI, "plot_info")
    
    ## for interaction effects like "bhistx(x) %X% bolsc(z)"
    if(!is.null(temp$numberLevels)){
      temp <- temp[[1]]
      warning("Of the composed base-learner ", l, " only the first effect is plotted.")
    }
    
    ## write the rows of the matrix into a list, 
    ## i.e., coefficients of each fold are one list entry
    if(!is.list(temp_CI)){
      
      if(temp$dim >= 2){
        temp$value <- split(temp_CI, seq(nrow(temp_CI)))
      }else{ 
        ## temp$dim == 1 like in scalar response with bsignal()
        ## put each fold into one list entry 
        if(length(x$yind) <= 1 & x$family != "Binomial Distribution (similar to glm)"){ 
          # scalar response and not Binomial
          temp$value <- split(temp_CI, rep(1:x$B_outer, each = length(temp_CI)/x$B_outer))
        }else{
          temp$value <- split(temp_CI, rep(1:x$B_outer, length(temp_CI)/x$B_outer))
        }

      }
      
    }else{
 
      if(is.null(temp$numberLevels) & is.factor(temp$x)){
        
        ## for time-varying factor effects
        temp$value <- temp_CI

      }else{
        
        ## for interaction effects like "bhistx(x) %X% bolsc(z)" 
        ## the values are already a list
        temp$value <- lapply(temp_CI, function(x) x[1, ])
      
      }
    }    
    
    if(!is.null(temp$dim) && temp$dim == 2 & !is.factor(temp$x)){
      temp$value <- lapply(temp$value, function(xx) 
        matrix(xx, ncol = sqrt(length(xx)), nrow = sqrt(length(xx)), byrow = FALSE) )
    }
    
    plot_bootstrapped_coef(temp = temp, l = l, 
                           offset = boot_offset, yind = x$yind, 
                           pers = pers, 
                           showNumbers = showNumbers, showQuantiles = showQuantiles, 
                           probs = probs, ylim = ylim, ...)
    
    
  } # end loop over effects
  
}



#' @rdname plot.bootstrapCI
#' @method print bootstrapCI
#' @export
#'
print.bootstrapCI <- function(x, ...)
{

  stopifnot(class(x)=="bootstrapCI")

  cat("\n")
     
  cat("\t Bootstrapped confidence interval object of FDboost fit\n")
  
  cat("\n")
  cat("Coefficients:\n\t", names(x$quantiles), sep="\t", fill = TRUE)
  cat("\n")
  cat("\n")
  cat("Summary for stopping iterations of inner validation:\n\n")
  print(summary(x$mstops))
  cat("\n")
  invisible(x)
    
}
