#' Prediction from a fitted pfr model
#'
#' Takes a fitted \code{pfr}-object produced by \code{\link{pfr}} and produces predictions given a
#' new set of values for the model covariates or the original values used for the model fit.
#' Predictions can be accompanied by standard errors, based on the posterior distribution of the
#' model coefficients. This is a wrapper function for \code{\link{predict.gam}}()
#' @param object a fitted \code{pfr} object as produced by \code{\link{pfr}}
#' @param newdata a named list containing the values of the model covariates at which predictions
#' are required. If this is not provided then predictions corresponding to the original data are
#' returned. All variables provided to newdata should be in the format supplied to \code{\link{pfr}},
#' i.e., functional predictors must be supplied as matrices with each row corresponding to one
#' observed function. Index variables for the functional covariates are reused from the fitted model
#' object or alternatively can be supplied as attributes of the matrix of functional predictor values.
#' Any variables in the model not specified in newdata are set to their average values from the data
#' supplied during fitting the model
#' 
#' @param type character; see \code{\link{predict.gam}} for details
#' @param se.fit logical; see \code{\link{predict.gam}} for details
#' @param terms character see \code{\link{predict.gam}} for details
#' @param PredOutOfRange logical; if this argument is true then any functional predictor values in
#' newdata corresponding to \code{pfr} terms that are greater[less] than the maximum[minimum] of the
#' domain of the marginal basis for the rows of the tensor product smooth are set to the maximum[minimum]
#' of the domain.  If this argument is false, attempting to predict a value of the functional predictor
#' outside the range of this basis produces an error
#' @param ... additional arguments passed on to \code{\link{predict.gam}}
#' 
#' @return If \code{type == "lpmatrix"}, the design matrix for the supplied covariate values in long
#' format. If \code{se == TRUE}, a list with entries fit and se.fit containing fits and standard errors,
#' respectively. If \code{type == "terms" or "iterms"} each of these lists is a list of matrices of the
#' same dimension as the response for newdata containing the linear predictor and its se for each term
#' 
#' @author Mathew W. McLean \email{mathew.w.mclean@@gmail.com} and Fabian Scheipl
#' @seealso \code{\link{pfr}}, \code{\link[mgcv]{predict.gam}}
#' @export
#' @examples
#' ######### Octane data example #########
#' data(gasoline)
#' N <- length(gasoline$octane)
#' wavelengths = 2*450:850
#' nir = matrix(NA, 60,401)
#' test <- sample(60,20)
#' for (i in 1:60) nir[i,] = gasoline$NIR[i, ] # changes class from AsIs to matrix
#' y <- gasoline$octane
#' fit <- pfr(y~af(nir,argvals=wavelengths,k=c(6,6), m=list(c(2,2),c(2,2))),
#'               subset=(1:N)[-test])
#' preds <- predict(fit,newdata=list(nir=nir[test,]),type='response')
#' plot(preds,y[test])
#' abline(a=0,b=1)
#' @importFrom mgcv predict.gam
#' @importFrom splines spline.des
predict.pfr <- function (object, newdata, type = "response", se.fit = FALSE,
                          terms = NULL, PredOutOfRange = FALSE, ...) {
  
  if (!("pfr"%in% class(object))) {
    # Call predict.pfr_old()
    call <- sys.call()
    call[[1]] <- as.symbol("predict.pfr_old")
    ret <- eval(call)
    return(ret)
  }
  
  call <- match.call()
  string <- NULL
  if (!missing(newdata)) {
    nobs <- nrow(as.matrix(newdata[[1]]))
    
    stopifnot(length(unique(sapply(newdata, function(x)
      ifelse(is.matrix(x), nrow(x), length(x))))) == 1)
    gamdata <- list()
    varmap <- sapply(object$pfr$termnames, function(x)
      all.vars(formula(paste("~", x))), simplify = FALSE)
    for (cov in names(newdata)) {
      trms <- which(sapply(varmap, function(x)
        any(grep(paste("^", cov, "$", sep = ""), x))))
      J <- ncol(as.matrix(newdata[[cov]]))
      if (length(trms) != 0) {
        for (trm in trms) {
          is.af <- object$pfr$termtype[trm]=="af"
          is.lf <- object$pfr$termtype[trm]=="lf"
          if (is.af) {
            af <- object$pfr$ft[[grep(paste(cov, "[,\\)]", sep = ""),
                                      names(object$pfr$ft))]]
            
            if(J!=length(af$xind) & type!='lpmatrix'){
              stop(paste('Provided data for functional covariate', cov,
                         'does not have same observation times as original data',
                         sep=''))
            }
            L <- matrix(af$L[1, ], nobs, J, byrow = T)
            tmat <- matrix(af$xind,nobs,J,byrow=TRUE)
            if (grepl(paste(cov, "\\.[ot]mat", sep = ""),
                      deparse(af$call$x))) {
              if (length(attr(newdata, "L"))) {
                if(type!='lpmatrix'){
                  warning('Supplying new L matrix of quadrature weights only implemented for type=\'lpmatrix\' and supplied L will be ignored')
                }else{
                  if (sum(dim(as.matrix(attr(newdata, "L")))==dim(as.matrix(newdata[[cov]])))!=2) {
                    warning(paste('Supplied L matrix for',cov,'is not the same dimension as the matrix of observations and will be ignored',sep=''))
                    
                  }else{
                    L <- as.vector(attr(newdata, "L"))
                  }
                }
              }
              if (length(attr(newdata, "tmat"))) {
                if(type!='lpmatrix'){
                  warning('Supplying new tmat matrix of observation times only implemented for type=\'lpmatrix\' and supplied L will be ignored')
                }else{
                  if (sum(dim(as.matrix(attr(newdata, "tmat")))==dim(as.matrix(newdata[[cov]])))!=2) {
                    warning(paste('Supplied tmat matrix for',cov,'is not the same dimension as the matrix of observations and will be ignored',sep=''))
                  }else{
                    tmat <- as.vector(attr(newdata, "tmat"))
                  }
                }
              }
              if (PredOutOfRange) {
                newdata[[cov]][newdata[[cov]]>af$Xrange[2]]  <- af$Xrange[2]
                newdata[[cov]][newdata[[cov]]<af$Xrange[1]]  <- af$Xrange[1]
              }
              if (!is.null(af$presmooth)) {
                if(type=='lpmatrix' & J!=length(af$xind)){
                  warning('Presmoothing of new functional covariates is only implemented for when new covariates observed at same time points as original data. No presmoothing of new covariates done.')
                } else if (is.logical(af$presmooth)) {
                  # af_old term
                  if (af$presmooth) {
                    newXfd <- fd(tcrossprod(af$Xfd$y2cMap,
                                            newdata[[cov]]), af$Xfd$basis)
                    newdata[[cov]] <- t(eval.fd(af$xind, newXfd))
                  }
                } else {
                  # af() term
                  newdata[[cov]] <- af$prep.func(newX = newdata[[cov]])  
                }
                
              }
              
              #if (af$Qtransform) {
              #  if(type=='lpmatrix' & J!=length(af$xind)){
              #    stop('Prediction with quantile transformation only implemented for when new data observation times match original data observation times')
              #  }
              #  for (i in 1:nobs) {
              #    newdata[[cov]][i, ] <- mapply(function(tecdf,
              #                                           x) {
              #      tecdf(x)
              #    }, tecdf = af$ecdflist, x = newdata[[cov]][i,
              #                                               ])
              #  }
              #}
              if(type=='lpmatrix') newdata[[cov]] <- as.vector(newdata[[cov]])
              gamdata[[paste(cov, ".omat", sep = "")]] <- newdata[[cov]]
              gamdata[[paste(cov, ".tmat", sep = "")]] <- tmat
              gamdata[[paste("L.", cov, sep = "")]] <- L
            }
          }
          if (is.lf) {
            lf <- object$pfr$ft[[grep(paste(cov, "[,\\)]",
                                             sep = ""), names(object$pfr$ft))]]
            if(J!=length(lf$xind) & type!='lpmatrix'){
              stop(paste('Provided data for functional covariate',cov,'does not have same observation times as original data',sep=''))
            }
            L <- matrix(lf$L[1, ], nobs, J, byrow = T)
            tmat <- matrix(lf$xind,nobs,J,byrow=TRUE)
            if (grepl(paste(cov, "\\.[t]mat", sep = ""),
                      deparse(lf$call$x))) {
              if (length(attr(newdata, "L"))) {
                if(type!='lpmatrix'){
                  warning('Supplying new L matrix of quadrature weights only implemented for type=\'lpmatrix\' and supplied L will be ignored')
                }else{
                  if (sum(dim(as.matrix(attr(newdata, "L")))==dim(as.matrix(newdata[[cov]])))!=2) {
                    warning(paste('Supplied L matrix for',cov,'is not the same dimension as the matrix of observations and will be ignored',sep=''))
                    
                  }else{
                    L <- as.vector(attr(newdata, "L"))
                  }
                }
              }
              if (length(attr(newdata, "tmat"))) {
                if(type!='lpmatrix'){
                  warning('Supplying new tmat matrix of observation times only implemented for type=\'lpmatrix\' and supplied L will be ignored')
                }else{
                  if (sum(dim(as.matrix(attr(newdata, "tmat")))==dim(as.matrix(newdata[[cov]])))!=2) {
                    warning(paste('Supplied tmat matrix for',cov,'is not the same dimension as the matrix of observations and will be ignored',sep=''))
                  }else{
                    tmat <- as.vector(attr(newdata, "tmat"))
                  }
                }
              }
              if (!is.null(lf$presmooth)) {
                if(type=='lpmatrix' & J!=length(lf$xind)){
                  warning('Presmoothing of new functional covariates is only implemented for when new covariates observed at same time points as original data. No presmoothing of new covariates done.')
                } else if (is.logical(lf$presmooth)) {
                  # lf_old term
                  if (lf$presmooth) {
                    newXfd <- fd(tcrossprod(lf$Xfd$y2cMap,
                                            newdata[[cov]]), lf$Xfd$basis)
                    newdata[[cov]] <- t(eval.fd(lf$xind, newXfd))
                  }
                } else {
                  # lf() term
                  newdata[[cov]] <- lf$prep.func(newX = newdata[[cov]])  
                }
                
              }
              if(type=='lpmatrix') newdata[[cov]] <- as.vector(newdata[[cov]])
              gamdata[[paste(cov, ".tmat", sep = "")]] <- tmat
              gamdata[[paste("L.", cov, sep = "")]] <- L *
                newdata[[cov]]
            }
          }
          if (!(is.af || is.lf)) {
            gamdata[[cov]] <- drop(newdata[[cov]])
          }
        }
      }
    }
    gamdata <- list2df(gamdata)
    call[["newdata"]] <- gamdata
  }
  else {
    call$newdata <- eval(call$newdata)
    nobs <- object$pfr$nobs
  }
  if (PredOutOfRange) {
    suppressMessages(trace(splines::spline.des, at = 2, quote({
      outer.ok <- TRUE
    }), print = FALSE))
    on.exit(suppressMessages(try(untrace(splines::spline.des),
                                 silent = TRUE)))
  }
  oterms <- terms
  if (type %in% c("terms", "iterms") & !all(terms %in% c(names(object$smooth),
                                                         attr(object$pterms, "term.labels")))) {
    if (terms %in% unlist(varmap)) {
      tnames <- c(names(object$smooth), attr(object$pterms,
                                             "term.labels"))
      terms <- tnames[grep(paste(string, "[,\\.\\)]", sep = ""),
                           tnames)]
    }
    else {
      stop("Invalid terms specified")
    }
  }
  call[[1]] <- mgcv::predict.gam
  call$object <- as.name("object")
  call$terms <- terms
  res <- eval(call)
  if (type %in% c("terms", "iterms"))
    colnames(res) <- oterms
  return(res)
}

