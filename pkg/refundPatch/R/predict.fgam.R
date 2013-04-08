predict.fgam <- function (object, newdata, type = "response", se.fit = FALSE, 
          terms = NULL, PredOutOfRange = FALSE, ...) 
{
  call <- match.call()
  string <- NULL
  if (!missing(newdata)) {
    nobs <- nrow(as.matrix(newdata[[1]]))

    stopifnot(length(unique(sapply(newdata, function(x) ifelse(is.matrix(x), 
                                                               nrow(x), length(x))))) == 1)
    gamdata <- list()
    varmap <- sapply(names(object$fgam$labelmap), function(x) all.vars(formula(paste("~", 
                                                                                     x))))
    for (cov in names(newdata)) {
      trms <- which(sapply(varmap, function(x) any(grep(paste("^", 
                                                              cov, "$", sep = ""), x))))
      J <- ncol(as.matrix(newdata[[cov]]))
      if (length(trms) != 0) {
        for (trm in trms) {
          is.af <- trm %in% object$fgam$where$where.af
          is.lf <- trm %in% object$fgam$where$where.lf
          if (is.af) {
            af <- object$fgam$ft[[grep(paste(cov, "[,\\)]", 
                                             sep = ""), names(object$fgam$ft))]]

            if(J!=length(af$xind) & type!='lpmatrix'){
              stop(paste('Provided data for functional covariate',cov,'does not have same observation times as original data',sep=''))
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
              if (af$presmooth) {
                if(type=='lpmatrix' & J!=length(af$xind)){
                  warning('Presmoothing of new functional covariates is only implemented for when new covariates observed at same time points as original data. No presmoothing of new covariates done.')
                }else{
                  newXfd <- fd(tcrossprod(af$Xfd$y2cMap, 
                                          newdata[[cov]]), af$Xfd$basis)
                  newdata[[cov]] <- t(eval.fd(af$xind, 
                                              newXfd))
                }
              }
              if (af$Qtransform) {
                if(type=='lpmatrix' & J!=length(af$xind)){
                  stop('Prediction with quantile transformation only implemented for when new data observation times match original data observation times')
                }
                for (i in 1:nobs) {
                  newdata[[cov]][i, ] <- mapply(function(tecdf, 
                                                         x) {
                    tecdf(x)
                  }, tecdf = af$ecdflist, x = newdata[[cov]][i, 
                                                             ])
                }
              }
              if(type=='lpmatrix') newdata[[cov]] <- as.vector(newdata[[cov]])
              gamdata[[paste(cov, ".omat", sep = "")]] <- newdata[[cov]]
              gamdata[[paste(cov, ".tmat", sep = "")]] <- tmat
              gamdata[[paste("L.", cov, sep = "")]] <- L
            }
          }
          if (is.lf) {
            lf <- object$fgam$ft[[grep(paste(cov, "[,\\)]", 
                                             sep = ""), names(object$fgam$ft))]]
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
              if (lf$presmooth) {
                if(type=='lpmatrix' & J!=length(lf$xind)){
                  warning('Presmoothing of new functional covariates is only implemented for when new covariates observed at same time points as original data. No presmoothing of new covariates done.')
                }else{
                  newXfd <- fd(tcrossprod(lf$Xfd$y2cMap, 
                                          newdata[[cov]]), lf$Xfd$basis)
                  newdata[[cov]] <- t(eval.fd(lf$xind, 
                                              newXfd))
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
    nobs <- object$fgam$nobs
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
