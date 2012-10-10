predict.fgam <- function(object,newdata,type='response',
                         se.fit=FALSE,terms=NULL,PredOutOfRange=FALSE,...){
  
  call <- match.call()
  
  if (!missing(newdata)) {
    nobs <- nrow(as.matrix(newdata[[1]]))
    
    #if (!(all(names(newdata) %in% names(object$model)))) {
      stopifnot(length(unique(sapply(newdata, function(x) ifelse(is.matrix(x),
                                                                 nrow(x), length(x))))) == 1)
      gamdata <- list()
      varmap <- sapply(names(object$fgam$labelmap), function(x) all.vars(formula(paste("~",
                                                                                      x))))

      for (cov in names(newdata)) {
        trms <- which(sapply(varmap, function(x) any(grep(paste("^",
                                                                cov, "$", sep = ""), x))))
        if (length(trms) != 0) {
          for (trm in trms) {
            is.af <- trm %in% object$fgam$where$where.af
            is.lf <- trm %in% object$fgam$where$where.lf
           
            if (is.af) {
              af <- object$fgam$ft[[grep(paste(cov, "[,\\)]",
                                               sep = ""), names(object$fgam$ft))]]
              if (grepl(paste(cov, "\\.[ot]mat", sep = ""),
                        deparse(af$call$x))) {
                
                if(length(attr(newdata[[cov]],"L"))){
                  if(nobs==nrow(attr(newdata[[cov]],"L"))){
                    L <- attr(newdata[[cov]],"L")
                  }else{
                    L <- matrix(af$L[1,],nobs,length(af$L[1,]),byrow=T)
                  }
                }else{
                  L <- matrix(af$L[1,],nobs,length(af$L[1,]),byrow=T)
                }
                if(length(attr(newdata[[cov]],"tmat"))){
                  if(nobs==nrow(attr(newdata[[cov]],"tmat"))){
                    tmat <- attr(newdata[[cov]],"tmat")
                  }else{
                    tmat <- matrix(af$xind, ncol = length(af$xind),nrow = nobs,byrow=T)  
                  }
                }else{
                  tmat <- matrix(af$xind, ncol = length(af$xind),nrow = nobs,byrow=T)
                }
#                 if (any(apply(L, 2, function(x) length(unique(x))) !=
#                   1)) {
#                   stop("Error for ", names(varmap)[trm],
#                        "-- Prediction for af-terms with varying rows in integration operator L not implememented yet.")
#                 }
                if(af$presmooth & !type=='lpmatrix'){
                  newXfd <- fd(tcrossprod(af$Xfd$y2cMap,newdata[[cov]]),af$Xfd$basis)
                  newdata[[cov]] <- t(eval.fd(af$xind,newXfd))
                }
                
                if(af$Qtransform & !type=='lpmatrix'){
                  for(i in 1:nobs){
                    newdata[[cov]][i,] <- mapply(function(tecdf,x){tecdf(x)},
                                                tecdf=af$ecdflist,x=newdata[[cov]][i,])
                  }
                }

                gamdata[[paste(cov, ".omat", sep = "")]] <- newdata[[cov]]
                gamdata[[paste(cov, ".tmat", sep = "")]] <- tmat
                gamdata[[paste("L.", cov, sep = "")]] <- L
     
              }
            }
            if (is.lf) {
              lf <- object$fgam$ft[[grep(paste(cov, "[,\\)]",
                                               sep = ""), names(object$fgam$ft))]]
      
              if (grepl(paste(cov, "\\.[t]mat", sep = ""),
                        deparse(lf$call$x))) {
              
                if(length(attr(newdata[[cov]],"L"))){
                  if(nobs==nrow(attr(newdata[[cov]],"L"))){
                    L <- attr(newdata[[cov]],"L")
                  }else{
                    L <- matrix(lf$L[1,],nobs,length(lf$L[1,]),byrow=T)
                  }
                }else{
                  L <- matrix(lf$L[1,],nobs,length(lf$L[1,]),byrow=T)
                }
                if(length(attr(newdata[[cov]],"tmat"))){
                  if(nobs==nrow(attr(newdata[[cov]],"tmat"))){
                    tmat <- attr(newdata[[cov]],"tmat")
                  }else{
                    tmat <- matrix(lf$xind, ncol = length(lf$xind),nrow = nobs,byrow=T)  
                  }
                }else{
                  tmat <- matrix(lf$xind, ncol = length(lf$xind),nrow = nobs,byrow=T)
                }
#                 if (any(apply(L, 2, function(x) length(unique(x))) !=
#                   1)) {
#                   stop("Error for ", names(varmap)[trm],
#                        "-- Prediction for lf-terms with varying rows in integration operator L not implememented yet.")
#                 }
#                 predL <- matrix(L[1, ], byrow = TRUE,
#                                 nrow = nrow(newdata[[cov]]), ncol = ncol(L))
                
                if(af$presmooth){
                  newXfd <- fd(tcrossprod(lf$Xfd$y2cMap,newdata[[cov]]),lf$Xfd$basis)
                  newdata[[cov]] <- t(eval.fd(lf$xind,newXfd))
                }

                gamdata[[paste(cov, ".tmat", sep = "")]] <- tmat
                gamdata[[paste("L.", cov, sep = "")]] <- L*newdata[[cov]]
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
    #}
  }
  else { #when newdata not specified
    call$newdata <- eval(call$newdata)
    nobs <- object$fgam$nobs
  }
  
  if(PredOutOfRange){
    suppressMessages(trace(splines::spline.des, at=2, quote({
      outer.ok <- TRUE
    }),print=FALSE))
    on.exit(suppressMessages(try(untrace(splines::spline.des),silent=TRUE)))
  }
  
  oterms <- terms
  if(type %in% c('terms','iterms') & !all(terms %in% 
              c(names(object$smooth),attr(object$pterms,'term.labels')))){
    if(terms %in% unlist(varmap)){
      tnames <- c(names(object$smooth),attr(object$pterms,'term.labels'))
      terms <- tnames[grep(paste(string, "[,\\.\\)]",sep = ""), tnames)]
    }else{
      stop('Invalid terms specified')
    }
  }

 
  call[[1]] <- mgcv::predict.gam
  call$object <- as.name("object")
  call$terms <- terms
  
  res <- eval(call)
  
  if(type %in% c('terms','iterms'))
    colnames(res) <- oterms
  
  return(res)
}