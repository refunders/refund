af <- function(X, xind = seq(0, 1, l = ncol(X)), basistype = c("te","t2", "s"), 
               integration = c("simpson", "trapezoidal", "riemann"), 
                L = NULL, splinepars = list(bs = "ps", 
                k= c(min(ceiling(nrow(X)/5),20),min(ceiling(ncol(X)/5),20)),
                m = list(c(2, 2), c(2, 2))), presmooth = TRUE,Xrange=range(X),Qtransform=FALSE) {
      
  n=nrow(X)
  nt=ncol(X)
  basistype <- match.arg(basistype)
  integration <- match.arg(integration)
  if(is.null(splinepars$bs)) splinepars$bs <- 'ps'
  if(is.null(splinepars$k)) splinepars$k <- c(min(ceiling(nrow(X)/5),20),
                                                min(ceiling(ncol(X)/5),20))
  if(is.null(splinepars$m)) splinepars$m = list(c(2, 2), c(2, 2))

  xindname <- paste(deparse(substitute(X)), ".omat", sep = "")
  tindname <- paste(deparse(substitute(X)), ".tmat", sep = "")
  Lname <- paste("L.", deparse(substitute(X)), sep = "")

  if (is.null(dim(xind))) {
    xind <- t(xind)
    stopifnot(ncol(xind) == nt)
    if (nrow(xind) == 1) {
      xind <- matrix(as.vector(xind), nrow = n, ncol = nt, 
                     byrow = T)
    }
    stopifnot(nrow(xind) == n)
  }
  
  Xfd=NULL
  if(presmooth){
    bbt=create.bspline.basis(rangeval=range(xind),nbasis=ceiling(nt/4),                                   
                             norder=splinepars$m[[2]][1]+2, breaks=NULL)
    
    # pre-smooth functional predictor
    temp <- smooth.basisPar(t(xind),t(X),bbt,int2Lfd(splinepars$m[[2]][2])) 
    Xfd <- temp$fd
    Xfd$y2cMap <-temp$y2cMap
    X <- t(sapply(1:n,function(i){eval.fd(xind[i,],Xfd[i])}))
    
    # need to check that smoothing didn't change range of data
    if(!Qtransform){
      if(max(X)>Xrange[2]){
        Xrange[2] <- max(X)
      } 
      if(min(X)<Xrange[1]){
        Xrange[1] <- min(X)
      }
    }
  }
  
  ecdf=NULL
  if(Qtransform){
    Xrange <- c(0,1)
    X <- apply(X, 2, function(x){ (rank(x)-1)/(length(x)-1)} )
    # need to keep ecdf's for prediction later
    ecdflist <- apply(X, 2, ecdf)     
  }

  if (!is.null(L)) {
    stopifnot(nrow(L) == n, ncol(L) == nt)
  }else {
    L <- switch(integration, simpson = {
      ((xind[, nt] - xind[, 1])/nt)/3 * matrix(c(1,rep(c(4, 2), length = nt - 2), 1), nrow = n, 
                                                       ncol = nt, byrow = T)
    }, trapezoidal = {
      diffs <- t(apply(xind, 1, diff))
      0.5 * cbind(diffs[, 1], t(apply(diffs, 1, filter,filter = c(1, 1)))[, -(nt - 1)], 
                     diffs[,(nt - 1)])
    }, riemann = {
      diffs <- t(apply(xind, 1, diff))
      cbind(rep(mean(diffs), n), diffs)
    })
  }

  data <- list(X, xind, L)
  names(data) <- c(xindname, tindname, Lname)
  splinefun <- as.symbol(basistype)
  frmls <- formals(getFromNamespace(deparse(splinefun), ns = "mgcv"))
  frmls <- modifyList(frmls[names(frmls) %in% names(splinepars)], 
                      splinepars)
  call <- as.call(c(list(splinefun, x = as.symbol(substitute(xindname)), 
                         z = as.symbol(substitute(tindname)), by = as.symbol(substitute(Lname))), 
                    frmls))
  res <-list(call = call, data = data, xind = xind[1,], L = L, xindname = xindname, tindname=tindname,
             Lname=Lname,Qtransform=Qtransform,presmooth=presmooth,Xrange=Xrange)
  if(Qtransform) res$ecdflist <- ecdflist
  if(presmooth) res$Xfd <- Xfd
  return(res)
}

      