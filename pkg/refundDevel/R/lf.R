lf <- function(X, xind = seq(0, 1, l = ncol(X)), 
               integration = c("simpson", "trapezoidal", "riemann"), 
               L = NULL, splinepars = list(bs = "ps", k= min(ceiling(n/4),40),
                                  m = c(2, 2)), presmooth = TRUE, Xrange=range(X)) {
  
  n=nrow(X)
  nt=ncol(X)
  integration <- match.arg(integration)
  if(is.null(splinepars$bs)) splinepars$bs <- 'ps'
  if(is.null(splinepars$k)) splinepars$k <- min(ceiling(n/4),40)
  if(is.null(splinepars$m)) splinepars$m = c(2, 2)

  tindname <- paste(deparse(substitute(X)), ".tmat", sep = "")
  LXname <- paste("L.", deparse(substitute(X)), sep = "")
  basistype = "s"
  
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
                             norder=splinepars$m[1]+2, breaks=NULL)
    
    # pre-smooth functional predictor
    temp <- smooth.basisPar(t(xind),t(X),bbt,int2Lfd(splinepars$m[2])) 
    Xfd <- temp$fd
    Xfd$y2cMap <-temp$y2cMap
    X <- t(sapply(1:n,function(i){eval.fd(xind[i,],Xfd[i])}))
    
    # need to check that smoothing didn't change range of data
    if(max(X)>Xrange[2]){
      Xranges[2] <- max(X)
    } 
    if(min(X)<Xrange[1]){
      Xranges[1] <- min(X)
    }
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
  LX <- L*X
  data <- list(xind, LX)
  names(data) <- c(tindname, LXname)
  splinefun <- as.symbol(basistype)
  frmls <- formals(getFromNamespace(deparse(splinefun), ns = "mgcv"))
  frmls <- modifyList(frmls[names(frmls) %in% names(splinepars)], 
                      splinepars)
  call <- as.call(c(list(splinefun, x = as.symbol(substitute(tindname)), 
                         by = as.symbol(substitute(LXname))),frmls))
  res <-list(call = call, data = data, xind = xind[1,], L = L, tindname=tindname,
             LXname=LXname,presmooth=presmooth,Xrange=Xrange)
  if(presmooth) res$Xfd <- Xfd
  return(res)
}

