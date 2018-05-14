#' @import fda

#library(fda)


#################################################################
constMin1 <- function(a, mu)
  # max  a'x such that ||x||_mu \e 1.
{
  nvar <- length(a)
  x <- rep(1,nvar)
  aa <- abs(a)
  aa.sort <- sort(aa, decreasing=T, index.return=T, method="quick")
  aas <- aa.sort$x
  aas.sum <- sum(aas)
  aas.cum <- cumsum(aas)
  flag <- c(aas<=mu*aas.cum/(1+(1:nvar-1)*mu), TRUE)
  mc <- which(flag)[1]-1

  tmp <- mu * aas.cum[mc] / (1+(mc-1)*mu)
  x <- c(aas[1:mc] - tmp, rep(0, nvar-mc))

  b <- sort(aa.sort$ix, index.return=T)
  x <- x[b$ix]

  norm= sqrt(mu*(sum(abs(x)))^2+(1-mu)*sum(x^2))
  x <- as.numeric(lapply(a,sign)) * x / norm

  return(list(x=x, norm.lambda=norm, m=mc,I=sort(aa.sort$ix[1:mc])))
}


#################################################################
maxLinearH <- function(b, LPhi, tau, lambda){
  # b: R^n, LPhi=L %*% phi.coef, b*p, Phi(t)=basis(t) %*% phi.coef, Phi(t): T*p, basis: T*b, phi.coef: b*p
  # max b'w+<<Phi, Psi>> s.t. ||w||_2^2+tau|||psi|||_lambda^2 \le 1
  # return w, Psi=Psi.coef

  a <- sqrt(crossprod(LPhi^2,rep(1,nrow(LPhi))))
  ret <- constMin1(a, lambda)  # u=L %*% Psistar.coef
  u <- ret$x
  I=which(a!=0)
  LPsistar=LPhi
  LPsistar[,I] <- scale(LPhi[,I], center=FALSE, scale=a[I]/u[I])
  PsiPhi <- sum(LPsistar * LPhi)
  norm <- sqrt(tau*sum(b^2)+ PsiPhi^2)
  coef <- sqrt(tau)/norm
  w.ret <- coef * b
  LPsi.ret <- PsiPhi * coef / tau * LPsistar #b*p
  x <- c(w.ret, as.vector(LPsi.ret))
  return(list(norm=norm, w=w.ret, LPsi=LPsi.ret, x=x, u=u, a=a))
}




#################################################################


constOptOrth.H <- function(gamma, orthConst.mtx, t0, Lambda.all, nobs, nbasis, tau, lambda, tol){
  #find t to minimize H, or find t so that f(gamma+orthConst.mtx %*% t) is orthogonal to the space spanned by columns in orthConst.mtx.
  # f(gamma + orthConst.mtx %*% t) is solution to problem (7.10) in paper.
  #gamma: (ncurve+nbasis*nvarX)*1, orthConst.mtx: (ncurve+kcomp)*(ncurve+nbasis*nvarX),  orthConst.mtx %*% t(orthConst.mtx)=I
  # assuming the number of basis for all alpha's are the same

  t <- -tcrossprod(orthConst.mtx, t(gamma))
  kcomp <- dim(orthConst.mtx)[1]
  orthConst.mtx.1=orthConst.mtx[,1:nobs]               #kcomp*ncurve

  gamma.orth <- gamma + crossprod(orthConst.mtx, t)   #check orhConst.mtx %*% t(orthConst.mtx)=I
  gammat=gamma.orth
  t <- rep(0, kcomp)
  Phi.mtx <- matrix(gammat[-(1:nobs)], nrow=nbasis)       #b*p
  ret <- maxLinearH(gammat[1:nobs], Phi.mtx, tau, lambda)
  w <- ret$w
  LPsi.mtx <- ret$LPsi      #b*p
  tau.sqrt <- sqrt(tau)
  h <- (tau.sqrt * ret$norm) * tcrossprod(orthConst.mtx, t(ret$x))
  H <- sum(h^2)

  if(sum(abs(t0))>0)
  {
    gammat=as.vector(gamma.orth + crossprod(orthConst.mtx, t0))   #gamma.orthogonal
    Phi.mtx.new <- matrix(gammat[-(1:nobs)], nrow=nbasis)       #b*p
    ret.new <- maxLinearH(gammat[1:nobs], Phi.mtx.new, tau, lambda)
    w.new <- ret.new$w
    LPsi.mtx.new <- ret.new$LPsi      #b*p
    h.new <- (tau.sqrt * ret.new$norm) * tcrossprod(orthConst.mtx,t(ret.new$x))   #c(w.new, as.vector(LPsi.mtx.new))))
    H.new <- sum(h.new^2)
    if(H.new<H)
    { # print("use new t")
      t <- t0
      w <- w.new
      LPsi.mtx <- LPsi.mtx.new
      Phi.mtx <- Phi.mtx.new
      h <- h.new
      H <- H.new
      ret <- ret.new
    }
  }

  E.mtx <- scale(Phi.mtx, center=FALSE, scale=ret$a)  #b*p (e1,...,ep), each ei corresponds to b basis
  I=which(ret$u!=0)
  I.len <- length(I)
  alpha <- lambda / (1+(I.len-1)*lambda)

  ones.nb <- rep(1,nbasis)
  count1 <- 0
  nonorth.ind=FALSE
  while(H>tol){
    #        print(c("count1", count1, H, I.len, I))
    count1 <- count1+1
    if(count1 > 20){
      nonorth.ind=TRUE
      break;
    }

    t1 <- rep((I-1)*nbasis, each=nbasis)
    t2 <- rep(1:(nbasis), I.len)
    col.id <- t1+t2
    orthConst.2.I <- orthConst.mtx[, col.id+nobs]      #kcomp*(length(I)*b)
    E.I.mtx <- as.matrix(E.mtx[,I])  #b*length(I)

    orthConst.2.I.arr <- array(orthConst.2.I, c(kcomp, nbasis, I.len))
    Xi <- sapply(1:I.len, function(i) crossprod(t(orthConst.2.I.arr[,,i]), E.I.mtx[,i]))


    nu2 <- (1-lambda)*sum(ret$a * ret$u)
    D.vec <- nu2 * ret$u[I] / ret$a[I]
    B2 <- tcrossprod(Xi * tcrossprod(rep(1, kcomp), 1-D.vec),Xi)
    tmp <- as.numeric(apply(Xi,1,sum))
    B3 <- tcrossprod(tmp) # %*% t(tmp)
    tmp=lapply(1:I.len, function(i) Lambda.all[[I[i]]]*D.vec[i])
    B4=Reduce("+",tmp)

    A <- tau * tcrossprod(orthConst.mtx.1) + 1/(1-lambda) * (B2 - alpha*B3 + B4)
    A <- (A+t(A))/2

    deriv.1 <- A %*% h

    t.inc=-solve(A, h)

    t.new <- t + t.inc

    gammat=as.vector(gamma.orth + crossprod(orthConst.mtx, t.new))   #gamma.orthogonal
    Phi.mtx.new <- matrix(gammat[-(1:nobs)], nrow=nbasis)       #b*p
    ret <- maxLinearH(gammat[1:nobs], Phi.mtx.new, tau, lambda)
    w.new <- ret$w
    LPsi.mtx.new <- ret$LPsi      #b*p
    h.new <- (tau.sqrt * ret$norm) * (tcrossprod(orthConst.mtx, t(ret$x))) #c(w.new, as.vector(LPsi.mtx.new)))))
    H.new <- sum(h.new^2)

    count2 <- 0
    #        print(c(H.new, H+(1e-4)*t(t.inc) %*% deriv.1))
    while(H.new>tol & H.new>H+(1e-4)*crossprod(t.inc, deriv.1)){
      count2 <- count2 + 1
      #           print(c("count.2",count2, H, H.new, sqrt(sum(t.inc^2))))
      if(count2 > 10){
        break;
      }

      t.inc <- 0.3 * t.inc
      t.new <- t + t.inc
      gammat=as.vector(gamma.orth + crossprod(orthConst.mtx, t.new))   #gamma.orthogonal
      Phi.mtx.new <- matrix(gammat[-(1:nobs)], nrow=nbasis)       #b*p
      ret <- maxLinearH(gammat[1:nobs], Phi.mtx.new, tau, lambda)
      w.new <- ret$w
      LPsi.mtx.new <- ret$LPsi      #b*p
      h.new <- (tau.sqrt * ret$norm) * (tcrossprod(orthConst.mtx, t(ret$x)))   #c(w.new, as.vector(LPsi.mtx.new)))))
      H.new <- sum(h.new^2)
    }

    t <- t.new
    w <- w.new
    LPsi.mtx <- LPsi.mtx.new
    Phi.mtx <- Phi.mtx.new
    E.mtx <- scale(Phi.mtx, center=FALSE, scale=ret$a)  #b*p (e1,...,ep), each ei corresponds to b basis
    I <- which(ret$u!=0)
    I.len <- length(I)
    alpha <- lambda/(1+(I.len-1)*lambda)
    h <- h.new
    H <- H.new
  } # end of while(H>tol)
  #print(c("count1", count1, H))
  return(list(gamma=gammat, t=t, count1=count1, w=w, LPsi=LPsi.mtx, H=H, I.len=I.len, x=ret$x, u=ret$u))
}


#################################################################

getComp.riemann.sof.int <- function(XbTransInv, Y, orthConst.mtx, nbasis, K.comp, tau, lambda, tol=1e-11, tol1=1e-11, repeat.solution=1){
  # internal function for getComp, not called by users
  # XbTransInv: a matrix
  time <- proc.time()

  nobs <- nrow(Y)
  nvarX <- ncol(XbTransInv) / nbasis
  Y.cent <- scale(Y, center=T, scale=FALSE)
  muy <- attributes(Y.cent)$`scaled:center`
  XY.prod <- crossprod(XbTransInv, Y.cent)
  Beta <- NULL
  for(ncomp in 1:K.comp){
    if(!is.null(Beta)){
      tmp <- c(rep(0, nobs), crossprod(XbTransInv %*% beta, XbTransInv))
      tmp <- tmp - crossprod(orthConst.mtx, orthConst.mtx %*% tmp)
      orthConst.mtx <- rbind(orthConst.mtx, matrix(tmp/sqrt(sum(tmp^2)), nrow=1))   #(nobs+kcomp) %*% (nobs+nbasis*nvarX)
    }
    nrow.orthconst <- nobs+ncomp-1
    orthConst.arr.2=array(orthConst.mtx[,-(1:nobs)], c(nobs+ncomp-1, nbasis, nvarX))
    orthConst.lst.2 <- lapply(1:nvarX, function(jx) orthConst.arr.2[,,jx])
    Lambda.all <- lapply(orthConst.lst.2, tcrossprod)


    qq <- rep(0, repeat.solution)
    mqq <- 0
    for(jj in 1:repeat.solution){
      #            pt <- proc.time()
      gamma <- rnorm(nobs+nbasis*nvarX, 0,10)
      gamma <- gamma/sqrt(gamma^2)
      gamma.old <- rep(0,nobs+nbasis*nvarX)
      count <- 0

      tmp <- crossprod(gamma[-(1:nobs)], XY.prod)
      value <- sum(tmp^2)
      temp <- tcrossprod(XY.prod, tmp)
      value.old <- 0
      t0 <- rep(0, nrow.orthconst)
      while((count<50)&(min(sum((gamma[-(1:nobs)]-gamma.old[-(1:nobs)])^2), sum((gamma[-(1:nobs)]+gamma.old[-(1:nobs)])^2))>sum(gamma[-(1:nobs)]^2)*1e-10)){
        if(count>0)
        {value.old <- value}

        count <- count+1
        ret <- constOptOrth.H(c(rep(0,nobs),temp), orthConst.mtx, t0, Lambda.all, nobs, nbasis, tau, lambda, tol)
        gamma.old=gamma
        t0 <- ret$t
        gamma=ret$x
        tmp <- crossprod(gamma[-(1:nobs)], XY.prod)
        value <- sum(tmp^2)
        #print(value)
        temp <- tcrossprod(XY.prod, tmp)
      }
      qq[jj] <- value
      if(qq[jj]>mqq){
        mqq <- qq[jj]
        beta=gamma[-(1:nobs)]
      }

    }
    #        print(proc.time()-time)
    Beta=rbind(Beta, beta)     #K*(bp)
  }
  return(list(Beta=Beta, Y=Y, Y.cent=Y.cent, muy=muy, XbTransInv=XbTransInv,XY.prod=XY.prod))
}



cv.folds <- function(n,nfolds=5)
  ## Randomly split the n samples into folds
  ## Returns a list of nfolds lists of indices, each corresponding to a fold
{
  return(split(sample(n),rep(1:nfolds,length=n)))
}


#################################################################
#' @export
cv.hd <- function(X, Y, t.x.list, t.y,   K.cv=5, s.n.basis=25, t.n.basis=50, thresh=0.01)
{
  if(!is.list(X))
  {stop("Error!!: X must be a list!")}
  if (sum(sapply(1:length(X),function(k){!is.matrix(X[[k]])})))
  {stop("Error!!: X must be a list and all its components must be matrix!")
  }
  if(!is.list(t.x.list))
  {stop("Error!!: t.x.list must be a list!")}
  if (length(X)!=length(t.x.list))
  {stop("Error!!: both X and t.x.list must be lists and they have the same numbers of  components!")
  }
  dim.1=sapply(1:length(X),function(k){dim(X[[k]])[1]})
  if((length(unique(dim.1))!=1))
  {stop("Error!!: all components of X must be matrix and have the same numbers of rows!")
  }
  if((dim(X[[1]])[1]!=dim(Y)[1]))
  {stop("Error!!: the number of observations of X (that is, the number of rows of each component of X) must be equal to the number of observations of Y (that is, the number of rows of Y)!")
  }
  if(sum(sapply(1:length(X), function(k){dim(X[[k]])[2]!=length(t.x.list[[k]])}))!=0)
  {stop("Error!!: The number of columns of each component of X must be equal to the length of the corresponsing component of t.x.list!")
  }
  if(dim(Y)[2]!=length(t.y))
  {stop("Error!!: the number of columns of Y must be equal to the length of the vector t.y of the observation points!")
  }

  tol=1e-11
  tol1=1e-11
  ytrange=c(0,1)
  length.Y=length(t.y)


  y.int.weights=(ytrange[2]-ytrange[1])/length(t.y)*diag(length(t.y))

  kappa=c(1e-10,1e-8,1e-6,1e-4,1e-2,1,100)
  y.basis<- create.bspline.basis(ytrange, t.n.basis, 4)
  y.params=list()
  y.params[[1]]=ytrange
  y.params[[2]]=t.y
  y.params[[3]]=y.basis
  y.params[[4]]=kappa
  y.params[[5]]= y.int.weights

  A=  cbind(c(0.1,1,10,100), c(0.1,0.2,0.3,0.4))
  params.set <- rbind( cbind(A, rep(1e-6,dim(A)[1])),cbind(A, rep(1e-4,dim(A)[1])),cbind(A, rep(1e-2,dim(A)[1])))

  maxK.ret <-determineMaxKcomp.hd(X, Y, t.x.list, t.y, s.n.basis, params.set,thresh)
  ind.vector=FALSE
  if(is.vector(Y))
  {Y=matrix(Y,length(Y),1)
  ind.vector=TRUE
  }
  max.K <- maxK.ret$max.K
  Xb.lst <- maxK.ret$Xb.lst
  normTransInv.lst <- maxK.ret$normTransInv.lst

  numParm <- dim(params.set)[1]
  errors <- list()
  for(j in 1:numParm)  #numParm need to be specified
    errors[[j]] <- matrix(0, length(y.params[[4]]), max.K[j])
  dims <- dim(Y)
  nsample <- dims[1]
  nvarY <- dims[2]
  all.folds <- cv.folds(nsample, K.cv)
  nvarX <- length(X)

  length.uniqeta <- length(unique(params.set[,3]))
  print("Starting the cross-validation procedure:")
  for(i in 1:K.cv){
    omit <- all.folds[[i]]
    Y.icv <- Y[-omit,]
    y.valid <- Y[omit,]

    if(ind.vector)
    {
      Y.icv <- matrix(Y[-omit,],ncol=1)
      y.valid <- matrix(Y[omit,],ncol=1)
    }
    nrow.icv <- nrow(Y.icv)
    Xb.icv <- lapply(1:nvarX, function(j) Xb.lst[[j]][-omit,])
    mu.xb <- lapply(1:nvarX, function(j) apply(Xb.icv[[j]],2,mean))
    Xb.icv.cent <- lapply(1:nvarX, function(j) scale(Xb.icv[[j]], center=mu.xb[[j]], scale=FALSE))
    Xb.valid.cent <- lapply(1:nvarX, function(j) scale(Xb.lst[[j]][omit,], center=mu.xb[[j]], scale=FALSE))

    XbTransInv.icv.lst <- XbTransInv.valid.lst <- Lambda.mtx.lst <- orthConst.mtx.lst <- orthConst.mtx.lst.2 <- list()
    for(j in 1:length.uniqeta){
      normTransInv.mtx <- normTransInv.lst[[j]]
      XbTransInv.icv.lst[[j]] <- do.call(cbind, lapply(1:nvarX, function(jx) crossprod(t(Xb.icv.cent[[jx]]), normTransInv.mtx)))
      XbTransInv.valid.lst[[j]] <- do.call(cbind, lapply(1:nvarX, function(jx) crossprod(t(Xb.valid.cent[[jx]]), normTransInv.mtx)))
      Lambda.mtx.lst[[j]] <- cbind(diag(nrow.icv), -XbTransInv.icv.lst[[j]])
      tmp <- qr(t(Lambda.mtx.lst[[j]]))
      orthConst.mtx.lst[[j]] <- t(qr.Q(tmp))
    }

    for(j in 1:numParm){
      j0 <- match(params.set[j,3], unique(params.set[,3]))
      XbTransInv.icv <- XbTransInv.icv.lst[[j0]]
      XbTransInv.valid <- XbTransInv.valid.lst[[j0]]
      ret <- getComp.riemann.sof.int(XbTransInv.icv, Y.icv, orthConst.mtx.lst[[j0]], s.n.basis, max.K[j], params.set[j,1], params.set[j,2], tol, tol1)  #!!! for the same number of basis
      errors[[j]] <- errors[[j]] + get.cv.errors.hd(XbTransInv.valid, y.valid, y.params, ret)
    }
    print(c("Fold", i, "is completed"))
  }
  min.error <- rep(0,numParm)
  for(j in 1:numParm)
    min.error[j] <- min(errors[[j]],na.rm = TRUE)
  temp <- which.min(min.error)
  opt.tau <- params.set[temp,1]
  opt.mu <- params.set[temp,2]
  opt.eta <- params.set[temp,3]
  tmp.id=which(errors[[temp]]==min.error[temp],  arr.ind =TRUE)
  opt.kappa=kappa[tmp.id[1,1]]
  opt.K=tmp.id[1,2]

  return(list( errors=errors,  min.error=min(min.error, na.rm=TRUE), opt.index=temp, params.set=params.set, opt.K=opt.K, opt.tau=opt.tau, opt.lambda=opt.mu, opt.eta=opt.eta,opt.kappa=opt.kappa, maxK.ret=maxK.ret, t.x.list=t.x.list,t.y=t.y, y.params=y.params))
}

##############################################################
get.cv.errors.hd <- function(xbTransInv.test, y.test,y.params, ret){
  #internal function called by the CV function
  #ret$Beta: Kcomp*bp

  t.y=y.params[[2]]
  y.basis=y.params[[3]]
  kappa=y.params[[4]]
  J.w= t(eval.basis(t.y,y.basis))
  K.w=getbasispenalty(y.basis, 2)
  y.int.weights=y.params[[5]]
  y.weights.aver=mean(diag(y.int.weights))
  q=length(kappa)

  Beta <- t(ret$Beta)
  ncol=dim(Beta)[2]
  T <- as.matrix(ret$XbTransInv %*% Beta)
  tmp.1 <- sqrt(as.numeric(apply(T^2,2,sum)))
  Beta <- scale(Beta, center=FALSE, scale=tmp.1)
  T <- as.matrix(ret$XbTransInv %*% Beta)

  t.mtx <- cbind(1/sqrt(dim(T)[1]),T)

  #print(t(t.mtx)%*%t.mtx)
  T.test=xbTransInv.test %*%  Beta
  t.test.mtx <- cbind(1/sqrt(dim(T)[1]),T.test)

  coef.w.0= J.w%*%y.int.weights%*%t(ret$Y)%*%t.mtx
  error.tmp = matrix(0,q, ncol)
  for(k in 1:q){
    lambda <- kappa[k]
    coef.w <- solve(J.w%*%y.int.weights%*%t(J.w)+lambda*K.w*y.weights.aver)%*%coef.w.0
    for(ncomp in 1:ncol){
      Y.pred <- t.test.mtx[,1:(1+ncomp)]%*%t(coef.w)[1:(1+ncomp),]%*%J.w
      error.tmp[k,ncomp] <- sum(diag((Y.pred-y.test)%*%y.int.weights%*%t(Y.pred-y.test)))
    }
  }

  return(error.tmp)
}



##############################################################

#' @export
pred.hd <- function(fit.cv, X.test, t.y.test=NULL)
{
  if (!is.list(X.test))
  {stop("Error!!: X.test must be a list!")
  }
  if (length(X.test)!=length(fit.cv$t.x.list))
  {stop("Error!!: X.test must be a list of length equal to the number of functional predcitors!")
  }
  if (sum(sapply(1:length(X.test),function(k){!is.matrix(X.test[[k]])})))
  {stop("Error!!: X.test must be a list and all elements must be matrices!")
  }else{
    tmp=fit.cv$t.x.list
    dim.1=sapply(1:length(X.test),function(k){dim(X.test[[k]])[1]})
    dim.2=sum(sapply(1:length(X.test),function(k){dim(X.test[[k]])[2]!=length(tmp[[k]])}))
    if((length(unique(dim.1))!=1)|(dim.2))
    {stop("Error!!: All elements of X.test must be matrices with the same numbers of rows and the number of columns equal to the number of observation time points for the corresponding functional predictor!")
    }
  }

  y.params=fit.cv$y.params
  maxK.ret=fit.cv$maxK.ret
  t.y=y.params[[2]]
  y.basis=y.params[[3]]
  kappa=y.params[[4]]
  J.w= t(eval.basis(t.y,y.basis))
  K.w=getbasispenalty(y.basis, 2)
  y.int.weights=y.params[[5]]
  y.weights.aver=mean(diag(y.int.weights))
  q=length(kappa)

  j=fit.cv$opt.index
  params.set=fit.cv$params.set
  nvarX <- length(X.test)
  eta <- params.set[j,3]
  j0 <- match(eta, unique(params.set[,3]))

  X.test.cent <- lapply(1:nvarX, function(irow) scale(X.test[[irow]], center=maxK.ret$mu.x[[irow]], scale=FALSE))

  xbTransInv.test <- do.call(cbind,lapply(1:nvarX, function(k){X.test.cent[[k]] %*%(maxK.ret$wb[[k]] %*% maxK.ret$normTransInv.lst[[j0]])})) #n*bp
  xbTransInv.test <- as.matrix(xbTransInv.test) / maxK.ret$scale
  Beta <- t(maxK.ret$Beta.max[[j]])
  Beta=Beta[,1:(fit.cv$opt.K)]
  T <- as.matrix(maxK.ret$XbTransInv.lst[[j0]] %*% Beta)
  tmp.1 <- sqrt(as.numeric(apply(T^2,2,sum)))
  Beta <- scale(Beta, center=FALSE, scale=tmp.1)
  T <- as.matrix(maxK.ret$XbTransInv.lst[[j0]]  %*% Beta)

  t.mtx <- cbind(1/sqrt(dim(T)[1]),T)

  #print(t(t.mtx)%*%t.mtx)
  T.test=xbTransInv.test %*%  Beta
  t.test.mtx <- cbind(1/sqrt(dim(T)[1]),T.test)

  coef.w.0= J.w%*%y.int.weights%*%t(maxK.ret$Y)%*%t.mtx
  opt.kappa <- fit.cv$opt.kappa
  coef.w <- solve(J.w%*%y.int.weights%*%t(J.w)+opt.kappa*K.w*y.weights.aver)%*%coef.w.0
  if(is.null(t.y.test)){
     Y.pred <- t.test.mtx%*%t(coef.w)%*%J.w
  }else{
     Y.pred <- t.test.mtx%*%t(coef.w)%*%t(eval.basis(t.y.test,y.basis))
  }
  return(Y.pred)
}



#################################################################

determineMaxKcomp.hd<- function(X, Y, t.x.list, t.y, s.n.basis, params.set,thresh)
{
  repeat.solution=1
  tol=1e-11
  tol1=1e-11
  x.basis<- create.bspline.basis(c(0,1), s.n.basis, 4)

  tmp=lapply(1:length(t.x.list), function(k){seq(0,1,length.out=length(t.x.list[[k]]))})
  t.x.list=tmp
  if(is.vector(Y))
  {Y=matrix(Y,length(Y),1)}

  numParm <- nrow(params.set)
  dims <- dim(Y)
  nobs <- dims[1]
  nvarY <- dims[2]
  nvarX <- length(X)
  K.max <- rep(8,numParm)
  K.comp <- rep(8,numParm)
  nbasis <- x.basis$nbasis

  X.cent <- mu.x <- Xb.lst <- list()
  norms.mtx <- matrix(0,nobs, nvarX)


  J.a <- getbasispenalty(x.basis,0)
  J2.a <- getbasispenalty(x.basis, 2)
  tmp <- x.basis$rangeval
  temp.fun=function(k)
  {
    basis.val <- eval.basis(t.x.list[[k]], x.basis)
    lengthT <- length(t.x.list[[k]])
    tmp.mtx <- diag(lengthT)
    x.int.wt <- (tmp[2]-tmp[1])/lengthT * tmp.mtx
    return(x.int.wt %*% basis.val)
  }
  wb1 <-lapply(1:nvarX, function(k){temp.fun(k)})

  X.cent <- lapply(1:nvarX, function(j) scale(X[[j]], center=TRUE, scale=FALSE))
  mu.x <- lapply(1:nvarX, function(j) attributes(X.cent[[j]])$`scaled:center`)
  Xb.lst <- lapply(1:nvarX, function(j) X.cent[[j]] %*% wb1[[j]])
  tmp <- rep(1,ncol(Xb.lst[[1]]))
  norms.mtx <- sapply(1:nvarX, function(j) sqrt(crossprod(t(Xb.lst[[j]]^2), tmp)))

  Xb.scale <- max(norms.mtx)
  Xb.lst <- lapply(1:nvarX, function(j) Xb.lst[[j]]/Xb.scale)
  total.basis <- nbasis * nvarX     #bp
  Y.cent <- scale(Y, center=T, scale=FALSE)
  muy <- attributes(Y.cent)$`scaled:center`

  zerosobs <- rep(0, nobs)
  obj.vals <- matrix(0,numParm,nvarY)

  XY.prod.lst <- normTrans.lst <- normTransInv.lst <- XbTransInv.lst <- Lambda.mtx.lst <- orthConst.mtx.lst <- orthConst.mtx.lst.2 <- list()
  etas.uniq <- unique(params.set[,3])
  length.uniqeta <- length(etas.uniq)
  for(j in 1:length.uniqeta){
    eta <- etas.uniq[j]
    K.a <- J.a + eta * J2.a      #b*b
    L.tmp <- chol(K.a)        #t(L) %*% L=K.a
    normTrans.lst[[j]] <- L.tmp
    normTransInv.lst[[j]] <- solve(L.tmp)        #bp*bp

    XbTransInv.lst[[j]] <- do.call(cbind, lapply(1:nvarX, function(jx) crossprod(t(Xb.lst[[jx]]), normTransInv.lst[[j]]) ))
    Lambda.mtx.lst[[j]] <- cbind(diag(nobs),-XbTransInv.lst[[j]])
    tmp <- qr(t(Lambda.mtx.lst[[j]]))
    orthConst.mtx.lst[[j]] <- t(qr.Q(tmp))
    orthConst.mtx.lst.2[[j]] <- (orthConst.mtx.lst[[j]])[,-(1:nobs)]
    XY.prod.lst[[j]] <- crossprod(XbTransInv.lst[[j]], Y.cent)
  }

  Beta.max=list()
  print("Calculate the maximum number of components for all tuning parameters")
  for(j in 1:numParm){
    tau <- params.set[j,1]
    mu <- params.set[j,2]
    eta <- params.set[j,3]
    j0 <- match(eta, etas.uniq)
    XbTransInv <- XbTransInv.lst[[j0]]
    orthConst.mtx <- orthConst.mtx.lst[[j0]]
    orthConst.mtx.2 <- orthConst.mtx.lst.2[[j0]]
    XY.prod <- XY.prod.lst[[j0]]

    Beta <- NULL
    for(ncomp in 1:K.max[j]){
      if(!is.null(Beta)){
        tmp <- c(zerosobs, crossprod(XbTransInv %*% beta, XbTransInv))
        tmp <- tmp - crossprod(orthConst.mtx, orthConst.mtx %*% tmp)
        tmp.mtx <-  matrix(tmp/sqrt(sum(tmp^2)), nrow=1)
        orthConst.mtx <- rbind(orthConst.mtx,tmp.mtx)   #(nobs+ncomp-1) %*% (nobs+nbasis*nvarX)
        orthConst.mtx.2 <- rbind(orthConst.mtx.2, tmp.mtx[1,-(1:nobs)])
      }

      nrow.orthconst <- nobs+ncomp-1 # q in paper
      orthConst.arr.2 <- array(orthConst.mtx.2, c(nrow.orthconst, nbasis, nvarX))
      orthConst.lst.2 <- lapply(1:nvarX, function(jx) orthConst.arr.2[,,jx])
      Lambda.all <- lapply(orthConst.lst.2, tcrossprod)


      qq <- rep(0, repeat.solution)
      mqq <- 0
      for(jj in 1:repeat.solution){
        #               pt <- proc.time()
        gamma <- rnorm(nobs+total.basis, 0,1)
        gamma <- gamma/sqrt(gamma^2)
        gamma.old <- rep(0,nobs+total.basis)
        count <- 0
        tmp <- crossprod(gamma[-(1:nobs)], XY.prod)
        value <- sum(tmp^2)
        temp <- tcrossprod(XY.prod, tmp)
        value.old <- 0
        t0 <- rep(0, nrow.orthconst)
        while((count<50)&(min(sum((gamma[-(1:nobs)]-gamma.old[-(1:nobs)])^2), sum((gamma[-(1:nobs)]+gamma.old[-(1:nobs)])^2))>sum(gamma[-(1:nobs)]^2)*1e-10)){
          if(count>0)
          {value.old <- value}
          count <- count+1
          ret <- constOptOrth.H(c(rep(0,nobs),temp), orthConst.mtx, t0, Lambda.all, nobs, nbasis, tau, mu, tol)
          gamma.old=gamma
          t0 <- ret$t
          gamma=ret$x
          tmp <- crossprod(gamma[-(1:nobs)], XY.prod)
          value <- sum(tmp^2)
          temp <- tcrossprod(XY.prod, tmp)

        }

        qq[jj] <- value
        if(qq[jj]>mqq){
          mqq <- qq[jj]
          beta=gamma[-(1:nobs)]
        }
      }

      Beta=rbind(Beta, beta)     #K*(bp)
      obj.vals[j, ncomp] <- mqq
      if(mqq/sum(obj.vals[j,])<thresh){
        K.comp[j] <- ncomp
        break
      }

    }#end for(ncomp)
    Beta.max[[j]]=Beta
  }#end for(j in 1:numParm)
  return(list(Beta.max=Beta.max, Y=Y, obj.vals=obj.vals, max.K=K.comp, Xb.lst=Xb.lst, Y.cent=Y.cent, mu.x=mu.x, mu.y=muy, wb=wb1, normTrans.lst=normTrans.lst, normTransInv.lst=normTransInv.lst, Lambda.mtx.lst=Lambda.mtx.lst, XbTransInv.lst=XbTransInv.lst, orthConst.mtx.lst=orthConst.mtx.lst, orthConst.mtx.lst.2=orthConst.mtx.lst.2, XY.prod.lst=XY.prod.lst, total.basis=total.basis, scale=Xb.scale))
}


######################################
#' @export
getcoef.hd=function(fit.cv)
{
  maxK.ret=fit.cv$maxK.ret
  k=length(fit.cv$t.x.list)
  m.y=length(fit.cv$t.y)
  m.list=sapply(1:k, function(j){length(fit.cv$t.x.list[[j]])})
  range.list=sapply(1:k, function(j){max(fit.cv$t.x.list[[j]])-min(fit.cv$t.x.list[[j]])})
  Beta=maxK.ret$Beta.max[[fit.cv$opt.index]]
  nbasis=dim(Beta)[2]/k
  nonzero=0
  for(i in 1:dim(Beta)[1])
  {
    D=matrix(Beta[i,],nbasis,k)
    L=apply(D^2,2,sum)
    nonzero=nonzero+1*(L!=0)
  }
  nonzero.id= which(nonzero!=0)

  X.test.0=lapply(1:k, function(j){matrix(0, 1, m.list[[j]])}) #T*S
  tmp.1=lapply(1:k, function(j){as.matrix(X.test.0[[j]], nrow=1)})
  mu=as.vector(pred.hd(X.test=tmp.1, fit.cv))

  beta=lapply(1:k, function(j){matrix(0, m.list[j], m.y)})
  for(i in 1:length(nonzero.id))
  {
    X.test.0=lapply(1:k, function(j){matrix(0, m.list[[nonzero.id[i]]], m.list[[j]])})
    X.test=X.test.0
    X.test[[nonzero.id[i]]]=diag(m.list[[nonzero.id[i]]])
    beta[[nonzero.id[i]]]=(pred.hd(fit.cv, X.test)-pred.hd(fit.cv, X.test.0))*m.list[[nonzero.id[i]]]/range.list[[nonzero.id[i]]]
  }


  return(list(mu=mu, beta=beta))

}

