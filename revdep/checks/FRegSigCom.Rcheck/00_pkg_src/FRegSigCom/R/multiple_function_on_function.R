#' @import fda
#' @importFrom stats rnorm
#require(fda)

#########################################

power.eig=function(Pi, M.basis, tol=1e-16)
{
  n=dim(Pi)[1]
  p=dim(M.basis)[1]
  tmp=rnorm(p)
  b=tmp/sqrt(sum(tmp^2))
  b.old=rep(0,p)
  count=0
  while((sum((b-b.old)^2)>tol)&(sum((b+b.old)^2)>tol))
  {
    count=count+1
    # print(c("count",count))
    if(count>50)
    {break;}
    #print(c(sum((b-b.old)^2), sum((b+b.old)^2)))
    b.old=b
    tmp=b-M.basis%*%(t(M.basis)%*%b)
    tmp1=rep(0,p)
    tmp1[1:n]=Pi%*%tmp[1:n]

    b=tmp1-M.basis%*%(t(M.basis)%*%tmp1)
    max.value=t(b.old)%*%b
    #print(c("maximum",max.value))
    b=b/sqrt(sum(b^2))
  }
  # print(c("count",count))
  return(list(max.value=max.value, beta=as.vector(b)))
}




#########################################

determine.max.comp.smooth=function(X,  Y, x.smooth.params, y.smooth.params, Z=NULL, tol=1e-12, thresh=0.01,upper.comp=10)
{
  n.curves=length(X)
  x.time.length=x.smooth.params[[3]]
  y.time.length=y.smooth.params[[3]]
  n.sample=dim(Y)[1]

  z.list=list()
  T.list=list()


  B.vals=x.smooth.params[[4]]
  tmp=lapply(1:length(B.vals), function(i){B.vals[[i]]%*%t(X[[i]])/x.time.length[i]})
  G.no.center.no.scale= Reduce(rbind,tmp)
  Q=dim(G.no.center.no.scale)[1]

  if(!is.null(Z))
  {G.no.center.no.scale=rbind(G.no.center.no.scale, t(Z))

  tmp=G.no.center.no.scale
  G=t(scale(t(tmp), scale=FALSE))

  G=G/sqrt(n.sample)
  G.1=G[1:Q,]
  G.2=G[-(1:Q),]
  Y.cent=scale(Y, scale=FALSE)
  Pi=tcrossprod(Y.cent)/n.sample/y.time.length

  Pi=(Pi+t(Pi))/2

  K=x.smooth.params[[5]]
  R= chol(K)
  R.inv=solve(R)
  B=diag(n.curves)
  R=B%x%R
  R.inv=B%x%R.inv
  R.inv.tran.G=t(R.inv)%*%G.1
  lambda.set=x.smooth.params[[2]]
  tuning.length=length(lambda.set)
  max.comp=rep(0, tuning.length)
  max.value.list=list()
  Z.1=t(G[-(1:Q),])
  tmp=svd(Z.1,nu=dim(Z.1)[1])
  U.2=tmp$u[,-(1:dim(Z.1)[2])]

  for(tuning.ind in 1:tuning.length)
  {
    lambda=lambda.set[tuning.ind]
    M=rbind(diag(n.sample)*sqrt(lambda), -R.inv.tran.G)
    M=tcrossprod(M, t(U.2))
    tmp <- qr(M)
    M.basis <-  qr.Q(tmp)
    # print(range(t(M.basis)%*%M.basis-diag(dim(M.basis)[2])))
    # print(range(M-M.basis%*%(t(M.basis)%*%M)))
    beta=NULL
    max.value=rep(0, upper.comp)

    for(n.comp in 1:upper.comp)
    {
      tmp=power.eig(Pi, M.basis)
      max.value[n.comp]=tmp$max.value
      beta=cbind(beta,tmp$beta)
      if(max.value[n.comp]/sum(max.value[1:n.comp])<thresh)
      {break}
      #print(range(beta[1:n.sample,n.comp]*sqrt(lambda)-t(R.inv.tran.G)%*%beta[-(1:n.sample),n.comp]))
      D=rep(0, dim(beta)[1])
      D[(1:n.sample)]=beta[1:n.sample,n.comp]
      D=D/sqrt(sum(D^2))
      tmp=D-M.basis%*%(t(M.basis)%*%D)
      M.basis=cbind(M.basis, tmp/sqrt(sum(tmp^2)))
      # print(range(t(M.basis)%*%M.basis-diag(dim(M.basis)[2])))
      # print(range(D-M.basis%*%(t(M.basis)%*%D)))
    }
    beta=as.matrix(beta)
    w=beta[1:n.sample,]
    v=beta[-(1:n.sample),]
    z.1=lambda^(-1/2)*R.inv%*%beta[-(1:n.sample),]
    w.2=w-crossprod(G.1,z.1)
    z.2=solve(tcrossprod(G.2), crossprod(t(G.2),w.2))
    # print(range(w-crossprod(G.1,z.1)-crossprod(G.2,z.2)))

    max.comp[tuning.ind]=n.comp
    max.value.list[[tuning.ind]]=max.value[1:n.comp]
    z.list[[tuning.ind]]= rbind(z.1,z.2)
    T.list[[tuning.ind]]=w
    # print(range(T.list[[tuning.ind]]-t(G)%*%z.list[[tuning.ind]]))
    # print(t(T.list[[tuning.ind]])%*%T.list[[tuning.ind]])
    # print(eigen(Pi)$values[1:n.comp])
    # print(max.value.list[[tuning.ind]])
  }
  }else{
    tmp=G.no.center.no.scale
    G=t(scale(t(tmp), scale=FALSE))

    G=G/sqrt(n.sample)

    Y.cent=scale(Y, scale=FALSE)
    Pi=tcrossprod(Y.cent)/n.sample/y.time.length

    Pi=(Pi+t(Pi))/2

    K=x.smooth.params[[5]]
    R= chol(K)
    R.inv=solve(R)
    B=diag(n.curves)
    R=B%x%R
    R.inv=B%x%R.inv
    R.inv.tran.G=t(R.inv)%*%G
    lambda.set=x.smooth.params[[2]]
    tuning.length=length(lambda.set)
    max.comp=rep(0, tuning.length)
    max.value.list=list()


    for(tuning.ind in 1:tuning.length)
    {
      lambda=lambda.set[tuning.ind]
      M=rbind(diag(n.sample)*sqrt(lambda), -R.inv.tran.G)
      tmp <- qr(M)
      M.basis <-  qr.Q(tmp)
      # print(range(t(M.basis)%*%M.basis-diag(dim(M.basis)[2])))
      # print(range(M-M.basis%*%(t(M.basis)%*%M)))
      beta=NULL
      max.value=rep(0, upper.comp)

      for(n.comp in 1:upper.comp)
      {
        tmp=power.eig(Pi, M.basis)
        max.value[n.comp]=tmp$max.value
        beta=cbind(beta,tmp$beta)
        if(max.value[n.comp]/sum(max.value[1:n.comp])<thresh)
        {break}
        #print(range(beta[1:n.sample,n.comp]*sqrt(lambda)-t(R.inv.tran.G)%*%beta[-(1:n.sample),n.comp]))
        D=rep(0, dim(beta)[1])
        D[(1:n.sample)]=beta[1:n.sample,n.comp]
        D=D/sqrt(sum(D^2))
        tmp=D-M.basis%*%(t(M.basis)%*%D)
        M.basis=cbind(M.basis, tmp/sqrt(sum(tmp^2)))
        # print(range(t(M.basis)%*%M.basis-diag(dim(M.basis)[2])))
        # print(range(D-M.basis%*%(t(M.basis)%*%D)))
      }
      beta=as.matrix(beta)
      w=beta[1:n.sample,]
      v=beta[-(1:n.sample),]
      z=lambda^(-1/2)*R.inv%*%beta[-(1:n.sample),]
      # print(range(w-crossprod(G,z)))

      max.comp[tuning.ind]=n.comp
      max.value.list[[tuning.ind]]=max.value[1:n.comp]
      z.list[[tuning.ind]]= z
      T.list[[tuning.ind]]=w
      # print(range(T.list[[tuning.ind]]-t(G)%*%z.list[[tuning.ind]]))
      # print(t(T.list[[tuning.ind]])%*%T.list[[tuning.ind]])
      # print(eigen(Pi)$values[1:n.comp])
      # print(max.value.list[[tuning.ind]])
    }
  }


  return(list(B.vals=B.vals, R=R, R.inv=R.inv, Q=Q,  Y=Y,  max.comp= max.comp, max.value.list=max.value.list, z.list=z.list, T.list=T.list, G.no.center.no.scale=G.no.center.no.scale))
}



#########################################


cv.folds <- function(n,nfolds=5)
  ## Randomly split the n samples into folds
  ## Returns a list of nfolds lists of indices, each corresponding to a fold
{
  return(split(sample(n),rep(1:nfolds,length=n)))
}


#########################################

cv.for.lambda=function(t.x, t.y, X, Y, Z=NULL, tau, s.n.basis, t.n.basis, all.folds, tol,upper.comp, thresh)
{

  K.cv=length(all.folds)
  tmp.1=sapply(1:length(t.x), function(i){length(t.x[[i]])})
  tmp.2=sapply(1:length(X), function(i){dim(X[[i]])[2]})
  x.smooth.params=list()
  x.smooth.params[[1]] =create.bspline.basis(c(0,1), s.n.basis)
  if(length(t.x)<=20){
     x.smooth.params[[2]]=c(1e-10, 1e-7, 1e-4,  1e-1)   #tuning parameter for the whole penalty on psi(s)
    #x.smooth.params[[2]]=c(1e-9, 1e-6, 1e-3,  1) #paper
  }else{
      x.smooth.params[[2]]=c(1,10, 10^3, 10^5, 10^7)
  }
  x.smooth.params[[3]]= tmp.1
  x.smooth.params[[4]]= lapply(1:length(tmp.1), function(i){t(eval.basis(seq(0,1,length=tmp.1[i]), x.smooth.params[[1]]))})

  tmp=getbasispenalty(x.smooth.params[[1]],0)+tau*getbasispenalty(x.smooth.params[[1]],2)
  x.smooth.params[[5]]=  (tmp+t(tmp))/2

  y.smooth.params=list()
  y.smooth.params[[1]]=create.bspline.basis(c(0,1), t.n.basis)
  #if(length(t.x)<10){
    y.smooth.params[[2]]= c(1e-8,1e-7,1e-6,1e-5,1e-4, 1e-3,1e-2,1,10)
  #}else{
   # y.smooth.params[[2]]= c(1e-10,1e-8,1e-6,1e-4, 1e-2,1)
  #}
  y.smooth.params[[3]]=  length(t.y)
  y.smooth.params[[4]]= t(eval.basis(seq(0,1,length=length(t.y)), y.smooth.params[[1]]))
  tmp=getbasispenalty(y.smooth.params[[1]],2)
  y.smooth.params[[5]]=  (tmp+t(tmp))/2

  fit.1=determine.max.comp.smooth(X,  Y, x.smooth.params, y.smooth.params,Z, upper.comp=upper.comp,thresh=thresh)
  # print("maximum number of components for lambdas:")
  # print(fit.1$max.comp)

  Q=fit.1$Q
  G.no.center.no.scale=fit.1$G.no.center.no.scale
  lambda.set=x.smooth.params[[2]]
  tuning.length=length(lambda.set)
  y.time.length=y.smooth.params[[3]]


  kappa.set=y.smooth.params[[2]]
  B.vals= y.smooth.params[[4]]
  K.w=y.smooth.params[[5]]
  y.weights.aver=1/y.smooth.params[[3]]
  q=length(kappa.set)
  B.vals.weig=B.vals*y.weights.aver
  y.penalty.inv=list()
  for(k in 1:q)
  {
    kappa=kappa.set[k]
    y.penalty.inv[[k]]=solve(B.vals.weig%*%t(B.vals)+kappa*K.w*y.weights.aver)
  }

  R=fit.1$R
  R.inv=fit.1$R.inv
  Y.all=fit.1$Y

  max.comp=fit.1$max.comp

  errors <- list()
  for(j in 1: length(max.comp))
  { errors[[j]] <- matrix(0, length(kappa.set), max.comp[j])}
  if(!is.null(Z))
  {
    for(fold.ind in 1:K.cv)
    {
      print(paste0("fold ",fold.ind))
      omit <- all.folds[[fold.ind]]
      Y=Y.all[-omit,]
      Y.valid=Y.all[omit,]
      tmp=t(G.no.center.no.scale[,-omit])
      G=scale(tmp, scale=FALSE)

      G.mean=attr(G,  "scaled:center")
      tmp=t(G.no.center.no.scale[,omit])
      G.valid=scale(tmp, scale=FALSE,center=G.mean)
      n.sample=dim(G)[1]

      G=G/sqrt(n.sample)
      G=t(G)
      G.1=G[1:Q,]
      G.2=G[-(1:Q),]

      G.valid=G.valid/sqrt(n.sample)
      G.valid=t(G.valid)
      G.valid.1=G.valid[1:Q,]
      G.valid.2=G.valid[-(1:Q),]
      Y.cent=scale(Y, scale=FALSE)
      t0 <- proc.time()[3]

      Pi=tcrossprod(Y.cent)/n.sample/y.time.length

      Pi=(Pi+t(Pi))/2

      R.inv.tran.G=t(R.inv)%*%G.1

      R.inv.tran.G.valid=t(R.inv)%*%G.valid.1

      Z.1=t(G[-(1:Q),])
      tmp=svd(Z.1,nu=dim(Z.1)[1])
      U.2=tmp$u[,-(1:dim(Z.1)[2])]


      for(tuning.ind in 1:tuning.length)
      {
        lambda=lambda.set[tuning.ind]
        M=rbind(diag(n.sample)*sqrt(lambda), -R.inv.tran.G)
        M=tcrossprod(M, t(U.2))
        tmp <- qr(M)
        M.basis <-  qr.Q(tmp)
        # print(range(t(M.basis)%*%M.basis-diag(dim(M.basis)[2])))
        # print(range(M-M.basis%*%(t(M.basis)%*%M)))
        beta=NULL

        for(n.comp in 1:max.comp[tuning.ind])
        {
          tmp=power.eig(Pi, M.basis)

          beta=cbind(beta,tmp$beta)
          D=rep(0, dim(beta)[1])
          D[(1:n.sample)]=beta[1:n.sample,n.comp]
          D=D/sqrt(sum(D^2))
          tmp=D-M.basis%*%(t(M.basis)%*%D)

          M.basis=cbind(M.basis, tmp/sqrt(sum(tmp^2)))
          # print(range(t(M.basis)%*%M.basis-diag(dim(M.basis)[2])))
          # print(range(D-M.basis%*%(t(M.basis)%*%D)))
        }
        beta=as.matrix(beta)
        w=beta[1:n.sample,]
        v=beta[-(1:n.sample),]
        z.1=lambda^(-1/2)*R.inv%*%beta[-(1:n.sample),]
        w.2=w-crossprod(G.1,z.1)
        z.2=solve(tcrossprod(G.2), crossprod(t(G.2),w.2))

        z=rbind(z.1,z.2)
        # print(range(w-crossprod(G,z)))


        T=w
        T.valid=crossprod(G.valid,z)
        tmp.1 <- sqrt(as.numeric(apply(T^2,2,sum)))
        T <- scale(T, center=FALSE, scale=tmp.1)
        T.valid <- scale(T.valid, center=FALSE, scale=tmp.1)

        errors[[tuning.ind]] <- errors[[tuning.ind]] + get.cv.errors.smooth(T, T.valid,Y, Y.valid, y.smooth.params,y.penalty.inv, B.vals.weig)

      }
    }
  }else{

    for(fold.ind in 1:K.cv)
    {
      print(paste0("fold ",fold.ind))
      omit <- all.folds[[fold.ind]]
      Y=Y.all[-omit,]
      Y.valid=Y.all[omit,]
      tmp=t(G.no.center.no.scale[,-omit])
      G=scale(tmp, scale=FALSE)

      G.mean=attr(G,  "scaled:center")
      tmp=t(G.no.center.no.scale[,omit])
      G.valid=scale(tmp, scale=FALSE,center=G.mean)
      n.sample=dim(G)[1]

      G=G/sqrt(n.sample)
      G=t(G)


      G.valid=G.valid/sqrt(n.sample)
      G.valid=t(G.valid)

      Y.cent=scale(Y, scale=FALSE)
      t0 <- proc.time()[3]

      Pi=tcrossprod(Y.cent)/n.sample/y.time.length

      Pi=(Pi+t(Pi))/2

      R.inv.tran.G=t(R.inv)%*%G

      R.inv.tran.G.valid=t(R.inv)%*%G.valid



      for(tuning.ind in 1:tuning.length)
      {
        lambda=lambda.set[tuning.ind]
        M=rbind(diag(n.sample)*sqrt(lambda), -R.inv.tran.G)
        tmp <- qr(M)
        M.basis <-  qr.Q(tmp)
        # print(range(t(M.basis)%*%M.basis-diag(dim(M.basis)[2])))
        # print(range(M-M.basis%*%(t(M.basis)%*%M)))
        beta=NULL

        for(n.comp in 1:max.comp[tuning.ind])
        {
          tmp=power.eig(Pi, M.basis)

          beta=cbind(beta,tmp$beta)
          D=rep(0, dim(beta)[1])
          D[(1:n.sample)]=beta[1:n.sample,n.comp]
          D=D/sqrt(sum(D^2))
          tmp=D-M.basis%*%(t(M.basis)%*%D)

          M.basis=cbind(M.basis, tmp/sqrt(sum(tmp^2)))
          # print(range(t(M.basis)%*%M.basis-diag(dim(M.basis)[2])))
          # print(range(D-M.basis%*%(t(M.basis)%*%D)))
        }
        beta=as.matrix(beta)
        w=beta[1:n.sample,]
        v=beta[-(1:n.sample),]
        z=lambda^(-1/2)*R.inv%*%beta[-(1:n.sample),]
        # print(range(w-crossprod(G,z)))
        T=w
        T.valid=crossprod(G.valid,z)
        tmp.1 <- sqrt(as.numeric(apply(T^2,2,sum)))
        T <- scale(T, center=FALSE, scale=tmp.1)
        T.valid <- scale(T.valid, center=FALSE, scale=tmp.1)
        errors[[tuning.ind]] <- errors[[tuning.ind]] + get.cv.errors.smooth(T, T.valid,Y, Y.valid, y.smooth.params,y.penalty.inv, B.vals.weig)

      }
    }

  }#end of if(!is.null(Z))else

  min.error <- rep(0,tuning.length)
  for(j in 1:tuning.length)
  {min.error[j] <- min(errors[[j]],na.rm = TRUE)}
  temp <- which.min(min.error)
  opt.lambda=lambda.set[temp]

  tmp.id=which(errors[[temp]]==min.error[temp],  arr.ind =TRUE)
  opt.kappa=kappa.set[tmp.id[1,1]]
  opt.K=tmp.id[1,2]
  return(list(opt.index=temp, opt.lambda=opt.lambda, opt.kappa=opt.kappa, min.error=min(min.error, na.rm=TRUE), errors=errors, opt.K=opt.K, x.smooth.params=x.smooth.params,  y.smooth.params=y.smooth.params, fit.1=fit.1, is.null.Z=is.null(Z)))
}

######################################
#' @export
cv.sigcom=function(X, Y, t.x, t.y, Z=NULL, s.n.basis=50, t.n.basis=50, K.cv=5, upper.comp=10, thresh=0.005, tol=1e-12)
{
  if(!is.list(X))
  {stop("Error!!: X must be a list!")}
  if (sum(sapply(1:length(X),function(k){!is.matrix(X[[k]])})))
  {stop("Error!!: X must be a list and all its components must be matrix!")
  }
  if(!is.list(t.x))
  {stop("Error!!: t.x must be a list!")}
  if (length(X)!=length(t.x))
  {stop("Error!!: both X and t.x must be lists and they have the same numbers of  components!")
  }
  dim.1=sapply(1:length(X),function(k){dim(X[[k]])[1]})
  if((length(unique(dim.1))!=1))
  {stop("Error!!: all components of X must be matrix and have the same numbers of rows!")
  }
  if((dim(X[[1]])[1]!=dim(Y)[1]))
  {stop("Error!!: the number of observations of X (that is, the number of rows of each component of X) must be equal to the number of observations of Y (that is, the number of rows of Y)!")
  }
  if(sum(sapply(1:length(X), function(k){dim(X[[k]])[2]!=length(t.x[[k]])}))!=0)
  {stop("Error!!: The number of columns of each component of X must be equal to the length of the corresponsing component of  t.x!")
  }
  if(dim(Y)[2]!=length(t.y))
  {stop("Error!!: the number of columns of Y must be equal to the length of the vector t.y of the observation points!")
  }

  tau.set=c(1e-4,1e-2,1, 10)     #tuning parameter for psi''(s)
  #tau.set=c(1e-2, 1e-1,1, 1e1)
  length.tau=length(tau.set)
  min.error=1e20
  opt.tau=tau.set[1]
  all.folds <- cv.folds(nrow(Y)[1], K.cv)

  for(i in 1:length.tau)
  {
    tau=tau.set[i]
    # print(t(c("Cross validation for tau=", tau, " and all lambda")))
    print(paste0("Cross validation for tau=", tau, " and all lambda"))
    tmp.cv=cv.for.lambda(t.x, t.y, X, Y, Z, tau, s.n.basis, t.n.basis, all.folds, tol, upper.comp, thresh=thresh)
    if(tmp.cv$min.error<min.error)
    {fit.cv=tmp.cv
    min.error=tmp.cv$min.error
    opt.tau=tau
    }
  }


  return(list(t.x.list=t.x, t.y=t.y, Z=Z, opt.index=fit.cv$opt.index, opt.lambda=fit.cv$opt.lambda, opt.tau=opt.tau,  opt.kappa=fit.cv$opt.kappa, opt.K=fit.cv$opt.K, min.error=fit.cv$min.error, errors=fit.cv$errors,  x.smooth.params=fit.cv$x.smooth.params,  y.smooth.params=fit.cv$y.smooth.params, fit.1=fit.cv$fit.1, is.null.Z=fit.cv$is.null.Z))

}




######################################
get.cv.errors.smooth<- function(T, T.test,Y, Y.test, y.smooth.params, y.penalty.inv, B.vals.weig){
  #internal function called by the CV function
  #ret$Beta: Kcomp*bp

  y.length=y.smooth.params[[3]]

  kappa.set=y.smooth.params[[2]]
  B.vals= y.smooth.params[[4]]
  q=length(kappa.set)

  ncol=dim(T)[2]

  t.mtx <- cbind(1/sqrt(dim(T)[1]),T)
  t.test.mtx <- cbind(1/sqrt(dim(T)[1]),T.test)

  coef.w.0= B.vals.weig%*%t(Y)%*%t.mtx
  error.tmp = matrix(0,q, ncol)
  for(k in 1:q)
  {
    coef.w <- y.penalty.inv[[k]]%*%coef.w.0
    V=crossprod(B.vals, coef.w)
    for(ncomp in 1:ncol)
    {
      Y.pred<- tcrossprod(t.test.mtx[,1:(1+ncomp)], V[, 1:(1+ncomp)])
      error.tmp[k,ncomp] <- sum((Y.pred-Y.test)^2/y.length)
    }
  }

  return(error.tmp)
}

######################################
#' @export
pred.sigcom <- function(fit.cv, X.test, Z.test=NULL, t.y.test=NULL){
  #internal function called by the CV function
  #ret$Beta: Kcomp*bp
  if (!is.list(X.test))
  {stop("Error!!: X.test must be a list!")
  }
  if (length(X.test)!=length(fit.cv$t.x.list))
  {stop("Error!!: X.test must be a list and the number of its components must be equal to the number of predcitor curves!")
  }
  if (sum(sapply(1:length(X.test),function(k){!is.matrix(X.test[[k]])})))
  {stop("Error!!: X.test must be a list and all its components must be matrix!")
  }else{
    tmp=fit.cv$t.x.list
    dim.1=sapply(1:length(X.test),function(k){dim(X.test[[k]])[1]})
    dim.2=sum(sapply(1:length(X.test),function(k){dim(X.test[[k]])[2]!=length(tmp[[k]])}))
    if((length(unique(dim.1))!=1)|(dim.2))
    {stop("Error!!: all components of X.test must be matrix, their numbers of rows are the same and the number of columns must be equal to the number of observation points of the corresponding predictor curve!")
    }
  }

  fit.1=fit.cv$fit.1
  Y=fit.1$Y
  n.sample=dim(Y)[1]
  opt.index=fit.cv$opt.index
  opt.K=fit.cv$opt.K
  lambda=fit.cv$opt.lambda
  kappa=fit.cv$opt.kappa

  x.time.length=fit.cv$x.smooth.params[[3]]

  B.vals.x=fit.cv$x.smooth.params[[4]]
#  y.smooth.params[[4]]= t(eval.basis(seq(0,1,length=length(t.y)), y.smooth.params[[1]]))
  if(is.null(t.y.test)){
    B.vals.y=fit.cv$y.smooth.params[[4]]
    y.weights.aver=1/fit.cv$y.smooth.params[[3]]
  }else{
    B.vals.y=eval.basis(t.y.test, y.smooth.params[[1]])
    y.weights.aver=1/length(t.y.test)
  }
  B.vals.weig.y=B.vals.y * y.weights.aver

  T=fit.1$T.list[[opt.index]]
  T=T[,1:opt.K]
  z=fit.1$z.list[[opt.index]]
  z=z[,1:opt.K]

  ncol=opt.K

  G.no.center.no.scale=fit.1$G.no.center.no.scale
  if(!fit.cv$is.null.Z)
  {
    if (!is.matrix(Z.test))
    {stop("Error!!: Z.test must be a matrix!")
    }
    if (dim(X.test[[1]])[1]!=dim(Z.test)[1])
    {stop("Error!!:  the number of the rows of each component in X.test must be equal to the number of rows of Z.test!")
    }
    if (dim(Z.test)[2]!=dim(fit.cv$Z)[2])
    {stop("Error!!:  the number of columns of Z.test must equal to the number of scalar predcitors!")
    }
    G.center=apply(G.no.center.no.scale,1,mean)
    tmp=lapply(1:length(B.vals.x), function(i){B.vals.x[[i]]%*%t(X.test[[i]])/x.time.length[i]})
    tmp.1=Reduce(rbind, tmp)
    tmp.1=rbind(tmp.1,t(Z.test))
    G.test=scale(t(tmp.1),center=G.center,scale=FALSE)/sqrt(n.sample)
    T.test=G.test%*%z
    T=as.matrix(T)
    T.test=as.matrix(T.test)
    tmp.1 <- sqrt(as.numeric(apply(T^2,2,sum)))
    T <- scale(T, center=FALSE, scale=tmp.1)
    T.test <- scale(T.test, center=FALSE, scale=tmp.1)
    t.mtx <- cbind(1/sqrt(dim(T)[1]),T)
    # print(t(t.mtx)%*%t.mtx)
    t.test.mtx <- cbind(1/sqrt(dim(T)[1]),T.test)
    y.smooth.params=fit.cv$y.smooth.params
    y.time.length=y.smooth.params[[3]]
    kappa.set=y.smooth.params[[2]]
    #B.vals= y.smooth.params[[4]]
    K.w=y.smooth.params[[5]]
    #y.weights.aver=1/y.smooth.params[[3]]
    #B.vals.weig=B.vals*y.weights.aver
    coef.w.0= B.vals.weig.y%*%t(Y)%*%t.mtx
    coef.w <- solve(B.vals.weig.y%*%t(B.vals.y)+kappa*K.w*y.weights.aver)%*%coef.w.0
    Y.pred <- t.test.mtx%*%(t(coef.w)%*%B.vals.y)

  }else
  {
    G.center=apply(G.no.center.no.scale,1,mean)
    tmp=lapply(1:length(B.vals.x), function(i){B.vals.x[[i]]%*%t(X.test[[i]])/x.time.length[i]})
    tmp.1=Reduce(rbind, tmp)
    G.test=scale(t(tmp.1),center=G.center,scale=FALSE)/sqrt(n.sample)
    T.test=G.test%*%z
    T=as.matrix(T)
    T.test=as.matrix(T.test)
    tmp.1 <- sqrt(as.numeric(apply(T^2,2,sum)))
    T <- scale(T, center=FALSE, scale=tmp.1)
    T.test <- scale(T.test, center=FALSE, scale=tmp.1)
    t.mtx <- cbind(1/sqrt(dim(T)[1]),T)
    # print(t(t.mtx)%*%t.mtx)
    t.test.mtx <- cbind(1/sqrt(dim(T)[1]),T.test)
    y.smooth.params=fit.cv$y.smooth.params
    y.time.length=y.smooth.params[[3]]
    kappa.set=y.smooth.params[[2]]
    #B.vals= y.smooth.params[[4]]
    K.w=y.smooth.params[[5]]
    #y.weights.aver=1/y.smooth.params[[3]]
    #B.vals.weig=B.vals*y.weights.aver
    coef.w.0= B.vals.weig.y%*%t(Y)%*%t.mtx
    coef.w <- solve(B.vals.weig.y%*%t(B.vals.y)+kappa*K.w*y.weights.aver)%*%coef.w.0
    Y.pred <- t.test.mtx%*%(t(coef.w)%*%B.vals.y)
  }

  return(Y.pred=Y.pred)
}


######################################
#' @export
getcoef.sigcom=function(fit.cv)
{

  k=length(fit.cv$t.x.list)
#  m.y=length(fit.cv$t.y)
  m.list=sapply(1:k, function(j){length(fit.cv$t.x.list[[j]])})
  range.list=sapply(1:k, function(j){max(fit.cv$t.x.list[[j]])-min(fit.cv$t.x.list[[j]])})
  X.test.0=lapply(1:k, function(j){matrix(0, 1, m.list[[j]])})

  beta=list()
  beta.scalar=list()
  if(fit.cv$is.null.Z)
  {
    Z.test=NULL
    Z.test.0=NULL
    tmp.1=lapply(1:k, function(j){as.matrix(X.test.0[[j]], nrow=1)})
    mu=as.vector(pred.sigcom(fit.cv, X.test=tmp.1))
    for(i in 1:k)
    {
      X.test.0=lapply(1:k, function(j){matrix(0, m.list[[i]], m.list[[j]])})
      X.test=X.test.0
      X.test[[i]]=diag(m.list[[i]])
      beta[[i]]=(pred.sigcom(fit.cv, X.test)-pred.sigcom(fit.cv, X.test.0))*m.list[[i]]/range.list[[i]]
    }
  }else{
    Z.test=diag(dim(fit.cv$Z)[2])
    Z.test.0=0*Z.test
    tmp.1=lapply(1:k, function(j){t(as.matrix(X.test.0[[j]][1,]))})
    tmp.2=t(as.matrix(Z.test.0[1,]))
    mu=as.vector(pred.sigcom(fit.cv, X.test=tmp.1, Z.test=tmp.2))
    Z.test.00=matrix(0, m.list, dim(Z.test)[2])
    X.test.00=lapply(1:k, function(j){matrix(0, dim(Z.test)[2], m.list[[j]])})
    for(i in 1:k)
    {
      X.test.0=lapply(1:k, function(j){matrix(0, m.list[[i]], m.list[[j]])})
      X.test=X.test.0
      X.test[[i]]=diag(m.list[[i]])
      Z.test.00=matrix(0, dim(X.test[[i]])[1], dim(Z.test)[2])
      beta[[i]]=(pred.sigcom(fit.cv, X.test,Z.test.00)-pred.sigcom(fit.cv, X.test.0,Z.test.00))*m.list[[i]]/range.list[[i]]
      beta.scalar=pred.sigcom(fit.cv, X.test.00,Z.test)-pred.sigcom(fit.cv, X.test.00,Z.test.0)
    }
  }
  return(list(mu=mu, beta=beta, beta.scalar=beta.scalar))

}
