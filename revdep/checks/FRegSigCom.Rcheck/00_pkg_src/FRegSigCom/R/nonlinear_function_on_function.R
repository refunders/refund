#' @import fda
#library(fda)

#########################################

eval.tensor.spline.basis=function(t.x, x, delta.x, bspline.x.obj, bspline.s.obj)
{

  lenghtX=length(t.x)
  x=as.vector(t(x))

  A=t(eval.basis(t.x,bspline.s.obj))
  B=t(eval.basis(x,bspline.x.obj))

  tmp=sapply(1:(length(x)/lenghtX), function(k){A%*%(diag(delta.x)%*%t(B[,(k-1)*lenghtX+1:lenghtX]))})


  return(tmp)
}


#########################################

eval.penalty.R.matrix=function(bspline.x.obj, bspline.s.obj)
{
  J=bspline.x.obj$nbasis
  L=bspline.s.obj$nbasis
  R=R.inv=matrix(0, J*L,J*L)

  K.x0=getbasispenalty(bspline.x.obj, 0)
  K.x1=getbasispenalty(bspline.x.obj, 1)
  K.x2=getbasispenalty(bspline.x.obj, 2)

  K.s0=getbasispenalty(bspline.s.obj, 0)
  K.s1=getbasispenalty(bspline.s.obj, 1)
  K.s2=getbasispenalty(bspline.s.obj, 2)

  K=K.x0%x%K.s0+K.x2%x%K.s0+K.x1%x%K.s1+K.x0%x%K.s2
  R=chol(K)
  R.inv=solve(R)

  return(list(R=R,R.inv=R.inv))
}



#########################################

power.eig=function(Pi, M.basis, tol=1e-12)
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

determine.max.comp=function(t.x.list, X, Y, lambda.set,  y.params, bspline.x.obj, bspline.s.obj, upper.com, thresh)
{
  n.curves=length(X)
  t.x.list.org=t.x.list
  t.x.list=lapply(1:n.curves, function(k){seq(0,1,length.out=length(t.x.list[[k]]))})

  delta.x.list=lapply(1:n.curves, function(k){rep(1/length(t.x.list[[k]]),length(t.x.list[[k]]))})
  shift.x.list=list()
  scale.x.list=list()


  for(i in 1:n.curves)
  {
    shift.x.list[[i]]=min(X[[i]])
    scale.x.list[[i]]=max(X[[i]])-min(X[[i]])
    X[[i]]=(X[[i]]-shift.x.list[[i]])/scale.x.list[[i]]
  }
  delta.y=y.params[[5]]
  z.list=list()
  T.list=list()
  G.list.no.center.no.scale=list()
  G.list=list()
  n.sample=dim(Y)[1]
  for(i in 1:n.curves)
  {
    G.list.no.center.no.scale[[i]]=t(eval.tensor.spline.basis(t.x.list[[i]], X[[i]], delta.x.list[[i]], bspline.x.obj, bspline.s.obj))
    tmp=G.list.no.center.no.scale[[i]]
    G=scale(tmp, scale=FALSE)
    G=G/sqrt(n.sample)
    G.list[[i]]=t(G)
  }
  Y.cent=scale(Y, scale=FALSE)
  Pi=Y.cent%*%(diag(delta.y)%*%t(Y.cent))/n.sample
  Pi=(Pi+t(Pi))/2
  tmp=eval.penalty.R.matrix(bspline.x.obj, bspline.s.obj)
  R=tmp$R
  Rinv=tmp$R.inv
  J=dim(R)[1]
  R.inv.tran.G.list=list()
  for(i in 1:n.curves)
  {
    R.inv.tran.G.list[[i]]=t(Rinv)%*%G.list[[i]]
  }
  R.inv.tran.G=do.call(rbind, R.inv.tran.G.list)
  max.comp=rep(0, length(lambda.set))
  max.value.list=list()

  for(lam.ind in 1:length(lambda.set))
  {
    lambda=lambda.set[lam.ind]
    M=rbind(diag(n.sample)*sqrt(lambda), -R.inv.tran.G)
    tmp <- qr(M)
    M.basis <-  qr.Q(tmp)
    # print(range(t(M.basis)%*%M.basis-diag(n.sample)))
    # print(range(M-M.basis%*%(t(M.basis)%*%M)))
    beta=NULL
    max.value=rep(0, upper.com)

    for(n.comp in 1:upper.com)
    {
      tmp=power.eig(Pi, M.basis)
      max.value[n.comp]=tmp$max.value
      beta=cbind(beta,tmp$beta)
      if(max.value[n.comp]/sum(max.value[1:n.comp])<thresh)
      {break}
      #print(range(beta[1:n.sample,n.comp]*sqrt(lambda)-t(G)%*%(Rinv%*%beta[-(1:n.sample),n.comp])))
      D=rep(0, dim(beta)[1])
      D[-(1:n.sample)]=R.inv.tran.G%*%(t(R.inv.tran.G)%*%beta[-(1:n.sample), n.comp])
      tmp=D-M.basis%*%(t(M.basis)%*%D)
      M.basis=cbind(M.basis, tmp/sqrt(sum(tmp^2)))
      # print(range(t(M.basis)%*%M.basis-diag(dim(M.basis)[2])))
      # print(range(D-M.basis%*%(t(M.basis)%*%D)))
    }
    max.comp[lam.ind]=n.comp
    max.value.list[[lam.ind]]=max.value[1:n.comp]
    tmp=lapply(1:n.curves, function(i){lambda^(-1/2)*Rinv%*%beta[n.sample+(i-1)*J+1:J,]})
    z.list[[lam.ind]]=do.call(rbind,tmp)
    T.list[[lam.ind]]=beta[1:n.sample,]
    # print(range(beta[(1:n.sample), ]-lambda^(-0.5)*t(R.inv.tran.G)%*%beta[-(1:n.sample), ]))
    # print(range(T.list[[lam.ind]]-t(G)%*%z.list[[lam.ind]]))
    # print(t(T.list[[lam.ind]])%*%T.list[[lam.ind]])
    # print(eigen(Pi)$values[1:n.comp])
    # print(max.value.list[[lam.ind]])
  }
  return(list(t.x.list.org=t.x.list.org, t.x.list=t.x.list, R=R,delta.x.list=delta.x.list, bspline.x.obj=bspline.x.obj, bspline.s.obj=bspline.s.obj, lambda.set=lambda.set, Rinv=Rinv, Y=Y, shift.x.list=shift.x.list, scale.x.list=scale.x.list, max.comp= max.comp, max.value.list=max.value.list, z.list=z.list, T.list=T.list, G.list.no.center.no.scale=G.list.no.center.no.scale))
}



#########################################


cv.folds <- function(n,nfolds=5)
  ## Randomly split the n samples into folds
  ## Returns a list of nfolds lists of indices, each corresponding to a fold
{
  return(split(sample(n),rep(1:nfolds,length=n)))
}


#########################################
#' @export
cv.nonlinear=function(X, Y, t.x.list, t.y, s.n.basis = 30, x.n.basis=30, t.n.basis = 30, K.cv=5, upper.com=15, thresh=0.005)
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
  {stop("Error!!: The number of columns of each component of X must be equal to the length of the corresponsing component of  t.x.list!")
  }
  if(dim(Y)[2]!=length(t.y))
  {stop("Error!!: the number of columns of Y must be equal to the length of the vector t.y of the observation points!")
  }
  tol=1e-12
  bspline.s.obj=create.bspline.basis(c(0,1), s.n.basis, 4)
  bspline.x.obj=create.bspline.basis(c(0, 1), x.n.basis, 4)
  lambda.set=c(1e-10,1e-8,1e-6, 1e-4,1e-2,1,100)
  t.y=seq(0,1,length.out=length(t.y))
  kappa=c(1e-10,1e-8,1e-6,1e-4,1e-2,1,100) #values for the tuning parameter kappa
  y.basis<- create.bspline.basis(c(0,1), t.n.basis, 4)
  y.params=list()
  y.params[[1]]=c(0,1)
  y.params[[2]]=t.y
  y.params[[3]]=y.basis
  y.params[[4]]=kappa
  y.params[[5]]= rep(1/length(t.y),length(t.y))


  fit.1=determine.max.comp(t.x.list, X, Y, lambda.set, y.params, bspline.x.obj, bspline.s.obj, upper.com=upper.com, thresh=thresh)

  delta.y=y.params[[5]]
  G.list.no.center.no.scale=fit.1$G.list.no.center.no.scale
  n.curves=length(G.list.no.center.no.scale)
  lambda.set=fit.1$lambda.set
  R=fit.1$R
  Rinv=fit.1$Rinv
  Y.all=fit.1$Y

  max.comp=fit.1$max.comp
  all.folds <- cv.folds(dim(Y.all)[1], K.cv)

  errors <- list()
  for(j in 1: length(max.comp))
  { errors[[j]] <- matrix(0, length(y.params[[4]]), max.comp[j])}
  print("Starting the cross-validation procedure:")
  for(i in 1:K.cv)
  {
    print(paste0("fold ",i))
    omit <- all.folds[[i]]
    Y=Y.all[-omit,]
    Y.valid=Y.all[omit,]
    n.sample=dim(Y)[1]
    G.list=list()
    G.valid.list=list()
    for(k in 1:n.curves)
    {
      tmp=G.list.no.center.no.scale[[k]]
      G=scale(tmp[-omit,] , scale=FALSE)
      G.mean=attr(G,  "scaled:center")
      G.valid=scale(tmp[omit,], scale=FALSE,center=G.mean)
      G=G/sqrt(n.sample)
      G.list[[k]]=t(G)
      G.valid=G.valid/sqrt(n.sample)
      G.valid.list[[k]]=t(G.valid)
    }


    Y.cent=scale(Y, scale=FALSE)
    Pi=Y.cent%*%(diag(delta.y)%*%t(Y.cent))/n.sample
    Pi=(Pi+t(Pi))/2
    R.inv.tran.G.list=list()
    for(k in 1:n.curves)
    {
      R.inv.tran.G.list[[k]]=t(Rinv)%*%G.list[[k]]
    }
    R.inv.tran.G=do.call(rbind, R.inv.tran.G.list)


    for(lam.ind in 1:length(lambda.set))
    {
      lambda=lambda.set[lam.ind]
      M=rbind(diag(n.sample)*sqrt(lambda), -R.inv.tran.G)
      tmp <- qr(M)
      M.basis <-  qr.Q(tmp)
      # print(range(t(M.basis)%*%M.basis-diag(n.sample)))
      # print(range(M-M.basis%*%(t(M.basis)%*%M)))
      beta=NULL

      for(n.comp in 1:max.comp[lam.ind])
      {
        tmp=power.eig(Pi, M.basis)
        beta=cbind(beta,tmp$beta)
        D=rep(0, dim(beta)[1])
        D[-(1:n.sample)]=R.inv.tran.G%*%(t(R.inv.tran.G)%*%beta[-(1:n.sample), n.comp])
        tmp=D-M.basis%*%(t(M.basis)%*%D)
        M.basis=cbind(M.basis, tmp/sqrt(sum(tmp^2)))
        # print(range(t(M.basis)%*%M.basis-diag(dim(M.basis)[2])))
        # print(range(D-M.basis%*%(t(M.basis)%*%D)))
      }
      T=beta[1:n.sample,]
      R.inv.tran.G.valid.list=list()
      for(i in 1:n.curves)
      {
        R.inv.tran.G.valid.list[[i]]=t(Rinv)%*%G.valid.list[[i]]
      }
      R.inv.tran.G.valid=do.call(rbind, R.inv.tran.G.valid.list)

      # print(range(T-lambda^(-0.5)*t(R.inv.tran.G)%*%beta[-(1:n.sample), ]))
      T.valid=lambda^(-0.5)*t(R.inv.tran.G.valid)%*%beta[-(1:n.sample), ]
      tmp.1 <- sqrt(as.numeric(apply(T^2,2,sum)))
      T <- scale(T, center=FALSE, scale=tmp.1)
      T.valid <- scale(T.valid, center=FALSE, scale=tmp.1)
      errors[[lam.ind]] <- errors[[lam.ind]] + get.cv.errors.nonlinear(T, T.valid,Y, Y.valid, y.params)
    }
  }
  min.error <- rep(0,length(lambda.set))
  for(j in 1:length(lambda.set)){
      min.error[j] <- min(errors[[j]],na.rm = TRUE)}
  temp <- which.min(min.error)
  opt.lambda <- lambda.set[temp]
  tmp.id=which(errors[[temp]]==min.error[temp],  arr.ind =TRUE)
#  kappa=y.params[[4]]
  opt.kappa=kappa[tmp.id[1,1]]
  opt.K=tmp.id[1,2]
  return(list(opt.index=temp, min.error=min(min.error, na.rm=TRUE),  errors=errors, lambda.set=lambda.set, kappa.set=kappa, opt.K=opt.K, opt.lambda=opt.lambda, opt.kappa=opt.kappa, fit.1=fit.1, y.params=y.params))
}

######################################
get.cv.errors.nonlinear <- function(T, T.test,Y, y.test, y.params){
  #internal function called by the CV function

  t.y=y.params[[2]]
  y.basis=y.params[[3]]
  kappa=y.params[[4]]
  J.w= t(eval.basis(t.y,y.basis))
  K.w=getbasispenalty(y.basis, 2)
  y.int.weights=diag(y.params[[5]])
  y.weights.aver=mean(diag(y.int.weights))
  q=length(kappa)

  ncol=dim(T)[2]

  t.mtx <- cbind(1/sqrt(dim(T)[1]),T)

  # print(t(t.mtx)%*%t.mtx)
  t.test.mtx <- cbind(1/sqrt(dim(T)[1]),T.test)

  coef.w.0= J.w%*%y.int.weights%*%t(Y)%*%t.mtx
  error.tmp = matrix(0,q, ncol)
  for(k in 1:q){
    coef.w <- solve(J.w%*%y.int.weights%*%t(J.w)+kappa[k]*K.w*y.weights.aver)%*%coef.w.0
    for(ncomp in 1:ncol){
      Y.pred <- t.test.mtx[,1:(1+ncomp)]%*%t(coef.w)[1:(1+ncomp),]%*%J.w
      error.tmp[k,ncomp] <- sum(diag((Y.pred-y.test)%*%y.int.weights%*%t(Y.pred-y.test)))
      #error.tmp[k,ncomp] <- mean(diag((Y.pred-y.test)%*%y.int.weights%*%t(Y.pred-y.test)))
    }
  }

  return(error.tmp)
}


######################################
#' @export
pred.nonlinear <- function(fit.cv, X.test, t.y.test=NULL){
  fit.1=fit.cv$fit.1
  t.x.list=fit.1$t.x.list
  y.params=fit.cv$y.params
  n.curves=length(X.test)
  bspline.x.obj=fit.1$bspline.x.obj
  bspline.s.obj=fit.1$bspline.s.obj
  delta.x.list=fit.1$delta.x.list
  shift.x.list=fit.1$shift.x.list
  scale.x.list=fit.1$scale.x.list
  for(k in 1:n.curves)
  {

    X.test[[k]]=(X.test[[k]]-shift.x.list[[k]])/scale.x.list[[k]]
    X.test[[k]]=((X.test[[k]]+1)-abs(X.test[[k]]-1))/2
    X.test[[k]]=((X.test[[k]])+abs(X.test[[k]]))/2
  }
  Y=fit.1$Y
  n.sample=dim(Y)[1]

  opt.index=fit.cv$opt.index
  opt.K=fit.cv$opt.K
  opt.lambda=fit.cv$opt.lambda
  kappa=fit.cv$opt.kappa

  t.y=y.params[[2]]
  y.basis=y.params[[3]]


  J.w= t(eval.basis(t.y,y.basis))
  K.w=getbasispenalty(y.basis, 2)
  y.int.weights=diag(y.params[[5]])
  y.weights.aver=mean(diag(y.int.weights))

  T=fit.1$T.list[[opt.index]]
  T=T[,1:opt.K]
  z=fit.1$z.list[[opt.index]]
  z=z[,1:opt.K]

  ncol=opt.K

  G.list.no.center.no.scale=fit.1$G.list.no.center.no.scale
  G.test.list=list()
  for(k in 1:n.curves)
  {
    G.center=apply(G.list.no.center.no.scale[[k]],2,mean)
    tmp=t(eval.tensor.spline.basis(t.x.list[[k]], X.test[[k]], delta.x.list[[k]], bspline.x.obj, bspline.s.obj))
    G.test.list[[k]]= scale(tmp,center=G.center,scale=FALSE)/sqrt(n.sample)
  }
  G.test=do.call(cbind,G.test.list)
  T.test=G.test%*%z

  T=as.matrix(T)
  T.test=as.matrix(T.test)
  tmp.1 <- sqrt(as.numeric(apply(T^2,2,sum)))
  T <- scale(T, center=FALSE, scale=tmp.1)
  T.test <- scale(T.test, center=FALSE, scale=tmp.1)


  t.mtx <- cbind(1/sqrt(dim(T)[1]),T)

  # print(t(t.mtx)%*%t.mtx)
  t.test.mtx <- cbind(1/sqrt(dim(T)[1]),T.test)

  coef.w.0= J.w%*%y.int.weights%*%t(Y)%*%t.mtx
  coef.w <- solve(J.w%*%y.int.weights%*%t(J.w)+kappa*K.w*y.weights.aver)%*%coef.w.0
  if(is.null(t.y.test)){
    Y.pred <- t.test.mtx%*%t(coef.w)%*%J.w
  }else{
    Y.pred <- t.test.mtx%*%t(coef.w)%*%t(eval.basis(t.y.test,y.basis))
  }
  return(Y.pred=Y.pred)
}

##########################################
#' @export
getcoef.nonlinear=function(fit.cv, n.x.grid=100)
{

  fit.1=fit.cv$fit.1
  t.x.list.org=fit.1$t.x.list.org
  t.x.list=fit.1$t.x.list
  n.curves=length(t.x.list.org)
  range.t.x.list=lapply(1:n.curves, function(k){max(t.x.list.org[[k]])-min(t.x.list.org[[k]])})
  X.grid=list()
  for(k in 1:n.curves)
  {
    a=(fit.1)$scale.x.list[[k]]
    b=(fit.1)$shift.x.list[[k]]
    X.grid[[k]]=seq(b,a+b, length.out=n.x.grid)
  }
  t.y=fit.cv$y.params[[2]]

  tmp.list=lapply(1:n.curves, function(k){matrix(0, 1,length(t.x.list[[k]]))})
  Y.pred.0=pred.nonlinear(fit.cv, tmp.list)
  mu=as.vector(Y.pred.0)

  F=list()
  X.new=lapply(1:n.curves, function(k){A=diag(length(t.x.list[[k]])); X.grid[[k]]%x%A})
    for(k in 1:n.curves)
  {
    X.0=lapply(1:n.curves, function(j){matrix(0, dim(X.new[[k]])[1], dim(X.new[[j]])[2])})
    X.0[[k]]=X.new[[k]]
    Y.pred=pred.nonlinear(fit.cv, X.0)
    Y.pred=t(sapply(1:dim(Y.pred)[1],function(k){Y.pred[k,]}))
    Y.pred=t(sapply(1:dim(Y.pred)[1],function(k){Y.pred[k,]-Y.pred.0}))
    F[[k]]=array(Y.pred, c(length(t.x.list[[k]]),length(X.grid[[k]]),length(t.y)))
    F[[k]]=aperm(F[[k]], c(2,1,3))*length(t.x.list[[k]])/(range.t.x.list[[k]])
  }
  return(list(mu=mu, F=F, X.grid=X.grid, t.x.list=t.x.list.org, t.y=t.y))
}



