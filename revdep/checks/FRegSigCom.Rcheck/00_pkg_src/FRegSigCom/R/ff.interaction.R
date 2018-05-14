#' @import fda
#' @import Matrix
#' @importFrom Rcpp evalCpp
#' @useDynLib FRegSigCom
#' @name FRegSigCom

#library(Rcpp)
#library(RcppEigen)
#library(fda)
#library(Matrix)

cv.folds <- function(n,nfolds=5)
  ## Randomly split the n samples into folds
  ## Returns a list of nfolds lists of indices, each corresponding to a fold
{
  return(split(sample(n),rep(1:nfolds,length=n)))
}


#########################################

###################################################################

#######################################################################
#' @export
step.ff.interaction=function(X, Y, t.x, t.y, s.n.basis=40, t.n.basis=40, inter.n.basis=20, basis.type="Bspline", K.cv=5, all.folds=NULL,  upper.comp=8, thresh=0.01)
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

  if(is.null(all.folds))
  {
    all.folds <- cv.folds(dim(Y)[1], K.cv)
  }
  K.cv=length(all.folds)
  n.sample=dim(Y)[1]
  x.smooth.params=list()
  n.curves=length(X)
  x.smooth.params[[1]]=n.curves

  if(basis.type=="Bspline")
  {
    basis.obj.x =create.bspline.basis(c(0,1), s.n.basis)
  }
  if(basis.type=="Fourier")
  {
    if(s.n.basis%%2==0)
    {
      s.n.basis=s.n.basis+1
      print("In **create.fourier.basis(c(0, 1), s.n.basis)** s.n.basis must be an odd integer; since s.n.basis is even now, it will be increased by 1")
    }
    basis.obj.x =create.fourier.basis(c(0,1), s.n.basis)
  }
  lambda.set=c(1e-8,1e-6, 1e-4, 1e-2, 1, 1e2)
  x.smooth.params[[2]]=lambda.set
  x.smooth.params[[3]]=n.sample
  tmp.1=sapply(1:length(t.x), function(i){length(t.x[[i]])})

  x.smooth.params[[4]]= lapply(1:length(tmp.1), function(i){t(eval.basis(seq(0,1,length=tmp.1[i]),    basis.obj.x))})



  tmp=list()

  tmp[[1]]=getbasispenalty(basis.obj.x,0)
  tmp[[2]]=getbasispenalty(basis.obj.x,2)

  x.smooth.params[[5]]= tmp
  tmp.1=sapply(1:length(t.x), function(i){length(t.x[[i]])})
  if(basis.type=="Bspline")
  {
    basis.obj=create.bspline.basis(c(0,1), inter.n.basis)
  }
  if(basis.type=="Fourier")
  {
    if(inter.n.basis%%2==0)
    {
      inter.n.basis=inter.n.basis+1
      print("In create.fourier.basis(c(0, 1), inter.n.basis) nbasis must be an odd integer; since inter.n.basisis is even now, it will be increased by 1")
    }
    basis.obj=create.fourier.basis(c(0,1), inter.n.basis)
  }

  x.smooth.params[[6]]=c(s.n.basis, inter.n.basis^2)
  x.smooth.params[[7]]= lapply(1:length(tmp.1), function(i){t(eval.basis(seq(0,1,length=tmp.1[i]), basis.obj))})

  tmp.3=list()

  tmp=getbasispenalty(basis.obj,0)
  tmp.0=(tmp+t(tmp))/2
  tmp=getbasispenalty(basis.obj,1)
  tmp.1=(tmp+t(tmp))/2
  tmp=getbasispenalty(basis.obj,2)
  tmp.2=(tmp+t(tmp))/2
  tmp.3[[1]]=tmp.0%x%tmp.0
  tmp.3[[2]]=tmp.2%x%tmp.0
  tmp.3[[3]]=tmp.1%x%tmp.1
  tmp.3[[4]]=tmp.0%x%tmp.2
  x.smooth.params[[8]]=  tmp.3

  tau.set=c(1e-3,1e-1,1e1, 1e3)
  x.smooth.params[[9]]= tau.set
  x.smooth.params[[10]]=t.x

  y.smooth.params = list()
  y.smooth.params[[1]] = create.bspline.basis(c(0, 1), t.n.basis)
  y.smooth.params[[2]] = c(1e-11, 1e-9,1e-7,1e-5,1e-3,1e-1,1e1, 1e3)
  y.smooth.params[[3]] =  length(t.y)
  y.smooth.params[[4]] = t(eval.basis(seq(0, 1, length = length(t.y)), y.smooth.params[[1]]))
  tmp = getbasispenalty(y.smooth.params[[1]], 2)
  y.smooth.params[[5]] =  (tmp + t(tmp)) / 2
  B.vals = y.smooth.params[[4]]
  K.w = y.smooth.params[[5]]
  y.weights.aver = 1 / y.smooth.params[[3]]

  B.vals.weig = B.vals * y.weights.aver
  y.penalty.inv = list()
  kappa.set = y.smooth.params[[2]]
  tmp=list()
  tmp[[1]]=B.vals.weig %*% t(B.vals)
  tmp[[2]]=K.w *y.weights.aver

  y.smooth.params[[6]] = B.vals.weig
  y.smooth.params[[7]] = tmp



  x.raw.params=x.smooth.params
  x.raw.params[[2]]=c(1e-8,1e-4, 1)
  x.raw.params[[9]]=c(1e-2,1,1e2)


  fit.step.c=C_stepwise_adaptive(t.x,  X, Y, x.raw.params, x.smooth.params,y.smooth.params, all.folds, upper.comp, thresh)
  tmp=(1:n.curves)
  opt.main.effects=tmp[fit.step.c$opt_main_index==1]
  opt.interaction.effects=NULL
  if(sum(fit.step.c$opt_inter_index)>0)
  {
    opt.interaction.effects=fit.step.c$inter_mat[fit.step.c$opt_inter_index==1,]+1
    if(sum(fit.step.c$opt_inter_index)==1)
    {
      opt.interaction.effects=matrix(opt.interaction.effects, 1, 2)
    }

  }
  fit.cv=C_cv_fix_effects(t.x,  X, Y, fit.step.c$opt_main_index, fit.step.c$opt_inter_index,  x.raw.params, x.smooth.params,y.smooth.params, all.folds, upper.comp, thresh)

  for(k in 1:length(x.smooth.params[[8]]))
  {
    tmp=x.smooth.params[[8]][[k]]
    x.smooth.params[[8]][[k]]=Matrix(tmp, sparse=TRUE)
  }


  return(list(opt.main.effects=opt.main.effects, opt.interaction.effects=opt.interaction.effects, fitted_model=fit.cv$fit_cv_fix_effects,  y_penalty_inv=fit.cv$y_penalty_inv, X=X, Y=Y, x.smooth.params=x.smooth.params, y.smooth.params=y.smooth.params))
}



#######################################################################
#' @export
cv.ff.interaction=function( X, Y, t.x, t.y, main.effect, interaction.effect=NULL, s.n.basis=40, t.n.basis=40, inter.n.basis=20, basis.type="Bspline",  K.cv=5, all.folds=NULL,  upper.comp=8, thresh=0.01)
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

  n.curves=length(X)
  if(sum(main.effect%in%(1:n.curves))!=length(main.effect))
  {
    stop("Error!!: the index in main.effect is not correct!")
  }

  if(!is.null(interaction.effect))
  {
    if(is.vector(interaction.effect))
    {
      if(length(interaction.effect)!=2)
      {
        stop("Error!!: interaction.effect must be a matrix with two columns or a vector of length 2!")
      }
      interaction.effect=matrix(interaction.effect, 1, 2)
    }

    if(is.matrix(interaction.effect))
    {
      if(dim(interaction.effect)[2]!=2)
      {
        stop("Error!!: interaction.effect must be a matrix with two columns or a vector of length 2!")
      }
      for(i in 1:nrow(interaction.effect))
      {
        if(sum(interaction.effect[i,]%in%(1:n.curves))!=length(interaction.effect[i,]))
        {
          stop("Error!!: the index in interaction.effect is not correct!")
        }
      }
    }else
    {
      stop("Error!!: interaction.effect must be a matrix with two columns or a vector of length 2!")
    }
  }

  if(is.null(all.folds))
  {
    all.folds <- cv.folds(dim(Y)[1], K.cv)
  }
  K.cv=length(all.folds)
  n.sample=dim(Y)[1]
  x.smooth.params=list()

  x.smooth.params[[1]]=n.curves

  if(basis.type=="Bspline")
  {basis.obj.x =create.bspline.basis(c(0,1), s.n.basis)
  }
  if(basis.type=="Fourier")
  {
    if(s.n.basis%%2==0)
    {
      s.n.basis=s.n.basis+1
      print("In **create.fourier.basis(c(0, 1), s.n.basis)** s.n.basis must be an odd integer; since s.n.basis is even now, it will be increased by 1")
    }
    basis.obj.x =create.fourier.basis(c(0,1), s.n.basis)
  }
  lambda.set=c(1e-8,1e-6, 1e-4, 1e-2, 1, 1e2)
  x.smooth.params[[2]]=lambda.set
  x.smooth.params[[3]]=n.sample
  tmp.1=sapply(1:length(t.x), function(i){length(t.x[[i]])})

  x.smooth.params[[4]]= lapply(1:length(tmp.1), function(i){t(eval.basis(seq(0,1,length=tmp.1[i]),    basis.obj.x))})



  tmp=list()

  tmp[[1]]=getbasispenalty(basis.obj.x,0)
  tmp[[2]]=getbasispenalty(basis.obj.x,2)

  x.smooth.params[[5]]= tmp
  tmp.1=sapply(1:length(t.x), function(i){length(t.x[[i]])})
  if(basis.type=="Bspline")
  {
    basis.obj=create.bspline.basis(c(0,1), inter.n.basis)
  }
  if(basis.type=="Fourier")
  {
    if(inter.n.basis%%2==0)
    {
      inter.n.basis=inter.n.basis+1
      print("In create.fourier.basis(c(0, 1), inter.n.basis) nbasis must be an odd integer; since inter.n.basisis is even now, it will be increased by 1")
    }
    basis.obj=create.fourier.basis(c(0,1), inter.n.basis)
  }

  x.smooth.params[[6]]=c(s.n.basis, inter.n.basis^2)
  x.smooth.params[[7]]= lapply(1:length(tmp.1), function(i){t(eval.basis(seq(0,1,length=tmp.1[i]), basis.obj))})

  tmp.3=list()

  tmp=getbasispenalty(basis.obj,0)
  tmp.0=(tmp+t(tmp))/2
  tmp=getbasispenalty(basis.obj,1)
  tmp.1=(tmp+t(tmp))/2
  tmp=getbasispenalty(basis.obj,2)
  tmp.2=(tmp+t(tmp))/2
  tmp.3[[1]]=tmp.0%x%tmp.0
  tmp.3[[2]]=tmp.2%x%tmp.0
  tmp.3[[3]]=tmp.1%x%tmp.1
  tmp.3[[4]]=tmp.0%x%tmp.2

  x.smooth.params[[8]]=  tmp.3
  tau.set=c(1e-3,1e-1,1e1, 1e3)
  x.smooth.params[[9]]= tau.set
  x.smooth.params[[10]]=t.x

  y.smooth.params = list()
  y.smooth.params[[1]] = create.bspline.basis(c(0, 1), t.n.basis)
  y.smooth.params[[2]] = c(1e-11, 1e-9,1e-7,1e-5,1e-3,1e-1,1e1, 1e3)
  y.smooth.params[[3]] =  length(t.y)
  y.smooth.params[[4]] = t(eval.basis(seq(0, 1, length = length(t.y)), y.smooth.params[[1]]))
  tmp = getbasispenalty(y.smooth.params[[1]], 2)
  y.smooth.params[[5]] =  (tmp + t(tmp)) / 2
  B.vals = y.smooth.params[[4]]
  K.w = y.smooth.params[[5]]
  y.weights.aver = 1 / y.smooth.params[[3]]

  B.vals.weig = B.vals * y.weights.aver
  y.penalty.inv = list()
  kappa.set = y.smooth.params[[2]]
  tmp=list()
  tmp[[1]]=B.vals.weig %*% t(B.vals)
  tmp[[2]]=K.w *y.weights.aver

  y.smooth.params[[6]] = B.vals.weig
  y.smooth.params[[7]] = tmp

  main.index=rep(0,n.curves)
  for(i in 1:length(main.effect))
  {
    main.index[main.effect[i]]=1
  }
  inter.index=rep(0, n.curves+n.curves*(n.curves-1)/2)
  if(!is.null(interaction.effect))
  {
    interaction.matrix=matrix(0, n.curves+n.curves*(n.curves-1)/2,2);
    k=1
    for(i in 1:n.curves)
    {
      for(j in i:n.curves)
      {
        interaction.matrix[k,1]=i;
        interaction.matrix[k,2]=j;
        k=k+1;
      }
    }
    for(i in 1:nrow(interaction.effect))
    {
      a=interaction.effect[i,1]
      b=interaction.effect[i,2]
      if(a>b)
      {
        tmp=a
        a=b
        b=tmp
      }
      for(j in 1:nrow(interaction.matrix))
      {
        if((interaction.matrix[j,1]==a)&(interaction.matrix[j,2]==b))
        {
          inter.index[j]=1
        }
      }
    }
  }


  x.raw.params=x.smooth.params
  x.raw.params[[2]]=c(1e-8,1e-4, 1)
  x.raw.params[[9]]=c(1e-2,1,1e2)

  fit.cv=C_cv_fix_effects(t.x,   X, Y, main.index, inter.index,  x.raw.params, x.smooth.params,y.smooth.params, all.folds, upper.comp, thresh)

  for(k in 1:length(x.smooth.params[[8]]))
  {
    tmp=x.smooth.params[[8]][[k]]
    x.smooth.params[[8]][[k]]=Matrix(tmp, sparse=TRUE)
  }
  return(list(fitted_model=fit.cv$fit_cv_fix_effects,  y_penalty_inv=fit.cv$y_penalty_inv, X=X,  Y=Y, x.smooth.params=x.smooth.params, y.smooth.params=y.smooth.params))
}



######################################
#' @export
pred.ff.interaction <- function(fit.obj,  X.test){
  fit.cv=fit.obj$fitted_model
  y_penalty_inv=fit.obj$y_penalty_inv
  x.smooth.params=fit.obj$x.smooth.params
  for(k in 1:length(x.smooth.params[[8]]))
  {
    tmp=x.smooth.params[[8]][[k]]
    x.smooth.params[[8]][[k]]=matrix(tmp, nrow(tmp), ncol(tmp))
  }
  y.smooth.params=fit.obj$y.smooth.params
  Y.train=fit.obj$Y
  Y.pred=C_pred_ff_inter(fit.cv, Y.train,  X.test, x.smooth.params, y.smooth.params, y_penalty_inv)

  return(Y.pred=Y.pred)
}


######################################
#' @export
getcoef.ff.interaction <- function(fit.obj){
  fit.cv=fit.obj$fitted_model
  x.smooth.params=fit.obj$x.smooth.params
  t.x=x.smooth.params[[10]]
  for(k in 1:length(x.smooth.params[[8]]))
  {
    tmp=x.smooth.params[[8]][[k]]
    x.smooth.params[[8]][[k]]=matrix(tmp, nrow(tmp), ncol(tmp))
  }
  y.smooth.params=fit.obj$y.smooth.params
  y_penalty_inv=fit.obj$y_penalty_inv
  Y.train=fit.obj$Y
  X.train=fit.obj$X
  coef.fit=C_find_coef_ff_interaction(fit.cv, X.train, Y.train, x.smooth.params, y.smooth.params, y_penalty_inv)
  intercept=coef.fit$intercept
  coef_main=coef.fit$coef_main
  for(i in 1:length(coef_main))
  {
    coef_main[[i]]=coef_main[[i]]/(max(t.x[[i]])-min(t.x[[i]]))
  }
  main_effects=coef.fit$main_effects+1
  inter_effects=coef.fit$inter_effects+1
  coef.inter.list=coef.fit$coef_inter
  coef_inter=list()
  if(length(coef.inter.list)>0)
  {
   for(i in 1:length(coef.inter.list))
   {
     coef_inter[[i]]=array(unlist(coef.inter.list[[i]]), c(dim(coef.inter.list[[i]][[1]]), length(coef.inter.list[[i]])))
     tmp=(max(t.x[[inter_effects[i,1]]])-min(t.x[[inter_effects[i,1]]]))*(max(t.x[[inter_effects[i,2]]])-min(t.x[[inter_effects[i,2]]]))
     coef_inter[[i]]=coef_inter[[i]]/tmp
   }
  }
  return(list(intercept=intercept, main_effects=main_effects, coef_main=coef_main, inter_effects=inter_effects, coef_inter=coef_inter))
}


