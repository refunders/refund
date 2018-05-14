face.sparse <- function(data, newdata = NULL,
                        center=TRUE,argvals.new=NULL,
                        knots=7, 
                        # knots.option="equally-spaced",
                        p=3,m=2,lambda=NULL,lambda_mean=NULL,
                        search.length=14,
                        lower=-3,upper=10, lower2 = -3, upper2=5,
                        calculate.scores=FALSE,pve=0.99,two_step=FALSE){

  if(!two_step){
    fit <- face.sparse.inner(data = data, newdata=newdata, W = NULL,
                             argvals.new = argvals.new,
                             center=center,
                             knots=knots,
                             # knots.option = knots.option,
                             p = p, m = m, lambda= lambda, lambda_mean = lambda_mean,
                             lower=lower,upper=upper,search.length=search.length,
                             calculate.scores=calculate.scores, pve = pve)
    
  }else{
    fit <- face.sparse.inner(data = data, newdata=NULL, W = NULL,
                             argvals.new = argvals.new,
                             center=center,
                             knots=knots,
                             # knots.option = knots.option,
                             p = p, m = m, lambda= lambda, lambda_mean = lambda_mean,
                             lower=lower,upper=upper,search.length=search.length,
                             pve = pve)
    
    W <- weight.construct(fit,data)
    
    
    fit <- face.sparse.inner(data = data, newdata=newdata, W = W,
                             argvals.new = argvals.new,
                             center=center,
                             knots=knots,
                             # knots.option = knots.option,
                             p = p, m = m, lambda= lambda, lambda_mean = lambda_mean,
                             lower=lower2,upper=upper2,search.length=length(lower2:upper2),
                             calculate.scores=calculate.scores,pve=pve)
  }
  

  return(fit)
}
  
