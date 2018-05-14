## About: 
## This function is to make subject-specific prediction based
## on an existing fit of class "face_sparse".

## Function arguments:

## "fit": an object of class "face_sparse"
## "Y.pred": a data frame with three arguments:
# (1) "argvals": observation times;
# (2) "subj": subject indices;
# (3) "y": values of observations;
## NA values are allowed in "y" but not in the other two.

predict.face.sparse <- function(object,newdata,...){
  
  ## check inputs
  if(class(object)!="face.sparse") stop("'fit' has to be a face.sparse object")
  check.data(newdata,type="predict") 

  subj.pred = newdata$subj
  subj_unique.pred = unique(subj.pred)
  y.pred = newdata$y
  
  Theta = object$Theta
  knots = object$knots
  npc = object$npc
  p = object$p
  
  center = object$center
  #const.var.error = object$const.var.error
  fit_mean = object$fit_mean
  #fit_var_error = object$fit_var_error
  sigma2 = object$sigma2
  
  mu.pred <- rep(0,nrow(newdata))
  var.error.pred <- rep(max(sigma2,0.000001),nrow(newdata))
  cov.pred <-  matrix(0,length(mu.pred),length(mu.pred))
    
  if(center) {mu.pred <- predict.pspline.face(fit_mean,newdata$argvals) }
  #if(!const.var.error){
  #  var.error.pred <- predict.pspline(fit_var_error,newdata$argvals)
  #  var.error.pred <-  sapply(var.error.pred,function(x) max(x,0))### added Sept 8, 2015 by Luo
  #}

  
  Bnew = spline.des(knots=knots, x=object$argvals.new, ord = p+1,outer.ok = TRUE,sparse=TRUE)$design
  scores = list(subj=subj_unique.pred,
                scores = matrix(NA,nrow=length(subj_unique.pred),ncol=npc),
                u = matrix(NA,nrow=length(subj_unique.pred),ncol=nrow(Theta))
                )
  
  for(i in 1:length(subj_unique.pred)){
    sel.pred = which(subj.pred==subj_unique.pred[i])
    pred.points <- newdata$argvals[sel.pred]
    mu.predi <- mu.pred[sel.pred]
    var.error.predi <- var.error.pred[sel.pred]
    
    y.predi = y.pred[sel.pred] - mu.predi
    sel.pred.obs = which(!is.na(y.predi))
    obs.points <- pred.points[sel.pred.obs]
     
    if(!is.null(obs.points)){
      var <- mean(var.error.predi[sel.pred.obs])
      if(var==0&length(sel.pred.obs) < npc)
        stop("Measurement error estimated to be zero and there are fewer observed points thans PCs; scores
               cannot be estimated.")
      B3i.pred = spline.des(knots=knots, x=pred.points, ord = p+1,outer.ok = TRUE,sparse=TRUE)$design
      B3i = spline.des(knots=knots, x=obs.points, ord = p+1,outer.ok = TRUE,sparse=TRUE)$design
      Chati = tcrossprod(B3i%*%Theta,B3i)
      
      if(length(sel.pred.obs)==1) Ri = var.error.predi[sel.pred.obs]
      if(length(sel.pred.obs)>1) Ri = diag(var.error.predi[sel.pred.obs])
      Vi.inv = as.matrix(solve(Chati + Ri))
      Vi.pred = as.matrix(tcrossprod(B3i.pred%*%Theta,B3i.pred))
      Hi = as.matrix(B3i.pred%*%tcrossprod(Theta,B3i)%*%Vi.inv)
      ui =tcrossprod(Theta,B3i)%*%Vi.inv %*%y.predi[sel.pred.obs]
      scores$u[i,] = as.vector(ui)
      y.pred[sel.pred] = as.numeric(Hi%*%y.predi[sel.pred.obs]) + mu.predi
      temp = as.matrix(B3i.pred%*%tcrossprod(Theta,B3i))
      if(length(sel.pred.obs) >1){
        cov.pred[sel.pred,sel.pred] = Vi.pred - temp%*%Vi.inv%*%t(temp)
      }
      if(length(sel.pred.obs) ==1){
        cov.pred[sel.pred,sel.pred] = Vi.pred - Vi.inv[1,1]*temp%*%t(temp)
      }
      ## predict scores
      if(object$calculate.scores==TRUE){ 
        temp = matrix(t(object$eigenfunctions),nrow=npc)%*%(as.matrix(Bnew)%*%ui)/sum(object$eigenfunctions[,1]^2)
        temp = as.matrix(temp)
        scores$scores[i,1:npc] = temp[,1]
      }
    }
  }
  
  B = spline.des(knots=knots, x=newdata$argvals, ord = p+1,outer.ok = TRUE,sparse=TRUE)$design
  Chat.pred = as.matrix(tcrossprod(B%*%Matrix(Theta),B))
  
  
  return(list(object=object,newdata=newdata,y.pred = y.pred,mu.pred=mu.pred,var.error.pred=var.error.pred,
              scores = scores, cov.pred = cov.pred, se.pred = sqrt(diag(cov.pred)),
              Chat.pred = Chat.pred,Chat.diag.pred = diag(Chat.pred)))  
}


