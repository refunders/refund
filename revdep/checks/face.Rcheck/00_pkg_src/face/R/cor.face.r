###############################################
## Extraction of correlation and mean via FACE
## Author: Luo Xiao
## Email: lxiao5@ncsu.edu
## Date: Sep 1, 2015
#############################################

cor.face <- function(object,argvals.new,option="raw"){
  
  ## arguments
  #  "object": fitted object from function ''face.sparse''
  #  "argvals.new": at which time points to evaluate correlations,
  ##  a correlation matrix evaulated at these time points will be the output
  ##
  
  ## check inputs
  if(class(object)!="face.sparse") stop("'fit' has to be a face_sparse object")
  if(is.null(argvals.new)) stop("'argvals.new' needs to be specified")
  if(length(argvals.new)<2) stop("'argvals.new' needs to contain at least two time points")
  if(!is.numeric(argvals.new[1])) stop("'argvals.new' need to be numeric")
  
  ## create a hypothetical data at time points 'argvals.new'
  ## and then make prediction and estimation
  k <- length(argvals.new)
  
  newdata_pred <- data.frame(
    subj=rep(1,k),
    argvals = argvals.new,
    y = rep(1,k)
  )
  
  yhat <- predict.face.sparse(object,newdata_pred)
  
  ## extract correlation and mean
  Chat <- yhat$Chat.pred
  if(option=="raw") Chat <- Chat + diag(yhat$var.error.pred)
  Cov_diag <- diag(Chat)
  Cor <- diag(sqrt(1/Cov_diag))%*%Chat%*%diag(sqrt(1/Cov_diag))
  

  
  mu <- yhat$mu.pred
 ## ouptut
  
  res <- list("argvals.new" = argvals.new,
              "option" = option,
              "Cor" = Cor, #estimated correlation matrix at argvals.new
              "mu" = mu #estimated group/population mean at argvals.new
  )
  return(res)
}