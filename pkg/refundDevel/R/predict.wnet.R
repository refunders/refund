predict.wnet <- function(object, newx, newcovt = NULL, ...){
	dim(newx) = c(nrow(newx), prod(dim(newx)[-1]))
	fhat = object$fhat
	dim(fhat) = c(prod(dim(fhat)), 1)
	coef.params = object$coef.params
	dim(coef.params) = c(length(coef.params), 1)
	yhat = cbind(rep(1, nrow(newx)), newcovt) %*% coef.params + newx %*% fhat
	if (object$family == 'binomial'){
		yhat = 1 / (1 + exp(-as.vector(yhat)))
	}
	return(yhat)
}