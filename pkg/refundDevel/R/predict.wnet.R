predict.wnet <- function(object, newx, covt = NULL){
	dim(newx) = c(nrow(newx), prod(dim(newx)[-1]))
	fhat = object$fhat
	dim(fhat) = c(prod(dim(fhat)), 1)
	const = object$const
	dim(const) = c(length(const), 1)
	yhat = cbind(rep(1, nrow(newx)), covt) %*% const + newx %*% fhat
	if (object$family == 'binomial'){
		yhat = 1 / (1 + exp(-yhat))
	}
	return(yhat)
}