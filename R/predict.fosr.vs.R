#' Prediction for Function-on Scalar Regression with variable selection
#' 
#' Given a "\code{\link{fosr.vs}}" object and new data, produces fitted values.
#'
#' @name predict.fosr.vs
#' @usage ## S3 method for class 'fosr.vs' 
#' predict.fosr.vs(object, newdata=NULL, ...)
#' @param object an object of class "\code{\link{fosr.vs}}".
#' @param newdata a data frame that contains the values of the model covariates at which predictors are required.
#' @param ... additional arguments.
#' 
#' @return fitted values.
#' 
#' @author Yakuan Chen \email{yc2641@@cumc.columbia.edu}
#' @seealso \code{\link{fosr.vs}}
#' @export
#'
#' @examples
#' I = 100
#' p = 20
#' D = 50
#' grid = seq(0, 1, length = D)
#' 
#' beta.true = matrix(0, p, D)
#' beta.true[1,] = sin(2*grid*pi)
#' beta.true[2,] = cos(2*grid*pi)
#' beta.true[3,] = 2
#' 
#' psi.true = matrix(NA, 2, D)
#' psi.true[1,] = sin(4*grid*pi)
#' psi.true[2,] = cos(4*grid*pi)
#' lambda = c(3,1)
#' 
#' set.seed(100)
#' 
#' X = matrix(rnorm(I*p), I, p)
#' C = cbind(rnorm(I, mean = 0, sd = lambda[1]), rnorm(I, mean = 0, sd = lambda[2]))
#' 
#' fixef = X%*%beta.true
#' pcaef = C %*% psi.true
#' error = matrix(rnorm(I*D), I, D)
#' 
#' Yi.true = fixef
#' Yi.pca = fixef + pcaef
#' Yi.obs = fixef + pcaef + error
#' 
#' data = as.data.frame(X)
#' data$Y = Yi.obs
#' fit.mcp = fosr.vs(Y~., data = data[1:80,], method="grMCP")
#' predicted.value = predict(fit.mcp, data[81:100,])
#'




predict.fosr.vs <- function(object, newdata=NULL, ...){
  x <- model.matrix(object$formula, newdata)
  y <- x %*% coef(object)
  return(y)
}