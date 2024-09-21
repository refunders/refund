#' Function-on Scalar Regression with variable selection
#'
#' Implements an iterative algorithm for function-on-scalar regression with variable selection
#' by alternatively updating the coefficients and covariance structure.
#'
#' @param formula an object of class "\code{{formula}}": an expression of the model to be fitted.
#' @param data a data frame that contains the variables in the model.
#' @param nbasis number of B-spline basis functions used.
#' @param method group variable selection method to be used ("grLasso", "grMCP", "grSCAD" refer to group Lasso, group MCP and group SCAD, respectively) or "\code{ls}" for least squares estimation.
#' @param epsilon the convergence criterion.
#' @param max.iter_num maximum number of iterations.
#'
#' @return A fitted fosr.vs-object, which is a list with the following elements:
#' \item{formula}{an object of class "\code{{formula}}": an expression of the model to be fitted.}
#' \item{coefficients}{the estimated coefficient functions.}
#' \item{fitted.values}{the fitted curves.}
#' \item{residuals}{the residual curves.}
#' \item{vcov}{the estimated variance-covariance matrix when convergence is achieved.}
#' \item{method}{group variable selection method to be used or "\code{ls}" for least squares estimation.}
#' @references
#' Chen, Y., Goldsmith, J., and Ogden, T. (2016).
#' Variable selection in function-on-scalar regression. \emph{Stat} 5 88-101
#'
#' @author Yakuan Chen \email{yc2641@@cumc.columbia.edu}
#' @seealso \code{{grpreg}}
#' @importFrom grpreg grpreg cv.grpreg
#' @importFrom splines bs
#' @importFrom stats terms.formula
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(100)
#'
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
#' fit.fosr.vs = fosr.vs(Y~., data = data, method="grMCP")
#' plot(fit.fosr.vs)
#' }
#'
#'
fosr.vs <- function(formula, data, nbasis=10, method=c("ls", "grLasso", "grMCP", "grSCAD"), epsilon=1e-5, max.iter_num = 100){
	cl = match.call()
	mf = match.call(expand.dots = FALSE)
	m = match(c("formula", "data"), names(mf), 0)
	mf = mf[c(1,m)]
	mf$drop.unused.levels <- TRUE

	mf[[1]] = as.name("model.frame")
	mf = eval.parent(mf)

	# mf = model.frame(formula, data = data)
	Y = model.response(mf)
	W.des = model.matrix(formula, data = mf)
	N = dim(Y)[1]
	D = dim(Y)[2]
	K = dim(W.des)[2]
	Theta = bs(1:D, df=nbasis, intercept=TRUE, degree=3)
	Z = kronecker(W.des,Theta)

	Y.vec = as.vector(t(Y))
	ls <- lm(Y.vec~Z+0)
	a <- matrix(ls$residuals,nrow=N,ncol=D,byrow=T)
	if(D < 15) {
		b = cov(a)
	} else {
		w <- fpca.sc(a, var=T)
		b <- Reduce("+", lapply(1:length(w$evalues), function(x) {w$evalues[x] * w$efunctions[,x] %*% t(w$efunctions[,x])})) + w$sigma2 * diag(1,D,D)
	}
	f <- coef(ls)
	f <- replace(f, which(is.na(f)), 0)
	d <- rep(0,nbasis*K)
	num.iter = 0

  cat("Beginning iterative algorithm \n")

	while(sum((f-d)^2)>epsilon & num.iter <= max.iter_num){
	  num.iter = num.iter + 1
	  d <- f
	  c <- chol(solve(b))
	  Y_new = Y %*% t(c)
	  Theta_new = c %*% Theta
	  Z_new = kronecker(W.des, Theta_new)
	  Y.vecnew = as.vector(t(Y_new))
	  if (method == "ls"){
	    ls <- lm(Y.vecnew~Z_new)
	    f <- as.vector(coef(ls))[-1]
	  } else {
	    cvlam <- cv.grpreg(Z_new, Y.vecnew, group=switch(colnames(W.des)[1], "(Intercept)" = c(rep(1:K,each=nbasis))-1, c(rep(1:K,each=nbasis))), penalty = method, max.iter=5000)
	    f <- as.vector(coef(cvlam))[-1]
	  }
	  a <- matrix(Y.vec-Z%*%f,nrow=N,ncol=D,byrow=T)

	  if(D < 15) {
	  	b = cov(a)
	  } else {
	  	w <- fpca.sc(a, var=T)
	  	b <- Reduce("+", lapply(1:length(w$evalues), function(x) {w$evalues[x] * w$efunctions[,x] %*% t(w$efunctions[,x])})) + w$sigma2 * diag(1,D,D)
	  }

    if(num.iter %% 10 == 1) cat(".")

	}

	B <- matrix(f, nrow=K, byrow=T)
	est <- B %*% t(Theta)

  rownames(est) <- colnames(W.des)
  fitted <- W.des %*% est
  res <- Y - fitted
	ret <- list(call = cl, formula = formula, coefficients = est, fitted.values = fitted, residuals = res, vcov = b, method = method)
  class(ret) <- "fosr.vs"
	return(ret)
}
