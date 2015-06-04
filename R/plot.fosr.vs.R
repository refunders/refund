#' Plot for Function-on Scalar Regression with variable selection
#' 
#' Given a "\code{\link{fosr.vs}}" object, produces a figure of estimated coefficient functions.
#'
#' @param x an object of class "\code{\link{fosr.vs}}".
#' @param ... additional arguments.
#' 
#' @return a figure of estimated coefficient functions.
#' 
#' @author Yakuan Chen \email{yc2641@@cumc.columbia.edu}
#' @seealso \code{\link{fosr.vs}}
#' @import ggplot2
#' @export
#'
#' @examples
#' \dontrun{
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
#' plot(fit.mcp)
#' }
#' 
#'


plot.fosr.vs <- function(x, ...){
  p <- dim(coef(x))[1]
  D <- dim(coef(x))[2]
  df = as.data.frame(cbind(as.vector(sapply(1:p, function(x) seq(1,D)/D)), as.vector(t(coef(x))), as.vector(sapply(1:p, function(x) rep(x,D)))))
  names(df) = c("grid", "value", "curve")
  ggplot(df, aes(x=grid, y=value, group=curve, colour = factor(curve))) + geom_path() + xlab("") + ylab("") + theme(legend.title = element_blank()) + scale_color_manual(values=1:p,labels=rownames(coef(object)))
}