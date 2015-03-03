##' Plotting of estimated regression functions obtained through \code{lpeer()}
##'
##' Plots the estimate of components of estimated regression function obtained
##' from an \code{\link{lpeer}} object along with pointwise confidence bands.
##'
##' Pointwise confidence interval is displayed only if the user set \code{se=T}
##' in the call to \code{\link{lpeer}}, and does not reflect any multiplicity
##' correction.
##'
##' @param x object of class \code{"\link{lpeer}"}.
##' @param conf pointwise confidence level.
##' @param ... additional arguments passed to \code{\link{plot}}.
##' @author Madan Gopal Kundu \email{mgkundu@@iupui.edu}
##' @seealso \code{peer}, \code{lpeer}, \code{plot.peer}
##' @export
##' @importFrom graphics matplot
##' @importFrom stats qnorm
##' @references Kundu, M. G., Harezlak, J., and Randolph, T. W. (2012).
##' Longitudinal functional models with structured penalties. (Please contact
##' J. Harezlak at \email{harezlak@iupui.edu}.)
##'
##' Randolph, T. W., Harezlak, J, and Feng, Z. (2012). Structured penalties for
##' functional linear models - partially empirical eigenvectors for regression.
##' \emph{Electronic Journal of Statistics}, 6, 323--353.
##' @examples
##'
##' data(DTI)
##' cca = DTI$cca[which(DTI$case == 1),]
##' DTI = DTI[which(DTI$case == 1),]
##' fit.cca.lpeer1 = lpeer(Y=DTI$pasat, t=DTI$visit, subj=DTI$ID, funcs = cca)
##' plot(fit.cca.lpeer1)
### Function to plot estimated regression function
plot.lpeer<- function(x, conf=0.95, ...){
  if(!class(x)=='lpeer') return (cat("Error: The object is not an lpeer object.\n"))
  if(conf>0.99 | conf<0.70) return (cat("Error: Confidence level should be within 0.70 and 0.99\n"))
  d<- x$d
  status<- x$status
  if(d==0) par(mfrow=c(1,1))
  if(d==1) par(mfrow=c(1,2))
  if(d>1) par(mfrow=c(2,2))
  for(i in 0:d)
  {
    est<- x$GammaHat[,(i+1)]
    if(status==0) matplot(est, type='l', main=paste('gamma', i, sep=''), ...)
    if(status==1){
      ll<- x$GammaHat[,(i+1)] - qnorm(0.5+conf/2)*x$se.Gamma[,(i+1)]
      ul<- x$GammaHat[,(i+1)] + qnorm(0.5+conf/2)*x$se.Gamma[,(i+1)]
      matplot(est, type='l', ylim=range(est, ll, ul),
              main=paste('gamma', i, sep=''), ...)
      matplot(ll, type='l', add=T, lty=2, col=2)
      matplot(ul, type='l', add=T, lty=2, col=2)
    }
    abline(h=0)
  }
}


