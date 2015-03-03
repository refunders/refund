##' Plotting of estimated regression functions obtained through \code{peer()}
##'
##' Plots the estimate of components of estimated regression function obtained
##' from a \code{\link{peer}} object along with pointwise confidence bands.
##'
##' Pointwise confidence interval is displayed only if the user set \code{se=T}
##' in the call to \code{\link{peer}}, and does not reflect any multiplicity
##' correction.
##'
##' @param x object of class \code{"\link{peer}"}.
##' @param conf pointwise confidence level.
##' @param ylab y-axis label.
##' @param main title for the plot.
##' @param ... additional arguments passed to \code{\link{plot}}.
##' @author Madan Gopal Kundu \email{mgkundu@@iupui.edu}
##' @seealso \code{peer}, \code{lpeer}, \code{plot.lpeer}
##' @importFrom graphics matplot
##' @importFrom stats qnorm
##' @export
##' @references Kundu, M. G., Harezlak, J., and Randolph, T. W. (2012).
##' Longitudinal functional models with structured penalties. (Please contact
##' J. Harezlak at \email{harezlak@@iupui.edu}.)
##'
##' Randolph, T. W., Harezlak, J, and Feng, Z. (2012). Structured penalties for
##' functional linear models - partially empirical eigenvectors for regression.
##' \emph{Electronic Journal of Statistics}, 6, 323--353.
##' @examples
##'
##' data(DTI)
##' cca = DTI$cca[which(DTI$case == 1),]
##' DTI = DTI[which(DTI$case == 1),]
##' fit.cca.peer1 = peer(Y=DTI$pasat, funcs = cca)
##' plot(fit.cca.peer1)
### Function to plot estimated regression function
plot.peer<- function(x, conf=0.95, ylab='Estimated regression function', main=expression(gamma),...){
  if(!class(x)=='peer') return (cat("Error: The object is not an peer object.\n"))
  if(conf>0.99 | conf<0.70) return (cat("Error: Confidence level should be within 0.70 and 0.99\n"))
  status<- x$status
  est<- x$GammaHat
  if(status==0) matplot(est, type='l', ylab=ylab,
                        main=main, ...)
  if(status==1){
    ll<- est - qnorm(0.5+conf/2)*x$se.Gamma
    ul<- est + qnorm(0.5+conf/2)*x$se.Gamma
    matplot(est, type='l', ylim=range(est, ll, ul), ylab=ylab,
            main=main, ...)
    matplot(ll, type='l', add=T, lty=2, col=2)
    matplot(ul, type='l', add=T, lty=2, col=2)
  }
  abline(h=0)
}
