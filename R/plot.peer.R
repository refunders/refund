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
