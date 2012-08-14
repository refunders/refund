### Function to plot estimated regression function
plot.lpeer<- function(x, conf=0.95, ...){
  if(!exists("x")) return (cat("Error: The value specified in lfit argument is not an lpeer object.\n"))
  if(!is.list(x)) return (cat("Error: The value specified in lfit argument is not an lpeer object.\n"))
  if(is.na(match("lpeerobj", names(x)))) return (cat("Error: The value specified in lfit argument is not an lpeer object.\n"))
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



