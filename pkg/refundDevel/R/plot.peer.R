### Function to plot estimated regression function
plot.peer<- function(x, conf=0.95, ...){
	if(!exists("x")) return (cat("Error: The value specified in fit argument is not an peer object.\n"))
	if(!is.list(x)) return (cat("Error: The value specified in fit argument is not an lpeer object.\n"))
	if(is.na(match("peerobj", names(x)))) return (cat("Error: The value specified in fit argument is not an lpeer object.\n"))
	if(conf>0.99 | conf<0.70) return (cat("Error: Confidence level should be within 0.70 and 0.99\n"))
	status<- x$status
	par(mfrow=c(1,1))
	est<- x$GammaHat
	if(status==0) matplot(est, type='l', main='gamma', ...)
	if(status==1){
			ll<- x$GammaHat - qnorm(0.5+conf/2)*x$se.Gamma
			ul<- x$GammaHat + qnorm(0.5+conf/2)*x$se.Gamma
			matplot(est, type='l', ylim=range(est, ll, ul), 
                    	main='gamma', ...)
			matplot(ll, type='l', add=T, lty=2, col=2)
			matplot(ul, type='l', add=T, lty=2, col=2)
			}
	}
