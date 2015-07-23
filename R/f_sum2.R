#' Sum computation 2
#' 
#' Internal function used compute a sum in FPCA-based covariance updates
#' 
#' @param y
#' @param fixef
#' @param mu.q.c
#' @param kt
#' @param theta
#' 
#' @author Jeff Goldsmith \email{ajg2202@@cumc.columbia.edu}
#' 
f_sum2 = function(y, fixef, mu.q.c, kt, theta){
  I = dim(mu.q.c)[1]
  kp = dim(mu.q.c)[2]
  ret.sum = matrix(0, nrow = kp*kt, ncol = 1)
  
  for(i in 1:I){
    obs.pts = !is.na(y[i,])
    ret.sum = ret.sum + kronecker((matrix(mu.q.c[i,])), theta[,obs.pts]) %*% matrix(y[i, obs.pts] - fixef[i,obs.pts])
  }
  return(ret.sum)
}