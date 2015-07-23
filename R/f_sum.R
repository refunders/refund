#' Sum computation 1
#' 
#' Internal function used compute a sum in FPCA-based covariance updates
#' 
#' @param mu.q.c
#' @param sig.q.c
#' @param theta
#' @param obspts.mat
#' 
#' @author Jeff Goldsmith \email{ajg2202@@cumc.columbia.edu}
#' 
f_sum = function(mu.q.c, sig.q.c, theta, obspts.mat){
  I = dim(mu.q.c)[1]
  kp = dim(mu.q.c)[2]
  kt = dim(theta)[1]
  ret.sum = matrix(0, kp*kt, kp*kt)
  
  for(i in 1:I){
    mu.mat = matrix(mu.q.c[i,], nrow = 1, ncol = kp)
    ret.sum = ret.sum + kronecker(t(mu.mat) %*% mu.mat + sig.q.c[[i]], (theta[,obspts.mat[i,]])%*%t(theta[,obspts.mat[i,]])) 
  }
  return(ret.sum)
}