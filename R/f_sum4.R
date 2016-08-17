#' Sum computation 2
#' 
#' Internal function used compute a sum in FPCA-based covariance updates
#' 
#' @param mu.q.c current value of mu.q.c
#' @param sig.q.c current value of sig.q.c
#' @param mu.q.bpsi current value of mu.q.bpsi
#' @param sig.q.bpsi current value of sig.q.bpsi
#' @param theta current value of theta
#' @param obspts.mat matrix indicating where curves are observed
#' 
#' @author Jeff Goldsmith \email{ajg2202@@cumc.columbia.edu}
#' 
f_sum4 = function(mu.q.c, sig.q.c, mu.q.bpsi, sig.q.bpsi, theta, obspts.mat){
  I = dim(mu.q.c)[1]
  kp = dim(mu.q.c)[2]
  kt = dim(theta)[2]
  ret.sum = matrix(0, 1, 1)
  
  for(i in 1:I){
    theta_i = t(theta)[,obspts.mat[i,]]
    temp = 
      f_trace(Theta_i = theta_i, Sig_q_Bpsi = sig.q.bpsi, Kp = kp, Kt = kt) %*% matrix(mu.q.c[i,], kp, 1) %*% matrix(mu.q.c[i,], 1, kp) +
      f_trace(Theta_i = theta_i, Sig_q_Bpsi = sig.q.bpsi, Kp = kp, Kt = kt) %*% sig.q.c[[i]] +
      t(mu.q.bpsi) %*% theta_i %*% t(theta_i) %*% mu.q.bpsi %*% sig.q.c[[i]]
    
    ret.sum = ret.sum + sum(diag(temp))   
  }
  return(ret.sum)
}


