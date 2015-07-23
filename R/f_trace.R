#' Trace computation
#' 
#' Internal function used compute a trace in FPCA-based covariance updates
#' 
#' @param Theta_i
#' @param Sig_q_Bpsi
#' @param Kp
#' @param Kt
#' 
#' @author Jeff Goldsmith \email{ajg2202@@cumc.columbia.edu}
f_trace = function(Theta_i, Sig_q_Bpsi, Kp, Kt){
  
  ret.mat = matrix(NA, nrow = Kp, ncol = Kp)
  A = Theta_i %*% t(Theta_i)
  
  for(i in 1:Kp){
    for(j in 1:Kp){
      ret.mat[i,j] = sum(diag(A %*% Sig_q_Bpsi[((-1 + i)*Kt + 1):(i*Kt), ((-1 + j)*Kt + 1):(j*Kt)]))
    }
  }
  
  return(ret.mat)
}
