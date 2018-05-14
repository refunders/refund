covZ <- function(t_seq, Y, C, sigma2)
{
  #########################################################################
  #     covZ: calculate the covariance matrix of paired-y product
  #
  #     Input                          
  #     t_seq: defined grid point m-by-1
  #     Y: response vector m-by-1    
  #     C: covariance matrix of X m-by-m         
  #     sigma2: variance of random error
  #                                    
  #     Output                          
  #     CovZZ: covariance matrix of paired-y product m^2-by-m^2      
  #     Z.info: vectorized Z.tri and  the corresponding grid points    
  #     Z.tri: upper triangle matrix of paired-y product m-by-m
  #                                
  #########################################################################
  temp <- expand.grid(Y,Y)
  t.pair <- expand.grid(t_seq,t_seq)
  d <- length(Y)
  id <- expand.grid(1:d,1:d)
  m <- d

  Z <- apply(temp, 1, function(z) prod(z))
  Z.info <- as.matrix(cbind(Z,t.pair,id))
  Z.info <- Z.info[Z.info[,5] <= Z.info[,4],]
  
#   Z.mat <- matrix(Z,nrow=length(Y),ncol=length(Y),byrow=T)
#   Z.tri <- Z.mat
#   Z.tri[lower.tri(Z.tri)] <- NA
  
  Z.tri <- matrix(NA,nrow=length(Y),ncol=length(Y))
  Z.tri[lower.tri(Z.tri, diag=T)] <- Z.info[,1]
  Z.tri <- t(Z.tri)

  # Quick-n-dirty implementation (for testing only)
#   CovZZ <- matrix(NA,nrow=length(Z),ncol=length(Z))
#   time.start <- proc.time()
#   for (i in 1:m)
#   {
#     for (j in 1:m)
#     {
#       row.id <- d * (i - 1) + j
#       for (k in 1:m)
#       {
#         for (l in 1:m)
#         {
#           col.id <- d * (k - 1) + l
#           if (row.id <= col.id)
#           {
#             CovZZ[row.id,col.id] = C[i,k] * C[j,l] + C[i,l] * C[j,k] + 
#               (i==k) * (j==l) *sigma2^2 + (i==l) * (j==k) *sigma2^2 +
#               (C[i,k] * (j==l) + C[i,l] * (j==k) + C[j,k] * (i==l) 
#                + C[j,l] * (i==k)) * sigma2
#           }
# 
#         }
#       }
#     }
#   }
#   print((proc.time()-time.start)[3]) 
#   
#   # Avoid loops (even slower, for testing only)
#   time.start <- proc.time()
#   id1 <- as.matrix(expand.grid(1:d,1:d,1:d,1:d))
#   id2 <- as.matrix(expand.grid(1:d^2,1:d^2))
#   id.tab <- as.matrix(cbind(id1,id2))
#   id.tab <- id.tab[id.tab[,5] <= id.tab[,6],]
#   CovZZ2.vec = apply(id.tab,1,function(x)
#     C[x[2],x[4]] * C[x[1],x[3]] + C[x[2],x[3]] * C[x[1],x[4]] + 
#     (x[2]==x[4]) * (x[1]==x[3]) *sigma2^2 + 
#     (x[2]==x[3]) * (x[1]==x[4]) *sigma2^2 +
#     (C[x[2],x[4]] * (x[1]==x[3]) + C[x[2],x[3]] * (x[1]==x[4]) + 
#        C[x[1],x[4]] * (x[2]==x[3]) + C[x[1],x[3]] * (x[2]==x[4])) * sigma2)
#   
#   CovZZ2 <- matrix(NA,nrow=length(Z),ncol=length(Z))
#   CovZZ2[upper.tri(CovZZ2, diag=T)] <- CovZZ2.vec
#   print((proc.time()-time.start)[3])
  
  # Alternative vectorizing
  #time.start <- proc.time()
  id1 <- as.matrix(expand.grid(1:d,1:d,1:d,1:d))
  id2 <- as.matrix(expand.grid(1:d^2,1:d^2))
  id.tab <- as.matrix(cbind(id1,id2))
  id.tab <- id.tab[id.tab[,5] <= id.tab[,6],] # symmetric
  # kronecker delta terms
  # \delta_{jk},\delta_{j'k'},\delta_{jk'},\delta_{j'k}
  k.del <- cbind(id.tab[,1]==id.tab[,3],id.tab[,2]==id.tab[,4],
                 id.tab[,1]==id.tab[,4],id.tab[,2]==id.tab[,3])
  # Cov operator terms
  # C_{jk},C_{j'k'},C_{jk'},C_{j'k}
  c.term <- cbind(C[cbind(id.tab[,1],id.tab[,3])],
                  C[cbind(id.tab[,2],id.tab[,4])],
                  C[cbind(id.tab[,1],id.tab[,4])],
                  C[cbind(id.tab[,2],id.tab[,3])])
  # Form a matrix with pre-calculated terms
  id.tab <- cbind(id.tab,k.del,c.term)
  # CovZZ = C_{jk}C_{j'k'} + C_{jk'}C_{j'k} 
  # + (\delta_{jk}\delta_{j'k'} + \delta_{jk'}\delta_{j'k})\sigma_2^2 + 
  # (C_{jk}\delta_{j'k'} + C_{jk'}\delta_{j'k} + C_{j'k}\delta_{jk'} +
  # C_{j'k'}\delta_{jk})\sigma_2
  CovZZ3.vec = id.tab[,12] * id.tab[,11] + id.tab[,14] * id.tab[,13] + 
      (id.tab[,8] * id.tab[,7] + id.tab[,10] * id.tab[,9]) *sigma2^2 +
      (id.tab[,12] * id.tab[,7] + id.tab[,14] * id.tab[,9] + 
         id.tab[,13] * id.tab[,10] + id.tab[,11] * id.tab[,8]) * sigma2
  # CovZZ.vec = vech(CovZZ)
  CovZZ3 <- matrix(NA,nrow=length(Z),ncol=length(Z))
  CovZZ3[upper.tri(CovZZ3, diag=T)] <- CovZZ3.vec
  #print((proc.time()-time.start)[3])
  
  return (list(CovZZ = CovZZ3,CovZZ.vech = CovZZ3.vec,
               Z.info = Z.info,Z.tri = Z.tri))

}
