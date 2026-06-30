# Software Name : PenFFR
# SPDX-FileCopyrightText: Copyright (c) Orange SA
# SPDX-License-Identifier:  GPL-2.0-or-later

# This software is distributed under the GNU General Public License v2.0 or later
# see the "LICENSE.md" file for more details or https://spdx.org/licenses/GPL-2.0-or-later.html

# Authors: see CONTRIBUTORS.md
# Software description: Linear regression for fully functional data.

#' Mean Square Error (MSE)
my.MSE <- function(betaact.vect, obs.grid, betahat.fd){
  b <- eval.fd(evalarg = obs.grid, fdobj = betahat.fd)
  tmp <- mean((betaact.vect - b)**2)
  return(tmp**(0.5))
}

#' Functionnal R2
fR2 <- function(Y, Yhat){
  Ybar <- matrix(rep(colMeans(Y),nrow(Y)),
                 nrow = nrow(Y), byrow = T)
  tmp1 <- colSums((Y - Yhat)**2)
  tmp2 <- colSums((Y - Ybar)**2)

  return(1 - mean(tmp1/tmp2))
}

#' Mean Square Prediction Error (MSPE)
my.MSPE <- function(Y.act, Y.hat, t.mat){
  Y.act <- as.vector(t(Y.act))
  Y.hat <- as.vector(t(Y.hat))
  m <- sapply(1:nrow(t.mat), function(u){
    length(as.vector(na.omit(as.vector(unlist(t.mat[u,])))))
  })
  m2 <- cumsum(m)
  tmp.1 <- sum((Y.act[1:m2[1]] - Y.hat[1:m2[1]])**2)/sum(Y.act[1:m2[1]]**2)
  if (nrow(t.mat)>1) {
    tmp.all <- sapply(2:nrow(t.mat), function(u){
      tmp1 <- Y.act[(m2[u-1]+1):(m2[u])]
      tmp2 <- Y.hat[(m2[u-1]+1):(m2[u])]

      sum((tmp1 - tmp2)**2)/sum(tmp1**2)
    })

    return(mean(c(tmp.1, tmp.all)))
  } else {
    return(tmp.1)
  }
}

#' Integrated Squared Error
my_ISE <- function(Y.act, Y.hat){

  # len <- sapply(1:nrow(t.mat), function(u){
  #   length(as.vector(na.omit(as.vector(unlist(t.mat[u,])))))
  # })
  #
  # if (nrow(t.mat)>1) {
  #   tmp.1 <- sum((Y.act[1:len[1]] - Y.hat[1:len[1]])**2)
  #   tmp.all <- sapply(1:(nrow(t.mat)-1), function(u){
  #     tmp1 <- Y.act[(sum(len[1:u])+1):(sum(len[1:(u+1)]))]
  #     tmp2 <- Y.hat[(sum(len[1:u])+1):(sum(len[1:(u+1)]))]
  #
  #     sum((tmp1 - tmp2)**2)
  #   })
  #
  #   return(c(tmp.1, tmp.all))
  #
  # } else {
  #   tmp.1 <- sum((Y.act[1:len] - Y.hat[1:len])**2)
  #   return(tmp.1)
  # }
  return(rowSums((Y.act - Y.hat)**2))
}

#' Predict in mixture of experts
predict.moe <- function(params, params.gated, data, data.gated, m){

  if (ncol(params.gated) != 1) {
    K <- ncol(params)

    ## Prediction in all classes
    pred <- sapply(1:K, function(k){
      as.matrix(data[,-c(1:3)]) %*% params[,k]
    })

    ## multinomial logit model
    tmp <- sapply(1:K, function(k){
      exp(cbind(1, as.matrix(data.gated)) %*% as.vector(params.gated[,k]))
    })
    pik <- tmp/rowSums(tmp)

    # ## posterior probabilities
    # pik <- sapply(1:K, function(k){
    #   sapply(1:nrow(pred), function(i){
    #     pi_ik[i,k] * dnorm(x = pred[i,k], mean = pred[i,k], sd =  sqrt(sig[k]))
    #   })
    # })
    # pik <- pik/rowSums(pik)

    ## Deduce the belonging class
    z <- apply(pik, 1, function(x) which.max(x))

    clust <- data.frame(x = z)
    tab <- table(clust$x)
    if (length(tab) == 1) {
      #print("one class")
      Y_hat <- matrix(pred[,clust$x[1]], ncol = m, byrow = T)
    } else {
      clust$x <- as.factor(clust$x)
      pik <- data.frame(model.matrix(~x - 1, data = clust))
      #Y_hat <- matrix(rowSums(pik * pred), ncol = m, byrow = T)
      Y_hat <- rowSums(pik * pred)
    }

    return(list(pred = Y_hat, class = z))
  } else {
    ## Prediction in all classes
    Y_hat <- as.matrix(data[,-c(1:3)]) %*% c(params)
    ## Prediction
    #Y_hat <- matrix(pred, ncol = m, byrow = T)

    return(list(pred = Y_hat, class = rep(1, length(pred))))
  }
}

#' Random effect matrix
Rdeff <- function(n,m){
  v <- matrix(1, nrow = m)
  tmp.mat <- v
  if (n>1) {
    for (i in 1:(n-1)) {
      tmp.mat <- rbind(cbind(tmp.mat, matrix(0, nrow = nrow(tmp.mat))),
                       cbind(matrix(0, nrow = m, ncol = ncol(tmp.mat)),v))
    }
  }
  return(tmp.mat)
}


#' Grid on multivariate unit ball
#' @import mvtnorm
generate_Grid<- function(nS,nR,n0,d){

  # inverse of the density of radius
  f <- function(x){
    exp(log(x)/d)
  }

  # unit d-dimensional vectors
  tmp_U <- rmvnorm(n=nS*nR, mean = rep(0,d), sigma = diag(d))
  rad <- f(runif(n=nR*nS))
  U_nS <- (tmp_U/rowNorms(tmp_U))*matrix(rep(rad, each = d),
                                         ncol=d,byrow=T)
  Y0=matrix(0,d,n0)
  return(cbind(Y0,t(U_nS)))
}


#' Multivariate quantile
#' @import transport
mvt.quant <- function(Y, Grid, tau = 0.05){

  opt_assign = function(Y, Grid){
    sizeY = length(Y[1, ])
    sizeG = length(Grid[1, ])
    distMat = as.matrix(dist(rbind(t(Y), t(Grid))))
    distMatsub = distMat[1:sizeY, (sizeY+1):(sizeY+sizeG)]^2
    return(transport(rep(1, sizeY), rep(1, sizeG), costm = distMatsub,
                     method = "networkflow"))
  }

  # set the level of quantile
  tau <- exp(log(1-tau)/nrow(Grid))

  # Optimal transport of data distribution to the closeball
  f <- opt_assign(Y = Y, Grid = Grid)
  f_sort <- f[order(f$from, f$mass),]
  f2 <- subset(f_sort, !duplicated(from, fromLast = T))
  ind <- which(colNorms(Grid) > tau)
  f3 <- unique(f2[(f2$to %in% ind),]$from)
  G2 <- Grid[,f3]
  Y2 <- Y[,f3]

  # compute the center-outward rank
  Rk <- 151*(colNorms(as.matrix(G2)))
  w <- 1/(Rk/sum(Rk))

  return(apply(Y2, 1, max))
}

##'Compute the outlier ratio (OR)
OR <- function(lwr, target, upr){
  s <- sum(sapply(1:length(lwr), function(j){
    1-as.numeric(between(x = target[j], left = lwr[j],#-abs(rnorm(1)),
                         right = upr[j]))# + abs(rnorm(1))))
  }))/length(lwr)
  return(s)
}
