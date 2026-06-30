# Software Name : PenFFR
# SPDX-FileCopyrightText: Copyright (c) Orange SA
# SPDX-License-Identifier:  GPL-2.0-or-later

# This software is distributed under the GNU General Public License v2.0 or later
# see the "LICENSE.md" file for more details or https://spdx.org/licenses/GPL-2.0-or-later.html

# Authors: see CONTRIBUTORS.md
# Software description: Linear regression for fully functional data.

#' Categorial variable for taking the repeated observation of the response
#' @param Y.mat the matrix response with of n rows and m = max m_i columns with na values after the end of the observation for each sample
#' @param t.mat the matrix of time with the corresponded values of Y.mat
my_funmat1 <- function(t.mat, Y.mat){

  time <- data <- ind <- c()
  for (u in 1:nrow(Y.mat)) {
    tmp.time <- as.vector(na.omit(as.vector(unlist(t.mat[u,]))))
    tmp.data <- as.vector(na.omit(as.vector(unlist(Y.mat[u,]))))
    if (length(tmp.data) != length(tmp.time)){
      k <- min(length(tmp.data), length(tmp.time))
      tmp.data <- tmp.data[1:k]
      tmp.time <- tmp.time[1:k]
    }
    ind <- c(ind, rep(u, length(tmp.data)))
    time <- c(time, tmp.time)
    data <- c(data, tmp.data)
  }

  return(cbind.data.frame(output = data,
                          id = as.factor(ind),
                          time = time))
}

#' Functional representation as matrix : Functionnal covariates
#' @param Y.mat the matrix response with of n rows and m = max m_i columns with na values after the end of the observation for each sample
#' @param t.mat the matrix of time with the corresponded values of Y.mat
#' @param X_fd.list the list of functional predictors
#' @param basis.beta the (B-splines) for functional parameters
#' @import fda
#' @return Matrix that can be use to build the design matrix of the integral model
#' @export
my_funmat2 <- function(X_fd.list, basis.beta, t.mat, Y.mat){

  d <- length(X_fd.list)
  m <- ncol(t.mat)
  n <- ncol(X_fd.list[[1]][["coefs"]])

  # Build the response ----
  data <- my_funmat1(t.mat = t.mat, Y.mat = Y.mat)

  tmp.basis <- function(x){
    d <- length(X_fd.list)
    tmp <- lapply(1:d, function(l){
      tmp_phi <- eval.basis(evalarg = x, basisobj = X_fd.list[[l]][["basis"]])
      tmp_psi <- eval.basis(evalarg = x, basisobj = basis.beta)
      list(phi = tmp_phi, psi = tmp_psi)
    })
    return(tmp)
  }

  ## Basis functions of parameters
  basis_psi <- function(t){
    tmp_basis <- tmp.basis(t)

    b <- matrix(tmp_basis[[1]][["psi"]], nrow = 1)
    if (d>1){
      for (l in 2:d) {
        tmp <- tmp_basis[[l]][["psi"]]
        cpt <- length(tmp)
        b <- rbind(cbind(b, matrix(0, nrow = nrow(b), ncol = cpt)),
                   c(rep(0,ncol(b)), tmp))
      }
    }

    return(b)
  }

  ## Basis functions of covariates
  basis_phi <- function(t){
    tmp_basis <- tmp.basis(t)

    b <- matrix(tmp_basis[[1]][["phi"]], nrow = 1)
    if (d>1){
      for (l in 2:d) {
        tmp <- tmp_basis[[l]][["phi"]]
        cpt <- length(tmp)
        b <- rbind(cbind(b, matrix(0, nrow = nrow(b), ncol = cpt)),
                   c(rep(0,ncol(b)), tmp))
      }
    }
    return(b)
  }


  ## Coefficients of basis expansion of covariates
  coefs <- c()
  for (l in 1:d) {
    coefs <- rbind(coefs, X_fd.list[[l]][["coefs"]])
  }
  coefs <- rbind(1, coefs)

  ## Design matrix
  R.mat <- c()
  for (i in 1:n) {
    tmp.time <- data$time[data$id == i]
    tmp <- c()
    for (tj in tmp.time) {
      tmp.phi <- basis_phi(tj)
      tmp.psi <- basis_psi(tj)
      tmp.b <- eval.basis(evalarg = tj, basisobj = basis.beta)

      tmp.phi <- rbind(cbind(matrix(1), matrix(0, ncol = ncol(tmp.phi))),
                       cbind(matrix(0, ncol = 1, nrow = nrow(tmp.phi)), tmp.phi))
      tmp.psi <- rbind(cbind(tmp.b, matrix(0, ncol = ncol(tmp.psi))),
                       cbind(matrix(0, ncol = ncol(tmp.b), nrow = nrow(tmp.psi)),
                             tmp.psi))

      tmp <- rbind(tmp, t(t(tmp.psi) %*% tmp.phi %*% coefs[,i]))

    }
    R.mat <- rbind(R.mat, tmp)
  }

  data <- cbind.data.frame(data, data.frame(R.mat))
  colnames(data) <- c("output", "id", "time",
                      paste("X", 1:(ncol(data)-3), sep = "."))
  data$id <- as.factor(data$id)

  return(data)
}


#' Function that build the design matrix of the linear concurrent model
#' @param Y.mat the matrix response with of n rows and m = max m_i columns with na values after the end of the observation for each sample
#' @param t.mat the matrix of time with the corresponded values of Y.mat
#' @param X_fd.list the list of functional predictors
#' @param nbasis the number of basis for functional parameter (an integer)
#' @param n.order the order of the splines basis
#' @import fda
#' @return Dataframe that containt the design matrix of the linear concurrent model
#' @export
get.data1 <- function(X_fd.list, Y.mat, t.mat, nbasis, n.order){

  basis.beta <- create.bspline.basis(rangeval = range(t.mat, na.rm = T),
                                     nbasis = nbasis, norder = n.order)
  data <- my_funmat2(X_fd.list = X_fd.list,
                     basis.beta = basis.beta,
                     Y.mat = Y.mat,
                     t.mat = t.mat)
  return(data)
}


#' Penalty matrix
#' @param obs.grid the time observation grid
#' @param nbasis the number of basis for functional parameter (an integer)
#' @param n.order the order of the splines basis
#' @param d the number of functional predictors
#' @param deg the derivative order to consider (default 2)
#' @import fda
#' @return Matrix that can be use to penalize the model
#' @export
my_penmat1 <- function(obs.grid, nbasis, n.order, d, deg=2){

  basis.beta <- create.bspline.basis(rangeval = range(obs.grid, na.rm = T),
                                     nbasis = nbasis, norder = n.order)
  t0 <- range(obs.grid, na.rm = T)[1]
  tf <- range(obs.grid, na.rm = T)[2]
  L.beta <- nbasis

  vals <- lapply(1:1, function(l){

    tmp.b1 <- diag(rep(1,L.beta))

    tmp.fd <- fd(coef = tmp.b1, basisobj = basis.beta)
    tmp.fd2 <- deriv.fd(tmp.fd, deg)

    sapply(1:L.beta, function(j){
      t(sapply(1:L.beta, function(i){
        f <- function(x){
          tmp.eval <- as.vector(t(eval.fd(evalarg = x, fdobj = tmp.fd2)))
          return(abs(tmp.eval[i]*tmp.eval[j]))
        }
        integrate(Vectorize(f), lower = t0, upper = tf,
                  subdivisions=2000)[["value"]]
      }))
    })
  })
  for (l in 2:d) {
    vals[[l]] <- vals[[1]]
  }

  return(vals)
}


