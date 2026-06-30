# Software Name : PenFFR
# SPDX-FileCopyrightText: Copyright (c) Orange SA
# SPDX-License-Identifier:  GPL-2.0-or-later

# This software is distributed under the GNU General Public License v2.0 or later
# see the "LICENSE.md" file for more details or https://spdx.org/licenses/GPL-2.0-or-later.html

# Authors: see CONTRIBUTORS.md
# Software description: Linear regression for fully functional data.

#' Calcul de l'intégrale des fonctions de base
#' @param t0 the initial time observation
#' @param tf the final time observation
#' @param basis.beta the order of the splines basis
#' @param basis.X the number of functional predictors
#' @import fda
#' @importFrom stats integrate
#' @return Matrix that compute the integral of basis functions
#' @export
basis.integral <- function(t0, tf, basis.beta, basis.X){
  "
  tf           :  Borne supérieure de l'intégrale
  t0           :  Borne inférieure de l'intégrale
  basis.beta   :  Vecteur de base de la décomposition du paramètre fonctionnel bivariée
  basis.X      :  Vecteur de base de la décomposition du prédicteur fonctionnel
  "

  L.X <- length(basis.X[["names"]])
  L.beta <- length(basis.beta[["names"]])

  tmp.mat <- sapply(1:L.X, function(i){
    sapply(1:L.beta, function(j){
      f <- function(x){
        b1 <- eval.basis(evalarg = x, basisobj = basis.beta)[j]
        b2 <- eval.basis(evalarg = x, basisobj = basis.X)[i]
        return(b1*b2)
      }
      stats::integrate(Vectorize(f), lower = t0, upper = tf)[["value"]]
    })
  })

  vals <- matrix(0,ncol = L.X, nrow = L.beta**2)
  for (j in 1:L.beta) {
    vals[(1:L.beta)+L.beta*(j-1),] <- matrix(rep(tmp.mat[j,], L.beta), nrow = L.beta, byrow = T)
  }

  return(vals)
}

#' Penalty matrix for the integral model
#' @param obs.grid the time observation grid
#' @param L.beta the number of basis functions of parameters (an integer)
#' @param n.order the order of the splines basis
#' @param d the number of functional predictors
#' @import fda
#' @importFrom stats integrate
#' @return Matrix that can be use to penalize the integral model
#' @export
my_penmat2 <- function(obs.grid, L.beta, n.order, d){

  basis.beta <- create.bspline.basis(rangeval = range(obs.grid, na.rm = T),
                                     nbasis = L.beta, norder = n.order)
  t0 <- min(obs.grid, na.rm = T)
  tf <- max(obs.grid, na.rm = T)

  vals <- list()

  source("R/func-to-mat1.R")
  vals[[1]] <- my_penmat1(obs.grid = obs.grid, nbasis = L.beta,
                          n.order = n.order, d = 1)[[1]]

  vals[[2]] <- lapply(1:1, function(l){

    tmp.fd <- fd(coef = diag(L.beta), basisobj = basis.beta)

    vect <- as.vector(t(sapply(1:L.beta, function(j){
      sapply(1:L.beta, function(i){
        f <- function(x){
          tmp.eval <- as.vector(t(eval.fd(evalarg = x, fdobj = tmp.fd)))
          tmp.vect <- abs(as.vector(t(tmp.eval %*% t(tmp.eval))))
          return(tmp.vect[i]*tmp.vect[j])
        }
        stats::integrate(Vectorize(f), lower = t0, upper = tf,
                  subdivisions=2000)[["value"]]
      })
    })))

    vect.1 <- as.vector(t(sapply(1:L.beta, function(j){
      sapply(1:L.beta, function(i){
        f <- function(x){
          tmp.eval <- as.vector(t(eval.fd(evalarg = x, fdobj = deriv.fd(tmp.fd, 1))))
          tmp.vect <- abs(as.vector(t(tmp.eval %*% t(tmp.eval))))
          return(tmp.vect[i]*tmp.vect[j])
        }
        stats::integrate(Vectorize(f), lower = t0, upper = tf,
                  subdivisions=2000)[["value"]]
      })
    })))

    vect.2 <- as.vector(t(sapply(1:L.beta, function(j){
      sapply(1:L.beta, function(i){
        f <- function(x){
          tmp.eval <- as.vector(t(eval.fd(evalarg = x, fdobj = deriv.fd(tmp.fd, 2))))
          tmp.vect <- abs(as.vector(t(tmp.eval %*% t(tmp.eval))))
          return(tmp.vect[i]*tmp.vect[j])
        }
        stats::integrate(Vectorize(f), lower = t0, upper = tf,
                  subdivisions=2000)[["value"]]
      })
    })))

    ## Now we build the three term of the penalization
    term.1 <- matrix(vect.2, ncol = 1) %*% matrix(vect, nrow = 1)
    term.2 <- 2 * (matrix(vect.1, ncol = 1) %*% matrix(vect.1, nrow = 1))
    term.3 <- matrix(vect, ncol = 1) %*% matrix(vect.2, nrow = 1)

    term.1 + term.2 + term.3
  })[[1]]

  if (d > 2) {
    for (l in 3:d) {
      vals[[l]] <- vals[[2]]
    }
  }

  return(vals)
}

#' Functional representation as matrix : Functionnal covariates
#' @param X_fd.list the list of functional predictors
#' @param basis.beta the (B-splines) for functional parameters
#' @param Y.mat the matrix response with of n rows and m = max m_i columns with na values after the end of the observation for each sample
#' @param t.mat the matrix of time with the corresponded values of Y.mat
#' @import fda
#' @return Matrix that can be use to build the design matrix of the integral model
#' @export
my_funmat3 <- function(X_fd.list, basis.beta, t.mat, Y.mat){

  source("R/func-to-mat1.R")

  d <- length(X_fd.list)
  n <- ncol(X_fd.list[[1]][["coefs"]])
  L.beta <- length(basis.beta[["names"]])

  ## Build the response ----
  data <- my_funmat1(t.mat = t.mat, Y.mat = Y.mat)

  if (!(file.exists(paste(paste("BasisIntegral", L.beta, sep = ""))))){

    dir.create(paste(paste("BasisIntegral", L.beta, sep = "")))
    ## Compute the basis function integral
    uniq <- unique(as.vector(t.mat))
    tmp.b <- lapply(1:d, function(l){

      L.X <- nrow(X_fd.list[[l]][["coefs"]])
      tmp <- matrix(0, ncol = L.X, nrow = L.beta**2)
      saveRDS(tmp, file = paste(paste("BasisIntegral", L.beta, sep = ""),
                                "/mat",l, 1,".rds", sep = ""))
      rm(tmp)

      for (j in 2:length(uniq)) {
        tmp <- basis.integral(t0 = uniq[j-1],
                              tf = uniq[j],
                              basis.beta = basis.beta,
                              basis.X = X_fd.list[[l]][["basis"]])
        saveRDS(tmp, file = paste(paste("BasisIntegral", L.beta, sep = ""),
                                  "/mat",l, j,".rds", sep = ""))
        rm(tmp)
      }
      l
    })
  }

  ## R.mat at a point t for a covariate l and subject i ----
  R.mat <- function(t, l, i){

    coefs.i <- as.vector(X_fd.list[[l]][["coefs"]][,i])

    b <- eval.basis(evalarg = t, basisobj = basis.beta)
    B1 <- t(diag(rep(b, length(b))))

    B <- eval.basis(evalarg = t, basisobj = X_fd.list[[l]][["basis"]])

    obs.grid <- as.vector(na.omit(as.vector(unlist(t.mat[i,]))))

    j <- which(obs.grid == t)
    B2 <- matrix(0, nrow = L.beta**2, ncol = length(coefs.i))
    for (k in 1:j) {
      tmp.B2 <- readRDS(file = paste(paste("BasisIntegral", L.beta, sep = ""),
                                     "/mat",l, k,".rds", sep = ""))
      B2 <- B2 + as.matrix(tmp.B2)
    }

    ### Vecteur R ----
    return(coefs.i %*% t(B2) %*% B1)
  }

  ## Build the design matrix ----
  R <- lapply(1:n, function(i){
    tmp.time <- data$time[data$id == i]

    t(sapply(tmp.time, function(tj) {
      B0 <- eval.basis(evalarg = tj, basisobj = basis.beta)
      c(B0, as.vector(sapply(1:d, function(l){
        R.mat(tj, l, i)
      })))
    }))

  })

  new.R <- c()
  for (i in 1:n) {
    new.R <- rbind(new.R, R[[i]])
  }
  rm(R)

  data <- cbind.data.frame(data, data.frame(new.R))
  colnames(data) <- c("output","id", "time",
                      paste("X", 1:(ncol(data)-3), sep = "."))
  data$id <- as.factor(data$id)

  unlink(paste(paste("BasisIntegral", L.beta, sep = "")),
         recursive = TRUE)

  return(data)
}

#' Function that build the design matrix of the linear integral model
#' @param X_fd.list the list of functional predictors
#' @param nbasis the number of basis for functional parameters
#' @param Y.mat the matrix response with of n rows and m = max m_i columns with na values after the end of the observation for each sample
#' @param t.mat the matrix of time with the corresponded values of Y.mat
#' @param n.order the order of the splines basis
#' @import fda
#' @return Matrix that can be use to build the design matrix of the integral model
#' @export
get.data2 <- function(X_fd.list, Y.mat, t.mat, nbasis, n.order = 4){
  source("R/func-to-mat1.R")
  basis.beta <- create.bspline.basis(rangeval = range(t.mat, na.rm = T),
                                     nbasis = nbasis, norder = n.order)
  data <- my_funmat3(X_fd.list = X_fd.list,
                     basis.beta = basis.beta,
                     Y.mat = Y.mat,
                     t.mat = t.mat)
  return(data)
}
