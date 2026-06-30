# Software Name : PenFFR
# SPDX-FileCopyrightText: Copyright (c) Orange SA
# SPDX-License-Identifier:  GPL-2.0-or-later

# This software is distributed under the GNU General Public License v2.0 or later
# see the "LICENSE.md" file for more details or https://spdx.org/licenses/GPL-2.0-or-later.html

# Authors: see CONTRIBUTORS.md
# Software description: Linear regression for fully functional data.

#' Functional integral model
#' @param Y.mat the matrix response with of n rows and m = max m_i columns with na values after the end of the observation for each sample
#' @param t.mat the matrix of time with the corresponded values of Y.mat
#' @param X_fd.list the list of functional predictors
#' @param X.scal the data frame of the non functional predictors
#' @param nbasis the number of basis for functional parameter (an integer)
#' @param n.order the order of the splines basis
#' @param pen the boolean value that indicates if we want to use the penalization
#' @import fda
#' @import latex2exp
#' @return a list that containt
#' \itemize{
#'   \item \code{model}: the fitting model object
#'   \item \code{beta.fd}: the estimated functional parameters.
#'   \item \code{beta.scal}: the estimated scalar parameters
#'   \item \code{lambda}: the optimal penalized parameter
#'   \item \code{n.grid}: the length of time observation (useful for the predict function)
#' }
#' @export
penffr2 <- function(Y.mat, t.mat = "all.ok", X_fd.list, X.scal = data.frame(), nbasis  = "all.ok",
                    n.order = 4, pen = T, weights = NULL){

  '
  Y.mat is matrix response with of n rows and m = max m_i columns with na values after the end of the observation for each sample
  t.mat is a matrix of time with the corresponded values of Y.mat
  X_fd.list is the list of functional predictors
  X.scal is the data frame of the non functional predictors
  nbasis is an integer for the number of basis for functional parameter
  Pen is the boolean value for the penalization
  '

  # Check if the t.mat is not provided ----
  n <- nrow(Y.mat)
  m <- ncol(Y.mat)
  d <- length(X_fd.list)
  if (all(t.mat == "all.ok")) {
    obs.grid <- ((1:m)-1)/(m-1)
    t.mat <- matrix(rep(obs.grid, n), nrow = n, byrow = T)
  }

  # Check if nbasis is provided ----
  if (all(nbasis == "all.ok")){
    nbasis = min(20, floor(((n*m)/(d+1))**0.5))
  }

  # Build the dataframe and penalty matrix for the model ----
  source("R/func-to-mat2.R")
  data <- get.data2(X_fd.list = X_fd.list,
                    Y.mat = Y.mat, nbasis = nbasis,
                    t.mat = t.mat, n.order = n.order)

  if (ncol(X.scal) > 0){
    # Formatting the non-functional data ----
    X.scal.new <- data.frame()
    for (i in 1:n) {
      ni <- length(data$time[data$id == i])
      tmp <- matrix(rep(as.vector(unlist(X.scal[i,])), each = ni),
                    nrow = ni, byrow = F)
      X.scal.new <- rbind(X.scal.new, tmp)
    }
    X.scal.new <- data.frame(X.scal.new)
    colnames(X.scal.new) <- paste("Z", 1:ncol(X.scal.new), sep = ".")
    data <- cbind.data.frame(data, X.scal.new)
  }

  if (pen == F){

    ## Run the model ----
    my_formula <- as.formula(paste("output ~ 0", paste(colnames(data)[-(1:3)],
                                                       collapse = " + "), sep = " + "))
    model <- lm(my_formula, data = cbind.data.frame(data))

    ## get the functional parameters ----
    params <- as.vector(model$coefficients)
    params[is.na(params)] <- 0
    beta.fd <- lapply(1:(d+1), function(l){
      if (l == 1){
        fd(coef = params[1:nbasis],
           basisobj = create.bspline.basis(rangeval = range(obs.grid),
                                           nbasis = nbasis, norder = n.order),
           fdnames = list(main = TeX("$\\beta_0(t)"),
                          xlab = "Time",
                          ylab = TeX("$\\beta_0(t)")))
      } else{
        bifd(coef = matrix(params[nbasis+(((nbasis**2)*(l-2)+1):((nbasis**2)*(l-1)))],
                           nbasis, nbasis, byrow = F),
             sbasisobj = create.bspline.basis(rangeval = range(obs.grid),
                                              nbasis = nbasis, norder = n.order),
             tbasisobj = create.bspline.basis(rangeval = range(obs.grid),
                                              nbasis = nbasis, norder = n.order))
      }
    })

    if (ncol(X.scal) > 0) {
      ## get the non functional parameters ----
      beta.scal <- params[-c(1:(nbasis + (nbasis**2)*d))]
      ## results ----
      return(list(model = model, beta.fd = beta.fd, beta.scal = beta.scal, n.grid = m))
    } else {
      ## results ----
      return(list(model = model, beta.fd = beta.fd, n.grid = m))
    }

  } else {

    source("R/func-to-mat2.R")
    ## Build the penalty matrix ----
    vals <- my_penmat2(obs.grid = t.mat, d = length(X_fd.list)+1,
                       L.beta = nbasis, n.order = n.order)

    # Run the functional model ----
    source("R/Pensim2.R")
    d <- length(X_fd.list)
    n.lam <- min(10, floor(exp(log(100)/(d+1))))
    fofreg <- Pensim2(data = data, lams = seq(0, 5, length.out = n.lam),
                      vals = vals, d = d, d.z = ncol(X.scal),
                      t.mat = t.mat, weights = weights)

    model <- fofreg$model

    ## get the functional parameters ----
    params <- model$coefficients
    params[is.na(params)] <- 0
    beta.fd <- lapply(1:(d+1), function(l){
      if (l == 1){
        fd(coef = params[1:nbasis],
           basisobj = create.bspline.basis(rangeval = range(obs.grid),
                                           nbasis = nbasis, norder = 4),
           fdnames = list(main = TeX("$\\beta_0(t)"),
                          xlab = "Time",
                          ylab = TeX("$\\beta_0(t)")))
      } else{
        bifd(coef = matrix(params[nbasis+(((nbasis**2)*(l-2)+1):((nbasis**2)*(l-1)))],
                           nbasis, nbasis, byrow = F),
             sbasisobj = create.bspline.basis(rangeval = range(obs.grid),
                                              nbasis = nbasis, norder = 4),
             tbasisobj = create.bspline.basis(rangeval = range(obs.grid),
                                              nbasis = nbasis, norder = 4))
      }
    })

    if (ncol(X.scal) > 0) {
      ## get the non functional parameters ----
      beta.scal <- params[-c(1:(nbasis + (nbasis**2)*d))]
      ## results ----
      return(list(model = model, beta.fd = beta.fd, n.grid = m,
                  beta.scal = beta.scal, lambda = fofreg$lambda))
    } else {
      ## results ----
      return(list(model = model, beta.fd = beta.fd,
                  lambda = fofreg$lambda, n.grid = m))
    }
  }
}

#' Prediction of Functional integral model
#' @param model the functional model
#' @param t.mat the matrix of time on where we want the prediction
#' @param newX_fd.list the list of functional predictors
#' @param newX.scal the data frame of the non functional predictors
#' @return the prediction of functional integral model
#' @export
pred.penffr2 <- function(model, newX_fd.list = list(), newX.scal = data.frame(), t.mat = "all.ok"){

  '
  t.mat must be the time matrix of predictions with of n rows and m columns
  newX_fd.list is the list of functional predictors
  newX.scal is the data frame of the non functional predictors
  model is the penffr model
  '

  # Check the providing of new functional predictors ----
  if (length(newX_fd.list) == 0) {
    if (nrow(newX.scal) == 0) {
      pred <- predict(model$model)
    } else {
      return(NA)
    }

  } else {

    # observation grid of predictions ----
    n <- ncol(newX_fd.list[[1]][["coefs"]])
    obs.grid <- seq(newX_fd.list[[1]][["basis"]][["rangeval"]][1],
                    newX_fd.list[[1]][["basis"]][["rangeval"]][2],
                    length.out = model$n.grid)
    t.mat <- matrix(rep(obs.grid, n), nrow = n, byrow = T)

    # Build the dataframe and penalty matrix for the model ----
    source("R/func-to-mat1.R")
    nbasis <- length(as.vector(model$beta.fd[[1]][["coefs"]]))
    data <- get.data2(X_fd.list = newX_fd.list,
                      Y.mat = t.mat, nbasis = nbasis,
                      t.mat = t.mat)

    if (ncol(newX.scal) > 0){
      # Formatting the non-functional data ----
      X.scal.new <- data.frame()
      for (i in 1:n) {
        ni <- length(data$time[data$id == i])
        tmp <- matrix(rep(as.vector(unlist(newX.scal[i,])), each = ni),
                      nrow = ni, byrow = F)
        X.scal.new <- rbind(X.scal.new, tmp)
      }
      X.scal.new <- data.frame(X.scal.new)
      colnames(X.scal.new) <- paste("Z", 1:ncol(X.scal.new), sep = ".")
      data <- cbind.data.frame(data, X.scal.new)
    }

    # prediction ----
    pred <- predict(model$model, newdata = data, type = "response")
  }

  # result ----
  return(pred)

}
