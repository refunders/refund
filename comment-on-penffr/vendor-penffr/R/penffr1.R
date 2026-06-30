# Software Name : PenFFR
# SPDX-FileCopyrightText: Copyright (c) Orange SA
# SPDX-License-Identifier:  GPL-2.0-or-later

# This software is distributed under the GNU General Public License v2.0 or later
# see the "LICENSE.md" file for more details or https://spdx.org/licenses/GPL-2.0-or-later.html

# Authors: see CONTRIBUTORS.md
# Software description: Linear regression for fully functional data.

#' Functional concurrent model
#' @param Y.mat the matrix response with of n rows and m = max m_i columns with na values after the end of the observation for each sample
#' @param t.mat the matrix of time with the corresponded values of Y.mat
#' @param X_fd.list the list of functional predictors
#' @param X.scal the data frame of the non functional predictors
#' @param nbasis the number of basis for functional parameter (an integer)
#' @param n.order the order of the splines basis
#' @param pen the boolean value that indicates if we want to use the penalization
#' @import fda
#' @return a list that containt
#' \itemize{
#'   \item \code{model}: the fitting model object
#'   \item \code{beta.fd}: the estimated functional parameters.
#'   \item \code{beta.scal}: the estimated scalar parameters
#'   \item \code{lambda}: the optimal penalized parameter
#'   \item \code{n.grid}: the length of time observation (useful for the predict function)
#' }
#' @export
penffr1 <- function(Y.mat, t.mat = "all.ok", X_fd.list, X.scal = data.frame(), nbasis = "all.ok",
                    n.order = 4, pen = T, weights = NULL){

  '
  Y.mat       : matrix response with of n rows and m = max m_i columns with na values after the end of the observation for each sample
  t.mat       : a matrix of time with the corresponded values of Y.mat
  X_fd.list   : list of functional predictors
  X.scal      : data frame of the non functional predictors
  nbasis      : integer for the number of basis for functional parameter
  Pen         : boolean value for the penalization
  '

  # Check if the t.mat is not provided ----
  n <- nrow(Y.mat)
  m <- ncol(Y.mat)
  d <- length(X_fd.list)
  if (all(t.mat == "all.ok")) {
    obs.grid <- seq(0,1,length.out=m)
    t.mat <- matrix(rep(obs.grid, n), nrow = n, byrow = T)
  }

  # Check if nbasis is provided ----
  if (all(nbasis == "all.ok")){
    nbasis = min(20, floor((n*m)/(d+1)))
  }

  # Build the dataframe and penalty matrix for the model ----
  source("R/func-to-mat1.R")
  data <- get.data1(X_fd.list = X_fd.list,
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
  #saveRDS(data, file = "data.rds")

  if (pen == F){

    ## Run the model ----
    my_formula <- as.formula(paste("output ~ 0", paste(colnames(data)[-(1:3)],
                                                       collapse = " + "), sep = " + "))
    model <- lm(my_formula, data = cbind.data.frame(data))

    ## get the functional parameters ----
    params <- model$coefficients
    params[is.na(params)] <- 0
    beta.fd <- lapply(1:(d+1), function(l){
      if (l != 1){
        fd(coef = params[(nbasis*(l-1)+1):(nbasis*l)],
           basisobj = create.bspline.basis(rangeval = range(t.mat, na.rm = T),
                                           nbasis = nbasis, norder = 4),
           fdnames = list(main = paste("beta_",(l-1), sep = ""),
                          xlab = "Time (millisec.)",
                          ylab = paste("beta_",(l-1),sep = "")))
      } else{
        fd(coef = params[(nbasis*(l-1)+1):(nbasis*l)],
           basisobj = create.bspline.basis(rangeval = range(t.mat, na.rm = T),
                                           nbasis = nbasis, norder = 4),
           fdnames = list(main = paste("beta_",(l-1),sep = ""),
                          xlab = "Time (millisec.)",
                          ylab = paste("beta_",(l-1),sep = "")))
      }
    })

    if (ncol(X.scal) > 0) {
      ## get the non functional parameters ----
      beta.scal <- model$coefficients[-c(1:(nbasis*(d+1)))]
      ## results ----
      return(list(model = model, beta.fd = beta.fd, beta.scal = beta.scal,
                  n.grid = m))
    } else {
      ## results ----
      return(list(model = model, beta.fd = beta.fd, n.grid = m))
    }

  } else {
    source("R/Pensim1.R")

    ## Build the penalty matrix ----
    vals <- my_penmat1(obs.grid = t.mat, d = length(X_fd.list)+1,
                       nbasis = nbasis, n.order = n.order, deg = 2)

    # Run the functional model ----
    n.lam <- min(10, floor(exp(log(100)/(d+1))))
    fofreg <- Pensim1(data = data, lams = seq(0.01, 5.0,length.out = n.lam),
                      vals = vals, d = d, d.z = ncol(X.scal),
                      t.mat = t.mat, weights = weights)

    model <- fofreg$model
    lams <- fofreg$lambda

    # ## Run the model ----
    # library(glmnet)
    # tmp.model <- cv.glmnet(x = as.matrix(data[,-c(1:3)]),
    #                        y = data$output,
    #                        alpha = 0,
    #                        intercept = FALSE)
    # #print(tmp.model$lambda.min)
    #
    # model <- glmnet(x = as.matrix(data[,-c(1:3)]),
    #                 y = data$output, intercept = FALSE,
    #                 alpha = 0, lambda = tmp.model$lambda.min)

    ## get the functional parameters ----
    params <- model$coefficients #as.vector(coef(model))[-1]
    params[is.na(params)] <- 0
    beta.fd <- lapply(1:(d+1), function(l){
      if (l != 1){
        fd(coef = params[(nbasis*(l-1)+1):(nbasis*l)],
           basisobj = create.bspline.basis(rangeval = range(t.mat, na.rm = T),
                                           nbasis = nbasis, norder = n.order),
           fdnames = list(main = paste("beta_",(l-1), sep = ""),
                          xlab = "Time (millisec.)",
                          ylab = paste("beta_",(l-1),sep = "")))
      } else{
        fd(coef = params[(nbasis*(l-1)+1):(nbasis*l)],
           basisobj = create.bspline.basis(rangeval = range(t.mat, na.rm = T),
                                           nbasis = nbasis, norder = n.order),
           fdnames = list(main = paste("beta_",(l-1),sep = ""),
                          xlab = "Time (millisec.)",
                          ylab = paste("beta_",(l-1),sep = "")))
      }
    })

    if (ncol(X.scal) > 0) {
      ## get the non functional parameters ----
      beta.scal <- params[-c(1:(nbasis*(d+1)))]
      ## results ----
      return(list(model = model, beta.fd = beta.fd,
                  beta.scal = beta.scal, lambda = lams, n.grid = m))
    } else {
      ## results ----
      return(list(model = model, beta.fd = beta.fd, lambda = lams,
                  n.grid = m))
    }
  }
}

#' Prediction of Functional concurrent model
#' @param model the functional model
#' @param t.mat the matrix of time on where we want the prediction
#' @param newX_fd.list the list of functional predictors
#' @param newX.scal the data frame of the non functional predictors
#' @importFrom stats predict
#' @return the prediction of functional concurrent model
#' @export
pred.penffr1 <- function(model, newX_fd.list = list(), newX.scal = data.frame(), t.mat = "all.ok"){

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
    if (all(t.mat == "all.ok")) {
      obs.grid <- seq(newX_fd.list[[1]][["basis"]][["rangeval"]][1],
                      newX_fd.list[[1]][["basis"]][["rangeval"]][2],
                      length.out = model$n.grid)
      t.mat <- matrix(rep(obs.grid, n), nrow = n, byrow = T)
    }

    d <- length(newX_fd.list)

    # Build the dataframe and penalty matrix for the model ----
    source("R/func-to-mat1.R")
    nbasis <- length(as.vector(model$beta.fd[[1]][["coefs"]]))
    n.order <- nbasis - length(model$beta.fd[[1]][["basis"]][["params"]])
    data <- get.data1(X_fd.list = newX_fd.list,
                      Y.mat = t.mat, nbasis = nbasis,
                      t.mat = t.mat, n.order = n.order)

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




