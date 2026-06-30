# Software Name : PenFFR
# SPDX-FileCopyrightText: Copyright (c) Orange SA
# SPDX-License-Identifier:  GPL-2.0-or-later

# This software is distributed under the GNU General Public License v2.0 or later
# see the "LICENSE.md" file for more details or https://spdx.org/licenses/GPL-2.0-or-later.html

# Authors: see CONTRIBUTORS.md
# Software description: Linear regression for fully functional data.

#' Functional Mixture-of-Experts model for concurrent model
#' @param Y.mat the matrix response with of n rows and m = max m_i columns with na values after the end of the observation for each sample
#' @param t.mat the matrix of time with the corresponded values of Y.mat
#' @param X_fd.list the list of functional predictors
#' @param X.scal the data frame of the non functional predictors
#' @param nbasis the number of basis for functional parameter (an integer)
#' @param n.order the order of the splines basis
#' @param K.min minimum number of components to consider
#' @param K.max maximum number of components to consider
#' @param pen the boolean value that indicates if we want to use the penalization
#' @param K.comp the number of fixed component to have or "BIC" to select the best model using BIC criterion
#' @param K.mixture the number of component mixture
#' @import fda
#' @import flexmix
#' @return a list that containt
#' \itemize{
#'   \item \code{model}: the fitting model object
#'   \item \code{beta.fd}: the estimated functional parameters.
#'   \item \code{beta.scal}: the estimated scalar parameters
#'   \item \code{lambda}: the optimal penalized parameter
#'   \item \code{n.grid}: the length of time observation (useful for the predict function)
#' }
#' @export
penffr1_mix <- function(Y.mat, t.mat = "all.ok", X_fd.list, X.scal = data.frame(), nbasis = "all.ok", n.order = 4, K.min = 1, K.max = 5, pen = F, K.comp = "BIC", seed = 2807){

  '
  Y.mat        : matrix response with of n rows and m = max m_i columns with na values after the end of the observation for each sample
  t.mat        : a matrix of time with the corresponded values of Y.mat
  X_fd.list    : list of functional predictors
  X.scal       : data frame of the non functional predictors
  nbasis       : an integer for the number of basis for functional parameter
  Pen          : boolean value for the penalization or not, default = True
  K.min        : Minimum number of mixture components
  K.max        : Maximum number of mixture components
  K.comp       : Number of components, default is the best according to BIC
  '

  # Check if the t.mat is not provided ----
  n <- nrow(Y.mat)
  d <- length(X_fd.list)
  if (all(t.mat == "all.ok")) {
    m <- ncol(Y.mat)
    obs.grid <- ((1:m)-1)/(m-1)
    t.mat <- matrix(rep(obs.grid, n), nrow = n, byrow = T)
  }
  m <- ncol(t.mat)

  # Check if nbasis is provided ----
  if (all(nbasis == "all.ok")){
    nbasis = min(20, floor((n*m)/(d+1)))
  }

  # Build the dataframe and penalty matrix for the model ----
  source("R/func-to-mat1.R")
  data <- get.data1(X_fd.list = X_fd.list,
                    Y.mat = Y.mat, nbasis = nbasis,
                    t.mat = t.mat, n.order = n.order)

  # Features of Gated Network function ----
  source("R/gated_features.R")
  data.gated <- gated.features(t.mat = t.mat, nbasis = nbasis,
                               X_fd.list = X_fd.list,
                               n.order = n.order)

  if (ncol(X.scal) > 0){
    ## Formatting the non-functional data ----
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
    data.gated <- cbind.data.frame(data.gated, X.scal.new)
  }


  if (pen == F){ # Non penalization case

    ## Model estimation ----
    Paste <- function(x, y = ""){paste(paste0(x, y), collapse = " + ")}
    tmp.formula <- sprintf("%s ~ 0 + %s | id",
                           "output",
                           Paste(colnames(data)[-(1:3)]))

    my_formula <- formula(tmp.formula)
    my_formula2 <- as.formula(paste("~", paste(colnames(data.gated),
                                               collapse = " + "),
                                    sep = ""))
    #tmp.model <- FLXMRglm(my_formula)
    tmp.model <- FLXMRlmm(random = ~ 1, lm.fit = "lm.wfit")
    control <- list(verbose = -0, iter.max = 2000, minprior = 0.1)

    if (K.comp == "BIC") {
      set.seed(seed)
      model <- initFlexmix(my_formula,
                           model = tmp.model,
                           concomitant = FLXPmultinom(my_formula2),
                           data = cbind(data, data.gated),
                           nrep = 5, k = K.min:K.max,
                           control = control, verbose = F)
      #print(model)
      m.BIC <- getModel(model, "BIC")

    } else {
      set.seed(seed)
      m.BIC <- flexmix(my_formula,
                       model = tmp.model,
                       k = K.comp,
                       concomitant = FLXPmultinom(my_formula2),
                       data = cbind(data, data.gated),
                       control = control)
    }


    ## get the functional parameters ----
    params <- parameters(m.BIC)
    params[is.na(params)] <- 0
    K <- ncol(params)
    #params <- params[-((nrow(params)-1):nrow(params)),]
    if (K == 1) {
      params <- matrix(params, ncol = 1)
    }
    beta.fd <- lapply(1:K, function(k){
      lapply(1:(d+1), function(l){
        fd(coef = params[(nbasis*(l-1)+1):(nbasis*l), k],
           basisobj = create.bspline.basis(rangeval = range(t.mat, na.rm = T),
                                           nbasis = nbasis, norder = 4),
           fdnames = list(main = paste("beta_",(l-1), sep = ""),
                          xlab = "Time",
                          ylab = paste("beta_",(l-1),sep = "")))
      })
    })

    if (ncol(X.scal) > 0) {
      ## get the non functional parameters ----
      beta.scal <- params[-c(1:(nbasis*(d+1))),]
      ## results ----
      return(list(model = m.BIC, beta.fd = beta.fd, beta.scal = beta.scal,
                  n.grid = m))
    } else {
      ## results ----
      return(list(model = m.BIC, beta.fd = beta.fd, n.grid = m))
    }

  } else { # Penalization case
    source("R/Pensim1_mix.R")

    ## Build the penalty matrix ----
    vals <- my_penmat1(obs.grid = t.mat, d = d+1,
                       nbasis = nbasis, deg = 2,
                       n.order = n.order)

    # Run the functional model ----
    n.lam <- min(10, floor(exp(log(100)/(d+1)))) + 2

    fofreg <- Pensim1_mix(data = data, lams = seq(0.01, 0.5, length.out = n.lam),
                          vals = vals, d = d, d.z = ncol(X.scal), seed = seed,
                          K.comp = K.comp, K.min = K.min, K.max = K.max, X.scal = X.scal,
                          t.mat = t.mat, data.gated = data.gated, n = n, m = ncol(t.mat))
    print(fofreg$lambda)
    m.BIC <- fofreg$model

    ## get the functional parameters ----
    params <- parameters(m.BIC)
    params[is.na(params)] <- 0
    K <- ncol(params)
    #params <- params[-((nrow(params)-1):nrow(params)),]
    if (K == 1) {
      params <- matrix(params, ncol = 1)
    }
    beta.fd <- lapply(1:K, function(k){
      lapply(1:(d+1), function(l){
        fd(coef = params[(nbasis*(l-1)+1):(nbasis*l), k],
           basisobj = create.bspline.basis(rangeval = range(t.mat, na.rm = T),
                                           nbasis = nbasis, norder = 4),
           fdnames = list(main = paste("beta_",(l-1), sep = ""),
                          xlab = "Time",
                          ylab = paste("beta_",(l-1),sep = "")))
      })
    })

    if (ncol(X.scal) > 0) {
      ## get the non functional parameters ----
      beta.scal <- params[-c(1:(nbasis*(d+1))),]
      ## results ----
      return(list(model = m.BIC, beta.fd = beta.fd, n.grid = m,
                  beta.scal = beta.scal, lambda = fofreg$lambda))
    } else {
      ## results ----
      return(list(model = m.BIC, beta.fd = beta.fd, n.grid = m,
                  lambda = fofreg$lambda))
    }
  }
}

#' Prediction of Functional Mixture-of-Experts model for concurrent model
#' @param model the functional Mixture-of-Experts model
#' @param t.mat the matrix of time on where we want the prediction
#' @param newX_fd.list the list of functional predictors
#' @param newX.scal the data frame of the non functional predictors
#' @importFrom stats predict
#' @return the prediction of functional Mixture-of-Experts model for concurrent model
#' @export
pred.penffr1_mix <- function(model, newX_fd.list = list(), newX.scal = data.frame(), t.mat = "all.ok"){

  # Check the providing of new functional predictors ----
  if (length(newX_fd.list) == 0) {
    if (nrow(newX.scal) == 0) {
      tmp.pred <- predict(model$model, type = "response")
      K <- length(tmp.pred)
      tmp <- matrix(0, ncol = K, nrow = length(unlist(tmp.pred[[1]])))
      for (k in 1:K) {
        tmp[,k] <- unlist(tmp.pred[[k]])
      }
      clust <- data.frame(x = clusters(model$model))
      clust$x <- as.factor(clust$x)
      pik <- data.frame(model.matrix(~x - 1, data = clust))
      pred <- rowSums(tmp * pik)

    } else {
      return(NA)
    }

  } else {

    # observation grid of predictions ----
    n <- ncol(newX_fd.list[[1]][["coefs"]])
    obs.grid <- seq(newX_fd.list[[1]][["basis"]][["rangeval"]][1],
                    newX_fd.list[[1]][["basis"]][["rangeval"]][2],
                    length.out = model$n.grid)
    if (all(t.mat == "all.ok")) {
      t.mat <- matrix(rep(obs.grid, n), nrow = n, byrow = T)
    }
    d <- length(newX_fd.list)

    source("R/func-to-mat1.R")
    nbasis <- length(as.vector(model$beta.fd[[1]][[1]][["coefs"]]))
    n.order <- nbasis - length(model$beta.fd[[1]][[1]][["basis"]][["params"]])
    data <- get.data1(X_fd.list = newX_fd.list,
                      Y.mat = t.mat, nbasis = nbasis,
                      t.mat = t.mat, n.order = n.order)

    # Features of Gated Network function ----
    source("R/gated_features.R")
    data.gated <- gated.features(t.mat = t.mat, nbasis = nbasis,
                                 X_fd.list = newX_fd.list,
                                 n.order = n.order)
    m <- ncol(t.mat)
    K <- ncol(parameters(model$model))

    if (ncol(newX.scal) > 0){
      ## Formatting the non-functional data ----
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
      data.gated <- cbind.data.frame(data.gated, X.scal.new)
    }

    # prediction ----
    params <- parameters(model$model)
    params <- params[-((nrow(params)-1):nrow(params)),]
    source("R/utils.R")
    params.gated <- model$model@concomitant@coef
    #params.gated <- params.gated/max(abs(params.gated))
    pred <- predict.moe(params = params, m=m,
                        params.gated = params.gated,
                        data = data, data.gated = data.gated)

    # result ----
    return(pred)
  }
}

