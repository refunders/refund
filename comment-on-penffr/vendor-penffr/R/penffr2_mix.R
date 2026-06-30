# Software Name : PenFFR
# SPDX-FileCopyrightText: Copyright (c) Orange SA
# SPDX-License-Identifier:  GPL-2.0-or-later

# This software is distributed under the GNU General Public License v2.0 or later
# see the "LICENSE.md" file for more details or https://spdx.org/licenses/GPL-2.0-or-later.html

# Authors: see CONTRIBUTORS.md
# Software description: Linear regression for fully functional data.

#' Functional Mixture-of-Experts model for integral model
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
#' @return a list that containt
#' \itemize{
#'   \item \code{model}: the fitting model object
#'   \item \code{beta.fd}: the estimated functional parameters.
#'   \item \code{beta.scal}: the estimated scalar parameters
#'   \item \code{lambda}: the optimal penalized parameter
#'   \item \code{n.grid}: the length of time observation (useful for the predict function)
#' }
#' @export
penffr2_mix <- function(Y.mat, t.mat = "all.ok", X_fd.list, X.scal = data.frame(), nbasis = "all.ok", K.comp = "BIC", pen = T){

  '
  Y.mat        : matrix response with of n rows and m = max m_i columns with na values after the end of the observation for each sample
  t.mat        : a matrix of time with the corresponded values of Y.mat
  X_fd.list    : list of functional predictors
  X.scal       : data frame of the non functional predictors
  nbasis       : an integer for the number of basis for functional parameter
  Pen          : boolean value for the penalization or not, default = True
  K.max        : Maximum number of mixture components, default = best based on BIC
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
                    t.mat = t.mat)

  # Features of Gated Network function ----
  source("R/gated_features.R")
  data.gated <- gated.features(t.mat = t.mat, nbasis = nbasis,
                               X_fd.list = X_fd.list)

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

  if (pen == F){

    K.max <- 5
    n.repeat = 3
    ## Run the model ----

    ### Run the non MoE model to detect no useful covariates
    my_formula <- as.formula(paste("output ~ 0", paste(colnames(data)[-c(1:3)],
                                                 collapse = " + "),
                                   sep = " + "))
    model <- lm(my_formula, data = data)
    params0 <- model$coefficients
    params0[is.na(params0)] <- 0
    ind = which(params0 == 0) + 3

    ### Now we run the MoE model with only useful covariates
    my_formula <- as.formula(paste(" ~ 0", paste(colnames(data)[-c(1:3, ind)],
                                                 collapse = " + "),
                                   sep = " + "))
    tmp.model <- FLXMRglmfix(my_formula)
    my_formula2 <- as.formula(paste(" ~ ", paste(colnames(data.gated),
                                                 collapse = " + "),
                                    sep = ""))

    if (K.comp == "BIC") {
      set.seed(999)
      model <- stepFlexmix(output ~ 0 | id,
                           model = tmp.model,
                           concomitant = FLXPmultinom(my_formula2),
                           data = cbind(data, data.gated),
                           nrep = n.repeat, k = 1:K.max)

      # select the best number of components by AIC
      new.d <- ncol(data) - 3
      loglik <- apply(model@logLiks, 1, max)
      my.BIC <- sapply(1:length(loglik), function(j){
        ((-2)*loglik[j]) + (2*(new.d*(model@k[j])*(model@k[j]+1))*(n/(n-(new.d*model@k[j])-1)))
      })


      m.BIC <- getModel(model, "BIC")#which = which.min(my.BIC))

    } else {
      set.seed(999)
      m.BIC <- flexmix(output ~ 0 | id,
                       model = tmp.model,
                       control = list(tolerance=1e100, classify="CEM"),
                       concomitant = FLXPmultinom(my_formula2),
                       data = cbind(data, data.gated),
                       k = K.comp)
    }



    ## get the functional parameters ----
    ### Parameters with the useful covariates
    params1 <- parameters(m.BIC)
    params1[is.na(params1)] <- 0
    params1 <- params1[-nrow(params1),]
    ### rearrange the parameter vector
    d.new <- ncol(data)-3
    K <- ncol(params1)
    params <- matrix(rep(0, d.new*K), ncol = K)
    params[-ind,] <- params1

    ## Recover the functional parameter
    beta.fd <- lapply(1:K, function(k){
      lapply(1:(d+1), function(l){
        if (l == 1){
          fd(coef = params[1:nbasis, k],
             basisobj = create.bspline.basis(rangeval = range(t.mat, na.rm = T),
                                             nbasis = nbasis, norder = 4),
             fdnames = list(main = "beta_0", xlab = "Time",
                            ylab = "beta_0"))
        } else{
          bifd(coef = matrix(params[nbasis+(((nbasis**2)*(l-2)+1):((nbasis**2)*(l-1))), k],
                             nbasis, nbasis, byrow = F),
               sbasisobj = create.bspline.basis(rangeval = range(obs.grid, na.rm = T),
                                                nbasis = nbasis, norder = 4),
               tbasisobj = create.bspline.basis(rangeval = range(obs.grid, na.rm = T),
                                                nbasis = nbasis, norder = 4))
        }
      })
    })

    if (ncol(X.scal) > 0) {
      ## get the non functional parameters ----
      beta.scal <- params[-c(1:(nbasis + (nbasis**2)*d)),]
      ## results ----
      return(list(model = m.BIC, beta.fd = beta.fd, beta.scal = beta.scal, n.grid = m))
    } else {
      ## results ----
      return(list(model = m.BIC, beta.fd = beta.fd, n.grid = m))
    }

  } else {

    source("R/func-to-mat2.R")
    ## Build the penalty matrix ----
    vals <- my_penmat2(obs.grid = t.mat, d = length(X_fd.list)+1,
                       L.beta = nbasis)

    # Run the penalized MoE functional model ----

    ### Run the non mixture model to get the useful covariates
    my_formula <- as.formula(paste("output ~ 0", paste(colnames(data)[-c(1:3)],
                                                 collapse = " + "),
                                   sep = " + "))
    model <- lm(my_formula, data = cbind.data.frame(data))
    params0 <- model$coefficients
    params0[is.na(params0)] <- 0
    ind = which(params0 == 0) + 3

    # Run the functional model ----
    source("R/Pensim2_mix.R")
    d <- length(X_fd.list)
    n.lam <- min(3, floor(exp(log(9)/(d+1))))
    fofreg <- Pensim2_mix(data = data, data.gated = data.gated, ind = ind,
                          lams = seq(10, 50, length.out = n.lam),
                          vals = vals, d = d, d.z = ncol(X.scal), K.comp = K.comp,
                          t.mat = t.mat)
    print(fofreg$lambda)
    m.BIC <- fofreg$model

    ## get the functional parameters ----
    ### Parameters on useful covariates
    params1 <- parameters(m.BIC)
    params1[is.na(params1)] <- 0
    params1 <- params1[-nrow(params1),]
    ### rearrange the parameter vector
    d.new <- ncol(data)-3
    K <- ncol(params1)
    params <- matrix(rep(0, d.new*K), ncol = K)
    params[-ind,] <- params1
    beta.fd <- lapply(1:K, function(k){
      lapply(1:(d+1), function(l){
        if (l == 1){
          fd(coef = params[1:nbasis, k],
             basisobj = create.bspline.basis(rangeval = range(t.mat, na.rm = T),
                                             nbasis = nbasis, norder = 4),
             fdnames = list(main = "beta_0", xlab = "Time",
                            ylab = "beta_0"))
        } else{
          bifd(coef = matrix(params[nbasis+(((nbasis**2)*(l-2)+1):((nbasis**2)*(l-1))), k],
                             nbasis, nbasis, byrow = F),
               sbasisobj = create.bspline.basis(rangeval = range(obs.grid, na.rm = T),
                                                nbasis = nbasis, norder = 4),
               tbasisobj = create.bspline.basis(rangeval = range(obs.grid, na.rm = T),
                                                nbasis = nbasis, norder = 4))
        }
      })
    })

    if (ncol(X.scal) > 0) {
      ## get the non functional parameters ----
      beta.scal <- params[-c(1:(nbasis + (nbasis**2)*d)),]
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

#' Prediction of Functional Mixture-of-Experts model for integral model
#' @param model the functional Mixture-of-Experts model
#' @param t.mat the matrix of time on where we want the prediction
#' @param newX_fd.list the list of functional predictors
#' @param newX.scal the data frame of the non functional predictors
#' @return the vector of predictions
#' @export
pred.penffr2_mix <- function(model, newX_fd.list = list(), newX.scal = data.frame(), t.mat = "all.ok"){

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

    nbasis <- length(as.vector(model$beta.fd[[1]][[1]][["coefs"]]))
    source("R/func-to-mat2.R")
    data <- get.data2(X_fd.list = newX_fd.list,
                      Y.mat = t.mat, nbasis = nbasis,
                      t.mat = t.mat)

    # Features of Gated Network function ----
    source("R/gated_features.R")
    data.gated <- gated.features(t.mat = t.mat, nbasis = nbasis,
                                 X_fd.list = newX_fd.list)
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
    params <- c()
    for (k in 1:K) {
      tmp <- c()
      for (l in 1:length(model$beta.fd[[k]])) {
        tmp <- c(tmp, as.vector(model$beta.fd[[k]][[l]][["coefs"]]))
      }
      params <- cbind(params, tmp)
    }

    if (ncol(newX.scal) > 0) {
      params <- rbind(params, model$beta.scal)
    }
    source("R/utils.R")
    pred <- predict.moe(params = params, m=m,
                        params.gated = model$model@concomitant@coef,
                        data = data, data.gated = data.gated)

    # result ----
    return(pred)
  }
}
