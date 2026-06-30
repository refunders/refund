# Software Name : PenFFR
# SPDX-FileCopyrightText: Copyright (c) Orange SA
# SPDX-License-Identifier:  GPL-2.0-or-later

# This software is distributed under the GNU General Public License v2.0 or later
# see the "LICENSE.md" file for more details or https://spdx.org/licenses/GPL-2.0-or-later.html

# Authors: see CONTRIBUTORS.md
# Software description: Linear regression for fully functional data.
Pensim2 <- function(data, lams, vals, t.mat, d, d.z, weights){

  '
  data is the data frame of functional data expand as matrix
  data.gated is the data frame for the gated features
  lams is a vector of possible value of penalty parameters
  vals is the list of penalty matrix for every parameters
  t.mat is a n*m matrix of observation
  d is the number of functional predictors
  d.z is the number of scalar predictors
  '

  ## grid search of the optimum lambdas ----
  grid_search <- as.matrix(expand.grid(rep(list(lams), d+1)))
  cvals <- lapply(1:length(vals), function(l){
    tryCatch({
      chol(vals[[l]])
      },
      error = function(e){
        tmp <- vals[[l]]
        tmp[lower.tri(tmp)] <- 0
        return(tmp)
      })
  })

  source("R/utils.R")
  mspe.vect <- sapply(1:nrow(grid_search), function(v){

    if (v==1) {
      ## Run the model ----
      my_formula <- as.formula(paste("output ~ 0", paste(colnames(data)[-(1:3)],
                                                         collapse = " + "), sep = " + "))
      tmp.model <- lm(my_formula, data = cbind.data.frame(data), weights = weights)

    } else{
      ### Build the covariates and response for penalization ----
      datapen <- (grid_search[v,1]**2) * cvals[[1]]
      for (l in 2:(d+1)) {
        tmp <- grid_search[v,l] * cvals[[l]]
        datapen <- rbind(cbind(datapen, matrix(0, nrow = nrow(datapen), ncol = ncol(tmp))),
                         cbind(matrix(0, nrow = nrow(tmp), ncol = ncol(datapen)), tmp))
      }
      for (i in 1:3) {
        datapen <- cbind.data.frame(0, datapen)
      }
      if (d.z > 0) {
        for (i in 1:d.z) {
          datapen <- cbind.data.frame(datapen, 0)
        }
      }
      colnames(datapen) <- colnames(data)
      datapen$id <- as.factor(datapen$id)
      newdat <- rbind.data.frame(data, datapen)

      ### model estimation ----
      my_formula <- as.formula(paste("output ~ 0", paste(colnames(newdat)[-(1:3)],
                                                         collapse = " + "), sep = " + "))
      tmp.model <- lm(my_formula, data = data.frame(newdat), weights = weights)
    }

    ### Compute the BIC
    BIC(tmp.model)

  })

  ind <- which.min(mspe.vect)
  lams <- as.vector(grid_search[ind,])

  ### Build the covariates and response for penalization ----
  datapen <- grid_search[ind,1] * cvals[[1]]
  for (l in 2:(d+1)) {
    tmp <- grid_search[ind,l] * cvals[[l]]
    datapen <- rbind(cbind(datapen, matrix(0, nrow = nrow(datapen), ncol = ncol(tmp))),
                     cbind(matrix(0, nrow = nrow(tmp), ncol = ncol(datapen)), tmp))
  }
  for (i in 1:3) {
    datapen <- cbind.data.frame(0, datapen)
  }
  if (d.z > 0) {
    for (i in 1:d.z) {
      datapen <- cbind.data.frame(datapen, 0)
    }
  }
  colnames(datapen) <- colnames(data)
  datapen$id <- as.factor(datapen$id)
  newdat <- rbind.data.frame(data, datapen)

  ### model estimation ----
  my_formula <- as.formula(paste("output ~ 0", paste(colnames(newdat)[-(1:3)],
                                                     collapse = " + "), sep = " + "))
  model <- lm(my_formula, data = data.frame(newdat), weights = weights)

  return(list(lambda = lams, model = model))
}
