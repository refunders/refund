# Software Name : PenFFR
# SPDX-FileCopyrightText: Copyright (c) Orange SA
# SPDX-License-Identifier:  GPL-2.0-or-later

# This software is distributed under the GNU General Public License v2.0 or later
# see the "LICENSE.md" file for more details or https://spdx.org/licenses/GPL-2.0-or-later.html

# Authors: see CONTRIBUTORS.md
# Software description: Linear regression and Mixture-of-Experts Linear Regression for fully functional data.

#' The functional model (PenFFR and PenFFMoE)
#' @param Y.mat the matrix response with of n rows and m = max m_i columns with na values after the end of the observation for each sample
#' @param t.mat the matrix of time with the corresponded values of Y.mat
#' @param X_fd.list the list of functional predictors
#' @param X.scal the data frame of the non functional predictors
#' @param nbasis the number of basis for functional parameter (an integer)
#' @param model.type the type of model we using and can be "concurrent" or "integral"
#' @param pen the boolean True/False value that indicates if we want to use the penalization
#' @param K.mixture the number of component mixture
#' @description
#' The generic function designed to perform the (penalised) function-on-function linear model in the concurrent and integral form named PenFFR.
#' And the (penalized) function-on-function Mixture-of-Experts linear regression named (PenFFMoE)
#'
#' @return a list that containt
#' \itemize{
#'   \item \code{model}: the fitting model object
#'   \item \code{beta.fd}: the list of estimated functional parameters.
#'   \item \code{beta.scal}: the vector/matrix of estimated scalar parameters
#'   \item \code{lambda}: the optimal penalized parameter
#'   \item \code{n.grid}: the length of time observation (useful for the predict function)
#' }
#' @export

library(fda)
library(flexmix, exclude = c("fitted", "predict", "logLik"))
library(refund)
library(latex2exp)
library(transport)
library(mvtnorm)

PenFFR <- function(Y.mat, t.mat = "all.ok", X_fd.list, X.scal = data.frame(), nbasis = "all.ok", model.type = "concurrent", K.mixture = "No", pen = T){


  if (K.mixture == "No") {

    if (model.type == "concurrent") {
      source("R/penffr1.R")
      fofreg <- penffr1(Y.mat = Y.mat,
                        t.mat = t.mat,
                        X_fd.list = X_fd.list,
                        X.scal = X.scal,
                        nbasis = nbasis,
                        pen = pen)
    } else {
      source("R/penffr2.R")
      fofreg <- penffr2(Y.mat = Y.mat,
                        t.mat = t.mat,
                        X_fd.list = X_fd.list,
                        X.scal = X.scal,
                        nbasis = nbasis,
                        pen = pen)
    }

  } else if (K.mixture == "BIC"){

    if (model.type == "concurrent") {
      source("R/penffr1_mix.R")
      fofreg <- penffr1_mix(Y.mat = Y.mat,
                            t.mat = t.mat,
                            X_fd.list = X_fd.list,
                            X.scal = X.scal,
                            nbasis = nbasis,
                            K.comp = "BIC",
                            pen = pen)
    } else {
      source("R/penffr2_mix.R")
      fofreg <- penffr2_mix(Y.mat = Y.mat,
                            t.mat = t.mat,
                            X_fd.list = X_fd.list,
                            X.scal = X.scal,
                            nbasis = nbasis,
                            K.comp = "BIC",
                            pen = pen)
    }

  } else {

    K.mixture <- as.integer(readline(prompt = "Enter the number of components:"))

    if (model.type == "concurrent") {
      source("R/penffr1_mix.R")
      fofreg <- penffr1_mix(Y.mat = Y.mat,
                            t.mat = t.mat,
                            X_fd.list = X_fd.list,
                            X.scal = X.scal,
                            nbasis = nbasis,
                            K.comp = K.mixture,
                            pen = pen)
    } else {
      source("R/penffr2_mix.R")
      fofreg <- penffr2_mix(Y.mat = Y.mat,
                            t.mat = t.mat,
                            X_fd.list = X_fd.list,
                            X.scal = X.scal,
                            nbasis = nbasis,
                            K.comp = K.mixture,
                            pen = pen)
    }

  }

}

#' Prediction of Functional model
#' @param model the functional model
#' @param t.mat the matrix of time on where we want the prediction
#' @param newX_fd.list the list of functional predictors
#' @param newX.scal the data frame of the non functional predictors
#' @return the prediction of functional model
#' @export
pred.PenFFR <- function(model, newX_fd.list = list(), newX.scal = data.frame(), t.mat = "all.ok",
                        K.mixture = "No", model.type = "concurrent"){

  '
  t.mat must be the time matrix of predictions with of n rows and m columns
  newX_fd.list is the list of functional predictors
  newX.scal is the data frame of the non functional predictors
  model is the penffr model
  '

  if (K.mixture == "No") {

    if (model.type == "concurrent") {
      source("R/penffr1.R")
      fofreg <- pred.penffr1(model = model,
                             newX_fd.list = newX_fd.list,
                             newX.scal = newX.scal,
                             t.mat = t.mat)
    } else {
      source("R/penffr2.R")
      fofreg <- pred.penffr2(model = model,
                             newX_fd.list = newX_fd.list,
                             newX.scal = newX.scal,
                             t.mat = t.mat)
    }

  } else if (K.mixture == "BIC"){

    if (model.type == "concurrent") {
      source("R/penffr1_mix.R")
      fofreg <- pred.penffr1_mix(model = model,
                                 newX_fd.list = newX_fd.list,
                                 newX.scal = newX.scal,
                                 t.mat = t.mat)
    } else {
      source("R/penffr2_mix.R")
      fofreg <- pred.penffr2_mix(model = model,
                                 newX_fd.list = newX_fd.list,
                                 newX.scal = newX.scal,
                                 t.mat = t.mat)
    }

  } else {

    K.mixture <- as.integer(readline(prompt = "Enter the number of components:"))

    if (model.type == "concurrent") {
      source("R/penffr1_mix.R")
      fofreg <- pred.penffr1_mix(model = model,
                                 newX_fd.list = newX_fd.list,
                                 newX.scal = newX.scal,
                                 t.mat = t.mat)
    } else {
      source("R/penffr2_mix.R")
      fofreg <- pred.penffr2_mix(model = model,
                                 newX_fd.list = newX_fd.list,
                                 newX.scal = newX.scal,
                                 t.mat = t.mat)
    }

  }

}
