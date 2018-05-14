#################################################################################
#' FDboost: Boosting Functional Regression Models
#' 
#' @description 
#' Regression models for functional data, i.e., scalar-on-function, 
#' function-on-scalar and function-on-function regression models, are fitted 
#' by a component-wise gradient boosting algorithm.
#' 
#' @details 
#' This package is intended to fit regression models with functional variables. 
#' It is possible to fit models with functional response and/or functional covariates, 
#' resulting in scalar-on-function, function-on-scalar and function-on-function regression. 
#' Details on the functional regression models that can be fitted with \pkg{FDboost}
#' can be found in Brockhaus et al. (2015, 2017, 2018) and Ruegamer et al. (2018).  
#' A hands-on tutorial for the package can be found 
#' in Brockhaus, Ruegamer and Greven (2017), see \url{https://arxiv.org/abs/1705.10662}.
#' 
#' Using component-wise gradient boosting as fitting procedure, \pkg{FDboost} relies on 
#' the R package \pkg{mboost} (Hothorn et al., 2017). 
#' A comprehensive tutorial to \pkg{mboost} is given in Hofner et al. (2014). 
#' 
#' The main fitting function is \code{\link{FDboost}}. 
#' The model complexity is controlled by the number of boosting iterations (mstop).   
#' Like the fitting procedures in \pkg{mboost}, the function \code{FDboost} DOES NOT 
#' select an appropriate stopping iteration. This must be chosen by the user.  
#' The user can determine an adequate stopping iteration by resampling methods like 
#' cross-validation or bootstrap. 
#' This can be done using the function \code{\link{applyFolds}}. 
#' 
#' @references 
#' Brockhaus, S., Ruegamer, D. and Greven, S. (2017):
#' Boosting Functional Regression Models with FDboost.
#' \url{https://arxiv.org/abs/1705.10662}
#' 
#' Brockhaus, S., Scheipl, F., Hothorn, T. and Greven, S. (2015): 
#' The functional linear array model. Statistical Modelling, 15(3), 279-300. 
#' 
#' Brockhaus, S., Melcher, M., Leisch, F. and Greven, S. (2017): 
#' Boosting flexible functional regression models with a high number of functional historical effects,  
#' Statistics and Computing, 27(4), 913-926.   
#' 
#' Brockhaus, S., Fuest, A., Mayr, A. and Greven, S. (2018): 
#' Signal regression models for location, scale and shape with an application to stock returns. 
#' Journal of the Royal Statistical Society: Series C (Applied Statistics), 67, 665-686.
#' 
#' Hothorn T., Buehlmann P., Kneib T., Schmid M., and Hofner B. (2017). mboost: Model-Based Boosting, 
#' R package version 2.8-1, \url{https://cran.r-project.org/package=mboost}
#' 
#' Hofner, B., Mayr, A., Robinzonov, N., Schmid, M. (2014). Model-based Boosting in R: 
#' A Hands-on Tutorial Using the R Package mboost. Computational Statistics, 29, 3-35. 
#' \url{https://cran.r-project.org/package=mboost/vignettes/mboost_tutorial.pdf}
#' 
#' Ruegamer D., Brockhaus, S., Gentsch K., Scherer, K., Greven, S. (2018). 
#' Boosting factor-specific functional historical models for the detection of synchronization in bioelectrical signals. 
#' Journal of the Royal Statistical Society: Series C (Applied Statistics), 67, 621-642.
#' 
#' @author 
#' Sarah Brockhaus and David Ruegamer 
#' 
#' @aliases FDboost_package package-FDboost FDboost-package
#' 
#' @seealso 
#' \code{\link{FDboost}} for the main fitting function and 
#' \code{\link{applyFolds}} for model tuning via resampling methods. 
#' 
"_PACKAGE"
#> [1] "_PACKAGE"