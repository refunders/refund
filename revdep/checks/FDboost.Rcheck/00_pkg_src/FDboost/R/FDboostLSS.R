#################################################################################
#' Model-based Gradient Boosting for Functional GAMLSS
#' 
#' Function for fitting generalized additive models for location, scale and shape (GAMLSS)  
#' with functional data using component-wise gradient boosting, for details see 
#' Brockhaus et al. (2018). 
#' 
#' @param formula a symbolic description of the model to be fit. 
#' If \code{formula} is a single formula, the same formula is used for all distribution parameters. 
#' \code{formula} can also be a (named) list, where each list element corresponds to one distribution 
#' parameter of the GAMLSS distribution. The names must be the same as in the \code{families}. 
#' @param timeformula one-sided formula for the expansion over the index of the response. 
#' For a functional response \eqn{Y_i(t)} typically \code{~bbs(t)} to obtain a smooth 
#' expansion of the effects along \code{t}. In the limiting case that \eqn{Y_i} is a scalar response
#' use \code{~bols(1)}, which sets up a base-learner for the scalar 1. 
#' Or you can use \code{timeformula=NULL}, then the scalar response is treated as scalar. 
#' Analogously to \code{formula}, \code{timeformula} can either be a one-sided formula or 
#' a named list of one-sided formulas. 
#' @param data a data frame or list containing the variables in the model.
#' @param families an object of class \code{families}. It can be either one of the pre-defined distributions 
#' that come along with the package \code{gamboostLSS} or a new distribution specified by the user 
#' (see \code{\link[gamboostLSS]{Families}} for details). 
#' Per default, the two-parametric \code{\link[gamboostLSS]{GaussianLSS}} family is used.
#' @param control  a list of parameters controlling the algorithm. 
#' For more details see \code{\link[mboost]{boost_control}}.  
#' @param weights does not work!
#' @param method fitting method, currently two methods are supported: 
#' \code{"cyclic"} (see Mayr et al., 2012) and \code{"noncyclic"} 
#' (algorithm with inner loss of Thomas et al., 2018).
#' @param ... additional arguments passed to \code{\link[FDboost]{FDboost}}, 
#' including, \code{family} and \code{control}.
#' 
#' @details For details on the theory of GAMLSS, see Rigby and Stasinopoulos (2005). 
#' \code{FDboostLSS} calls \code{FDboost} to fit the distribution parameters of a GAMLSS - 
#' a functional boosting model is fitted for each parameter of the response distribution.  
#' In \code{\link[gamboostLSS]{mboostLSS}}, details on boosting of GAMLSS based on 
#' Mayr et al. (2012) and Thomas et al. (2018) are given.   
#' In \code{\link{FDboost}}, details on boosting regression models with functional variables 
#' are given (Brockhaus et al., 2015, Brockhaus et al., 2017). 
#' 
#' @return An object of class \code{FDboostLSS} that inherits from \code{mboostLSS}. 
#' The \code{FDboostLSS}-object is a named list containing one list entry per distribution parameter
#' and some attributes. The list is named like the parameters, e.g. mu and sigma, 
#' if the parameters mu and sigma are modeled. Each list-element is an object of class \code{FDboost}.  
#' 
#' @author Sarah Brockhaus 
#' 
#' @seealso Note that \code{FDboostLSS} calls \code{\link{FDboost}} directly.  
#' 
#' @keywords models, nonlinear 
#' 
#' @references 
#' Brockhaus, S., Scheipl, F., Hothorn, T. and Greven, S. (2015). 
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
#' Mayr, A., Fenske, N., Hofner, B., Kneib, T. and Schmid, M. (2012): 
#' Generalized additive models for location, scale and shape for high-dimensional 
#' data - a flexible approach based on boosting. 
#' Journal of the Royal Statistical Society: Series C (Applied Statistics), 61(3), 403-427. 
#' 
#' Rigby, R. A. and D. M. Stasinopoulos (2005):  
#' Generalized additive models for location, scale and shape (with discussion). 
#' Journal of the Royal Statistical Society: Series C (Applied Statistics), 54(3), 507-554. 
#' 
#' Thomas, J., Mayr, A., Bischl, B., Schmid, M., Smith, A., and Hofner, B. (2018), 
#' Gradient boosting for distributional regression - faster tuning and improved 
#' variable selection via noncyclical updates. 
#' Statistics and Computing, 28, 673-687. 
#' 
#' @examples 
#' ########### simulate Gaussian scalar-on-function data
#' n <- 500 ## number of observations
#' G <- 120 ## number of observations per functional covariate
#' set.seed(123) ## ensure reproducibility
#' z <- runif(n) ## scalar covariate
#' z <- z - mean(z)
#' s <- seq(0, 1, l=G) ## index of functional covariate
#' ## generate functional covariate
#' if(require(splines)){
#'    x <- t(replicate(n, drop(bs(s, df = 5, int = TRUE) %*% runif(5, min = -1, max = 1))))
#' }else{
#'   x <- matrix(rnorm(n*G), ncol = G, nrow = n)
#' }
#' x <- scale(x, center = TRUE, scale = FALSE) ## center x per observation point
#' 
#' mu <- 2 + 0.5*z + (1/G*x) %*% sin(s*pi)*5 ## true functions for expectation
#' sigma <- exp(0.5*z - (1/G*x) %*% cos(s*pi)*2) ## for standard deviation
#' 
#' y <- rnorm(mean = mu, sd = sigma, n = n) ## draw respone y_i ~ N(mu_i, sigma_i)
#' 
#' ## save data as list containing s as well 
#' dat_list <- list(y = y, z = z, x = I(x), s = s)
#' 
#' ## model fit with noncyclic algorithm assuming Gaussian location scale model 
#' m_boost <- FDboostLSS(list(mu = y ~ bols(z, df = 2) + bsignal(x, s, df = 2, knots = 16), 
#'                            sigma = y ~ bols(z, df = 2) + bsignal(x, s, df = 2, knots = 16)), 
#'                            timeformula = NULL, data = dat_list, method = "noncyclic")
#' summary(m_boost)
#' 
#' \dontrun{
#'  if(require(gamboostLSS)){
#'   ## find optimal number of boosting iterations on a grid in 1:1000
#'   ## using 5-fold bootstrap
#'   ## takes some time, easy to parallelize on Linux
#'   set.seed(123) 
#'   cvr <- cvrisk(m_boost, folds = cv(model.weights(m_boost[[1]]), B = 5),
#'                 grid = 1:1000, trace = FALSE)
#'   ## use model at optimal stopping iterations 
#'   m_boost <- m_boost[mstop(cvr)] ## 832
#'    
#'   ## plot smooth effects of functional covariates for mu and sigma
#'   par(mfrow = c(1,2))
#'   plot(m_boost$mu, which = 2, ylim = c(0,5))
#'   lines(s, sin(s*pi)*5, col = 3, lwd = 2)
#'   plot(m_boost$sigma, which = 2, ylim = c(-2.5,2.5))
#'   lines(s, -cos(s*pi)*2, col = 3, lwd = 2)
#'  }
#' }
#' @export
## function that calls FDboost for each distribution parameter
FDboostLSS <- function(formula, timeformula, data = list(), families = GaussianLSS(),
                       control = boost_control(), weights = NULL, 
                       method = c("cyclic", "noncyclic"), ...){
  
  cl <- match.call()
  if(is.null(cl$families)) cl$families <- families
  
  ## warnings for functional response are irrelevant for scalar response  
  if( !is.null(timeformula) && timeformula != ~bols(1) ){
    message("No smooth offsets over time are used, just global scalar offsets.")
    message("No integration weights are used to compute the loss for the functional response.")
  }
  method <- match.arg(method)
  
  fit <- mboostLSS_fit(formula = formula, timeformula = timeformula, 
                       data = data, families = families,
                       control = control, weights = weights, ...,
                       fun = FDboost, funchar = "FDboost", call = cl, 
                       method = method)

  ## make sure that the first class of the model object is 'FDboostLSS'
  class(fit) <- class(fit)[class(fit) != "FDboostLSS"]
  class(fit) <- c("FDboostLSS", class(fit))
  
  return(fit)
}


#################################################################################
#' Cross-validation for FDboostLSS
#' 
#' Multidimensional cross-validated estimation of the empirical risk for hyper-parameter selection, 
#' for an object of class \code{FDboostLSS} setting the folds per default to resampling curves.  
#' 
#' @param object an object of class \code{FDboostLSS}. 
#' @param folds a weight matrix a weight matrix with number of rows equal to the number of observations. 
#' The number of columns corresponds to the number of cross-validation runs, 
#' defaults to 25 bootstrap samples, resampling whole curves  
#' @param grid defaults to a grid up to the current number of boosting iterations. 
#' The default generates the grid according to the defaults of 
#' \code{\link[gamboostLSS]{cvrisk.mboostLSS}} and \code{\link[gamboostLSS]{cvrisk.nc_mboostLSS}} for
#' models with cyclic or noncyclic fitting.  
#' @param papply (parallel) apply function, defaults to \code{\link[parallel]{mclapply}}, 
#' see \code{\link[gamboostLSS]{cvrisk.mboostLSS}} for details 
#' @param trace print status information during cross-validation? Defaults to \code{TRUE}.
#' @param fun if \code{fun} is \code{NULL}, the out-of-sample risk is returned. 
#' \code{fun}, as a function of \code{object}, 
#' may extract any other characteristic of the cross-validated models. These are returned as is.
#' @param ... additional arguments passed to \code{\link[parallel]{mclapply}}.
#' 
#' @details The function \code{cvrisk.FDboostLSS} is a wrapper for 
#' \code{\link[gamboostLSS]{cvrisk.mboostLSS}} in package gamboostLSS.  
#' It overrieds the default for the folds, so that the folds are sampled on the level of curves 
#' (not on the level of single observations, which does not make sense for functional response).  
#' 
#' @return An object of class \code{cvriskLSS} (when \code{fun} was not specified), 
#' basically a matrix containing estimates of the empirical risk for a varying number 
#' of bootstrap iterations. \code{plot} and \code{print} methods are available as well as an 
#' \code{mstop} method, see \code{\link[gamboostLSS]{cvrisk.mboostLSS}}.
#' 
#' @seealso \code{\link[gamboostLSS]{cvrisk.mboostLSS}} in packge gamboostLSS. 
#' 
#' @export
## wrapper for cvrisk of gamboostLSS, specifying folds on the level of curves
cvrisk.FDboostLSS <- function(object, folds = cvLong(id = object[[1]]$id, 
                                                     weights = model.weights(object[[1]])),
                              grid = NULL,
                              papply = mclapply, trace = TRUE, 
                              fun = NULL, ...){
  
  ## message not necessary as currently only a scalar offset is possible for FDboostLSS-models
  ## if(!length(unique(object$offset)) == 1) message("The smooth offset is fixed over all folds.")
  
  class(object) <- class(object)[class(object) != "FDboostLSS"]
  
  ## set up grid according to defaults of cvrisk.nc_mboostLSS and cvrisk.mboostLSS
  if(is.null(grid)){
    
    if(any(class(object) == "nc_mboostLSS")){
      grid <- 1:sum(mstop(object))
    }else{
      grid <- make.grid(mstop(object))
    }
    
  }

  ## call cvrisk.nc_mboostLSS or cvrisk.mboostLSS
  ret <- cvrisk(object = object, folds = folds,
                grid = grid,
                papply = papply, trace = trace, 
                fun = fun, ...)
  return(ret) 
}

