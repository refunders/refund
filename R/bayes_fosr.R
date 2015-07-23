#' Bayesian Function-on-scalar regression
#' 
#' Wrapper function that implements several approaches to Bayesian function-
#' on-scalar regression. Currently handles real-valued response curves; models
#' can include subject-level random effects in a multilevel framework. The 
#' residual curve error structure can be estimated using Bayesian FPCA or a 
#' Wishart prior. Model parameters can be estimated using a Gibbs sampler
#' or variational Bayes.
#' 
#' @param formula a formula indicating the structure of the proposed model. 
#' Random intercepts are designated using \code{\link{re}}().
#' @param data an optional data frame, list or environment containing the 
#' variables in the model. If not found in data, the variables are taken from 
#' environment(formula), typically the environment from which the function is 
#' called.
#' @param est.method method used to estimate model parameters. Options are "VB",
#' "Gibbs", and "GLS" with "VB" as default. Variational Bayes is a fast approximation to
#' the full posterior and often provides good point estimates, but may be 
#' unreliable for inference. "GLS" doesn't do anything Bayesian -- just fits an
#' unpenalized GLS estimator for the specified model.
#' @param cov.method method used to estimate the residual covariance structure.
#' Options are "FPCA" and "Wishart", with default "FPCA"
#' @param ... additional arguments that are passed to individual fitting functions.
#' 
#' @references
#' Goldsmith, J., Kitago, T. (Under Review).
#' Assessing Systematic Effects of Stroke on Motor Control using Hierarchical 
#' Function-on-Scalar Regression.
#' 
#' @author Jeff Goldsmith \email{ajg2202@@cumc.columbia.edu}
#' @export
#' 
bayes_fosr = function(formula, data=NULL, est.method = "VB", cov.method = "FPCA", ...){
  
  ranef = sum(grepl("re", formula))

  if(ranef == 0 & est.method == "GLS"){
    ret = gls_cs(formula = formula, data = data, ...)
  } else if(ranef == 0 & est.method == "OLS"){
    ret = ols_cs(formula = formula, data = data, ...)
  } else if(ranef == 0 & est.method == "VB" & cov.method == "FPCA"){
    ret = vb_cs_fpca(formula = formula, data = data, ...)
  } else if(ranef == 0 & est.method == "VB" & cov.method == "Wishart"){
    ret = vb_cs_wish(formula = formula, data = data, ...)
  } else if(ranef == 0 & est.method == "Gibbs" & cov.method == "FPCA"){
    ret = gibbs_cs_fpca(formula = formula, data = data, ...)
  } else if(ranef == 0 & est.method == "Gibbs" & cov.method == "Wishart"){
    ret = gibbs_cs_wish(formula = formula, data = data, ...)
  } else if(ranef == 1 & est.method == "VB" & cov.method == "FPCA"){
    ret = vb_mult_fpca(formula = formula, data = data, ...)
  } else if(ranef == 1 & est.method == "VB" & cov.method == "Wishart"){
    ret = vb_mult_wish(formula = formula, data = data, ...)
  } else if(ranef == 1 & est.method == "Gibbs" & cov.method == "FPCA"){
    ret = gibbs_mult_fpca(formula = formula, data = data, ...)
  } else if(ranef == 1 & est.method == "Gibbs" & cov.method == "Wishart"){
    ret = gibbs_mult_wish(formula = formula, data = data, ...)
  } else {
    error("The combination of estimation method and covariance structure you specified is not yet implemented.")
  }
  
  ret

}

###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################