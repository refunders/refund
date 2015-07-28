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
#' Goldsmith, J., Kitago, T. (Accepted).
#' Assessing Systematic Effects of Stroke on Motor Control using Hierarchical 
#' Function-on-Scalar Regression.
#' Journal of the Royal Statistical Society: Series C.
#' 
#' @author Jeff Goldsmith \email{ajg2202@@cumc.columbia.edu}
#' @export
#' 
#' @examples
#' 
#' \dontrun{
#' 
#' library(reshape2)
#' library(dplyr)
#' library(ggplot2)
#' 
#' ##### Cross-sectional real-data examples #####
#' 
#' ## organize data
#' data(DTI)
#' DTI = subset(DTI, select = c(cca, case, pasat))
#' DTI = DTI[complete.cases(DTI),]
#' DTI$gender = factor(sample(c("male","female"), dim(DTI)[1], replace = TRUE))
#' DTI$status = factor(sample(c("RRMS", "SPMS", "PPMS"), dim(DTI)[1], replace = TRUE))
#' 
#' ## fit models
#' default = bayes_fosr(cca ~ pasat, data = DTI)
#' VB = bayes_fosr(cca ~ pasat, data = DTI, Kp = 4, Kt = 10)
#' Gibbs = bayes_fosr(cca ~ pasat, data = DTI, Kt = 10, est.method = "Gibbs", cov.method = "Wishart",
#'                    N.iter = 500, N.burn = 200)
#' OLS = bayes_fosr(cca ~ pasat, data = DTI, Kt = 10, est.method = "OLS")
#' GLS = bayes_fosr(cca ~ pasat, data = DTI, Kt = 10, est.method = "GLS")
#' 
#' ## plot results
#' models = c("default", "VB", "Gibbs", "OLS", "GLS")
#' intercepts = sapply(models, function(u) get(u)$beta.hat[1,])
#' slopes = sapply(models, function(u) get(u)$beta.hat[2,])
#' 
#' plot.dat = melt(intercepts); colnames(plot.dat) = c("grid", "method", "value")
#' ggplot(plot.dat, aes(x = grid, y = value, group = method, color = method)) + geom_path() + theme_bw()
#' 
#' plot.dat = melt(slopes); colnames(plot.dat) = c("grid", "method", "value")
#' ggplot(plot.dat, aes(x = grid, y = value, group = method, color = method)) + geom_path() + theme_bw()
#' 
#' ## fit a model with an interaction
#' fosr.dti.interaction = bayes_fosr(cca ~ pasat*gender, data = DTI, Kp = 4, Kt = 10)
#' 
#' 
#' ##### Longitudinal real-data examples #####
#' 
#' data(DTI2)
#' DTI2 = subset(DTI2, select = c(cca, id, pasat))
#' DTI2 = DTI2[complete.cases(DTI2),]
#' 
#' default = bayes_fosr(cca ~ pasat + re(id), data = DTI2)
#' VB = bayes_fosr(cca ~ pasat + re(id), data = DTI2, Kt = 10, cov.method = "Wishart")
#' 
#' }
#' 
bayes_fosr = function(formula, data=NULL, est.method = "VB", cov.method = "FPCA", ...){
  
  tf <- terms.formula(formula, specials = "re")
  trmstrings <- attr(tf, "term.labels")
  specials <- attr(tf, "specials")    # if there are no random effects this will be NULL
  where.re <-specials$re - 1
  ranef = length(where.re)
  
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
    stop("The combination of estimation method and covariance structure you specified is not yet implemented.")
  }
  
  ret

}

###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################