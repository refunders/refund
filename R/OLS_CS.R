#' Cross-sectional FoSR using GLS
#' 
#' Fitting function for function-on-scalar regression for cross-sectional data.
#' This function estimates model parameters using GLS: first, an OLS estimate of 
#' spline coefficients is estimated; second, the residual covariance is estimated
#' using an FPC decomposition of the OLS residual curves; finally, a GLS estimate
#' of spline coefficients is estimated. Although this is in the `BayesFoSR` package,
#' there is nothing Bayesian about this FoSR.
#' 
#' @param formula a formula indicating the structure of the proposed model. 
#' @param Kt number of spline basis functions used to estimate coefficient functions
#' @param data an optional data frame, list or environment containing the 
#' variables in the model. If not found in data, the variables are taken from 
#' environment(formula), typically the environment from which the function is 
#' called.
#' @param basis basis type; options are "bs" for b-splines and "pbs" for periodic
#' b-splines
#' @param verbose logical defaulting to \code{TRUE} -- should updates on progress be printed?
#'  
#' @references
#' Goldsmith, J., Kitago, T. (2016).
#' Assessing Systematic Effects of Stroke on Motor Control using Hierarchical 
#' Function-on-Scalar Regression. \emph{Journal of the Royal Statistical Society:
#' Series C}, 65 215-236.
#' 
#' @author Jeff Goldsmith \email{ajg2202@@cumc.columbia.edu}
#' @importFrom splines bs
#' @importFrom pbs pbs
#' @export
#' 
ols_cs = function(formula, data=NULL, Kt=5, basis = "bs", verbose = TRUE){
  
  call <- match.call()
  tf <- terms.formula(formula, specials = "re")
  trmstrings <- attr(tf, "term.labels")
  specials <- attr(tf, "specials")
  where.re <-specials$re - 1
  if (length(where.re) != 0) {
    mf_fixed <- model.frame(tf[-where.re], data = data)
    formula = tf[-where.re]
    responsename <- attr(tf, "variables")[2][[1]]
    ###
    REs = list(NA, NA)
    REs[[1]] = names(eval(parse(text=attr(tf[where.re], "term.labels")), envir=data)$data)
    REs[[2]]=paste0("(1|",REs[[1]],")")
    ###
    formula2 <- paste(responsename, "~", REs[[1]], sep = "")
    newfrml <- paste(responsename, "~", REs[[2]], sep = "")
    newtrmstrings <- attr(tf[-where.re], "term.labels")
    formula2 <- formula(paste(c(formula2, newtrmstrings), 
                              collapse = "+"))
    newfrml <- formula(paste(c(newfrml, newtrmstrings), collapse = "+"))
    mf <- model.frame(formula2, data = data)
    if (length(data) == 0) {
      Z = lme4::mkReTrms(lme4::findbars(newfrml), fr = mf)$Zt
    }
    else {
      Z = lme4::mkReTrms(lme4::findbars(newfrml), fr = data)$Zt
    }
  }
  else {
    mf_fixed <- model.frame(tf, data = data)
  }
  mt_fixed <- attr(mf_fixed, "terms")
  
  # get response (Y)
  Y <- model.response(mf_fixed, "numeric")
  
  # x is a matrix of fixed effects
  # automatically adds in intercept
  X <- model.matrix(mt_fixed, mf_fixed, contrasts)
  
  ### model organization ###
  D = dim(Y)[2]
  I = dim(X)[1]
  p = dim(X)[2]
  
  if(basis == "bs"){
    Theta = bs(1:D, df = Kt, intercept=TRUE, degree=3)
  } else if(basis == "pbs"){
    Theta = pbs(1:D, df = Kt, intercept=TRUE, degree=3)
  }

  X.des = X
  Y.vec = as.vector(t(Y)) 
  X = kronecker(X.des, Theta)
  n.coef = dim(X.des)[2]
  
  ## OLS model fitting and processing results
  if(verbose) { cat("Using OLS to estimate model parameters \n") }
  model.ols = lm(Y.vec ~ -1 + X)
  Bx.ols = matrix(model.ols$coef, nrow = Kt, ncol = n.coef)  
  beta.hat.ols = t(Bx.ols) %*% t(Theta)
    
  resid.mat = matrix(resid(model.ols), I, D, byrow = TRUE)
    
  ## Get Residual Structure using FPCA
  ## note: this is commented out because, in simulations based on the headstart data, 
  ## using FPCA lead to higher-than-nominal sizes for tests of nested models. 
  ## using the raw covariance worked better. using FPCA is possible, but relies
  ## on some case-specific choices.
  # raw.resid.cov = cov(resid.mat)
  # fpca.resid = fpca.sc(resid.mat, pve = .9995, nbasis = 20)
  # resid.cov = with(fpca.resid, efunctions %*% diag(evalues) %*% t(efunctions))
    
  ## account for (possibly non-constant) ME nugget effect
  # sm.diag = Theta %*% solve(crossprod(Theta)) %*% t(Theta) %*% (diag(raw.resid.cov) - diag(resid.cov))
  # if(sum( sm.diag < 0 ) >0) { sm.diag[ sm.diag < 0] = min((diag(raw.resid.cov) - diag(resid.cov))[ sm.diag < 0])}
  # diag(resid.cov) = diag(resid.cov) + sm.diag
    
  sigma = cov(resid.mat) * (I - 1) / (I - p)
  
  ## get confidence intervals
  beta.UB = beta.LB = matrix(NA, p, D)
  for(p.cur in 1:p){
    ## confidence intervals for this model shouldn't be trusted
  }
  
  Yhat = X.des %*% beta.hat.ols
  
  ret = list(beta.hat.ols, beta.UB, beta.LB, Yhat, mt_fixed, data, model.ols, sigma)
  names(ret) = c("beta.hat", "beta.UB", "beta.LB", "Yhat", "terms", "data", "model.ols", "sigma")
  class(ret) = "fosr"
  ret
    
}

###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################