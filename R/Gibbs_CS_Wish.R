#' Cross-sectional FoSR using a Gibbs sampler and Wishart prior
#' 
#' Fitting function for function-on-scalar regression for cross-sectional data.
#' This function estimates model parameters using a Gibbs sampler and estimates
#' the residual covariance surface using a Wishart prior.
#' 
#' @param formula a formula indicating the structure of the proposed model. 
#' @param Kt number of spline basis functions used to estimate coefficient functions
#' @param data an optional data frame, list or environment containing the 
#' variables in the model. If not found in data, the variables are taken from 
#' environment(formula), typically the environment from which the function is 
#' called.
#' @param N.iter number of iterations used in the Gibbs sampler
#' @param N.burn number of iterations discarded as burn-in
#' @param alpha tuning parameter balancing second-derivative penalty and
#' zeroth-derivative penalty (alpha = 0 is all second-derivative penalty)
#' @param Aw hyperparameter for inverse gamma controlling variance of spline terms
#' for population-level effects
#' @param Bw hyperparameter for inverse gamma controlling variance of spline terms
#' for population-level effects
#' @param v hyperparameter for inverse Wishart prior on residual covariance
#' @param SEED seed value to start the sampler; ensures reproducibility
#' @param verbose logical defaulting to \code{TRUE} -- should updates on progress be printed?
#' @param min.iter minimum number of iterations
#' @param max.iter maximum number of iterations
#' 
#' @references
#' Goldsmith, J., Kitago, T. (Under Review).
#' Assessing Systematic Effects of Stroke on Motor Control using Hierarchical 
#' Function-on-Scalar Regression.
#' 
#' @author Jeff Goldsmith \email{ajg2202@@cumc.columbia.edu}
#' @importFrom splines bs
#' @importFrom MASS mvrnorm
#' @export
#' 
gibbs_cs_wish = function(formula, Kt=5, data=NULL, verbose = TRUE, N.iter = 5000, N.burn = 1000, alpha = .1, 
                         min.iter = 10, max.iter = 50, Aw = NULL, Bw = NULL, v = NULL, SEED = NULL){

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
  
  if(!is.null(SEED)) { set.seed(SEED) }
  
  ## fixed effect design matrix
  W.des = X

  I = dim(Y)[1]
  p = dim(W.des)[2]
  D = dim(Y)[2]
  
  ## bspline basis and penalty matrix
  Theta = bs(1:D, df=Kt, intercept=TRUE, degree=3)

  diff0 = diag(1, D, D)
  diff2 = matrix(rep(c(1,-2,1, rep(0, D-2)), D-2)[1:((D-2)*D)], D-2, D, byrow = TRUE)
  P0 = t(Theta) %*% t(diff0) %*% diff0 %*% Theta
  P2 = t(Theta) %*% t(diff2) %*% diff2 %*% Theta
  P.mat = alpha * P0 + (1-alpha) * P2

  ## data organization; these computations only need to be done once
  Y.vec = as.vector(t(Y))
  t.designmat.X = t(kronecker(W.des, Theta))
  sig.X = kronecker(t(W.des) %*% W.des, t(Theta)%*% Theta)

  ## initial estimation and hyperparameter choice
  vec.BW = solve(kronecker(t(W.des)%*% W.des, t(Theta) %*% Theta)) %*% t(kronecker(W.des, Theta)) %*% Y.vec
  mu.q.BW = matrix(vec.BW, Kt, p)
  
  Yhat = as.matrix(W.des %*% t(mu.q.BW) %*% t(Theta))
  
  if(is.null(v)){
    fpca.temp = fpca.sc(Y = Y - Yhat, pve = .95, var = TRUE)
    cov.hat = fpca.temp$efunctions %*% tcrossprod(diag(fpca.temp$evalues, nrow = length(fpca.temp$evalues), 
                                                       ncol = length(fpca.temp$evalues)), fpca.temp$efunctions)    
    cov.hat = cov.hat + diag(fpca.temp$sigma2, D, D)
    Psi = cov.hat * I
  } else {
    Psi = diag(v, D, D)
  }
  
  v = ifelse(is.null(v), I, v)
  inv.sig = solve(Psi/v)
  
  Aw = ifelse(is.null(Aw), Kt/2, Aw)
  if(is.null(Bw)){
    Bw = b.q.lambda.BW = sapply(1:p, function(u) max(1, .5*sum(diag( t(mu.q.BW[,u]) %*% P.mat %*% (mu.q.BW[,u])))))
  } else {
    Bw = b.q.lambda.BW = rep(Bw, p)
  }
  
  ## matrices to store within-iteration estimates 
  BW = array(NA, c(Kt, p, N.iter))
    BW[,,1] = bw = matrix(rnorm(Kt * p, 0, 10), Kt, p)
  INV.SIG = array(NA, c(D, D, N.iter))
    INV.SIG[,,1] = inv.sig = diag(10, D, D)
  LAMBDA.BW = matrix(NA, nrow = N.iter, ncol = p)
    LAMBDA.BW[1,] = lambda.bw = runif(p, .1, 10)
  
  y.post = array(NA, dim = c(I, D, (N.iter - N.burn)))

  if(verbose) { cat("Beginning Sampler \n") }

  for(i in 1:N.iter){

    ###############################################################
    ## update b-spline parameters for fixed effects
    ###############################################################
    
    sigma = solve(Xt_siginv_X(tx = t.designmat.X, siginv = inv.sig) + 
                  kronecker(diag(lambda.bw), P.mat ))
    mu = sigma %*% Xt_siginv_X(tx = t.designmat.X, siginv = inv.sig, y = Y.vec)
      
    bw = matrix(mvrnorm(1, mu = mu, Sigma = sigma), nrow = Kt, ncol = p)

    beta.cur = t(bw) %*% t(Theta)

    ###############################################################
    ## update inverse covariance matrix
    ###############################################################

    resid.cur = Y - W.des %*% beta.cur
    inv.sig = rWishart(1, v + I, solve(Psi + t(resid.cur) %*% resid.cur))
    
    ###############################################################
    ## update variance components
    ###############################################################

    ## lambda for beta's
    for(term in 1:p){
      a.post = Aw + Kt/2
      b.post = Bw[term] + 1/2 * bw[,term] %*% P.mat %*% bw[,term]
      lambda.bw[term] = rgamma(1, a.post, b.post)
    }
      
    ###############################################################
    ## save this iteration's parameters
    ###############################################################

    BW[,,i] = as.matrix(bw)
    
    INV.SIG[,,i] = inv.sig
    LAMBDA.BW[i,] = lambda.bw
    
    if(verbose) { if(round(i %% (N.iter/10)) == 0) {cat(".")} }
    
  }

  ###############################################################
  ## compute posteriors for this dataset
  ###############################################################

  ## main effects
  beta.post = array(NA, dim = c(p, D, (N.iter - N.burn)))
  for(n in 1:(N.iter - N.burn)){
    beta.post[,,n] = t(BW[,, n + N.burn]) %*% t(Theta)
  }
  beta.pm = apply(beta.post, c(1,2), mean)
  beta.LB = apply(beta.post, c(1,2), quantile, c(.025))
  beta.UB = apply(beta.post, c(1,2), quantile, c(.975))


  ## covariance matrix
  sig.pm = solve(apply(INV.SIG, c(1,2), mean))
  
  ## export fitted values
  fixef.pm = W.des %*% beta.pm
  
  ret = list(beta.pm, beta.UB, beta.LB, fixef.pm, mt_fixed, data)
  names(ret) = c("beta.hat", "beta.UB", "beta.LB", "Yhat", "terms", "data")
  class(ret) = "fosr"
  ret
  
}


###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################