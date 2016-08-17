#' Cross-sectional FoSR using Variational Bayes and Wishart prior
#' 
#' Fitting function for function-on-scalar regression for cross-sectional data.
#' This function estimates model parameters using VB and estimates
#' the residual covariance surface using a Wishart prior.
#' 
#' @param formula a formula indicating the structure of the proposed model. 
#' @param Kt number of spline basis functions used to estimate coefficient functions
#' @param data an optional data frame, list or environment containing the 
#' variables in the model. If not found in data, the variables are taken from 
#' environment(formula), typically the environment from which the function is 
#' called.
#' @param alpha tuning parameter balancing second-derivative penalty and
#' zeroth-derivative penalty (alpha = 0 is all second-derivative penalty)
#' @param min.iter minimum number of interations of VB algorithm
#' @param max.iter maximum number of interations of VB algorithm
#' @param Aw hyperparameter for inverse gamma controlling variance of spline terms
#' for population-level effects; if \code{NULL}, defaults to \code{Kt/2}.
#' @param Bw hyperparameter for inverse gamma controlling variance of spline terms
#' for population-level effects; if \code{NULL}, defaults to 
#' 1/2 tr(mu.q.beta %*% P.inv %*% mu.q.beta) where mu.q.beta is based on a OLS fit 
#' of the model
#' @param v hyperparameter for inverse Wishart prior on residual covariance; if \code{NULL},
#' Psi defaults to an FPCA decomposition of the residual covariance in which residuals are 
#' estimated based on an OLS fit of the model (note the "nugget effect" on this covariance
#' is assumed to be constant over the time domain).
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
#' @export
#' 
vb_cs_wish = function(formula, data=NULL, verbose = TRUE, Kt=5, alpha = .1, min.iter = 10, max.iter = 50,
                      Aw = NULL, Bw = NULL, v = NULL){
  
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
  
  ## fixed effect design matrix
  W.des = X
  
  I = dim(Y)[1]
  D = dim(Y)[2]
  p = dim(W.des)[2]

  ## bspline basis and penalty matrix
  Theta = bs(1:D, df = Kt, intercept=TRUE, degree=3)
  
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
    fpca.temp = fpca.sc(Y = Y - Yhat, pve = .995, var = TRUE)
    cov.hat = fpca.temp$efunctions %*% tcrossprod(diag(fpca.temp$evalues, nrow = length(fpca.temp$evalues), 
                                                       ncol = length(fpca.temp$evalues)), fpca.temp$efunctions)    
    cov.hat = cov.hat + diag(fpca.temp$sigma2, D, D)
#    cov.hat = cov(Y - Yhat)
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
  
  lpxq=c(0,1)
  j=2
  
  if(verbose) { cat("Beginning Algorithm \n") }
  
  while((j < (min.iter + 2) | (lpxq[j]-lpxq[j-1])>1.0E-1) & (j < max.iter)){
    
    ###############################################################
    ## update b-spline parameters for fixed effects
    ###############################################################
    
    sigma.q.BW = solve( Xt_siginv_X(tx = t.designmat.X, siginv = inv.sig) +
                         kronecker(diag((Aw+Kt/2)/b.q.lambda.BW), P.mat ))
    mu.q.BW = matrix( sigma.q.BW %*% Xt_siginv_X(tx = t.designmat.X, siginv = inv.sig, y = Y.vec), nrow = Kt, ncol = p)
    
    beta.cur = t(mu.q.BW) %*% t(Theta)
    
    ###############################################################
    ## update inverse covariance matrix
    ###############################################################
    
    T.BW.W = Theta %*% mu.q.BW %*% t(W.des)
    
    mu.q.v = v + I
    mu.q.Psi = Psi + t(Y) %*% Y - 
      t(Y) %*% t(T.BW.W) -
      T.BW.W %*% Y +
      T.BW.W %*% t(T.BW.W) + 
      matrix(apply(sapply(1:p, function(u) I * Theta %*% sigma.q.BW[(Kt * (u-1) +1):(Kt * u),(Kt * (u-1) +1):(Kt * u)] %*% t(Theta)), 1, sum), D, D)
    
    inv.sig = solve(mu.q.Psi/mu.q.v)
    
    ###############################################################
    ## update variance components
    ###############################################################
    
    ## lambda for fixed effects
    for(term in 1:dim(W.des)[2]){
      b.q.lambda.BW[term] = Bw[term] + .5 * (t(mu.q.BW[,term]) %*% P.mat %*% mu.q.BW[,term] + 
                                        sum(diag(P.mat %*% sigma.q.BW[(Kt*(term-1)+1):(Kt*term),(Kt*(term-1)+1):(Kt*term)])))
    }
    
    ###############################################################
    ## lower bound
    ###############################################################
    
    curlpxq = 10
    lpxq = c(lpxq, curlpxq)
    j=j+1
    
    if(verbose) { cat(".") }
    
  }
  
  lpxq=lpxq[-(1:2)]
  
  ## compute CI for fixed effects
  beta.sd = beta.LB = beta.UB = matrix(NA, nrow = p, ncol = D)
  for(i in 1:p){
    beta.sd[i,] = sqrt(diag((Theta) %*% sigma.q.BW[(Kt*(i-1)+1):(Kt*i),(Kt*(i-1)+1):(Kt*i)] %*% t(Theta)))
    beta.LB[i,] = beta.cur[i,]-1.96*beta.sd[i,]
    beta.UB[i,] = beta.cur[i,]+1.96*beta.sd[i,]
  }
  
  ## convert objects from spam to matrix
  beta.cur = as.matrix(beta.cur)
  
  ## export fitted values
  Yhat.fixed = W.des %*% beta.cur
  Yhat = Yhat.fixed
  
  ## export various r2 values
  #r2.f = 1 - (sum((Y - Yhat.fixed)^2)/(IJ*D)) / (sum((Y)^2)/(IJ*D))
  #r2.fr = 1 - (sum((Y - Yhat)^2)/(IJ*D)) / (sum((Y)^2)/(IJ*D))
  
  data = if(is.null(data)) { mf_fixed }  else { data }
  
  ret = list(beta.cur, beta.UB, beta.LB, Yhat.fixed, mt_fixed, data)
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