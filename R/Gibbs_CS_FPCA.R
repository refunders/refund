#' Cross-sectional FoSR using a Gibbs sampler and FPCA
#' 
#' Fitting function for function-on-scalar regression for cross-sectional data.
#' This function estimates model parameters using a Gibbs sampler and estimates
#' the residual covariance surface using FPCA.
#' 
#' @param formula a formula indicating the structure of the proposed model. 
#' @param Kt number of spline basis functions used to estimate coefficient functions
#' @param Kp number of FPCA basis functions to be estimated
#' @param data an optional data frame, list or environment containing the 
#' variables in the model. If not found in data, the variables are taken from 
#' environment(formula), typically the environment from which the function is 
#' called.
#' @param N.iter number of iterations used in the Gibbs sampler
#' @param N.burn number of iterations discarded as burn-in
#' @param sig2.me starting value for measurement error variance
#' @param alpha tuning parameter balancing second-derivative penalty and
#' zeroth-derivative penalty (alpha = 0 is all second-derivative penalty)
#' @param Aw hyperparameter for inverse gamma controlling variance of spline terms
#' for population-level effects
#' @param Bw hyperparameter for inverse gamma controlling variance of spline terms
#' for population-level effects
#' @param Apsi hyperparameter for inverse gamma controlling variance of spline terms
#' for FPC effects
#' @param Bpsi hyperparameter for inverse gamma controlling variance of spline terms
#' for FPC effects
#' @param seed seed value to start the sampler; ensures reproducibility
#' @param verbose logical defaulting to \code{TRUE} -- should updates on progress be printed?
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
gibbs_cs_fpca = function(formula, Kt=5, Kp=2, data=NULL, verbose = TRUE, N.iter = 5000, N.burn = 1000, SEED = NULL, 
                         sig2.me = .01, alpha = .1, Aw = NULL, Bw = NULL, Apsi = NULL, Bpsi = NULL){

  # not used now but may need this later
  call <- match.call()
  
  tf <- terms.formula(formula, specials = "re")
  trmstrings <- attr(tf, "term.labels")
  specials <- attr(tf, "specials")    # if there are no random effects this will be NULL
  where.re <-specials$re - 1
  
  # gets matrix of fixed and random effects
  if(length(where.re)!=0){
    mf_fixed <- model.frame(tf[-where.re], data = data)
    formula = tf[-where.re]
    
    # get random effects matrix
    responsename <- attr(tf, "variables")[2][[1]]
    REs = eval(parse(text=attr(tf[where.re], "term.labels")))
    
    # set up dataframe if data = NULL
    formula2 <- paste(responsename, "~", REs[[1]],sep = "")
    newfrml <- paste(responsename, "~", REs[[2]],sep = "")
    newtrmstrings <- attr(tf[-where.re], "term.labels")
    
    formula2 <- formula(paste(c(formula2, newtrmstrings), collapse = "+"))
    newfrml <- formula(paste(c(newfrml, newtrmstrings), collapse = "+"))
    mf <- model.frame(formula2, data = data)
    
    # creates the Z matrix. get rid of $zt if you want a list with more stuff.
    if(length(data)==0){Z = lme4::mkReTrms(lme4::findbars(newfrml),fr=mf)$Zt
    }else
    {Z = lme4::mkReTrms(lme4::findbars(newfrml),fr=data)$Zt}
    
    
  } else {
    mf_fixed <- model.frame(tf, data = data)
  }
  mt_fixed <- attr(mf_fixed, "terms")
  
  # get response (Y)
  Y <- model.response(mf_fixed, "numeric")
  
  if(!is.null(SEED)) { set.seed(SEED) }
  
  # x is a matrix of fixed effects
  # automatically adds in intercept
  W.des = X <- model.matrix(mt_fixed, mf_fixed, contrasts)
  
  ## subject covariates
  I = dim(X)[1]
  D = dim(Y)[2]
  p = dim(X)[2]
  
  ## bspline basis and penalty matrix
  Theta = bs(1:D, df=Kt, intercept=TRUE, degree=3)
  
  diff0 = diag(1, D, D)
  diff2 = matrix(rep(c(1,-2,1, rep(0, D-2)), D-2)[1:((D-2)*D)], D-2, D, byrow = TRUE)
  P0 = t(Theta) %*% t(diff0) %*% diff0 %*% Theta
  P2 = t(Theta) %*% t(diff2) %*% diff2 %*% Theta
  P.mat = alpha * P0 + (1-alpha) * P2
  
  ## data organization; these computations only need to be done once   
  Y.vec = as.vector(t(Y))
  t.designmat.W = t(kronecker(W.des, Theta))
  sig.W = kronecker(t(W.des) %*% W.des, t(Theta)%*% Theta)
  
  ## initial estimation and hyperparameter choice
  vec.BW = solve(kronecker(t(W.des)%*% W.des, t(Theta) %*% Theta)) %*% t(kronecker(W.des, Theta)) %*% Y.vec
  mu.q.BW = matrix(vec.BW, Kt, p)
  
  Yhat = as.matrix(W.des %*% t(mu.q.BW) %*% t(Theta))
  
  Aw = ifelse(is.null(Aw), Kt/2, Aw)
  if(is.null(Bw)){
    Bw = b.q.lambda.BW = sapply(1:p, function(u) max(1, .5*sum(diag( t(mu.q.BW[,u]) %*% P.mat %*% (mu.q.BW[,u])))))
  } else {
    Bw = b.q.lambda.BW = rep(Bw, p)
  }
  
  Apsi = ifelse(is.null(Apsi), Kt/2, Apsi)
  Bpsi = ifelse(is.null(Bpsi), Kt/2, Bpsi)
  Asig = 1; Bsig = 1
  
  ## matrices to store within-iteration estimates 
  BW = array(NA, c(Kt, p, N.iter))
      BW[,,1] = bw = matrix(0, Kt, p)
  BPSI = array(NA, c(Kt, Kp, N.iter))
      BPSI[,,1] = bpsi = matrix(0, Kt, Kp)
  C = array(NA, c(I, Kp, N.iter))
      C[,,1] = c.mat = matrix(rnorm(I*Kp, 0, .01), I, Kp)
  SIGMA = rep(NA, N.iter)
      SIGMA[1] = sig2.me = sig2.me
  LAMBDA.BW = matrix(NA, nrow = N.iter, ncol = p)
      LAMBDA.BW[1,] = lambda.bw = rep(1, p)
  LAMBDA.PSI = matrix(NA, nrow = N.iter, ncol = Kp)
      LAMBDA.PSI[1,] = lambda.psi = rep(1,Kp)
  y.post = array(NA, dim = c(I, D, (N.iter - N.burn))) 


  ## initialize estimates of fixed, random and pca effects
  beta.cur = t(bw) %*% t(Theta)
  fixef.cur = W.des %*% beta.cur
  psi.cur = t(bpsi) %*% t(Theta)
  pcaef.cur = c.mat %*% psi.cur

  if(verbose) { cat("Beginning Sampler \n") }
  
  for(i in 1:N.iter){
    
    ###############################################################
    ## update b-spline parameters for fixed effects
    ###############################################################

    mean.cur = as.vector(t(pcaef.cur))

    sigma = solve( (1/sig2.me) * sig.W + kronecker(diag(1/lambda.bw, length(lambda.bw), length(lambda.bw)), P.mat)  )
    mu = (1/sig2.me) * sigma %*% (t.designmat.W %*%  (Y.vec - mean.cur))
      
    bw = matrix(mvrnorm(1, mu = mu, Sigma = sigma), nrow = Kt, ncol = p)
    beta.cur = t(bw) %*% t(Theta)
    fixef.cur = W.des %*% beta.cur
   
    ###############################################################
    ## update b-spline parameters for PC basis functions
    ###############################################################

    mean.cur = as.vector(t(fixef.cur))

    sigma = solve( (1/sig2.me) * kronecker(t(c.mat) %*% c.mat, t(Theta)%*% Theta) + kronecker(diag(1/lambda.psi), P.mat  ))
    mu = (1/sig2.me) * sigma %*% t(kronecker(c.mat, Theta)) %*%  (Y.vec - mean.cur)

    bpsi = matrix(mvrnorm(1, mu = mu, Sigma = sigma), nrow = Kt, ncol = Kp)
    psi.cur = t(bpsi) %*% t(Theta)
    
    ppT = psi.cur %*% t(psi.cur)

    ###############################################################
    ## scores for each individual
    ###############################################################

    for(c in 1:I){
      sigma = solve( (1/sig2.me)* ppT + diag(1, Kp, Kp)  )
      mu = (1/sig2.me) * sigma %*% psi.cur %*%  (Y[c,] - fixef.cur[c,] )     
      c.mat[c,] = mvrnorm(1, mu = mu, Sigma = sigma)
    }
    
    pcaef.cur = c.mat %*% psi.cur

    ###############################################################
    ## update variance components
    ###############################################################
      
    ## sigma.me
    Y.cur = fixef.cur +  pcaef.cur
    a.post = Asig + I*D/2
    b.post = Bsig + 1/2*crossprod(as.vector(Y - Y.cur))
    sig2.me = 1/rgamma(1, a.post, b.post)

    ## lambda for beta's
    for(term in 1:p){
      a.post = Aw + Kt/2
      b.post = Bw[term] + 1/2*bw[,term] %*% P.mat %*% bw[,term]
      lambda.bw[term] = 1/rgamma(1, a.post, b.post)
    }

    ## lambda for psi's
    for(K in 1:Kp){
      a.post = Apsi + Kt/2
      b.post = Bpsi + 1/2*bpsi[,K] %*% P.mat %*% bpsi[,K]
      lambda.psi[K] = 1/rgamma(1, a.post, b.post)
    }

    ###############################################################
    ## save this iteration's parameters
    ###############################################################

    BW[,,i] = as.matrix(bw)
    BPSI[,,i] = as.matrix(bpsi)
    C[,,i] = as.matrix(c.mat)
      
    SIGMA[i] = sig2.me
    LAMBDA.BW[i,] = lambda.bw
    LAMBDA.PSI[i,] = lambda.psi

    if(i > N.burn){
      y.post[,,i - N.burn] = fixef.cur+pcaef.cur
    }
  
    if(verbose) { if(round(i %% (N.iter/10)) == 0) {cat(".")} }
    
  }
  
  ###############################################################
  ## compute posteriors for this dataset
  ###############################################################

  ## main effects
  beta.pm = beta.LB = beta.UB = matrix(NA, nrow = p, ncol = D)
  for(i in 1:p){
  	beta.post = matrix(NA, nrow = (N.iter - N.burn), ncol = D)
  	for(n in 1:(N.iter - N.burn)){
  	  beta.post[n,] = BW[,i, n + N.burn] %*% t(Theta)
  	}
  	beta.pm[i,] = apply(beta.post, 2, mean)
    beta.LB[i,] = apply(beta.post, 2, quantile, c(.025))
    beta.UB[i,] = apply(beta.post, 2, quantile, c(.975))
  }


  ## FPCA basis functions -- OUT OF DATE
  psi.pm = matrix(NA, nrow = Kp, ncol = D)
  for(i in 1:Kp){
  	psi.post = matrix(NA, nrow = (N.iter - N.burn), ncol = D)
  	for(n in 1:(N.iter - N.burn)){
  	  psi.post[n,] = BPSI[,i, n + N.burn] %*% t(Theta)
  	}
  	psi.pm[i,] = apply(psi.post, 2, mean)
  }

  ## export fitted values
  Yhat = apply(y.post, c(1,2), mean)
  y.LB = apply(y.post, c(1,2), quantile, c(.025))
  y.UB = apply(y.post, c(1,2), quantile, c(.975))
  
  ret = list(beta.pm, beta.UB, beta.LB, Yhat, mt_fixed, data, psi.pm)
  names(ret) = c("beta.hat", "beta.UB", "beta.LB", "Yhat", "terms", "data", "psi.pm")
  class(ret) = "fosr"
  ret
  
}


###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################