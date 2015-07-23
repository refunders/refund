#' Multilevel FoSR using a Gibbs sampler and FPCA
#' 
#' Fitting function for function-on-scalar regression for longitudinal data.
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
gibbs_mult_fpca = function(formula, Kt=5, Kp = 2, data=NULL, N.iter = 5000, N.burn = 1000, sig2.me = .01, alpha = .1, SEED = NULL){

  # not used now but may need this later
  call <- match.call()
  
  tf <- terms.formula(formula, specials = "re.fosr")
  trmstrings <- attr(tf, "term.labels")
  specials <- attr(tf, "specials")    # if there are no random effects this will be NULL
  where.re.fosr <-specials$re.fosr - 1
  
  # gets matrix of fixed and random effects
  if(length(where.re.fosr)!=0){
    mf_fixed <- model.frame(tf[-where.re.fosr], data = data)
    formula = tf[-where.re.fosr]
    
    # get random effects matrix
    responsename <- attr(tf, "variables")[2][[1]]
    REs = eval(parse(text=attr(tf[where.re.fosr], "term.labels")))
    
    # set up dataframe if data = NULL
    formula2 <- paste(responsename, "~", REs[[1]],sep = "")
    newfrml <- paste(responsename, "~", REs[[2]],sep = "")
    newtrmstrings <- attr(tf[-where.re.fosr], "term.labels")
    
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
  
  # x is a matrix of fixed effects
  # automatically adds in intercept
  X <- model.matrix(mt_fixed, mf_fixed, contrasts)
  
  if(!is.null(SEED)) { set.seed(SEED) }
    
  ## fixed and random effect design matrices
  W.des = X
  Z.des = t(as.matrix(Z))

  I = dim(Z.des)[2]
  D = dim(Y)[2]
  Ji = as.numeric(apply(Z.des, 2, sum))
  IJ = sum(Ji)

  SUBJ = factor(apply(Z.des %*% 1:dim(Z.des)[2], 1, sum))

  ## find first observation
  firstobs = rep(NA, length(unique(SUBJ)))
  for(i in 1:length(unique(SUBJ))){
    firstobs[i] = which(SUBJ %in% unique(SUBJ)[i])[1]
  }
  Wi = W.des[firstobs,]
  
  ## bspline basis and penalty matrix
  Theta = bs(1:D, df = Kt, intercept=TRUE, degree=3)

  diff0 = diag(1, D, D)
  diff2 = matrix(rep(c(1,-2,1, rep(0, D-2)), D-2)[1:((D-2)*D)], D-2, D, byrow = TRUE)
  P0 = t(Theta) %*% t(diff0) %*% diff0 %*% Theta
  P2 = t(Theta) %*% t(diff2) %*% diff2 %*% Theta
  P.mat = alpha * P0 + (1-alpha) * P2
  
  ## hyper parameters for inverse gaussian
  A = .01
  B = .01

  ## number of fixed effects
  p = dim(W.des)[2]


  ## matrices to store within-iteration estimates 
  BW = array(NA, c(Kt, p, N.iter))
      BW[,,1] = matrix(0, Kt, p)
      bw = BW[,,1]
  BZ = array(NA, c(Kt, I, N.iter))
      BZ[,,1] = matrix(0, Kt, I)
      bz = BZ[,,1]
  Bpsi = array(NA, c(Kt, Kp, N.iter))
      Bpsi[,,1] = matrix(0, Kt, Kp)
      bpsi = Bpsi[,,1]
  C = array(NA, c(IJ, Kp, N.iter))
      C[,,1] = matrix(rnorm(IJ*Kp, 0, .01), IJ, Kp)
      c.mat = C[,,1]
  SIGMA = rep(NA, N.iter)
      SIGMA[1] = sig2.me
      sig2.me = SIGMA[1]
  LAMBDA.BW = matrix(NA, nrow = N.iter, ncol = p)
      LAMBDA.BW[1,] = rep(1, p)
      lambda.bw = LAMBDA.BW[1,]
  LAMBDA.BZ = rep(NA, N.iter)
      LAMBDA.BZ[1] = 1
      lambda.ranef = LAMBDA.BZ[1]
  LAMBDA.PSI = matrix(NA, nrow = N.iter, ncol = Kp)
      LAMBDA.PSI[1,] = rep(1,Kp)
      lambda.psi = LAMBDA.PSI[1,]

  ## data organization; these computations only need to be done once   
  Y.vec = as.vector(t(Y))
  t.designmat.W = t(kronecker(W.des, Theta))
  sig.W = kronecker(t(W.des) %*% W.des, t(Theta)%*% Theta)

  ## initialize estimates of fixed, random and pca effects
  beta.cur = t(bw) %*% t(Theta)
  fixef.cur = W.des %*% beta.cur
  ranef.cur = Z.des %*% t(bz) %*% t(Theta)
  psi.cur = t(bpsi) %*% t(Theta)
  pcaef.cur = c.mat %*% psi.cur

  if(verbose) { cat("Beginning Sampler \n") }
  
  for(i in 1:N.iter){
    
    ###############################################################
    ## update b-spline parameters for fixed effects
    ###############################################################

    mean.cur = as.vector(t(ranef.cur + pcaef.cur))

    sigma = solve( (1/sig2.me) * sig.W + kronecker(diag(1/lambda.bw), P.mat)  )
    mu = (1/sig2.me) * sigma %*% (t.designmat.W %*%  (Y.vec - mean.cur))
      
    bw = matrix(mvrnorm(1, mu = mu, Sigma = sigma), nrow = Kt, ncol = p)
    beta.cur = t(bw) %*% t(Theta)
    fixef.cur = W.des %*% beta.cur

    ###############################################################
    ## update b-spline parameters for subject random effects
    ###############################################################
      
    for(subj in 1:length(unique(SUBJ))){

      t.designmat.Z = t(kronecker(Theta, rep(1, Ji[subj])))
      sig.Z = kronecker(Ji[subj], t(Theta)%*% Theta)

      mean.cur = as.vector((fixef.cur[which(SUBJ == unique(SUBJ)[subj]), ] + 
                            pcaef.cur[which(SUBJ == unique(SUBJ)[subj]), ])) 

      sigma = solve( (1/sig2.me) * sig.Z + (1/lambda.ranef) * P.mat )
      mu = (1/sig2.me) * sigma %*% (t.designmat.Z %*%  (as.vector(Y[which(SUBJ == unique(SUBJ)[subj]),]) - mean.cur))

      bz[,subj] = mvrnorm(1, mu = mu, Sigma = sigma)

    }
    
    ranef.cur = Z.des %*% t(bz) %*% t(Theta)
    
    ###############################################################
    ## update b-spline parameters for PC basis functions
    ###############################################################

    mean.cur = as.vector(t(fixef.cur + ranef.cur))

    sigma = solve( (1/sig2.me) * kronecker(t(c.mat) %*% c.mat, t(Theta)%*% Theta) + kronecker(diag(1/lambda.psi), P.mat  ))
    mu = (1/sig2.me) * sigma %*% t(kronecker(c.mat, Theta)) %*%  (Y.vec - mean.cur)

    bpsi = matrix(mvrnorm(1, mu = mu, Sigma = sigma), nrow = Kt, ncol = Kp)
    psi.cur = t(bpsi) %*% t(Theta)
    ppT = psi.cur %*% t(psi.cur)

    ###############################################################
    ## scores for each individual
    ###############################################################

    for(subj.vis in 1:(IJ)){
      sigma = solve( (1/sig2.me)* ppT + diag(1, Kp, Kp)  )
      mu = (1/sig2.me) * sigma %*% psi.cur %*%  (Y[subj.vis,] - fixef.cur[subj.vis,] - ranef.cur[subj.vis,]  )
      
      c.mat[subj.vis,] = mvrnorm(1, mu = mu, Sigma = sigma)
    }
    
    pcaef.cur = c.mat %*% psi.cur

    ###############################################################
    ## update variance components
    ###############################################################
      
    ## sigma.me
    Y.cur = fixef.cur + ranef.cur +  pcaef.cur
    a.post = A + IJ*D/2
    b.post = B + 1/2*crossprod(as.vector(Y - Y.cur))

    sig2.me = 1/rgamma(1, a.post, b.post)

    ## lambda for beta's
    for(term in 1:p){
      a.post = A + Kt/2
      b.post = B + 1/2*bw[,term] %*% P.mat %*% bw[,term]
      lambda.bw[term] = 1/rgamma(1, a.post, b.post)
    }
      
    ## lambda for random effects
    a.post = A + I*Kt/2
    b.post = B + 1/2*sum(diag(t(bz) %*% P.mat %*% bz))
    lambda.ranef = 1/rgamma(1, a.post, b.post)
            
    ## lambda for psi's
    for(K in 1:Kp){
      a.post = A + Kt/2
      b.post = B + 1/2*bpsi[,K] %*% P.mat %*% bpsi[,K]
      lambda.psi[K] = 1/rgamma(1, a.post, b.post)
    }

    ###############################################################
    ## save this iteration's parameters
    ###############################################################

    BW[,,i] = as.matrix(bw)
    BZ[,,i] = as.matrix(bz)
    Bpsi[,,i] = as.matrix(bpsi)
    C[,,i] = as.matrix(c.mat)
      
    SIGMA[i] = sig2.me
    LAMBDA.BW[i,] = lambda.bw
    LAMBDA.BZ[i] = lambda.ranef
    LAMBDA.PSI[i,] = lambda.psi

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

  ## random effects
  b.pm = matrix(NA, nrow = I, ncol = D)
  for(i in 1:I){
  	b.post = matrix(NA, nrow = (N.iter - N.burn), ncol = D)
  	for(n in 1:(N.iter - N.burn)){
  	  b.post[n,] = BZ[,i, n + N.burn] %*% t(Theta)
  	}
  	b.pm[i,] = apply(b.post, 2, mean)
  }

  ## FPCA basis functions -- OUT OF DATE
  psi.pm = matrix(NA, nrow = Kp, ncol = D)
  for(i in 1:Kp){
  	psi.post = matrix(NA, nrow = (N.iter - N.burn), ncol = D)
  	for(n in 1:(N.iter - N.burn)){
  	  psi.post[n,] = Bpsi[,i, n + N.burn] %*% t(Theta)
  	}
  	psi.pm[i,] = apply(psi.post, 2, mean)
  }

  ## export fitted values
  fixef.pm = W.des %*% beta.pm
  ranef.pm = Z.des %*% b.pm
  Yhat = matrix(NA, nrow = IJ, ncol = D)
#  for(i in 1:IJ){
#  	y.post = matrix(NA, nrow = (N.iter - N.burn), ncol = D)
#  	for(n in 1:(N.iter - N.burn)){
#  	  y.post[n,] = W.des[i,] %*% (t(BW[,,n + N.burn]) %*% t(Theta)) +
#  	               Z.des[i,] %*% (t(BZ[,,n + N.burn]) %*% t(Theta)) +
#  	               C[i,, n + N.burn] %*% (t(Bpsi[,,n + N.burn]) %*% t(Theta))
#  	}
#  	Yhat[i,] = apply(y.post, 2, mean)
#  }
  
  ## rotate scores to be orthonormal
  psi.pm = t(svd(t(psi.pm))$u)
  
  ret = list(beta.pm, beta.UB, beta.LB, fixef.pm, mt_fixed, data, psi.pm)
  names(ret) = c("beta.hat", "beta.UB", "beta.LB", "Yhat", "terms", "data", "psi.pm")
  class(ret) = "fosr"

}


###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################