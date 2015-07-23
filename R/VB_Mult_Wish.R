#' Multilevel FoSR using Variational Bayes and Wishart prior
#' 
#' Fitting function for function-on-scalar regression for cross-sectional data.
#' This function estimates model parameters using VB and estimates
#' the residual covariance surface using a Wishart prior. If prior hyperparameters
#' are \code{NULL} they are estimated using the data.
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
#' @param Az hyperparameter for inverse gamma controlling variance of spline terms
#' for subject-level effects
#' @param Bz hyperparameter for inverse gamma controlling variance of spline terms
#' for subject-level effects
#' @param Aw hyperparameter for inverse gamma controlling variance of spline terms
#' for population-level effects
#' @param Bw hyperparameter for inverse gamma controlling variance of spline terms
#' for population-level effects
#' @param v hyperparameter for inverse Wishart prior on residual covariance
#' @param verbose logical defaulting to \code{TRUE} -- should updates on progress be printed?
#' 
#' @references
#' Goldsmith, J., Kitago, T. (Under Review).
#' Assessing Systematic Effects of Stroke on Motor Control using Hierarchical 
#' Function-on-Scalar Regression.
#' 
#' @author Jeff Goldsmith \email{ajg2202@@cumc.columbia.edu}
#' @importFrom splines bs
#' @export
#' 
vb_mult_wish = function(formula, data=NULL, verbose = TRUE, Kt = 5, alpha = .1, min.iter = 10, max.iter = 50,
                        Az = NULL, Bz = NULL, Aw = NULL, Bw = NULL, v = NULL){

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
  
  ## fixed and random effect design matrices
  W.des = X
  Z.des = t(as.matrix(Z))
  
  I = dim(Z.des)[2]
  D = dim(Y)[2]
  Ji = as.numeric(apply(Z.des, 2, sum))
  IJ = sum(Ji)
  p = dim(W.des)[2]

  ## bspline basis and penalty matrix
  Theta = bs(1:D, df=Kt, intercept=TRUE, degree=3)

  diff0 = diag(1, D, D)
  diff2 = matrix(rep(c(1,-2,1, rep(0, D-2)), D-2)[1:((D-2)*D)], D-2, D, byrow = TRUE)
  P0 = t(Theta) %*% t(diff0) %*% diff0 %*% Theta
  P2 = t(Theta) %*% t(diff2) %*% diff2 %*% Theta
  P.mat = alpha * P0 + (1-alpha) * P2

  SUBJ = factor(apply(Z.des %*% 1:dim(Z.des)[2], 1, sum))

  ## find first observation
  firstobs = rep(NA, length(unique(SUBJ)))
  for(i in 1:length(unique(SUBJ))){
    firstobs[i] = which(SUBJ %in% unique(SUBJ)[i])[1]
  }
  Wi = W.des[firstobs,]

  ## data organization; these computations only need to be done once
  Y.vec = as.vector(t(Y))
  IIP = kronecker(diag(1, I, I), P.mat)
  WIk = kronecker(Wi, diag(1, Kt, Kt))
  tWIW = t(WIk) %*% IIP %*% WIk
  tWI = t(WIk) %*% IIP

  ## initial estimation and hyperparameter choice
  mu.q.BZ = matrix(NA, nrow = Kt, ncol = I)
  for(subj in 1:length(unique(SUBJ))){
    mu.q.BZ[,subj] = solve(kronecker(Ji[subj], t(Theta) %*% Theta)) %*% t(kronecker(rep(1, Ji[subj]), Theta)) %*% 
                     (as.vector(t(Y[which(SUBJ == unique(SUBJ)[subj]),])))
    
  }
  vec.BZ = as.vector(mu.q.BZ)
  
  vec.BW = solve(tWIW) %*% tWI %*% vec.BZ
  mu.q.BW = matrix(vec.BW, Kt, p)

  Yhat = as.matrix(Z.des %*% t(mu.q.BZ) %*% t(Theta))

  if(is.null(v)){
    fpca.temp = fpca.sc(Y = Y - Yhat, pve = .95, var = TRUE)
    cov.hat = fpca.temp$efunctions %*% tcrossprod(diag(fpca.temp$evalues, nrow = length(fpca.temp$evalues), 
                                                       ncol = length(fpca.temp$evalues)), fpca.temp$efunctions)    
    cov.hat = cov.hat + diag(fpca.temp$sigma2, D, D)
    Psi = cov.hat * IJ
  } else {
    Psi = diag(v, D, D)
  }
  
  v = ifelse(is.null(v), IJ, v)
  inv.sig = solve(Psi/v)
  
  Az = ifelse(is.null(Az), I*Kt/2, Az)
  Bz = b.q.lambda.BZ = ifelse(is.null(Bz), .5*sum(diag((t(mu.q.BZ) - Wi %*% t(mu.q.BW)) %*% P.mat %*% t(t(mu.q.BZ) - Wi %*% t(mu.q.BW)))), Bz)
  
  Aw = ifelse(is.null(Aw), Kt/2, Aw)
  if(is.null(Bw)){
    Bw = b.q.lambda.BW = sapply(1:p, function(u) max(1, .5*sum(diag( t(mu.q.BW[,u]) %*% P.mat %*% (mu.q.BW[,u])))))
  } else {
    Bw = b.q.lambda.BW = rep(Bw, p)
  }

  sigma.q.BZ = vector("list", I)

  lpxq=c(0,1)
  j=2

  if(verbose) { cat("Beginning Algorithm \n") }

  while((j < (min.iter + 2) | (lpxq[j]-lpxq[j-1])>1.0E-1) & (j < max.iter)){
  
    ###############################################################
    ## update b-spline parameters for subject random effects
    ###############################################################

    for(subj in 1:length(unique(SUBJ))){
      
      t.designmat.Z = t(kronecker(rep(1, Ji[subj]), Theta))

      sigma.q.BZ[[subj]] = solve(t.designmat.Z %*% kronecker(diag(1, Ji[subj], Ji[subj]), inv.sig) %*% t(t.designmat.Z) + 
                         ((Az+I*Kt/2)/b.q.lambda.BZ) * P.mat)
      mu.q.BZ[,subj] = sigma.q.BZ[[subj]] %*% (t.designmat.Z %*% kronecker(diag(1, Ji[subj], Ji[subj]), inv.sig) %*% 
                                               (as.vector(t(Y[which(SUBJ == unique(SUBJ)[subj]),]))) + 
                          (((Az+I*Kt/2)/b.q.lambda.BZ) * P.mat) %*% mu.q.BW %*% (Wi[subj,]) )
    
    }

    ranef.cur = as.matrix(Z.des %*% t(mu.q.BZ) %*% t(Theta))

    ###############################################################
    ## update b-spline parameters for fixed effects
    ###############################################################
    
    sigma.q.BW = solve(  (((Az+I*Kt/2)/b.q.lambda.BZ) * tWIW) + kronecker(diag((Aw+Kt/2)/b.q.lambda.BW), P.mat ))
    mu.q.BW = matrix( sigma.q.BW %*% (((Az+I*Kt/2)/b.q.lambda.BZ) * tWI %*% as.vector(mu.q.BZ)), nrow = Kt, ncol = p)

    beta.cur = t(mu.q.BW) %*% t(Theta)
 
    ###############################################################
    ## update inverse covariance matrix
    ###############################################################
    
    T.BZ.Z = Theta %*% mu.q.BZ %*% t(Z.des)

    mu.q.v = v + IJ
    mu.q.Psi = Psi + t(Y) %*% Y - 
               t(Y) %*% t(T.BZ.Z) -
               T.BZ.Z %*% Y +
               T.BZ.Z %*% t(T.BZ.Z) + 
               matrix(apply(sapply(1:I, function(u) Ji[u] * Theta %*% sigma.q.BZ[[u]] %*% t(Theta)), 1, sum), D, D)

    inv.sig = solve(mu.q.Psi/mu.q.v)

    ###############################################################
    ## update variance components
    ###############################################################

    ## lambda for fixed effects
    for(term in 1:dim(W.des)[2]){
      b.q.lambda.BW[term] = Bw[term] + .5 * (t(mu.q.BW[,term]) %*% P.mat %*% mu.q.BW[,term] + 
                                        sum(diag(P.mat %*% sigma.q.BW[(Kt*(term-1)+1):(Kt*term),(Kt*(term-1)+1):(Kt*term)])))
    }

    ## lambda for random effects
    trace.temp = sapply(1:dim(W.des)[2], function(p) sum(diag(P.mat %*% sigma.q.BW[(Kt*(p-1)+1):(Kt*p),(Kt*(p-1)+1):(Kt*p)])) )
    b.q.lambda.BZ = as.numeric(Bz + .5 * sum(sapply(1:I, function(u) {
    	                                  mu.q.BZ[,u] %*% P.mat %*% mu.q.BZ[,u] + sum(diag(P.mat %*% sigma.q.BZ[[u]])) -
    	                                  2 * Wi[u,] %*% t(mu.q.BW) %*% P.mat %*% mu.q.BZ[,u] +
    	                                  Wi[u,] %*% t(mu.q.BW) %*% P.mat %*% mu.q.BW %*% (Wi[u,]) + 
    	                                  Wi[u,] %*% diag(trace.temp) %*% (Wi[u,])
                                        })))
    
    ###############################################################
    ## lower bound
    ###############################################################

    curlpxq = .5 * sum(sapply(1:I, function(u) det(sigma.q.BZ[[u]]))) + 
              .5 * sum(sapply(1:p, function(u) det(sigma.q.BW[(Kt*(p-1)+1):(Kt*p),(Kt*(p-1)+1):(Kt*p)]))) -
              (Az + (I*Kt)/2) * log(b.q.lambda.BZ) - 
              (Aw + Kt/2) * sum(b.q.lambda.BW) -
              (v + IJ)/2 * log(det(mu.q.Psi))
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
  Yhat = ranef.cur
  
  ## export various r2 values
  r2.f = 1 - (sum((Y - Yhat.fixed)^2)/(IJ*D)) / (sum((Y)^2)/(IJ*D))
  r2.fr = 1 - (sum((Y - Yhat)^2)/(IJ*D)) / (sum((Y)^2)/(IJ*D))
    
  ret = list(beta.cur, beta.UB, beta.LB, Yhat)
  names(ret) = c("beta.pm", "beta.UB", "beta.LB", "Yhat")
  
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