#' Multilevel FoSR using Variational Bayes and FPCA
#' 
#' Fitting function for function-on-scalar regression for multilevel data.
#' This function estimates model parameters using a VB and estimates
#' the residual covariance surface using FPCA.
#' 
#' @param formula a formula indicating the structure of the proposed model. 
#' @param Kt number of spline basis functions used to estimate coefficient functions
#' @param Kp number of FPCA basis functions to be estimated
#' @param data an optional data frame, list or environment containing the 
#' variables in the model. If not found in data, the variables are taken from 
#' environment(formula), typically the environment from which the function is 
#' called.
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
vb_mult_fpca = function(formula, data=NULL, verbose = TRUE, Kt=5, Kp=2, alpha = .1){

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
  
  # x is a matrix of fixed effects
  # automatically adds in intercept
  X <- model.matrix(mt_fixed, mf_fixed, contrasts)

  ## fixed and random effect design matrices
  W.des = X
  Z.des = t(as.matrix(Z))

  ## subject covariates
  I = dim(Z.des)[2]
  D = dim(Y)[2]
  Ji = as.numeric(apply(Z.des, 2, sum))
  IJ = sum(Ji)
  p = dim(W.des)[2]

  SUBJ = factor(apply(Z.des %*% 1:dim(Z.des)[2], 1, sum))

  ## find first observation
  firstobs = rep(NA, length(unique(SUBJ)))
  for(i in 1:length(unique(SUBJ))){
    firstobs[i] = which(SUBJ %in% unique(SUBJ)[i])[1]
  }
  Wi = W.des[firstobs,]

  ## bspline basis and penalty matrix
  Theta = bs(1:D, df=Kt, intercept=TRUE, degree=3)

  diff0 = diag(1, D, D)
  diff2 = matrix(rep(c(1,-2,1, rep(0, D-2)), D-2)[1:((D-2)*D)], D-2, D, byrow = TRUE)
  P0 = t(Theta) %*% t(diff0) %*% diff0 %*% Theta
  P2 = t(Theta) %*% t(diff2) %*% diff2 %*% Theta
  P.mat = alpha * P0 + (1-alpha) * P2

  ## hyper parameters for inverse gaussians
  A = .01
  B = .01

  ## matrices to to approximate paramater values
  sigma.q.BW = vector("list", p)
  for(k in 1:p){
    sigma.q.BW[[k]] = diag(1, Kt)
  }
  mu.q.BW = matrix(0, nrow = Kt, ncol = p)  

  sigma.q.BZ = diag(1, Kt)
  mu.q.BZ = matrix(0, nrow = Kt, ncol = length(unique(SUBJ)))  

  sigma.q.Bpsi = vector("list", Kp)
  for(k in 1:Kp){
    sigma.q.Bpsi[[k]] = diag(1, Kt)
  }
  mu.q.Bpsi = matrix(0, nrow = Kt, ncol = Kp)  

  sigma.q.C = vector("list", IJ)
  for(k in 1:IJ){
    sigma.q.C[[k]] = diag(1, Kp)
  }
  mu.q.C = matrix(rnorm(IJ*Kp, 0, .01), IJ, Kp)
    
  b.q.lambda.BW = rep(1, p)
  b.q.lambda.BZ = 1
  b.q.lambda.Bpsi = rep(1, Kp)
  b.q.sigma.me = 1

  ## data organization; these computations only need to be done once
  Y.vec = as.vector(t(Y))
  t.designmat.X = t(kronecker(W.des, Theta))
  sig.X = kronecker(t(W.des) %*% W.des, t(Theta)%*% Theta)

  ## initialize estimates of fixed, random and pca effects
  fixef.cur = matrix(0, nrow = IJ, ncol = D)
  ranef.cur = matrix(0, nrow = IJ, ncol = D)
  pcaef.cur = matrix(0, IJ, D)
  
  lpxq=c(0,1)
  j=2

  if(verbose) { cat("Beginning Algorithm \n") }

#  while(j<4 | (lpxq[j]-lpxq[j-1])>1.0E-1){
  while(j<11){
  
    ###############################################################
    ## update b-spline parameters for fixed effects
    ###############################################################

    mean.cur = as.vector(t(ranef.cur + pcaef.cur))
    
    sigma.q.BW = solve(as.numeric((A + IJ*D/2)/(b.q.sigma.me)) * sig.X + kronecker(diag((A+Kt/2)/b.q.lambda.BW), P.mat ))
    mu.q.BW = matrix(sigma.q.BW %*% (as.numeric((A + IJ*D/2)/(b.q.sigma.me)) * t.designmat.X %*%  (Y.vec - mean.cur)), nrow = Kt, ncol = p)

    beta.cur = t(mu.q.BW) %*% t(Theta)
    fixef.cur = as.matrix(W.des %*% beta.cur)

    ###############################################################
    ## update b-spline parameters for subject random effects
    ###############################################################

    for(subj in 1:length(unique(SUBJ))){
      
      t.designmat.Z = t(kronecker(Theta, rep(1, Ji[subj])))
      sig.Z = kronecker(t(rep(1, Ji[subj])) %*% rep(1, Ji[subj]), t(Theta)%*% Theta)

      mean.cur = as.vector(fixef.cur[which(SUBJ == unique(SUBJ)[subj]),] + pcaef.cur[which(SUBJ == unique(SUBJ)[subj]),])

      sigma.q.BZ = solve( as.numeric((A + IJ*D/2)/(b.q.sigma.me)) * sig.Z + ((A+I*Kt/2)/b.q.lambda.BZ) * P.mat  )
      mu.q.BZ[,subj] = sigma.q.BZ %*% (as.numeric((A + IJ*D/2)/(b.q.sigma.me)) * t.designmat.Z %*%  (as.vector(Y[which(SUBJ == unique(SUBJ)[subj]),]) - mean.cur))
    
    }

    ranef.cur = as.matrix(Z.des %*% t(mu.q.BZ) %*% t(Theta))

    ###############################################################
    ## update b-spline parameters for PC basis functions
    ###############################################################

    mean.cur = as.vector(t(fixef.cur + ranef.cur))
    designmat = kronecker(mu.q.C, Theta)
    
    sigma.q.Bpsi = solve( as.numeric((A + IJ*D/2)/(b.q.sigma.me)) * kronecker(t(mu.q.C)%*%mu.q.C + diag(IJ, Kp, Kp), t(Theta)%*%Theta) + 
                                                                   kronecker(diag((A+Kt/2)/b.q.lambda.Bpsi), P.mat  ))
    mu.q.Bpsi = matrix(((A + IJ*D/2)/(b.q.sigma.me)) * sigma.q.Bpsi %*% (t(designmat) %*%  (Y.vec - mean.cur)), nrow = Kt, ncol = Kp)

    psi.cur = t(mu.q.Bpsi) %*% t(Theta)
    ppT = (psi.cur) %*% t(psi.cur)

    ###############################################################
    ## scores for each individual
    ###############################################################

    for(subj in 1:IJ){

      sigma.q.C[[subj]] = solve( ((A + IJ*D/2)/(b.q.sigma.me))* ppT + 
                                     diag(sapply(1:Kp, function(u) sum(diag(sigma.q.Bpsi[(Kt*(u-1)+1):(Kt*u),(Kt*(u-1)+1):(Kt*u)] %*% t(Theta) %*% Theta  )))) + ##### note -- double check this
                                     diag(1, Kp, Kp)  )

      mu.q.C[subj,] = ((A + IJ*D/2)/(b.q.sigma.me)) * sigma.q.C[[subj]] %*% as.matrix(psi.cur) %*%  (Y[subj,] - fixef.cur[subj,] - ranef.cur[subj,]  )
    }

  	pcaef.cur =  as.matrix(mu.q.C %*% psi.cur)
 
    ###############################################################
    ## update variance components
    ###############################################################

    ## measurement error variance
    b.q.sigma.me = as.numeric(B + .5 *crossprod(as.vector(Y - (fixef.cur + ranef.cur + pcaef.cur))))

    ## lambda for fixed effects
    for(term in 1:dim(W.des)[2]){
      b.q.lambda.BW[term] = B + .5 * (t(mu.q.BW[,term]) %*% P.mat %*% mu.q.BW[,term] + 
                                        sum(diag(P.mat %*% sigma.q.BW[(Kt*(term-1)+1):(Kt*term),(Kt*(term-1)+1):(Kt*term)])))
    }

    ## lambda for random effects
    vec.BZ = as.vector(mu.q.BZ)
    b.q.lambda.BZ = as.numeric(B + .5 * (t(vec.BZ) %*% kronecker(diag(1, length(Ji)), P.mat) %*% vec.BZ +
                   I * sum(diag(P.mat %*% sigma.q.BZ))))

    ## lambda for FPCA basis functions
    for(K in 1:Kp){
      b.q.lambda.Bpsi[K] = B + .5 * (t(mu.q.Bpsi[,K]) %*% P.mat %*% mu.q.Bpsi[,K] + 
                                     sum(diag(P.mat %*% sigma.q.Bpsi[(Kt*(K-1)+1):(Kt*K),(Kt*(K-1)+1):(Kt*K)])))
    }
    
    ###############################################################
    ## lower bound
    ###############################################################
    
    curlpxq = 10
    lpxq = c(lpxq, curlpxq)
    j=j+1

    if(verbose) { cat(".") }

  }

  ## compute CI for fixed effects
  beta.sd = beta.LB = beta.UB = matrix(NA, nrow = p, ncol = D)
  for(i in 1:p){
    beta.sd[i,] = sqrt(diag((Theta) %*% sigma.q.BW[(Kt*(i-1)+1):(Kt*i),(Kt*(i-1)+1):(Kt*i)] %*% t(Theta)))
    beta.LB[i,] = beta.cur[i,]-1.96*beta.sd[i,]
    beta.UB[i,] = beta.cur[i,]+1.96*beta.sd[i,]
  }

  ## convert objects from spam to matrix
  beta.cur = as.matrix(beta.cur)
  psi.cur = as.matrix(psi.cur)
    
  ## effective degrees of freedom for fixed effects
  edf =  10 #sum(diag( (as.numeric((A + IJ*D/2)/(b.q.sigma.me))) * sigma.q.BW %*% t.designmat.X %*% t(t.designmat.X)  ))

  ## subj level random effects
  ranef.subj = as.matrix(t(mu.q.BZ) %*% t(Theta))

  ## do svd to get rotated fpca basis
  temp = svd(t(psi.cur))
  psi.cur = t(temp$u)
  lambda.pm = temp$d
  
  ## export fitted values
  Yhat.fixed = fixef.cur 
  Yhat.subj = fixef.cur + ranef.cur
  Yhat = fixef.cur + ranef.cur + pcaef.cur
  
  ## export variance components
  sigeps.pm = 1 / as.numeric((A + IJ*D/2)/(b.q.sigma.me))

  ## export various r2 values
  r2.f = 1 - (sum((Y - Yhat.fixed)^2)/(IJ*D)) / (sum((Y)^2)/(IJ*D))
  r2.fr = 1 - (sum((Y - Yhat.subj)^2)/(IJ*D)) / (sum((Y)^2)/(IJ*D))
  r2.frp = 1 - (sum((Y - Yhat)^2)/(IJ*D)) / (sum((Y)^2)/(IJ*D))
    
#  ret = list(beta.cur, beta.UB, beta.LB, psi.cur, pcaef.cur, ranef.cur, ranef.subj, Yhat, Yhat.subj, edf, 
#             sigeps.pm, r2.f, r2.fr, r2.frp, mu.q.BW, mu.q.BZ, mu.q.Bpsi, mu.q.C, 
#             1/((A + IJ*D/2)/(b.q.sigma.me)), 1/((A+I*Kt/2)/b.q.lambda.BZ), 1/((A+Kt/2)/b.q.lambda.BW), 1/((A+Kt/2)/b.q.lambda.Bpsi))
#  names(ret) = c("beta.pm", "beta.UB", "beta.LB", "psi.pm", "pcaef.cur", "ranef.pm", "ranef.subj", "Yhat", "Yhat.subj", "edf", 
#             "sigeps.pm", "r2.f", "r2.fr", "r2.frp", "BW.pm", "BZ.pm", "Bpsi.pm", "BC.pm", "sige.pm", "sigz.pm", "sigw.pm", "sigpsi.pm")

  ret = list(beta.cur, beta.UB, beta.LB, fixef.cur, mt_fixed, data, psi.cur, sigeps.pm, lambda.pm)
  names(ret) = c("beta.hat", "beta.UB", "beta.LB", "Yhat", "terms", "data", "psi.hat", "sigeps.pm", "lambda.pm")
  class(ret) = "fosr"
  ret

}
  

###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################