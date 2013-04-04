# Change lambdahat to eigenvals, eta.hat to mean, phiU.hat to eigenfuncs, etc. ???
# npc=1 seems to give error
fpca.sc <-
function(Y, Y.pred=NULL, nbasis = 10, pve = .99, npc = NULL, var = FALSE, simul = FALSE, sim.alpha = .95,
        useSymm = FALSE, makePD = FALSE){
  ## if Y.pred is not provided, use Y
  if (is.null(Y.pred)) Y.pred = Y
  
  D = NCOL(Y)                 # size of grid
  I = NROW(Y)                 # number of curves
  I.pred = NROW(Y.pred)       # number of curves for prediction
  d.vec = rep(1:D, each = I)  # grid on which curves are observed

  ## estimate the mean function
  gam0 = gam(as.vector(Y) ~ s(d.vec, k = nbasis))
  mu = predict(gam0, newdata = data.frame(d.vec = 1:D))

  ## de-mean bootstrap curves
  Y.tilde = Y - matrix(mu, I, D, byrow=TRUE)

  ## estimate covariance functions using method of moments
  cov.sum = cov.count = cov.mean = matrix(0, D, D)
  for (i in 1:I){
    obs.points = which(!is.na(Y[i,]))
    cov.count[obs.points, obs.points] = cov.count[obs.points, obs.points] + 1
    cov.sum[obs.points, obs.points] = cov.sum[obs.points, obs.points] + tcrossprod(Y.tilde[i, obs.points])
                                      # Y.tilde[i, obs.points ] %*% t(Y.tilde[i, obs.points ])  # modified by PR
  } 
  G.0 = ifelse(cov.count==0, NA, cov.sum/cov.count)  
  
  ## save raw diagonal elements, remove from diagonal from estimate for smoothing
  diag.G0 = diag(G.0)
  diag(G.0) = NA
  
  ## smooth the off-diagonal covariance elements:
  if(!useSymm){
      row.vec = rep(1:D, each = D)
      col.vec = rep(1:D, D)
      npc.0 = matrix(predict(gam(as.vector(G.0) ~ te(row.vec, col.vec, k = nbasis)),
                      newdata = data.frame(row.vec = row.vec, col.vec = col.vec)), D, D)
      npc.0 = (npc.0 + t(npc.0)) / 2        ## bivariate smoothing using the gam function 
  } else {
      #browser()
      
      #use upper tri to reduce computation (and allow kink at diagonal)
      use <- upper.tri(G.0, diag=TRUE)
      #use elements in upper and lower corners s.t. predict works for whole range
      use[2,1] <- use[ncol(G.0),ncol(G.0)-1] <- TRUE
      #modify weights accordingly s.t. these points don't affect fit:
      usecov.count <- cov.count
      usecov.count[2,1] <- usecov.count[ncol(G.0),ncol(G.0)-1] <- 0
      usecov.count <- as.vector(usecov.count)[use]
      
      use <- as.vector(use) 
      vG.0 <- as.vector(G.0)[use]
      row.vec <- rep(1:D, each = D)[use]
      col.vec <- rep(1:D, times = D)[use]
      
      mCov <- gam(vG.0 ~ te(row.vec, col.vec, k=nbasis), weights=usecov.count)    
      npc.0 <- matrix(NA, D, D)
      spred <- rep(1:D, each = D)[upper.tri(npc.0, diag=TRUE)]
      tpred <- rep(1:D, times = D)[upper.tri(npc.0, diag=TRUE)]
      smVCov <- predict(mCov, newdata=data.frame(row.vec=spred, col.vec=tpred))
      npc.0[upper.tri(npc.0, diag=TRUE)] <- smVCov 
      npc.0[lower.tri(npc.0)] <- t(npc.0)[lower.tri(npc.0)]
  }
  if(makePD){
      # project smoothed cov. surface into pos-semidef. cone
      npc.0 <- {
          tmp <- Matrix:::nearPD(npc.0,
                  corr = FALSE, 
                  keepDiag = FALSE, 
                  do2eigen = TRUE, # enforce positive def!
                  trace=TRUE
          )
          as.matrix(tmp$mat)                     
      }
  }
  
  ## estimate the score variances, truncation lag, and basis functions
  evalues = eigen(npc.0, symmetric = TRUE, only.values = TRUE)$values
  evalues = replace(evalues, which(evalues<=0), 0)
  npc = ifelse(is.null(npc), min(which(cumsum(evalues)/sum(evalues)>pve)), npc)
  
  efunctions = matrix(eigen(npc.0, symmetric = TRUE)$vectors[, seq(len=npc)], nrow = D, ncol = npc)
  evalues = eigen(npc.0, symmetric = TRUE, only.values = TRUE)$values[1:npc]
    
  ## estimate covariance matrix using retained basis functions
  cov.hat = efunctions %*% tcrossprod(diag(evalues, nrow = npc, ncol = npc), efunctions)
  # cov.hat = efunctions %*% tcrossprod(diag(evalues), efunctions)
  DIAG = (diag.G0 - diag(cov.hat))[floor(D*0.2):ceiling(D*0.8)]
  
  ## estimate measurement error variance; stop if the estimate is 0
  sigma2 = max(mean(DIAG, na.rm=TRUE), 0)  

  ## construct design matrices for mixed model estimation
  D.inv = diag(1/evalues, nrow = npc, ncol = npc)  # commented out by PR
  Z = efunctions
  Y.tilde = Y.pred - matrix(mu, I.pred, D, byrow=TRUE)
  
  ## create matrices to store results
  Yhat = matrix(0, nrow = I.pred, ncol = D)
  scores = matrix(NA, nrow = I.pred, ncol = npc)
  VarMats = vector("list", I.pred)
  for (i in 1:I.pred) VarMats[[i]] = matrix(NA, nrow = D, ncol = D)
  diag.var = matrix(NA, nrow = I.pred, ncol = D)
  crit.val = rep(0, I.pred)

  ## estimate scores, functions and variances for all curves
  for (i.subj in 1:I.pred) {
    obs.points = which(!is.na(Y.pred[i.subj,]))
    if(sigma2 ==0 & length(obs.points) < npc){
      stop("Measurement error estimated to be zero and there are fewer observed points than PCs; scores cannot be estimated.")
    }
    Zcur = matrix(Z[obs.points,], nrow = length(obs.points), ncol = dim(Z)[2])  # PR: Is this ncol just npc??
    ZtZ_sD.inv = solve(crossprod(Zcur) + sigma2 * D.inv)

    scores[i.subj,] = ZtZ_sD.inv %*% t(Zcur) %*% (Y.tilde[i.subj,obs.points])
    Yhat[i.subj,] = t(as.matrix(mu)) + scores[i.subj,]%*%t(efunctions)
  
    if (var) {
      VarMats[[i.subj]] = sigma2 * Z %*% ZtZ_sD.inv %*%t(Z)
      diag.var[i.subj,] = diag(VarMats[[i.subj]])
  
      ## estimate critical values for simultaneous intervals
      if (simul & sigma2 != 0) {
        norm.samp = mvrnorm(2500, mu = rep(0, D), Sigma = VarMats[[i.subj]])/
          matrix(sqrt(diag(VarMats[[i.subj]])), nrow = 2500, ncol = D, byrow = TRUE)
        crit.val[i.subj] = quantile(apply(abs(norm.samp), 1, max), sim.alpha)      
      }
    }
  }

  ret.objects = c("Yhat", "scores", "mu", "efunctions", "evalues", "npc")
  if (var) {
      ret.objects = c(ret.objects, "sigma2", "diag.var", "VarMats")
      if (simul) ret.objects = c(ret.objects, "crit.val") 
  }
  
  ret = lapply(1:length(ret.objects), function(u) get(ret.objects[u]))
  names(ret) = ret.objects
  
  return(ret)
}

