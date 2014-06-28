rlrt.pfr <- function (pfr.obj=pfr.obj, test=NULL, ...)
{
  if(is.null(test) || !(test %in% c("constancy","inclusion")) ){
    print("test must be 'constancy' or 'inclusion'")
    break;
  }
  ## parse pfr.obj and pfr.obj$fit slots so that code works
  Y <- pfr.obj$Y
  fit <- pfr.obj$fit
  fitted.vals <- pfr.obj$fitted.vals
  beta.covariates <- pfr.obj$beta.covariates
  BetaHat <- BetaHat.ma <- pfr.obj$BetaHat
  totD <- pfr.obj$totD
  X <- pfr.obj$X
  D <- pfr.obj$D
  kb <- pfr.obj$kb
  p <- pfr.obj$p
  N_subj <- pfr.obj$N_subj
  subj <- pfr.obj$subj
  CJ <- pfr.obj$CJ
  N.Pred <- pfr.obj$N.Pred
  phi <- pfr.obj$phi
  fixed.mat <- pfr.obj$fixed.mat
  rand.mat <- pfr.obj$rand.mat


  family <- fit$family

  ## if user-specified, one of four tests
  if(test=="constancy"){
    ## constancy is a RLRT.  Refit full model under method=REML if necessary
    ##        BetaHat.ma <- varBeta <- varBetaHat <- Bounds <- list()
    BetaHat.ma <- list()
    if(fit$method!="REML"){

      print("warning: test 'constancy' is a RLRT for bsplines and a refit of the full (alternative) model")
      print("was completed with method = 'REML'")

      ## (re)fit alternative under method="ML"
      ma = gam(Y ~ X - 1, paraPen = list(X = D), method = "REML", family = family)#, ...)
      logLik.ma <- -summary(ma)$sp.criterion
      coefs.ma = ma$coef
      fitted.vals.ma <- as.matrix(X[, 1:length(coefs.ma)]) %*% coefs.ma
      beta.covariates.ma = coefs.ma[1:(p + 1)]
      for(i in 1:N.Pred){
        BetaHat.ma[[i]] = phi[[i]] %*% coefs.ma[-1*(1:(N_subj+p+1))][((i-1)*kb+1):(kb*i)]
      }
    }else{## fit=REML is specified by user for test="constancy"
      ma                 <- fit
      logLik.ma          <- -summary(ma)$sp.criterion
      coefs.ma           =  ma$coef
      fitted.vals.ma     <- fitted.vals
      beta.covariates.ma =  beta.covariates
      ##for(i in 1:N.Pred){ ## taken care of above in parsing pfr.obj
      ##  BetaHat.ma[[i]] = BetaHat[[i]]
      ##}

    }
    ## .m0 is the null model; that is the variance component
    ## for CJ[[N.Pred]] (random effs, but not fixed effect are 0) which
    ## correspondingly eliminates ALL  random CJ[[N.Pred]] coeffs.  So
    ## we shave off the end of X for X.m0
    X.m0 <- X[, -1*((ncol(X)-kb+2):ncol(X))]
    D.m0 <- list()
    ##for (i in 1:(totD-1)) {
    for(i in 1:(totD + is.null(subj)-1)){## adjustment to totD for subj==NULL
    D.m0[[i]] <- D[[i]][-1*((ncol(X)-kb+2):ncol(X)), -1*((ncol(X)-kb+2):ncol(X))]
    }
    ## detect if any of penalty matrices are all 0.
    all0 <- rep(NA, length(D.m0))
    for(i in 1:length(D.m0)){ all0[i] <- sum(D.m0[[i]]==0) == nrow(D.m0[[i]])*ncol(D.m0[[i]])}
    if(any(all0)){## if all 0 penalty, don't specify paraPen
      m0 <- gam(Y ~ X.m0 - 1, method = "REML", family = family)#, ...)
    }else{m0 <- gam(Y ~ X.m0 - 1, paraPen = list(X.m0=D.m0), method = "REML", family = family)#, ...)
        }
    logLik.m0 <- -summary(m0)$sp.criterion
    coefs.m0 = m0$coef
    fitted.vals.m0 <- as.matrix(X[, 1:length(coefs.m0)]) %*% coefs.m0
    beta.covariates.m0 = coefs.m0[1:(p + 1)]
    BetaHat.m0 <- list()
    if(N.Pred > 1){
      for(i in 1:(N.Pred-1)){
        BetaHat.m0[[i]] = phi[[i]] %*% coefs.m0[-1*(1:(N_subj+p+1))][((i-1)*kb+1):(kb*i)]
      }
    }else{BetaHat.m0 <- phi[[i]][,1] *coefs.m0[-1*(1:(N_subj+p+1))][1]}
    BetaHat.null <- BetaHat.m0 ## cheat for now...; doesn't affect testing, just returned coeff. func.
    ## .m is the model that only contains the variance components
    ## being tested; that is the CJ[[N.Pred]] which correspondingly
    ## eliminates the other CJ[[]] and Z1
    X.m <- cbind(X[,1:(1+p)], CJ[[N.Pred]])
    D.m <- list(length=1)
    D.m[[1]] <- as.matrix(D[[totD]][c(1:(1+p),(ncol(X)-kb+1):ncol(X)),c(1:(1+p),(ncol(X)-kb+1):ncol(X)) ])
    m <- gam(Y ~ X.m - 1, paraPen = list(X.m=D.m), method = "REML", family = family)##, ...)
    logLik.m <--summary(m)$sp.criterion
    ##calculate the observed RLRT statistic
    (rlrt.obs = mean(max(0, 2*logLik.ma - 2*logLik.m0 )))
    ## Subtle point:  need the X.m and other .m info for LRT sampling.  This is
    ## an approximation to the exact, but correctly follows Schiepl and Greven
    fixed.mat <- fixed.mat[,c(1:(p+1),ncol(fixed.mat))]
    rand.mat  <- rand.mat[,(ncol(rand.mat)-kb+2):ncol(rand.mat)]
    rand.pen  <- D.m[[1]][-1*c(1:(2+p)),-1*c(1:(2+p))]
    sample = RLRTSim(fixed.mat, rand.mat, qrX=qr(X), sqrt.Sigma = chol(cov2cor(rand.pen)),
      seed = NA, nsim = 10000, log.grid.hi = 8, log.grid.lo = -10, gridlength = 200)
    test.stat <- rlrt.obs
    (p.val = mean(rlrt.obs < sample))
    (p.val.sl = mean(rlrt.obs < cbind(rchisq(5000,0), rchisq(5000,1))))

  }
  if(test=="inclusion"){
    ## inclusion is a LRT.  Refit full model under method=ML if necessary
    if(fit$method=="REML"){

      print("warning: test 'confounding' is a LRT and a refit of the full (alternative) model")
      print("was completed with method = 'ML'")

      ## (re)fit alternative under method="ML"
      ma = gam(Y ~ X - 1, paraPen = list(X = D), method = "ML", family = family)#, ...)
      logLik.ma <- -summary(ma)$sp.criterion
      coefs.ma = ma$coef
      fitted.vals.ma <- as.matrix(X[, 1:length(coefs.ma)]) %*% coefs.ma
      beta.covariates.ma = coefs.ma[1:(p + 1)]
      BetaHat.ma <- list()
      for(i in 1:N.Pred){
        BetaHat.ma[[i]] = phi[[i]] %*% coefs.ma[-1*(1:(N_subj+p+1))][((i-1)*kb+1):(kb*i)]
      }
    }else{## fit=ML is specified by user along with test="confounding"
      ma                 <- fit
      logLik.ma          <- -summary(ma)$sp.criterion
      coefs.ma           =  ma$coef
      fitted.vals.ma     <- fitted.vals
      beta.covariates.ma =  beta.covariates
      ##for(i in 1:N.Pred){ ## taken care of above in parsing pfr.obj
      ##  BetaHat.ma[[i]] = BetaHat[[i]]
      ##}
    }
    ## .m0 is the null model; that is the variance component
    ## for CJ[[N.Pred]] (random effs as well as fixed) which
    ## correspondingly eliminates ALL CJ[[N.Pred]] coeffs.  So
    ## we shave off the end of X for X.m0
    X.m0 <- X[, -1*((ncol(X)-kb+1):ncol(X))]
    D.m0 <- list()
    ##for (i in 1:(totD-1)) {
    for(i in 1:(totD + is.null(subj)-1)){## adjustment to totD for subj==NULL
    D.m0[[i]] <- D[[i]][-1*((ncol(X)-kb+1):ncol(X)), -1*((ncol(X)-kb+1):ncol(X))]
    }
    ## detect if any of penalty matrices are all 0.
    all0 <- rep(NA, length(D.m0))
    for(i in 1:length(D.m0)){ all0[i] <- sum(D.m0[[i]]==0) == nrow(D.m0[[i]])*ncol(D.m0[[i]])}
    if(any(all0)){## if all 0 penalty, don't specify paraPen
      m0 <- gam(Y ~ X.m0 - 1, method = "ML", family = family)#, ...)
    }else{m0 <- gam(Y ~ X.m0 - 1, paraPen = list(X.m0=D.m0), method = "ML", family = family)#, ...)
        }
    logLik.m0 <- -summary(m0)$sp.criterion
    coefs.m0 = m0$coef
    fitted.vals.m0 <- as.matrix(X[, 1:length(coefs.m0)]) %*% coefs.m0
    beta.covariates.m0 = coefs.m0[1:(p + 1)]
    BetaHat.m0 <- list()
    if(N.Pred > 1){
      for(i in 1:(N.Pred-1)){
        BetaHat.m0[[i]] = phi[[i]] %*% coefs.m0[-1*(1:(N_subj+p+1))][((i-1)*kb+1):(kb*i)]
      }
    }else{BetaHat.m0 <- NULL}
    ## .m is the model that only contains the variance components
    ## being tested; that is the CJ[[N.Pred]] which correspondingly
    ## eliminates the other CJ[[]] and Z1
    X.m <- cbind(X[,1:(1+p)], CJ[[N.Pred]])
    D.m <- list(length=1)
    D.m[[1]] <- as.matrix(D[[totD]][c(1:(1+p),(ncol(X)-kb+1):ncol(X)),c(1:(1+p),(ncol(X)-kb+1):ncol(X)) ])
    m <- gam(Y ~ X.m - 1, paraPen = list(X.m=D.m), method = "ML", family = family)##, ...)
    logLik.m <--summary(m)$sp.criterion
    ##calculate the observed LRT statistic
    (lrt.obs = mean(max(0, 2*logLik.ma - 2*logLik.m0 )))
    ## Subtle point:  need the X.m and other .m info for LRT sampling.  This is
    ## an approximation to the exact, but correctly follows Schiepl and Greven
    fixed.mat <- fixed.mat[,c(1:(p+1),ncol(fixed.mat))]
    rand.mat  <- rand.mat[,(ncol(rand.mat)-kb+2):ncol(rand.mat)]
    rand.pen  <- D.m[[1]][-1*c(1:(2+p)),-1*c(1:(2+p))]
    sample = LRTSim(fixed.mat, rand.mat, q=1, sqrt.Sigma = chol(cov2cor(rand.pen)),
      seed = NA, nsim = 10000, log.grid.hi = 8, log.grid.lo = -10, gridlength = 200)
    test.stat <- lrt.obs
    (p.val = mean(lrt.obs < sample))
  }
  ## return pertinent results
  ## ret <-     list( BetaHat ,  Bounds ,  BetaHat.null ,  p.val ,  test.stat,   ma ,  m0 ,  m)
  ## names(ret) <- c("BetaHat", "Bounds", "BetaHat.null", "p.val", "test.stat", "ma", "m0", "m")
  ## ret
  ret <-     list( p.val ,  test.stat,   ma ,  m0 ,  m)
  names(ret) <- c("p.val", "test.stat", "ma", "m0", "m")
  ret
}




