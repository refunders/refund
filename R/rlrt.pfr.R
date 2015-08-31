##' Likelihood Ratio Test and Restricted Likelihood Ratio Test for inference of
##' functional predictors
##'
##' NOTE: this function is designed to work with pfr_old() rather than pfr().
##' Given a pfr object of family="gaussian", tests whether the function is
##' identically equal to its mean (constancy), or whether the functional
##' predictor significantly improves the model (inclusion).  Based on
##' zero-variance-component work of Crainiceanu et al. (2004), Scheipl et al.
##' (2008), and Swihart et al. (2012).
##'
##' A Penalized Functional Regression of family="gaussian" can be represented
##' as a linear mixed model dependent on variance components. Testing whether
##' certain variance components and (potentially) fixed effect coefficients are
##' 0 correspond to tests of constancy and inclusion of functional predictors.
##'
##' For rlrt.pfr, Restricted Likelihood Ratio Test is preferred for the
##' constancy test as under the special B-splines implementation of pfr for the
##' coefficient function basis the test involves only the variance component.
##' Therefore, the constancy test is best for pfr objects with method="REML";
##' if the method was something else, a warning is printed and the model refit
##' with "REML" and a test is then conducted.
##'
##' For rlrt.pfr, the Likelihood Ratio Test is preferred for the inclusion test
##' as under the special B-splines implementation of pfr for the coefficient
##' function basis the test involves both the variance component and a fixed
##' effect coefficient in the linear mixed model representation. Therefore, the
##' inclusion test is best for pfr objects with method="ML"; if the method was
##' something else, a warning is printed and the model refit with "ML" and a
##' test is then conducted.
##'
##' @param pfr.obj an object returned by pfr_old()
##' @param test "constancy" will test functional form of the coefficient
##' function of the last function listed in funcs in pfr.obj against the null
##' of a constant line: the average of the functional predictor.  "inclusion"
##' will test functional form of the coefficient function of the last function
##' listed in funcs in pfr.obj against the null of 0: that is, whether the
##' functional predictor should be included in the model.
##' @param ... additional arguments
##' @return \item{p.val }{the p-value for the full model (alternative) against
##' the null specified by the test} \item{test.stat }{the test statistic, see
##' Scheipl et al. 2008 and Swihart et al 2012}
##'
##' \item{ma }{the alternative model as fit with mgcv::gam} \item{m0 }{the null
##' model as fit with mgcv::gam} \item{m }{the model containing only the
##' parameters being tested as fit with mgcv::gam}
##' @author Jeff Goldsmith <jeff.goldsmith@@columbia.edu> and Bruce Swihart
##' <bswihart@@jhsph.edu>
##' @seealso \code{\link{pfr}}, \code{\link{predict.pfr}}, package
##' \code{RLRsim}
##' @export
##' @importFrom RLRsim RLRTSim LRTSim
##' @importFrom mgcv summary.gam gam
##' @references Goldsmith, J., Bobb, J., Crainiceanu, C., Caffo, B., and Reich,
##' D. (2011). Penalized functional regression. \emph{Journal of Computational
##' and Graphical Statistics}, 20(4), 830--851.
##'
##' Goldsmith, J., Crainiceanu, C., Caffo, B., and Reich, D. (2012).
##' Longitudinal penalized functional regression for cognitive outcomes on
##' neuronal tract measurements. \emph{Journal of the Royal Statistical
##' Society: Series C}, 61(3), 453--469.
##'
##' Crainiceanu, C. and Ruppert, D. (2004) Likelihood ratio tests in linear
##' mixed models with one variance component. \emph{Journal of the Royal
##' Statistical Society: Series B}, 66, 165--185.
##'
##' Scheipl, F. (2007) Testing for nonparametric terms and random effects in
##' structured additive regression. Diploma thesis.\
##' http://www.statistik.lmu.de/~scheipl/downloads/DIPLOM.zip.
##'
##' Scheipl, F., Greven, S. and Kuechenhoff, H (2008) Size and power of tests
##' for a zero random effect variance or polynomial regression in additive and
##' linear mixed models.  \emph{Computational Statistics & Data Analysis},
##' 52(7), 3283--3299.
##'
##' Swihart, Bruce J., Goldsmith, Jeff; and Crainiceanu, Ciprian M. (2012).
##' Testing for functional effects. Johns Hopkins University Dept. of
##' Biostatistics Working Paper 247. Available at
##' \url{http://biostats.bepress.com/jhubiostat/paper247}
##' @examples
##'
##' \dontrun{
##' ##################################################################
##' #########               DTI Data Example                 #########
##' ##################################################################
##'
##' ##################################################################
##' # For more about this example, see Swihart et al. 2012
##' # Testing for Functional Effects
##' ##################################################################
##'
##' ## load and reassign the data;
##' data(DTI2)
##' O  <- DTI2$pasat ## PASAT outcome
##' id <- DTI2$id    ## subject id
##' W1 <- DTI2$cca   ## Corpus Callosum
##' W2 <- DTI2$rcst  ## Right corticospinal
##' V  <- DTI2$visit ## visit
##'
##' ## prep scalar covariate
##' visit.1.rest <- matrix(as.numeric(V > 1), ncol=1)
##' covar.in <- visit.1.rest
##'
##'
##' ## note there is missingness in the functional predictors
##' apply(is.na(W1), 2, mean)
##' apply(is.na(W2), 2, mean)
##'
##' ## fit two univariate models, then one model with both functional predictors
##' pfr.obj.t1 <- pfr_old(Y = O, covariates=covar.in, funcs = list(W1),     subj = id, kz = 10, kb = 50)
##' pfr.obj.t2 <- pfr_old(Y = O, covariates=covar.in, funcs = list(W2),     subj = id, kz = 10, kb = 50)
##' pfr.obj.t3 <- pfr_old(Y = O, covariates=covar.in, funcs = list(W1, W2), subj = id, kz = 10, kb = 50)
##'
##' ## plot the coefficient function and bounds
##' dev.new()
##' par(mfrow=c(2,2))
##' ran <- c(-2,.5)
##' matplot(cbind(pfr.obj.t1$BetaHat[[1]], pfr.obj.t1$Bounds[[1]]),
##'   type = 'l', lty = c(1,2,2), col = c(1,2,2), ylab = "BetaHat",
##'   main = "CCA", xlab="Location", ylim=ran)
##' abline(h=0, col="blue")
##' matplot(cbind(pfr.obj.t2$BetaHat[[1]], pfr.obj.t2$Bounds[[1]]),
##'   type = 'l', lty = c(1,2,2), col = c(1,2,2), ylab = "BetaHat",
##'   main = "RCST", xlab="Location", ylim=ran)
##' abline(h=0, col="blue")
##' matplot(cbind(pfr.obj.t3$BetaHat[[1]], pfr.obj.t3$Bounds[[1]]),
##'   type = 'l', lty = c(1,2,2), col = c(1,2,2), ylab = "BetaHat",
##'   main = "CCA  - mult.", xlab="Location", ylim=ran)
##' abline(h=0, col="blue")
##' matplot(cbind(pfr.obj.t3$BetaHat[[2]], pfr.obj.t3$Bounds[[2]]),
##'   type = 'l', lty = c(1,2,2), col = c(1,2,2), ylab = "BetaHat",
##'   main = "RCST - mult.", xlab="Location", ylim=ran)
##' abline(h=0, col="blue")
##'
##' ## do some testing
##' t1 <- rlrt.pfr(pfr.obj.t1, "constancy")
##' t2 <- rlrt.pfr(pfr.obj.t2, "constancy")
##' t3 <- rlrt.pfr(pfr.obj.t3, "inclusion")
##'
##' t1$test.stat
##' t1$p.val
##'
##' t2$test.stat
##' t2$p.val
##'
##' t3$test.stat
##' t3$p.val
##'
##'
##' ## do some testing with rlrt.pfr(); same as above but subj = NULL
##' pfr.obj.t1 <- pfr(Y = O, covariates=covar.in, funcs = list(W1),     subj = NULL, kz = 10, kb = 50)
##' pfr.obj.t2 <- pfr(Y = O, covariates=covar.in, funcs = list(W2),     subj = NULL, kz = 10, kb = 50)
##' pfr.obj.t3 <- pfr(Y = O, covariates=covar.in, funcs = list(W1, W2), subj = NULL, kz = 10, kb = 50)
##'
##' t1 <- rlrt.pfr(pfr.obj.t1, "constancy")
##' t2 <- rlrt.pfr(pfr.obj.t2, "constancy")
##' t3 <- rlrt.pfr(pfr.obj.t3, "inclusion")
##'
##' t1$test.stat
##' t1$p.val
##'
##' t2$test.stat
##' t2$p.val
##'
##' t3$test.stat
##' t3$p.val
##' }
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

      warning("warning: test 'constancy' is a RLRT for bsplines and a refit of the full (alternative) model was completed with method = 'REML'")

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

      warning("test 'confounding' is a LRT and a refit of the full (alternative) model was completed with method = 'ML'")

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




