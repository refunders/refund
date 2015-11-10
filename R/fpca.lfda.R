#' Longtiudinal Functional Data Analysis using FPCA
#'
#' Implements longitudinal functional data analysis (Park and Staicu, 2015).
#' It decomposes longitudinally-observed functional observations in two steps.
#' It first applies FPCA on a properly defined marginal covariance function and obtain estimated scores (mFPCA step).
#' Then it further models the underlying process dynamics by applying another FPCA on a covariance of the estimated scores
#' obtained in the mFPCA step. The function also allows to use a random effects model to study the underlying process dynamics
#' instead of a KL expansion model in the second step. Scores in mFPCA step are estimated
#' using numerical integration. Scores in sFPCA step are estimated under a mixed model framework.
#' 
#' @section Details:
#' A random effects model is recommended when a set of visit times for all subjects and visits is not dense in its range.
#' 
#' @param Y a matrix of which each row corresponds to one curve observed on a regular and dense grid 
#'          (dimension of N by m; N = total number of observed functions; m = number of grid points)
#' @param subject.index subject id; vector of length N with each element corresponding a row of Y
#' @param visit.index index for visits (repeated measures); vector of length N with each element corresponding a row of Y
#' @param obsT actual time of visits at which a function is observed; vector of length N with each element corresponding a row of Y
#' @param funcArg numeric; function argument
#' @param numTEvalPoint total number of evaluation time points for visits; used for pre-binning in sFPCA step; defaults to 41
#' @param newdata an optional data frame providing predictors (i for subject id / Ltime for visit time) with which prediction is desired; defaults to NULL
#' 
#' @param fbps.knots list of two vectors of knots or number of equidistanct knots for all dimensions
#'                   for a fast bivariate \emph{P}-spline smoothing (fbps) method used to estimate a bivariate, smooth mean function; defauls to c(5,10); see \command{fbps} 
#' @param fbps.p integer;degrees of B-spline functions to use for a fbps method; defaults to 3; see \command{fbps} 
#' @param fbps.m integer;order of differencing penalty to use for a fbps method; defaults to 2; see \command{fbps} 
#' 
#' @param mFPCA.pve proportion of variance explained for a mFPCA step; used to choose the number of principal components (PCs); defaults to 0.95; see \command{fpca.face} 
#' @param mFPCA.knots number of knots to use or the vectors of knots in a mFPCA step; used for obtain a smooth estimate of a covarance function; defaults to 35; see \command{fpca.face} 
#' @param mFPCA.p integer; the degree of B-spline functions to use in a mFPCA step; defaults to 3; see \command{fpca.face} 
#' @param mFPCA.m integer;order of differencing penalty to use in a mFPCA step; defaults to 2; see \command{fpca.face} 
#' 
#' @param LongiModel.method model and estimation method for estimating covariance of estimated scores from a mFPCA step; 
#'                          either KL expansion model or random effects model; defaults to fpca.sc
#' @param sFPCA.pve proportion of variance explained for sFPCA step; used to choose the number of principal components; defaults to 0.95; see \command{fpca.sc} 
#' @param sFPCA.nbasis number of B-spline basis functions used in sFPCA step for estimation of the mean function and bivariate smoothing of the covariance surface; defaults to 10; see \command{fpca.sc} 
#' @param sFPCA.npc pre-specified value for the number of principal components; if given, it overrides \code{pve}; defaults to NULL; see \command{fpca.sc} 
#' 
#' @param gam.method smoothing parameter estimation method when \command{gam} is used for predicting score functions at unobserved visit time, T; defaults to \code{REML}; see \command{gam} 
#' @param gam.kT dimension of basis functions to use; see \command{gam} 
#' 
#' @return A list with components \item{i }{subject id} \item{funcArg }{function argument} \item{visitTime }{visit times}
#' \item{fitted.values }{fitted values (in-sample); of the same dimension as Y} \item{fitted.values.all }{a list of which each component consists of a subject's fitted values at all pairs of evaluation points (s and T)} 
#' \item{predicted.values }{predicted values for variables provided in newdata}
#' \item{bivariateSmoothMeanFunc }{estimated bivariate smooth mean function}
#' \item{mFPCA.efunctions }{estimated eigenfunction in a mFPCA step} \item{mFPCA.evalues }{estimated eigenvalues in a mFPCA step}
#' \item{mFPCA.npc }{number of principal components selected with pre-specified pve in a mFPCA step} \item{mFPCA.scree.eval }{estiamted eigenvalues obtained with pre-specified pve = 0.9999; for scree plot}
#' \item{sFPCA.xiHat.bySubj }{a list of which each component consists of a subject's predicted score functions evaluated at equidistanced grid in direction of visit time, T} 
#' \item{sFPCA.npc}{a vector of numbers of principal components selected in a sFPCA step with pre-specified pve; length of mFPCA.npc}
#' \item{mFPCA.covar }{estimated marginal covariance} \item{sFPCA.longDynCov.k }{a list of estimated covariance of score function; length of mFPCA.npc}
#' 
#' @author So Young Park \email{spark13@@ncsu.edu}, Ana-Maria Staicu
#' @export
#' @references Park, S.Y. and Staicu, A.M. (2015). Longitudinal functional data analysis. Stat 4 212-226.
#'
#' 
#' @importFrom lme4 lmer
#' @importFrom mgcv gam s 
#' @importFrom splines spline.des
#' @importFrom Matrix kronecker as.matrix
#'  
#' @examples 
#'
#'   ###########################################################################################
#'   # data generation
#'   ###########################################################################################
#'   library(rgl)
#'   set.seed(1)
#'   n <- 100 # number of subjects
#'   ss <- seq(0,1,length.out=101) 
#'   TT <- seq(0, 1, length.out=41)
#'   mi <- runif(n, min=6, max=15)
#'   ij <- sapply(mi, function(a) sort(sample(1:41, size=a, replace=FALSE)))
#'   
#'   # error variances
#'   sigma <- 0.1 
#'   sigma_wn <- 0.2
#'
#'   lambdaTrue <- c(1,0.5) # True eigenvalues
#'   eta1True <- c(0.5, 0.5^2, 0.5^3) # True eigenvalues
#'   eta2True <- c(0.5^2, 0.5^3) # True eigenvalues
#'   
#'   phi <- sqrt(2)*cbind(sin(2*pi*ss),cos(2*pi*ss))
#'   psi1 <- cbind(rep(1,length(TT)), sqrt(3)*(2*TT-1), sqrt(5)*(6*TT^2-6*TT+1))
#'   psi2 <- sqrt(2)*cbind(sin(2*pi*TT),cos(2*pi*TT))
#'   
#'   zeta1 <- sapply(eta1True, function(a) rnorm(n = n, mean = 0, sd = a))
#'   zeta2 <- sapply(eta2True, function(a) rnorm(n = n, mean = 0, sd = a))
#'   
#'   xi1 <- unlist(lapply(1:n, function(a) (zeta1 %*% t(psi1))[a,ij[[a]]] ))
#'   xi2 <- unlist(lapply(1:n, function(a) (zeta2 %*% t(psi2))[a,ij[[a]]] ))
#'   xi <- cbind(xi1, xi2)
#'   
#'   Tij <- unlist(lapply(1:n, function(i) TT[ij[[i]]] ))
#'   i <- unlist(lapply(1:n, function(i) rep(i, length(ij[[i]]))))
#'   j <- unlist(lapply(1:n, function(i) 1:length(ij[[i]])))
#'   
#'   X <- xi %*% t(phi)
#'   meanFn <- function(s,t){ 0.5*t + 1.5*s + 1.3*s*t}
#'   mu <- matrix(meanFn(s = rep(ss, each=length(Tij)), t=rep(Tij, length(ss)) ) , nrow=nrow(X))
#'
#'   Y <- mu +  X + 
#'      matrix(rnorm(nrow(X)*ncol(phi), 0, sigma), nrow=nrow(X)) %*% t(phi) + # correlated error process
#'      matrix(rnorm(length(X), 0, sigma_wn), nrow=nrow(X)) # white noise
#'
#'   matplot(ss, t(Y[which(i==2),]), type='l', ylab="", xlab="functional argument", main="observations from subject i = 2")
#'   # END: data generation
#'   
#'   ###########################################################################################
#'   # Illustration I : when covariance of scores from a mFPCA step is estimated using fpca.sc
#'   ###########################################################################################
#'   est <- fpca.lfda(Y = Y, 
#'                    subject.index = i,
#'                    visit.index = j,
#'                    obsT = Tij,
#'                    funcArg = ss,
#'                    numTEvalPoints = length(TT), newdata = data.frame(i = c(1:3), Ltime = c(Tij[1], 0.2, 0.5)), 
#'                    fbps.knots = 35, fbps.p = 3, fbps.m = 2,
#'                    LongiModel.method='fpca.sc',
#'                    mFPCA.pve = 0.95, mFPCA.knots = 35, mFPCA.p = 3, mFPCA.m = 2, mFPCA.npc = NULL,
#'                    sFPCA.pve = 0.95, sFPCA.nbasis = 10, sFPCA.npc = NULL,
#'                    gam.method = 'REML', gam.kT = 10)
#'   
#'   
#'   # mean function (true vs. estimated)
#'   persp3d(x=TT, y = ss, z= t(sapply(TT, function(a) meanFn(s=ss, t = a))),xlab="visit times", ylab="s", zlab="estimated mean fn")
#'   persp3d(x = TT, y = ss, est$bivariateSmoothMeanFunc, add = TRUE, col='light blue')
#'   
#'   ################   mFPCA step   ################
#'   par(mfrow=c(1,2))
#'   
#'   # marginal covariance fn (true vs. estimated)
#'   image(phi%*%diag(lambdaTrue)%*%t(phi))
#'   image(est$mFPCA.covar) 
#'   
#'   # eigenfunctions (true vs. estimated)
#'   matplot(ss, phi, type='l') 
#'   matlines(ss, cbind(est$mFPCA.efunctions[,1], est$mFPCA.efunctions[,2])  , type='l', lwd=2)
#'   
#'   # scree plot
#'   plot(cumsum(est$mFPCA.scree.eval)/sum(est$mFPCA.scree.eval), type='l', ylab = "Percentage of variance explained")
#'   points(cumsum(est$mFPCA.scree.eval)/sum(est$mFPCA.scree.eval), pch=16)
#'   
#'   ################   sFPCA step   ################
#'   par(mfrow=c(1,2))
#'   print(est$mFPCA.npc)  # k = 2
#'   
#'   # covariance of score functions for k = 1 (true vs. estimated)
#'   image(psi1%*%diag(eta1True)%*%t(psi1), main='TRUE')
#'   image(est$sFPCA.longDynCov.k[[1]], main='ESTIMATED')
#'   
#'   # covariance of score functions for k = 2 (true vs. estimated)
#'   image(psi2%*%diag(eta2True)%*%t(psi2))
#'   image(est$sFPCA.longDynCov.k[[2]])
#'   
#'   # estimated scores functions
#'   matplot(TT, do.call(cbind,lapply(est$sFPCA.xiHat.bySubj, function(a) a[,1])), xlab="visit time", main="k=1", type='l', ylab="", col=rainbow(100, alpha = 1), lwd=1, lty=1)
#'   matplot(TT, do.call(cbind,lapply(est$sFPCA.xiHat.bySubj, function(a) a[,2])), xlab="visit time", main="k=2",type='l', ylab="", col=rainbow(100, alpha = 1), lwd=1, lty=1)
#'   
#'   ################   In-sample and Out-of-sample Prediction   ################
#'   par(mfrow=c(1,2))
#'   # fitted
#'   matplot(ss, t(Y[which(i==1),]), type='l', ylab="", xlab="functional argument")
#'   matlines(ss, t(est$fitted.values[which(i==1),]), type='l', lwd=2)
#'   
#'  # sanity check : expect fitted and predicted (obtained using info from newdata) values to be the same
#'   plot(ss, est$fitted.values[1,], type='p', xlab="", ylab="", pch = 1, cex=1)
#'   lines(ss, est$predicted.values[1,], type='l', lwd=2, col='blue')
#'   all.equal(est$predicted.values[1,], est$fitted.values[1,])
#'   
#'   ###########################################################################################
#'   # Illustration II : when covariance of scores from a mFPCA step is estimated using lmer
#'   ###########################################################################################
#'   est.lme <- fpca.lfda(Y = Y, 
#'                        subject.index = i,
#'                        visit.index = j,
#'                        obsT = Tij,
#'                        funcArg = ss,
#'                        numTEvalPoints = length(TT), newdata = data.frame(i = c(1:3), Ltime = c(Tij[1], 0.2, 0.5)), 
#'                        fbps.knots = 35, fbps.p = 3, fbps.m = 2,
#'                        LongiModel.method='lme',
#'                        mFPCA.pve = 0.95, mFPCA.knots = 35, mFPCA.p = 3, mFPCA.m = 2, mFPCA.npc = NULL,
#'                        gam.method = 'REML', gam.kT = 10)
#'   
#'   par(mfrow=c(2,2))
#'   
#'   # fpca.sc vs. lme (assumes linearity)
#'   matplot(TT, do.call(cbind,lapply(est$sFPCA.xiHat.bySubj, function(a) a[,1])), xlab="visit time", main="k=1", type='l', ylab="", col=rainbow(100, alpha = 1), lwd=1, lty=1)
#'   matplot(TT, do.call(cbind,lapply(est$sFPCA.xiHat.bySubj, function(a) a[,2])), xlab="visit time", main="k=2",type='l', ylab="", col=rainbow(100, alpha = 1), lwd=1, lty=1)
#'   matplot(TT, do.call(cbind,lapply(est.lme$sFPCA.xiHat.bySubj, function(a) a[,1])), xlab="visit time", main="k=1", type='l', ylab="", col=rainbow(100, alpha = 1), lwd=1, lty=1)
#'   matplot(TT, do.call(cbind,lapply(est.lme$sFPCA.xiHat.bySubj, function(a) a[,2])), xlab="visit time", main="k=2", type='l', ylab="", col=rainbow(100, alpha = 1), lwd=1, lty=1)

#####################################################################################################################
#####################################################################################################################

fpca.lfda <- function(Y, subject.index, visit.index, obsT = NULL, funcArg = NULL, numTEvalPoints = 41, newdata = NULL, 
                      fbps.knots = c(5,10), fbps.p = 3, fbps.m = 2,
                      mFPCA.pve = 0.95, mFPCA.knots = 35, mFPCA.p = 3, mFPCA.m = 2, mFPCA.npc = NULL,
                      LongiModel.method = c('fpca.sc', 'lme'),
                      sFPCA.pve = 0.95, sFPCA.nbasis = 10, sFPCA.npc = NULL,
                      gam.method = 'REML', gam.kT = 10){
  
  
  #######################################
  # Retriving info
  #######################################
  y <- as.matrix(Y) ; colnames(y)<-NULL 
  if( is.null(funcArg) ){
    ss <- seq(0, 1, length.out = ncol(y))
  }else{
    if(length(funcArg)==ncol(y)){
      ss <- funcArg  
    }else{
      warning('length(funcArg) and ncol(Y) do not match. funcArg is re-defined: funcArg = ncol(Y) equally spaced grid points in [0,1].')
      ss <- seq(0, 1, length.out = ncol(y))
    }
  }
  
  Tij <- obsT
  TT <- seq(min(Tij), max(Tij), length.out=numTEvalPoints)
  
  n <- length(unique(subject.index))   # number of subject
  M <- length(ss)   # number of eval.points
  J <- length(TT)
  Ncurves <- nrow(y)  # sum of J_i
  uTij <- unique(Tij)
  
  #######################################
  # bivariate smooth mean function
  #######################################
  fit.fbps <- fbps(data=y, covariates=list(Tij, ss), knots=fbps.knots)
  
  mu.hat <- fit.fbps$Yhat
  
  #######################################
  # 1. marginal FPCA
  #######################################  
  new.y <- y-mu.hat
  
  fit1 <- fpca.face(Y=new.y, pve=mFPCA.pve, knots=mFPCA.knots, p = mFPCA.p, m = mFPCA.m, npc = mFPCA.npc)
  phi.hat <- fit1$efunctions*sqrt(M)   # estimate eigenfunctions
  K.hat <- fit1$npc 
  
  # marginal covariance estimate: sum_k lambda_hat_k * phi_hat_k(s) * phi_hat_k(s')
  lambda.hat <- fit1$evalues/M
  marCovar.hat <- Reduce('+',lapply(seq_len(length(lambda.hat)), function(a) lambda.hat[a]*t(t(phi.hat[,a]))%*%t(phi.hat[,a]) ))
  
  # for scree plot
  fitfit <- fpca.face(Y=new.y, pve=0.9999)
  scree.PVE <- cumsum(fitfit$evalues)/ sum(fitfit$evalues)
  
  # binning
  v <- unlist(lapply(Tij, function(a) which(abs(TT-a) == min(abs(TT-a))))) # Tij close to TT[v]
  
  if(LongiModel.method == 'fpca.sc'){
    
    #######################################
    # 2.1 2nd FPCA for modeling basis coefs
    #######################################   
    
    # create 'sparse functional data' of estimated basis coefs (xi_tilde's) (k=1, 2, ..., npc)
    xi.hat0 <- list()
    for(k in seq_len(K.hat)){
      xihat0.vec <- fit1$scores[,k] / sqrt(M)
      xi.hat0[[k]]<-t(sapply(seq_len(n), function(i) {xi.subj <- matrix(nrow=1, ncol=J) 
      xi.subj[v[which(subject.index==i)]]<-xihat0.vec[which(subject.index==i)]
      return(xi.subj)}))
    }  
    
    # 2nd fpca, and get fitted values for basis coefficients (for all T)
    fit2 <- lapply(xi.hat0, function(a) fpca.sc(Y=a, pve=sFPCA.pve, var=TRUE))
    xi.hat<- lapply(fit2, function(a) a$Yhat)                     # xi_hat(\cdot)
    
    longDynamicsCov.hat.k <- lapply(seq_len(fit1$npc), function(a) Reduce('+', lapply(seq_len(length(fit2[[a]]$evalues)), 
                                                                                      function(b) fit2[[a]]$evalues[b]*t(t(fit2[[a]]$efunctions[,b]))%*%t(fit2[[a]]$efunctions[,b]))))
    
  } else if(LongiModel.method == 'lme'){
    
    #######################################
    # 2.2 linear random effects model (lme)
    #######################################  
    
    lme.coef <- lme.pred <- lme.full.pred <- lme.fit <- list()
    for(k in seq_len(K.hat)){
      lme.dat <- data.frame(Y = fit1$scores[,k]/sqrt(M), X = Tij, subj = subject.index)
      lme.fit0 <- lmer(Y ~ (X|subj), data=lme.dat, REML = TRUE)
      lme.fit[[k]] <- lme.fit0
      lme.coef[[k]] <- coef(lme.fit0)
      lme.pred[[k]] <- fitted(lme.fit0)
      lme.full.pred[[k]] <- t(matrix(predict(lme.fit0, newdata= data.frame(X= rep(TT, n), subj = rep(unique(subject.index), each=length(TT)) )),
                                     length(TT)))
    }
    
    xi.hat <- lme.full.pred
    longDynamicsCov.hat.k <- lapply(seq_len(fit1$npc), function(a) cbind(rep(1, numTEvalPoints), TT) %*% 
                                      as.matrix(as.data.frame( VarCorr(lme.fit[[a]])[[1]]))%*%t(cbind(rep(1, numTEvalPoints), TT))  )
  }
  
  #######################################
  # 3. Fitted values
  #######################################   
  
  fitted <- mu.hat +  t(matrix(rep(fit1$mu, Ncurves), nrow=length(ss))) + 
    do.call(rbind,lapply(seq_len(Ncurves), function(icv) unlist(lapply(xi.hat, function(a) a[subject.index[icv], v[icv]]))))%*% t(phi.hat)
  
  #######################################
  # 4. Yhat for all s and TT
  #######################################
  
  muHat <- matrix(predict.fbps(object=fit.fbps, 
                               newdata=data.frame( x = rep(TT, M),
                                                   z = rep(ss, each=numTEvalPoints) ))$fitted.values, nrow = numTEvalPoints )
  xi.hat.bySubj <- lapply(1:n, function(i) sapply(xi.hat, function(a) a[i,]))
  xi.hat.phi.hat.bySubj <- lapply(xi.hat.bySubj, function(a) t( apply(a, 1, function(b) matrix(b, nrow=1) %*%t(phi.hat) )  ))
  Yhat.all <-  lapply(xi.hat.phi.hat.bySubj , function(a) a + muHat + t(matrix(rep(fit1$mu, J), nrow=length(ss))) )
  
  #######################################
  # 5. Prediction (if specified in arg)
  #######################################
  
  if(!is.null(newdata)){
    
    Jpred <- nrow(newdata)
    
    randDev <- function(row){ 
      i <- newdata$i[row]
      Tpred <- newdata$Ltime[row]
      temp.xi.hat.bySubj <- xi.hat.bySubj[[i]]
      temp.fit <- apply(xi.hat.bySubj[[i]], 2, function(a){ 
        temp.data <- data.frame(y=a, x=TT)
        gam(y~s(x, k=gam.kT, bs='cr'), data=temp.data, method = gam.method)})
      xi.hat.atT <- sapply(temp.fit, function(a) predict.gam(a, newdata=data.frame(x=Tpred)))
      return(xi.hat.atT %*% t(phi.hat))
    }
    
    muHat.pred <- matrix(predict.fbps(object=fit.fbps, 
                                      newdata=data.frame( x = rep(newdata$Ltime, M),
                                                          z = rep(ss, each=Jpred) ))$fitted.values, nrow = Jpred )
    
    predicted <- muHat.pred  + t(matrix(rep(fit1$mu, Jpred), nrow=length(ss))) + do.call(rbind,lapply(seq_len(Jpred), function(icv) randDev(icv)))
    
  }else{
    
    predicted <- NULL
    newdata.hat <- NULL
    
  }
  
  #######################################
  # OUTPUT
  #######################################
  
  if(LongiModel.method == 'fpca.sc'){
    return(list(i = subject.index, # index for subject
                funcArg = ss,   # eval points in s direction
                visitTime = TT, # eval points in T direction
                
                fitted.values = fitted,     # fitted values
                fitted.values.all = Yhat.all, # fitted and predicted values
                predicted.values = predicted, # predicted values for specified subject's future observation
                
                bivariateSmoothMeanFunc = muHat, # estimated mean surface
                
                mFPCA.efunctions = phi.hat, # estimated marginal eigenfunctions
                mFPCA.evalues = lambda.hat, # estimated marginal eigenvalues
                mFPCA.npc = K.hat, # number of marginal PCs selected at pre-specified PVE
                mFPCA.scree.eval = fitfit$evalues / M , # estimated eigenvalues for PVE=0.9999
                
                sFPCA.xiHat.bySubj = xi.hat.bySubj, # estimated basis coef functiosn by subject
                sFPCA.npc = unlist(lapply(fit2, function(a) a$npc)), # number of PCs selected in the second FPCA step for each k
                
                mFPCA.covar = marCovar.hat,   # estimated marginal covariance
                sFPCA.longDynCov.k = longDynamicsCov.hat.k) # estimated covariance of longitudinal dynamics
    )  
  }else if(LongiModel.method == 'lme'){
    
    return(list(i = subject.index, # index for subject
                funcArg = ss,   # eval points in s direction
                visitTime = TT, # eval points in T direction
                
                fitted.values = fitted,     # fitted values
                fitted.values.all = Yhat.all, # fitted and predicted values
                predicted.values = predicted, # predicted values for specified subject's future observation
                
                bivariateSmoothMeanFunc = muHat, # estimated mean surface
                
                mFPCA.efunctions = phi.hat, # estimated marginal eigenfunctions
                mFPCA.evalues = lambda.hat, # estimated marginal eigenvalues
                mFPCA.npc = K.hat, # number of marginal PCs selected at pre-specified PVE
                mFPCA.scree.eval = fitfit$evalues / M , #estimated eigenvalues for PVE=0.9999
                
                sFPCA.xiHat.bySubj = xi.hat.bySubj, # estimated basis coef functiosn by subject
                
                mFPCA.covar = marCovar.hat,   # estimated marginal covariance
                sFPCA.longDynCov.k = longDynamicsCov.hat.k) # estimated covariance of longitudinal dynamics
    )   
    
  }
  
  
}



