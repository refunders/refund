##' Multilevel functional principal components analysis by smoothed covariance
##'
##' Decomposes functional observations using functional principal components
##' analysis. A mixed model framework is used to estimate scores and obtain
##' variance estimates.
##'
##' This function computes a multilevel FPC decomposition for a set of observed curves,
##' which may be sparsely observed and/or measured with error. A mixed model
##' framework is used to estimate level 1 and level 2  scores.
##'
##' MFPCA was proposed in Di et al. (2009), with variations for 
##' MFPCA with sparse data in Di et al. (2014). 
##' \code{mfpca.sc} uses penalized splines to smooth the covariance functions, as
##' Described in Di et al. (2009) and Goldsmith et al. (2013).
##'
##' @param Y, The user must supply a matrix of functions on a regular grid
##' @param id Must be supplied, a vector containing the id information used to identify clusters
##' @param visit A vector containing information used to identify visits. Defaults to \code{NULL}.
##' @param twoway logical, indicating whether to carry out twoway ANOVA and calculate visit-specific means. Defaults to \code{FALSE}.
##' @param argvals function argument.
##' @param nbasis number of B-spline basis functions used for estimation of the
##' mean function and bivariate smoothing of the covariance surface.
##' @param pve proportion of variance explained: used to choose the number of
##' principal components.
##' @param npc prespecified value for the number of principal components (if
##' given, this overrides \code{pve}).
##'
##' @param makePD logical: should positive definiteness be enforced for the
##' covariance surface estimate? Defaults to \code{FALSE} Only \code{FALSE} is currently supported.
##' @param center logical: should an estimated mean function be subtracted from
##' \code{Y}? Set to \code{FALSE} if you have already demeaned the data using
##' your favorite mean function estimate.
##'
##'  @param cov.est.method covariance estimation method. If set to \code{1}, a
##' one-step method that applies a bivariate smooth to the \eqn{y(s_1)y(s_2)}
##' values. This can be very slow. If set to \code{2} (the default), a two-step
##' method that obtains a naive covariance estimate which is then smoothed. \code{2} is currently supported.
##' @param integration quadrature method for numerical integration; only
##' \code{"trapezoidal"} is currently supported.
##' 
##' @return An object of class \code{mfpca} containing:
##' \item{Yhat}{FPC approximation (projection onto leading components)
##' of \code{Y}, estimated curves for all subjects and visits}
##' \item{Yhat.subject}{estimated subject specific curves for all subjects} 
##' \item{Y}{the observed data}\item{scores}{\eqn{n
##' \times npc} matrix of estimated FPC scores for level1 and level2.} \item{mu}{estimated mean
##' function (or a vector of zeroes if \code{center==FALSE}).} \item{efunctions
##' }{\eqn{d \times npc} matrix of estimated eigenfunctions of the functional
##' covariance, i.e., the FPC basis functions for levels 1 and 2.} \item{evalues}{estimated
##' eigenvalues of the covariance operator, i.e., variances of FPC scores for levels 1 and 2.}
##' \item{npc }{number of FPCs: either the supplied \code{npc}, or the minimum
##' number of basis functions needed to explain proportion \code{pve} of the
##' variance in the observed curves for levels 1 and 2.} \item{sigma2}{estimated measurement error
##' variance.} \item{eta}{the estimated visit specific shifts from overall mean.}
##' @author Julia Wrobel \email{jw3134@@cumc.columbia.edu}, Jeff Goldsmith \email{jeff.goldsmith@@columbia.edu}, and Chongzhi Di
##' @references Di, C., Crainiceanu, C., Caffo, B., and Punjabi, N. (2009).
##' Multilevel functional principal component analysis. \emph{Annals of Applied
##' Statistics}, 3, 458--488.
##' 
##' Di, C., Crainiceanu, C., Caffo, B., and Punjabi, N. (2014).
##' Multilevel sparse functional principal component analysis. \emph{Stat}, 3, 126--143.
##'
##' Goldsmith, J., Greven, S., and Crainiceanu, C. (2013). Corrected confidence
##' bands for functional data using principal components. \emph{Biometrics},
##' 69(1), 41--51.
##'
##' @export
##' @importFrom Matrix nearPD Matrix t as.matrix
##' @importFrom mgcv gam predict.gam
##' @importFrom gamm4 gamm4
##' 
##' @examples 
##'  \dontrun{
##'  data(DTI)
##'  DTI = subset(DTI, Nscans < 6)  ## example where all subjects have 6 or fewer visits
##'  id  = DTI$ID
##'  Y = DTI$cca
##'  mfpca.DTI =  mfpca.sc(Y=Y, id = id, twoway = TRUE)
##'  }
## npc=1 seems to give error


###############################################################################################
#
# mfpca.sc : multilevel functional principal components analysis by smoothed covariance
#
#################################################################################################
mfpca.sc <- function(Y = NULL, id=NULL, visit=NULL, twoway = FALSE,
                     argvals = NULL, nbasis = 10, pve = 0.99, npc = NULL,
                     #####################################################################
                     ## the following arguments from fpca.sc have not yet been integrated. should they be?
                     makePD = FALSE, center = TRUE, cov.est.method = 2, integration = "trapezoidal") {
  
  stopifnot((!is.null(Y) && !is.null(id)))
  
  ## Y.df keeps visit and id values so that we can sort more easily later on
  # make sure to separate out Y from the data so you can use it for what you need later
  
  if (!is.null(visit)){
    visit = as.integer(factor(visit))
  }else{ visit = ave(id, id, FUN=seq_along)}
  
  Y.df <- data.frame(id=id, visit = visit)
  Y.df$Y <- Y
  J = length(unique(visit))  ## gets max number of visits
  M = length(unique(id)) ## number of subjects
  D = NCOL(Y)  
  I = NROW(Y)  
  nVisits <- data.frame(table(id))  ## calculate number of visitis for each subject
  colnames(nVisits) <- c("id", "numVisits")
  
  if (is.null(argvals))  argvals = seq(0, 1, , D)  
  d.vec = rep(argvals, each = I) 
  
  if (center) {
    gam0 = gam(as.vector(Y) ~ s(d.vec, k = nbasis))
    mu = predict(gam0, newdata = data.frame(d.vec = argvals))  
  } else{
    Y.tilde = Y
    mu = rep(0, D)
  }
  
  mueta = matrix(0, D, J) ## matrix to store visit-specific means
  eta = matrix(0, D, J)
  
  if (center && twoway) {
    Y.split = as.list(rep(NA, J))   
    for(j in 1:J) {
      Y.split[[j]] = subset(Y.df, visit == j)
      d.vec.split = rep(argvals, each = NROW(Y.split[[j]]))
      fit.eta = gam(as.vector(Y.split[[j]]$Y) ~ s(d.vec.split, k = nbasis))
      mueta[, j] = predict(fit.eta, newdata = data.frame(d.vec.split = argvals)) 
      eta[,j] <- mueta[,j] - mu
      
      Y.split[[j]]$Y.tilde = Y.split[[j]]$Y - matrix(mueta[,j], NROW(Y.split[[j]]), D, byrow = TRUE)      
    }
    ## merge split datasets to get Y.tilde matrix
    Y.df.new = Reduce(function(...) merge(..., by = c("id", "visit",  "Y", "Y.tilde"), all = TRUE, sort = FALSE), Y.split)
    Y.df.new = Y.df.new[order(Y.df.new$id, Y.df.new$visit),]
    Y.tilde = Y.df.new$Y.tilde
  } else{
    Y.tilde = Y - matrix(mu, I, D, byrow = TRUE)  ## subtract the mean function from Y
  }
  

  ####################################################################################
  ####################################################################################################
  # cov.est.method == 2 obtains a naive covariance estimate which is then smoothed
  ####################################################################################################
  
  ##################################################################################
  #      CALCULATE Kt, THE TOTAL COVARIANCE
  # 
  ##################################################################################
  
  if (cov.est.method == 2) {  ## make sure to exclude this if you don't incorporate cov.est.method == 1
    cov.sum = cov.count = cov.mean = matrix(0, D, D)
    row.ind = 0  ## get row index of data frame 
    for (m in 1:M) {
      for(j in 1:nVisits[m, 2]) {
        row.ind.temp <- row.ind  + j
        obs.points = which(!is.na(Y[row.ind.temp,]))
        cov.count[ obs.points, obs.points ] <- cov.count[ obs.points, obs.points ] + 1
        cov.sum[ obs.points, obs.points ] <- cov.sum[ obs.points, obs.points ] + tcrossprod(Y.tilde[row.ind.temp, obs.points])
      }
      row.ind = row.ind + nVisits[m, 2]
    } 
    
    
    G.0 = ifelse(cov.count == 0, NA, cov.sum/cov.count)  
    diag.G0 = diag(G.0) 
    diag(G.0) = NA
    row.vec = rep(argvals, each = D)
    col.vec = rep(argvals, D)
    s.npc.0 =predict(gam(as.vector(G.0) ~ te(row.vec, col.vec, k = nbasis), weights = as.vector(cov.count)), 
                     newdata = data.frame(row.vec = row.vec, col.vec = col.vec))
    npc.0 = matrix(s.npc.0, D, D)
    npc.0 = (npc.0 + t(npc.0))/2  ## the smoothed (total) covariance matrix
    
  } ## end of if(cov.est.method == 2)
  
  ###
  ###
  ###
  ###
  #  next add option cov.est.method == 1 and makePD argument
  ###
  ###
  ###
  ###
  
  
  ##################################################################################
  #      CALCULATE Kb, THE BETWEEN COVARIANCE
  #
  ###     Calculate the pairs to calculate the between covariance function
  ###     first, calculate the possible pairs of same-subject visits
  ##################################################################################
  
  if (cov.est.method == 2) {  ## make sure to exclude this if you don't incorporate cov.est.method == 1
    cov.sum = cov.count = cov.mean = matrix(0, D, D)
    row.ind = 0
    ids.KB = nVisits[nVisits$numVisits > 1, c("id")]
    
    for(m in 1:M) {
      if (Y.df$id[m] %in% ids.KB){ ## check if mth subject has at least 2 visits
        for(j in 1:(nVisits[m, 2]-1) ) 
          #row.ind1 <- (m-1)*J + j 
          row.ind1 <- row.ind + j
        obs.points1 = which(!is.na(Y[row.ind1,]))
        for(k in (j+1):nVisits[m, 2] ) {
          #row.ind2 <- (m-1)*J + k
          row.ind2 <- row.ind + k
          obs.points2 = which(!is.na(Y[row.ind2,]))
          cov.count[ obs.points1, obs.points2 ] <- cov.count[ obs.points1, obs.points2 ] + 1
          cov.sum[ obs.points1, obs.points2 ] <- cov.sum[ obs.points1, obs.points2 ] + tcrossprod(Y.tilde[row.ind1, obs.points1], Y.tilde[row.ind2, obs.points2]) 
          cov.count[ obs.points2, obs.points1 ] <- cov.count[ obs.points2, obs.points1 ] + 1
          cov.sum[ obs.points2, obs.points1 ] <- cov.sum[ obs.points2, obs.points1 ] +  tcrossprod(Y.tilde[row.ind2, obs.points2], Y.tilde[row.ind1, obs.points1])                                                        
        }
      }
      row.ind = row.ind + nVisits[m, 2]
    }
    
    G.0b <- ifelse(cov.count==0, NA,  cov.sum/cov.count)  ## between covariance
    
    row.vec = rep(argvals, each = D)  
    col.vec = rep(argvals, D)        
    
    s.npc.0b =predict(gam(as.vector(G.0b) ~ te(row.vec, col.vec, k = nbasis), weights = as.vector(cov.count)), 
                      newdata = data.frame(row.vec = row.vec, col.vec = col.vec))
    npc.0b = matrix(s.npc.0b, D, D)
    npc.0b = (npc.0b + t(npc.0b))/2  ##  smoothed (between) covariance matrix
    
  }
  
  
  ###
  #
  #  add option for covariance method == 1 and makePD arguments
  #
  ###
  
  ###########################################################################################
  #      CALCULATE Kw, THE WITHIN COVARIANCE
  #
  ### use smoothed covariance total and covariance between to calculated covariance within
  ############################################################################################
  
  s.npc.0w = s.npc.0 - s.npc.0b
  npc.0w = matrix(s.npc.0w, D, D)
  npc.0w = (npc.0w + t(npc.0w))/2  ## smoothed within covariance matrix
  
  ###################################################################################
  ###     Estimate eigen values and eigen functions at two levels by calling the 
  ###     eigen function (in R "base" package) on discretized covariance matrices.
  
  w <- quadWeights(argvals, method = integration)
  Wsqrt <- diag(sqrt(w))
  Winvsqrt <- diag(1/(sqrt(w)))
  
  npc.0wb = list(level1 = npc.0b, level2 = npc.0w)   
  V <- lapply(npc.0wb, function(x) Wsqrt %*% x %*% Wsqrt)
  evalues = lapply(V, function(x) eigen(x, symmetric = TRUE, only.values = TRUE)$values)
  evalues = lapply(evalues, function(x) replace(x, which(x <= 0), 0))
  npc = lapply(evalues, function(x) ifelse(is.null(npc), min(which(cumsum(x)/sum(x) >  pve)), npc ))
  efunctions = lapply(names(V), function(x) matrix(Winvsqrt %*% eigen(V[[x]], symmetric = TRUE)$vectors[, seq(len = npc[[x]])], nrow = D, ncol = npc[[x]] ))
  evalues = lapply(names(V), function(x) eigen(V[[x]], symmetric = TRUE, only.values = TRUE)$values[1:npc[[x]]])
  
  names(efunctions) = names(evalues) = names(npc) = c("level1", "level2")
  
  ##################################################################################
  ###
  ###     Estimate principal component scores and predict various functions
  ###
  ##################################################################################
  
  ###################################################################
  # calculate the measurement error variance
  ###################################################################
  
  cov.hat = lapply(c("level1", "level2"), function(x) efunctions[[x]] %*% diag(evalues[[x]]) %*% t(efunctions[[x]]))
  
  T.len <- argvals[D] - argvals[1]
  T1.min <- min(which(argvals >= argvals[1] + 0.25 * T.len))
  T1.max <- max(which(argvals <= argvals[D] - 0.25 * T.len))
  
  DIAG = (diag.G0 - diag(cov.hat[[1]])-diag(cov.hat[[2]]))[T1.min:T1.max]
  w2 <- quadWeights(argvals[T1.min:T1.max], method = integration)
  sigma2 <- max(weighted.mean(DIAG, w = w2, na.rm = TRUE), 0) ### estimated measurement error variance
  
  ################################################################################
  Yhat = Yhat.subject =  matrix(0, I, D)
  
  score1 <- matrix(0, nrow=M, ncol=npc[[1]])    
  score2 <- matrix(0, nrow=I, ncol=npc[[2]])    
  Z1 = efunctions[[1]]
  Z2 = efunctions[[2]]
    
  row.ind = 0
  for(m in 1:M) {
    
    Jm = nVisits[m, 2]  ## number of visits for mth subject
    obs.points = lapply(1:Jm, function(j) which(!is.na(Y[row.ind + j, ])) )
    
    ## Ti is number of non-NA observations across all visits for a subject
    numObs = lapply(obs.points, length)
    Ti = Reduce("+", numObs)
    
    ## clean this up. the idea is to be able to index along observed values for each visit j
    Ti.index = data.frame(Ti = 1:Ti, Jm = unlist(lapply(1:Jm, function(j) rep(j, numObs[[j]]))))
     
    
    cov.y = matrix(0, Ti, Ti)
    Ai = matrix(0, npc[[1]], Ti)
    Bi = matrix(0, npc[[2]] * Jm, Ti)
    
    Yij.center = lapply(1:Jm, function(j) {
      observed = obs.points[[j]]  
      Y[row.ind + j, observed]-mu[observed]-eta[observed, j]
    })
    
    Yi.center = matrix(unlist(Yij.center))
    
    for (j1 in 1:Jm){     
      ## sets diagonal of cov.y matrix
      indices1 = Ti.index[Ti.index$Jm == j1, c("Ti")]
      cov.y [indices1 , indices1] <- Z1[obs.points[[j1]],] %*% diag( evalues[[1]]) %*%  t( Z1[obs.points[[j1]], ] ) +  Z2[obs.points[[j1]], ] %*% diag( evalues[[2]]) %*%  t( Z2[obs.points[[j1]], ] ) + diag( rep(sigma2, numObs[[j1]]) )
      
      if(j1 < Jm && nVisits$numVisits[m] > 1){
        for (j2 in (j1+1):Jm){
          #indices2 = 1:numObs[[j2]] + (j2-1)*numObs[[j2-1]]
          indices2 = Ti.index[Ti.index$Jm == j2, c("Ti")] 
          cov.y[indices2, indices1] <- Z1[obs.points[[j2]], ] %*% diag( evalues[[1]]) %*%  t( Z1[obs.points[[j1]], ] )
          cov.y[indices1, indices2] <- t(cov.y[indices2, indices1])
        }
      }      
    }
    
    ##the following defines Ai and Bi
    for(j in 1:Jm) { 
      
      Ai[1:npc[[1]], Ti.index[Ti.index$Jm == j, c("Ti")]] <- diag(evalues[[1]]) %*% t(Z1[obs.points[[j]],])
      Bi[1:npc[[2]] + npc[[2]] * (j-1), Ti.index[Ti.index$Jm == j, c("Ti")]] <- diag(evalues[[2]]) %*% t(Z2[obs.points[[j]],])
    }
        
    subj.indices = row.ind + 1:Jm
    score1[m ,] <- Ai %*% ginv(cov.y) %*% Yi.center
    score2[subj.indices,] <- t( matrix( Bi %*% ginv(cov.y) %*% Yi.center, ncol=Jm) )    
    
    for(j in 1:Jm) {
      Yhat[row.ind + j,] <- as.matrix(mu) + eta[,j] + Z1 %*% score1[m,] + Z2 %*% score2[row.ind + j ,] 
      Yhat.subject[row.ind + j, ] <- as.matrix(mu) + eta[,j] + as.vector(Z1 %*% score1[m,])
    }
    row.ind = row.ind + nVisits[m, 2]
  }
  
  scores <- list(level1 = score1, level2 = score2)
  
  ############################################################################
  # Returns an object of class 'mfpca'
  ##############################################################################
  
  ret.objects = c("Yhat", "Yhat.subject","Y.df", "scores", "mu", "efunctions", "evalues", "npc", "sigma2", "eta")
  
  ret = lapply(1:length(ret.objects), function(u) get(ret.objects[u]))
  names(ret) = ret.objects
  class(ret) = "mfpca"
  return(ret)
  
} ## end of mfpca function