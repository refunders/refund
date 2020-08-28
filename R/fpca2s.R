##' Functional principal component analysis by a two-stage method
##'
##' This function performs functional PCA by performing an ordinary singular
##' value decomposition on the functional data matrix, then smoothing the right
##' singular vectors by smoothing splines.
##'
##' Note that \code{fpca2s} computes smoothed orthonormal eigenvectors
#'  of the supplied function evaluations (and associated scores), not (!)
#'  evaluations of the smoothed orthonormal eigenfunctions. The smoothed
#'  orthonormal eigenvectors are then rescaled by the length of the domain
#'  defined by \code{argvals} to have a quadratic integral approximately equal
#'  to one (instead of crossproduct equal to one), so they approximate the behavior
#'  of smooth eigenfunctions. If \code{argvals} is not equidistant,
#'  \code{fpca2s} will simply return the smoothed eigenvectors without rescaling,
#'  with a warning.
##'
#'@param Y data matrix (rows: observations; columns: grid of eval. points)
#'@param ydata a data frame \code{ydata} representing
#'  irregularly observed functions. NOT IMPLEMENTED for this method.
#'@param argvals the argument values of the function evaluations in \code{Y},
#'  defaults to a equidistant grid from 0 to 1. See Details.
#'@param npc how many smooth SVs to try to extract, if \code{NA} (the default)
#'  the hard thresholding rule of Donoho, Gavish (2013) is used (see Details,
#'  References).
#'@param center center \code{Y} so that its column-means are 0? Defaults to
#'  \code{TRUE}
##' @param smooth logical; defaults to TRUE, if NULL, no smoothing of
##' eigenvectors.
#'@return an \code{fpca} object like that returned from \code{\link{fpca.sc}},
#'  with entries \code{Yhat}, the smoothed trajectories, \code{Y}, the observed
#'  data, \code{scores}, the estimated FPC loadings, \code{mu}, the column means
#'  of \code{Y} (or a vector of zeroes if \code{!center}),  \code{efunctions},
#'  the estimated smooth FPCs (note that these are orthonormal vectors, not
#'  evaluations of orthonormal functions if \code{argvals} is not equidistant),
#'  \code{evalues}, their associated eigenvalues, and \code{npc}, the number of
#'  smooth components that were extracted.
##' @author Luo Xiao \email{lxiao@@jhsph.edu}, Fabian Scheipl
##' @export
##' @importFrom stats smooth.spline
##' @importFrom stats lm.fit
##' @seealso \code{\link{fpca.sc}} and \code{\link{fpca.face}} for FPCA based
##' on smoothing a covariance estimate; \code{\link{fpca.ssvd}} for another
##' SVD-based approach.
##' @references
##'
##' Xiao, L., Ruppert, D., Zipunnikov, V., and Crainiceanu, C., (2013), Fast
##' covariance estimation for high-dimensional functional data. (submitted)
##' \url{https://arxiv.org/abs/1306.5718}.
##'
##' Gavish, M., and Donoho, D. L.  (2014). The optimal hard threshold for
##' singular values is 4/sqrt(3).  \emph{IEEE Transactions on Information Theory}, 60(8), 5040--5053.
##' @examples
##'
##'   #### settings
##'   I <- 50 # number of subjects
##'   J <- 3000 # dimension of the data
##'   t <- (1:J)/J # a regular grid on [0,1]
##'   N <- 4 #number of eigenfunctions
##'   sigma <- 2 ##standard deviation of random noises
##'   lambdaTrue <- c(1,0.5,0.5^2,0.5^3) # True eigenvalues
##'
##'   case = 1
##'   ### True Eigenfunctions
##'
##'   if(case==1) phi <- sqrt(2)*cbind(sin(2*pi*t),cos(2*pi*t),
##'                                    sin(4*pi*t),cos(4*pi*t))
##'   if(case==2) phi <- cbind(rep(1,J),sqrt(3)*(2*t-1),
##'                            sqrt(5)*(6*t^2-6*t+1),
##'                            sqrt(7)*(20*t^3-30*t^2+12*t-1))
##'
##'   ###################################################
##'   ########     Generate Data            #############
##'   ###################################################
##'   xi <- matrix(rnorm(I*N),I,N);
##'   xi <- xi%*%diag(sqrt(lambdaTrue))
##'   X <- xi%*%t(phi); # of size I by J
##'   Y <- X + sigma*matrix(rnorm(I*J),I,J)
##'
##'   results <- fpca2s(Y,npc=4,argvals=t)
##'   ###################################################
##'   ####               SVDS               ########
##'   ###################################################
##'   Phi <- results$efunctions
##'   eigenvalues <- results$evalues
##'
##'   for(k in 1:N){
##'     if(Phi[,k]%*%phi[,k]< 0)
##'       Phi[,k] <- - Phi[,k]
##'   }
##'
##'  ### plot eigenfunctions
##'  par(mfrow=c(N/2,2))
##'  seq <- (1:(J/10))*10
##'  for(k in 1:N){
##'       plot(t[seq],Phi[seq,k]*sqrt(J),type='l',lwd = 3,
##'            ylim = c(-2,2),col = 'red',
##'            ylab = paste('Eigenfunction ',k,sep=''),
##'            xlab='t',main='SVDS')
##'
##'       lines(t[seq],phi[seq,k],lwd = 2, col = 'black')
##'       }
fpca2s <- function(Y = NULL, ydata = NULL, argvals = NULL, npc = NA, center = TRUE,
  smooth = TRUE) {

  ## data: Y, I by J data matrix argvals: vector of J
  stopifnot(!is.null(Y))
  if (any(is.na(Y)))
    stop("No missing values in <Y> allowed.")
  if (!is.null(ydata)) {
    stop(paste("<ydata> argument for irregular data is not supported,", "please use fpca.sc instead."))
  }

  X <- Y
  data_dim <- dim(X)
  I <- data_dim[1]
  J <- data_dim[2]

  if (is.na(npc)) {
    npc <- getNPC.DonohoGavish(X)
  }

  irregular <- FALSE
  if (!is.null(argvals)) {
    stopifnot(is.numeric(argvals), length(argvals) == J, all(!is.na(argvals)))
    if (any(diff(argvals)/mean(diff(argvals)) > 1.05 | diff(argvals)/mean(diff(argvals)) <
      0.95)) {
      warning(paste("non-equidistant <argvals>-grid detected:", "fpca2s will return orthonormal eigenvectors of the function evaluations",
        "not evaluations of the orthonormal eigenvectors.", "Use fpca.sc() for the latter instead."))
      irregular <- TRUE
    }
  } else {
    argvals <- seq(0, 1, length = J)
  }

  meanX <- rep(0, J)
  if (center) {
    meanX <- apply(X, 2, function(x) mean(x, na.rm = TRUE))
    meanX <- smooth.spline(argvals, meanX, all.knots = TRUE)$y
    X <- t(t(X) - meanX)
  }
  ### SVD decomposition
  if (J > I) {
    VV <- X %*% t(X)
    Eigen <- eigen(VV)
    D <- Eigen$values
    sel <- (D > 0)
    V <- Eigen$vectors[, sel == 1]
    D <- D[sel == 1]
    D <- sqrt(D)
    U <- t(X) %*% V %*% diag(1/D)
  }

  if (J <= I) {
    UU <- t(X) %*% X
    Eigen <- eigen(UU)
    D <- Eigen$values
    U <- Eigen$vectors[, D > 0]
    D <- D[D > 0]
    D <- sqrt(D)
  }

  lambda <- D^2/(I - 1)/J

  if (!is.numeric(npc))
    stop("Invalid <npc>.")
  if (npc < 1 | npc > min(I, J))
    stop("Invalid <npc>.")
  #### end: borrowed from Fabian's code
  message("Extracted ", npc, " smooth components.")

  if (smooth == TRUE) {
    #### smoothing
    for (j in 1:npc) {
      temp = smooth.spline(argvals, U[, j], all.knots = TRUE)$y
      U[, j] = temp
    }
  }
  if (!irregular) {
    # scale smooth eigenvectors so they're scaled as realizations of orthonormal
    # eigenfunctions i.e. so that colSums(diff(argvals) * U^2) == 1 instead of
    # crossprod(U) == diag(npc)
    scale <- sqrt(mean(diff(argvals)))
  } else {
    scale <- 1
  }

  eigenvectors = U[, 1:npc, drop = FALSE]/scale
  scores = unname(t(lm.fit(x = eigenvectors, y = t(X))$coefficients))
  eigenvalues = diag(var(scores))


  Yhat = t(eigenvectors %*% t(scores) + meanX)
  ret = list(Yhat = Yhat, Y = Y, scores = scores, mu = meanX, efunctions = eigenvectors,
    evalues = eigenvalues, npc = npc, argvals = argvals)
  class(ret) = "fpca"
  return(ret)
}
