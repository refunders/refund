##' Functional principal component analysis by a two-stage method
##'
##' This function performs functional PCA by performing an ordinary singular
##' value decomposition on the functional data matrix, then smoothing the right
##' singular vectors by smoothing splines.
##'
##' The eigenvalues are the same as those obtained from eigendecomposition of
##' the sample covariance matrix. Please note that we expect to merge this
##' function into \code{\link{fpca.ssvd}} in future versions of the package.
##'
##' @param Y data matrix (rows: observations; columns: grid of eval. points)
##' @param npc how many smooth SVs to try to extract. If \code{NA} (the
##' default), the hard thresholding rule of Donoho and Gavish (2013) is used.
##' Application of this rule to functional PCA should be regarded as
##' experimental.
##' @param center center \code{Y} so that its column-means are 0? Defaults to
##' \code{TRUE}
##' @param argvals index vector of \eqn{J} entries for data in X; defaults to a
##' sequence from 0 to 1.
##' @param smooth logical; defaults to TRUE, if NULL, no smoothing of
##' eigenvectors.
##' @return A list with components \item{Yhat}{predicted data matrix}
##' \item{scores}{matrix of scores} \item{mu}{mean function} \item{npc}{number
##' of principal components} \item{eigenvectors}{matrix of eigenvectors}
##' \item{eigenvalues}{vector of eigenvalues}
##' @author Luo Xiao \email{lxiao@@jhsph.edu}
##' @export
##' @importFrom stats smooth.spline
##' @seealso \code{\link{fpca.sc}} and \code{\link{fpca.face}} for FPCA based
##' on smoothing a covariance estimate; \code{\link{fpca.ssvd}} for another
##' SVD-based approach.
##' @references
##'
##' Xiao, L., Ruppert, D., Zipunnikov, V., and Crainiceanu, C., (2013), Fast
##' covariance estimation for high-dimensional functional data. (submitted)
##' \url{http://arxiv.org/abs/1306.5718}.
##'
##' Donoho, D.L. and Gavish, M. (2013). The optimal hard threshold for singular
##' values is 4/sqrt(3).  eprint arXiv:1305.5870. Available from
##' \url{http://arxiv.org/abs/1305.5870}.
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
##'   Phi <- results$eigenvectors
##'   eigenvalues <- results$eigenvalues
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
##'       plot(t[seq],Phi[seq,k]*sqrt(J),type="l",lwd = 3,
##'            ylim = c(-2,2),col = "red",
##'            ylab = paste("Eigenfunction ",k,sep=""),
##'            xlab="t",main="SVDS")
##'
##'       lines(t[seq],phi[seq,k],lwd = 2, col = "black")
##'       }
fpca2s <-
function(Y,npc=NA,center = TRUE, argvals = NULL,smooth=TRUE){

  ## data: X, I by J data matrix
  ## argvals: vector of J
  X <- Y
  data_dim <- dim(X)
  I <- data_dim[1]
  J <- data_dim[2]

  if(is.na(npc)){
    npc <- getNPC.DonohoGavish(X)
  }

  if(is.null(argvals)) argvals <- seq(0, 1, length=J)

  meanX <- rep(0,J)
  if(center) {
    meanX <- apply(X,2,function(x) mean(x,na.rm=TRUE))
    meanX <- smooth.spline(argvals,meanX,all.knots =TRUE)$y
    X <- t(t(X)-meanX)
  }
  ### SVD decomposition
  if(J>I){
    VV <- X%*%t(X)
    Eigen <- eigen(VV)
    D <- Eigen$values
    sel <- (D>0)
    V <- Eigen$vectors[,sel==1]
    D <- D[sel==1]
    D <- sqrt(D)
    U <- t(X)%*%V%*%diag(1/D)
  }

  if(J<=I){
    UU <- t(X)%*%X
    Eigen <- eigen(UU)
    D <- Eigen$values
    U <- Eigen$vectors[,D>0]
    D <- D[D>0]
    D <- sqrt(D)
  }

  lambda <- D^2/(I-1)/J

  if(!is.numeric(npc)) stop("Invalid <npc>.")
  if(npc<1 | npc>min(I,J)) stop("Invalid <npc>.")
  #### end: borrowed from Fabian's code
  message("Extracted ", npc, " smooth components.")

  if(smooth==TRUE){
    #### smoothing
    for(j in 1:npc){
      #temp = pspline(U[,j],argvals,knots=knots,p=p,m=m)$fitted.values
      temp = smooth.spline(argvals,U[,j],all.knots =TRUE)$y
      U[,j] = temp/sqrt(sum(temp^2))
    }
  }
  eigenvectors=U[,1:npc]
  eigenvalues = lambda[1:npc]
  scores = X%*%U[,1:npc]/sqrt(J)

  Yhat = t(eigenvectors%*%t(scores)*sqrt(J) + meanX)
  return(list(Yhat = Yhat,scores = scores, mu = meanX,eigenvectors=eigenvectors,eigenvalues=eigenvalues, npc=npc))
}
