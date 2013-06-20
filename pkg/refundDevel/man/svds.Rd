\name{svds}
\alias{svds}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Singular Value Decomposition plus Smoothing
}
\description{
Smoothing by \emph{P}-splines the eigenvectors obtained from the singular value decomposition of the functional data matrix.
}
\usage{
svds(Y, t = NULL, smooth = TRUE, N = NULL)
}

\arguments{
  \item{Y}{
J by I data matrix
}
  \item{t}{
vector of J
}
  \item{smooth}{
logical; defaults to TRUE, if NULL, no smoothing of eigenvectors
}
  \item{N}{
first N eigenvectors to smooth
}
}

\details{
The eigenvalues are the same as those obtained from eigendecomposition of sample covariance matrix.
}
\value{
\item{N}{number of principal components}
\item{eigenvectors}{matrix of eigenvectors}
\item{eigenvalues}{vector of eigenvalues}
\item{scores}{matrix of scores}
}
\references{

Xiao, L., Ruppert, D., Zipunnikov, V., and Crainiceanu, C., (2013),
Fast covariance estimation for high-dimensional functional data, \emph{manuscript}.
}
\author{
Luo Xiao \email{lxiao@jhsph.edu}
}

\seealso{
\code{\link{face}}
}
\examples{
  #### settings
  I <- 50 # number of subjects
  J <- 3000 # dimension of the data
  t <- (1:J)/J # a regular grid on [0,1]
  N <- 4 #number of eigenfunctions
  sigma <- 2 ##standard deviation of random noises
  lambdaTrue <- c(1,0.5,0.5^2,0.5^3) # True eigenvalues
  
  case = 1
  ### True Eigenfunctions
  
  if(case==1) phi <- sqrt(2)*cbind(sin(2*pi*t),cos(2*pi*t),
                                   sin(4*pi*t),cos(4*pi*t))
  if(case==2) phi <- cbind(rep(1,J),sqrt(3)*(2*t-1),
                           sqrt(5)*(6*t^2-6*t+1),
                           sqrt(7)*(20*t^3-30*t^2+12*t-1))

  ###################################################
  ########     Generate Data            #############
  ###################################################
  xi <- matrix(rnorm(I*N),I,N);
  xi <- xi\%*\%diag(sqrt(lambdaTrue))
  X <- xi\%*\%t(phi); # of size I by J
  Y <- X + sigma*matrix(rnorm(I*J),I,J)
  
  results <- svds(Y,t=t,N=4)
  ###################################################
  ####               SVDS               ########
  ###################################################  
  Phi <- results$eigenvectors
  eigenvalues <- results$eigenvalues
  
  for(k in 1:N){
    if(Phi[,k]\%*\%phi[,k]< 0) 
      Phi[,k] <- - Phi[,k]
  }
 
 ### plot eigenfunctions
 par(mfrow=c(N/2,2))
 seq <- (1:(J/10))*10
 for(k in 1:N){
      plot(t[seq],Phi[seq,k]*sqrt(J),type="l",lwd = 3, 
           ylim = c(-2,2),col = "red",
           ylab = paste("Eigenfunction ",k,sep=""),
           xlab="t",main="SVDS")
      
      lines(t[seq],phi[seq,k],lwd = 2, col = "black")
      }
     
}
