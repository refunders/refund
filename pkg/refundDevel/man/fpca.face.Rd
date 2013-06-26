\name{fpca.face}
\alias{fpca.face}
\title{
Fast Covariance Estimation
}
\description{
A fast implementation of the sandwich smoother (Xiao et al., 2013) for covariance matrix smoothing. Pooled generalized cross validation
at the data level is used for selecting the smoothing parameter. 
}
\usage{
fpca.face(X, center = TRUE, t = NULL, knots = 35, p = 3, m = 2, lambda = NULL, 
     pve = 0.99, score.method = "int",
     search.grid = TRUE, search.length = 100, method = "L-BFGS-B",
     lower = -20, upper = 20, control = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
I by J data matrix; missing data allowed
}
  \item{center}{
logical; defaults to TRUE and data will be centered across subjects  	
  }
  
  \item{t}{
vector of J
}
  \item{knots}{
vector of knots or number of (interior) knots; defaults to 35
}
  \item{p}{
degrees of B-splines; defaults to 3
}
  \item{m}{
order of differencing penalty; defaults to 2
}
  \item{lambda}{
user-specified smoothing parameter; defaults to NULL and uses grid-searching method or \command{optim}
}
  \item{pve}{
percentage of variance explained for selecting the number of principal components; defaults to 0.99  	
}
  \item{score.method}{
method for predicting scores; defaults to "int" for numerical integration,
"blup" for best linear unbiased prediction
}
  \item{search.grid}{
logicial; defaults to TRUE, if FALSE, uses \command{optim}
}
  \item{search.length}{
number of equidistant (log scale) smoothing parameters to search; defaults to 100
}
  \item{method}{
see \command{optim}; defaults to \command{L-BFGS-B}
}
  \item{lower}{
see \command{optim}; defaults to -20 for the log of the smoothing paramter
}
  \item{upper}{
see \command{optim}; defaults to 20 for the log of the smoothing paramter
}
  \item{control}{
see \command{optim}
}
}


\value{
\item{npc}{number of principal components}
\item{eigenvectors}{matrix of eigenvectors}
\item{eigenvalues}{vector of eigenvalues}
\item{scores}{matrix of scores}
\item{lambda}{value of smoothing parameter}
}
\references{
Xiao, L., Li, Y., and Ruppert, D. (2013),
Fast bivariate \emph{P}-splines: the sandwich smoother,
\emph{Journal of the Royal Statistical Society: Series B}, 75(3), 577-599.

Xiao, L., Ruppert, D., Zipunnikov, V., and Crainiceanu, C., (2013),
Fast covariance estimation for high-dimensional functional data, \emph{manuscript}.
}
\author{
Luo Xiao \email{lxiao@jhsph.edu}
}

%\seealso{
%\code{\link{fpca.sc}}
%}
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

results <- fpca.face(Y,center = TRUE, t=t,knots=100,pve=0.99)
###################################################
####               FACE                ########
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
       xlab="t",main="FACE")
  
  lines(t[seq],phi[seq,k],lwd = 2, col = "black")
}

     
}

%\keyword{ ~FACE}
%\keyword{ ~covariance}
%\keyword{ ~smoothing}
