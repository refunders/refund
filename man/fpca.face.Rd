\name{fpca.face}
\alias{fpca.face}
\title{
Functional principal component analysis with fast covariance estimation
}
\description{
A fast implementation of the sandwich smoother (Xiao et al., 2013) for covariance matrix smoothing. Pooled generalized cross validation
at the data level is used for selecting the smoothing parameter. 
}
\usage{
fpca.face(Y=NULL, ydata=NULL, Y.pred=NULL, argvals = NULL,pve = 0.99,npc=NULL,
          var=FALSE,simul = FALSE, sim.alpha = 0.95,
          center=TRUE,knots = 35, p = 3, m = 2, lambda = NULL, alpha = 1,
          search.grid = TRUE,search.length = 100,
          method = "L-BFGS-B", lower = -20,upper = 20, control = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{Y, ydata}{data matrix (currently only support \code{Y}; rows: observations; columns: grid
  of eval. points); missing data allowed but caution is needed for very sparse data}
  \item{Y.pred }{if desired, a matrix of observed functions to be approximated using the FPC decomposition.}
  \item{argvals}{function argument.}
  \item{pve}{percentage of variance explained for selecting the number of principal components; defaults to 0.99.}
  \item{npc}{number of extracted principal components; defaults to NULL and use \code{pve}.}
  \item{var}{\code{TRUE} or \code{FALSE} indicating whether model-based estimates for the variance of FPCA expansions should be computed.}
\item{simul}{logical: should critical values be estimated for simultaneous confidence intervals?}
\item{sim.alpha}{1 - coverage probability of the simultaneous intervals.}
\item{center}{center \code{Y} so that its column-means are 0? Defaults to \code{TRUE}}
  \item{knots}{vector of knots or number of (interior) knots; defaults to 35.}
  \item{p}{degrees of B-splines; defaults to 3.}
  \item{m}{order of differencing penalty; defaults to 2.}
  \item{lambda}{user-specified smoothing parameter; defaults to NULL and uses grid-searching method or \command{optim}.}

  \item{alpha}{
  the coefficient of the penalty for degrees of freedom in the GCV criterion; defaults to 1
  }
  \item{search.grid}{
logical; defaults to TRUE, if FALSE, uses \command{optima}.
}
  \item{search.length}{
number of equidistant (log scale) smoothing parameters to search; defaults to 100.
}
  \item{method}{
see \command{optim}; defaults to \command{L-BFGS-B}.
}
  \item{lower, upper}{
bounds for log smoothing parameter, passed to \command{optim}; defaults are -20 and 20.
}
  \item{control}{
passed to \command{optim}.
}
}


\value{A list with components
\item{Yhat}{If Y.pred is specified, the smooth version of Y.pred.   Otherwise, if Y.pred=NULL, the smooth version of Y.}
\item{scores}{ \eqn{ n \itmes npc} matrix of scores.}
\item{mu}{estimated mean function (or a vector of zeroes if \code{center==FALSE}).}
\item{efunctions}{matrix of estimated eigenvectors.}
\item{evalues}{vector of estimated eigenvalues.}
\item{npc}{number of principal components.}
\item{sigma2}{estimated measurement error variance.}
\item{diag.var}{diagonal elements of the covariance matrices for each estimated curve.}
\item{VarMats}{a list containing the estimated covariance matrices for each curve in \code{Y}.}
\item{crit.val}{estimated critical values for constructing simultaneous confidence intervals.}
}

\note{We expect to merge this function into \code{\link{fpca.sc}} in future releases of the package; see the details of \code{\link{fpca.sc}} for explanation of some of the inputs and outpus.}

\references{
Xiao, L., Li, Y., and Ruppert, D. (2013).
Fast bivariate \emph{P}-splines: the sandwich smoother,
\emph{Journal of the Royal Statistical Society: Series B}, 75(3), 577-599.

Xiao, L., Ruppert, D., Zipunnikov, V., and Crainiceanu, C., (2013).
Fast covariance estimation for high-dimensional functional data. Submitted). To appear in \emph{Statistics and Computing}. Available at \url{http://arxiv.org/abs/1306.5718}.

Goldsmith, J., Greven, S., and Crainiceanu, C. (2013). Corrected confidence bands for functional data using principal components. \emph{Biometrics}, 69(1), 41--51.
}
\author{
Luo Xiao \email{lxiao@jhsph.edu}
}

\seealso{
  \code{\link{fpca.sc}}  for
  another covariance-estimate based smoothing of \code{Y};
\code{\link{fpca2s}} and \code{\link{fpca.ssvd}} for two SVD-based smoothings.
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

results <- fpca.face(Y,center = TRUE, argvals=t,knots=100,pve=0.99)
###################################################
####               FACE                ########
###################################################  
Phi <- results$efunctions
eigenvalues <- results$evalues

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


