##' Functional Models with Structured Penalties
##'
##' Implements functional model with structured penalties (Randloph et al.,
##' 2012) with scalar outcome and single functional predictor through mixed
##' model equivalence.
##'
##' If there are any missing or infinite values in \code{Y}, and \code{funcs},
##' the corresponding row (or observation) will be dropped. Neither \code{Q}
##' nor \code{L} may contain missing or infinite values.
##'
##' \code{peer_old()} fits the following model:
##'
##' \eqn{y_i=\int {W_i(s)\gamma(s) ds} + \epsilon_i}
##'
##' where \eqn{\epsilon_i ~ N(0,\sigma^2)}.  For all the observations,
##' predictor function \eqn{W_i(s)} is evaluated at K sampling points. Here,
##' \eqn{\gamma (s)} denotes the regression function.
##'
##' Values of \eqn{y_i} and \eqn{W_i(s)}are passed through argument Y and
##' funcs, respectively. Number of elements or rows in \code{Y} and
##' \code{funcs} need to be equal.
##'
##' The estimate of regression functions \eqn{\gamma(s)} is obtained as
##' penalized estimated. Following 3 types of penalties can be used:
##'
##' i.  Ridge: \eqn{I_K}
##'
##' ii.  Second-order difference: [\eqn{d_{i,j}}] with \eqn{d_{i,i} = d_{i,i+2}
##' = 1, d_{i,i+1} = -2}, otherwise \eqn{d_{i,i} =0}
##'
##' iii. Decomposition based penalty: \eqn{bP_Q+a(I-P_Q)} where \eqn{P_Q=
##' Q^T(QQ^T)^{-1}Q}
##'
##' For Decomposition based penalty user need to specify
##' \code{pentype='DECOMP'} and associated Q matrix need to be passed through
##' \code{Q} argument.
##'
##' Alternatively, user can pass directly penalty matrix through argument L.
##' For this user need to specify \code{pentype='USER'} and associated L matrix
##' need to be passed through \code{L} argument.
##'
##' Default penalty is Ridge penalty and user needs to specify \code{RIDGE}.
##' For second-order difference penalty, user needs to specify \code{D2}.
##'
##' @param Y vector of all outcomes
##' @param funcs matrix containing observed functional predictors as rows. Rows
##' with \code{NA} and \code{Inf} values will be deleted.
##' @param argvals matrix (or vector) of indices of evaluations of \eqn{X_i(t)}; i.e. a matrix with
##' \emph{i}th row \eqn{(t_{i1},.,t_{iJ})}
##' @param pentype type of penalty. It can be either decomposition based
##' penalty (\code{DECOMP}) or ridge (\code{RIDGE}) or second-order difference
##' penalty (\code{D2}) or any user defined penalty (\code{USER}). For
##' decomposition based penalty user need to specify Q matrix in Q argument
##' (see details). For user defined penalty user need to specify L matrix in L
##' argument (see details). For Ridge and second-order difference penalty,
##' specification for arguments L and Q will be ignored. Default is
##' \code{RIDGE}.
##' @param L.user penalty matrix. Need to be specified with
##' \code{pentype='USER'}. Number of columns need to be equal with number of
##' columns of matrix specified to \code{funcs}. Each row represents a
##' constraint on functional predictor. This argument will be ignored when
##' value of \code{pentype} is other than \code{USER}.
##' @param Q Q matrix to derive decomposition based penalty. Need to be
##' specified with \code{pentype='DECOMP'}. Number of columns need to be equal
##' with number of columns of matrix specified to \code{funcs}. Each row
##' represents a basis function where functional predictor is expected lie
##' according to prior belief. This argument will be ignored when value of
##' \code{pentype} is other than \code{DECOMP}.
##' @param phia Scalar value of a in decomposition based penalty. Need to be
##' specified with \code{pentype='DECOMP'}.
##' @param se logical; calculate standard error when \code{TRUE}.
##' @param ... additional arguments passed to the \code{\link{lme}} function.
##' @return a list containing: \item{fit}{result of the call to \code{lme}}
##' \item{fitted.vals }{predicted outcomes} \item{Gamma}{estimates with
##' standard error for regression function} \item{GammaHat}{estimates of
##' regression function} \item{se.Gamma}{standard error associated with
##' \code{GammaHat}} \item{AIC }{AIC value of fit (smaller is better)}
##' \item{BIC }{BIC value of fit (smaller is better)}
##' \item{logLik}{(restricted) log-likelihood at convergence}
##' \item{lambda}{estimates of smoothing parameter} \item{N}{number of
##' subjects} \item{K}{number of Sampling points in functional predictor}
##' \item{sigma}{estimated within-group error standard deviation.}
##' @author Madan Gopal Kundu \email{mgkundu@@iupui.edu}
##' @seealso \code{lpeer}, \code{plot.peer}
##' @references Kundu, M. G., Harezlak, J., and Randolph, T. W. (2012).
##' Longitudinal functional models with structured penalties (arXiv:1211.4763
##' [stat.AP]).
##'
##' Randolph, T. W., Harezlak, J, and Feng, Z. (2012). Structured penalties for
##' functional linear models - partially empirical eigenvectors for regression.
##' \emph{Electronic Journal of Statistics}, 6, 323--353.
##' @examples
##'
##' \dontrun{
##' #------------------------------------------------------------------------
##' # Example 1: Estimation with D2 penalty
##' #------------------------------------------------------------------------
##'
##' ## Load Data
##' data(DTI)
##'
##' ## Extract values for arguments for peer() from given data
##' cca = DTI$cca[which(DTI$case == 1),]
##' DTI = DTI[which(DTI$case == 1),]
##'
##' ##1.1 Fit the model
##' fit.cca.peer1 = peer(Y=DTI$pasat, funcs = cca, pentype='D2', se=TRUE)
##' plot(fit.cca.peer1)
##'
##' #------------------------------------------------------------------------
##' # Example 2: Estimation with structured penalty (need structural
##' #            information about regression function or predictor function)
##' #------------------------------------------------------------------------
##'
##' ## Load Data
##' data(PEER.Sim)
##'
##' ## Extract values for arguments for peer() from given data
##' PEER.Sim1<- subset(PEER.Sim, t==0)
##' K<- 100
##' W<- PEER.Sim1[,c(3:(K+2))]
##' Y<- PEER.Sim1[,K+3]
##'
##' ##Load Q matrix containing structural information
##' data(Q)
##'
##' ##2.1 Fit the model
##' Fit1<- peer(Y=Y, funcs=W, pentype='Decomp', Q=Q, se=TRUE)
##' plot(Fit1)
##' }
##'
##' @export
##' @importFrom nlme VarCorr lme pdIdent
peer_old<- function(Y, funcs, argvals=NULL, pentype='Ridge', L.user=NULL, Q=NULL, 
                phia=10^3, se=FALSE, ...)
{
  if(!is.null(argvals))
    stop("argvals is not supported in the current version of refund.")
  
  #Determining K, converting W and Y to matrix
  W<- as.matrix(funcs)
  K<- ncol(W)
  Y<- as.matrix(Y)
  
  #Check 1:Making sure Y has only 1 column
  if(dim(Y)[2]>1) stop("No. of column for Y cannot be greater than 1. \nThe peer() will not proceed further.")
  
  #Check 2: Check the dimension of Y, id, t, W and X
  Yl<- dim(Y)[1]
  chk.eq<- ifelse(Yl==nrow(W), 0 ,1)
  if(chk.eq==1) stop("Length of Y and number of rows of funcs are not equal.\n The peer() will not proceed further.\n")
  
  #Removal of missing and infinite observations
  tdata<- data.frame(Y, W)
  tdata<- tdata[which(apply(is.na(tdata), 1, sum)==0),]
  tdata<- tdata[which(apply(is.infinite(as.matrix(tdata)), 1, sum)==0),]
  tdata<- tdata[!is.na(tdata$Y),]
  Y<- tdata$Y
  W.e<- dim(W)[2]+1; W<- as.matrix(tdata[,c(2:W.e)])
  
  #Determine N
  N<- length(Y)
  
  #Check 3: Checking entry for pentype
  pentypecheck<- toupper(pentype) %in% c('DECOMP', 'DECOMPOSITION', 'RIDGE', 'D2', 'USER')
  if(!pentypecheck)  stop("Error: Specify valid object for argument PENTYPE.")
  
  #Check 4: Some checking/processing for decomposition type of penalty
  if(toupper(pentype)=='DECOMP' | toupper(pentype)=='DECOMPOSITION'){
    
    #4.1: Removing rows containing missing and infinite values
    Q<- Q[which(apply(is.na(Q), 1, sum)==0),]
    Q<- Q[which(apply(is.infinite(Q), 1, sum)==0),]
    Q<- matrix(Q, ncol=K)
    
    #4.2: Compatibility of Q and W matrix
    if(ncol(Q)!=ncol(W)) stop('number of columns of Q need to be equal with number of columns of funcs.\nThe peer() will not proceed further.')
    
    #4.3: Singularity of Q matrix
    Q.eig<- abs(eigen(Q %*% t(Q))$values)
    if(any(Q.eig<1e-12)) stop('Q matrix is singular or near singular.\nThe peer() will not proceed further.\n')
    
    #4.4: Checking for phia
    if(!exists("phia")) stop("Error: Specify valid object for argument PHIA.")
    if(!is.numeric(phia)|is.matrix(phia)|is.matrix(phia)) stop("Specify valid object for argument PHIA.")
  }
  
  #Check 5: Some checking/processing for user type of penalty
  if(toupper(pentype)=='USER'){
    L<- L.user
    
    #Check 5.1: Removing rows containing missing and infinite values
    L<- L[which(apply(is.na(L), 1, sum)==0),]
    L<- L[which(apply(is.infinite(L), 1, sum)==0),]
    L<- matrix(L, ncol=K)
    
    #Check 5.2: Dimension of L matrix
    if(ncol(L)!=ncol(W)) stop('number of columns of L need to be equal with number of columns of funcs.\nThe peer() will not proceed further.')
    
    #Check 5.3: Singularity of L'L matrix
    LL<- t(L)%*%L
    LL.eig<- abs(eigen(LL %*% t(LL))$values)
    if(any(LL.eig<1e-12)) stop("Error: L'L matrix is singular or near singular.\nThe peer() will not proceed further.")
  }
  
  #Generate L matrix for D2 penalty
  if(toupper(pentype)=='D2'){
    Left<- cbind(diag(rep(1,K-2)),rep(0,K-2),rep(0,K-2))
    Middle<- cbind(rep(0,K-2),diag(rep(-2,K-2)),rep(0,K-2))
    Right<- cbind(rep(0,K-2),rep(0,K-2),diag(rep(1,K-2)))
    D.2<- rbind(Left+Middle+Right, c(rep(0, K-2), 1, -2), c(rep(0, K-2), 0, 1))
  }
  
  #Generate W* matrix
  if(toupper(pentype)=='DECOMP' | toupper(pentype)=='DECOMPOSITION'){
    P_Q <- t(Q) %*% solve(Q %*% t(Q)) %*% Q
    L_PEER<- phia*(diag(K)- P_Q) + 1*P_Q
  } else
    if(toupper(pentype)=='RIDGE'){
      L_PEER<- diag(K)
    } else
      if(toupper(pentype)=='D2'){
        L_PEER<- D.2
      } else
        if(toupper(pentype)=='USER'){
          L_PEER<- L
        }
  
  v<- diag(K)
  if(K>N) v<-  svd((data.matrix(W))%*% solve(L_PEER))$v
  W1_PEER<- (data.matrix(W))%*% solve(L_PEER) %*% v
  
  #Fitting the model
  id.bd1<- factor(rep(1, length(Y)))
  out_PEER<- nlme::lme(fixed=Y~1,
                       random=list(id.bd1=nlme::pdIdent(~W1_PEER-1)),
                       ...
  )
  message('The fit is successful.\n')
  
  #Extracting estimates
  Gamma.PEER.hat<-matrix(out_PEER$coeff$random$id.bd1, ncol=1)
  GammaHat <- solve(L_PEER) %*% v %*%Gamma.PEER.hat
  colnames(GammaHat)<- c('Gamma')
  fitted.vals<- summary(out_PEER)$fitted[,2]
  
  #Extracting model diagnostics
  logLik<- summary(out_PEER)$logLik
  AIC<- summary(out_PEER)$AIC
  BIC<- summary(out_PEER)$BIC
  
  #Extracting lambda and variances
  tVarCorr<- nlme::VarCorr(out_PEER, rdig=4)[,2]
  r<- ncol(v)
  lambda<- 1/ as.numeric(unique(tVarCorr[1:r]))
  message(paste(lambda, sep = " = ", collapse = ", "))
  sigma<- out_PEER$sigma
  sigma.e<- sigma
  
  #Returning output when se=F
  if(!se){
    status<- 0
    ret <- list(out_PEER, fitted.vals, GammaHat,
                AIC, BIC, logLik, Y, W, L_PEER,
                lambda, N, K, sigma.e, status)
    names(ret)<- c("fit", "fitted.vals", "GammaHat",
                   "AIC", "BIC", "logLik", "Y", "W", "L",
                   "lambda", "N", "K", "sigma", "status")
    
    class(ret) <- "peer"
    return(ret)
  }
  
  ###---- Standard Error
  tsigma<- as.numeric(unique(tVarCorr[1:K]))
  LL.inv<- lambda^(-2)*solve(t(L_PEER)%*%L_PEER)
  
  V.1<- W%*%LL.inv%*%t(W)+sigma.e^2*diag(rep(1, N))
  V<- sigma.e^2*diag(rep(1, N))
  
  X<- as.matrix(rep(1,N))
  X.Vinv.X.inv<- solve(t(X)%*%solve(V.1)%*%X)
  X.Vinv.Y<- t(X)%*%solve(V.1)%*%Y
  
  p1<- LL.inv%*%t(W)%*%solve(V.1)
  p2<- V.1 - X%*%X.Vinv.X.inv%*%t(X)
  p3<- solve(V.1)
  p<- p1%*%p2%*%p3
  vcov.gamma<- p%*%V%*%t(p)
  
  se.Gamma<- sqrt(diag(vcov.gamma))
  Gamma<- cbind(GammaHat, se.Gamma)
  colnames(Gamma)<- c('Estimate', 'SE')
  
  #Returning output when se=T
  status<- 1
  ret <- list(out_PEER, fitted.vals, Gamma,  V, V.1,
              GammaHat, vcov.gamma, se.Gamma, AIC, BIC, logLik, Y, W, L_PEER,
              lambda, N, K, sigma.e, status)
  names(ret)<- c("fit", "fitted.vals", "Gamma", "V", "V.1",
                 "GammaHat", "vcov.gamma", "se.Gamma", "AIC", "BIC", "logLik", "Y", "W", "L",
                 "lambda", "N", "K", "sigma", "status")
  
  class(ret) <- "peer"
  ret
}

