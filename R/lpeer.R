##' Longitudinal Functional Models with Structured Penalties
##'
##' Implements longitudinal functional model with structured penalties (Kundu
##' et al., 2012) with scalar outcome, single functional predictor, one or more
##' scalar covariates and subject-specific random intercepts through mixed
##' model equivalence.
##'
##' If there are any missing or infinite values in \code{Y}, \code{subj},
##' \code{t}, \code{covariates}, \code{funcs} and \code{f_t}, the corresponding
##' row (or observation) will be dropped, and infinite values are not allowed
##' for these arguments. Neither \code{Q} nor \code{L} may contain missing or
##' infinite values.  \code{lpeer()} fits the following model:
##'
##' \eqn{y_{i(t)}=X_{i(t)}^T \beta+\int {W_{i(t)}(s)\gamma(t,s) ds}
##' +Z_{i(t)}u_i + \epsilon_{i(t)}}
##'
##' where \eqn{\epsilon_{i(t)} ~ N(0,\sigma ^2)} and \eqn{u_i ~ N(0,
##' \sigma_u^2)}.  For all the observations, predictor function
##' \eqn{W_{i(t)}(s)} is evaluated at K sampling points. Here, regression
##' function \eqn{\gamma (t,s)} is represented in terms of (d+1) component
##' functions \eqn{\gamma_0(s)},..., \eqn{\gamma_d(s)} as follows
##'
##' \eqn{\gamma (t,s)= \gamma_0(s)+f_1(t) \gamma_1(s) + f_d(t) \gamma_d(s)}
##'
##' Values of \eqn{y_{i(t)} , X_{i(t)}} and \eqn{W_{i(t)}(s)} are passed
##' through argument \code{Y}, \code{covariates} and \code{funcs},
##' respectively. Number of elements or rows in \code{Y}, \code{t},
##' \code{subj}, \code{covariates} (if not \code{NULL}) and \code{funcs} need
##' to be equal.
##'
##' Values of \eqn{f_1(t),...,f_d(t)} are passed through f_t argument. The
##' matrix passed through \code{f_t} argument should have d columns where each
##' column represents one and only one of \eqn{f_1(t),..., f_d(t)}.
##'
##' The estimate of (d+1) component functions \eqn{\gamma_0(s)},...,
##' \eqn{\gamma_d(s)} is obtained as penalized estimated. The following 3 types
##' of penalties can be used for a component function:
##'
##' i.  Ridge: \eqn{I_K}
##'
##' ii.  Second-order difference: [\eqn{d_{i,j}}] with \eqn{d_{i,i} = d_{i,i+2}
##' = 1, d_{i,i+1} = -2}, otherwise \eqn{d_{i,j} =0}
##'
##' iii. Decomposition based penalty: \eqn{bP_Q+a(I-P_Q)} where \eqn{P_Q= Q^T
##' (QQ^T)^{-1}Q}
##'
##' For Decomposition based penalty the user must specify \code{pentype=
##' 'DECOMP'} and the associated Q matrix must be passed through the \code{Q}
##' argument. Alternatively, one can directly specify the penalty matrix by
##' setting \code{pentype= 'USER'} and using the \code{L} argument to supply
##' the associated L matrix.
##'
##' If Q (or L) matrix is similar for all the component functions then argument
##' \code{comm.pen} should have value \code{TRUE} and in that case specified
##' matrix to argument \code{Q} (or \code{L}) should have K columns. When Q (or
##' L) matrix is different for all the component functions then argument
##' \code{comm.pen} should have value \code{FALSE} and in that case specified
##' matrix to argument \code{Q} (or \code{L}) should have K(d+1) columns. Here
##' first K columns pertains to first component function, second K columns
##' pertains to second component functions, and so on.
##'
##' Default penalty is Ridge penalty for all the component functions and user
##' needs to specify \code{'RIDGE'}. For second-order difference penalty, user
##' needs to specify \code{'D2'}. When pentype is \code{'RIDGE'} or \code{'D2'}
##' the value of \code{comm.pen} is always \code{TRUE} and
##' \code{comm.pen=FALSE} will be ignored.
##'
##' @param Y vector of all outcomes over all visits or timepoints
##' @param subj vector containing the subject number for each observation
##' @param t vector containing the time information when the observation are
##' taken
##' @param covariates matrix of scalar covariates.
##' @param funcs matrix containing observed functional predictors as rows. Rows
##' with \code{NA} and \code{Inf} values will be deleted.
##' @param argvals matrix (or vector) of indices of evaluations of \eqn{X_i(t)}; i.e. a matrix with
##' \emph{i}th row \eqn{(t_{i1},.,t_{iJ})}
##' @param comm.pen logical value indicating whether common penalty for all the
##' components of regression function. Default is TRUE.
##' @param pentype type of penalty: either decomposition based penalty
##' (\code{'DECOMP'}) or ridge (\code{'RIDGE'}) or second-order difference
##' penalty (\code{'D2'}) or any user defined penalty (\code{'USER'}). For
##' decomposition based penalty user need to specify Q matrix in \code{Q}
##' argument (see details). For user defined penalty user need to specify L
##' matrix in \code{L} argument (see details). For Ridge and second-order
##' difference penalty, specification for arguments \code{L} and \code{Q} will
##' be ignored. Default is \code{'RIDGE'}.
##' @param f_t vector or matrix with number of rows equal to number of total
##' observations and number of columns equal to d (see details). If matrix then
##' each column pertains to single function of time and the value in the column
##' represents the realization corresponding to time vector t. The column with
##' intercept or multiple of intercept will be dropped. A \code{NULL} value
##' refers to time-invariant regression function. Default value is \code{NULL}.
##' @param Q Q matrix to derive decomposition based penalty. Need to be
##' specified with \code{pentype='DECOMP'}. When \code{comm.pen=TRUE}, number
##' of columns must equal number of columns of matrix specified to
##' \code{funcs}. When \code{comm.pen=FALSE}, Number of columns need to be
##' equal with the number of columns of matrix specified to \code{funcs} times
##' the number of components of regression function. Each row represents a
##' basis function where functional predictor is expected lie according to
##' prior belief. This argument will be ignored when value of \code{pentype} is
##' other than \code{'DECOMP'}.
##' @param L.user penalty matrix. Need to be specified with
##' \code{pentype='USER'}. When \code{comm.pen=TRUE}, Number of columns need to
##' be equal with number of columns of matrix specified to \code{funcs}. When
##' \code{comm.pen=FALSE}, Number of columns need to be equal with the number
##' of columns of matrix specified to \code{funcs} times the number of
##' components of regression function. Each row represents a constraint on
##' functional predictor. This argument will be ignored when value of
##' \code{pentype} is other than \code{'USER'}.
##' @param phia scalar value of a in decomposition based penalty. Needs to be
##' specified with \code{pentype='DECOMP'}.
##' @param se logical; calculate standard error when \code{TRUE}.
##' @param ... additional arguments passed to \code{\link[nlme]{lme}}.
##' @return A list containing: \item{fit}{result of the call to \code{lme}}
##' \item{fitted.vals}{predicted outcomes} \item{BetaHat}{parameter estimates
##' for scalar covariates including intercept} \item{se.Beta}{standard error of
##' parameter estimates for scalar covariates including intercept}
##' \item{Beta}{parameter estimates with standard error for scalar covariates
##' including intercept} \item{GammaHat}{estimates of components of regression
##' functions. Each column represents one component function. }
##' \item{Se.Gamma}{standard error associated with \code{GammaHat}}
##' \item{AIC}{AIC value of fit (smaller is better) } \item{BIC}{BIC value of
##' fit (smaller is better) } \item{logLik}{(restricted) log-likelihood at
##' convergence} \item{lambda}{list of estimated smoothing parameters
##' associated with each component function} \item{V}{conditional variance of Y
##' treating only random intercept as random one. } \item{V1}{unconditional
##' variance of Y } \item{N}{number of subjects} \item{K}{number of Sampling
##' points in functional predictor} \item{TotalObs}{total number of
##' observations over all subjects} \item{Sigma.u}{estimated sd of random
##' intercept. } \item{sigma}{estimated within-group error standard deviation.}
##' @author Madan Gopal Kundu \email{mgkundu@@iupui.edu}
##' @seealso \code{peer}, \code{plot.lpeer}
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
##' # Example 1: Estimation with Ridge penalty
##' #------------------------------------------------------------------------
##'
##' ##Load Data
##' data(DTI)
##'
##' ## Extract values for arguments for lpeer() from given data
##' cca = DTI$cca[which(DTI$case == 1),]
##' DTI = DTI[which(DTI$case == 1),]
##'
##' ##1.1 Fit the model with single component function
##' ##    gamma(t,s)=gamm0(s)
##' t<- DTI$visit
##' fit.cca.lpeer1 = lpeer(Y=DTI$pasat, t=t, subj=DTI$ID, funcs = cca)
##' plot(fit.cca.lpeer1)
##'
##' ##1.2 Fit the model with two component function
##' ##    gamma(t,s)=gamm0(s) + t*gamma1(s)
##' fit.cca.lpeer2 = lpeer(Y=DTI$pasat, t=t, subj=DTI$ID, funcs = cca,
##'                       f_t=t, se=TRUE)
##' plot(fit.cca.lpeer2)
##'
##' #------------------------------------------------------------------------
##' # Example 2: Estimation with structured penalty (need structural
##' #            information about regression function or predictor function)
##' #------------------------------------------------------------------------
##'
##' ##Load Data
##' data(PEER.Sim)
##'
##' ## Extract values for arguments for lpeer() from given data
##' K<- 100
##' W<- PEER.Sim[,c(3:(K+2))]
##' Y<- PEER.Sim[,K+3]
##' t<- PEER.Sim[,2]
##' id<- PEER.Sim[,1]
##'
##' ##Load Q matrix containing structural information
##' data(Q)
##'
##' ##2.1 Fit the model with two component function
##' ##    gamma(t,s)=gamm0(s) + t*gamma1(s)
##' Fit1<- lpeer(Y=Y, subj=id, t=t, covariates=cbind(t), funcs=W,
##' 	    pentype='DECOMP', f_t=cbind(1,t), Q=Q, se=TRUE)
##'
##' Fit1$Beta
##' plot(Fit1)
##'
##' ##2.2 Fit the model with three component function
##' ##    gamma(t,s)=gamm0(s) + t*gamma1(s) + t^2*gamma1(s)
##' Fit2<- lpeer(Y=Y, subj=id, t=t, covariates=cbind(t), funcs=W,
##' 		     pentype='DECOMP', f_t=cbind(1,t, t^2), Q=Q, se=TRUE)
##'
##' Fit2$Beta
##' plot(Fit2)
##'
##' ##2.3 Fit the model with two component function with different penalties
##' ##    gamma(t,s)=gamm0(s) + t*gamma1(s)
##' Q1<- cbind(Q, Q)
##' Fit3<- lpeer(Y=Y, subj=id, t=t, covariates=cbind(t), comm.pen=FALSE, funcs=W,
##' 		     pentype='DECOMP', f_t=cbind(1,t), Q=Q1, se=TRUE)
##'
##' ##2.4 Fit the model with two component function with user defined penalties
##' ##    gamma(t,s)=gamm0(s) + t*gamma1(s)
##' phia<- 10^3
##' P_Q <- t(Q)%*%solve(Q%*%t(Q))%*% Q
##' L<- phia*(diag(K)- P_Q) + 1*P_Q
##' Fit4<- lpeer(Y=Y, subj=id, t=t, covariates=cbind(t), funcs=W,
##' 		     pentype='USER', f_t=cbind(1,t), L=L, se=TRUE)
##'
##' L1<- adiag(L, L)
##' Fit5<- lpeer(Y=Y, subj=id, t=t, covariates=cbind(t), comm.pen=FALSE, funcs=W,
##' 		     pentype='USER', f_t=cbind(1,t), L=L1, se=TRUE)
##' }
##'
##' @export
##' @importFrom nlme pdIdent pdBlocked VarCorr
##' @importFrom magic adiag
lpeer<- function(Y, subj, t, funcs, argvals=NULL, covariates=NULL, comm.pen=TRUE,  
                 pentype='Ridge', L.user=NULL, f_t=NULL, Q=NULL, phia=10^3,
                 se=FALSE, ...)
{
  if (!is.null(argvals)) 
    stop("argvals is not implemented in this version of refund")
  
  pd1 = pd2 = pd3 = pd4 = pd5 = pd6 = NULL

  #Determining K, converting W, Y, id and t to matrix
  W<- as.matrix(funcs)
  K<- ncol(W)
  Y<- as.matrix(Y)
  id<- as.matrix(subj)
  t<-as.matrix(t)

  #Check 1:Making sure Y, subj and t have only 1 column
  if(dim(Y)[2]>1) return(cat("Error: No. of column for Y cannot be greater than 1. \nThe lpeer() will not proceed further.\n"))
  if(dim(id)[2]>1) return(cat("Error: No. of column for subj cannot be greater than 1. \nThe lpeer() will not proceed further.\n"))
  if(dim(t)[2]>1) return(cat("Error: No. of column for t cannot be greater than 1. \nThe lpeer() will not proceed further.\n"))

  #Check 2: Do check for intercept in X matrix
  if (!is.null(covariates)){
    covariates<- as.matrix(covariates)
    X.chk<- apply(covariates, 2, sd)
    if(any(X.chk==0)) return(cat("Error: Drop intercept or equivalent to intercept term from covariate. \nThe lpeer() will not proceed further.\n"))
  }
  X<-  cbind(1, covariates)

  #Check 3: Check the dimension of Y, id, t, W and X
  Yl<- dim(Y)[1]
  idl<- dim(id)[1]
  tl<- dim(t)[1]
  chk.eq<- ifelse(Yl==idl & idl==tl & tl==nrow(W), 0 ,1)
  if(chk.eq==1) return(cat("Error: At least one of (1) length of Y, (2) lenght of subj, (3) length of \nt, and (4) number of row of funcs are not equal.\n The lpeer() will not proceed further.\n"))
  if(!is.null(covariates) & Yl!=nrow(cbind(X,X))) return(cat("Error: length of Y and number of rows of X is not equal.\n The lpeer() will not proceed further.\n"))

  #Organizing f(t)
  if(length(dim(f_t))>0 ){
    f_t<- f_t
  } else
    if(length(f_t)>0){
      f_t<- matrix(f_t, ncol=1)
    } else
    {
      f_t<- matrix(rep(1, Yl), ncol=1)
    }
  f_t<- f_t[,which(apply(f_t, 2, sd)>0)]
  f_t<- cbind(1, f_t)
  d=ncol(f_t)-1
  if(d>5) warning("Only first 5 time components will be used", call. = FALSE)
  d=min(d,5)

  #Check 4: check in f(t)
  if(dim(f_t)[1]!=Yl) return(cat("Error: f_t and Y are not compatible in dimension. \nThe lpeer() will not proceed further.\n"))

  #Sort the data by id and t & removal of missing and infinite observations
  tdata<- data.frame(id, t, Y, W, X, f_t)
  tdata<- tdata[which(apply(is.na(tdata), 1, sum)==0),]
  tdata<- tdata[which(apply(is.infinite(as.matrix(tdata)), 1, sum)==0),]
  tdata<- tdata[order(id, t), ]
  tdata<- tdata[!is.na(tdata$id) & !is.na(tdata$t) & !is.na(tdata$Y),]
  id<- tdata$id
  t<- tdata$t
  Y<- tdata$Y
  W.e<- dim(W)[2]+3; W<- as.matrix(tdata[,c(4:W.e)])
  X.s<- W.e + 1; X.e<- dim(X)[2]+W.e; X<- as.matrix(tdata[,c(X.s:X.e)])
  f_t.s<- X.e + 1; f_t.e<- dim(f_t)[2]+X.e; f_t<- as.matrix(tdata[,c(f_t.s:f_t.e)])

  #Determining N and NT
  N<- length(unique(id))
  NT<- length(Y)

  #Checking entry for pentype
  pentypecheck<- toupper(pentype) %in% c('DECOMP', 'DECOMPOSITION', 'RIDGE', 'D2', 'USER')
  if(!pentypecheck)  return (cat("Error: Specify valid object for argument PENTYPE.\n"))

  #Check 5: Some checking/processing for decomposition type of penalty
  if(toupper(pentype)=='DECOMP' | toupper(pentype)=='DECOMPOSITION'){


    #5.1: Removing rows containing missing and infinite values
    Q<- Q[which(apply(is.na(Q), 1, sum)==0),]
    Q<- Q[which(apply(is.infinite(Q), 1, sum)==0),]

    #5.2: Compatibility of Q and W matrix
    if (!comm.pen)
    {
      if(ncol(Q)!=(d+1)*ncol(W)) return(cat('Error: For different penalty, number of columns of Q need to be (d+1) \ntimes of number of columns of funcs.\nThe lpeer() will not proceed further.\n'))
    }
    if (comm.pen)
    {
      if(ncol(Q)!=ncol(W)) return(cat('Error: For common penalty, number of columns of func and Q need to be equal.\nThe lpeer() will not proceed further.\n'))
      Q1<- Q
      for(i in 1:d) Q1<- cbind(Q1, Q)
      Q<- Q1
    }

    #5.3: Singularity of Q matrix
    Q.eig<- abs(eigen(Q %*% t(Q))$values)
    if(any(Q.eig<1e-12)) return(cat('Error: Q matrix is singular or near singular.\nThe lpeer() will not proceed further.\n'))

    #5.4: Checking for phia
    if(!exists("phia")) return (cat("Error: Specify valid object for argument PHIA.\n"))
    if(!is.numeric(phia)|is.matrix(phia)|is.matrix(phia)) return (cat("Error: Specify valid object for argument PHIA.\n"))
  }

  #Check 6: Some checking/processing for user type of penalty
  if(toupper(pentype)=='USER'){
    L<- L.user

    #6.1: Removing rows containing missing and infinite values
    L<- L[which(apply(is.na(L), 1, sum)==0),]
    L<- L[which(apply(is.infinite(L), 1, sum)==0),]

    #6.2: Dimension of L matrix
    if (!comm.pen)
    {
      if(ncol(L)!=(d+1)*ncol(W)) return(cat('Error: For different penalty, number of columns of L need to be (d+1) \ntimes of number of columns of func.\nThe lpeer() will not proceed further.\n'))
    }
    if (comm.pen)
    {
      if(ncol(L)!=ncol(W)) return(cat('Error: For common penalty, number of columns of func and L.user need to be equal.\nThe lpeer() will not proceed further.\n'))
      L1<- L
      for(i in 1:d) L1<- magic::adiag(L1, L)
      L<- L1
    }

    #6.3: Singularity of L'L matrix
    LL<- t(L)%*%L
    LL.eig<- abs(eigen(LL %*% t(LL))$values)
    if(any(LL.eig<1e-12)) return(cat("Error: L'L matrix is singular or near singular.\nThe lpeer() will not proceed further.\n"))
  }

  #Generate L matrix for D2 penalty
  if(toupper(pentype)=='D2'){
    Left<- cbind(diag(rep(1,K-2)),rep(0,K-2),rep(0,K-2))
    Middle<- cbind(rep(0,K-2),diag(rep(-2,K-2)),rep(0,K-2))
    Right<- cbind(rep(0,K-2),rep(0,K-2),diag(rep(1,K-2)))
    D.2<- rbind(Left+Middle+Right, c(rep(0, K-2), 1, -2), c(rep(0, K-2), 0, 1))
  }

  #Generate W1 matrix (This is referred as in the paper W matrix)
  for(i in 0:d){
    if (i==0) W1<-data.matrix(W)*f_t[,(i+1)]
    if (i>0) W1<- cbind(W1, data.matrix(W)*f_t[,(i+1)])
  }

  #Generate W* matrix
  for(i in 0:d){
    if(toupper(pentype)=='DECOMP' | toupper(pentype)=='DECOMPOSITION'){
      tQ<- Q[,(i*K+1):((i+1)*K)]
      tP_Q <- t(tQ) %*% solve(tQ %*% t(tQ)) %*% tQ
      tL_PEER<- phia*(diag(K)- tP_Q) + 1*tP_Q
      rm(tQ); rm(tP_Q)
    } else
      if(toupper(pentype)=='RIDGE'){
        tL_PEER<- diag(K)
      } else
        if(toupper(pentype)=='D2'){
          tL_PEER<- D.2
        } else
          if(toupper(pentype)=='USER'){
            tL_PEER<- L[(i*K+1):((i+1)*K),(i*K+1):((i+1)*K)]
          }

    v <- diag(K)
    if (K > N)
      v <- svd((data.matrix(W) * f_t[, (i + 1)]) %*% solve(tL_PEER))$v
    assign(paste("W", i + 1, "_PEER", sep = ""), (data.matrix(W) *
      f_t[, (i + 1)]) %*% solve(tL_PEER) %*% v)
    rm(tL_PEER)
    rm(v)
  }

  #Generate Z
  id.bd1<- factor(rep(1, NT))
  ni<- tapply(id, id, length)
  for(i in 1:N){
    if (i==1) Z<- matrix(1, nrow=ni[i])
    if (i>1) Z<- magic::adiag(Z, matrix(1, nrow=ni[i]))
  }

  #Input for random argument of lme function
  for(i in 0:d) assign(paste('pd', i+1, sep=''),
                       nlme::pdIdent(form=as.formula(paste('~W', i+1, '_PEER -1', sep=''))))
  pdid<- nlme::pdIdent(~Z-1)

  if(d==0) tXX<- nlme::pdBlocked(list(pd1, pdid))
  if(d==1) tXX<- nlme::pdBlocked(list(pd1, pd2, pdid))
  if(d==2) tXX<- nlme::pdBlocked(list(pd1, pd2, pd3, pdid))
  if(d==3) tXX<- nlme::pdBlocked(list(pd1, pd2, pd3, pd4, pdid))
  if(d==4) tXX<- nlme::pdBlocked(list(pd1, pd2, pd3, pd4, pd5, pdid))
  if(d==5) tXX<- nlme::pdBlocked(list(pd1, pd2, pd3, pd4, pd5, pd6, pdid))

  #Fitting the model
  out_PEER<- nlme::lme(fixed=Y~X-1, random=list(id.bd1=tXX), ... )
  cat('The fit is successful.\n')

  #Extracting the estimates
  Gamma_PEER<-matrix(out_PEER$coeff$random$id.bd1, ncol=1)
  for(i in 0:d) {
    if(toupper(pentype)=='DECOMP' | toupper(pentype)=='DECOMPOSITION'){
      tQ<- Q[,(i*K+1):((i+1)*K)]
      tP_Q <- t(tQ) %*% solve(tQ %*% t(tQ)) %*% tQ
      tL_PEER<- phia*(diag(K)- tP_Q) + 1*tP_Q
      rm(tQ); rm(tP_Q)
    } else
      if(toupper(pentype)=='RIDGE'){
        tL_PEER<- diag(K)
      } else
        if(toupper(pentype)=='D2'){
          tL_PEER<- D.2
        } else
          if(toupper(pentype)=='USER'){
            tL_PEER<- L[(i*K+1):((i+1)*K),(i*K+1):((i+1)*K)]
          }
    v <- diag(K)
    if (K > N)
      v <- svd((data.matrix(W) * f_t[, (i + 1)]) %*% solve(tL_PEER))$v
    r <- ncol(v)
    tGamma.PEER.hat <- matrix(Gamma_PEER[(i * r + 1):((i +
      1) * r)], ncol = 1)
    tGammaHat <- solve(tL_PEER) %*% v %*% tGamma.PEER.hat
    if (i == 0)
      GammaHat <- matrix(tGammaHat, ncol = 1)
    if (i > 0)
      GammaHat <- cbind(GammaHat, tGammaHat)
    rm(tL_PEER)
    rm(v)
    rm(tGammaHat)
  }
  colnames(GammaHat)<- paste('Gamma', 0:d, sep='')
  BetaHat<- summary(out_PEER)$tTable[,1]
  names(BetaHat)<- c('Intercept', colnames(covariates))
  fitted.vals<- summary(out_PEER)$fitted

  #Extracting model diagnostics
  logLik<- summary(out_PEER)$logLik
  AIC<- summary(out_PEER)$AIC
  BIC<- summary(out_PEER)$BIC

  #Extracting lambda and variance
  tVarCorr<- nlme::VarCorr(out_PEER, rdig=4)[,2]
  for(i in 0:d) assign(paste('lambda', i, sep=''),
                       1/ as.numeric(unique(tVarCorr[(i*r+1):((i+1)*r)])))
  for(i in 0:d)
  {
    tLambda<- 1/ as.numeric(unique(tVarCorr[(i*r+1):((i+1)*r)]))
    if(i==0) lambda<- tLambda
    if(i>0) lambda<- c(lambda, tLambda)
  }
  names(lambda)<- paste('lambda', 0:d, sep='')
  print(lambda)
  sd_int.est<- as.numeric(unique(tVarCorr[((d+1)*r+1):((d+1)*r+N)]))
  sigma<- out_PEER$sigma
  Sigma.u<- sd_int.est
  sigma.e<- sigma

  #Returning output when se=F
  if(!se){
    status<- 0
    ret <- list(out_PEER, BetaHat,  fitted.vals,
                GammaHat, AIC, BIC, logLik,
                lambda, N, K, NT, Sigma.u, sigma.e, d, status)
    names(ret)<- c("fit", "BetaHat", "fitted.vals",
                   "GammaHat", "AIC", "BIC", "logLik",
                   "lambda", "N", "K", "TotalObs", "Sigma.u", "sigma", "d", "status")
    class(ret) <- "lpeer"
    return(ret)
  }

  ###---- Standard Error
  for(i in 0:d)
  {
    if(toupper(pentype)=='DECOMP' | toupper(pentype)=='DECOMPOSITION'){
      tQ<- Q[,(i*K+1):((i+1)*K)]
      tP_Q <- t(tQ) %*% solve(tQ %*% t(tQ)) %*% tQ
      tL_PEER<- phia*(diag(K)- tP_Q) + 1*tP_Q
      rm(tQ); rm(tP_Q)
    } else
      if(toupper(pentype)=='RIDGE'){
        tL_PEER<- diag(K)
      } else
        if(toupper(pentype)=='D2'){
          tL_PEER<- D.2
        } else
          if(toupper(pentype)=='USER'){
            tL_PEER<- L[(i*K+1):((i+1)*K),(i*K+1):((i+1)*K)]
          }
    tsigma<- as.numeric(unique(tVarCorr[(i*r+1):((i+1)*r)]))
    if(i==0) LL.inv<- tsigma^2*solve(t(tL_PEER)%*%tL_PEER)
    if(i>0) LL.inv<- magic::adiag(LL.inv, tsigma^2*solve(t(tL_PEER)%*%tL_PEER))

    rm(tsigma); rm(tL_PEER)
  }


  rand.var<- Sigma.u^2*(Z%*%t(Z))

  V.1<- W1%*%LL.inv%*%t(W1)+rand.var+sigma.e^2*diag(rep(1, NT))
  V<- rand.var+sigma.e^2*diag(rep(1, NT))

  X.Vinv.X.inv<- solve(t(X)%*%solve(V.1)%*%X)
  X.Vinv.Y<- t(X)%*%solve(V.1)%*%Y

  se.Beta<- sqrt(diag(X.Vinv.X.inv%*%t(X)%*%solve(V.1)%*%V%*%solve(V.1)%*%X%*%X.Vinv.X.inv))
  names(se.Beta)<- c('Intercept', colnames(covariates))
  Beta<- cbind(BetaHat, se.Beta)

  p1<- LL.inv%*%t(W1)%*%solve(V.1)
  p2<- V.1 - X%*%X.Vinv.X.inv%*%t(X)
  p3<- solve(V.1)
  p<- p1%*%p2%*%p3

  SE.gamma<- sqrt(diag(p%*%V%*%t(p)))
  for(i in 0:d)
  {
    if(i==0) se.Gamma<- matrix(SE.gamma[(i*K+1):((i+1)*K)], ncol=1)
    if(i>0) se.Gamma<- cbind(se.Gamma, SE.gamma[(i*K+1):((i+1)*K)])
  }
  colnames(se.Gamma)<- paste('Gamma', 0:d, sep='')

  #Returning output when se=T
  status<- 1
  ret <- list(out_PEER, BetaHat,  se.Beta, Beta, fitted.vals,
              GammaHat, se.Gamma, AIC, BIC, logLik, V.1,
              V, lambda, N, K, NT, Sigma.u, sigma.e, d, status)
  names(ret)<- c("fit", "BetaHat", "se.Beta", "Beta", "fitted.vals",
                 "GammaHat", "se.Gamma", "AIC", "BIC", "logLik", "V1",
                 "V", "lambda", "N", "K", "TotalObs", "Sigma.u", "sigma", "d", "status")

  class(ret) <- "lpeer"
  ret
}
