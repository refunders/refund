\name{cv.sigcom}
\alias{cv.sigcom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Cross-validation for linear function-on-function regression
%%  ~~function to do ... ~~
}
\description{This function is used to perform cross-validation and build the final model using the signal compression approach for the following linear function-on-function regression model that includes functional (and scalar) predictors:
\deqn{Y(t)= \mu(t)+Z^T\alpha(t)+\sum_{i=1}^p\int_{a_i}^{b_i}X_i(s)\beta_i(s,t)ds+\varepsilon(t)}{Y(t)= \mu(t)+Z'\alpha(t)+\int X_1(s)\beta_1(s,t)ds+...+\int X_p(s)\beta_p(s,t)ds+\epsilon(t),}
where \eqn{\mu(t)}{\mu(t)} is the intercept function. \eqn{Z}{Z} is a
multivariate predictor and \eqn{\alpha(t)}{\alpha(t)} is the vector of
corresponding coefficient functions. \eqn{Z=}NULL means that there is no
scalar predictor. The \eqn{\{X_i(s),1\le i\le p\}}{{X_i(s),1\le i\le p}}
are \eqn{p} functional predictors and \eqn{\{\beta_i(s,t),1\le i\le
  p\}}{{\beta_i(s,t),1\le i\le p}} are their corresponding coefficient
functions, where \eqn{p} is a positive integer. When \eqn{p} is large
(e.g., greater than 6), one may want to consider sparsity penalty (see
\code{\link{cv.hd}}()). The \eqn{\epsilon(t)}{\epsilon(t)} is the noise function.

We require that all the sample curves of each functional predictor are observed in a common dense grid of time points, but the grid can be different for different predictors. All the sample curves of the functional response are observed in a common dense grid.

This method estimates a special decomposition of the coefficient functions:
\deqn{\beta_i(s,t)=\sum_{k=1}^\infty\psi_{ik}(s)w_k(t),  1\le i\le p}{\beta_i(s,t)=\psi_{i1}(s)w_1(t)+...+\psi_{iK}(s)w_K(t),   1\le i\le p.}
We first estimate \eqn{\{\psi_{ik}(s), 1\le i\le p\}}{{\psi_{ik}(s), 1\le i\le p}} for each \eqn{k>0} by solving a generalized functional eigenvalue problem with penalty
\deqn{\lambda\sum_{i=1}^p\{|| \psi_{ik}||^2+ \tau ||\psi''_{ik}||^2\}}{\lambda{|| \psi_{1k}||^2+...+  ||\psi_{pk}^2||^2} + \lambda\tau{|| \psi''_{1k}||^2+...+||\psi''_{pk}||^2},}
%\deqn{\lambda\sum_{i=1}^p[\int \psi_{ik}^2(s)ds+ \tau\int \psi''_{ik}^2(s)ds]}{\lambda[\int \psi_{1k}^2(s)ds+...+ \int \psi_{pk}^2(s)ds] + \lambda\tau[\int \psi''_{1k}^2(s)ds+...+\int \psi''_{pk}^2(s)ds],}
and then estimate \eqn{w_{k}(t), k>0}{w_{k}(t), k>0} by penalized least square with smoothness tuning parameter \eqn{\kappa}.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cv.sigcom(X, Y, t.x, t.y, Z = NULL, s.n.basis = 50, t.n.basis = 50,
           K.cv = 5, upper.comp = 10, thresh = 0.005, tol = 1e-12)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a list of length \eqn{p}, the number of functional predictors. Its \eqn{i}-th element is the \eqn{n\times m_i}{n*m_i} data matrix for the \eqn{i}-th functional predictor \eqn{X_i(s)}{X_i(s)}, where \eqn{n} is the sample size and \eqn{m_i} is the number of observation time points for  \eqn{X_i(s)}{X_i(s)}.
%%     ~~Describe \code{X} here~~
}
  \item{Y}{the \eqn{n\times m}{n*m} data matrix for the functional response \eqn{Y(t)}{Y(t)}, where  \eqn{n}  is the sample size and  \eqn{m}  is the number of the observation time points for \eqn{Y(t)}{Y(t)}.
%%     ~~Describe \code{Y} here~~
}
  \item{t.x}{a list of length \eqn{p}. Its \eqn{i}-th element is the vector of obesrvation time points of the \eqn{i}-th functional predictor \eqn{X_i(s)}{X_i(s)}, \eqn{1\le i\le p}{1\le i\le p}.
%%     ~~Describe \code{t.x} here~~
}
  \item{t.y}{the vector of obesrvation time points of the functional response \eqn{Y(t)}{Y(t)}.
%%     ~~Describe \code{t.y} here~~
}
  \item{Z}{the \eqn{n\times q}{n*q} data matrix for multivariate scalar predictors, where  \eqn{n}  is the sample size and  \eqn{q}  is the number of the scalar predictors. Default is NULL, indicating no scalar predictors.
%%     ~~Describe \code{Z} here~~
}
  \item{s.n.basis}{the number of B-spline basis functions used for estimating the functions \eqn{\psi_{ik}(s)}{\psi_{ik}(s)'s}. Default is 50.
%%     ~~Describe \code{s.n.basis} here~~
}
  \item{t.n.basis}{the number of B-spline basis functions used for estimating the functions \eqn{w_{k}(t)}{w_{k}(t)'s}. Default is 50.
%%     ~~Describe \code{t.n.basis} here~~
}
  \item{K.cv}{the number of CV folds. Default is 5.
%%     ~~Describe \code{K.cv} here~~
}

  \item{upper.comp}{the upper bound for the maximum number of components to be calculated. Default is 10.
%%     ~~Describe \code{upper.comp} here~~
}
  \item{thresh}{a number between 0 and 1 specifying the minimum proportion of variation to be explained by each selected component relative to all the selected components. This will help determining the upper bound of the number of components to be used in CV, and the optimal number of components will be determined from 1,2,..., to the minimum of this uppber bound and upper.com by CV. A smaller thresh value leads to more components to be selected and a longer running time. A larger thresh value needs less running time, but may miss some important components and lead to a larger prediction error. Default is 0.005.
%%     ~~Describe \code{thresh} here~~
}
  \item{tol}{the iteration convergence tolerance.
%%     ~~Describe \code{tol} here~~
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{An object of the ``cv.sigcom'' class, which is used in the function \code{\link{pred.sigcom}} for prediction and \code{\link{getcoef.sigcom}} for extracting the estimated coefficient functions.
\item{t.x.list}{the input t.x.}
\item{t.y}{the input t.y.}
\item{Z}{the input Z}
\item{opt.index}{index of the optimal \eqn{lambda}.}
\item{opt.lambda}{optimal value for \eqn{lambda}, the parameter tuning the penalty on \eqn{{\psi_{ik}(s), 1\le i \le p}}.}
\item{opt.tau}{optimal value for \eqn{tau}, the smoothness parameter for \eqn{{\psi_{ik}(s), 1\le i \le p}}.}
\item{opt.kappa}{optimal value for \eqn{kappa}, the smoothness parameter for \eqn{w_k(t)}{w_k(t)}.}
\item{opt.K}{optimal number of components to be selected.}
\item{min.error}{minimum CV error.}
\item{errors}{list for CV errors. }
\item{x.smooth.params}{a list for internal use.}
\item{y.smooth.params}{a list for internal use.}
\item{fit.1}{a list for internal use.}
\item{is.null.Z}{a logic value indicating whether \eqn{Z} is null.}
}
\references{Ruiyan Luo and Xin Qi (2017) Function-on-Function Linear Regression by Signal Compression, Journal of the American Statistical Association.  112(518), 690-705. \url{ https://doi.org/10.1080/01621459.2016.1164053}
%% ~put references to the literature/web site here ~
}
\author{Ruiyan Luo and Xin Qi
%%  ~~who you are~~
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
#################################################################
# Example 1: function-on-function model without scalar predictors
#########################################################

ptm <- proc.time()
library(MASS)

n.curves <- 2 #number of predictor curves
nnew <- 50 # number of observations in new data
ntrain <- 10 # number of observations in training data
ntot <- nnew+ntrain
t.x <- seq(0,2,length=50) # all the four predictor curves are observed
              # at 50 equally spaced points in [0,2].
t.y <- seq(0,1,length=50) # the response curve is observed at 50
              # equally spaced points in [0,1].

# functions for mu(t) and beta_i(s,t)
mu.fun <- function(t){
    2*exp(-(t-1)^2)*sin(2*pi*t)
}
beta.fun.1 <- function(s,t)
{
    sin(1.5*pi*s)*sin( pi*t)
}
beta.fun.2 <- function(s,t)
{
    exp(-5*(s-0.5)^2-3*(t-0.5))+2*exp(-5*(s-1.5)^2-3*(t-0.5))
}

#generate the predictor curves using sin() and cos() basis functions.
 X=list()
 Sigma=matrix(0.5, 2, 2)
 diag(Sigma)=1
 sins=lapply(1:5, function(i){t(sin(i*pi*t.x))/i})
 coss=lapply(1:5, function(i){t(cos(i*pi*t.x))/i})
 coefs1=lapply(1:5, function(i){mvrnorm(ntot,rep(0,2),Sigma)})
 coefs2=lapply(1:5, function(i){mvrnorm(ntot,rep(0,2),Sigma)})
 X[[1]]=Reduce("+",lapply(1:5, function(i){coefs1[[i]][,1] \%*\% sins[[i]]
             +coefs2[[i]][,1] \%*\% coss[[i]]}))
 X[[2]]=Reduce("+",lapply(1:5, function(i){coefs1[[i]][,2] \%*\% sins[[i]]
             +coefs2[[i]][,2] \%*\% coss[[i]]}))
 mu.val=mu.fun(t.y)
 beta.val=list()
 beta.val[[1]]<- outer(t.x,t.y,beta.fun.1)
 beta.val[[2]]<- outer(t.x,t.y,beta.fun.2)
# generate sample curves for the noise and the respose function.
 E<-matrix(rnorm(ntot*length(t.y)),ntot, length(t.y))
 delta<-(max(t.x)-min(t.x))/length(t.x)
 Y<-t(sapply(1:ntot, function(i){mu.val+(X[[1]][i,]\%*\% beta.val[[1]]
   +X[[2]][i,]\%*\% beta.val[[2]])*(max(t.x)-min(t.x))/length(t.x)+ E[i,]}))

# gereate the training data and perform CV
 X.train=lapply(1:n.curves, function(j){X[[j]][1:ntrain,]})
 t.x.list=lapply(1:n.curves,function(j){t.x})
 Y.train <- Y[1:ntrain, ]
 fit.cv=cv.sigcom(X.train, Y.train, t.x.list, t.y, s.n.basis=20, t.n.basis=20)

# prediction and estimation error
 X.new=lapply(1:n.curves, function(j){X[[j]][-(1:ntrain),]})
 Y.new <- Y[-(1:ntrain), ]
 E.new=E[-(1:ntrain), ]
 Y.pred=pred.sigcom(fit.cv, X.new)
 error <- mean((Y.pred-Y.new)^2)
 print(c(" prediction error=", error))
 est.error <- mean((Y.pred-Y.new+E.new)^2)
# extract the estimated intercept and coefficient functions
 print(c("estimation error for regression function (or signal function)=", est.error))
 coef.est.all=getcoef.sigcom(fit.cv)
 mu.est=coef.est.all[[1]]
 beta.est=coef.est.all[[2]]
print(proc.time()-ptm)

#################################################################
# Example 2: function-on-function model with scalar predictors
#########################################################

\donttest{
ptm <- proc.time()
library(MASS)

n.curves <- 2 # there are two predictor curves. They have different range
              # and observation points
ntest <- 500
ntrain <- 100
ntot <- ntest+ntrain
t.x.1=seq(0,2,length=100) # the vector of observation points for X_1(s)
t.x.2=seq(0,1,length=50) # the vector of observation points for X_1(s)
t.y=seq(0,1,length=60) # the vector of observation points for Y(t)

#########################################
# functions for mu(t) and beta_i(s,t)
##########################################
mu.fun <- function(t){
    cos(4*pi*t)
}

beta.fun.1 <- function(s,t)
{
  1/(1+s^2+t^2)
}
beta.fun.2 <- function(s,t)
{
   log(1+s*t)
}
##########################################
X=list()
X[[1]]=X[[2]]=0
Sigma=matrix(0.3, 2, 2)
diag(Sigma)=1
for(k in 1:20){
   coef.mtx.1=mvrnorm(ntot, rep(0,2), Sigma)
   coef.mtx.2=mvrnorm(ntot, rep(0,2), Sigma)
   X[[1]]=X[[1]]+coef.mtx.1[,1]\%*\% t(sin(k*pi*t.x.1))/k^2
                +coef.mtx.2[,1]\%*\% t(cos(k*pi*t.x.1))/k^2
   X[[2]]=X[[2]]+coef.mtx.1[,2]\%*\% t(sin(k*pi*t.x.2))/k^2
                +coef.mtx.2[,2]\%*\% t(cos(k*pi*t.x.2))/k^2
  }
 mu.val=mu.fun(t.y)
 beta.val=list()
 beta.val[[1]]<- outer(t.x.1,t.y,beta.fun.1)
 beta.val[[2]]<- outer(t.x.2,t.y,beta.fun.2)

#  generate the scalar predictors and the corresponding coefficient functions
 Sigma.scalar=matrix(0.6, 3, 3)
 diag(Sigma.scalar)=1
 Z=mvrnorm(ntot, rep(0,3), Sigma.scalar)
 beta.scalar=rbind(exp(-t.y), 1+t.y^2, sin(3*pi*t.y))
 E<-matrix(rnorm(ntot*length(t.y)),ntot, length(t.y))
 Y<- t(sapply(1:ntot, function(i){mu.val
     +X[[1]][i,] \%*\% beta.val[[1]]*(max(t.x.1)-min(t.x.1))/length(t.x.1)
     +X[[2]][i,] \%*\% beta.val[[2]]*(max(t.x.2)-min(t.x.2))/length(t.x.2)+ E[i,]}))
 Y=Y+Z\%*\%beta.scalar
 X.train=lapply(1:n.curves, function(j){X[[j]][1:ntrain,]})
 X.test=lapply(1:n.curves, function(j){X[[j]][-(1:ntrain),]})
 t.x.list=list()
 t.x.list[[1]]=t.x.1
 t.x.list[[2]]=t.x.2
 Z.train<- Z[1:ntrain, ]
 Z.test <- Z[-(1:ntrain), ]
 Y.train <- Y[1:ntrain, ]
 Y.test <- Y[-(1:ntrain), ]
 E.test=E[-(1:ntrain), ]
 fit.cv=cv.sigcom(X.train, Y.train, t.x.list, t.y, Z.train, s.n.basis=40,
                   t.n.basis=40)
 Y.pred=pred.sigcom(fit.cv, X.test, Z.test)
 error<- mean((Y.pred-Y.test)^2)
 print(c(" prediction error=", error))
 est.error<- mean((Y.pred-Y.test+E.test)^2)
 print(c("estimation error for regression function (or signal function)=", est.error))
 coef.est.all=getcoef.sigcom(fit.cv)
 mu.est=coef.est.all[[1]]
 beta.est=coef.est.all[[2]]
 beta.scalar.est=coef.est.all[[3]]
 print(proc.time()-ptm)
#   user  system elapsed
#   2.90    0.69    3.96
}

#################################################################
#Example 3: application to the DTI data in 'refund' package
#########################################################

\donttest{
ptm <- proc.time()
library(refund)
data(DTI)
I=which(is.na(apply(DTI$cca,1,mean)))
Y=DTI$cca[-I,]
X=list(DTI$rcst[-I,-(1:12)])
ntot=dim(Y)[1]

lengthX=dim(X[[1]])[2]
lengthY=dim(Y)[2]
t.x <- list(seq(0,1,length=lengthX))
t.y <- seq(0,1,length=lengthY)
# randomly split all the observations into a training set with 200 observations
# and a test set.
ntrain=200
I=sample(1:ntot, ntrain)
X.train <- list(X[[1]][I,])
Y.train <- Y[I, ]
X.test <- list(X[[1]][-I,])
Y.test <- Y[-I, ]
fit.cv=cv.sigcom(X.train, Y.train, t.x, t.y, s.n.basis=40, t.n.basis=40)
Y.pred=pred.sigcom(fit.cv, X.test)
error<- mean((Y.pred-Y.test)^2)
print(c(" prediction error=", error))
coef.est.all=getcoef.sigcom(fit.cv)
mu.est=coef.est.all[[1]]
beta.est=coef.est.all[[2]]
print(proc.time()-ptm)
#   user  system elapsed
#   6.09    0.10    6.17
}
}