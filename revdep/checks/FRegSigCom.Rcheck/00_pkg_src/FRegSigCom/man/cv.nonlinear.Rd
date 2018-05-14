\name{cv.nonlinear}
\alias{cv.nonlinear}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Cross-validation for nonlinear function-on-function regression
%%  ~~function to do ... ~~
}
\description{This function is used to perform cross-validation and build the final model using signal compression approach for the following nonlinear function-on-function regression model:
\deqn{Y(t)= \mu(t)+\sum_{i=1}^p\int_{a_i}^{b_i}F_i(X_i(s), s,t)ds+\varepsilon(t)}{Y(t)= \mu(t)+\int_{a_1}^{b_1}F_1(X_1(s), s,t)ds+...+\int_{a_p}^{b_p}F_p(X_p(s), s,t)ds+\epsilon(t),}
where \eqn{\mu(t)}{\mu(t)} is the intercept function, \eqn{\{F_i(x,
  s,t),1\le i\le p\}}{{F_i(x,s,t), 1\le i\le p}} are all unspecified
nonlinear functions of \eqn{x, s, t},  \eqn{\{X_i(s),1\le i\le p\}}{{X_i(s), 1\le i\le p}} are functional predictors, and \eqn{\epsilon(t)}{\epsilon(t)} is the noise function.

In this method, we require that all the sample curves of each functional predictor be observed in a common dense set, but the observation points can be different for different functional predictors. All the sample curves of the functional response are observed in a common dense set.

This method estimates a special decomposition:
\deqn{F_i(x,s,t)=\sum_{k=1}^\infty G_{ik}(x,s)w_k(t),  1\le i\le p.}{F_i(x,s,t)=G_{i1}(x,s)w_1(t)+...+G_{iK}(x,s)w_K(t),   1\le i\le p.}
We first estimate \eqn{\{G_{ik}(x,s), 1\le i\le p\}}{{G_{ik}(x,s), 1\le
    i\le p}} for each \eqn{k>0}  by solving a generalized functional
eigenvalue problem with smoothness tuning parameter \eqn{\lambda} on its
Sobolev norm, and then estimate \eqn{\{w_{k}(t), k>0\}}{{w_{k}(t), k>0}} by penalized least square with smoothness tuning parameter \eqn{\kappa}.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cv.nonlinear(X, Y, t.x.list, t.y, s.n.basis = 30, x.n.basis = 30,
        t.n.basis = 30, K.cv = 5, upper.com = 15, thresh = 0.005)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a list of length \eqn{p}, the number of functional predictors. Its \eqn{i}-th element is the \eqn{n\times m_i}{n*m_i} data matrix for the \eqn{i}-th functional predictor \eqn{X_i(s)}{X_i(s)}, where \eqn{n} is the sample size and \eqn{m_i} is the number of observation time points for  \eqn{X_i(s)}{X_i(s)}.
%%     ~~Describe \code{X} here~~
}
  \item{Y}{the \eqn{n\times m}{n*m} data matrix for the functional response \eqn{Y(t)}{Y(t)}, where  \eqn{n}  is the sample size and  \eqn{m}  is the number of the observation time points for \eqn{Y(t)}{Y(t)}.
%%     ~~Describe \code{Y} here~~
}
  \item{t.x.list}{a list of length \eqn{p}. Its \eqn{i}-th element is the vector of obesrvation time points of the \eqn{i}-th functional predictor \eqn{X_i(s)}{X_i(s)}, \eqn{1\le i\le p}{1\le i\le p}.
%%     ~~Describe \code{t.x} here~~
}
  \item{t.y}{the vector of obesrvation time points of the functional response \eqn{Y(t)}{Y(t)}.
%%     ~~Describe \code{t.y} here~~
}
  \item{s.n.basis}{the number of B-spline basis functions for the argument \eqn{s}. Default  is 30.
%%     ~~Describe \code{s.n.basis} here~~
}
  \item{x.n.basis}{the number of B-spline basis functions for \eqn{x}. Default is 30.
%%     ~~Describe \code{x.n.basis} here~~
}
  \item{t.n.basis}{the number of B-spline basis functions for \eqn{t}. Default is 30.
%%     ~~Describe \code{t.n.basis} here~~
}
  \item{K.cv}{the number of CV folds. Default is 5.
%%     ~~Describe \code{K.cv} here~~
}
  \item{upper.com}{the upper bound for the maximum number of components need to calculate. Default is 15.
%%     ~~Describe \code{upper.com} here~~
}
  \item{thresh}{a number between 0 and 1 specifying the minimum proportion of variation to be explained by each selected component relative to all the selected components. This will help determining the upper bound of the number of components to be used in CV, and the optimal number of components will be determined from 1,2,..., to the minimum of this uppber bound and upper.com by CV. A smaller thresh value leads to more components to be selected and a longer running time. A larger thresh value needs less running time, but may miss some important components and lead to a larger prediction error. Default is 0.005.
%%     ~~Describe \code{thresh} here~~
}
}
%\details{}

\value{A fitted CV-object, which is used in the function \code{\link{pred.nonlinear}} for prediction and \code{\link{getcoef.nonlinear}} for  extracting the estimated coefficient functions.
\item{opt.index}{index of the optimal lambda.}
\item{min.error}{minimum CV error.}
\item{errors}{list for CV errors. }
\item{lambda.set}{the set of \eqn{\lambda} values used in CV.}
\item{kappa.set}{the set of \eqn{\kappa} values used in CV.}
\item{opt.K}{optimal number of components to be selected.}
\item{opt.lambda}{optimal value for \eqn{\lambda}, the tuning parameter
  for \eqn{\{G_{ik}(x,s), 1 \le i \le p\}}{{G_{ik}(x,s), 1 \le i \le p}}.}
\item{opt.kappa}{optimal value for \eqn{\kappa}, the smoothness tuning parameter for \eqn{w_k(t)}{w_k(t)}.}
\item{fit.1}{a list for internal use.}
\item{y.params}{a list for internal use.}
}
\references{
   Xin Qi and Ruiyan Luo. (Accepted) Nonlinear function on function additive model
   with multiple predictor curves. Stistica Sinica.
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
#########################################
# toy example
# fit a nonlinear regression model with p=1 in the air quality data
#########################################

ptm <- proc.time()
data(air)
t.x=seq(0,1,length=24)
t.y=seq(0,1,length=24)
air.nonlinear.cv=cv.nonlinear(X=list(air[[2]][1:60,]), Y=air[[1]][1:60,], list(t.x),
           t.y, K.cv=2, s.n.basis = 8, x.n.basis = 10, t.n.basis = 8)
air.pred=pred.nonlinear(air.nonlinear.cv, list(air[[2]][1:2,]))
est.coefficient=getcoef.nonlinear(air.nonlinear.cv)
print(proc.time()-ptm)



\donttest{
#########################################
# example with simulated data and p=3
#########################################
rm(list=ls())
ptm <- proc.time()
library(MASS)

n.curves=3 # the number of predictor curves.
t.x=seq(0,1,length.out=80) # observation points for all six predictor curves.
t.y=seq(0,1,length.out=100) # observation points for response curve.
nnew <- 50
ntrain <-100
ntot <- nnew+ntrain

######################################################
#generate the predictor curves using sin() and cos() basis functions.
#########################################################
X=list()
X[[1]]=X[[2]]=X[[3]]=0
Sigma=matrix(0.5, 3, 3)
diag(Sigma)=1
for(k in 1:20){
  coef.mtx.1=mvrnorm(ntot, rep(0,3), Sigma)
  coef.mtx.2=mvrnorm(ntot, rep(0,3), Sigma)
  X[[1]]=X[[1]]+coef.mtx.1[,1]\%*\%t(sin(k*pi*t.x))/k+coef.mtx.2[,1]\%*\%t(cos(k*pi*t.x))/k
  X[[2]]=X[[2]]+coef.mtx.1[,2]\%*\%t(sin(k*pi*t.x))/k+coef.mtx.2[,2]\%*\%t(cos(k*pi*t.x))/k
  X[[3]]=X[[3]]+coef.mtx.1[,3]\%*\%t(sin(k*pi*t.x))/k+coef.mtx.2[,3]\%*\%t(cos(k*pi*t.x))/k
}

######################################################
#generate the nonlinear function and calculate the \int_0^1F_j(X_j(s),s,t)ds, j=1,..,3.
#########################################################
A.1=sapply(1:ntot, function(i){t.x})
A.1=t(A.1)
S.1.mtx=sapply(1:length(t.y), function(i){
         tmp= (A.1+X[[1]]-3*t.y[i])^2; apply(tmp,1,mean)})
S.2.mtx= sapply(1:length(t.y), function(i){
         tmp=cos(2*X[[2]]*A.1)+5*sin(X[[2]]*t.y[i]^(0.5)) ; apply(tmp,1,mean)})
S.3.mtx=sapply(1:length(t.y), function(i){
         tmp= 1/(1+X[[3]]^2+A.1*t.y[i]); apply(tmp,1,mean)})

######################################################
#generate the noise function, and response function
#########################################################
S.mtx=0.4*(S.1.mtx+S.2.mtx+S.3.mtx)
E<-matrix(rnorm(ntot*length(t.y)),ntot, length(t.y))
mu.val<- cos(4*pi*t.y)
Y<-t(sapply(1:ntot, function(i){mu.val+S.mtx[i,] + E[i,]}))
X.train=list()
X.new=list()
for(j in 1:n.curves)
{  X.train[[j]]=X[[j]][1:ntrain,]
   X.new[[j]]=X[[j]][-(1:ntrain),]
}
Y.train <- Y[1:ntrain, ]
Y.new <- Y[-(1:ntrain), ]
E.new <- E[-(1:ntrain),]

######################################################
#fit the nonlinear function-on-function model and make prediction
#########################################################
t.x.list=lapply(1:n.curves, function(k){t.x})
fit.cv=cv.nonlinear(X.train, Y.train, t.x.list, t.y)
Y.pred=pred.nonlinear(fit.cv, X.new)
nonlinear.error=mean((Y.pred-Y.new)^2)
print(c("nonlinear prediction error=",nonlinear.error))
nonlinear.esterror=mean((Y.pred-Y.new+E.new)^2)
print(c("nonlinear estimation.error=",nonlinear.esterror))

######################################################
#extract the estimated mu(t), F(x,s,t)
#########################################################
F.list=getcoef.nonlinear(fit.cv)
mu=F.list$mu
F=F.list$F
X.grid=F.list$X.grid
print(proc.time()-ptm)
#   user  system elapsed
#  59.30    7.63   66.94
}
}