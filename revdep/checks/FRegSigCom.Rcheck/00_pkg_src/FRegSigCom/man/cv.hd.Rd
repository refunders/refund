\name{cv.hd}
\alias{cv.hd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Cross-validation for sparse linear function-on-function regression with a large number of functional predictors
%%  ~~function to do ... ~~
}
\description{Conduct cross-validation and build the final model for the following function-on-function regression model:
\deqn{Y(t)= \mu(t)+\sum_{i=1}^p\int_{a_i}^{b_i}X_i(s)\beta_i(s,t)ds+\varepsilon(t)}{Y(t)= \mu(t)+\int X_1(s)\beta_1(s,t)ds+...+\int X_p(s)\beta_p(s,t)ds+\varepsilon(t),}
where \eqn{\mu(t)}{\mu(t)} is the intercept
function. The \eqn{\{X_i(s),1\le i\le p\}}{{X_i(s), 1\le i\le p}} are
\eqn{p} functional predictors and \eqn{\{\beta_i(s,t),1\le i\le p\}}{{\beta_i(s,t),1\le i\le p}}
 are corresponding coefficient functions.
The \eqn{\epsilon(t)}{\epsilon(t)} is the noise function. The \eqn{p} can be much larger than the sample size.

We require that all the sample curves of each functional predictor are observed in a common dense set of time points, but the set can be different for different functional predictor. All the sample curves of the response are observed in a common dense set.

This method estimates a special decomposition of the coefficient functions induced by the KL expansion of the signal function:
\deqn{\beta_i(s,t)=\sum_{k=1}^\infty\psi_{ik}(s)w_k(t),  1\le i\le p.}{\beta_i(s,t)=\psi_{i1}(s)w_1(t)+\psi_{i2}(s)w_2(t)+...,   1\le i\le p.}
We first estimate \eqn{\{\psi_{ik}(s), 1\le i\le p\}}{{\psi_{ik}(s),
    1\le i\le p}} for each \eqn{0<k<K}{0<k<K}, with a simultaneous sparse and smooth penalty imposed,
\deqn{P(\psi_1(s),...,\psi_p(s))=\tau \{(1-\lambda) \sum_{i=1}^p||\psi_{i}||_\eta^2 + \lambda (\sum_{i=1}^p ||\psi_{i}||_\eta)^2\},}{P(\psi_1(s),...,\psi_p(s))=\tau {(1-\lambda) (||\psi_{1}||_\eta^2+...+||\psi_{p}||_\eta^2)+\lambda(||\psi_{1}||_\eta+||\psi_{p}||_\eta)^2},}
where \eqn{||\psi_{i}||_\eta^2=||\psi_i||^2+\eta||\psi_i''||^2}.
Then we estimate \eqn{\{w_{k}(t), 0<k<K_0\}}{{w_{k}(t), 0<k<K}} by a
penalized least square problem with the smoothness of \eqn{w_k(t)} tuned by \eqn{\kappa}, where \eqn{K} is the number of componets to be selected and will be chosen by cross-validation. The simultaneous sparse and smooth penalty will make most estimates of  \eqn{\psi_{ik}(s)}{\psi_{ik}(s)}'s, and hence \eqn{\beta_i(s,t)}{\beta_i(s,t)}'s, exactly equal to zero.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}

\usage{
cv.hd(X, Y, t.x.list, t.y, K.cv = 5, s.n.basis = 25, t.n.basis = 50,
       thresh = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a list of length \eqn{p}. Its \eqn{i}-th element is the \eqn{n\times m_i}{n*m_i} data matrix for the \eqn{i}-th funtional predictor \eqn{X_i(s)}{X_i(s)}, where \eqn{n} is the sample size and \eqn{m_i} is the number of the observation points for  \eqn{X_i(s)}{X_i(s)}.
%%     ~~Describe \code{X} here~~
}
  \item{Y}{the \eqn{n\times m}{n*m} data matrix for the functional response \eqn{Y(t)}{Y(t)}, where \eqn{n}{n} is the sample size and \eqn{m}{m} is the number of the observation points for \eqn{Y(t)}{Y(t)}.
%%     ~~Describe \code{Y} here~~
}
  \item{t.x.list}{a list of length \eqn{p}. Its \eqn{i}-th element is the vector of observation points of the \eqn{i}-th functional predictor \eqn{X_i(s)}{X_i(s)}, \eqn{1\le i\le p}{1\le i\le p}.
%%     ~~Describe \code{t.x.list} here~~
}
  \item{t.y}{the vector of obesrvation points of the functional response \eqn{Y(t)}{Y(t)}.
%%     ~~Describe \code{t.y} here~~
}
  \item{K.cv}{the number of CV folds. Default is 5.
%%     ~~Describe \code{K.cv} here~~
}
  \item{s.n.basis}{the number of B-spline basis functions used for estimating the functions \eqn{\psi_{ik}(s)}{\psi_{ik}(s)'s}. Default  is 25.
%%     ~~Describe \code{s.n.basis} here~~
}
  \item{t.n.basis}{the number of B-spline basis functions used for estimating the functions \eqn{w_{k}(t)}{w_{k}(t)'s}. Default is 50.
%%     ~~Describe \code{t.n.basis} here~~
}
  \item{thresh}{a number between 0 and 1 specifying the minimum proportion of variation to be explained by each selected component relative to all the selected components. This will determine the upper bound of the number of components for CV, and the optimal number of components will be determined from 1,2,..., to this uppber bound by CV. A smaller thresh value leads to more components to be selected and a longer running time. A larger thresh value needs less running time, but may miss some important components and lead to a larger prediction error. Default is 0.01.
  }
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{An object of the ``cv.hd'' class, which is used in the function \code{\link{pred.hd}} for prediction and \code{\link{getcoef.hd}} for extracting the estimated coefficient functions.
\item{errors}{list for CV errors. }
\item{min.error}{minimum CV error.}
\item{opt.index}{index of the optimal tunning parameters.}
\item{params.set}{the set of candidate values for tuning parameters used in CV. It's a matrix with 4 columns  corresponding to the tuning parameters \eqn{\tau}, \eqn{\lambda}, \eqn{\eta}, \eqn{\kappa}, respectively.}
\item{opt.K}{optimal number of components to be selected.}
\item{opt.tau}{optimal value for \eqn{\tau} tuning the simultaenous
  sparse-smooth penalty on \eqn{\{\psi_{ik}(s), 1\le i \le p
    \}}{{\psi_{ik}(s), 1 \le i \le p}}.}
\item{opt.lambda}{optimal value for \eqn{\lambda}, the sparsity parameter for   \eqn{\{\psi_{ik}(s), 1\le i \le p
    \}}{{\psi_{ik}(s), 1 \le i \le p}}.}
\item{opt.eta}{optimal value for \eqn{\eta}, the smoothness parameter for \eqn{\psi_{ik}(s)}{\psi_{ik}(s)}.}
\item{opt.kappa}{optimal value for \eqn{\kappa}, the smoothness parameter for \eqn{w_k(t)}{w_k(t)}.}
\item{maxK.ret}{a list for internal use.}
\item{t.y}{the input \code{t.y}.}
\item{y.params}{a list for internal use.}
}
\references{Xin Qi and Ruiyan Luo (2018) Function-on-Function
  Regression with thousands of predictive curves, Jounal of Multivariate
  Analysis 163:51-66. \url{https://doi.org/10.1016/j.jmva.2017.10.002}}
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
#toy example using the air quality data with p=1
#########################################
data(air)
t.x=seq(0,1,length=24)
t.y=seq(0,1,length=24)
air.cv=cv.hd(X=list(air[[2]][1:20,]), Y=air[[1]][1:20,], list(t.x), t.y,
             K.cv=2, s.n.basis = 8, t.n.basis = 8)
air.pred=pred.hd(air.cv, list(air[[2]][1:2,]))
predict.error=mean((air.pred-air[[1]][1:2,])^2)
print(c("predict error=", predict.error))


\donttest{
#########################################
#example with simulated data and p=200
#########################################

rm(list=ls())
ptm <- proc.time()
library(MASS)

n.curves=200 # total number of predictor curves.

t.x=seq(0,1,length.out=80) # observation points for X_i(s).
t.y=seq(0,1,length.out=100) # observation points for Y(t).
ntrain <- 100 # number of observations in training data
nnew <-50 # number of new observations
ntot <- ntrain+nnew
########################################
##generate the predictor curves using cos() basis functions.
############################################
W <- list()
for(j in 1:(n.curves+5))
{
  W[[j]] <- 0
  for(k in 1:30)
  W[[j]] <- W[[j]]+rnorm(ntot) \%*\% t(cos(k*pi*t.x))/k
}
X=lapply(1:n.curves,
function(k){(W[[k]]+W[[k+1]]+W[[k+2]]+W[[k+3]]+W[[k+4]]
             +W[[k+5]])/k})
########################################
##generate the coefficient functions. Only the first five are nonzero.
############################################
mu=sin(2*pi*t.y)
B.coef=list()
B.coef[[1]]=sapply(1:length(t.y), function(k){4*t.x^2*t.y[k]})
B.coef[[2]]=sapply(1:length(t.y), function(k){4*sin(pi*t.x)*cos(pi*t.y[k])})
B.coef[[3]]=sapply(1:length(t.y), function(k){4*exp(-((t.x-0.5)^2+(t.y[k]-0.5)^2)/2)})
B.coef[[4]]=sapply(1:length(t.y), function(k){4*t.x/(1+t.y[k]^2)})
B.coef[[5]]=sapply(1:length(t.y), function(k){4*log(1+t.x)*sqrt(t.y[k])})
########################################
##generate the response curves
############################################
Y=0
for(j in 1:5){
Y<- Y+ (X[[j]] \%*\% B.coef[[j]])/length(t.x)
}
E=sqrt(0.01)*matrix(rnorm(ntot*length(t.y)), ntot, length(t.y))
Y=rep(1,ntot)\%*\%t(mu)+Y+E
X.train=lapply(1:n.curves, function(k){X[[k]][(1:ntrain),]})
X.new=lapply(1:n.curves, function(k){X[[k]][-(1:ntrain),]})
Y.train <- Y[1:ntrain,]
Y.new <- Y[-(1:ntrain),]
E.new=E[-(1:ntrain),]

########################################
##fit the model using the sparse function-on-function method
############################################
t.x.list=lapply(1:n.curves, function(k){t.x})
fit.cv <- cv.hd(X.train, Y.train, t.x.list, t.y,   K.cv=5, s.n.basis=25, t.n.basis=40)
Y.pred=pred.hd(fit.cv, X.new)
predict.error=mean((Y.pred-Y.new)^2)
est.error=mean((Y.pred-Y.new+E.new)^2)
print(c("predict error=", predict.error))
print(c("estimation error=", est.error))
est.coefficient=getcoef.hd(fit.cv)
mu.est=est.coefficient[[1]]
beta.est=est.coefficient[[2]]
}
}
