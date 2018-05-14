\name{cv.ff.interaction}
\alias{cv.ff.interaction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Cross-validation for function-on-function regression models with specified main effects and two-way
  interaction terms
%%  ~~function to do ... ~~
}
\description{This function is used to perform cross-validation and build
  the final model using the signal compression approach for the
  following linear function-on-function regression model with main
  effects, two-way interaction between them, and quadratic terms. Let
  \eqn{\{X_i(s), 1\le i\le p\}}{{X_i(s), 1\le i \le p}} be \eqn{p} potential functional predictors. The model is given by
\deqn{Y(t)= \mu(t)+\sum_{i \in M}\int_{a_i}^{b_i} X_i(s)\beta_i(s,t)ds+\sum_{(i,j) \in I}\int_{a_i}^{b_i}\int_{a_j}^{b_j} X_i(u)X_j(v)\gamma_{ij}(u,v,t)dudv+\epsilon(t),}{Y(t)= \mu(t)+\sum_{i \in M}\int_{a_i}^{b_i} X_i(s)\beta_i(s,t)ds+\sum_{(i,j) \in I}\int_{a_i}^{b_i}\int_{a_j}^{b_j} X_i(u)X_j(v)\gamma_{ij}(u,v,t)dudv+\epsilon(t),}
where \eqn{\mu(t)}{\mu(t)} is the intercept function. The index set
\eqn{M} of main effects is a subset of \eqn{\{1,...,p\}}{{1,...,p}}, and the index
set \eqn{I} of interactions and quadratic effects is a subset of the
collection of all possible pairs \eqn{\{(i,j), 1\le i\le j\le
    p\}}{{(i,j), 1\le i\le j\le p}}. The \eqn{\{\beta_i(s,t), i \in
  M\}}{{\beta_i(s,t), i\in M}} and  \eqn{\{\gamma_{ij}(u,v,t),(i,j)\in I\}}{{\gamma_{ij}(u,v,t),(i,j)\in I}} are the corresponding coefficient functions. The \eqn{\epsilon(t)}{\epsilon(t)} is the noise function.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cv.ff.interaction(X, Y, t.x, t.y, main.effect, interaction.effect=NULL,
       s.n.basis=40, t.n.basis=40, inter.n.basis=20, basis.type="Bspline",
       K.cv=5, all.folds=NULL, upper.comp=8, thresh=0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a list of \eqn{p} potential functional predictors. Its \eqn{i}-th element is the \eqn{n\times m_i}{n*m_i} data matrix for the \eqn{i}-th potential functional predictor \eqn{X_i(s)}{X_i(s)}, where \eqn{n} is the sample size and \eqn{m_i} is the number of observation time points for  \eqn{X_i(s)}{X_i(s)}.
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
  \item{main.effect}{a vector of indices for main effects. It
    is a subset of \eqn{\{1,2,...,p\}}{{1,2,...,p}}.
}
  \item{interaction.effect}{a matrix of two columns.  Each row of this matrix specifies the index of a two-way interaction or a quadratic effect. Default is NULL.
}
  \item{s.n.basis}{the number of basis functions used for estimating the functions \eqn{\psi_{ik}(s)}{\psi_{ik}(s)'s} (see the reference for details). Default is 40.
%%     ~~Describe \code{s.n.basis} here~~
}
  \item{t.n.basis}{the number of basis functions used for estimating the functions \eqn{w_{k}(t)}{w_{k}(t)'s}. Default is 40.
%%     ~~Describe \code{t.n.basis} here~~
}
  \item{inter.n.basis}{the number of one-dimensional basis functions
    used to construct the tensor product basis functions for estimating the functions
    \eqn{\phi_{ijk}(u,v)}{\phi_{ijk}(u,v)'s}. Default is 20.
%%     ~~Describe \code{t.n.basis} here~~
}
   \item{basis.type}{the type of basis functions. Only "BSpline" (default) and "Fourier" are supported.
}
  \item{K.cv}{the number of CV folds. Default is 5.
%%     ~~Describe \code{K.cv} here~~
}
  \item{all.folds}{a list of cross validation sets. The length of the
    list is equal to \code{K.cv}. The \eqn{i}-th element of the list is
    the collection of indices of observations assigned to the \eqn{i}-th
    cross validation set. Default is NULL, and the cross validation
    sets are generated in function cv.ff.interaction.
%%     ~~Describe \code{K.cv} here~~
}
  \item{upper.comp}{the upper bound for the maximum number of components to be calculated. Default is 8.
%%     ~~Describe \code{upper.comp} here~~
}
  \item{thresh}{a number between 0 and 1 specifying the minimum proportion of variation to be explained by each selected component relative to all the selected components. This will help determining the upper bound of the number of components to be used in CV, and the optimal number of components will be determined from 1,2,..., to the minimum of this uppber bound and upper.com by CV. A smaller thresh value leads to more components to be selected and a longer running time. A larger thresh value needs less running time, but may miss some important components and lead to a larger prediction error. Default is 0.01.
%%     ~~Describe \code{thresh} here~~
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{An object of the ``cv.ff.interaction'' class, which is used in the function \code{\link{pred.ff.interaction}} for prediction and \code{\link{getcoef.ff.interaction}} for extracting the estimated coefficient functions.
\item{fitted_model}{a list for interval use.}
\item{y_penalty_inv}{a list for internal use.}
\item{X}{the input \code{X}.}
\item{Y}{the input \code{Y}.}
\item{x.smooth.params}{a list for internal use.}
\item{y.smooth.params}{a list for internal use.}
}
%\value{Y}
\references{Ruiyan Luo and Xin Qi (under revision) Interaction model and model selection for function-on-function regression.
%% ~put references to the literature/web site here ~
}
\author{Xin Qi and Ruiyan Luo}
\seealso{ \code{\link{pred.ff.interaction}}, \code{\link{getcoef.ff.interaction}}}
\examples{

library(FRegSigCom)
data(ocean)

Y=ocean$Salinity
X=list()
X[[1]]=ocean$Potential.density
X[[2]]=ocean$Temperature
X[[3]]=ocean$Oxygen
X[[4]]=ocean$Chloropigment
n.curves=length(X)
ntot=dim(Y)[1]
ntrain=50
ntest=ntot-ntrain
X.uncent=X
for(i in 1:n.curves){
  X[[i]]=scale(X.uncent[[i]],center=TRUE, scale=FALSE)
}
lengthX=dim(X[[1]])[2]
lengthY=dim(Y)[2]
t.x=seq(0,1,length=lengthX)
t.y=seq(0,1,length=lengthY)
I.train=sample(1:ntot, ntrain)
X.train=list()
X.test=list()
t.x.all=list()
for(i in 1:n.curves)
{
  X.train[[i]]=X[[i]][I.train,]
  X.test[[i]]=X[[i]][-I.train,]
  t.x.all[[i]]=t.x
}
Y.train=Y[I.train, ]
Y.test=Y[-I.train, ]


#########################################################################
# an interaction model with given main effects and two-way interactions
#########################################################################
main.effect=c(1,2,3)
inter.effect=rbind(c(1,1),c(1,2))
fit.inter=cv.ff.interaction(X.train, Y.train, t.x.all, t.y,  main.effect, inter.effect)
Y.pred=pred.ff.interaction(fit.inter,  X.test)
error.inter=mean((Y.pred-Y.test)^2)
print(c("error.inter=", error.inter))
#coef.obj=getcoef.ff.interaction(fit.inter)
#str(coef.obj)
}
