\name{fgam}
\Rdversion{1.0}
\alias{fgam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functional Generalized Additive Models}
\description{
Implements functional generalized additive models for functional and scalar covariates and scalar responses. Additionally implements functional linear models.  This function is a wrapper for mgcv's gam and its siblings to fit models of the general form 
\deqn{g(E(Y_i)) = \beta_0 + \int_{T_1} F(X_{i1},t)dt+ \int_{T_2} \beta(t)X_{i2}dt + f(z_{i1}) + f(z_{i2}, z_{i3}) + \ldots)}
with a scalar (but not necessarily continuous) response Y, and link function g}
\usage{
fgam(formula,fitter=NA,tensortype=c('te','t2'),...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
a formula with special terms as for gam, with additional special terms \code{\link{af}}() and \code{\link{lf}}().}
  \item{fitter}{
the name of the function used to estimate the model. Defaults to \code{\link{gam}} if the matrix of functional responses has less than 2e5 data points and to \code{\link{bam}} if not. "gamm" (see \code{\link{gamm}}) and "gamm4" (see \code{\link{gamm4}}) are valid options as well.}
\item{tensortype}{which type of tenor product splines to use. One of "te" or "t2", defaults to te}
  \item{...}{
   additional arguments that are valid for gam or bam; for example, specify a \code{gamma} > 1 to increase amount of smoothing when using GCV to choose smoothing parameters or \code{method="REML"} to change to REML for estimation of smoothing parameters (default is GCV).}
}

\section{Warning}{
Binomial responses should be specified as a numeric vector rather than as a matrix or a factor.}

\value{
a fitted fgam-object, which is a \code{\link{gam}}-object with some additional information in a fgam-entry. If fitter is "gamm" or "gamm4", only the $gam part of the returned list is modified in this way..}
\references{
McLean, M. W., Hooker, G., Staicu, A.-M., Scheipl, F., and Ruppert, D. (2013). Functional generalized additive models. \emph{Journal of Computational and Graphical Statistics}, to appear. Available at \url{http://courses2.cit.cornell.edu/mwmclean}}

\author{
Mathew W. McLean \email{mwm79@cornell.edu} and Fabian Scheipl}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{af}}, \code{\link{lf}}, \code{\link{predict.fgam}}, \code{\link{vis.fgam}}}

\examples{
################# DTI Example #####################
data(DTI)
## only consider first visit and cases (no PASAT scores for controls)
y <- DTI$pasat[DTI$visit==1 & DTI$case==1]
X <- DTI$cca[DTI$visit==1 & DTI$case==1,]
X_2 <- DTI$rcst[DTI$visit==1 & DTI$case==1,]

## remove samples containing missing data
ind <- rowSums(is.na(X))>0
ind2 <- rowSums(is.na(X_2))>0

y <- y[!(ind | ind2)]
X <- X[!(ind | ind2),]
X_2 <- X_2[!(ind | ind2),]

N <- length(y)

## fit fgam using FA measurements along corpus callosum 
## as functional predictor with PASAT as response
## using 8 cubic B-splines for marginal bases with third 
## order marginal difference penalties
## specifying gamma>1 enforces more smoothing when using 
## GCV to choose smoothing parameters
fit <- fgam(y~af(X,splinepars=list(k=c(8,8),m=list(c(2,3),c(2,3)))),gamma=1.2)
vis.fgam(fit)


## fgam term for the cca measurements plus an flm term for the rcst measurements
## leave out 10 samples for prediction
test <- sample(N,10)
fit <- fgam(y~af(X,splinepars=list(k=c(7,7),m=list(c(2,2),c(2,2))))+
  lf(X_2,splinepars=list(k=7,m=c(2,2))),subset=(1:N)[-test])
## plot the fits
plot(fit)
## vis.fgam(fit,af.term='X')
vis.fgam(fit,af.term='X',xval=.6)
## predict the ten left outs samples
pred <- predict(fit,newdata=list(X=X[test,],X_2=X_2[test,]),type='response',PredOutOfRange=TRUE)
sqrt(mean((y[test]-pred)^2))

## Try to predict the binary response disease status (case or control)
##   using the quantile transformed measurements from the rcst tract
##   with a smooth component for a scalar covariate that is pure noise
y <- DTI$case[DTI$visit==1]
X <- DTI$cca[DTI$visit==1,]
X_2 <- DTI$rcst[DTI$visit==1,]

ind <- rowSums(is.na(X))>0
ind2 <- rowSums(is.na(X_2))>0

y <- y[!(ind | ind2)]
X <- X[!(ind | ind2),]
X_2 <- X_2[!(ind | ind2),]
z1 <- rnorm(length(y))

# select=TRUE allows terms to be zeroed out of model completely
fit <- fgam(y~s(z1,k=10)+af(X_2,splinepars =list(k=c(7,7),m=list(c(2,1),c(2,1))),Qtransform=TRUE),
            family=binomial(),select=TRUE)
plot(fit)
vis.fgam(fit,af.term='X_2',plot.type='contour')
}
