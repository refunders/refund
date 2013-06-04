\name{test.pfr}
\alias{test.pfr}
\title{
Likelihood Ratio Test and Restricted Likelihood Ratio Test for inference of functional predictors
}
\description{
Given a pfr object of family="gaussian", conducts either a test to see if the function is
just a mean (constancy) or a test to see if the functional predictor
significantly improves the model (inclusion).  Based on work 0-variance
component work of Crainiceanu et al. (2004) and Scheipl et al. (2008)
and Swihart et al. (2012).}
\usage{
test.pfr(pfr.obj=pfr.obj, test=NULL, ...)
}
\arguments{
  \item{pfr.obj}{an object returned by pfr}
  \item{test }{"constancy" will test functional form of the coefficient
    function of the last function listed in funcs in pfr.obj against the
    null of a constant line: the average of the functional predictor.
    "inclusion"  will test functional form of the coefficient
    function of the last function listed in funcs in pfr.obj against the
    null of 0: that is, whether the functional predictor should be
    included in the model.}
  \item{...}{additional arguments}
}
\details{
A Penalized Functional Regression of family="gaussian" can be
represented as a linear mixed model dependent on variance components.
Testing whether certain variance components and (potentially) fixed
effect coefficients are 0 correspond to tests of constancy and inclusion
of functional predictors.

For test.pfr, Restricted Likelihood Ratio Test is preferred for the
constancy test as under the special B-splines implementation of pfr for
the coefficient function basis the test involves only the variance
component.  Therefore, the constancy test is best for pfr objects with
method="REML"; if the method was something else, a warning is printed
and the model refit with "REML" and a test is then conducted.

For test.pfr, the Likelihood Ratio Test is preferred for the
inclusion test as under the special B-splines implementation of pfr for
the coefficient function basis the test involves both the variance
component and a fixed effect coefficient in the linear mixed model
representation.
Therefore, the inclusion test is best for pfr objects with
method="ML"; if the method was something else, a warning is printed
and the model refit with "ML" and a test is then conducted.
}

\value{
  \item{p.val }{the p-value for the full model (alternative) against the
    null specified by the test}
  \item{test.stat }{the test statistic, see Scheipl et al. 2008 and
    Swihart et al 2012}

  \item{ma }{the alternative model as fit with mgcv::gam}
  \item{m0 }{the null model as fit with mgcv::gam}
  \item{m }{the model containing only the parameters being tested as fit
  with mgcv::gam}

}

\references{
Goldsmith, J., Bobb, J., Crainiceanu, C., Caffo, B., and Reich, D. (2011).
Penalized functional regression. \emph{Journal of Computational and Graphical Statistics}, 20(4), 830-851.

Goldsmith, J., Crainiceanu, C., Caffo, B., and Reich, D. (2012).
Longitudinal penalized functional regression for cognitive outcomes on
neuronal tract measurements. \emph{Journal of the Royal Statistical Society: Series C}, 61(3), 453--469.

Crainiceanu, C. and Ruppert, D. (2004) Likelihood ratio tests in linear
mixed models with one variance component. \emph{Journal of the Royal
Statistical Society: Series B},66,165--185.

Scheipl, F. (2007) Testing for nonparametric terms and random effects in structured additive regression. Diploma thesis.\ http://www.statistik.lmu.de/~scheipl/downloads/DIPLOM.zip.

Scheipl, F., Greven, S. and Kuechenhoff, H (2008)
Size and power of tests for a zero random effect variance or polynomial regression in additive and linear mixed models.             \emph{Computational Statistics & Data Analysis}, 52(7):3283-3299    

Swihart, Bruce J., Goldsmith, Jeff; and Crainiceanu, Ciprian M. (2012).
Testing For Functional Effects 
\emph{Johns Hopkins University, Dept. of Biostatistics Working Papers. Working Paper 247.}
http://biostats.bepress.com/jhubiostat/paper247
}

\seealso{\code{\link{pfr}}, \code{\link{predict.pfr}}, package \code{RLRsim}}
\author{
  Jeff Goldsmith <jeff.goldsmith@columbia.edu>
  Bruce Swihart <bswihart@jhsph.edu>
}
\examples{
##################################################################
#########               DTI Data Example                 #########
##################################################################

##################################################################
# For more about this example, see Swihart et al. 2012
# Testing for Functional Effects
##################################################################

## load and reassign the data;
data(DTI2)
Y  <- DTI2$pasat ## PASAT outcome
id <- DTI2$id    ## subject id
W1 <- DTI2$cca   ## Corpus Callosum
W2 <- DTI2$rcst  ## Right corticospinal
V  <- DTI2$visit ## visit

## prep scalar covariate
visit.1.rest <- matrix(as.numeric(V > 1), ncol=1)
covar.in <- visit.1.rest 


## note there is missingness in the functional predictors
apply(is.na(W1), 2, mean)
apply(is.na(W2), 2, mean)

## fit two univariate models, then one model with both functional predictors
pfr.obj.t1 <- pfr(Y = Y, covariates=covar.in, funcs = list(W1),     subj = id, kz = 10, kb = 50)
pfr.obj.t2 <- pfr(Y = Y, covariates=covar.in, funcs = list(W2),     subj = id, kz = 10, kb = 50)
pfr.obj.t3 <- pfr(Y = Y, covariates=covar.in, funcs = list(W1, W2), subj = id, kz = 10, kb = 50)

## plot the coefficient function and bounds
dev.new()
par(mfrow=c(2,2))
ran <- c(-2,.5)
matplot(cbind(pfr.obj.t1$BetaHat[[1]], pfr.obj.t1$Bounds[[1]]),
  type = 'l', lty = c(1,2,2), col = c(1,2,2), ylab = "BetaHat", 
  main = "CCA", xlab="Location", ylim=ran)
abline(h=0, col="blue")
matplot(cbind(pfr.obj.t2$BetaHat[[1]], pfr.obj.t2$Bounds[[1]]),
  type = 'l', lty = c(1,2,2), col = c(1,2,2), ylab = "BetaHat", 
  main = "RCST", xlab="Location", ylim=ran)
abline(h=0, col="blue")
matplot(cbind(pfr.obj.t3$BetaHat[[1]], pfr.obj.t3$Bounds[[1]]),
  type = 'l', lty = c(1,2,2), col = c(1,2,2), ylab = "BetaHat", 
  main = "CCA  - mult.", xlab="Location", ylim=ran)
abline(h=0, col="blue")
matplot(cbind(pfr.obj.t3$BetaHat[[2]], pfr.obj.t3$Bounds[[2]]),
  type = 'l', lty = c(1,2,2), col = c(1,2,2), ylab = "BetaHat", 
  main = "RCST - mult.", xlab="Location", ylim=ran)
abline(h=0, col="blue")

## do some testing
t1 <- test.pfr(pfr.obj.t1, "constancy")
t2 <- test.pfr(pfr.obj.t2, "constancy")
t3 <- test.pfr(pfr.obj.t3, "inclusion")

\dontrun{

t1$test.stat
t1$p.val

t2$test.stat
t2$p.val

t3$test.stat
t3$p.val


## do some testing with test.pfr(); same as above but subj = NULL
pfr.obj.t1 <- pfr(Y = Y, covariates=covar.in, funcs = list(W1),     subj = NULL, kz = 10, kb = 50)
pfr.obj.t2 <- pfr(Y = Y, covariates=covar.in, funcs = list(W2),     subj = NULL, kz = 10, kb = 50)
pfr.obj.t3 <- pfr(Y = Y, covariates=covar.in, funcs = list(W1, W2), subj = NULL, kz = 10, kb = 50)

t1 <- test.pfr(pfr.obj.t1, "constancy")
t2 <- test.pfr(pfr.obj.t2, "constancy")
t3 <- test.pfr(pfr.obj.t3, "inclusion")

t1$test.stat
t1$p.val

t2$test.stat
t2$p.val

t3$test.stat
t3$p.val
}
}