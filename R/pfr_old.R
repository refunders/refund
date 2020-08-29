#' Penalized Functional Regression (old version)
#'
#' This code implements the function pfr() available in refund 0.1-11. It is included
#' to maintain backwards compatibility. 
#'   
#' Functional predictors are entered as a matrix or, in the case of
#' multiple functional predictors, as a list of matrices using the
#' \code{funcs} argument. Missing values are allowed in the functional
#' predictors, but it is assumed that they are observed over the same
#' grid. Functional coefficients and confidence bounds are returned as
#' lists in the same order as provided in the \code{funcs} argument, as
#' are principal component and spline bases.  Increasing values of
#' \code{nbasis} will increase computational time and the values of
#' \code{nbasis}, \code{kz}, and \code{kb} in relation to each other may
#' need to be adjusted in application-specific ways.
#'
#' @param Y vector of all outcomes over all visits
#' @param subj vector containing the subject number for each observation
#' @param covariates matrix of scalar covariates
#' @param funcs matrix, or list of matrices, containing observed functional 
#'    predictors as rows. NA values are allowed.
#' @param kz can be NULL; can be a scalar, in which case this will be the 
#'    dimension of principal components basis for each and every observed 
#'    functional predictors; can be a vector of length equal to the number 
#'    of functional predictors, in which case each element will correspond 
#'    to the dimension of principal components basis for the corresponding 
#'    observed functional predictors
#' @param kb dimension of the B-spline basis for the coefficient function 
#'    (note: this is a change from versions 0.1-7 and previous) 
#' @param nbasis passed to refund::fpca.sc (note: using fpca.sc is a change 
#'    from versions 0.1-7 and previous)
#' @param family generalized linear model family
#' @param method method for estimating the smoothing parameters; defaults 
#'    to REML
#' @param smooth.option method to do FPC decomposition on the predictors. 
#'    Two options available -- "fpca.sc" or "fpca.face". If using "fpca.sc", 
#'    a number less than 35 for \code{nbasis} should be used while if using 
#'    "fpca.face",35 or more is recommended. 
#' @param pve proportion of variance explained used to choose the number of 
#'    principal components to be included in the expansion.
#' @param ... additional arguments passed to \code{\link[mgcv]{gam}} to 
#'    fit the regression model.
#'
#' @section Warning:
#'   Binomial responses should be specified as a numeric vector rather than as a
#'   matrix or a factor.
#' @return
#'  \item{fit }{result of the call to \code{gam}}
#' \item{fitted.vals }{predicted outcomes}
#' \item{fitted.vals.level.0 }{predicted outcomes at population level}
#' \item{fitted.vals.level.1 }{predicted outcomes at subject-specific level (if applicable)}
#' \item{betaHat }{list of estimated coefficient functions}
#' \item{beta.covariates }{parameter estimates for scalar covariates}
#' \item{varBetaHat }{list containing covariance matrices for the estimated coefficient functions}
#' \item{Bounds }{list of bounds of a pointwise 95\% confidence interval for the estimated coefficient functions}
#' \item{X }{design matrix used in the model fit}
#' \item{D }{penalty matrix used in the model fit}
#' \item{phi }{list of B-spline bases for the coefficient functions}
#' \item{psi }{list of principal components basis for the functional predictors}
#' \item{C }{stacked row-specific principal component scores}
#' \item{J }{transpose of psi matrix multiplied by phi}
#' \item{CJ }{C matrix multiplied J}
#' \item{Z1 }{design matrix of random intercepts}
#' \item{subj }{subject identifiers as specified by user}
#' \item{fixed.mat }{the fixed effects design matrix of the pfr as a mixed model}
#' \item{rand.mat }{the fixed effects design matrix of the pfr as a mixed model}
#' \item{N_subj }{the number of unique subjects, if subj is specified}
#' \item{p }{number of scalar covariates}
#' \item{N.pred }{number of functional covariates}
#' \item{kz }{as specified}
#' \item{kz.adj}{For smooth.option="fpca.sc", will be same as kz (or a vector of repeated values of the specified scalar kz).  For smooth.option="fpca.face", will be the corresponding number of principal components for each functional predictor as determined by fpca.face; will be less than or equal to kz on an elemental-wise level.}
#' \item{kb }{as specified}
#' \item{nbasis }{as specified}
#' \item{totD }{number of penalty matrices created for mgcv::gam}
#' \item{funcs }{as specified}
#' \item{covariates }{as specified}
#' \item{smooth.option}{as specified}
#' 
#' @references
#' Goldsmith, J., Bobb, J., Crainiceanu, C., Caffo, B., and Reich, D. (2011).
#' Penalized functional regression. \emph{Journal of Computational and Graphical
#' Statistics}, 20(4), 830-851.
#'
#' Goldsmith, J., Crainiceanu, C., Caffo, B., and Reich, D. (2012). Longitudinal
#' penalized functional regression for cognitive outcomes on neuronal tract
#' measurements. \emph{Journal of the Royal Statistical Society: Series C},
#' 61(3), 453-469.
#'
#' Swihart, Bruce J., Goldsmith, Jeff; and Crainiceanu, Ciprian M. (July 2012). 
#' Testing for functional effects. Johns Hopkins University Dept. of Biostatistics 
#' Working Paper 247, available at \url{https://biostats.bepress.com/jhubiostat/paper247/}
#' American Statistical Association, 109(508): 1425-1439.
#'
#' @author Bruce Swihart \email{bruce.swihart@@gmail.com} and 
#' Jeff Goldsmith \email{jeff.goldsmith@@columbia.edu}
#' 
#' @seealso \code{\link{rlrt.pfr}}, \code{\link{predict.pfr}}.
#' @importFrom mgcv gam
#' @export
#'
#' @examples
#' \dontrun{
#' ##################################################################
#' #########               DTI Data Example                 #########
#' ##################################################################
#' 
#' ##################################################################
#' # For more about this example, see Swihart et al. 2013
#' ##################################################################
#' 
#' ## load and reassign the data;
#' data(DTI2)
#' Y  <- DTI2$pasat ## PASAT outcome
#' id <- DTI2$id    ## subject id
#' W1 <- DTI2$cca   ## Corpus Callosum
#' W2 <- DTI2$rcst  ## Right corticospinal
#' V  <- DTI2$visit ## visit
#' 
#' ## prep scalar covariate
#' visit.1.rest <- matrix(as.numeric(V > 1), ncol=1)
#' covar.in <- visit.1.rest 
#' 
#' 
#' ## note there is missingness in the functional predictors
#' apply(is.na(W1), 2, mean)
#' apply(is.na(W2), 2, mean)
#' 
#' ## fit two univariate models
#' pfr.obj.t1 <- pfr(Y = Y, covariates=covar.in, funcs = list(W1),     subj = id, kz = 10, kb = 50)
#' pfr.obj.t2 <- pfr(Y = Y, covariates=covar.in, funcs = list(W2),     subj = id, kz = 10, kb = 50)
#' 
#' ### one model with two functional predictors using "smooth.face"
#' ###  for smoothing predictors
#' pfr.obj.t3 <- pfr(Y = Y, covariates=covar.in, funcs = list(W1, W2), 
#'                   subj = id, kz = 10, kb = 50, nbasis=35,smooth.option="fpca.face")
#' 
#' ## plot the coefficient function and bounds
#' dev.new()
#' par(mfrow=c(2,2))
#' ran <- c(-2,.5)
#' matplot(cbind(pfr.obj.t1$BetaHat[[1]], pfr.obj.t1$Bounds[[1]]),
#'         type = 'l', lty = c(1,2,2), col = c(1,2,2), ylab = "BetaHat", 
#'         main = "CCA", xlab="Location", ylim=ran)
#' abline(h=0, col="blue")
#' matplot(cbind(pfr.obj.t2$BetaHat[[1]], pfr.obj.t2$Bounds[[1]]),
#'         type = 'l', lty = c(1,2,2), col = c(1,2,2), ylab = "BetaHat", 
#'         main = "RCST", xlab="Location", ylim=ran)
#' abline(h=0, col="blue")
#' matplot(cbind(pfr.obj.t3$BetaHat[[1]], pfr.obj.t3$Bounds[[1]]),
#'         type = 'l', lty = c(1,2,2), col = c(1,2,2), ylab = "BetaHat", 
#'         main = "CCA  - mult.", xlab="Location", ylim=ran)
#' abline(h=0, col="blue")
#' matplot(cbind(pfr.obj.t3$BetaHat[[2]], pfr.obj.t3$Bounds[[2]]),
#'         type = 'l', lty = c(1,2,2), col = c(1,2,2), ylab = "BetaHat", 
#'         main = "RCST - mult.", xlab="Location", ylim=ran)
#' abline(h=0, col="blue")
#' 
#' 
#' ##################################################################
#' # use baseline data to regress continuous outcomes on functional 
#' # predictors (continuous outcomes only recorded for case == 1)
#' ##################################################################
#' 
#' data(DTI)
#' 
#' # subset data as needed for this example
#' cca = DTI$cca[which(DTI$visit ==1 & DTI$case == 1),]
#' rcst = DTI$rcst[which(DTI$visit ==1 & DTI$case == 1),]
#' DTI = DTI[which(DTI$visit ==1 & DTI$case == 1),]

#' # note there is missingness in the functional predictors
#' apply(is.na(cca), 2, mean)
#' apply(is.na(rcst), 2, mean)
#' 
#' # fit two models with single functional predictors and plot the results
#' fit.cca = pfr(Y=DTI$pasat, funcs = cca, kz=10, kb=50, nbasis=20)
#' fit.rcst = pfr(Y=DTI$pasat, funcs = rcst, kz=10, kb=50, nbasis=20)
#' 
#' par(mfrow = c(1,2))
#' matplot(cbind(fit.cca$BetaHat[[1]], fit.cca$Bounds[[1]]),
#'         type = 'l', lty = c(1,2,2), col = c(1,2,2), ylab = "BetaHat", 
#'         main = "CCA")
#' matplot(cbind(fit.rcst$BetaHat[[1]], fit.rcst$Bounds[[1]]),
#'         type = 'l', lty = c(1,2,2), col = c(1,2,2), ylab = "BetaHat", 
#'         main = "RCST")
#' 
#' # fit a model with two functional predictors and plot the results
#' fit.cca.rcst = pfr(Y=DTI$pasat, funcs = list(cca, rcst), kz=10, kb=30, nbasis=20)
#' 
#' par(mfrow = c(1,2))
#' matplot(cbind(fit.cca.rcst$BetaHat[[1]], fit.cca.rcst$Bounds[[1]]),
#'         type = 'l', lty = c(1,2,2), col = c(1,2,2), ylab = "BetaHat", 
#'         main = "CCA")
#' matplot(cbind(fit.cca.rcst$BetaHat[[2]], fit.cca.rcst$Bounds[[2]]),
#'         type = 'l', lty = c(1,2,2), col = c(1,2,2), ylab = "BetaHat", 
#'         main = "RCST")
#' 
#' ##################################################################
#' # use baseline data to regress binary case-status outcomes on 
#' # functional predictors
#' ##################################################################
#' 
#' data(DTI)
#' 
#' # subset data as needed for this example
#' cca = DTI$cca[which(DTI$visit == 1),]
#' rcst = DTI$rcst[which(DTI$visit == 1),]
#' DTI = DTI[which(DTI$visit == 1),]
#' 
#' # fit two models with single functional predictors and plot the results
#' fit.cca = pfr(Y=DTI$case, funcs = cca, family = "binomial")
#' fit.rcst = pfr(Y=DTI$case, funcs = rcst, family = "binomial")
#' 
#' par(mfrow = c(1,2))
#' matplot(cbind(fit.cca$BetaHat[[1]], fit.cca$Bounds[[1]]),
#'         type = 'l', lty = c(1,2,2), col = c(1,2,2), ylab = "BetaHat", 
#'         main = "CCA")
#' matplot(cbind(fit.rcst$BetaHat[[1]], fit.rcst$Bounds[[1]]),
#'         type = 'l', lty = c(1,2,2), col = c(1,2,2), ylab = "BetaHat", 
#'         main = "RCST")
#' 
#' ##################################################################
#' #########              Octane Data Example               #########
#' ##################################################################
#' 
#' data(gasoline)
#' Y = gasoline$octane
#' funcs = gasoline$NIR
#' wavelengths = as.matrix(2*450:850)
#' 
#' # fit the model using pfr and the smoothing option "fpca.face"
#' fit = pfr(Y=Y, funcs=funcs, kz=15, kb=50,nbasis=35,smooth.option="fpca.face")
#' 
# plot the estimated coefficient function
#' matplot(wavelengths, cbind(fit$BetaHat[[1]], fit$Bounds[[1]]), 
#'         type='l', lwd=c(2,1,1), lty=c(1,2,2), xlab = "Wavelengths", 
#'         ylab = "Coefficient Function", col=1)
#' }
pfr_old  <- function (Y, subj=NULL, covariates = NULL, funcs, kz = 10, kb = 30, nbasis=10,
                            family = "gaussian", method="REML", smooth.option="fpca.sc", pve=0.99,...)
{
  ## Step 1:
  ## parse formulae, etc.
  ## parse.pfr() In progress.

  ## Step 2:
  ## Preprocess in prep for gam() fit.
  pre <- preprocess.pfr(subj=subj, covariates=covariates, funcs=funcs, kz=kz, kb=kb, nbasis=nbasis, smooth.option=smooth.option, pve=pve)

  ## Step 3:
  ## gam() fit.
  fit = with(pre, gam(Y ~ X - 1, paraPen = list(X = D), method = method, family = family, ...))

  ## Step 4:
  ## Postprocess objects within "fit" to be of use to user and rlrt.pfr(), predict.pfr(), plot.pfr()
  pos <- postprocess.pfr(fit=fit, X=pre$X, p=pre$p, N_subj=pre$N_subj, phi=pre$phi, subj=subj,
                         N.Pred=pre$N.Pred, kb=kb)  

  ## Step 5:  return everything.
  ret <- c(pos, pre, list(Y=Y))
  ret
}
