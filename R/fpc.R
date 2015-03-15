##' Construct a FPC regression term
##' 
##' Constructs a functional principal component regression (Reiss and Ogden, 
##' 2007, 2010) term for inclusion in an \code{mgcv::gam}-formula (or
##' \code{\link{bam}} or \code{\link{gamm}} or \code{gamm4:::gamm}) as
##' constructed by \code{\link{pfr}}. The predictor \code{X} can either be a
##' one-dimensional function, or a two-dimensional image.
##' 
##' @param X functional predictors. For 1D predictors, this is expressed as an
##'   \code{N} by \code{J} matrix, where \code{N} is the number of subjects
##'   and \code{J} is the number of evaluation points. For 2D predictors, this
##'   is expressed as an \code{N} by {J1} by {J2} array, representing \code{N}
##'   images of dimension \code{J1} by \code{J2}.
##' @param argvals indices of evaluation of \code{X}. 
##' @param integration method used for numerical integration
##' @param ncomp number of principal components. if \code{NULL}, chosen by \code{pve}
##' @param pve proportion of variance explained; used to choose the number of
##'   principal components
##' @param basistype for 2D smooths, selects the function for defining the
##'   two-dimensional basis.
##' @param ... additional options to be passed to \code{mgcv}'s \code{\link{s}}
##'   function for defining the pre-smoothing basis
##' 
##' @details We implement the FPCR-R method of Reiss and Ogden (2007). This
##'   method first performs a penalized basis expansion of the functional
##'   predictors, and then calculates the principal components of the resulting
##'   smooth functions. The loadings for each predictor function are used as the
##'   covariates in a linear model.
##' 
##' 

fpc <- function(X, argvals=seq(0, 1, l=ncol(X)),
                integration = c("simpson", "trapezoidal", "riemann"),
                ncomp=NULL, pve=0.99, basistype=c("te", "s", "t2"),
                ...) {
  
  
  
}

# fpcr <- function(y, xfuncs = NULL, fdobj = NULL, ncomp=NULL, pve = 0.99, nbasis = NULL, basismat = NULL, 
#                  penmat = NULL, argvals = NULL, covt = NULL, mean.signal.term = FALSE, 
#                  spline.order = NULL, family = "gaussian", method = "REML", sp = NULL, 
#                  pen.order = 2, cv1 = FALSE, nfold = 5, store.cv = FALSE, store.gam = TRUE, ...){
#   # Error handling
# 
# lf <- function(X, argvals = seq(0, 1, l = ncol(X)), xind = NULL,
#                integration = c("simpson", "trapezoidal", "riemann"),
#                L = NULL, presmooth = NULL, presmooth.opts = NULL, ...)
#   
#   