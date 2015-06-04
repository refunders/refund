#' Extract coefficient functions from a fitted pfr-object
#' 
#' @param object return object from \code{\link{pfr}}
#' @param which which terms to extract coefficients from (integer vector), defaults
#'    to all terms
#' @param n see \code{\link[mgcv]{plot.gam}}
#' @param n2 see \code{\link[mgcv]{plot.gam}}
#' @param se compute standard errors? defaults to TRUE. See \code{\link[mgcv]{plot.gam}}
#' @param seWithMean see \code{\link[mgcv]{plot.gam}}
#' @param ... further arguments passed on to \pkg{mgcv}'s \code{plot.mgcv.smooth}
#' @return a list of data.frames containing the evaluation points, 
#'    coefficient function values and optionally their se's for each term in \code{which}
#' @author Fabian Scheipl, se-computation parts adapted from Simon Wood's \code{plot.gam}.
#' @examples 
#' #TODO: see pfr 
#' @export
coefficients.pfr <- function(object, which, n=100, n2=40, se=TRUE, 
                             seWithMean=TRUE, exclude=FALSE, ...){
  if (missing(which)) {
    which <- seq_along(object$smooth)
  } 
  coef <- vector(length(which), mode="list")
  for (i in which) {
    is.re <- "random.effect" %in% class(object$smooth[[i]])
    plotf <- if(is.re) {
      mgcv:::plot.random.effect
    } else {
      mgcv:::plot.mgcv.smooth
    }
    plotdata <- plotf(object$smooth[[i]], P=NULL, n=n, 
                                         n2=n2, data = object$model, ...)
    if(is.re) plotdata$x <- levels(object$model[[object$smooth[[i]]$term]])
    first <- object$smooth[[i]]$first.para
    last <- object$smooth[[i]]$last.para
    coef[[i]] <- list()
    coef[[i]]$value <- drop(plotdata$X %*% object$coefficients[first:last])
    
    # ctrl-c-v from plot.mgcv.smooth :
    if (exclude & !is.null(plotdata$exclude))
      coef[[i]]$value[plotdata$exclude] <- NA
    if (se && plotdata$se) { ## get standard errors for fit
      ## test whether mean variability to be added to variability (only for centred terms)
      if (seWithMean && attr(object$smooth[[i]], "nCons") > 0) {
        if (length(object$cmX) < ncol(object$Vp)){
          object$cmX <- c(object$cmX,rep(0,ncol(object$Vp)-length(object$cmX)))
        } 
        X1 <- matrix(object$cmX, nrow(plotdata$X), ncol(object$Vp), byrow=TRUE)
        meanL1 <- object$smooth[[i]]$meanL1
        if (!is.null(meanL1)) {
          X1 <- X1 / meanL1
        } 
        X1[,first:last] <- plotdata$X
        coef[[i]]$se <- sqrt(pmax(0,rowSums((X1 %*% object$Vp) * X1)))
      } else {
        coef[[i]]$se <- ## se in centred (or anyway unconstained) space only
          sqrt(pmax(0, 
                    rowSums((plotdata$X %*% 
                               object$Vp[first:last, first:last, drop=FALSE]) * 
                              plotdata$X)))
      }
      if (exclude & !is.null(plotdata$exclude)) {
        coef[[i]]$se[plotdata$exclude] <- NA
      }
    }
    if(object$smooth[[i]]$dim == 1) {
      if(!is.re) {
        coef[[i]][[gsub("tmat", "argvals", plotdata$xlab)]] <- plotdata$x
      } else {
        coef[[i]][[object$smooth[[i]]$term]] <- plotdata$x
      }
      
    } else {
      grid <- expand.grid(x=plotdata$x, y=plotdata$y)
      coef[[i]][[gsub("tmat", "argvals", plotdata$ylab)]] <- grid$y
      coef[[i]][[gsub("\\.omat", "", plotdata$xlab)]] <- grid$x
    }
    coef[[i]] <- as.data.frame(coef[[i]])
  } #end i in which
  if(length(which) == 1) {
    #special extra wish from Jon:
    coef[[1]]
  } else {
    coef
  }
}
#' @rdname coefficients.pfr
coef.pfr <- coefficients.pfr