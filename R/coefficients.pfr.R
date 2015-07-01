#' Extract coefficient functions from a fitted pfr-object
#' 
#' @param object return object from \code{\link{pfr}}
#' @param which which of the smooth terms to extract coefficients from
#'    (integer vector). If \code{NULL}, returns all smooth terms.
#' @param se compute standard errors? defaults to TRUE. See \code{\link[mgcv]{plot.gam}}
#' @param exclude if \code{TRUE}, excludes reporting of the estimate at coordinates that are
#'   "too far" from data used to fit the model, as determined by
#'   \code{mgcv::plot.mgcv.smooth}, by setting the estimate to \code{NA}.
#' @param ... further arguments passed on to \pkg{mgcv}'s \code{plot.gam}. Common
#'   arguments include \code{n} and \code{n2} to set the number of coordinates
#'   to estimate for 1-D and 2-D coefficient functions, and \code{seWithMean}
#'   if the standard error should include uncertainty about the overall mean.
#'   See \code{\link[mgcv]{plot.gam}}.
#' @return a list of data.frames containing the evaluation points, 
#'    coefficient function values and optionally their se's for each term in \code{which}.
#'    If only one term is selected, the one data frame is unlisted.
#' @author Fabian Scheipl and Jonathan Gellar, se-computation parts adapted from
#'    Simon Wood's \code{plot.gam}.
#' @examples 
#' #TODO: see ?pfr 
#' @export

coefficients.pfr <- function(object, which=NULL, se=TRUE, exclude=FALSE, ...) {
  
  # If (se==TRUE), replace with 1 (indicating we want 1 SE returned by plot.gam)
  if (is.logical(se)) if (se) se <- 1
  
  # exclude sets the too.far value appropriately
  too.far <- list(...)$too.far
  too.far <- ifelse(is.null(too.far), ifelse(exclude, 0.1, 0), too.far)
  
  plotdata <- mgcv::plot.gam(object, select=which, se=se, too.far=too.far, ...)
  
  coef <- lapply(1:length(plotdata), function(i) {
    pd <- plotdata[[i]]
    
    # Create coef.i with coordinates
    is.re <- "random.effect" %in% class(object$smooth[[i]])
    if(is.re) pd$x <- levels(object$model[[object$smooth[[i]]$term]])
    
    coef.i <- if(object$smooth[[i]]$dim == 1) {
      setNames(data.frame(pd$x),
               ifelse(is.re, object$smooth[[i]]$term,
                      gsub("tmat", "argvals", pd$xlab)))
    } else {
      grid <- expand.grid(x=pd$x, y=pd$y)
      setNames(data.frame(grid$x, grid$y), 
               c(gsub("\\.omat", "", pd$xlab),
                 gsub("tmat", "argvals", pd$ylab)))
    }
    
    # Add estimate & SE
    coef.i$value <- pd$fit
    if (se)
      coef.i$se  <- pd$se
    
    coef.i
  })
  
  #special extra wish from Jon:
  if(length(coef) == 1) {
    coef[[1]]
  } else {
    names(coef) <- sapply(object$smooth, function(x) x$label)
    coef
  }
}
#' @rdname coefficients.pfr
coef.pfr <- coefficients.pfr