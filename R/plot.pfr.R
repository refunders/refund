#' Plot a pfr object
#' 
#' This function plots the smooth coefficients of a pfr object. These include
#' functional coefficients as well as any smooths of scalar covariates. The
#' function simply dispatches to \code{plot.gam}.
#' 
#' @param x a fitted \code{pfr}-object
#' @param Qscale For additive functional terms fit with
#'   \code{af(Qtransform=TRUE)}, should the coefficient be plotted on the
#'   raw scale or the quantile-transformed scale?
#' @param ... arguments handed over to \code{\link[mgcv]{plot.gam}}
#' 
#' @return This function's main purpose is its side effect of generating plots.
#' It also silently returns a list of the data used to produce the plots, which
#' can be used to generate customized plots.
#' 
#' @author Jonathan Gellar
#' @seealso \code{\link{af}}, \code{\link{pfr}}
#' @export

plot.pfr <- function(x, Qscale=FALSE, ...) {
  class(x) <- class(x)[-1]
  
  if (Qscale) {
    
    
    # Inject code into plot.gam to rescale
    suppressMessages(
      trace(plot.gam,
            at=which(sapply(as.list(body(plot.gam)), function(x)
              
              
              
              any(grepl(x, pattern="mf[[timetrans$var[i]]]", fixed=TRUE)))) + 1,
            print=FALSE,
            tracer = quote({
              # Inserted Code
              
              
            })
      ))
    on.exit({
      suppressMessages(try(untrace(plot.gam), silent = TRUE))
    })
    
    
    
  }
  
  
  plot(x, ...)
  
  
}