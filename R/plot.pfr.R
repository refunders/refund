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
#' @importFrom mgcv plot.gam
#' @export

plot.pfr <- function(x, Qscale=FALSE, ...) {
  class(x) <- class(x)[-1]
  if (Qscale) {
    stop("Qscale not yet implemented")
    # Inject code into plot.gam to rescale
    suppressMessages(
      trace(plot.gam,
            at=which(sapply(as.list(body(plot.gam)), function(x)
              any(grepl(x, pattern="is.null(P)", fixed=TRUE)))) + 1,
            print=FALSE,
            tracer = quote({
              # Inserted Code
              
              
            })
      ))
    on.exit({
      suppressMessages(try(untrace(plot.gam), silent = TRUE))
    })
    
    
  } else {
    # Inject code into plot.gam to exclude coordinates outside range
    bod <- as.list(body(plot.gam))
    locvec <- which(sapply(bod, function(x)
      any(grepl(x, pattern="is.null(P)", fixed=TRUE)))) + 1
    #locvec <- locfcn(bod)
    
    suppressMessages(
      trace(plot.gam,
            at=list(locvec), print=FALSE, tracer = quote({
              # Inserted Code
              for (i in 1:m) if (!is.null(x$smooth[[i]]$QTransform)) {
                tf <- x$smooth[[i]]$tf[[1]]
                if (!is.null(tf)) {
                  rna <- environment(tf)$retNA
                  environment(tf)$retNA <- TRUE
                  cgrid <- expand.grid(pd[[i]]$x, pd[[i]]$y)
                  new.exclude <- is.na(tf(cgrid[,1], cgrid[,2]))
                  environment(tf)$retNA <- rna
                  if (!is.null(pd[[i]]$exclude))
                    pd[[i]]$exclude[new.exclude] <- TRUE
                  pd[[i]]$fit[new.exclude] <- NA
                  if (se && pd[[i]]$se)
                    pd[[i]]$se.fit[new.exclude] <- NA
                }}})))
    on.exit({
      suppressMessages(try(untrace(plot.gam), silent = TRUE))
    })
    
  }
  
  plot(x, ...)
}

#locfcn <- function(bod) {
#  loc <- which(sapply(bod, function(x)
#    any(grepl(x, pattern="is.null(P)", fixed=TRUE))))
#  ret <- if (length(loc)) {
#    newbod <- bod[[loc]]
#    c(loc, locfcn(newbod))
#  } else c()
#  ret
#}
