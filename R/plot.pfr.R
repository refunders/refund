#' Plot a pfr object
#' 
#' This function plots the smooth coefficients of a pfr object. These include
#' functional coefficients as well as any smooths of scalar covariates. The
#' function dispatches to \code{plot.gam}, after injecting a small amount of
#' code into the function to handle special cases.
#' 
#' @param x a fitted \code{pfr}-object
#' @param Qtransform For additive functional terms fit with
#'   \code{af(Qtransform=TRUE)}, \code{TRUE} indicates the coefficient should be
#'   plotted on the quantile-transformed scale, whereas \code{FALSE} indicates
#'   the scale of the original data
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

plot.pfr <- function(x, Qtransform=FALSE, ...) {
  class(x) <- class(x)[-1]
  
  if (Qtransform) {
    # Inject code into plot.gam to rescale
    bod <- as.list(body(mgcv:::plot.mgcv.smooth))
    loc <- locfcn(bod, "rep(ym, rep(n2, n2))")
    loc[length(loc)] <- locs2[length(loc)] + 2
    
    suppressMessages(
      trace(mgcv:::plot.mgcv.smooth,
            at=list(loc), print=FALSE, tracer = quote({
              # Inserted Code
              if (!is.null(x$QT)) {
                tf <- x$tf[[1]]
                raw$y <- tf(raw$x, raw$y)
                ym <- seq(min(raw$y), max(raw$y), length = n2)
                yy <- rep(ym, rep(n2, n2))
                if (too.far > 0) 
                  exclude <- exclude.too.far(xx, yy, raw$x, raw$y, dist = too.far)
                x0 <- environment(tf)$x0
                t0 <- environment(tf)$t0
                idx <- factor(t0)
                newidx <- factor(xx)
                tmp <- tapply(x0, t0, function(y) y,
                              simplify=F)
                for (lev in levels(newidx)) {
                  yy[newidx==lev] <- if (lev %in% levels(idx)) {
                    quantile(tmp[[which(levels(idx)==lev)]], yy[newidx==lev])
                  } else {
                    u1 <- as.numeric(levels(idx))
                    idx1 <- which(u1 == max(u1[u1<as.numeric(lev)]))
                    idx2 <- which(u1 == min(u1[u1>as.numeric(lev)]))
                    bounds <- sapply(c(idx1, idx2), function(i) {
                      quantile(tmp[[i]], yy[newidx==lev])
                    })
                    apply(bounds, 1, function(y) {
                      approx(as.numeric(levels(idx)[idx1:idx2]), c(y[1], y[2]),
                             xout = as.numeric(lev), rule=2)$y
                    })
                  }
                }                
              }
            })))
    on.exit({
      suppressMessages(try(untrace(mgcv:::plot.mgcv.smooth), silent = TRUE))
    })

  } else {
    # Inject code into plot.gam to exclude coordinates outside range
    bod <- as.list(body(plot.gam))
    loc <- which(sapply(bod, function(x)
      any(grepl(x, pattern="is.null(P)", fixed=TRUE)))) + 1
    suppressMessages(
      trace(plot.gam,
            at=loc, print=FALSE, tracer = quote({
              # Inserted Code
              for (i in 1:m) if (!is.null(x$smooth[[i]]$QT)) {
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

locfcn <- function(bod, txt) {
  loc <- which(sapply(bod, function(x)
    any(grepl(x, pattern=txt, fixed=TRUE))))
  ret <- if (length(loc)) {
    newbod <- bod[[loc]]
    c(loc, locfcn(newbod, txt))
  } else c()
  ret
}
