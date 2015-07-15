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
    loc[length(loc)] <- loc[length(loc)] + 2
    
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
                  exclude <- rep(FALSE, length(yy))
                  #exclude <- exclude.too.far(xx, yy, raw$x, raw$y, dist = too.far)
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

# # @export
# plot.mgcv.smooth <- function (x, P = NULL, data = NULL, label = "", se1.mult = 1, 
#                               se2.mult = 2, partial.resids = FALSE, rug = TRUE, se = TRUE, 
#                               scale = -1, n = 100, n2 = 40, pers = FALSE, theta = 30, phi = 30, 
#                               jit = FALSE, xlab = NULL, ylab = NULL, main = NULL, ylim = NULL, 
#                               xlim = NULL, too.far = 0.1, shade = FALSE, shade.col = "gray80", 
#                               shift = 0, trans = I, by.resids = FALSE, scheme = 0, ...) 
# {
#   sp.contour <- function(x, y, z, zse, xlab = "", ylab = "", 
#                          zlab = "", titleOnly = FALSE, se.plot = TRUE, se.mult = 1, 
#                          trans = I, shift = 0, ...) {
#     gap <- median(zse, na.rm = TRUE)
#     zr <- max(trans(z + zse + shift), na.rm = TRUE) - min(trans(z - 
#                                                                   zse + shift), na.rm = TRUE)
#     n <- 10
#     while (n > 1 && zr/n < 2.5 * gap) n <- n - 1
#     zrange <- c(min(trans(z - zse + shift), na.rm = TRUE), 
#                 max(trans(z + zse + shift), na.rm = TRUE))
#     zlev <- pretty(zrange, n)
#     yrange <- range(y)
#     yr <- yrange[2] - yrange[1]
#     xrange <- range(x)
#     xr <- xrange[2] - xrange[1]
#     ypos <- yrange[2] + yr/10
#     args <- as.list(substitute(list(...)))[-1]
#     args$x <- substitute(x)
#     args$y <- substitute(y)
#     args$type = "n"
#     args$xlab <- args$ylab <- ""
#     args$axes <- FALSE
#     do.call("plot", args)
#     cs <- (yr/10)/strheight(zlab)
#     if (cs > 1) 
#       cs <- 1
#     tl <- strwidth(zlab)
#     if (tl * cs > 3 * xr/10) 
#       cs <- (3 * xr/10)/tl
#     args <- as.list(substitute(list(...)))[-1]
#     n.args <- names(args)
#     zz <- trans(z + shift)
#     args$x <- substitute(x)
#     args$y <- substitute(y)
#     args$z <- substitute(zz)
#     if (!"levels" %in% n.args) 
#       args$levels <- substitute(zlev)
#     if (!"lwd" %in% n.args) 
#       args$lwd <- 2
#     if (!"labcex" %in% n.args) 
#       args$labcex <- cs * 0.65
#     if (!"axes" %in% n.args) 
#       args$axes <- FALSE
#     if (!"add" %in% n.args) 
#       args$add <- TRUE
#     do.call("contour", args)
#     if (is.null(args$cex.main)) 
#       cm <- 1
#     else cm <- args$cex.main
#     if (titleOnly) 
#       title(zlab, cex.main = cm)
#     else {
#       xpos <- xrange[1] + 3 * xr/10
#       xl <- c(xpos, xpos + xr/10)
#       yl <- c(ypos, ypos)
#       lines(xl, yl, xpd = TRUE, lwd = args$lwd)
#       text(xpos + xr/10, ypos, zlab, xpd = TRUE, pos = 4, 
#            cex = cs * cm, off = 0.5 * cs * cm)
#     }
#     if (is.null(args$cex.axis)) 
#       cma <- 1
#     else cma <- args$cex.axis
#     axis(1, cex.axis = cs * cma)
#     axis(2, cex.axis = cs * cma)
#     box()
#     if (is.null(args$cex.lab)) 
#       cma <- 1
#     else cma <- args$cex.lab
#     mtext(xlab, 1, 2.5, cex = cs * cma)
#     mtext(ylab, 2, 2.5, cex = cs * cma)
#     if (!"lwd" %in% n.args) 
#       args$lwd <- 1
#     if (!"lty" %in% n.args) 
#       args$lty <- 2
#     if (!"col" %in% n.args) 
#       args$col <- 2
#     if (!"labcex" %in% n.args) 
#       args$labcex <- cs * 0.5
#     zz <- trans(z + zse + shift)
#     args$z <- substitute(zz)
#     do.call("contour", args)
#     if (!titleOnly) {
#       xpos <- xrange[1]
#       xl <- c(xpos, xpos + xr/10)
#       lines(xl, yl, xpd = TRUE, lty = args$lty, col = args$col)
#       text(xpos + xr/10, ypos, paste("-", round(se.mult), 
#                                      "se", sep = ""), xpd = TRUE, pos = 4, cex = cs * 
#              cm, off = 0.5 * cs * cm)
#     }
#     if (!"lty" %in% n.args) 
#       args$lty <- 3
#     if (!"col" %in% n.args) 
#       args$col <- 3
#     zz <- trans(z - zse + shift)
#     args$z <- substitute(zz)
#     do.call("contour", args)
#     if (!titleOnly) {
#       xpos <- xrange[2] - xr/5
#       xl <- c(xpos, xpos + xr/10)
#       lines(xl, yl, xpd = TRUE, lty = args$lty, col = args$col)
#       text(xpos + xr/10, ypos, paste("+", round(se.mult), 
#                                      "se", sep = ""), xpd = TRUE, pos = 4, cex = cs * 
#              cm, off = 0.5 * cs * cm)
#     }
#   }
#   if (is.null(P)) {
#     if (!x$plot.me || x$dim > 2) 
#       return(NULL)
#     if (x$dim == 1) {
#       raw <- data[x$term][[1]]
#       if (is.null(xlim)) 
#         xx <- seq(min(raw), max(raw), length = n)
#       else xx <- seq(xlim[1], xlim[2], length = n)
#       if (x$by != "NA") {
#         by <- rep(1, n)
#         dat <- data.frame(x = xx, by = by)
#         names(dat) <- c(x$term, x$by)
#       }
#       else {
#         dat <- data.frame(x = xx)
#         names(dat) <- x$term
#       }
#       X <- PredictMat(x, dat)
#       if (is.null(xlab)) 
#         xlabel <- x$term
#       else xlabel <- xlab
#       if (is.null(ylab)) 
#         ylabel <- label
#       else ylabel <- ylab
#       if (is.null(xlim)) 
#         xlim <- range(xx)
#       return(list(X = X, x = xx, scale = TRUE, se = TRUE, 
#                   raw = raw, xlab = xlabel, ylab = ylabel, main = main, 
#                   se.mult = se1.mult, xlim = xlim))
#     }
#     else {
#       xterm <- x$term[1]
#       if (is.null(xlab)) 
#         xlabel <- xterm
#       else xlabel <- xlab
#       yterm <- x$term[2]
#       if (is.null(ylab)) 
#         ylabel <- yterm
#       else ylabel <- ylab
#       raw <- data.frame(x = as.numeric(data[xterm][[1]]), 
#                         y = as.numeric(data[yterm][[1]]))
#       n2 <- max(10, n2)
#       if (is.null(xlim)) 
#         xm <- seq(min(raw$x), max(raw$x), length = n2)
#       else xm <- seq(xlim[1], xlim[2], length = n2)
#       if (is.null(ylim)) 
#         ym <- seq(min(raw$y), max(raw$y), length = n2)
#       else ym <- seq(ylim[1], ylim[2], length = n2)
#       xx <- rep(xm, n2)
#       yy <- rep(ym, rep(n2, n2))
#       if (too.far > 0) 
#         #exclude <- exclude.too.far(xx, yy, raw$x, raw$y, 
#         #                           dist = too.far)
#         exclude <- rep(NA, length(yy))
#       else exclude <- rep(FALSE, n2 * n2)
#       if (x$by != "NA") {
#         by <- rep(1, n2^2)
#         dat <- data.frame(x = xx, y = yy, by = by)
#         names(dat) <- c(xterm, yterm, x$by)
#       }
#       else {
#         dat <- data.frame(x = xx, y = yy)
#         names(dat) <- c(xterm, yterm)
#       }
#       X <- PredictMat(x, dat)
#       if (is.null(main)) {
#         main <- label
#       }
#       if (is.null(ylim)) 
#         ylim <- range(ym)
#       if (is.null(xlim)) 
#         xlim <- range(xm)
#       return(list(X = X, x = xm, y = ym, scale = FALSE, 
#                   se = TRUE, raw = raw, xlab = xlabel, ylab = ylabel, 
#                   main = main, se.mult = se2.mult, ylim = ylim, 
#                   xlim = xlim, exclude = exclude))
#     }
#   }
#   else {
#     if (se) {
#       if (x$dim == 1) {
#         if (scheme == 1) 
#           shade <- TRUE
#         ul <- P$fit + P$se
#         ll <- P$fit - P$se
#         if (scale == 0 && is.null(ylim)) {
#           ylimit <- c(min(ll), max(ul))
#           if (partial.resids) {
#             max.r <- max(P$p.resid, na.rm = TRUE)
#             if (max.r > ylimit[2]) 
#               ylimit[2] <- max.r
#             min.r <- min(P$p.resid, na.rm = TRUE)
#             if (min.r < ylimit[1]) 
#               ylimit[1] <- min.r
#           }
#         }
#         if (!is.null(ylim)) 
#           ylimit <- ylim
#         if (shade) {
#           plot(P$x, trans(P$fit + shift), type = "n", 
#                xlab = P$xlab, ylim = trans(ylimit + shift), 
#                xlim = P$xlim, ylab = P$ylab, main = P$main, 
#                ...)
#           polygon(c(P$x, P$x[n:1], P$x[1]), trans(c(ul, 
#                                                     ll[n:1], ul[1]) + shift), col = shade.col, 
#                   border = NA)
#           lines(P$x, trans(P$fit + shift), ...)
#         }
#         else {
#           plot(P$x, trans(P$fit + shift), type = "l", 
#                xlab = P$xlab, ylim = trans(ylimit + shift), 
#                xlim = P$xlim, ylab = P$ylab, main = P$main, 
#                ...)
#           if (is.null(list(...)[["lty"]])) {
#             lines(P$x, trans(ul + shift), lty = 2, ...)
#             lines(P$x, trans(ll + shift), lty = 2, ...)
#           }
#           else {
#             lines(P$x, trans(ul + shift), ...)
#             lines(P$x, trans(ll + shift), ...)
#           }
#         }
#         if (partial.resids && (by.resids || x$by == "NA")) {
#           if (length(P$raw) == length(P$p.resid)) {
#             if (is.null(list(...)[["pch"]])) 
#               points(P$raw, trans(P$p.resid + shift), 
#                      pch = ".", ...)
#             else points(P$raw, trans(P$p.resid + shift), 
#                         ...)
#           }
#           else {
#             warning("Partial residuals do not have a natural x-axis location for linear functional terms")
#           }
#         }
#         if (rug) {
#           if (jit) 
#             rug(jitter(as.numeric(P$raw)), ...)
#           else rug(as.numeric(P$raw), ...)
#         }
#       }
#       else if (x$dim == 2) {
#         P$fit[P$exclude] <- NA
#         if (pers) 
#           scheme <- 1
#         if (scheme == 1) {
#           persp(P$x, P$y, matrix(trans(P$fit + shift), 
#                                  n2, n2), xlab = P$xlab, ylab = P$ylab, zlab = P$main, 
#                 ylim = P$ylim, xlim = P$xlim, theta = theta, 
#                 phi = phi, ...)
#         }
#         else if (scheme == 2) {
#           image(P$x, P$y, matrix(trans(P$fit + shift), 
#                                  n2, n2), xlab = P$xlab, ylab = P$ylab, main = P$main, 
#                 xlim = P$xlim, ylim = P$ylim, col = heat.colors(50), 
#                 ...)
#           contour(P$x, P$y, matrix(trans(P$fit + shift), 
#                                    n2, n2), add = TRUE, col = 3, ...)
#           if (rug) {
#             if (is.null(list(...)[["pch"]])) 
#               points(P$raw$x, P$raw$y, pch = ".", ...)
#             else points(P$raw$x, P$raw$y, ...)
#           }
#         }
#         else {
#           sp.contour(P$x, P$y, matrix(P$fit, n2, n2), 
#                      matrix(P$se, n2, n2), xlab = P$xlab, ylab = P$ylab, 
#                      zlab = P$main, titleOnly = !is.null(main), 
#                      se.mult = 1, trans = trans, shift = shift, 
#                      ...)
#           if (rug) {
#             if (is.null(list(...)[["pch"]])) 
#               points(P$raw$x, P$raw$y, pch = ".", ...)
#             else points(P$raw$x, P$raw$y, ...)
#           }
#         }
#       }
#       else {
#         warning("no automatic plotting for smooths of more than two variables")
#       }
#     }
#     else {
#       if (x$dim == 1) {
#         if (scale == 0 && is.null(ylim)) {
#           if (partial.resids) 
#             ylimit <- range(P$p.resid, na.rm = TRUE)
#           else ylimit <- range(P$fit)
#         }
#         if (!is.null(ylim)) 
#           ylimit <- ylim
#         plot(P$x, trans(P$fit + shift), type = "l", xlab = P$xlab, 
#              ylab = P$ylab, ylim = trans(ylimit + shift), 
#              xlim = P$xlim, main = P$main, ...)
#         if (rug) {
#           if (jit) 
#             rug(jitter(as.numeric(P$raw)), ...)
#           else rug(as.numeric(P$raw), ...)
#         }
#         if (partial.resids && (by.resids || x$by == "NA")) {
#           if (is.null(list(...)[["pch"]])) 
#             points(P$raw, trans(P$p.resid + shift), pch = ".", 
#                    ...)
#           else points(P$raw, trans(P$p.resid + shift), 
#                       ...)
#         }
#       }
#       else if (x$dim == 2) {
#         P$fit[P$exclude] <- NA
#         if (!is.null(main)) 
#           P$title <- main
#         if (pers) 
#           scheme <- 1
#         if (scheme == 1) {
#           persp(P$x, P$y, matrix(trans(P$fit + shift), 
#                                  n2, n2), xlab = P$xlab, ylab = P$ylab, zlab = P$main, 
#                 theta = theta, phi = phi, xlim = P$xlim, 
#                 ylim = P$ylim, ...)
#         }
#         else if (scheme == 2) {
#           image(P$x, P$y, matrix(trans(P$fit + shift), 
#                                  n2, n2), xlab = P$xlab, ylab = P$ylab, main = P$main, 
#                 xlim = P$xlim, ylim = P$ylim, col = heat.colors(50), 
#                 ...)
#           contour(P$x, P$y, matrix(trans(P$fit + shift), 
#                                    n2, n2), add = TRUE, col = 3, ...)
#           if (rug) {
#             if (is.null(list(...)[["pch"]])) 
#               points(P$raw$x, P$raw$y, pch = ".", ...)
#             else points(P$raw$x, P$raw$y, ...)
#           }
#         }
#         else {
#           contour(P$x, P$y, matrix(trans(P$fit + shift), 
#                                    n2, n2), xlab = P$xlab, ylab = P$ylab, main = P$main, 
#                   xlim = P$xlim, ylim = P$ylim, ...)
#           if (rug) {
#             if (is.null(list(...)[["pch"]])) 
#               points(P$raw$x, P$raw$y, pch = ".", ...)
#             else points(P$raw$x, P$raw$y, ...)
#           }
#         }
#       }
#       else {
#         warning("no automatic plotting for smooths of more than one variable")
#       }
#     }
#   }
# }
