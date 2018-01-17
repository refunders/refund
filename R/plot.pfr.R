#' Plot a pfr object
#'
#' This function plots the smooth coefficients of a pfr object. These include
#' functional coefficients as well as any smooths of scalar covariates. The
#' function dispatches to \code{pfr_plot.gam}, which is our local copy of
#' \code{\link[mgcv]{plot.gam}} with some minor changes.
#'
#' @param x a fitted \code{pfr}-object
#' @param Qtransform For additive functional terms, \code{TRUE} indicates the
#'   coefficient should be plotted on the quantile-transformed scale, whereas
#'   \code{FALSE} indicates the scale of the original data. Note this is
#'   different from the \code{Qtransform} arguemnt of \code{af}, which specifies
#'   the scale on which the term is fit.
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
  smooth.types <- x$pfr$termtype[x$pfr$termtype != "par"]

  if (Qtransform) {
    # Flag the smooths
    for (i in 1:length(x$smooth))
      if (smooth.types[i]=="af")
        x$smooth[[i]]$Qplot <- TRUE
  }
  pfr_plot.gam(x, ...)
}

# locfcn <- function(bod, txt) {
#   loc <- which(sapply(bod, function(x)
#     any(grepl(x, pattern=txt, fixed=TRUE))))
#   ret <- if (length(loc)) {
#     newbod <- bod[[loc]]
#     c(loc, locfcn(newbod, txt))
#   } else c()
#   ret
# }


#' Local version of \code{plot.gam}
#'
#' These internal functions were copied from Simon Wood's \code{mgcv} package,
#' with some minor changes to allow for plotting \code{pfr} objects.
#'
#' @seealso \code{\link[mgcv]{plot.gam}}
#' @keywords internal
pfr_plot.gam <- function (x, residuals = FALSE, rug = TRUE, se = TRUE, pages = 0,
                         select = NULL, scale = -1, n = 100, n2 = 40, n3=3, pers = FALSE,
                         theta = 30, phi = 30, jit = FALSE, xlab = NULL, ylab = NULL,
                         main = NULL, ylim = NULL, xlim = NULL, too.far = 0.1, all.terms = FALSE,
                         shade = FALSE, shade.col = "gray80", shift = 0, trans = I,
                         seWithMean = FALSE, unconditional = FALSE, by.resids = FALSE,
                         scheme = 0, ...) {
  sub.edf <- function(lab, edf) {
    pos <- regexpr(":", lab)[1]
    if (pos < 0) {
      pos <- nchar(lab) - 1
      lab <- paste(substr(lab, start = 1, stop = pos),
                   ",", round(edf, digits = 2), ")", sep = "")
    }
    else {
      lab1 <- substr(lab, start = 1, stop = pos - 2)
      lab2 <- substr(lab, start = pos - 1, stop = nchar(lab))
      lab <- paste(lab1, ",", round(edf, digits = 2), lab2,
                   sep = "")
    }
    lab
  }
  if (unconditional) {
    if (is.null(x$Vc))
      warning("Smoothness uncertainty corrected covariance not available")
    else x$Vp <- x$Vc
  }
  w.resid <- NULL
  if (length(residuals) > 1) {
    if (length(residuals) == length(x$residuals))
      w.resid <- residuals
    else warning("residuals argument to plot.gam is wrong length: ignored")
    partial.resids <- TRUE
  }
  else partial.resids <- residuals
  m <- length(x$smooth)
  if (length(scheme) == 1)
    scheme <- rep(scheme, m)
  if (length(scheme) != m) {
    warn <- paste("scheme should be a single number, or a vector with",
                  m, "elements")
    warning(warn)
    scheme <- rep(scheme[1], m)
  }
  order <- if (is.list(x$pterms))
    unlist(lapply(x$pterms, attr, "order"))
  else attr(x$pterms, "order")
  if (all.terms)
    n.para <- sum(order == 1)
  else n.para <- 0
  if (se) {
    if (is.numeric(se))
      se2.mult <- se1.mult <- se
    else {
      se1.mult <- 2
      se2.mult <- 1
    }
    if (se1.mult < 0)
      se1.mult <- 0
    if (se2.mult < 0)
      se2.mult <- 0
  }
  else se1.mult <- se2.mult <- 1
  if (se && x$Vp[1, 1] < 0) {
    se <- FALSE
    warning("No variance estimates available")
  }
  if (partial.resids) {
    if (is.null(w.resid)) {
      if (is.null(x$residuals) || is.null(x$weights))
        partial.resids <- FALSE
      else {
        wr <- sqrt(x$weights)
        w.resid <- x$residuals * wr/mean(wr)
      }
    }
    if (partial.resids)
      fv.terms <- predict(x, type = "terms")
  }
  pd <- list()
  i <- 1
  if (m > 0)
    for (i in 1:m) {
      first <- x$smooth[[i]]$first.para
      last <- x$smooth[[i]]$last.para
      edf <- sum(x$edf[first:last])
      term.lab <- sub.edf(x$smooth[[i]]$label, edf)
      attr(x$smooth[[i]], "coefficients") <- x$coefficients[first:last]
      P <- plot(x$smooth[[i]], P = NULL, data = x$model,
                partial.resids = partial.resids, rug = rug, se = se,
                scale = scale, n = n, n2 = n2, pers = pers, theta = theta,
                phi = phi, jit = jit, xlab = xlab, ylab = ylab,
                main = main, label = term.lab, ylim = ylim, xlim = xlim,
                too.far = too.far, shade = shade, shade.col = shade.col,
                se1.mult = se1.mult, se2.mult = se2.mult, shift = shift,
                trans = trans, by.resids = by.resids, scheme = scheme[i],
                ...)
      if (is.null(P))
        pd[[i]] <- list(plot.me = FALSE)
      else if (is.null(P$fit)) {
        if (!is.null(x$smooth[[i]]$QT) & is.null(x$smooth[[i]]$Qplot)) {
          # af term that was fit on quantile scale, but plotting on raw scale
          # need to get rid of coordinates outside original data range
          tf <- x$smooth[[i]]$tf[[1]]
          if (!is.null(tf)) {
            rna <- environment(tf)$retNA
            environment(tf)$retNA <- TRUE
            cgrid <- expand.grid(P$x, P$y)
            new.exclude <- is.na(tf(cgrid[,1], cgrid[,2]))
            environment(tf)$retNA <- rna
            if (!is.null(P$exclude))
              P$exclude[new.exclude] <- TRUE
          }
        }

        p <- x$coefficients[first:last]
        offset <- attr(P$X, "offset")
        if (is.null(offset))
          P$fit <- P$X %*% p
        else P$fit <- P$X %*% p + offset
        if (!is.null(P$exclude))
          P$fit[P$exclude] <- NA
        if (se && P$se) {
          if (seWithMean && attr(x$smooth[[i]], "nCons") >
                0) {
            if (length(x$cmX) < ncol(x$Vp))
              x$cmX <- c(x$cmX, rep(0, ncol(x$Vp) - length(x$cmX)))
            X1 <- matrix(x$cmX, nrow(P$X), ncol(x$Vp),
                         byrow = TRUE)
            meanL1 <- x$smooth[[i]]$meanL1
            if (!is.null(meanL1))
              X1 <- X1/meanL1
            X1[, first:last] <- P$X
            se.fit <- sqrt(pmax(0, rowSums((X1 %*% x$Vp) *
                                             X1)))
          }
          else se.fit <- sqrt(pmax(0, rowSums((P$X %*%
                                                 x$Vp[first:last, first:last, drop = FALSE]) *
                                                P$X)))
          if (!is.null(P$exclude))
            P$se.fit[P$exclude] <- NA
        }
        if (partial.resids) {
          P$p.resid <- fv.terms[, length(order) + i] +
            w.resid
        }
        if (se && P$se)
          P$se <- se.fit * P$se.mult
        P$X <- NULL
        P$plot.me <- TRUE
        pd[[i]] <- P
        rm(P)
      }
      else {
        if (partial.resids) {
          P$p.resid <- fv.terms[, length(order) + i] +
            w.resid
        }
        P$plot.me <- TRUE
        pd[[i]] <- P
        rm(P)
      }
    }
  n.plots <- n.para
  if (m > 0)
    for (i in 1:m) n.plots <- n.plots + as.numeric(pd[[i]]$plot.me)
  if (n.plots == 0)
    stop("No terms to plot - nothing for plot.gam() to do.")
  if (pages > n.plots)
    pages <- n.plots
  if (pages < 0)
    pages <- 0
  if (pages != 0) {
    ppp <- n.plots%/%pages
    if (n.plots%%pages != 0) {
      ppp <- ppp + 1
      while (ppp * (pages - 1) >= n.plots) pages <- pages -
        1
    }
    c <- r <- trunc(sqrt(ppp))
    if (c < 1)
      r <- c <- 1
    if (c * r < ppp)
      c <- c + 1
    if (c * r < ppp)
      r <- r + 1
    oldpar <- par(mfrow = c(r, c))
  }
  else {
    ppp <- 1
    oldpar <- par()
  }
  if ((pages == 0 && prod(par("mfcol")) < n.plots && dev.interactive()) ||
        pages > 1 && dev.interactive())
    ask <- TRUE
  else ask <- FALSE
  if (!is.null(select)) {
    ask <- FALSE
  }
  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }
  if (scale == -1 && is.null(ylim)) {
    k <- 0
    if (m > 0)
      for (i in 1:m) if (pd[[i]]$plot.me && pd[[i]]$scale) {
        if (se && length(pd[[i]]$se) > 1) {
          ul <- pd[[i]]$fit + pd[[i]]$se
          ll <- pd[[i]]$fit - pd[[i]]$se
          if (k == 0) {
            ylim <- c(min(ll, na.rm = TRUE), max(ul,
                                                 na.rm = TRUE))
            k <- 1
          }
          else {
            if (min(ll, na.rm = TRUE) < ylim[1])
              ylim[1] <- min(ll, na.rm = TRUE)
            if (max(ul, na.rm = TRUE) > ylim[2])
              ylim[2] <- max(ul, na.rm = TRUE)
          }
        }
        else {
          if (k == 0) {
            ylim <- range(pd[[i]]$fit, na.rm = TRUE)
            k <- 1
          }
          else {
            if (min(pd[[i]]$fit, na.rm = TRUE) < ylim[1])
              ylim[1] <- min(pd[[i]]$fit, na.rm = TRUE)
            if (max(pd[[i]]$fit, na.rm = TRUE) > ylim[2])
              ylim[2] <- max(pd[[i]]$fit, na.rm = TRUE)
          }
        }
        if (partial.resids) {
          ul <- max(pd[[i]]$p.resid, na.rm = TRUE)
          if (ul > ylim[2])
            ylim[2] <- ul
          ll <- min(pd[[i]]$p.resid, na.rm = TRUE)
          if (ll < ylim[1])
            ylim[1] <- ll
        }
      }
  }
  if (m > 0)
    for (i in 1:m) if (pd[[i]]$plot.me && (is.null(select) ||
                                             i == select)) {
      plot(x$smooth[[i]], P = pd[[i]], partial.resids = partial.resids,
           rug = rug, se = se, scale = scale, n = n, n2 = n2,
           pers = pers, theta = theta, phi = phi, jit = jit,
           xlab = xlab, ylab = ylab, main = main, ylim = ylim,
           xlim = xlim, too.far = too.far, shade = shade,
           shade.col = shade.col, shift = shift, trans = trans,
           by.resids = by.resids, scheme = scheme[i], ...)
    }
  if (n.para > 0) {
    class(x) <- c("gam", "glm", "lm")
    if (is.null(select)) {
      attr(x, "para.only") <- TRUE
      termplot(x, se = se, rug = rug, col.se = 1, col.term = 1,
               main = attr(x$pterms, "term.labels"), ...)
    }
    else {
      if (select > m) {
        select <- select - m
        term.labels <- attr(x$pterms, "term.labels")
        term.labels <- term.labels[order == 1]
        if (select <= length(term.labels)) {
          termplot(x, terms = term.labels[select], se = se,
                   rug = rug, col.se = 1, col.term = 1, ...)
        }
      }
    }
  }
  if (pages > 0)
    par(oldpar)
  invisible(pd)
}


#' @rdname pfr_plot.gam
#' @keywords internal
plot.mgcv.smooth <- function (x, P = NULL, data = NULL, label = "", se1.mult = 1,
                              se2.mult = 2, partial.resids = FALSE, rug = TRUE, se = TRUE,
                              scale = -1, n = 100, n2 = 40, pers = FALSE, theta = 30, phi = 30,
                              jit = FALSE, xlab = NULL, ylab = NULL, main = NULL, ylim = NULL,
                              xlim = NULL, too.far = 0.1, shade = FALSE, shade.col = "gray80",
                              shift = 0, trans = I, by.resids = FALSE, scheme = 0, ...) {
  sp.contour <- function(x, y, z, zse, xlab = "", ylab = "",
                         zlab = "", titleOnly = FALSE, se.plot = TRUE, se.mult = 1,
                         trans = I, shift = 0, ...) {
    gap <- median(zse, na.rm = TRUE)
    zr <- max(trans(z + zse + shift), na.rm = TRUE) - min(trans(z -
                                                                  zse + shift), na.rm = TRUE)
    n <- 10
    while (n > 1 && zr/n < 2.5 * gap) n <- n - 1
    zrange <- c(min(trans(z - zse + shift), na.rm = TRUE),
                max(trans(z + zse + shift), na.rm = TRUE))
    zlev <- pretty(zrange, n)
    yrange <- range(y)
    yr <- yrange[2] - yrange[1]
    xrange <- range(x)
    xr <- xrange[2] - xrange[1]
    ypos <- yrange[2] + yr/10
    args <- as.list(substitute(list(...)))[-1]
    args$x <- substitute(x)
    args$y <- substitute(y)
    args$type = "n"
    args$xlab <- args$ylab <- ""
    args$axes <- FALSE
    do.call("plot", args)
    cs <- (yr/10)/strheight(zlab)
    if (cs > 1)
      cs <- 1
    tl <- strwidth(zlab)
    if (tl * cs > 3 * xr/10)
      cs <- (3 * xr/10)/tl
    args <- as.list(substitute(list(...)))[-1]
    n.args <- names(args)
    zz <- trans(z + shift)
    args$x <- substitute(x)
    args$y <- substitute(y)
    args$z <- substitute(zz)
    if (!"levels" %in% n.args)
      args$levels <- substitute(zlev)
    if (!"lwd" %in% n.args)
      args$lwd <- 2
    if (!"labcex" %in% n.args)
      args$labcex <- cs * 0.65
    if (!"axes" %in% n.args)
      args$axes <- FALSE
    if (!"add" %in% n.args)
      args$add <- TRUE
    do.call("contour", args)
    if (is.null(args$cex.main))
      cm <- 1
    else cm <- args$cex.main
    if (titleOnly)
      title(zlab, cex.main = cm)
    else {
      xpos <- xrange[1] + 3 * xr/10
      xl <- c(xpos, xpos + xr/10)
      yl <- c(ypos, ypos)
      lines(xl, yl, xpd = TRUE, lwd = args$lwd)
      text(xpos + xr/10, ypos, zlab, xpd = TRUE, pos = 4,
           cex = cs * cm, off = 0.5 * cs * cm)
    }
    if (is.null(args$cex.axis))
      cma <- 1
    else cma <- args$cex.axis
    axis(1, cex.axis = cs * cma)
    axis(2, cex.axis = cs * cma)
    box()
    if (is.null(args$cex.lab))
      cma <- 1
    else cma <- args$cex.lab
    mtext(xlab, 1, 2.5, cex = cs * cma)
    mtext(ylab, 2, 2.5, cex = cs * cma)
    if (!"lwd" %in% n.args)
      args$lwd <- 1
    if (!"lty" %in% n.args)
      args$lty <- 2
    if (!"col" %in% n.args)
      args$col <- 2
    if (!"labcex" %in% n.args)
      args$labcex <- cs * 0.5
    zz <- trans(z + zse + shift)
    args$z <- substitute(zz)
    do.call("contour", args)
    if (!titleOnly) {
      xpos <- xrange[1]
      xl <- c(xpos, xpos + xr/10)
      lines(xl, yl, xpd = TRUE, lty = args$lty, col = args$col)
      text(xpos + xr/10, ypos, paste("-", round(se.mult),
                                     "se", sep = ""), xpd = TRUE, pos = 4, cex = cs *
             cm, off = 0.5 * cs * cm)
    }
    if (!"lty" %in% n.args)
      args$lty <- 3
    if (!"col" %in% n.args)
      args$col <- 3
    zz <- trans(z - zse + shift)
    args$z <- substitute(zz)
    do.call("contour", args)
    if (!titleOnly) {
      xpos <- xrange[2] - xr/5
      xl <- c(xpos, xpos + xr/10)
      lines(xl, yl, xpd = TRUE, lty = args$lty, col = args$col)
      text(xpos + xr/10, ypos, paste("+", round(se.mult),
                                     "se", sep = ""), xpd = TRUE, pos = 4, cex = cs *
             cm, off = 0.5 * cs * cm)
    }
  }
  if (is.null(P)) {
    if (!x$plot.me || x$dim > 2)
      return(NULL)
    if (x$dim == 1) {
      raw <- data[x$term][[1]]
      if (is.null(xlim))
        xx <- seq(min(raw), max(raw), length = n)
      else xx <- seq(xlim[1], xlim[2], length = n)
      if (x$by != "NA") {
        by <- rep(1, n)
        dat <- data.frame(x = xx, by = by)
        names(dat) <- c(x$term, x$by)
      }
      else {
        dat <- data.frame(x = xx)
        names(dat) <- x$term
      }
      X <- PredictMat(x, dat)
      if (is.null(xlab))
        xlabel <- x$term
      else xlabel <- xlab
      if (is.null(ylab))
        ylabel <- label
      else ylabel <- ylab
      if (is.null(xlim))
        xlim <- range(xx)
      return(list(X = X, x = xx, scale = TRUE, se = TRUE,
                  raw = raw, xlab = xlabel, ylab = ylabel, main = main,
                  se.mult = se1.mult, xlim = xlim))
    }
    else {
      xterm <- x$term[1]
      if (is.null(xlab))
        xlabel <- xterm
      else xlabel <- xlab
      yterm <- x$term[2]
      if (is.null(ylab))
        ylabel <- yterm
      else ylabel <- ylab
      raw <- data.frame(x = as.numeric(data[xterm][[1]]),
                        y = as.numeric(data[yterm][[1]]))
      n2 <- max(10, n2)
      if (is.null(xlim))
        xm <- seq(min(raw$x), max(raw$x), length = n2)
      else xm <- seq(xlim[1], xlim[2], length = n2)
      if (is.null(ylim))
        ym <- seq(min(raw$y), max(raw$y), length = n2)
      else ym <- seq(ylim[1], ylim[2], length = n2)
      xx <- rep(xm, n2)
      yy <- rep(ym, rep(n2, n2))
      if (too.far > 0)
        exclude <- mgcv::exclude.too.far(xx, yy, raw$x, raw$y, dist = too.far)
      else exclude <- rep(FALSE, n2 * n2)

      if (!is.null(x$Qplot)) {
        # Transform raw data to quantile scale
        x0 <- raw$y
        t0 <- raw$x
        raw$y <- QTFunc(t0, x0, retNA=FALSE)

        # Create grid on quantile scale
        ym <- seq(min(raw$y), max(raw$y), length = n2)
        yy <- rep(ym, rep(n2, n2))
        if (too.far > 0)
          exclude <- mgcv::exclude.too.far(xx, yy, raw$x, raw$y, dist = too.far)

        # Transform back to data scale for prediction
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

      if (x$by != "NA") {
        by <- rep(1, n2^2)
        dat <- data.frame(x = xx, y = yy, by = by)
        names(dat) <- c(xterm, yterm, x$by)
      }
      else {
        dat <- data.frame(x = xx, y = yy)
        names(dat) <- c(xterm, yterm)
      }

      X <- PredictMat(x, dat)

      #if (!is.null(x$Qplot) & is.null())

      if (is.null(main)) {
        main <- label
      }
      if (is.null(ylim))
        ylim <- range(ym)
      if (is.null(xlim))
        xlim <- range(xm)
      return(list(X = X, x = xm, y = ym, scale = FALSE,
                  se = TRUE, raw = raw, xlab = xlabel, ylab = ylabel,
                  main = main, se.mult = se2.mult, ylim = ylim,
                  xlim = xlim, exclude = exclude))
    }
  }
  else {
    if (se) {
      if (x$dim == 1) {
        if (scheme == 1)
          shade <- TRUE
        ul <- P$fit + P$se
        ll <- P$fit - P$se
        if (scale == 0 && is.null(ylim)) {
          ylimit <- c(min(ll), max(ul))
          if (partial.resids) {
            max.r <- max(P$p.resid, na.rm = TRUE)
            if (max.r > ylimit[2])
              ylimit[2] <- max.r
            min.r <- min(P$p.resid, na.rm = TRUE)
            if (min.r < ylimit[1])
              ylimit[1] <- min.r
          }
        }
        if (!is.null(ylim))
          ylimit <- ylim
        if (shade) {
          plot(P$x, trans(P$fit + shift), type = "n",
               xlab = P$xlab, ylim = trans(ylimit + shift),
               xlim = P$xlim, ylab = P$ylab, main = P$main,
               ...)
          polygon(c(P$x, P$x[n:1], P$x[1]), trans(c(ul,
                                                    ll[n:1], ul[1]) + shift), col = shade.col,
                  border = NA)
          lines(P$x, trans(P$fit + shift), ...)
        }
        else {
          plot(P$x, trans(P$fit + shift), type = "l",
               xlab = P$xlab, ylim = trans(ylimit + shift),
               xlim = P$xlim, ylab = P$ylab, main = P$main,
               ...)
          if (is.null(list(...)[["lty"]])) {
            lines(P$x, trans(ul + shift), lty = 2, ...)
            lines(P$x, trans(ll + shift), lty = 2, ...)
          }
          else {
            lines(P$x, trans(ul + shift), ...)
            lines(P$x, trans(ll + shift), ...)
          }
        }
        if (partial.resids && (by.resids || x$by == "NA")) {
          if (length(P$raw) == length(P$p.resid)) {
            if (is.null(list(...)[["pch"]]))
              points(P$raw, trans(P$p.resid + shift),
                     pch = ".", ...)
            else points(P$raw, trans(P$p.resid + shift),
                        ...)
          }
          else {
            warning("Partial residuals do not have a natural x-axis location for linear functional terms")
          }
        }
        if (rug) {
          if (jit)
            rug(jitter(as.numeric(P$raw)), ...)
          else rug(as.numeric(P$raw), ...)
        }
      }
      else if (x$dim == 2) {
        P$fit[P$exclude] <- NA
        if (pers)
          scheme <- 1
        if (scheme == 1) {
          persp(P$x, P$y, matrix(trans(P$fit + shift),
                                 n2, n2), xlab = P$xlab, ylab = P$ylab, zlab = P$main,
                ylim = P$ylim, xlim = P$xlim, theta = theta,
                phi = phi, ...)
        }
        else if (scheme == 2) {
          image(P$x, P$y, matrix(trans(P$fit + shift),
                                 n2, n2), xlab = P$xlab, ylab = P$ylab, main = P$main,
                xlim = P$xlim, ylim = P$ylim, col = heat.colors(50),
                ...)
          contour(P$x, P$y, matrix(trans(P$fit + shift),
                                   n2, n2), add = TRUE, col = 3, ...)
          if (rug) {
            if (is.null(list(...)[["pch"]]))
              points(P$raw$x, P$raw$y, pch = ".")
            else points(P$raw$x, P$raw$y)
          }
        }
        else {
          sp.contour(P$x, P$y, matrix(P$fit, n2, n2),
                     matrix(P$se, n2, n2), xlab = P$xlab, ylab = P$ylab,
                     zlab = P$main, titleOnly = !is.null(main),
                     se.mult = 1, trans = trans, shift = shift,
                     ...)
          if (rug) {
            if (is.null(list(...)[["pch"]]))
              points(P$raw$x, P$raw$y, pch = ".", ...)
            else points(P$raw$x, P$raw$y, ...)
          }
        }
      }
      else {
        warning("no automatic plotting for smooths of more than two variables")
      }
    }
    else {
      if (x$dim == 1) {
        if (scale == 0 && is.null(ylim)) {
          if (partial.resids)
            ylimit <- range(P$p.resid, na.rm = TRUE)
          else ylimit <- range(P$fit)
        }
        if (!is.null(ylim))
          ylimit <- ylim
        plot(P$x, trans(P$fit + shift), type = "l", xlab = P$xlab,
             ylab = P$ylab, ylim = trans(ylimit + shift),
             xlim = P$xlim, main = P$main, ...)
        if (rug) {
          if (jit)
            rug(jitter(as.numeric(P$raw)), ...)
          else rug(as.numeric(P$raw), ...)
        }
        if (partial.resids && (by.resids || x$by == "NA")) {
          if (is.null(list(...)[["pch"]]))
            points(P$raw, trans(P$p.resid + shift), pch = ".",
                   ...)
          else points(P$raw, trans(P$p.resid + shift),
                      ...)
        }
      }
      else if (x$dim == 2) {
        P$fit[P$exclude] <- NA
        if (!is.null(main))
          P$title <- main
        if (pers)
          scheme <- 1
        if (scheme == 1) {
          persp(P$x, P$y, matrix(trans(P$fit + shift),
                                 n2, n2), xlab = P$xlab, ylab = P$ylab, zlab = P$main,
                theta = theta, phi = phi, xlim = P$xlim,
                ylim = P$ylim, ...)
        }
        else if (scheme == 2) {
          image(P$x, P$y, matrix(trans(P$fit + shift),
                                 n2, n2), xlab = P$xlab, ylab = P$ylab, main = P$main,
                xlim = P$xlim, ylim = P$ylim, col = heat.colors(50),
                ...)
          contour(P$x, P$y, matrix(trans(P$fit + shift),
                                   n2, n2), add = TRUE, col = 3, ...)
          if (rug) {
            if (is.null(list(...)[["pch"]]))
              points(P$raw$x, P$raw$y, pch = ".", ...)
            else points(P$raw$x, P$raw$y, ...)
          }
        }
        else {
          contour(P$x, P$y, matrix(trans(P$fit + shift),
                                   n2, n2), xlab = P$xlab, ylab = P$ylab, main = P$main,
                  xlim = P$xlim, ylim = P$ylim, ...)
          if (rug) {
            if (is.null(list(...)[["pch"]]))
              points(P$raw$x, P$raw$y, pch = ".", ...)
            else points(P$raw$x, P$raw$y, ...)
          }
        }
      }
      else {
        warning("no automatic plotting for smooths of more than one variable")
      }
    }
  }
}


#' @rdname pfr_plot.gam
#' @keywords internal
plot.random.effect <- getFromNamespace("plot.random.effect", "mgcv")

