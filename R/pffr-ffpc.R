#' Construct a PC-based function-on-function regression term
#'
#' Defines a term \eqn{\int X_i(s)\beta(t,s)ds}
#' for inclusion in an \code{mgcv::gam}-formula (or \code{bam} or \code{gamm} or \code{gamm4:::gamm4}) as constructed
#' by \code{\link{pffr}}.
#'
#' In contrast to \code{\link{ff}}, \code{ffpc}
#' does an FPCA decomposition \eqn{X(s) \approx \sum^K_{k=1} \xi_{ik} \Phi_k(s)} using \code{\link{fpca.sc}} and
#' represents \eqn{\beta(t,s)} in the function space spanned by these \eqn{\Phi_k(s)}.
#' That is, since
#' \deqn{\int X_i(s)\beta(t,s)ds = \sum^K_{k=1} \xi_{ik} \int \Phi_k(s) \beta(s,t) ds = \sum^K_{k=1} \xi_{ik} \tilde \beta_k(t),}
#' the function-on-function term can be represented as a sum of \eqn{K} univariate functions \eqn{\tilde \beta_k(t)} in \eqn{t} each multiplied by the FPC
#' scores \eqn{\xi_{ik}}. The truncation parameter \eqn{K} is chosen as described in \code{\link{fpca.sc}}.
#' Using this instead of \code{ff()} can be beneficial if the covariance operator of the \eqn{X_i(s)}
#' has low effective rank (i.e., if \eqn{K} is small). If the covariance operator of the \eqn{X_i(s)}
#' is of (very) high rank, i.e., if \eqn{K} is large, \code{ffpc()} will not be very efficient.
#'
#' To reduce model complexity, the \eqn{\tilde \beta_k(t)} all have a single joint smoothing parameter
#' (in \code{mgcv}, they get the same \code{id}, see \code{\link[mgcv]{s}}).\cr
#'
#' Please see \code{\link[refund]{pffr}} for details on model specification and
#' implementation.
#'
#' @param X an n by \code{ncol(xind)} matrix of function evaluations \eqn{X_i(s_{i1}),\dots, X_i(s_{iS})}; \eqn{i=1,\dots,n}.
#' @param yind \emph{DEPRECATED} used to supply matrix (or vector) of indices of evaluations of \eqn{Y_i(t)}, no longer used.
#' @param xind matrix (or vector) of indices of evaluations of \eqn{X_i(t)}, defaults to \code{seq(0, 1, length=ncol(X))}.
#' @param splinepars optional arguments supplied to the \code{basistype}-term. Defaults to a cubic
#' 	B-spline with first difference penalties and 8 basis functions for each \eqn{\tilde \beta_k(t)}.
#' @param decomppars  parameters for the FPCA performed with \code{\link{fpca.sc}}.
#'   Unless they include \code{argvals}, the FPCA uses \code{argvals = xind},
#'   so the FPCs are orthonormal and the scores are integrals over the domain
#'   of \code{xind}. The coefficient surface implied by the fit is returned by
#'   \code{\link{coef.pffr}} and plotted by \code{\link{ffpcplot}}.
#' @param npc.max maximal number \eqn{K} of FPCs to use, regardless of \code{decomppars}; defaults to 15
#' @return A list containing the necessary information to construct a term to be included in a \code{mgcv::gam}-formula.
#'
#' @author Fabian Scheipl
#' @export
#' @examples \dontrun{
#' set.seed(1122)
#' n <- 55
#' S <- 60
#' T <- 50
#' s <- seq(0,1, l=S)
#' t <- seq(0,1, l=T)
#'
#' #generate X from a polynomial FPC-basis:
#' rankX <- 5
#' Phi <- cbind(1/sqrt(S), poly(s, degree=rankX-1))
#' lambda <- rankX:1
#' Xi <- sapply(lambda, function(l)
#'             scale(rnorm(n, sd=sqrt(l)), scale=FALSE))
#' X <- Xi %*% t(Phi)
#'
#' beta.st <- outer(s, t, function(s, t) cos(2 * pi * s * t))
#'
#' y <- (1/S*X) %*% beta.st + 0.1 * matrix(rnorm(n * T), nrow=n, ncol=T)
#'
#' data <- list(y=y, X=X)
#' # set number of FPCs to true rank of process for this example:
#' m.pc <- pffr(y ~ c(1) + 0 + ffpc(X, yind=t, decomppars=list(npc=rankX)),
#'         data=data, yind=t)
#' summary(m.pc)
#' m.ff <- pffr(y ~ c(1) + 0 + ff(X, yind=t), data=data, yind=t)
#' summary(m.ff)
#'
#' # fits are very similar:
#' all.equal(fitted(m.pc), fitted(m.ff))
#'
#' # plot implied coefficient surfaces:
#' layout(t(1:3))
#' persp(t, s, t(beta.st), theta=50, phi=40, main="Truth",
#'     ticktype="detailed")
#' plot(m.ff, select=1, zlim=range(beta.st), theta=50, phi=40,
#'     ticktype="detailed")
#' title(main="ff()")
#' ffpcplot(m.pc, type="surf", auto.layout=FALSE, theta = 50, phi = 40)
#' title(main="ffpc()")
#'
#' # show default ffpcplot:
#' ffpcplot(m.pc)
#' }
ffpc <- function(
  X,
  yind = NULL,
  xind = seq(0, 1, length = ncol(X)),
  splinepars = list(bs = "ps", m = c(2, 1), k = 8),
  decomppars = list(pve = .99, useSymm = TRUE),
  npc.max = 15
) {
  # Deprecation warning for yind
  if (!is.null(yind)) {
    .Deprecated(
      msg = paste0(
        "The 'yind' argument in ffpc() is deprecated and ignored. ",
        "The y-index is now obtained automatically from pffr()."
      )
    )
  }

  nxgrid <- length(xind)

  # check & format index for X
  stopifnot(length(xind) == ncol(X))
  stopifnot(all.equal(order(xind), 1:nxgrid))

  # Run the FPCA on the covariate's own index, so that the eigenfunctions are
  # orthonormal and the scores are integrals over the domain of xind.
  if (is.null(decomppars$argvals) && is.null(decomppars$ydata)) {
    decomppars$argvals <- xind
  }
  decomppars$Y <- X
  klX <- do.call(fpca.sc, decomppars)
  npc <- min(ncol(klX$scores), npc.max)
  xiMat <- klX$scores[, 1:npc, drop = FALSE]

  #assign unique names based on the given args
  colnames(xiMat) <- paste(
    make.names(deparse(substitute(X))),
    ".PC",
    1:ncol(xiMat),
    sep = ""
  )
  id <- paste(make.names(deparse(substitute(X))), ".ffpc", sep = "")
  return(list(
    data = xiMat,
    PCMat = klX$efunctions[, 1:npc, drop = FALSE],
    meanX = klX$mu,
    eigenvalues = klX$evalues,
    xind = xind,
    id = id,
    splinepars = splinepars,
    # grid on which fpca.sc normalised the eigenfunctions
    argvals = klX$argvals,
    # what fpca.sc's score BLUPs need, to compute scores for new data
    score_pars = list(
      efunctions = klX$efunctions,
      evalues = klX$evalues,
      sigma2 = klX$sigma2
    )
  ))
} #end ffpc()

# Helpers for fitted ffpc terms ----------------------------------------------

#' Scale factor from FPCs to the coefficient surface of an ffpc term
#'
#' \code{fpca.sc()} normalises the eigenfunctions \eqn{\phi_k} on its grid
#' \code{argvals}, and its scores are \eqn{\xi_k = \int X^c(u)\phi_k(u)du} on that
#' grid. With \eqn{\psi_k = c\,\phi_k} and \eqn{c} the ratio of the domain
#' lengths of \code{argvals} and \code{xind}, \eqn{\int X^c(s)\psi_k(s)ds = \xi_k} over
#' the domain of \code{xind}, so the implied surface is
#' \eqn{\beta(s,t) = c\sum_k \phi_k(s)\tilde\beta_k(t)}. \code{c = 1} for terms
#' from this version of \code{\link{ffpc}}, which runs the FPCA on \code{xind}; earlier
#' versions ran it on \code{seq(0, 1)}.
#'
#' @param trm An element of \code{object$pffr$ffpc}.
#' @returns A positive scalar.
#' @keywords internal
ffpc_beta_scale <- function(trm) {
  argvals <- trm$argvals %||% c(0, 1)
  diff(range(argvals)) / diff(range(trm$xind))
}

#' FPC scores of new covariate curves for an ffpc term
#'
#' Computes the scores as \code{fpca.sc()} computes them for the curves the term was
#' fitted on: the BLUPs \eqn{(Z^\top Z + \sigma^2 \Lambda^{-1})^{-1} Z^\top
#' (x - \mu)} over the observed points of each curve, using all estimated FPCs,
#' then truncated to the FPCs in the model. Terms from earlier versions of
#' \code{\link{ffpc}} lack \eqn{\sigma^2} and fall back to least-squares projection on
#' the FPCs in the model.
#'
#' @param trm An element of \code{object$pffr$ffpc}.
#' @param X Matrix of covariate curves evaluated on \code{trm$xind}.
#' @returns Matrix of scores, one row per curve.
#' @keywords internal
ffpc_scores <- function(trm, X) {
  X <- as.matrix(X)
  npc <- ncol(trm$PCMat)
  Xc <- sweep(X, 2, as.vector(trm$meanX))
  pars <- trm$score_pars
  if (is.null(pars)) {
    return(t(qr.coef(qr(trm$PCMat), t(Xc))))
  }
  Z <- pars$efunctions
  D_inv <- diag(1 / pars$evalues, nrow = ncol(Z))
  blup <- function(x) {
    obs <- which(!is.na(x))
    Zo <- Z[obs, , drop = FALSE]
    solve(crossprod(Zo) + pars$sigma2 * D_inv, crossprod(Zo, x[obs]))
  }
  scores <- matrix(
    unlist(lapply(seq_len(nrow(Xc)), \(i) blup(Xc[i, ]))),
    nrow = nrow(Xc),
    byrow = TRUE
  )
  scores[, seq_len(npc), drop = FALSE]
}

#' Smooth indices of the FPC-specific coefficient functions of ffpc terms
#'
#' @param object A fitted \code{pffr} object.
#' @returns A list with one integer vector per ffpc term: the indices in
#'   \code{object$smooth} of \eqn{\tilde\beta_1(t), \dots, \tilde\beta_K(t)}, in
#'   the order of the FPCs.
#' @keywords internal
ffpc_smooth_indices <- function(object) {
  lapply(names(object$pffr$ffpc), \(nm) {
    idx <- match(object$pffr$label_map[[nm]], names(object$smooth))
    pc <- vapply(
      object$smooth[idx],
      \(sm) as.integer(sub(".*\\.PC([0-9]+)$", "\\1", sm$by)),
      integer(1)
    )
    idx <- idx[order(pc)]
    stopifnot(
      !anyNA(idx),
      length(idx) == ncol(object$pffr$ffpc[[nm]]$PCMat)
    )
    idx
  })
}

#' Linear map from the model coefficients to the surface of an ffpc term
#'
#' Rows are the points of the grid \code{s} x \code{t} (\code{s} varying fastest), columns
#' the model coefficients, so that \code{L \%*\% coef(object, raw = TRUE)} is the
#' surface \eqn{\beta(s,t) = \sum_k \psi_k(s)\tilde\beta_k(t)} at these points
#' (see \code{\link{ffpc_beta_scale}}). The estimated FPCs are treated as fixed.
#'
#' @param object A fitted \code{pffr} object.
#' @param which Index of the ffpc term in \code{object$pffr$ffpc}.
#' @param t Evaluation points along the response index.
#' @returns A list with the matrix \code{L} and the grid values \code{s} and \code{t}.
#' @keywords internal
ffpc_surface_map <- function(object, which, t) {
  trm <- object$pffr$ffpc[[which]]
  idx <- ffpc_smooth_indices(object)[[which]]
  psi <- ffpc_beta_scale(trm) * trm$PCMat
  L <- matrix(0, length(trm$xind) * length(t), length(object$coefficients))
  for (k in seq_along(idx)) {
    sm <- object$smooth[[idx[k]]]
    nd <- data.frame(t, 1)
    names(nd) <- c(sm$term, sm$by)
    B <- mgcv::PredictMat(sm, nd)
    L[, sm$first.para:sm$last.para] <- kronecker(B, psi[, k, drop = FALSE])
  }
  list(L = L, s = trm$xind, t = t)
}

#' Plot PC-based function-on-function regression terms
#'
#' Convenience function for graphical summaries of \code{ffpc}-terms from a
#' \code{pffr} fit.
#'
#' @param object a fitted \code{pffr}-model
#' @param type one of "fpc+surf", "surf" or "fpc": "surf" shows a perspective plot of the coefficient surface implied
#'          by the estimated effect functions of the FPC scores, "fpc" shows three plots:
#'          1) a scree-type plot of the estimated eigenvalues of the functional covariate, 2) the estimated eigenfunctions,
#'          and 3) the estimated coefficient functions associated with the FPC scores. Defaults to showing both.
#' @param se.mult display estimated coefficient functions associated with the FPC scores with plus/minus this number time the estimated standard error.
#' Defaults to 2.
#' @param pages  the number of pages over which to spread the output. Defaults to 1. (Irrelevant if \code{auto.layout=FALSE}.)
#' @param ticktype see \code{\link[graphics]{persp}}.
#' @param theta see \code{\link[graphics]{persp}}.
#' @param phi see \code{\link[graphics]{persp}}.
#' @param plot produce plots or only return plotting data? Defaults to \code{TRUE}.
#' @param auto.layout should the the function set a suitable layout automatically? Defaults to TRUE
#' @return primarily produces plots, invisibly returns a list containing
#' the data used for the plots: \code{betatilde} and \code{betatilde.se}, the
#' estimated coefficient functions of the FPC scores on the response index and
#' their standard errors (from \code{object$Vp}), and \code{phibeta}, the
#' implied coefficient surfaces \eqn{\beta(s,t)} (rows: the covariate's index
#' \code{xind}, columns: the response index), as returned by
#' \code{\link{coef.pffr}}.
#'
#' @author Fabian Scheipl
#' @export
#' @importFrom graphics persp layout polygon matplot
#' @importFrom mgcv gam
#' @examples \dontrun{
#'  #see ?ffpc
#' }
ffpcplot <- function(
  object,
  type = c("fpc+surf", "surf", "fpc"),
  pages = 1,
  se.mult = 2,
  ticktype = "detailed",
  theta = 30,
  phi = 30,
  plot = TRUE,
  auto.layout = TRUE
) {
  type <- match.arg(type)
  T <- object$pffr$nyindex
  nterms <- length(object$pffr$ffpc)
  ffpcnames <- names(object$pffr$ffpc)

  # one row per response index value, every FPC score set to 1
  betadata <- object$model[rep(1, T), ]
  betadata[, paste0(object$pffr$yind_name, ".vec")] <- object$pffr$yind
  betadata[, grep(".PC[[:digit:]]+$", colnames(betadata))] <- 1
  termsffpc <- predict.gam(
    object,
    newdata = betadata,
    type = "iterms",
    se.fit = TRUE
  )

  betatilde <- termsffpc$fit[,
    grep(".PC[[:digit:]]+$", colnames(termsffpc$fit)),
    drop = FALSE
  ]
  betatilde.se <- termsffpc$se.fit[,
    grep(".PC[[:digit:]]+$", colnames(termsffpc$se.fit)),
    drop = FALSE
  ]
  betatilde.up <- betatilde + se.mult * betatilde.se
  betatilde.lo <- betatilde - se.mult * betatilde.se
  betatildemap <- lapply(
    ffpcnames,
    function(n) which(colnames(betatilde) %in% object$pffr$label_map[[n]])
  )

  phibeta <- vector(length = nterms, mode = "list")
  names(phibeta) <- sapply(object$pffr$ffpc, "[[", "id")
  for (i in 1:nterms) {
    ### betatilde_k(t) = int phi_k(s) beta(s,t) ds for orthonormal phi_k, so
    ### beta(s,t) = sum_k phi_k(s) betatilde_k(t) (see ffpc_beta_scale()).
    trm <- object$pffr$ffpc[[i]]
    phibeta[[i]] <- ffpc_beta_scale(trm) *
      trm$PCMat %*%
        t(betatilde[, betatildemap[[i]], drop = FALSE])
  }

  if (plot) {
    if (auto.layout) {
      nplots <- switch(
        type,
        "surf" = nterms,
        "fpc" = 3 * nterms,
        "fpc+surf" = 4 * nterms
      )

      #define layout
      plotsperpage <- ceiling(nplots / pages)
      columns <- switch(type, "surf" = plotsperpage, "fpc" = 3, "fpc+surf" = 4)
      layout(matrix(
        1:plotsperpage,
        ncol = columns,
        nrow = ceiling(plotsperpage / columns),
        byrow = TRUE
      ))
    }

    for (i in 1:nterms) {
      trm <- object$pffr$ffpc[[i]]
      if (type == "fpc+surf" | type == "fpc") {
        #
        npc <- ncol(trm$PCMat)
        plot(
          1:npc,
          trm$eigenvalues[1:npc],
          col = 1:npc,
          xlab = "FPC",
          type = "b",
          pch = 19,
          ylab = paste0("Estimated eigenvalues: ", trm$id),
          bty = "n"
        )
        matplot(
          trm$xind,
          trm$PCMat[, npc:1],
          type = "l",
          lty = 1,
          col = npc:1,
          xlab = "",
          ylab = paste0("Estimated FPCs: ", trm$id),
          bty = "n"
        )
        matplot(
          object$pffr$yind,
          betatilde[, rev(betatildemap[[i]])],
          type = "l",
          lty = 1,
          ylim = range(
            betatilde.up[, betatildemap[[i]]],
            betatilde.lo[, betatildemap[[i]]]
          ),
          col = npc:1,
          xlab = "",
          ylab = paste0("Effects of FPC scores"),
          bty = "n"
        )
        abline(h = 0, col = "grey", lwd = .5)
        secol <- length(betatildemap[[i]])
        for (j in rev(betatildemap[[i]])) {
          polygon(
            cbind(x = c(object$pffr$yind, rev(object$pffr$yind))),
            y = c(betatilde.up[, j], rev(betatilde.lo[, j])),
            col = do.call(rgb, as.list(c(col2rgb(secol) / 255, .1))),
            border = NA
          )
          secol <- secol - 1
        }
      }
      if (type == "fpc+surf" | type == "surf") {
        persp(
          trm$xind,
          object$pffr$yind,
          z = phibeta[[i]],
          theta = theta,
          phi = phi,
          ticktype = ticktype,
          xlab = "x.index",
          ylab = "y.index",
          zlim = range(as.vector(phibeta[[i]])),
          zlab = trm$id
        )
      }
    } #end for(i)
  }
  invisible(list(
    betatilde = betatilde,
    betatilde.se = betatilde.se,
    phibeta = phibeta
  ))
} #end ffpcplot()
