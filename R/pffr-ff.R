#' Construct a function-on-function regression term
#'
#' Defines a term \eqn{\int^{s_{hi, i}}_{s_{lo, i}} X_i(s)\beta(t,s)ds} for
#' inclusion in an \code{mgcv::gam}-formula (or \code{bam} or \code{gamm} or
#' \code{gamm4:::gamm4}) as constructed by \code{\link{pffr}}. \cr Defaults to a
#' cubic tensor product B-spline with marginal first order differences penalties
#' for \eqn{\beta(t,s)} and numerical integration over the entire range
#' \eqn{[s_{lo, i}, s_{hi, i}] = [\min(s_i), \max(s_i)]} by using Simpson
#' weights. Can't deal with any missing \eqn{X(s)}, unequal lengths of
#' \eqn{X_i(s)} not (yet?) possible. Unequal integration ranges for different
#' \eqn{X_i(s)} should work. \eqn{X_i(s)} is assumed to be numeric (duh...).
#'
#' If \code{check.ident==TRUE} and \code{basistype!="s"}  (the default), the
#' routine checks conditions for non-identifiability of the effect.  This occurs
#' if a) the marginal basis for the functional covariate is rank-deficient
#' (typically because the functional covariate has lower rank than the spline
#' basis along its index) and simultaneously b) the kernel of Cov\eqn{(X(s))} is
#' not disjunct from the kernel of the marginal penalty over \code{s}. In
#' practice, a) occurs quite frequently, and b) occurs usually because
#' curve-wise mean centering has removed all constant components from the
#' functional covariate. \cr If there is kernel overlap, \eqn{\beta(t,s)} is
#' constrained to be orthogonal to functions in that overlap space (e.g., if the
#' overlap contains constant functions, constraints "\eqn{\int \beta(t,s) ds =
#' 0} for all t" are enforced). See reference for details.\cr A warning is
#' always given if the effective rank of Cov\eqn{(X(s))} (defined as the number
#' of eigenvalues accounting for at least 0.995 of the total variance in
#' \eqn{X_i(s)}) is lower than 4. If \eqn{X_i(s)} is of very low rank,
#' \code{\link{ffpc}}-term may be preferable.
#'
#' @section Effective rank and weak identifiability:
#'
#' \code{check.ident} also compares the effective rank of Cov\eqn{(X(s))}
#' against the marginal basis dimension \eqn{k_s} used for \eqn{\beta(t,s)}
#' along \eqn{s} (i.e. \code{splinepars$k[1]}), and warns whenever the
#' effective rank is below \eqn{1.5 k_s}. The hard case \eqn{\mathrm{rank} <
#' k_s} (the surface is identified through the penalty alone) is a special case
#' of this warning and is flagged explicitly in its message.
#'
#' The reason for the margin is that identifiability is not a yes/no property
#' here: the data only inform \eqn{\beta(t,s)} in directions that the observed
#' curves span. Once \eqn{k_s} approaches the effective rank, an appreciable
#' part of a rough true surface lies outside that span, and this component is
#' returned by the penalty rather than estimated. It is a shared bias floor: it
#' affects every estimator and every kind of standard error, so both point
#' estimates and interval coverage for the \code{ff} term become hard to
#' interpret, and it does not show up as a fitting failure. A factor of about
#' 1.5 was the smallest margin at which this contamination was negligible in
#' the simulation studies behind the cluster-robust \code{\link{pffr}}
#' intervals.
#'
#' The effective rank can never exceed \code{min(nrow(X), ncol(X))}, so with
#' few curves the warning may be impossible to satisfy at the requested
#' \eqn{k_s}. Remedies, in order of preference: lower \eqn{k_s} (via
#' \code{splinepars = list(k = c(k_s, k_t))}); use more or more varied curves;
#' switch to \code{\link{ffpc}}, which parameterizes the effect in the leading
#' functional principal components of \eqn{X} and is designed for the low-rank
#' case. If none of these is possible, the fit is still usable, but conclusions
#' about \eqn{\beta(t,s)} -- including the width and coverage of its confidence
#' bands -- should be drawn with that caveat in mind. The check can be switched
#' off with \code{check.ident = FALSE}.
#'
#' @param X an n by \code{ncol(xind)} matrix of function evaluations
#'   \eqn{X_i(s_{i1}),\dots, X_i(s_{iS})}; \eqn{i=1,\dots,n}.
#' @param yind \emph{DEPRECATED} used to supply matrix (or vector) of indices of
#'   evaluations of \eqn{Y_i(t)}, no longer used.
#' @param xind vector of indices of evaluations of \eqn{X_i(s)},
#'   i.e, \eqn{(s_{1},\dots,s_{S})}
#' @param basistype defaults to "\code{\link[mgcv]{te}}", i.e. a tensor product
#'   spline to represent \eqn{\beta(t,s)}. Alternatively, use \code{"s"} for
#'   bivariate basis functions (see \code{mgcv}'s \code{\link[mgcv]{s}}) or
#'   \code{"t2"} for an alternative parameterization of tensor product splines
#'   (see \code{mgcv}'s \code{\link[mgcv]{t2}}).
#' @param integration method used for numerical integration. Defaults to
#'   \code{"simpson"}'s rule for calculating entries in \code{L}. Alternatively
#'   and for non-equidistant grids, \code{"trapezoidal"} or \code{"riemann"}.
#'   \code{"riemann"} integration is always used if \code{limits} is specified
#' @param L optional: an n by \code{ncol(xind)} matrix giving the weights for
#'   the numerical integration over \eqn{s}.
#' @param limits defaults to NULL for integration across the entire range of
#'   \eqn{X(s)}, otherwise specifies the integration limits \eqn{s_{hi}(t),
#'   s_{lo}(t)}: either one of \code{"s<t"} or \code{"s<=t"} for
#'   \eqn{(s_{hi}(t), s_{lo}(t)) = (t, 0]} or \eqn{[t, 0]}, respectively, or a
#'   function that takes \code{s} as the first and \code{t} as the second
#'   argument and returns TRUE for combinations of values \code{(s,t)} if
#'   \code{s} falls into the integration range for the given \code{t}. This is
#'   an experimental feature and not well tested yet; use at your own risk.
#' @param splinepars optional arguments supplied to the \code{basistype}-term.
#'   Defaults to a cubic tensor product B-spline with marginal first difference
#'   penalties, i.e. \code{list(bs="ps", m=list(c(2, 1), c(2,1)))}. See
#'   \code{\link[mgcv]{te}} or \code{\link[mgcv]{s}} in \pkg{mgcv} for details
#' @param check.ident check identifiability of the model spec. See Details and
#'   References. Defaults to \code{TRUE}.
#'
#' @seealso \code{mgcv}'s \code{\link[mgcv]{linear.functional.terms}}
#' @return A list containing \item{call}{a "call" to
#'   \code{\link[mgcv]{te}} (or \code{\link[mgcv]{s}} or \code{\link[mgcv]{t2}})
#'   using the appropriately constructed covariate and weight matrices}
#'   \item{data}{a list containing the necessary covariate and weight matrices}
#'
#' @author Fabian Scheipl, Sonja Greven
#' @references For background on \code{check.ident}:\cr Scheipl, F., Greven,
#'   S. (2016). Identifiability in penalized function-on-function regression
#'   models. Electronic Journal of Statistics, 10(1), 495--526.
#'   \url{https://projecteuclid.org/journals/electronic-journal-of-statistics/volume-10/issue-1/Identifiability-in-penalized-function-on-function-regression-models/10.1214/16-EJS1123.full}
#' @export
#' @importFrom MASS Null
# FIXME: weights for Simpson's rule on non-equidistant grids
# TODO: allow X to be of class fd (?)
# TODO: allow X to be a factor -- would result in one beta(s,t) surface for each level? (?)
# TODO: by variables
# TODO: add FAME penalty?
ff <- function(
  X,
  yind = NULL,
  xind = seq(0, 1, l = ncol(X)),
  basistype = c("te", "t2", "ti", "s", "tes"),
  integration = c("simpson", "trapezoidal", "riemann"),
  L = NULL,
  limits = NULL,
  splinepars = if (basistype != "s") {
    list(bs = "ps", m = list(c(2, 1), c(2, 1)), k = c(5, 5))
  } else {
    list(bs = "tp", m = NA)
  },
  check.ident = TRUE
) {
  # Deprecation warning for yind
  if (!is.null(yind)) {
    .Deprecated(
      msg = paste0(
        "The 'yind' argument in ff() is deprecated and ignored. ",
        "The y-index is now obtained automatically from pffr()."
      )
    )
  }

  n <- nrow(X)
  nxgrid <- ncol(X)

  # Validate X has no NA values
  if (anyNA(X)) {
    stop("`X` must not contain NA values.")
  }

  # Validate and expand xind to matrix form
  xind <- validate_and_expand_xind(xind, n, nxgrid, arg_name = "xind")

  basistype <- match.arg(basistype)
  integration <- match.arg(integration)

  # Check for non-equidistant grid and adjust integration method
  xind_sc <- xind - min(xind)
  xind_sc <- xind_sc / max(xind_sc)
  diff_xind <- t(round(apply(xind_sc, 1, diff), 3))

  if (
    is.null(L) &&
      any(apply(diff_xind, 1, \(x) length(unique(x))) != 1) &&
      integration == "simpson"
  ) {
    message(
      "Non-equidistant grid detected for ",
      deparse(substitute(X)),
      ".\n Changing to trapezoidal rule for integration."
    )
    integration <- "trapezoidal"
  }

  if (!is.null(limits) && integration != "riemann") {
    integration <- "riemann"
    message(
      "<limits>-argument detected. ",
      "Changing to Riemann sums for numerical integration."
    )
  }

  # Compute integration weights
  if (!is.null(L)) {
    if (nrow(L) != n || ncol(L) != nxgrid) {
      stop("`L` must be a ", n, " x ", nxgrid, " matrix.")
    }
  } else {
    L <- compute_integration_weights(xind, integration)
  }
  LX <- L * X

  # Parse limits argument
  limits <- build_limits_function(limits)

  # assign unique names based on the given args
  xindname <- paste(deparse(substitute(X)), ".smat", sep = "")
  yindname <- paste(deparse(substitute(X)), ".tmat", sep = "")
  LXname <- paste("L.", deparse(substitute(X)), sep = "")

  # make call
  splinefun <- as.symbol(basistype) # if(basistype=="te") quote(te) else quote(s)
  frmls <- if (exists(basistype, asNamespace("mgcv"), inherits = FALSE)) {
    formals(getFromNamespace(basistype, ns = "mgcv"))
  } else {
    formals(basistype)
  }
  frmls <- modifyList(frmls[names(frmls) %in% names(splinepars)], splinepars)
  call <- as.call(c(
    list(
      splinefun,
      x = as.symbol(substitute(xindname)),
      z = as.symbol(substitute(yindname)),
      by = as.symbol(substitute(LXname))
    ),
    frmls
  ))

  if (check.ident) {
    ## check whether (number of basis functions) < (number of relevant eigenfunctions of X)
    evX <- svd(X, nu = 0, nv = 0)$d^2
    maxK <- max(1, min(which((cumsum(evX) / sum(evX)) >= .995)))
    term_spec <- eval(call)
    bsdim <- if (!is.null(term_spec$margin)) {
      term_spec$margin[[1]]$bs.dim
    } else {
      term_spec$bs.dim
    }
    if (maxK <= 4)
      warning(
        "Very low effective rank of <",
        deparse(match.call()$X),
        "> detected. ",
        maxK,
        " largest eigenvalues of its covariance alone account for >99.5% of ",
        "variability. <ffpc> might be a better choice here."
      )
    ## Weak-identifiability guard. The effective rank of Cov(X(s)) has to
    ## exceed the marginal basis dimension along s comfortably, not merely
    ## match it: maxK < bsdim is the hard case (the surface is pinned down by
    ## the penalty alone), while maxK < 1.5 * bsdim is the practical margin
    ## below which part of beta(t, s) lies outside the span of the observed
    ## curves and is therefore not estimable from the data.
    if (length(bsdim) == 1 && is.finite(bsdim) && bsdim > 0) {
      if (maxK < 1.5 * bsdim) {
        warning(
          "Effective rank of <",
          deparse(match.call()$X),
          "> is ",
          maxK,
          ", below 1.5 * k = ",
          format(1.5 * bsdim),
          " for the k = ",
          bsdim,
          " basis functions along <s>",
          if (maxK < bsdim) {
            paste0(
              "; <k> is larger than the effective rank, so the model is ",
              "identifiable only through the penalty"
            )
          } else {
            ""
          },
          ". The coefficient surface is only weakly identified: components of ",
          "beta(t, s) in directions the observed curves do not span are ",
          "determined by the penalty alone, so estimates can be biased and ",
          "interval coverage unreliable in those directions. Reduce k along ",
          "<s>, or use more / richer curves -- the effective rank cannot ",
          "exceed min(nrow(X), ncol(X)) = ",
          min(n, nxgrid),
          ". See ?ff (Details) and Scheipl & Greven (2016).",
          call. = FALSE
        )
      }
    }
    if (basistype != "s") {
      # check whether span(Null(X)), span(L * B_s%*%Null(penalty)) are disjunct:
      # set up marginal spline basis:
      smConstr <- get(paste0(
        "smooth.construct.",
        attr(eval(call)$margin[[1]], "class")
      ))
      basisdata <- list(sort(unique(xind)))
      names(basisdata) <- xindname
      basis <- smConstr(
        object = list(
          term = xindname,
          bs.dim = ifelse(!is.null(call$k[1]), call$k[1], -1),
          fixed = FALSE,
          dim = 1,
          p.order = if (!is.null(call$m)) call$m[[1]] else NA,
          by = NA
        ),
        data = basisdata,
        knots = list()
      )

      # get condition number of marginal design matrix
      evDs <- svd(LX %*% basis$X, nu = 0, nv = 0)$d^2
      logCondDs <- log10(max(evDs)) - log10(min(evDs))

      N.X <- Null(t(X))
      ## for the artificial examples in the paper below produces surprising
      ## results if integration != "riemann":
      ## -- usually more constraints than expected, e.g. constraints on "linearish"
      ## even though nullspace of X only contains constants by construction --
      ## unless diag(L[1,]) is rm'ed from N.pen: non-constant integration wts
      ## seem to implicate higher order eigenfunctions in the kernel as well, e.g.
      ## sv's of N.pen have high frequency oscillations, etc (...waves hands...)
      N.pen <- diag(L[1, ]) %*% basis$X %*% Null(basis$S[[1]])
      if (any(c(NCOL(N.X) == 0, NCOL(N.pen) == 0))) {
        nullOverlap <- 0
      } else {
        nullOverlap <- trace_lv(svd(N.X)$u, svd(N.pen)$u)
      }
      if (nullOverlap > 0.95 & logCondDs > 6) {
        warning(
          "Found badly conditioned design matrix for the functional effect",
          " and kernel overlap for <",
          deparse(match.call()$X),
          "> and the specified basis and penalty. ",
          "Enforcing constraint to force function components in this overlap to 0 ",
          "since coefficient surface is not identifiable in that function space.",
          "See Scheipl/Greven (2016) for details & alternatives."
        )

        C_overlap <- {
          tmp <- svd(qr.fitted(qr(N.X), N.pen))
          t(tmp$u[,
            which(tmp$d > max(tmp$d) * .Machine$double.eps^.66),
            drop = FALSE
          ]) %*%
            basis$X
        }
        if (is.null(call$xt)) {
          call$xt <- list(C1 = C_overlap)
        } else {
          call$xt <- c(call$xt, C1 = C_overlap)
        }

        call$bs <- c(
          "ps_c",
          sub(
            ".smooth.spec",
            "",
            attr(eval(call)$margin[[2]], "class"),
            fixed = TRUE
          )
        )
        call[[1]] <- as.symbol("ti")
        call$mc <- c(TRUE, FALSE)
      } else {
      }
    }
  }
  return(list(
    call = call,
    xind = xind[1, ],
    LX = LX,
    L = L,
    xindname = xindname,
    yindname = yindname,
    LXname = LXname,
    limits = limits
  ))
} #end ff()
