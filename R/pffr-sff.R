#' Construct a smooth function-on-function regression term
#'
#' Defines a term \eqn{\int^{s_{hi, i}}_{s_{lo, i}} f(X_i(s), s, t) ds} for
#' inclusion in an \code{mgcv::gam}-formula (or \code{bam} or \code{gamm} or
#' \code{gamm4:::gamm}) as constructed by \code{\link{pffr}}. Defaults to a
#' cubic tensor product B-spline with marginal second differences penalties for
#' \eqn{f(X_i(s), s, t)} and integration over the entire range \eqn{[s_{lo, i},
#' s_{hi, i}] = [\min(s_i), \max(s_i)]}. Can't deal with any missing \eqn{X(s)},
#' unequal lengths of \eqn{X_i(s)} not (yet?) possible. Unequal ranges for
#' different \eqn{X_i(s)} should work. \eqn{X_i(s)} is assumed to be numeric.\cr
#' \code{sff()} IS AN EXPERIMENTAL FEATURE AND NOT WELL TESTED YET -- USE AT
#' YOUR OWN RISK.
#'
#' @param X an n by \code{ncol(xind)} matrix of function evaluations
#'   \eqn{X_i(s_{i1}),\dots, X_i(s_{iS})}; \eqn{i=1,\dots,n}.
#' @param yind \emph{DEPRECATED} matrix (or vector) of indices of evaluations of
#'   \eqn{Y_i(t)}; i.e. matrix with rows \eqn{(t_{i1},\dots,t_{iT})}; no longer
#'   used.
#' @param xind vector of indices of evaluations of \eqn{X_i(s)},
#'   i.e, \eqn{(s_{1},\dots,s_{S})}
#' @param basistype defaults to "\code{\link[mgcv]{te}}", i.e. a tensor product
#'   spline to represent \eqn{f(X_i(s), t)}. Alternatively, use \code{"s"} for
#'   bivariate basis functions (see \code{\link[mgcv]{s}}) or \code{"t2"} for an
#'   alternative parameterization of tensor product splines (see
#'   \code{\link[mgcv]{t2}}).
#' @param integration method used for numerical integration. Defaults to
#'   \code{"simpson"}'s rule. Alternatively and for non-equidistant grids,
#'   \code{"trapezoidal"}.
#' @param L optional: an n by \code{ncol(xind)} giving the weights for the
#'   numerical integration over \eqn{s}.
#' @param limits defaults to NULL for integration across the entire range of
#'   \eqn{X(s)}, otherwise specifies the integration limits \eqn{s_{hi, i},
#'   s_{lo, i}}: either one of \code{"s<t"} or \code{"s<=t"} for \eqn{(s_{hi,
#'   i}, s_{lo, i}) = (0, t)} or a function that takes \code{s} as the first and
#'   \code{t} as the second argument and returns TRUE for combinations of values
#'   \code{(s,t)} if \code{s} falls into the integration range for the given
#'   \code{t}. This is an experimental feature and not well tested yet; use at
#'   your own risk.
#' @param splinepars optional arguments supplied to the \code{basistype}-term.
#'   Defaults to a cubic tensor product B-spline with marginal second
#'   differences, i.e. \code{list(bs="ps", m=c(2,2,2))}. See
#'   \code{\link[mgcv]{te}} or \code{\link[mgcv]{s}} for details
#'
#' @return a list containing \itemize{ \item \code{call} a "call" to
#'   \code{\link[mgcv]{te}} (or \code{\link[mgcv]{s}}, \code{\link[mgcv]{t2}})
#'   using the appropriately constructed covariate and weight matrices (see
#'   \code{\link[mgcv]{linear.functional.terms}}) \item \code{data} a list
#'   containing the necessary covariate and weight matrices }
#'
#' @author Fabian Scheipl, based on Sonja Greven's trick for fitting functional
#'   responses.
#' @export
#' @importFrom utils getFromNamespace modifyList
# FIXME: figure out weights for simpson's rule on non-equidistant grids
# TODO: allow X to be of class fd (?)
# TODO: by variables (?)
sff <- function(
  X,
  yind = NULL,
  xind = seq(0, 1, length.out = ncol(X)),
  basistype = c("te", "t2", "s"),
  integration = c("simpson", "trapezoidal"),
  L = NULL,
  limits = NULL,
  splinepars = list(bs = "ps", m = c(2, 2, 2))
) {
  # Deprecation warning for yind
  if (!is.null(yind)) {
    .Deprecated(
      msg = paste0(
        "The 'yind' argument in sff() is deprecated and ignored. ",
        "The y-index is now obtained automatically from pffr()."
      )
    )
  }

  n <- nrow(X)
  nxgrid <- ncol(X)

  # Validate X has no NA values
  if (anyNA(X)) {
    cli::cli_abort("{.arg X} must not contain NA values.")
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
    warning(
      "Non-equidistant grid detected for ",
      deparse(substitute(X)),
      ".\n Changing to trapezoidal rule for integration."
    )
    integration <- "trapezoidal"
  }

  # Compute integration weights
  if (!is.null(L)) {
    if (nrow(L) != n || ncol(L) != nxgrid) {
      cli::cli_abort("{.arg L} must be a {n} x {nxgrid} matrix.")
    }
  } else {
    L <- compute_integration_weights(xind, integration)
  }

  # Parse limits argument
  limits <- build_limits_function(limits)
  #assign unique names based on the given args
  xname <- paste(deparse(substitute(X)), ".mat", sep = "")
  xindname <- paste(deparse(substitute(X)), ".smat", sep = "")
  yindname <- paste(deparse(substitute(X)), ".tmat", sep = "")
  LXname <- paste("L.", deparse(substitute(X)), sep = "")

  # make call
  splinefun <- as.symbol(basistype) # if(basistype=="te") quote(te) else quote(s)
  frmls <- formals(getFromNamespace(deparse(splinefun), ns = "mgcv"))
  frmls <- modifyList(frmls[names(frmls) %in% names(splinepars)], splinepars)
  call <- as.call(c(
    list(
      splinefun,
      x = as.symbol(substitute(yindname)),
      y = as.symbol(substitute(xindname)),
      z = as.symbol(substitute(xname)),
      by = as.symbol(substitute(LXname))
    ),
    frmls
  ))

  return(list(
    call = call,
    xind = xind[1, ],
    L = L,
    X = X,
    xname = xname,
    xindname = xindname,
    yindname = yindname,
    LXname = LXname
  ))
} #end sff()
