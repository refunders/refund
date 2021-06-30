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
ff <- function(X,
  yind = NULL,
  xind=seq(0, 1, l=ncol(X)),
  basistype= c("te", "t2", "ti", "s", "tes"),
  integration=c("simpson", "trapezoidal", "riemann"),
  L=NULL,
  limits=NULL,
  splinepars=if(basistype != "s") {
    list(bs="ps", m=list(c(2, 1), c(2,1)), k=c(5, 5))
  } else {
    list(bs="tp", m=NA)
  },
  check.ident=TRUE
){
  n <- nrow(X)
  nxgrid <- ncol(X)
  stopifnot(all(!is.na(X)))


  # check & format index for X
  if(is.null(dim(xind))){
    xind <- t(as.matrix(xind))
  }
  stopifnot(ncol(xind) == nxgrid)
  if(nrow(xind)== 1){
    xind <- matrix(as.vector(xind), nrow=n, ncol=nxgrid, byrow=T)
  } else {
    stop("<xind> has to be supplied as a vector or matrix with a single row.")
  }
  stopifnot(nrow(xind) == n)
  stopifnot(all.equal(order(xind[1,]), 1:nxgrid))

  basistype <- match.arg(basistype)
  integration <- match.arg(integration)

  # scale xind to [0, 1] and check for reasonably equidistant gridpoints
  xind.sc <- xind - min(xind)
  xind.sc <- xind.sc/max(xind.sc)
  diffXind <- t(round(apply(xind.sc, 1, diff), 3))
  if(is.null(L) & any(apply(diffXind, 1, function(x) length(unique(x))) != 1) &&
      # gridpoints for any  X_i(s) not equidistant?
      integration=="simpson"){
    message("Non-equidistant grid detected for ", deparse(substitute(X)),
      ".\n Changing to trapezoidal rule for integration.")
    integration <- "trapezoidal"
  }
  if(!is.null(limits) && integration != "riemann"){
    integration <- "riemann"
    message("<limits>-argument detected. ",
      "Changing to Riemann sums for numerical integration.")
  }
  # FIXME: figure out weights for simpson's rule on non-equidistant grids instead of all this...

  #make weight matrix for by-term
  if(!is.null(L)){
    stopifnot(nrow(L) == n, ncol(L) == nxgrid)
    #TODO: check whether supplied L is compatibel with limits argument
  } else {

    L <- switch(integration,
      "simpson" = {
        # \int^b_a f(t) dt = (b-a)/gridlength/3 * [f(a) + 4*f(t_1) + 2*f(t_2) + 4*f(t_3) + 2*f(t_3) +...+ f(b)]
        ((xind[,nxgrid]-xind[,1])/nxgrid)/3 *
          matrix(c(1, rep(c(4, 2), length=nxgrid-2), 1), nrow=n, ncol=nxgrid, byrow=T)
      },
      "trapezoidal" = {
        # \int^b_a f(t) dt = .5* sum_i (t_i - t_{i-1}) f(t_i) + f(t_{i-1}) =
        #	(t_2 - t_1)/2 * f(a=t_1) + sum^{nx-1}_{i=2} ((t_i - t_i-1)/2 + (t_i+1 - t_i)/2) * f(t_i) + ... +
        #			+ (t_nx - t_{nx-1})/2 * f(b=t_n)
        diffs <- t(apply(xind, 1, diff))
        .5 * cbind(diffs[,1],
          t(apply(diffs, 1, filter, filter=c(1,1)))[,-(nxgrid-1)],
          diffs[,(nxgrid-1)])
      },
      "riemann" = {
        # simple quadrature rule:
        # \int^b_a f(t) dt = sum_i (t_i-t_{i-1})*(f(t_i))
        diffs <- t(apply(xind, 1, diff))
        #assume delta(t_0=a, t_1) = avg. delta
        cbind(rep(mean(diffs),n), diffs)
      }
    )

  }
  LX <- L*X

  if(!is.null(limits)){
    if(!is.function(limits)){
      if(!(limits %in% c("s<t","s<=t"))){
        stop("supplied <limits> argument unknown")
      }
      if(limits=="s<t"){
        limits <- function(s, t){
          s < t
        }
      } else {
        if(limits=="s<=t"){
          limits <- function(s, t){
            (s < t) | (s == t)
          }
        }
      }
    }
  }


  # assign unique names based on the given args
  xindname <- paste(deparse(substitute(X)), ".smat", sep="")
  yindname <- paste(deparse(substitute(X)), ".tmat", sep="")
  LXname <- paste("L.", deparse(substitute(X)), sep="")

  # make call
  splinefun <- as.symbol(basistype) # if(basistype=="te") quote(te) else quote(s)
  frmls <- if(exists(basistype, asNamespace("mgcv"),  inherits = FALSE)) {
    formals(getFromNamespace(basistype, ns="mgcv"))
  } else {
    formals(basistype)
  }
  frmls <- modifyList(frmls[names(frmls) %in% names(splinepars)], splinepars)
  call <- as.call(c(
    list(splinefun,
      x = as.symbol(substitute(xindname)),
      z = as.symbol(substitute(yindname)),
      by =as.symbol(substitute(LXname))),
    frmls))

  if(check.ident){
    ## check whether (number of basis functions) < (number of relevant eigenfunctions of X)
    evX <- svd(X, nu = 0, nv = 0)$d^2
    maxK <- max(1, min(which((cumsum( evX)/sum(evX)) >= .995)))
    bsdim <- eval(call)$margin[[1]]$bs.dim
    if(maxK <= 4)
      warning("Very low effective rank of <", deparse(match.call()$X),
        "> detected. ", maxK,
        " largest eigenvalues of its covariance alone account for >99.5% of ",
        "variability. <ffpc> might be a better choice here.")
    if(maxK < bsdim){
      warning("<k> larger than effective rank of <",deparse(match.call()$X),
        ">. Model identifiable only through penalty.")
    }
    if(basistype!="s"){
      # check whether span(Null(X)), span(L * B_s%*%Null(penalty)) are disjunct:
      # set up marginal spline basis:
      smConstr <- get(paste0("smooth.construct.",
        attr(eval(call)$margin[[1]], "class")))
      basisdata <- list(sort(unique(xind)))
      names(basisdata) <- xindname
      basis <- smConstr(object=list(term=xindname,
          bs.dim=ifelse(!is.null(call$k[1]), call$k[1], -1),
          fixed=FALSE, dim=1,
          p.order=if(!is.null(call$m)) call$m[[1]] else NA,
          by=NA),
        data=basisdata, knots=list())

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
      N.pen <- diag(L[1,]) %*% basis$X %*% Null(basis$S[[1]])
      if(any(c(NCOL(N.X) == 0, NCOL(N.pen) == 0))) {
        nullOverlap <- 0
      } else {
        nullOverlap <- trace_lv(svd(N.X)$u, svd(N.pen)$u)
      }
      if(nullOverlap > 0.95 & logCondDs > 6){
        warning("Found badly conditioned design matrix for the functional effect",
          " and kernel overlap for <", deparse(match.call()$X),
          "> and the specified basis and penalty. ",
          "Enforcing constraint to force function components in this overlap to 0 ",
          "since coefficient surface is not identifiable in that function space.",
          "See Scheipl/Greven (2016) for details & alternatives.")

        C_overlap <- {
          tmp <- svd(qr.fitted(qr(N.X), N.pen))
          t(tmp$u[, which(tmp$d > max(tmp$d)*.Machine$double.eps^.66), drop=FALSE]) %*%
            basis$X
        }
        if(is.null(call$xt)) {
          call$xt <- list(C1 = C_overlap)
        } else {
          call$xt <- c(call$xt, C1 = C_overlap)
        }

        call$bs <- c("ps_c",
          sub(".smooth.spec", "", attr(eval(call)$margin[[2]], "class"), fixed=TRUE))
        call[[1]] <- as.symbol("ti")
        call$mc <- c(TRUE, FALSE)
      } else {

      }
    }
  }
  return(list(call=call, xind=xind[1,], LX=LX, L=L,
    xindname=xindname, yindname=yindname,
    LXname=LXname, limits=limits))
}#end ff()
