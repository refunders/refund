#' Construct a VDFR regression term
#' 
#' This function defines the a variable-domain functional regression term
#' for inclusion in an \code{mgcv::gam-formula} (or \code{bam} or
#' \code{gamm} or \code{gamm4:::gamm}) as constructed by \code{\link{pfr}}.
#' These are functional predictors for which each function is observed over a
#' domain of different width.
#' The default is the ``untransformed" term,\eqn{1/T_i\int_0^{T_i}X_i(t)\beta(t,T_i)dt},
#' where \eqn{X_i(t)} is a functional predictor of length \eqn{T_i} and \eqn{\beta(t,T_i)}
#' is an unknown bivariate coefficient function. Lagged-domain and
#' standardized-domain models are also allowed. For the standardized-domain
#' model, the interaction between \eqn{t} and \eqn{T_i} could be nonparametric, linear,
#' quadratic, or not present at all. Basis choice is fully custiomizable using
#' the options of \code{mgcv::s} and \code{mgcv::te}, though tensor-product
#' bases are not allowed except in the standardized-domain case.
#' 
#' @param X matrix containing variable-domain functions. Should be \eqn{N x J},
#'    where \eqn{N} is the number of subjects and \eqn{J} is the maximum number of time
#'    points per subject. Most rows will have \code{NA} values in the right-most
#'    columns, corresponding to unobserved time points.
#' @param argvals matrix (or vector) containing the time indices of evaluations of
#'    \eqn{X_i(t)}. If a matrix, it must be the same dimensionality as \code{X}; if a
#'    vector, must be of length \code{ncol(X)}.
#' @param Tind vector of values of \eqn{T_i}. Defaults to the \code{argvals} value
#'    corresponding to the last observation of \eqn{X_i(t)}.
#' @param T.trans optional function applied to \code{Tind} to allow the interaction
#'    to occur on a transformed scale, e.g. the log or quantile scale.
#' @param domain defines the domain for each function \eqn{X_i(t)}; see Details.
#' @param interaction defines the type of interaction between \eqn{t} and \eqn{T_i};
#'    see Details. Must be nonparametric if \code{domain} is not ``standardized."
#' @param integration method used for numerical integration. Defaults to
#'    ``\code{simpson}"'s rule. Alternatively and for non-equidistant grids,
#'    ``\code{trapezoidal}" or ``\code{riemann}".
#' @param basistype type of bivariate basis used. Corresponds to either \code{mgcv::s}
#'    or \code{mgcv::te}. ``\code{te}" option is only allowed when
#'    \code{domain="standardized"} and \code{interaction="nonparametric"}.
#' @param rescale.unit logical, indicating whether the \code{argvals} and {Tind}
#'    indices should be rescaled to go from 0 to 1. Rescaling occurs after
#'    \code{T.trans} is applied.
#' @param splinepars optional arguments specifying options for representing
#'    and penalizing the functional coefficient. These are passed directly to
#'    \code{mgcv::s} or \code{mgcv::te}. Defaults to the default choices
#'    of the associated \code{basistype}. See \code{\link{s}} or \code{\link{te}}
#'    in \code{mgcv} for details.
#' @details The default (``untransformed") variable-domain functional regression
#'    model uses the term \eqn{\frac1{T_i}\int_0^{T_i}X_i(t)\beta(t,T_i)dt} to
#'    incorporate the functional predictor. This term imposes
#'    a smooth (nonparametric) interaction between \eqn{t} and \eqn{T_i}. The domain
#'    of the coefficient function is the triangular (or trapezoidal) surface
#'    defined by \eqn{{t,T_i: 0\le t\le T_i}}. Tensor product smooths
#'    (\code{basistype="te"}) are not allowed, but any bivariate basis
#'    allowed by \code{mgcv::s} is supported, such as thin-plate regression splines
#'    (the default).
#'    
#'    The ``lagged" VDFR model is similar to the ``untransformed" model, but
#'    each function is aligned according to their final measurement of
#'    \eqn{X_i(t)} as opposed to their first. This model is equivalent to
#'    the untransformed model, up to the assumptions imposed by the basis
#'    choice and smoothness.
#'    
#'    The ``standardized" models apply the subject-specific domain transformation
#'    \eqn{s = t/T_i}, which linearly stretches (or compresses) each function to
#'    the domain \eqn{[0,1]}. For nonparametric interactions, the functional
#'    predictor is incorporated into the model by the term
#'    \eqn{\int_0^1 X*_i(s)\beta*(s,T_i) dt}, where \eqn{X*_i(s) = X_i(sT_i)}
#'    are the rescaled functional predictors. Because we still allow the coefficient
#'    function \eqn{\beta*(s,T_i)} to change with \eqn{T_i}, this model is
#'    equivalent to the untransformed model, up to the assumptions imposed by the
#'    numerical integration, basis choice, and smoothness. Practically, results differ
#'    primarily due to the smoothness assumptions. Because smoothness is imposed on the
#'    rescaled domain \eqn{{s,T_i: 0\le s\le 1, 0\le T_i\le max_i(T_i)}}, when transformed
#'    back to the original time domain \eqn{t}, the resulting coefficient function can be
#'    less smooth (across \eqn{t}) for smaller values of \eqn{T_i} than for larger values.
#'    
#'    Since the domain of the standardized coefficient function is rectangular, tensor product
#'    bases are allowed. This form also allows us to easily parameterize the interaction
#'    between \eqn{t} and \eqn{T_i}. The software supports three different parameterizations
#'    of this interaction: ``linear" implies \eqn{\beta*(s,T_i) = \beta_1(s) + T_i\beta_2(s)},
#'    ``quadratic" implies \eqn{\beta*(s,T_i) = \beta_1(s) + T_i\beta_2(s) + T_i^2\beta_3(s)},
#'    and ``none" implies \eqn{\beta*(s,T_i) = \beta(s)}. Note that this last parameterization
#'    implies no interaction at all, and is equivalent to \code{lf()} using the standardized
#'    functional predictors.
#' @return a list with the following entries
#'    \item{call}{a \code{call} to \code{s} or \code{te}, using the appropriately constructed
#'      weight matrices}
#'    \item{data}{data used by the \code{call}, which includes the matrices indicated
#'      by \code{tindname}, \code{Tindname}, and \code{LXname}}
#'    \item{L}{the matrix of weights used for the integration}
#'    \item{tindname}{the name used for the \code{argvals} variable in the \code{formula}
#'      used by \code{mgcv::gam}}
#'    \item{Tindname}{the name used for the \code{Tind} variable in the \code{formula}
#'      used by \code{mgcv::gam}}
#'    \item{LXname}{the name of the \code{by} variable used by \code{s} or \code{te}
#'      in the \code{formula} for \code{mgcv::gam}}
#' @export
#' @author Jonathan E. Gellar <JGellar@@mathematica-mpr.com>
#' @references Gellar, Jonathan E., Elizabeth Colantuoni, Dale M. Needham, and
#'    Ciprian M. Crainiceanu. Variable-Domain Functional Regression for Modeling
#'    ICU Data. Journal of the American Statistical Association,
#'    109(508):1425-1439, 2014.
#' @seealso \code{\link{pfr}}, \code{\link{lf}}, mgcv's
#'    \code{\link{linear.functional.terms}}.
#' @examples
#'   data(sofa)
#'   fit.vd1 <- pfr(death ~ lf.vd(SOFA) + age + los, family="binomial",
#'                  data=sofa)
#'   fit.vd2 <- pfr(death ~ lf.vd(SOFA, domain="lagged") + age + los,
#'                  family="binomial", data=sofa)
#'   fit.vd3 <- pfr(death ~ lf.vd(SOFA, domain="standardized") + age + los,
#'                  family="binomial", data=sofa)
#'   fit.vd4 <- pfr(death ~ lf.vd(SOFA, domain="standardized", basistype="te")
#'                                + age + los,
#'                  family="binomial", data=sofa)
#'   ests <- lapply(1:4, function(i) {
#'     est.i <- coef(get(paste0("fit.vd", i)))
#'     est.i[est.i$SOFA.tmat<=est.i$SOFA.Tmat,]
#'   })
#'   
#' 

lf.vd <- function(X, argvals = seq(0, 1, l = ncol(X)), Tind=NULL,
                  T.trans=identity, domain=c("untransformed", "lagged", "standardized"),
                  interaction=c("nonparametric", "none", "linear", "quadratic"),
                  integration = c("simpson", "trapezoidal", "riemann"),
                  basistype=c("s","te","t2"), rescale.unit = TRUE, splinepars = NULL) {
  
  n = nrow(X)
  J = ncol(X)
  J.i <- apply(X, 1, function(x) max(which(!is.na(x))))
  domain <- match.arg(domain)
  interaction <- match.arg(interaction)
  integration <- match.arg(integration)
  basistype <- match.arg(basistype)
  tindname <- paste(deparse(substitute(X)), ".tmat", sep = "")
  Tindname <- paste(deparse(substitute(X)), ".Tmat", sep = "")
  LXname <- paste("L.", deparse(substitute(X)), sep = "")
  splinefun <- as.symbol(basistype)
  frmls <- if (is.null(splinepars)) {
    NULL
  } else {
    frmls <- formals(getFromNamespace(deparse(splinefun), ns = "mgcv"))
    modifyList(frmls[names(frmls) %in% names(splinepars)], 
               splinepars)
  }
  
  # Check domain/interaction/basis compatability
  if (domain %in% c("untransformed","lagged")) {
    if (interaction!="nonparametric") {
      stop("Untransformed and lagged domains require nonparametric interactions.")
    } else if (basistype %in% c("te","ti","t2")) {
      stop("Tensor product smooths are not supported for non-rectangular domains.")
    } else if (!is.null(splinepars)) {
      if (!is.null(splinepars$bs)) {
        if (!(splinepars$bs %in% c("ts","tp"))) {
          stop("Basis not supported for non-rectangular domains.")
        }
      }
    }
  } else if (interaction!="nonparametric" & basistype %in% c("te","ti","t2")) {
    stop("Tensor product smooths are not allowed for parametric interactions.")
  }
  
  # Create index matrices
  if (is.null(dim(argvals))) {
    argvals <- t(argvals)
    stopifnot(ncol(argvals) == J)
    if (nrow(argvals) == 1) {
      argvals <- matrix(as.vector(argvals), nrow = n, ncol = J, 
                        byrow = T)
    }
    stopifnot(nrow(argvals) == n)
  }
  if (is.null(Tind)) 
    Tind <- sapply(1:nrow(X), function(i) argvals[i,max(which(!is.na(X[i,])))])
  Tind <- T.trans(Tind)
  if (is.null(dim(Tind))) {
    Tind <- t(Tind)
    stopifnot(ncol(Tind) == n)
    if (nrow(Tind) == 1) {
      Tind <- matrix(as.vector(Tind), nrow = n, ncol = J)
    }
    stopifnot(nrow(argvals) == n)
  }
  if (rescale.unit) {
    argvals <- (argvals-min(argvals))/(max(argvals)-min(argvals))
    Tind <- (Tind-min(Tind))/(max(Tind)-min(Tind))
  }
  
  # Process functional predictor
  if (domain=="standardized") {
    X <- t(apply(X, 1, function(x) {
      J.i <- sum(!is.na(x))
      if (J.i==1) {
        rep(x[1], J)
      } else {
        approx(x=seq(0,1,length=J.i), y=x[1:J.i],
               xout=seq(0,1,length=J))$y
      }
    }))
    L <- getL(argvals, integration=integration)
    LX <- L*X
  } else {
    L <- getL(argvals, integration=integration, n.int=J.i)
    LX <- L*X
    if (domain=="lagged") {
      LX <- t(apply(LX, 1, function(x) {
        c(rep(NA,sum(is.na(x))), x[!is.na(x)])
      }))
    }
    LX[is.na(LX)] <- 0
  }
  
  if (interaction=="nonparametric") {
    data <- list(argvals, Tind, LX)
    names(data) <- c(tindname, Tindname, LXname)
    call <- as.call(c(list(splinefun),
                      as.symbol(substitute(tindname)),
                      as.symbol(substitute(Tindname)),
                      by=as.symbol(substitute(LXname)), frmls))
  } else {
    data <- list(argvals, LX)
    names(data) <- c(tindname, LXname)
    call <- as.call(c(list(splinefun), as.symbol(substitute(tindname)),
                      by=as.symbol(substitute(LXname)), frmls))    
    if (interaction %in% c("linear","quadratic")) {
      stop("Linear and quadratic interactions are not currently supported... stay tuned")
      LX.lin <- LX * Tind
      LXname.lin <- paste0(LXname, ".lin")
      data[[LXname.lin]] <- LX.lin
      call <- call("+", call, as.call(c(list(splinefun), as.symbol(substitute(tindname)),
                                        by=as.symbol(substitute(LXname.lin)), frmls)))
    }
    if (interaction == "quadratic") {
      LX.qud <- LX * Tind^2
      LXname.qud <- paste0(LXname, ".qud")
      data[[LXname.qud]] <- LX.qud
      call <- call("+", call, as.call(c(list(splinefun), as.symbol(substitute(tindname)),
                                        by=as.symbol(substitute(LXname.qud)), frmls)))
    }
  }
  
  res <- list(call = call, data = data, L = L,
              tindname = tindname, Tindname=Tindname, LXname = LXname)
  return(res)
}


#' Get the weight matrix for a linear functional term
#' 
#' @keywords internal

getL <- function(argvals, integration, n.int=NULL) {
  nt <- ncol(argvals)
  if (is.null(n.int)) {n.int=rep(nt,nrow(argvals))}
  L <- t(sapply(1:nrow(argvals), function(i) {
    # L <- t(sapply(1:2, function(i) {
    nt.i <- n.int[i]
    if (nt.i==1) {
      c(1,rep(0,nt-1))
    } else {
      argvals.i <- argvals[i,1:nt.i]
      L.i <- switch(integration, simpson = {
        ((argvals.i[nt.i] - argvals.i[1])/nt.i)/3 * c(1, rep(c(4, 
                                                               2), length = nt.i - 2), 1)
      }, trapezoidal = {
        diffs <- diff(argvals.i)
        if (length(diffs)>1) {
          0.5 * c(diffs[1], filter(diffs, filter=c(1,1))[-(nt.i-1)],
                  diffs[(nt.i-1)])
        } else {
          rep(0.5*diffs,2)
        }
      }, riemann = {
        diffs <- diff(argvals.i)
        c(mean(diffs), diffs)
      })
      c(L.i, rep(0,nt-nt.i))/sum(L.i)
    }
  }))
  L
}
