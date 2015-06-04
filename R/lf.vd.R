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
#' @param vd vector of values of containing the variable-domain width (\eqn{T_i}
#'    above). Defaults to the \code{argvals} value corresponding to the last
#'    non-\code{NA} element of \eqn{X_i(t)}.
#' @param integration method used for numerical integration. Defaults to
#'    ``\code{simpson}"'s rule. Alternatively and for non-equidistant grids,
#'    ``\code{trapezoidal}" or ``\code{riemann}".
#' @param basistype type of bivariate basis used. Corresponds to either \code{mgcv::s}
#'    or \code{mgcv::te}. ``\code{te}" option is only allowed when
#'    \code{domain="standardized"} and \code{interaction="nonparametric"}.
#' @param transform character string indicating an optional basis transformation;
#'    see Details for options.
#' @param ... optional arguments for basis and penalization to be passed to the
#'    function indicated by \code{basistype}. These could include, for example,
#'    \code{"bs"}, \code{"k"}, \code{"m"}, etc. See \code{\link{s}} or
#'    \code{\link{te}} for details.
#'    
#' @details The variable-domain functional regression model uses the term
#'    \eqn{\frac1{T_i}\int_0^{T_i}X_i(t)\beta(t,T_i)dt} to incorporate a
#'    functional predictor with subject-specific domain width. This term imposes
#'    a smooth (nonparametric) interaction between \eqn{t} and \eqn{T_i}. The domain
#'    of the coefficient function is the triangular (or trapezoidal) surface
#'    defined by \eqn{{t,T_i: 0\le t\le T_i}}. The default basis uses
#'    bivariate thin-plate regression splines.
#'    
#'    Different basis tranformations can result in different properties; see
#'    Gellar, et al. (2014) for a more complete description. We make five basis
#'    transformations easily accessable using the \code{transform} argument.
#'    This argument is a character string that can take one of the following
#'    values:
#'    \enumerate{
#'      \item \code{"lagged"}: transforms \code{argvals} to \code{argvals - vd}
#'      \item \code{"standardized"}: transforms \code{argvals} to \code{argvals/vd},
#'        and then rescales \code{vd} linearly so it ranges from 0 to 1
#'      \item \code{"linear"}: first transforms the domain as in
#'        \code{"standardized"}, then parameterizes the interaction with
#'        \code{"vd"} to be linear
#'      \item \code{"quadratic"}: first transforms the domain as in
#'        \code{"standardized"}, then parameterizes the interaction with
#'        \code{"vd"} to be quadratic
#'      \item \code{"noInteraction"}: first transforms the domain as in
#'        \code{"standardized"}, then reduces the bivariate basis to univariate
#'        with no effect of \code{vd}
#'    }
#'    
#'    These basis transformations rely on the basis constructors available in
#'    the \code{mgcvTrans} package. For more specific control over the
#'    transformations, you can use 
#'    
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
#'      by \code{argname}, \code{Tindname}, and \code{LXname}}
#'    \item{L}{the matrix of weights used for the integration}
#'    \item{argname}{the name used for the \code{argvals} variable in the \code{formula}
#'      used by \code{mgcv::gam}}
#'    \item{Tindname}{the name used for the \code{Tind} variable in the \code{formula}
#'      used by \code{mgcv::gam}}
#'    \item{LXname}{the name of the \code{by} variable used by \code{s} or \code{te}
#'      in the \code{formula} for \code{mgcv::gam}}
#' @export
#' @author Jonathan E. Gellar <jgellar1@@jhu.edu>
#' @references Gellar, Jonathan E., Elizabeth Colantuoni, Dale M. Needham, and
#'    Ciprian M. Crainiceanu. Variable-Domain Functional Regression for Modeling
#'    ICU Data. Journal of the American Statistical Association,
#'    109(508):1425-1439, 2014.
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
#'     coef(get(paste0("fit.vd", i)))
#'   })
#' @seealso \code{\link{pfr}}, \code{\link{lf}}, mgcv's
#'    \code{\link{linear.functional.terms}}.

lf.vd <- function(X, argvals = seq(0, 1, l = ncol(X)), vd=NULL,
                  integration = c("simpson", "trapezoidal", "riemann"),
                  basistype=c("s","te","t2"),
                  transform=NULL, ...
) {
  
  integration <- match.arg(integration)
  basistype <- match.arg(basistype)
  
  # Set up functions
  n = nrow(X)
  J = ncol(X)
  J.i <- apply(X, 1, function(x) max(which(!is.na(x))))
  
  # Create coordinate matrices
  if (is.null(dim(argvals))) {
    argvals <- t(argvals)
    stopifnot(ncol(argvals) == J)
    if (nrow(argvals) == 1) {
      argvals <- matrix(as.vector(argvals), nrow = n, ncol = J, byrow = T)
    }
    stopifnot(nrow(argvals) == n)
  }
  if (is.null(vd)) 
    # Defaults to the last non-NA value of argvals for that function
    vd <- sapply(1:nrow(X), function(i) argvals[i,J.i[i]])
  if (is.null(dim(vd))) {
    vd <- t(vd)
    stopifnot(ncol(vd) == n)
    if (nrow(vd) == 1) {
      vd <- matrix(as.vector(vd), nrow = n, ncol = J)
    }
    stopifnot(nrow(vd) == n)
  }
  
  # Process Functional Predictor
  L <- getL(argvals, integration=integration, n.int=J.i)
  LX <- L*X
  
  # Zero-out unused coordinates. For argvals and vd, this means we set them to
  # any other (used) coordinate combination. This ensures they won't affect the
  # range of values used to set up the basis.
  # For LX, we set the weight to 0.
  argvals[is.na(LX)] <- argvals[!is.na(LX)][1]
  vd[is.na(LX)]      <- vd[!is.na(LX)][1]
  LX[is.na(LX)]      <- 0
  
  # Term names for basis construction
  argname <- paste(deparse(substitute(X)), ".arg", sep = "")
  vdname  <- paste(deparse(substitute(X)), ".vd",  sep = "")
  LXname <- paste("L.", deparse(substitute(X)), sep = "")
  splinefun <- as.symbol(basistype)
  
  # Set up transformations
  # TRANSFORMS
  dots <- list(...)
  bs0 <- dots$bs
  xt0 <- dots$xt
  if (!is.null(transform)) {
    if (transform=="lagged") {
      dots$bs <- "dt"
      dots$xt <- list(tf=list("s-t"), bs=bs0, xt=xt0)
    } else if (transform=="standardized") {
      dots$bs <- "dt"
      dots$xt <- list(tf=list("s/t", "linear01"), bs=bs0, xt=xt0)
    } else if (transform=="noInteraction") {
      dots$bs <- "dt"
      dots$xt <- list(tf="s/t", bs="pi", xt=list(g="none", bs=bs0, xt=xt0))
    } else if (transform=="linear") {
      dots$bs <- "dt"
      dots$xt <- list(tf="s/t", bs="pi", xt=list(g="linear", bs=bs0, xt=xt0))
    } else if (transform=="quadratic") {
      dots$bs <- "dt"
      dots$xt <- list(tf="s/t", bs="pi", xt=list(g="quadratic", bs=bs0, xt=xt0))
    }
  }
  
  # Set up basis
  data <- list(argvals, vd, LX)
  names(data) <- c(argname, vdname, LXname)
  call <- as.call(c(list(splinefun),
                    as.symbol(substitute(argname)),
                    as.symbol(substitute(vdname)),
                    by=as.symbol(substitute(LXname)), 
                    dots
  ))
  
  res <- list(call = call, data = data, L = L,
              argname = argname, vdname=vdname, LXname = LXname)
  return(res)
}


#' Get the weight matrix for a linear functional term
#' 
#' @keywords internal

getL <- function(tind, integration, n.int=NULL) {
  nt <- ncol(tind)
  if (is.null(n.int)) {n.int=rep(nt,nrow(tind))}
  L <- t(sapply(1:nrow(tind), function(i) {
    # L <- t(sapply(1:2, function(i) {
    nt.i <- n.int[i]
    if (nt.i==1) {
      c(1,rep(0,nt-1))
    } else {
      tind.i <- tind[i,1:nt.i]
      L.i <- switch(integration, simpson = {
        ((tind.i[nt.i] - tind.i[1])/nt.i)/3 * c(1, rep(c(4, 
                                                         2), length = nt.i - 2), 1)
      }, trapezoidal = {
        diffs <- diff(tind.i)
        if (length(diffs)>1) {
          0.5 * c(diffs[1], filter(diffs, filter=c(1,1))[-(nt.i-1)],
                  diffs[(nt.i-1)])
        } else {
          rep(0.5*diffs,2)
        }
      }, riemann = {
        diffs <- diff(tind.i)
        c(mean(diffs), diffs)
      })
      c(L.i, rep(0,nt-nt.i))/sum(L.i)
    }
  }))
  L
}
