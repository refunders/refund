#' Construct a VDFR regression term
#' 
#' This function defines the a variable-domain functional regression term
#' for inclusion in an \code{\link[mgcv]{gam}}-formula (or \code{\link[mgcv]{bam}} or
#' \code{\link[mgcv]{gamm}} or \code{gamm4::gamm} as constructed by
#' \code{\link{pfr}}. These are functional predictors for which each function is
#' observed over a domain of different width.
#' The default is the term \eqn{1/T_i\int_0^{T_i}X_i(t)\beta(t,T_i)dt},
#' where \eqn{X_i(t)} is a functional predictor of length \eqn{T_i} and \eqn{\beta(t,T_i)}
#' is an unknown bivariate coefficient function. Various domain transformations
#' are available, such as lagging or domain-standardizing the coordinates, or
#' parameterizing the interactions; these often result in improved model fit.
#' Basis choice is fully custiomizable using the options of
#' \code{\link[mgcv]{s}} and \code{\link[mgcv]{te}}.
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
#' @param basistype character string indicating type of bivariate basis used.
#'    Options include \code{"s"} (the default), \code{"te"}, and \code{"t2"},
#'    which correspond to \code{mgcv::s}, \code{mgcv::te}, and \code{mgcv::t2}.
#' @param transform character string indicating an optional basis transformation;
#'    see Details for options.
#' @param mp for \code{transform=="linear"} or \code{transform=="quadratic"},
#'    \code{TRUE} to use multiple penalties for the smooth (one for each marginal
#'    basis). If \code{FALSE}, penalties are concatonated into a single
#'    block-diagonal penalty matrix (with one smoothing parameter).
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
#'        with no effect of \code{vd}. This would be equivalent to using
#'        \code{\link{lf}} on the domain-standardized predictor functions.
#'    }
#'    
#'    The practical effect of using the \code{"lagged"} basis is to increase
#'    smoothness along the right (diagonal) edge of the resultant estimate.
#'    The practical effect of using a \code{"standardized"} basis is to allow
#'    for greater smoothness at high values of \eqn{T_i} compared to lower
#'    values.
#'    
#'    These basis transformations rely on the basis constructors
#'    available in the \code{mgcvTrans} package. For more specific control over
#'    the transformations, you can use \code{bs="dt"} and/or \code{bs="pi"};
#'    see \code{\link{smooth.construct.dt.smooth.spec}} or
#'    \code{\link{smooth.construct.pi.smooth.spec}} for an explanation of the
#'    options (entered through the \code{xt} argument of \code{lf.vd}/\code{s}).
#'    
#'    Note that tensor product bases are only recommended when a standardized
#'    transformation is used. Without this transformation, just under half of
#'    the "knots" used to define the basis will fall outside the range of the
#'    data and have no data available to estimate them. The penalty allows
#'    the corresponding coefficients to be estiamted, but results may be
#'    unstable.
#'    
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
#' @author Jonathan E. Gellar <JGellar@@mathematica-mpr.com>
#' @references Gellar, Jonathan E., Elizabeth Colantuoni, Dale M. Needham, and
#'    Ciprian M. Crainiceanu. Variable-Domain Functional Regression for Modeling
#'    ICU Data. Journal of the American Statistical Association,
#'    109(508):1425-1439, 2014.
#' @examples
#'   data(sofa)
#'   fit.vd1 <- pfr(death ~ lf.vd(SOFA) + age + los, family="binomial",
#'                  data=sofa)
#'   fit.vd2 <- pfr(death ~ lf.vd(SOFA, transform="lagged") + age + los,
#'                  family="binomial", data=sofa)
#'   fit.vd3 <- pfr(death ~ lf.vd(SOFA, transform="standardized") + age + los,
#'                  family="binomial", data=sofa)
#'   fit.vd4 <- pfr(death ~ lf.vd(SOFA, transform="standardized",
#'                                basistype="te") + age + los,
#'                  family="binomial", data=sofa)
#'   fit.vd5 <- pfr(death ~ lf.vd(SOFA, transform="linear", bs="ps") + age + los,
#'                  family="binomial", data=sofa)
#'   fit.vd6 <- pfr(death ~ lf.vd(SOFA, transform="quadratic", bs="ps") + age + los,
#'                  family="binomial", data=sofa)
#'   fit.vd7 <- pfr(death ~ lf.vd(SOFA, transform="noInteraction", bs="ps") + age + los,
#'                  family="binomial", data=sofa)
#'   
#'   ests <- lapply(1:7, function(i) {
#'     c.i <- coef(get(paste0("fit.vd", i)), n=173, n2=173) 
#'     c.i[(c.i$SOFA.arg <= c.i$SOFA.vd),]
#'   })
#'   
#'   # Try plotting for each i
#'   i <- 1
#'   lims <- c(-2,8)
#'   if (requireNamespace("ggplot2", quietly = TRUE) &
#'       requireNamespace("RColorBrewer", quietly = TRUE)) {
#'         est <- ests[[i]]
#'         est$value[est$value<lims[1]] <- lims[1]
#'         est$value[est$value>lims[2]] <- lims[2]
#'         ggplot2::ggplot(est, ggplot2::aes(SOFA.arg, SOFA.vd)) +
#'           ggplot2::geom_tile(ggplot2::aes(colour=value, fill=value)) +
#'           ggplot2::scale_fill_gradientn(  name="", limits=lims,
#'                     colours=rev(RColorBrewer::brewer.pal(11,"Spectral"))) +
#'           ggplot2::scale_colour_gradientn(name="", limits=lims,
#'                     colours=rev(RColorBrewer::brewer.pal(11,"Spectral"))) +
#'           ggplot2::scale_y_continuous(expand = c(0,0)) +
#'           ggplot2::scale_x_continuous(expand = c(0,0)) +
#'           ggplot2::theme_bw()
#'   }
#'   
#' @seealso \code{\link{pfr}}, \code{\link{lf}}, mgcv's
#'    \code{\link{linear.functional.terms}}.

lf.vd <- function(X, argvals = seq(0, 1, l = ncol(X)), vd=NULL,
                  integration = c("simpson", "trapezoidal", "riemann"),
                  basistype=c("s","te","t2"),
                  transform=NULL, mp=TRUE, ...
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
  
  # Set up transformations
  dots <- list(...)
  bs0 <- dots$bs
  xt0 <- dots$xt
  if (!is.null(transform)) {
    # Set up dt basis call
    dots$bs <- "dt"
    if (transform=="lagged") {
      dots$xt <- list(tf=list("s-t"))
      if (!is.null(bs0)) dots$xt$bs=bs0
      if (!is.null(xt0)) dots$xt$xt=xt0
    } else if (transform=="standardized") {
      dots$xt <- list(tf=list("s/t", "linear01"))
      if (!is.null(bs0)) dots$xt$bs=bs0
      if (!is.null(xt0)) dots$xt$xt=xt0
    } else if (transform=="noInteraction") {
      dots$xt <- list(tf="s/t", bs="pi", xt=list(g="none"))
      if (!is.null(bs0)) dots$xt$xt$bs=bs0
      if (!is.null(xt0)) dots$xt$xt$xt=xt0
    } else if (transform=="linear") {
      dots$xt <- list(tf="s/t", bs="pi", xt=list(g="linear", mp=mp))
      if (!is.null(bs0)) dots$xt$xt$bs=bs0
      if (!is.null(xt0)) dots$xt$xt$xt=xt0
    } else if (transform=="quadratic") {
      dots$xt <- list(tf="s/t", bs="pi", xt=list(g="quadratic", mp=mp))
      if (!is.null(bs0)) dots$xt$xt$bs=bs0
      if (!is.null(xt0)) dots$xt$xt$xt=xt0
    }
    if (basistype!="s") {
      # dt basis call must go through s to allow bivariate transformations
      # (te would split the coordinates up)
      dots$xt$basistype <- basistype
      basistype <- "s"
    }
  }
  
  # Set up basis
  data <- list(argvals, vd, LX)
  names(data) <- c(argname, vdname, LXname)
  #splinefun <- as.symbol(basistype)
  call <- as.call(c(list(as.symbol(basistype)),
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

#' SOFA (Sequential Organ Failure Assessment) Data
#' 
#' A dataset containing the SOFA scores (Vincent et al, 1996). for 520 patients,
#' hospitalized in the intensive care unit (ICU) with Acute Lung Inury. Daily
#' measurements are available for as long as each one remains in the ICU. This is an
#' example of variable-domain functional data, as described by Gellar et al. (2014).
#' 
#' The data was collected as part of the Improving Care of ALI Patients (ICAP)
#' study (Needham et al., 2006). If you use this dataset as an example in
#' written work, please cite the study protocol.
#' 
#' @format A data frame with 520 rows (subjects) and 7 variables:
#' \describe{
#'   \item{death}{binary indicator that the subject died in the ICU}
#'   \item{SOFA}{520 x 173 matrix in variable-domain format (a ragged array). 
#'     Each column represents an ICU day. Each row contains the SOFA scores for
#'     a subject, one per day, for as long as the subject remained in the ICU.
#'     The remaining cells of each row are padded with \code{NA}s. SOFA scores
#'     range from 0 to 24, increasing with severity of organ failure. Missing
#'     values during one's ICU stay have been imputed using LOCF.}
#'   \item{SOFA_raw}{Identical to the \code{SOFA} element, except that it contains
#'     some missing values during one's hospitalization. These missing values
#'     arise when a subject leaves the ICU temporarily, only to be re-admitted.
#'     SOFA scores are not monitored outside the ICU.}
#'   \item{los}{ICU length of stay, i.e., the number of days the patient remained
#'     in the ICU prior to death or final discharge.}
#'   \item{age}{Patient age}
#'   \item{male}{Binay indicator for male gender}
#'   \item{Charlson}{Charlson co-morbidity index, a measure of baseline health
#'     status (before hospitalization and ALI).}
#' }
#' 
#' @references 
#'    Vincent, JL, Moreno, R, Takala, J, Willatts, S, De Mendonca, A,
#'    Bruining, H, Reinhart, CK, Suter, PM, Thijs, LG (1996). The SOFA ( Sepsis
#'    related Organ Failure Assessment) score to describe organ
#'    dysfunction/failure. Intensive Care Medicine, 22(7): 707-710.
#'    
#'    Needham, D. M., Dennison, C. R., Dowdy, D. W., Mendez-Tellez, P. A.,
#'    Ciesla, N., Desai, S. V., Sevransky, J., Shanholtz, C., Scharfstein, D.,
#'    Herridge, M. S., and Pronovost, P. J. (2006). Study protocol: The
#'    Improving Care of Acute Lung Injury Patients (ICAP) study. Critical Care
#'    (London, England), 10(1), R9.
#'    
#'    Gellar, Jonathan E., Elizabeth Colantuoni, Dale M. Needham, and
#'    Ciprian M. Crainiceanu. Variable-Domain Functional Regression for Modeling
#'    ICU Data. Journal of the American Statistical Association,
#'    109(508):1425-1439, 2014.
#'    
"sofa"
