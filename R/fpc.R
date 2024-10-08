#' Construct a FPC regression term
#'
#' Constructs a functional principal component regression (Reiss and Ogden,
#' 2007, 2010) term for inclusion in an \code{mgcv::gam}-formula (or
#' \code{{bam}} or \code{{gamm}} or \code{gamm4:::gamm}) as
#' constructed by \code{{pfr}}. Currently only one-dimensional functions
#' are allowed.
#'
#' @param X functional predictors, typically expressed as an \code{N} by \code{J} matrix,
#'   where \code{N} is the number of columns and \code{J} is the number of
#'   evaluation points. May include missing/sparse functions, which are
#'   indicated by \code{NA} values. Alternatively, can be an object of class
#'   \code{"fd"}; see \code{[fda]{fd}}.
#' @param argvals indices of evaluation of \code{X}, i.e. \eqn{(t_{i1},.,t_{iJ})} for
#'   subject \eqn{i}. May be entered as either a length-\code{J} vector, or as
#'   an \code{N} by \code{J} matrix. Indices may be unequally spaced. Entering
#'   as a matrix allows for different observations times for each subject. If
#'   \code{NULL}, defaults to an equally-spaced grid between 0 or 1 (or within
#'   \code{X$basis$rangeval} if \code{X} is a \code{fd} object.)
#' @param method the method used for finding principal components. The default
#'   is an unconstrained SVD of the \eqn{XB} matrix. Alternatives include
#'   constrained (functional) principal components approaches
#' @param ncomp number of principal components. if \code{NULL}, chosen by \code{pve}
#' @param pve proportion of variance explained; used to choose the number of
#'   principal components
#' @param penalize if \code{TRUE}, a roughness penalty is applied to the
#'   functional estimate. Defaults to \code{TRUE} if \code{method=="svd"}
#'   (corresponding to the FPCR_R method of Reiss and Ogden (2007)), and
#'   \code{FALSE} if \code{method!="svd"} (corresponding to FPCR_C).
#' @param bs two letter character string indicating the \code{mgcv}-style basis
#'   to use for pre-smoothing \code{X}
#' @param k the dimension of the pre-smoothing basis
#' @param ... additional options to be passed to \code{{lf}}. These include
#'   \code{argvals}, \code{integration}, and any additional options for the
#'   pre-smoothing basis (as constructed by \code{mgcv::s}), such as \code{m}.
#'
#' @details
#' \code{fpc} is a wrapper for \code{{lf}}, which defines linear
#' functional predictors for any type of basis for inclusion in a \code{pfr}
#' formula. \code{fpc} simply calls \code{lf} with the appropriate options for
#' the \code{fpc} basis and penalty construction.
#'
#' This function implements both the FPCR-R and FPCR-C methods of Reiss and
#' Ogden (2007). Both methods consist of the following steps:
#' \enumerate{
#'  \item project \eqn{X} onto a spline basis \eqn{B}
#'  \item perform a principal components decomposition of \eqn{XB}
#'  \item use those PC's as the basis in fitting a (generalized) functional
#'    linear model
#' }
#'
#' This implementation provides options for each of these steps. The basis
#' for in step 1 can be specified using the arguments \code{bs} and \code{k},
#' as well as other options via \code{...}; see \code{[mgcv]{s}} for
#' these options. The type of PC-decomposition is specified with \code{method}.
#' And the FLM can be fit either penalized or unpenalized via \code{penalize}.
#'
#' The default is FPCR-R, which uses a b-spline basis, an unconstrained
#' principal components decomposition using \code{{svd}}, and the FLM
#' fit with a second-order difference penalty. FPCR-C can be selected by
#' using a different option for \code{method}, indicating a constrained
#' ("functional") PC decomposition, and by default an unpenalized fit of the
#' FLM.
#'
#' FPCR-R is also implemented in \code{{fpcr}}; here we implement the
#'   method for inclusion in a \code{pfr} formula.
#'
#' @section NOTE:
#' Unlike \code{{fpcr}}, \code{fpc} within a \code{pfr} formula does
#' not automatically decorrelate the functional predictors from additional
#' scalar covariates.
#'
#' @return The result of a call to \code{{lf}}.
#'
#' @references
#' Reiss, P. T. (2006). Regression with signals and images as predictors. Ph.D.
#' dissertation, Department of Biostatistics, Columbia University. Available
#' at http://works.bepress.com/phil_reiss/11/.
#'
#' Reiss, P. T., and Ogden, R. T. (2007). Functional principal component
#' regression and functional partial least squares. \emph{Journal of the
#' American Statistical Association}, 102, 984-996.
#'
#' Reiss, P. T., and Ogden, R. T. (2010). Functional generalized linear models
#' with images as predictors. \emph{Biometrics}, 66, 61-69.
#'
#' @author Jonathan Gellar \email{JGellar@@mathematica-mpr.com}, Phil Reiss
#'   \email{phil.reiss@@nyumc.org}, Lan Huo \email{lan.huo@@nyumc.org}, and
#'   Lei Huang \email{huangracer@@gmail.com}
#'
#'
#' @examples
#' data(gasoline)
#' par(mfrow=c(3,1))
#'
#' # Fit PFCR_R
#' gasmod1 <- pfr(octane ~ fpc(NIR, ncomp=30), data=gasoline)
#' plot(gasmod1, rug=FALSE)
#' est1 <- coef(gasmod1)
#'
#' # Fit FPCR_C with fpca.sc
#' gasmod2 <- pfr(octane ~ fpc(NIR, method="fpca.sc", ncomp=6), data=gasoline)
#' plot(gasmod2, se=FALSE)
#' est2 <- coef(gasmod2)
#'
#' # Fit penalized model with fpca.face
#' gasmod3 <- pfr(octane ~ fpc(NIR, method="fpca.face", penalize=TRUE), data=gasoline)
#' plot(gasmod3, rug=FALSE)
#' est3 <- coef(gasmod3)
#'
#' par(mfrow=c(1,1))
#' ylm <- range(est1$value)*1.35
#' plot(value ~ X.argvals, type="l", data=est1, ylim=ylm)
#' lines(value ~ X.argvals, col=2, data=est2)
#' lines(value ~ X.argvals, col=3, data=est3)
#'
#' @seealso \code{{lf}}, \code{{smooth.construct.fpc.smooth.spec}}
#' @importFrom methods is
#' @export

fpc <- function(X, argvals=NULL,
                method=c("svd", "fpca.sc", "fpca.face", "fpca.ssvd"),
                ncomp=NULL, pve=0.99, penalize=(method=="svd"),
                bs="ps", k=40, ...) {
  method <- match.arg(method)

  if (is(X, "fd")) {
    # If X is an fd object, turn it back into a (possibly pre-smoothed) matrix
    if (is.null(argvals))
      argvals <- seq(X$basis$rangeval[1], X$basis$rangeval[2],
                                length = length(X$fdnames[[1]]))
    X <- t(eval.fd(argvals, X))
  } else if (is.null(argvals))
    argvals <- seq(0, 1, l = ncol(X))

  if("xt" %in% names(list(...)))
    stop("fpc() does not accept an xt-argument.")
  xt <- call("list", X=substitute(X), method=method, npc=ncomp, pve=pve,
             penalize=penalize, bs=bs)
  lf (X=X, argvals=argvals, bs="fpc", k=k, xt=xt, ...)
}



#' Basis constructor for FPC terms
#'
#' @param object a \code{fpc.smooth.spec} object, usually generated by a
#'   term \code{s(x, bs="fpc")}; see Details.
#' @param data a list containing the data (including any \code{by} variable)
#'   required by this term, with names corresponding to \code{object$term}
#'   (and \code{object$by}). Only the first element of this list is used.
#' @param knots not used, but required by the generic \code{smooth.construct}.
#'
#' @details
#' \code{object} must contain an \code{xt} element. This is a list that can
#'   contain the following elements:
#' \describe{
#'   \item{X}{(required) matrix of functional predictors}
#'   \item{method}{(required) the method of finding principal components;
#'     options include \code{"svd"} (unconstrained), \code{"fpca.sc"},
#'     \code{"fpca.face"}, or \code{"fpca.ssvd"}}
#'   \item{npc}{(optional) the number of PC's to retain}
#'   \item{pve}{(only needed if \code{npc} not supplied) the percent variance
#'     explained used to determine \code{npc}}
#'   \item{penalize}{(required) if \code{FALSE}, the smoothing parameter is
#'     set to 0}
#'   \item{bs}{the basis class used to pre-smooth \code{X}; default is \code{"ps"}}
#' }
#'
#' Any additional options for the pre-smoothing basis (e.g. \code{k}, \code{m},
#'   etc.) can be supplied in the corresponding elements of \code{object}.
#'   See \code{[mgcv]{s}} for a full list of options.
#'
#' @return An object of class \code{"fpc.smooth"}. In addtional to the elements
#'   listed in \code{{smooth.construct}}, the object will contain
#'   \item{sm}{the smooth that is fit in order to generate the basis matrix
#'     over \code{object$term}}
#'   \item{V.A}{the matrix of principal components}
#'
#' @references
#' Reiss, P. T., and Ogden, R. T. (2007). Functional principal component
#' regression and functional partial least squares. \emph{Journal of the
#' American Statistical Association}, 102, 984-996.
#'
#' @author Jonathan Gellar \email{JGellar@@mathematica-mpr.com}
#' @seealso \code{{fpcr}}
#' @export

smooth.construct.fpc.smooth.spec <- function(object, data, knots) {
  xt <- object$xt
  if (is.null(xt$npc) & is.null(xt$pve))
    stop("Need either the number of PC's or the PVE used to determine this number")
  if (is.null(xt$X)) stop("X matrix must be supplied as part of xt")
  if (is.null(xt$bs)) xt$bs <- "ps"

  # Create mini-basis
  obj <- object
  obj$by <- "NA"
  class(obj) <- sub("fpc", xt$bs, class(obj))
  obj$xt <- obj$xt[!(names(obj$xt) %in% c("bs", "X", "npc", "pve"))]
  obj$sp <- 0
  sm <- mgcv::smooth.construct(obj, data=data[obj$term], knots=NULL)

  # Create PC Basis: need X, B, and V.A
  XB <- xt$X %*% sm$X
  res <- switch(xt$method,
                svd = {
                  XB.svd <- svd(XB)
                  pve <- XB.svd$d/sum(XB.svd$d)
                  npc <- ifelse(is.null(xt$npc),
                                min(which(cumsum(pve) > xt$pve)),
                                #min(which(cumsum(XB.svd$d) > xt$pve * sum(XB.svd$d))),
                                xt$npc)
                  list(V.A = XB.svd$v[,1:npc], pve=pve[1:npc])
                }, fpca.sc = {
                  fp <- fpca.sc(XB, npc=xt$npc, pve=xt$pve)
                  list(V.A=fp$efunctions, pve=cumsum(fp$evalues)/sum(fp$evalues))
                }, fpca.face = {
                  if (!is.null(xt$npc)) xt$pve <- 1
                  fp <- fpca.face(XB, npc=xt$npc, pve=xt$pve)
                  list(V.A=fp$efunctions, pve=cumsum(fp$evalues)/sum(fp$evalues))
                }, fpca.ssvd = {
                  fp <- fpca.ssvd(XB, npc=ifelse(is.null(xt$npc), NA, xt$npc))
                  list(V.A=fp$efunctions, pve=cumsum(fp$evalues)/sum(fp$evalues))
                })
  V.A <- res$V.A
  npc <- ncol(V.A)

  ## Limit to the first npc PC's
  #npc <- ifelse(is.null(xt$npc),
  #            min(which(cumsum(XB.svd$d) > xt$pve * sum(XB.svd$d))), xt$npc)
  #V.A <- XB.svd$v[,1:npc]

  # Return Object
  object$X <- sm$X %*% V.A
  object$S <- lapply(sm$S, function(smat) crossprod(V.A, smat) %*% V.A)
  if (!xt$penalize) object$sp <- 0

  # Other stuff... need to check these
  object$null.space.dim <- npc - qr(object$S[[1]])$rank
  object$rank <- npc - object$null.space.dim
  object$df   <- npc #need this for gamm
  object$sm   <- sm
  object$V.A  <- V.A
  object$pve  <- res$pve
  class(object) <- "fpc.smooth"
  object
}


#' mgcv-style constructor for prediction of FPC terms
#'
#' @param object a \code{fpc.smooth} object created by
#'   \code{{smooth.construct.fpc.smooth.spec}}, see
#'   \code{[mgcv]{smooth.construct}}
#' @param data  see \code{[mgcv]{smooth.construct}}
#' @return design matrix for FPC terms
#' @author Jonathan Gellar
#' @export

Predict.matrix.fpc.smooth <- function(object, data) {
  mgcv::Predict.matrix(object$sm, data) %*% object$V.A
}



