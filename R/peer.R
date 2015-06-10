##' Construct a PEER regression term in a \code{pfr} formula
##' 
##' Defines a term \eqn{\int_{T}\beta(t)X_i(t)dt} for inclusion in a
##' \code{\link{pfr}} formula, where \eqn{\beta(t)} is estimated with
##' structured penalties (Randloph et al., 2012).
##' 
##' @param X the functional predictors, expressed as an \code{N} by \code{J}
##'   matrix
#' @param argvals vector (or matrix) of indices of evaluations of \eqn{X_i(t)}; i.e. a matrix with
#' \emph{i}th row \eqn{(t_{i1},.,t_{iJ})}.
##' @param pentype the type of penalty to apply; see Details.
##' @param Q matrix \eqn{Q} used for \code{pentype="DECOMP"}; see Details.
##' @param phia scalar \eqn{a} used for \code{pentype="DECOMP"}; see Details.
##' @param L user-supplied penalty matrix for \code{pentype="USER"}; see
##'   Details.
##' @param ... additional arguments to be passed to \code{lf} (and then
##'   possibly \code{s}). Arguments processed by \code{lf} include, for example,
##'   \code{argvals} and \code{integration} for specifying the functional
##'   argument and numerical integration. Arguments processed by \code{s}
##'   include information related to basis and penalization, such as \code{m}
##'   for specifying the order of the difference penalty; See Details.
##' 
##' @details
##' \code{peer} is a wrapper for \code{\link{lf}}, which defines linear
##' functional predictors for any type of basis. It simply calls \code{lf}
##' with the appropriate options for basis and penalty construction.
##' The type of penalty is determined by the \code{pentype} argument. There
##' are four types of penalties available:
##' \enumerate{
##'   \item \code{pentype=="RIDGE"} for a ridge penalty, the default
##'   \item \code{pentype=="D"} for a difference penalty. The order of the
##'     difference penalty may be specified by supplying an \code{m} argument
##'     (default is 2).
##'   \item \code{pentype=="DECOMP"} for a decomposition-based penalty,
##'     \eqn{bP_Q + a(I-P_Q)}, where \eqn{P_Q = Q^t(QQ^t)^{-1}Q}. The \eqn{Q}
##'     matrix must be specified by \code{Q}, and the scalar \eqn{a} by
##'     \code{phia}. The number of columns of \code{Q} must be equal to the
##'     length of the data. Each row represents a basis function where the
##'     functional predictor is expected to lie, according to prior belief.
##'   \item \code{pentype=="USER"} for a user-specified penalty matrix,
##'     supplied by the \code{L} argument.
##' }
##' 
##' The original stand-alone implementation by Madan Gopal Kundu is available in 
##' \code{\link{peer_old}}. 
##' 
##' 
##' @author Jonathan Gellar \email{JGellar@@mathematica-mpr.com} and
##'         Madan Gopal Kundu \email{mgkundu@@iupui.edu}
##' 
##' @references
##' Randolph, T. W., Harezlak, J, and Feng, Z. (2012). Structured penalties for
##' functional linear models - partially empirical eigenvectors for regression.
##' \emph{Electronic Journal of Statistics}, 6, 323-353.
##' 
##' Kundu, M. G., Harezlak, J., and Randolph, T. W. (2012). Longitudinal
##' functional models with structured penalties (arXiv:1211.4763 [stat.AP]).
##' 
##' @seealso \code{\link{pfr}}, \code{\link{smooth.construct.peer.smooth.spec}}
##'
##' @examples
##'
##' \dontrun{
##' #------------------------------------------------------------------------
##' # Example 1: Estimation with D2 penalty
##' #------------------------------------------------------------------------
##'
##' data(DTI)
##' DTI = DTI[which(DTI$case == 1),]
##' fit.D2 = pfr(pasat ~ peer(cca, pentype="D"), data=DTI)
##' plot(fit.D2)
##'
##' #------------------------------------------------------------------------
##' # Example 2: Estimation with structured penalty (need structural
##' #            information about regression function or predictor function)
##' #------------------------------------------------------------------------
##'
##' data(PEER.Sim)
##' data(Q)
##' PEER.Sim1<- subset(PEER.Sim, t==0)
##' 
##' # Setting k to max possible value
##' fit.decomp <- pfr(Y ~ peer(W, pentype="Decomp", Q=Q, k=99), data=PEER.Sim1)
##' plot(fit.decomp)
##' }
##'
##'

peer <- function(X, argvals=seq(0, 1, l=ncol(X)), pentype="RIDGE",
                 Q=NULL, phia=10^3, L=NULL,  ...) {
  
  # Catch if peer_old syntax is used
  dots <- list(...)
  dots.unmatched <- names(dots)[!(names(dots) %in% c(names(formals(lf)),
                                                     names(formals(s))))]
  if (any(dots.unmatched %in% names(formals(peer_old)))) {
    warning(paste0("The interface for peer() has changed, see ?peer and ?pfr ",
                   "for details. This interface will not be supported in the ",
                   "next refund release."))
    # Call peer_old()
    call <- sys.call()
    call[[1]] <- as.symbol("peer_old")
    ret <- eval(call, envir=parent.frame())
    return(ret)
  }
  
  # Extract xt if in ...
  xt <- if ("xt" %in% names(dots)) dots$xt else list()
  dots$xt <- NULL
  
  # Update xt
  xt$pentype <- pentype
  xt$W <- X
  # FIXME: this evaluates X, leading to ugly model output -- would want to 
  # hand over the unevaluated symbol to lf instead, but this fails as the call stack
  # is so complicated in pfr.R#249: eval(newcall)
  # substitute(X) does not work as inputs get renamed before handing over to fitter (?) 
  xt$Q <- Q
  xt$phia <- phia
  xt$L <- L
  
  lf(X=X, argvals=argvals, bs="peer", xt=xt, ...)
}

