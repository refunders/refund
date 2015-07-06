#' Construct a PEER regression term in a \code{pfr} formula
#'
#' Defines a term \eqn{\int_{T}\beta(t)X_i(t)dt} for inclusion in a
#' \code{\link{pfr}} formula, where \eqn{\beta(t)} is estimated with
#' structured penalties (Randloph et al., 2012).
#'
#' @param X functional predictors, expressed as an \code{N} by \code{J} matrix,
#'   where \code{N} is the number of columns and \code{J} is the number of
#'   evaluation points. May include missing/sparse functions, which are
#'   indicated by \code{NA} values.
#' @param argvals indices of evaluation of \code{X}, i.e. \eqn{(t_{i1},.,t_{iJ})} for
#'   subject \eqn{i}. May be entered as either a length-\code{J} vector, or as
#'   an \code{N} by \code{J} matrix. Indices may be unequally spaced. Entering
#'   as a matrix allows for different observations times for each subject.
#' @param pentype the type of penalty to apply, one of \code{"RIDGE"}, \code{"D"},
#'  \code{"DECOMP"}, or \code{"USER"}; see Details.
#' @param Q matrix \eqn{Q} used for \code{pentype="DECOMP"}; see Details.
#' @param phia scalar \eqn{a} used for \code{pentype="DECOMP"}; see Details.
#' @param L user-supplied penalty matrix for \code{pentype="USER"}; see
#'   Details.
#' @param ... additional arguments to be passed to \code{lf} (and then
#'   possibly \code{s}). Arguments processed by \code{lf} include, for example,
#'   \code{integration} for specifying the method of numerical integration.
#'   Arguments processed by \code{s}
#'   include information related to basis and penalization, such as \code{m}
#'   for specifying the order of the difference penalty; See Details.
#'   \code{xt}-argument is not allowed for \code{peer}-terms and will cause
#'   an error.
#'
#' @details
#' \code{peer} is a wrapper for \code{\link{lf}}, which defines linear
#' functional predictors for any type of basis. It simply calls \code{lf}
#' with the appropriate options for the \code{peer} basis and penalty construction.
#' The type of penalty is determined by the \code{pentype} argument. There
#' are four types of penalties available:
#' \enumerate{
#'   \item \code{pentype=="RIDGE"} for a ridge penalty, the default
#'   \item \code{pentype=="D"} for a difference penalty. The order of the
#'     difference penalty may be specified by supplying an \code{m} argument
#'     (default is 2).
#'   \item \code{pentype=="DECOMP"} for a decomposition-based penalty,
#'     \eqn{bP_Q + a(I-P_Q)}, where \eqn{P_Q = Q^t(QQ^t)^{-1}Q}. The \eqn{Q}
#'     matrix must be specified by \code{Q}, and the scalar \eqn{a} by
#'     \code{phia}. The number of columns of \code{Q} must be equal to the
#'     length of the data. Each row represents a basis function where the
#'     functional predictor is expected to lie, according to prior belief.
#'   \item \code{pentype=="USER"} for a user-specified penalty matrix,
#'     supplied by the \code{L} argument.
#' }
#'
#' The original stand-alone implementation by Madan Gopal Kundu is available in
#' \code{\link{peer_old}}.
#'
#'
#' @author Jonathan Gellar \email{JGellar@@mathematica-mpr.com} and
#'         Madan Gopal Kundu \email{mgkundu@@iupui.edu}
#'
#' @references
#' Randolph, T. W., Harezlak, J, and Feng, Z. (2012). Structured penalties for
#' functional linear models - partially empirical eigenvectors for regression.
#' \emph{Electronic Journal of Statistics}, 6, 323-353.
#'
#' Kundu, M. G., Harezlak, J., and Randolph, T. W. (2012). Longitudinal
#' functional models with structured penalties (arXiv:1211.4763 [stat.AP]).
#'
#' @seealso \code{\link{pfr}}, \code{\link{smooth.construct.peer.smooth.spec}}
#'
#' @examples
#'
#' \dontrun{
#' #------------------------------------------------------------------------
#' # Example 1: Estimation with D2 penalty
#' #------------------------------------------------------------------------
#'
#' data(DTI)
#' DTI = DTI[which(DTI$case == 1),]
#' fit.D2 = pfr(pasat ~ peer(cca, pentype="D"), data=DTI)
#' plot(fit.D2)
#'
#' #------------------------------------------------------------------------
#' # Example 2: Estimation with structured penalty (need structural
#' #            information about regression function or predictor function)
#' #------------------------------------------------------------------------
#'
#' data(PEER.Sim)
#' data(Q)
#' PEER.Sim1<- subset(PEER.Sim, t==0)
#'
#' # Setting k to max possible value
#' fit.decomp <- pfr(Y ~ peer(W, pentype="Decomp", Q=Q, k=99), data=PEER.Sim1)
#' plot(fit.decomp)
#' }
#'
#'

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
  if("xt" %in% names(dots)) stop("peer() does not accept an xt-argument.")
  xt <- call("list", pentype=pentype, W=substitute(X), phia=phia, L=L,
    Q=substitute(Q))
  lf(X=X, argvals=argvals, bs="peer", xt=xt, ...)
}
