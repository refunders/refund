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



#' Basis constructor for PEER terms
#' 
#' Smooth basis constructor to define structured penalties (Randolph et al.,
#' 2012) for smooth terms.
#' 
#' @param object a \code{peer.smooth.spec} object, usually generated by a 
#'   term \code{s(x, bs="peer")}; see Details.
#' @param data a list containing the data (including any \code{by} variable)
#'   required by this term, with names corresponding to \code{object$term}
#'   (and \code{object$by}). Only the first element of this list is used.
#' @param knots not used, but required by the generic \code{smooth.construct}.
#' 
#' @details The smooth specification object, defined using \code{s()}, should
#'   contain an \code{xt} element. \code{xt} will be a list that contains
#'   additional information needed to specify the penalty. The type of penalty
#'   is indicated by \code{xt$pentype}. There are four types of penalties
#'   available:
#' \enumerate{
#'   \item \code{xt$pentype=="RIDGE"} for a ridge penalty, the default
#'   \item \code{xt$pentype=="D"} for a difference penalty. The order of the
#'     difference penalty is specified by the \code{m} argument of
#'     \code{s()}.
#'   \item \code{xt$pentype=="DECOMP"} for a decomposition-based penalty,
#'     \eqn{bP_Q + a(I-P_Q)}, where \eqn{P_Q = Q^t(QQ^t)^{-1}Q}. The \eqn{Q}
#'     matrix must be specified by \code{xt$Q}, and the scalar \eqn{a} by
#'     \code{xt$phia}. The number of columns of \code{Q} must be equal to the
#'     length of the data. Each row represents a basis function where the
#'     functional predictor is expected to lie, according to prior belief.
#'   \item \code{xt$pentype=="USER"} for a user-specified penalty matrix
#'     \eqn{L}, supplied by \code{xt$L}.
#' }
#' 
#' @return An object of class \code{"peer.smooth"}. See
#'   \code{\link{smooth.construct}} for the elements that this object will
#'   contain.
#' @author Madan Gopal Kundu \email{mgkundu@@iupui.edu} and Jonathan Gellar
#' @seealso \code{\link{peer}}
#' @references
#' Randolph, T. W., Harezlak, J, and Feng, Z. (2012). Structured penalties for
#' functional linear models - partially empirical eigenvectors for regression.
#' \emph{Electronic Journal of Statistics}, 6, 323-353.
#' 
#' 

smooth.construct.peer.smooth.spec <- function(object, data, knots) {
  # Constuctor Method for PEER basis and penalization
  
  #if(!is.null(argvals))
  #  stop("argvals is not supported in the current version of refund.")
  
  K <- length(data[[1]])
  k <- object$bs.dim
  m <- object$p.order
  xt <- object$xt
  pentype <- xt$pentype
  if (is.null(pentype)) pentype="RIDGE"
  
  L <- if (toupper(pentype)=='DECOMP' | toupper(pentype)=='DECOMPOSITION') {
    # Decomposition Penalty
    
    Q <- xt$Q
    phia <- xt$phia
    if (is.null(Q)) stop("Must enter a non-null Q matrix for DECOMP penalty")
    if (is.null(phia)) phia <- 10^3
    else if (!is.numeric(phia)|is.matrix(phia))
      stop("Invalid entry for phia")
    if (ncol(Q) != K) stop("Width of Q matrix must match width of functions")
    
    # Check singularity of Q matrix
    Q <- Q[complete.cases(Q),]    
    Q.eig<- abs(eigen(Q %*% t(Q))$values)
    if(any(Q.eig<1e-12)) stop("Q matrix is singular or near singular")
    
    P_Q <- t(Q) %*% solve(Q %*% t(Q)) %*% Q
    phia*(diag(K)- P_Q) + 1*P_Q
  } else if (toupper(pentype) %in% c("RIDGE", "NONE")) {
    diag(K)
  } else if (toupper(pentype)=="D") {
    # Difference Penalty
    if (is.na(m)) m <- 2
    L <- diag(K)
    for (i in 1:m) L <- diff(L)
    L1 <- L[nrow(L),]
    for (i in 1:m) {
      L1 <- c(0, L1[-K])
      L <- rbind(L, L1)
    }
    rownames(L) <- NULL
    L
  } else if (toupper(pentype)=="USER") {
    L <- xt$L
    if (is.null(L)) stop("Must enter a non-null L matrix for USER penalty")
    if (ncol(L) != K) stop("Width of L matrix must match width of functions")
    
    # Check singularity of L matrix
    L <- L[complete.cases(L),]
    LL<- t(L)%*%L
    LL.eig<- abs(eigen(LL %*% t(LL))$values)
    if(any(LL.eig<1e-12)) stop("L'L matrix is singular or near singular")
    L
  } else {
    stop("Invalid pentype entry for PEER smooth")
  }
  
  # Default k
  if (k<0) k <- K
  if (k==K) {
    v <- diag(K)
  } else {
    W <- xt$W
    v <- svd(data.matrix(W) %*% solve(L))$v[,1:k]
  }
  
  D <- L %*% v
  D <- t(D) %*% D
  D <- (D + t(D))/2
  
  # Return object
  object$X <- v
  object$S <- if (toupper(pentype)=="NONE") list() else list(D)
  ## numerically determine null space dim -- seems safer than relying on m, 
  ## for pentype DECOMP and USER
  object$null.space.dim <- k - qr(D)$rank
  object$rank <- k - object$null.space.dim
  object$df <- k #need this for gamm
  object$argvals <- data[[1]]
  object$v <- v
  class(object) <- "peer.smooth"
  object
}

#' mgcv-style constructor for prediction of PEER terms
#' 
#' @param object a \code{peer.smooth} object created by
#'   \code{\link{smooth.construct.peer.smooth.spec}}, see
#'   \code{\link[mgcv]{smooth.construct}}
#' @param data  see \code{\link[mgcv]{smooth.construct}}
#' @return design matrix for PEER terms
#' @author Jonathan Gellar

Predict.matrix.peer.smooth <- function(object, data) {
  apply(object$v, 2, function(x)
    approx(object$argvals, x, data[[object$term]])$y)
}
