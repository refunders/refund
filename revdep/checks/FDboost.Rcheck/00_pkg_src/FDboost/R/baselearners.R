
#' Functions to compute integration weights
#' 
#' Computes trapezoidal integration weights (Riemann sums) for a functional variable 
#' \code{X1} that has evaluation points \code{xind}.
#' 
#' @param X1 for functional data that is observed on one common grid, 
#' a matrix containing the observations of the functional variable. 
#' For a functional variable that is observed on curve specific grids, a long vector.  
#' @param xind evaluation points (index) of functional variable
#' @param id defaults to \code{NULL}. Only necessary for response in long format. 
#' In this case \code{id} specifies which curves belong together. 
#' @param leftWeight one of \code{c("mean", "first", "zero")}. With left Riemann sums 
#' different assumptions for the weight of the first observation are possible. 
#' The default is to use the mean over all integration weights, \code{"mean"}. 
#' Alternatively one can use the first integration weight, \code{"first"}, or 
#' use the distance to zero, \code{"zero"}. 
#' 
#' @aliases integrationWeightsLeft
#' 
#' @details The function \code{integrationWeights()} computes trapezoidal integration weights, 
#' that are symmetric. Per default those weights are used in the \code{\link{bsignal}}-base-learner. 
#' In the special case of evaluation points (\code{xind}) with equal distances, 
#' all integration weights are equal.
#'   
#' The function \code{integrationWeightsLeft()} computes weights,
#' that take into account only the distance to the prior observation point. 
#' Thus one has to decide what to do with the first observation. 
#' The left weights are adequate for historical effects like in \code{\link{bhist}}.
#'  
#' @seealso \code{\link{bsignal}} and \code{\link{bhist}} for the base-learners. 
#'
#' @examples 
#' ## Example for trapezoidal integration weights
#' xind0 <- seq(0,1,l = 5)
#' xind <- c(0, 0.1, 0.3, 0.7, 1)
#' X1 <- matrix(xind^2, ncol = length(xind0), nrow = 2)
#' 
#' # Regualar observation points
#' integrationWeights(X1, xind0)
#' # Irregular observation points
#' integrationWeights(X1, xind)
#' 
#' # with missing value
#' X1[1,2] <- NA
#' integrationWeights(X1, xind0)
#' integrationWeights(X1, xind)
#' 
#' ## Example for left integration weights
#' xind0 <- seq(0,1,l = 5)
#' xind <- c(0, 0.1, 0.3, 0.7, 1)
#' X1 <- matrix(xind^2, ncol = length(xind0), nrow = 2)
#' 
#' # Regular observation points
#' integrationWeightsLeft(X1, xind0, leftWeight = "mean") 
#' integrationWeightsLeft(X1, xind0, leftWeight = "first") 
#' integrationWeightsLeft(X1, xind0, leftWeight = "zero")
#' 
#' # Irregular observation points
#' integrationWeightsLeft(X1, xind, leftWeight = "mean") 
#' integrationWeightsLeft(X1, xind, leftWeight = "first") 
#' integrationWeightsLeft(X1, xind, leftWeight = "zero")
#' 
#' # obervation points that do not start with 0
#' xind2 <- xind + 0.5
#' integrationWeightsLeft(X1, xind2, leftWeight = "zero")
#'  
#' @export
################################# 
# Trapezoidal integration weights for a functional variable X1 on grid xind
# corresponds to mean of left and right Riemann integration sum
integrationWeights <- function(X1, xind, id = NULL){
  
  if(is.null(id)) if(ncol(X1) != length(xind) ) stop("Dimension of xind and X1 do not match")
  
  # compute integraion weights for irregular data in long format
  if(!is.null(id)){
    Lneu <- tapply(xind, id, FUN = function(x) colMeans(rbind(c(0, diff(x)), c(diff(x), 0))) )
    Lneu <- unlist(Lneu)
    names(Lneu) <- NULL    
    return(Lneu) 
  }  
  
  ## special case that grid has equal distances = regular grid
  if(all( abs(diff(xind) - mean(diff(xind))) < .Machine$double.eps *10^10 )){
    ## use the first difference 
    L <- matrix(diff(xind)[1], nrow=nrow(X1), ncol=ncol(X1)) 
    ##L <- matrix( 1/length(xind)*( max(xind)-min(xind) ), nrow=nrow(X1), ncol=ncol(X1))
    
  }else{ ## case with irregular grid
    
    #Li <- c(diff(xind)[1]/2, diff(xind)[2:(length(xind)-1)], diff(xind)[length(xind)-1]/2)
    
    # Riemann integration weights: mean weights between left and right sum
    # \int^b_a f(t) dt = sum_i (t_i-t_{i-1})*(f(t_i)) (left sum)
    # Li <- c(0,diff(xind))
    
    # trapezoidal
    # \int^b_a f(t) dt = .5* sum_i (t_i - t_{i-1}) f(t_i) + f(t_{i-1}) = 
    #  (t_2 - t_1)/2 * f(a=t_1) + sum^{nx-1}_{i=2} ((t_i - t_{i-1})/2 + (t_{i+1} - t_i)/2) * f(t_i) 
    #      + (t_{nx} - t_{nx-1})/2 * f(b=t_{nx})
    Li <- colMeans(rbind(c(0,diff(xind)), c(diff(xind), 0))) 
    
    # alternative calculation of trapezoidal weights
    #diffs <- diff(xind)
    #nxgrid <- length(xind)
    #Li <- 0.5 * c(diffs[1],  filter(diffs, filter=c(1,1))[-(nxgrid-1)], diffs[(nxgrid-1)] )
    
    L <- matrix(Li, nrow=nrow(X1), ncol=ncol(X1), byrow=TRUE)
    
  }
  
  # taking into account missing values
  if(any(is.na(X1))){
    Lneu <- sapply(1:nrow(X1), function(i){
      x <- X1[i,]
      
      if(!any(is.na(x))){
        l <- L[i, ] # no missing values in curve i
      }else{
        xindL <- xind # lower
        xindL[is.na(x)] <- NA
        xindU <- xindL # upper
        
        if(is.na(xindL[1])){ # first observation is missing
          xindL[1] <- xind[1] - diff(c(xind[1], xind[2])) 
        } 
        if(is.na(xindU[length(xind)])){ # last observation is missing
          xindU[length(xind)] <- xind[length(xind)] + diff(c(xind[length(xind)-1], xind[length(xind)])) 
        }

        xindL <- na.locf(xindL, na.rm=FALSE) # index for lower sum
        xindU <- na.locf(xindU, fromLast=TRUE, na.rm=FALSE) # index for upper sum
        
        l <- colMeans(rbind(c(0,diff(xindL)), c(diff(xindU), 0))) # weight is 0 for missing values
      }
      return(l)
    }
    )
    
    return(t(Lneu))  
    
  }else{ 
    return(L)
  }
}


#### Computes Riemann-weights that only take into account the distance to the previous 
# observation point
# important for bhist()

#' @rdname integrationWeights
#' @export
integrationWeightsLeft <- function(X1, xind, leftWeight = c("first", "mean", "zero")){
  
  if(ncol(X1) != length(xind)) stop("Dimension of xind and X1 do not match")
  # !is.unsorted(xind, strictly = FALSE) # is xind sorted?
  
  leftWeight <- match.arg(leftWeight)
  
  # use lower/left Riemann sum
  Li <- diff(xind)
  # assume  delta(xind_0) = avg. delta
  Li <- switch(leftWeight,
               mean = c(mean(Li), Li), 
               first = c(Li[1], Li), 
               zero = c(xind[1], Li)
  )

  L <- matrix(Li, nrow=nrow(X1), ncol=ncol(X1), byrow=TRUE)
  
  return(L)
}


################################################################################
### syntax for base learners is modified code of the package mboost, see bl.R



################################################################################
################################################################################
# Base-learners for functional covariates 

### hyper parameters for signal baselearner with P-splines
hyper_signal <- function(mf, vary, inS="smooth", knots = 10, boundary.knots = NULL, degree = 3,
                      differences = 1, df = 4, lambda = NULL, center = FALSE,
                      cyclic = FALSE, constraint = "none", deriv = 0L, 
                      Z=NULL, penalty = "ps", check.ident = FALSE,
                      s=NULL) {
  
  knotf <- function(x, knots, boundary.knots) {
    if (is.null(boundary.knots))
      boundary.knots <- range(x, na.rm = TRUE)
    ## <fixme> At the moment only NULL or 2 boundary knots can be specified.
    ## Knot expansion is done automatically on an equidistand grid.</fixme>
    if ((length(boundary.knots) != 2) || !boundary.knots[1] < boundary.knots[2])
      stop("boundary.knots must be a vector (or a list of vectors) ",
           "of length 2 in increasing order")
    if (length(knots) == 1) {
      knots <- seq(from = boundary.knots[1],
                   to = boundary.knots[2], length = knots + 2)
      knots <- knots[2:(length(knots) - 1)]
    }
    list(knots = knots, boundary.knots = boundary.knots)
  }
  
  #   nm <- colnames(mf)[colnames(mf) != vary]
  #   if (is.list(knots)) if(!all(names(knots) %in% nm))
  #     stop("variable names and knot names must be the same")
  #   if (is.list(boundary.knots)) if(!all(names(boundary.knots) %in% nm))
  #     stop("variable names and boundary.knot names must be the same")
  if (!identical(center, FALSE) && cyclic)
    stop("centering of cyclic covariates not yet implemented")
  #    ret <- vector(mode = "list", length = length(nm))
  #    names(ret) <- nm
   
  ret <- knotf(s, knots, boundary.knots)
  
  if (cyclic & constraint != "none")
    stop("constraints not implemented for cyclic B-splines")
  stopifnot(is.numeric(deriv) & length(deriv) == 1)
  
  ## prediction is usually set in/by newX()
  list(knots = ret, degree = degree, differences = differences,
       df = df, lambda = lambda, center = center, cyclic = cyclic,
       Ts_constraint = constraint, deriv = deriv, prediction = FALSE, 
       Z = Z, penalty = penalty, check.ident = check.ident, s=s, inS=inS)
}


### model.matrix for P-splines base-learner of signal matrix mf
### with index/time as attribute
X_bsignal <- function(mf, vary, args) {
  
  stopifnot(is.data.frame(mf))
  xname <- names(mf)
  X1 <- as.matrix(mf)
  xind <- attr(mf[[1]], "signalIndex")
  if(is.null(xind)) xind <- args$s # if the attribute is NULL use the s of the model fit
  
  if(ncol(X1)!=length(xind)) stop(xname, ": Dimension of signal matrix and its index do not match.")
  
  # compute design-matrix in s-direction
  Bs <- switch(args$inS, 
               # B-spline basis of specified degree 
               # "smooth" = bsplines(xind, knots=args$knots$knots, 
               #                   boundary.knots=args$knots$boundary.knots, 
               #                   degree=args$degree),
               "smooth" = mboost_intern(xind, knots = args$knots$knots, 
                                        boundary.knots = args$knots$boundary.knots, 
                                        degree = args$degree, 
                                        fun = "bsplines"),
               "linear" = matrix(c(rep(1, length(xind)), xind), ncol=2),
               "constant"=  matrix(c(rep(1, length(xind))), ncol=1))
  
  colnames(Bs) <- paste(xname, 1:ncol(Bs), sep="")
  
  
  # use cyclic splines
  if (args$cyclic) {
    if(args$inS != "smooth") stop("Cyclic splines are only meaningful for a smooth effect.")
    # Bs <- cbs(xind, knots = args$knots$knots,
    #          boundary.knots = args$knots$boundary.knots,
    #          degree = args$degree)
    Bs <- mboost_intern(xind, knots = args$knots$knots,
                        boundary.knots = args$knots$boundary.knots,
                        degree = args$degree, 
                        fun = "cbs")
  }
  
  colnames(Bs) <- paste(xname, 1:ncol(Bs), sep="")  
  
  ### Penalty matrix: product differences matrix
  if (args$differences > 0){
    if (!args$cyclic) {
      K <- diff(diag(ncol(Bs)), differences = args$differences)
    } else {
      ## cyclic P-splines
      differences <- args$differences
      K <- diff(diag(ncol(Bs) + differences),
                differences = differences)
      tmp <- K[,(1:differences)]   # save first "differences" columns
      K <- K[,-(1:differences)]    # drop first "differences" columns
      indx <- (ncol(Bs) - differences + 1):(ncol(Bs))
      K[,indx] <- K[,indx] + tmp   # add first "differences" columns
    }
  } else {
    if (args$differences != 0)
      stop(sQuote("differences"), " must be an non-negative integer")
    K <- diag(ncol(Bs))
  }
  
  ### penalty matrix is squared difference matrix
  K <- crossprod(K)
  
  if(args$inS != "smooth"){ 
    K <- diag(ncol(Bs))
  }
  
  #----------------------------------
  ### <SB> Calculate constraints if necessary
  ### use the transformation matrix Z if necessary
  ### Check whether integral over trajectories is different, then centering is advisable
  ## as arbitrary constants can be added to the coefficient surface
  if(is.null(args$Z) && 
       all( abs(rowMeans(X1, na.rm = TRUE)-mean(rowMeans(X1, na.rm = TRUE))) < .Machine$double.eps *10^10)){
    C <- t(Bs) %*% rep(1, nrow(Bs))
    Q <- qr.Q(qr(C), complete=TRUE) # orthonormal matrix of QR decomposition
    args$Z <- Q[  , 2:ncol(Q)] # only keep last columns  
  }else{ # nicer solution that Z not produced for prediction with new data with mean 0?
    args$Z <- diag(x=1, ncol=ncol(Bs), nrow=ncol(Bs))
  }
  
  if(!is.null(args$Z)){ 
    ### Transform design and penalty matrix
    Bs <- Bs %*% args$Z
    K <- t(args$Z) %*% K %*% args$Z
  }
  #---------------------------------- 
  
  ### Weighting with matrix of functional covariate
  L <- integrationWeights(X1=X1, xind=xind)
  # Design matrix is product of weighted X1 and basis expansion over xind 
  X <- (L*X1) %*% Bs
  
  colnames(X) <- paste0(xname, 1:ncol(X))
 
  ## see Scheipl and Greven (2016): 
  ## Identifiability in penalized function-on-function regression models  
  if(args$check.ident){
    res_check <- check_ident(X1=X1, L=L, Bs=Bs, K=K, xname=xname, 
                             penalty=args$penalty)
    args$penalty <- res_check$penalty
    args$logCondDs <- res_check$logCondDs
    args$overlapKe <- res_check$overlapKe
    args$maxK <- res_check$maxK
  }
  
  if(args$penalty == "pss"){
    # instead of using 0.1, allow for flexible shrinkage parameter in penalty_pss()?
    K <- penalty_pss(K = K, difference = args$difference, shrink = 0.1)
  }
  
  #####################################################
  ####### K <- crossprod(K) has been computed before!
  if (!identical(args$center, FALSE)) {

    ### L = \Gamma \Omega^1/2 in Section 2.3. of
    ### Fahrmeir et al. (2004, Stat Sinica); "spectralDecomp"
    SVD <- eigen(K, symmetric = TRUE)
    ev <- SVD$vector[, 1:(ncol(X) - args$differences), drop = FALSE]
    ew <- SVD$values[1:(ncol(X) - args$differences), drop = FALSE]
    # penalized part of X: X L (L^t L)^-1
    X <- X %*% ev %*% diag(1/sqrt(ew))
    
    ## unpenalized part of X: 
    ## for differences = 2 gives equivalent results to specifying inS='linear' 
    # X <- X %*% Null(K)
    
    # attributes(X)[c("degree", "knots", "Boundary.knots")] <- tmp
    K <- diag(ncol(X))
  } 
  #####################################################
  
  ## compare specified degrees of freedom to dimension of null space
  if (!is.null(args$df)){
    rns <- ncol(K) - qr(as.matrix(K))$rank # compute rank of null space
    if (rns == args$df)
      warning( sQuote("df"), " equal to rank of null space ",
               "(unpenalized part of P-spline);\n  ",
               "Consider larger value for ", sQuote("df"),
               ## " or set ", sQuote("center = TRUE"), 
               ".", immediate.=TRUE)
    if (rns > args$df)
      stop("not possible to specify ", sQuote("df"),
           " smaller than the rank of the null space\n  ",
           "(unpenalized part of P-spline). Use larger value for ",
           sQuote("df"), 
           ## " or set ", sQuote("center = TRUE"), 
           ".")
  }
  return(list(X = X, K = K, args=args))
}

###############################################################################

#' Base-learners for Functional Covariates
#' 
#' Base-learners that fit effects of functional covariates.  
#' 
#' @param x matrix of functional variable x(s). The functional covariate has to be 
#' supplied as n by <no. of evaluations> matrix, i.e., each row is one functional observation. 
#' @param s vector for the index of the functional variable x(s) giving the 
#' measurement points of the functional covariate. 
#' @param time vector for the index of the functional response y(time) 
#' giving the measurement points of the functional response. 
#' @param index a vector of integers for expanding the covariate in \code{x} 
#' For example, \code{bsignal(X, s, index = index)} is equal to \code{bsignal(X[index,], s)}, 
#' where index is an integer of length greater or equal to \code{NROW(x)}.
#' @param knots either the number of knots or a vector of the positions 
#' of the interior knots (for more details see \code{\link[mboost]{bbs}}).
#' @param boundary.knots boundary points at which to anchor the B-spline basis 
#' (default the range of the data). A vector (of length 2) 
#' for the lower and the upper boundary knot can be specified.
#' @param degree degree of the regression spline.
#' @param differences a non-negative integer, typically 1, 2 or 3. Defaults to 1.  
#' If \code{differences} = \emph{k}, \emph{k}-th-order differences are used as 
#' a penalty (\emph{0}-th order differences specify a ridge penalty).
#' @param df trace of the hat matrix for the base-learner defining the 
#' base-learner complexity. Low values of \code{df} correspond to a 
#' large amount of smoothing and thus to "weaker" base-learners.
#' @param lambda smoothing parameter of the penalty, computed from \code{df} when \code{df} is specified. 
#' @param center See \code{\link[mboost]{bbs}}. 
#' The effect is re-parameterized such that the unpenalized part of the fit is subtracted and only 
#' the penalized effect is fitted, using a spectral decomposition of the penalty matrix.  
#' The unpenalized, parametric part has then to be included in separate 
#' base-learners using \code{bsignal(..., inS = 'constant')} or \code{bsignal(..., inS = 'linear')} 
#' for first (\code{difference = 1}) and second (\code{difference = 2}) order difference penalty respectively. 
#' See the help on the argument \code{center} of \code{\link[mboost]{bbs}}.   
#' @param cyclic if \code{cyclic = TRUE} the fitted coefficient function coincides at the boundaries 
#' (useful for cyclic covariates such as day time etc.).
#' @param Z a transformation matrix for the design-matrix over the index of the covariate.
#' \code{Z} can be calculated as the transformation matrix for a sum-to-zero constraint in the case
#' that all trajectories have the same mean 
#' (then a shift in the coefficient function is not identifiable).
#' @param penalty for \code{bsignal}, by default, \code{penalty = "ps"}, the difference penalty for P-splines is used, 
#' for \code{penalty = "pss"} the penalty matrix is transformed to have full rank, 
#' so called shrinkage approach by Marra and Wood (2011). 
#' For \code{bfpc} the penalty can be either \code{"identity"} for a ridge penalty 
#' (the default) or \code{"inverse"} to use the matrix with the inverse eigenvalues 
#' on the diagonal as penalty matrix or \code{"no"} for no penalty. 
#' @param check.ident use checks for identifiability of the effect, based on Scheipl and Greven (2016) 
#' for linear functional effect using \code{bsignal} and 
#' based on Brockhaus et al. (2017) for historical effects using \code{bhist}
#' @param standard the historical effect can be standardized with a factor. 
#' "no" means no standardization, "time" standardizes with the current value of time and 
#' "length" standardizes with the length of the integral 
#' @param intFun specify the function that is used to compute integration weights in \code{s} 
#' over the functional covariate \eqn{x(s)}
#' @param inS the functional effect can be smooth, linear or constant in s, 
#' which is the index of the functional covariates x(s). 
#' @param inTime the historical effect can be smooth, linear or constant in time, 
#' which is the index of the functional response y(time). 
#' @param limits defaults to \code{"s<=t"} for an historical effect with s<=t;  
#' either one of \code{"s<t"} or \code{"s<=t"} for [l(t), u(t)] = [T1, t]; 
#' otherwise specify limits as a function for integration limits [l(t), u(t)]: 
#' function that takes \eqn{s} as the first and \code{t} as the second argument and returns 
#' \code{TRUE} for combinations of values (s,t) if \eqn{s} falls into the integration range for 
#' the given \eqn{t}.  
#' @param pve proportion of variance explained by the first K functional principal components (FPCs): 
#' used to choose the number of functional principal components (FPCs).
#' @param npc prespecified value for the number K of FPCs (if given, this overrides \code{pve}).
#' @param npc.max maximal number K of FPCs to use; defaults to 15. 
#' @param getEigen save the eigenvalues and eigenvectors, defaults to \code{TRUE}. 
#' 
#' @aliases bconcurrent bhist bfpc 
#' 
#' @details \code{bsignal()} implements a base-learner for functional covariates to  
#' estimate an effect of the form \eqn{\int x_i(s)\beta(s)ds}. Defaults to a cubic  
#' B-spline basis with first difference penalties for \eqn{\beta(s)} and numerical 
#' integration over the entire range by using trapezoidal Riemann weights. 
#' If \code{bsignal()} is used within \code{FDboost()}, the base-learner of 
#' \code{timeformula} is attached, resulting in an effect varying over the index
#' of the response \eqn{\int x_i(s)\beta(s, t)ds} if \code{timeformula = bbs(t)}. 
#' The functional variable must be observed on one common grid \code{s}.  
#' 
#' \code{bconcurrent()} implements a concurrent effect for a functional covariate
#' on a functional response, i.e., an effect of the form \eqn{x_i(t)\beta(t)} for
#' a functional response \eqn{Y_i(t)} and concurrently observed covariate \eqn{x_i(t)}. 
#' \code{bconcurrent()} can only be used if \eqn{Y(t)} and \eqn{x(s)} are observed over
#' the same domain \eqn{s,t \in [T1, T2]}.  
#'
#' \code{bhist()} implements a base-learner for functional covariates with 
#' flexible integration limits \code{l(t)}, \code{r(t)} and the possibility to
#' standardize the effect by \code{1/t} or the length of the integration interval. 
#' The effect is \eqn{stand * \int_{l(t)}^{r_{t}} x(s)\beta(t,s)ds}, where \eqn{stand} is 
#' the chosen standardization which defaults to 1. 
#' The base-learner defaults to a historical effect of the form 
#' \eqn{\int_{T1}^{t} x_i(s)\beta(t,s)ds}, 
#' where \eqn{T1} is the minimal index of \eqn{t} of the response \eqn{Y(t)}. 
#' The functional covariate must be observed on one common grid \code{s}.  
#' See Brockhaus et al. (2017) for details on historical effects.   
#' 
#' \code{bfpc()} is a base-learner for a linear effect of functional covariates based on 
#' functional principal component analysis (FPCA). 
#' For the functional linear effect \eqn{\int x_i(s)\beta(s)ds} the functional covariate 
#' and the coefficient function are both represented by a FPC basis. 
#' The functional covariate
#' \eqn{x(s)} is decomposed into \eqn{x(s) \approx \sum_{k=1}^K \xi_{ik} \Phi_k(s)} using 
#' \code{\link[refund]{fpca.sc}} for the truncated Karhunen-Loeve decomposition. 
#' Then \eqn{\beta(s)} is represented in the function
#' space spanned by \eqn{\Phi_k(s)}, k=1,...,K, see Scheipl et al. (2015) for details. 
#' As penalty matrix, the identity matrix is used. 
#' The implementation is similar to \code{\link[refund]{ffpc}}.  
#' 
#' It is recommended to use centered functional covariates with 
#' \eqn{\sum_i x_i(s) = 0} for all \eqn{s} in \code{bsignal()}-, 
#' \code{bhist()}- and \code{bconcurrent()}-terms. 
#' For centered covariates, the effects are centered per time-point of the response. 
#' If all effects are centered, the functional intercept 
#' can be interpreted as the global mean function. 
#' 
#' The base-learners for functional covariates cannot deal with any missing 
#' values in the covariates.
#' 
#' @return Equally to the base-learners of package \code{mboost}: 
#' 
#' An object of class \code{blg} (base-learner generator) with a 
#' \code{dpp()} function (dpp, data pre-processing). 
#' 
#' The call of \code{dpp()} returns an object of class 
#' \code{bl} (base-learner) with a \code{fit()} function. The call to 
#' \code{fit()} finally returns an object of class \code{bm} (base-model).
#' 
#' @seealso \code{\link{FDboost}} for the model fit. 
#' @keywords models
#' 
#' @references 
#' Brockhaus, S., Scheipl, F., Hothorn, T. and Greven, S. (2015): 
#' The functional linear array model. Statistical Modelling, 15(3), 279-300.
#' 
#' Brockhaus, S., Melcher, M., Leisch, F. and Greven, S. (2017): 
#' Boosting flexible functional regression models with a high number of functional historical effects,  
#' Statistics and Computing, 27(4), 913-926.   
#' 
#' Marra, G. and Wood, S.N. (2011): Practical variable selection for generalized additive models. 
#' Computational Statistics & Data Analysis, 55, 2372-2387.
#' 
#' Scheipl, F., Staicu, A.-M. and Greven, S. (2015): 
#' Functional Additive Mixed Models, Journal of Computational and Graphical Statistics, 24(2), 477-501. 
#' 
#' Scheipl, F. and Greven, S. (2016): Identifiability in penalized function-on-function regression models. 
#' Electronic Journal of Statistics, 10(1), 495-526. 
#'  
#' @examples 
#' ######## Example for scalar-on-function-regression with bsignal()  
#' data("fuelSubset", package = "FDboost")
#' 
#' ## center the functional covariates per observed wavelength
#' fuelSubset$UVVIS <- scale(fuelSubset$UVVIS, scale = FALSE)
#' fuelSubset$NIR <- scale(fuelSubset$NIR, scale = FALSE)
#' 
#' ## to make mboost:::df2lambda() happy (all design matrix entries < 10)
#' ## reduce range of argvals to [0,1] to get smaller integration weights
#' fuelSubset$uvvis.lambda <- with(fuelSubset, (uvvis.lambda - min(uvvis.lambda)) /
#'                                   (max(uvvis.lambda) - min(uvvis.lambda) ))
#' fuelSubset$nir.lambda <- with(fuelSubset, (nir.lambda - min(nir.lambda)) /
#'                                 (max(nir.lambda) - min(nir.lambda) ))
#'
#' ## model fit with scalar response and two functional linear effects 
#' ## include no intercept 
#' ## as all base-learners are centered around 0 
#' mod2 <- FDboost(heatan ~ bsignal(UVVIS, uvvis.lambda, knots = 40, df = 4, check.ident = FALSE) 
#'                + bsignal(NIR, nir.lambda, knots = 40, df=4, check.ident = FALSE), 
#'                timeformula = NULL, data = fuelSubset) 
#' summary(mod2) 
#' ## plot(mod2)
#' 
#' 
#' ###############################################
#' ### data simulation like in manual of pffr::ff
#' 
#' if(require(refund)){
#' 
#' #########
#' # model with linear functional effect, use bsignal()
#' # Y(t) = f(t) + \int X1(s)\beta(s,t)ds + eps
#' set.seed(2121)
#' data1 <- pffrSim(scenario = "ff", n = 40)
#' data1$X1 <- scale(data1$X1, scale = FALSE)
#' dat_list <- as.list(data1)
#' dat_list$t <- attr(data1, "yindex")
#' dat_list$s <- attr(data1, "xindex")
#' 
#' ## model fit by FDboost 
#' m1 <- FDboost(Y ~ 1 + bsignal(x = X1, s = s, knots = 5), 
#'               timeformula = ~ bbs(t, knots = 5), data = dat_list, 
#'               control = boost_control(mstop = 21))
#' 
#' ## search optimal mSTOP
#' \dontrun{
#'   set.seed(123)
#'   cv <- validateFDboost(m1, grid = 1:100) # 21 iterations
#' }
#' 
#' ## model fit by pffr
#' t <- attr(data1, "yindex")
#' s <- attr(data1, "xindex")
#' m1_pffr <- pffr(Y ~ ff(X1, xind = s), yind = t, data = data1)
#' 
#' \dontrun{
#'   par(mfrow = c(2, 2))
#'   plot(m1, which = 1); plot(m1, which = 2) 
#'   plot(m1_pffr, select = 1, shift = m1_pffr$coefficients["(Intercept)"]) 
#'   plot(m1_pffr, select = 2)
#' }
#' 
#' 
#' ############################################
#' # model with functional historical effect, use bhist() 
#' # Y(t) = f(t)  + \int_0^t X1(s)\beta(s,t)ds + eps
#' set.seed(2121)
#' mylimits <- function(s, t){
#'   (s < t) | (s == t)
#' }
#' data2 <- pffrSim(scenario = "ff", n = 40, limits = mylimits)
#' data2$X1 <- scale(data2$X1, scale = FALSE)
#' dat2_list <- as.list(data2)
#' dat2_list$t <- attr(data2, "yindex")
#' dat2_list$s <- attr(data2, "xindex")
#' 
#' ## model fit by FDboost 
#' m2 <- FDboost(Y ~ 1 + bhist(x = X1, s = s, time = t, knots = 5), 
#'               timeformula = ~ bbs(t, knots = 5), data = dat2_list, 
#'               control = boost_control(mstop = 40))
#'               
#' ## search optimal mSTOP
#' \dontrun{
#'   set.seed(123)
#'   cv2 <- validateFDboost(m2, grid = 1:100) # 40 iterations
#' }               
#' 
#' ## model fit by pffr
#' t <- attr(data2, "yindex")
#' s <- attr(data2, "xindex")
#' m2_pffr <- pffr(Y ~ ff(X1, xind = s, limits = "s<=t"), yind = t, data = data2)
#' 
#' \dontrun{
#' par(mfrow = c(2, 2))
#' plot(m2, which = 1); plot(m2, which = 2)
#' ## plot of smooth intercept does not contain m1_pffr$coefficients["(Intercept)"]
#' plot(m2_pffr, select = 1, shift = m2_pffr$coefficients["(Intercept)"]) 
#' plot(m2_pffr, select = 2) 
#' 
#' }
#' 
#' 
#' }
#' 
#' 
#' @export
### P-spline base-learner for signal matrix with index vector
bsignal <- function(x, s, index = NULL, inS = c("smooth", "linear", "constant"), #by = NULL,
                    knots = 10, boundary.knots = NULL, degree = 3, differences = 1, df = 4, 
                    lambda = NULL, center = FALSE, 
                    cyclic = FALSE, Z = NULL, 
                    penalty = c("ps","pss"), check.ident = FALSE
){
  
  if (!is.null(lambda)) df <- NULL
  
  cll <- match.call()
  cll[[1]] <- as.name("bsignal")
  penalty <- match.arg(penalty)
  inS <- match.arg(inS)
  #print(cll)
  
  # if(!isMATRIX(x)) stop("signal has to be a matrix")
  if( ! mboost_intern(x, fun = "isMATRIX")) stop("signal has to be a matrix")
  
  varnames <- all.vars(cll)
  #   if(length(mfL)==1){ 
  #     mfL[[2]] <- 1:ncol(mfL[[1]]); cll[[3]] <- "xind" 
  #     varnames <- c(all.vars(cll), "xindDefault")
  #   }
  
  # Reshape mfL so that it is the dataframe of the signal with the index as attribute
  xname <- varnames[1]
  indname <- varnames[2] 
  if(is.null(colnames(x))) colnames(x) <- paste(xname, 1:ncol(x), sep="_")
  attr(x, "signalIndex") <- s
  attr(x, "xname") <- xname
  attr(x, "indname") <- indname
  
  mf <- data.frame("z"=I(x))
  names(mf) <- xname
  
  if(!all( abs(colMeans(x, na.rm = TRUE)) < .Machine$double.eps*10^10)){
    message(xname, " is not centered per column, inducing a non-centered effect.")
  }
  
  if(is.null(Z) && 
       all( abs(rowMeans(x, na.rm = TRUE)-mean(rowMeans(x, na.rm = TRUE))) < .Machine$double.eps *10^10)){
    message(paste("All trajectories in ", xname, " have the same mean. Coefficient function is centered.", sep=""))
  }
  
  #   mf <- mfL
  #   names(mf) <- varnames
  
  vary <- ""
  
  # CC <- all(Complete.cases(mf))
  CC <- all(mboost_intern(mf, fun = "Complete.cases"))
  if (!CC)
    warning("base-learner contains missing values;\n",
            "missing values are excluded per base-learner, ",
            "i.e., base-learners may depend on different",
            " numbers of observations.")
  
  #index <- NULL
  
  ## call X_bsignal in oder to compute parameter settings, e.g. 
  ## the transformation matrix Z, shrinkage penalty, identifiability problems...
  temp <- X_bsignal(mf, vary, 
                    args = hyper_signal(mf, vary, inS=inS, knots = knots, 
                                        boundary.knots = boundary.knots, degree = degree, 
                                        differences = differences,
                                        df = df, lambda = lambda, center = center, cyclic = cyclic,
                                        Z = Z, penalty = penalty, check.ident = check.ident,
                                        s = s))
  temp$args$check.ident <- FALSE
  
  ret <- list(model.frame = function() 
    if (is.null(index)) return(mf) else{
      mftemp <- mf
      mf <- mftemp[index,,drop = FALSE] # this is necessary to pass the attributes
      attributes(mftemp[,xname])
      attr(mf[,xname], "signalIndex") <- attr(mftemp[,xname], "signalIndex")
      attr(mf[,xname], "xname") <- attr(mftemp[,xname], "xname")
      attr(mf[,xname], "indname") <- attr(mftemp[,xname], "indname")
      return(mf)
    } ,
    get_call = function(){
      cll <- deparse(cll, width.cutoff=500L)
      if (length(cll) > 1)
        cll <- paste(cll, collapse="")
      cll
    },
    get_data = function() mf,
    get_index = function() index,
    get_vary = function() vary,
    get_names = function(){
      attr(xname, "indname") <- indname 
      xname
      #c(xname, indname)
    }, #colnames(mf),
    set_names = function(value) {
      #if(length(value) != length(colnames(mf)))
      if(length(value) != names(mf[1]))
        stop(sQuote("value"), " must have same length as ",
             sQuote("names(mf[1])"))
      for (i in 1:length(value)){
        cll[[i+1]] <<- as.name(value[i])
      }
      attr(mf, "names") <<- value
    })
  class(ret) <- "blg"
  
  #print("bsignal")
  #print(Z)
  
  # ret$dpp <- bl_lin(ret, Xfun = X_bsignal, args = temp$args)
  ret$dpp <- mboost_intern(ret, Xfun = X_bsignal, args = temp$args, fun = "bl_lin")
  
  ## function that comoutes a design matrix such that new_des %*% hat{theta} = beta(s)
  ## use ng equally spaced observation points 
  ret$get_des <- function(ng = 40){

    ## use a new grid of s-values with ng grid points 
    s_grid <- seq(min(s), max(s), l = ng)
    ## matrix with inverse integraion weights
    dummyX <- I(diag(length(s_grid)) /  integrationWeights(diag(length(s_grid)), s_grid ))

    ## setup for X_signal 
    attr(dummyX, "signalIndex") <- s_grid
    attr(dummyX, "xname") <- xname
    attr(dummyX, "indname") <- indname
    
    new_mf <- data.frame("z" = I(dummyX))
    names(new_mf) <- xname
     
    ## use vary and args like in bl  
    new_des <- X_bsignal(mf = new_mf, vary = vary, args = temp$args)$X
    
    ## give arguments to new_des for easier use in the plot-function 
    attr(new_des, "x") <- s_grid
    attr(new_des, "xlab") <- indname
    attr(new_des, "varname") <- xname
    
    return(new_des)
  }

  # rm(temp)
  temp$X <- NULL
  temp$K <- NULL
  
  return(ret)
}

# testX <- I(matrix(rnorm(40), ncol=5))
# s <- seq(0,1,l=5)
# test <- bsignal(testX, s)
# test$get_names()
# test$get_data()
# names(test$dpp(rep(1,nrow(testX))))
########################


#################################
# Base-learner for concurrent effect of functional covariate

### model.matrix for P-splines base-learner of signal matrix mf
X_conc <- function(mf, vary, args) {
  
  stopifnot(is.data.frame(mf))
  xname <- names(mf)
  X1 <- as.matrix(mf)
  class(X1) <- "matrix"
  xind <- attr(mf[[1]], "signalIndex")
  yind <- attr(mf[[1]], "indexY")
  if(is.null(xind)) xind <- args$s # if the attribute is NULL use the s of the model fit 
  if(is.null(yind)) yind <- args$time # if the attribute is NULL use the time of the model fit  
  
  nobs <- nrow(X1)
  
  # get id-variable
  id <- attr(mf[,xname], "id") # for data in long format
  # id is NULL for regular response
  # id has values 1, 2, 3, ... for response in long format
  
  ## <FIXME> is that line still necessary? 
  ## important for prediction, otherwise id=NULL and yind is multiplied accordingly
  if(is.null(id)) id <- 1:nrow(X1)
  
  ## check yind 
  if(args$format=="long" && length(yind)!=length(id)) stop(xname, ": Index of response and id do not have the same length")
  ## check dimensions of s and x(s)
  if(args$format=="wide" && ncol(X1)!=length(xind)) stop(xname, ": Dimension of signal matrix and its index do not match.")
  if(args$format=="long" && nrow(X1)!=length(xind)) stop(xname, ": Dimension of signal matrix and its index do not match.")
  
  # compute design-matrix in s-direction
  Bs <- switch(args$inS, 
               # B-spline basis of specified degree 
               # "smooth" = bsplines(xind, knots=args$knots$s$knots, 
               #                   boundary.knots=args$knots$s$boundary.knots, 
               #                   degree=args$degree),
               "smooth" = mboost_intern(xind, knots = args$knots$s$knots, 
                                        boundary.knots = args$knots$s$boundary.knots, 
                                        degree = args$degree, 
                                        fun = "bsplines"),
               "linear" = matrix(c(rep(1, length(xind)), xind), ncol = 2),
               "constant"=  matrix(c(rep(1, length(xind))), ncol = 1))
  
  # use cyclic splines
  if (args$cyclic) {
    if(args$inS != "smooth") stop("Cyclic splines are only meaningful for a smooth effect.")
    Bs <- mboost_intern(xind, knots = args$knots$s$knots,
                        boundary.knots = args$knots$s$boundary.knots,
                        degree = args$degree, 
                        fun = "cbs")
  }
  
  colnames(Bs) <- paste(xname, 1:ncol(Bs), sep="")
    
  # set up design matrix for concurrent model
  if(args$format=="wide"){
    listCol <- list()
    for(i in 1:ncol(X1)){
      listCol[[i]] <- X1[,i]
    }
    X1des <- as.matrix(bdiag(listCol))
    # Design matrix is product of expanded X1 and basis expansion over xind 
    X <- (X1des) %*% Bs
    rm(X1des, listCol)
  }else{
    # Design matrix contains rows x_i(t_{ig_i})*Bs[at row t_{ig_i},]
    X <- X1[,1]*Bs
  }
  
  ### Penalty matrix: product differences matrix
  differenceMatrix <- diff(diag(ncol(X)), differences = args$differences)
  K <- crossprod(differenceMatrix)
  
  ## compare specified degrees of freedom to dimension of null space
  if (!is.null(args$df)){
    rns <- ncol(K) - qr(as.matrix(K))$rank # compute rank of null space
    if (rns == args$df)
      warning( sQuote("df"), " equal to rank of null space ",
               "(unpenalized part of P-spline);\n  ",
               "Consider larger value for ", sQuote("df"),
               ## " or set ", sQuote("center = TRUE"), 
               ".", immediate.=TRUE)
    if (rns > args$df)
      stop("not possible to specify ", sQuote("df"),
           " smaller than the rank of the null space\n  ",
           "(unpenalized part of P-spline). Use larger value for ",
           sQuote("df"), 
           ## " or set ", sQuote("center = TRUE"), 
           ".")
  }
  
  # tidy up workspace 
  rm(X1)
  
  return(list(X = X, K = K, args = args))
}

#' @rdname bsignal
#' @export
### P-spline base learner for signal matrix with index vector
bconcurrent <- function(x, s, time, index = NULL, #by = NULL, 
                        knots = 10, boundary.knots = NULL, degree = 3, differences = 1, df = 4, 
                        lambda = NULL, #center = FALSE, 
                        cyclic = FALSE
){
  
  if (!is.null(lambda)) df <- NULL
  
  cll <- match.call()
  cll[[1]] <- as.name("bconcurrent")
  #print(cll)
  
  if(!mboost_intern(x, fun = "isMATRIX") && is.null(index)) stop("signal has to be a matrix for regular response")
  if( mboost_intern(x, fun = "isMATRIX") && NCOL(x)!=length(s)) stop("Dimension of x and s do not match.")
  if(!mboost_intern(x, fun = "isMATRIX") && length(x)!=length(s)) stop("Dimension of x and s do not match.")
  
  varnames <- all.vars(cll)
  
  if(!is.atomic(s)) stop("index of signal has to be a vector")
  if(!is.atomic(time)) stop("index of response has to be a vector")
  
  if(substitute(time)==substitute(s)){
    #warning("Do not use the same variable t as time-variable in y(t) and in x(t).")
    warning("Do not use the same variable for s and time in bconcurrent().")
  }
  
  # compare range of index signal and index response
  # the index of the signal s, has to contain all values of time
  if( !all(s %in% time) ) stop("Index s of functional variable has to contain all values of time.")
  
  # Reshape mfL so that it is the dataframe of the signal with the index as attribute
  xname <- varnames[1]
  indname <- varnames[2]
  indnameY <- varnames[3]
  if(length(varnames)==2) indnameY <- varnames[2]
  attr(x, "indexY") <- time
  attr(x, "indnameY") <- indnameY
  attr(x, "id") <- index
  
  if(mboost_intern(x, fun = "isMATRIX") && 
     is.null(colnames(x))) colnames(x) <- paste(xname, 1:ncol(x), sep="_")
  attr(x, "signalIndex") <- s
  attr(x, "xname") <- xname
  attr(x, "indname") <- indname 
  
  mf <- data.frame("z"=I(x))
  names(mf) <- xname 
  
  #   if(all(round(colSums(mf, na.rm = TRUE), 4)!=0)){
  #     warning(xname, " is not centered. 
  #     Functional covariates should be mean-centered in each measurement point.")
  #   }
  
  #   mf <- mfL
  #   names(mf) <- varnames
  
  vary <- ""
  
  # CC <- all(Complete.cases(mf))
  CC <- all(mboost_intern(mf, fun = "Complete.cases"))
  if (!CC)
    warning("base-learner contains missing values;\n",
            "missing values are excluded per base-learner, ",
            "i.e., base-learners may depend on different",
            " numbers of observations.")
  
  #index <- NULL 
  
  if(is.null(index)){
    ### X_conc for data in wide format with regular response
    temp <- X_conc(mf, vary, 
                   args = hyper_hist(mf, vary, knots = knots, boundary.knots = boundary.knots, 
                                     degree = degree, differences = differences,
                                     df = df, lambda = lambda, center = FALSE, cyclic = cyclic,
                                     s = s, time=time, limits = NULL, 
                                     inS = "smooth", inTime = "smooth", 
                                     penalty = "ps", check.ident = FALSE, 
                                     format="wide"))
  }else{
    ### X_conc for data in long format with irregular response
    temp <- X_conc(mf, vary, 
                   args = hyper_hist(mf, vary, knots = knots, boundary.knots = boundary.knots, 
                                     degree = degree, differences = differences,
                                     df = df, lambda = lambda, center = FALSE, cyclic = cyclic,
                                     s = s, time=time, limits = NULL, 
                                     inS = "smooth", inTime = "smooth", 
                                     penalty = "ps", check.ident = FALSE, 
                                     format="long"))
  }
  
  ret <- list(model.frame = function() 
    if (is.null(index)) return(mf) else{
      mftemp <- mf
      mf <- mftemp[index,,drop = FALSE] # this is necessary to pass the attributes
      attributes(mftemp[,xname])
      attr(mf[,xname], "signalIndex") <- attr(mftemp[,xname], "signalIndex")
      attr(mf[,xname], "xname") <- attr(mftemp[,xname], "xname")
      attr(mf[,xname], "indname") <- attr(mftemp[,xname], "indname")
      attr(mf[,xname], "indexY") <- attr(mftemp[,xname], "indexY")
      attr(mf[,xname], "indnameY") <- attr(mftemp[,xname], "indnameY")
      attr(mf[,xname], "id") <- attr(mftemp[,xname], "id")
      return(mf)
    },
    get_call = function(){
      cll <- deparse(cll, width.cutoff=500L)
      if (length(cll) > 1)
        cll <- paste(cll, collapse="")
      cll
    },
    get_data = function() mf,
    ## set index to NULL, as the index is treated within X_conc()
    ##get_index = function() index, 
    get_index = function() NULL,
    get_vary = function() vary,
    get_names = function(){
      attr(xname, "indname") <- indname 
      attr(xname, "indnameY") <- indnameY
      xname 
    }, #colnames(mf),
    set_names = function(value) {
      #if(length(value) != length(colnames(mf)))
      if(length(value) != names(mf[1]))
        stop(sQuote("value"), " must have same length as ",
             sQuote("names(mf[1])"))
      for (i in 1:length(value)){
        cll[[i+1]] <<- as.name(value[i])
      }
      attr(mf, "names") <<- value
    })
  class(ret) <- "blg"
  
  # ret$dpp <- bl_lin(ret, Xfun = X_conc, args = temp$args)
  ret$dpp <- mboost_intern(ret, Xfun = X_conc, args = temp$args, fun = "bl_lin")
  
  return(ret)
}

# testX <- I(matrix(rnorm(40), ncol=5))
# s <- seq(0,1,l=5)
# test <- bconcurrent(testX, s)
# test$get_names()
# test$get_data()
# names(test$dpp(rep(1,nrow(testX))))






#################################
#### Base-learner for historic effect of functional covariate
### with integral over specific limits, e.g. s<=t

### hyper parameters for signal baselearner with P-splines
hyper_hist <- function(mf, vary, knots = 10, boundary.knots = NULL, degree = 3,
                         differences = 1, df = 4, lambda = NULL, center = FALSE,
                         cyclic = FALSE, constraint = "none", deriv = 0L, 
                         Z=NULL, s=NULL, time=NULL, limits=NULL, 
                         standard="no", intFun=integrationWeightsLeft,  
                         inS="smooth", inTime="smooth", 
                         penalty = "ps", check.ident = FALSE, 
                         format="long") {
  
  knotf <- function(x, knots, boundary.knots) {
    if (is.null(boundary.knots))
      boundary.knots <- range(x, na.rm = TRUE)
    ## <fixme> At the moment only NULL or 2 boundary knots can be specified.
    ## Knot expansion is done automatically on an equidistand grid.</fixme>
    if ((length(boundary.knots) != 2) || !boundary.knots[1] < boundary.knots[2])
      stop("boundary.knots must be a vector (or a list of vectors) ",
           "of length 2 in increasing order")
    if (length(knots) == 1) {
      knots <- seq(from = boundary.knots[1],
                   to = boundary.knots[2], length = knots + 2)
      knots <- knots[2:(length(knots) - 1)]
    }
    list(knots = knots, boundary.knots = boundary.knots)
  }
  
  #   nm <- colnames(mf)[colnames(mf) != vary]
  #   if (is.list(knots)) if(!all(names(knots) %in% nm))
  #     stop("variable names and knot names must be the same")
  #   if (is.list(boundary.knots)) if(!all(names(boundary.knots) %in% nm))
  #     stop("variable names and boundary.knot names must be the same")
  if (!identical(center, FALSE) && cyclic)
    stop("centering of cyclic covariates not yet implemented")
  #    ret <- vector(mode = "list", length = length(nm))
  #    names(ret) <- nm
  
  ret <- vector(mode = "list", length = 2)
  indVars <- c("s","time")
  names(ret) <- indVars
  for (n in 1:2) ret[[n]] <- knotf(get(indVars[[n]]), if (is.list(knots)) 
    knots[[n]]
    else knots, if (is.list(boundary.knots)) 
      boundary.knots[[n]]
    else boundary.knots)
  
  if (cyclic & constraint != "none")
    stop("constraints not implemented for cyclic B-splines")
  stopifnot(is.numeric(deriv) & length(deriv) == 1)
  
  ## prediction is usually set in/by newX()
  list(knots = ret, degree = degree, differences = differences,
       df = df, lambda = lambda, center = center, cyclic = cyclic,
       Ts_constraint = constraint, deriv = deriv, prediction = FALSE, 
       Z = Z, s = s, time = time, limits = limits, 
       standard = standard, intFun = intFun, 
       inS = inS, inTime = inTime, 
       penalty = penalty, check.ident = check.ident, format = format)
}


### model.matrix for P-splines base-learner of signal matrix mf
### for response observed over a common grid, args$format="wide"
### or irregularly observed reponse, args$format="long" 
X_hist <- function(mf, vary, args) {
  
  stopifnot(is.data.frame(mf))
  xname <- names(mf)
  X1 <- as.matrix(mf)
  class(X1) <- "matrix"
  xind <- attr(mf[[1]], "signalIndex")
  yind <- attr(mf[[1]], "indexY")
  
  if(is.null(xind)) xind <- args$s # if the attribute is NULL use the s of the model fit
  if(is.null(yind)) yind <- args$time # if the attribute is NULL use the time of the model fit  
  
  nobs <- nrow(X1)
  
  # get id-variable
  id <- attr(mf[,xname], "id") # for data in long format
  # id is NULL for regular response
  # id has values 1, 2, 3, ... for response in long format
  
  ## <FIXME> is that line still necessary? should it be there in long and wide format?
  ###### EXTRA LINE in comparison to X_hist
  ## important for prediction, otherwise id=NULL and yind is multiplied accordingly
  if(is.null(id)) id <- 1:nrow(X1)
  
  ## check yind 
  if(args$format=="long" && length(yind)!=length(id)) stop(xname, ": Index of response and id do not have the same length")
  ## check dimensions of s and x(s)
  if(ncol(X1)!=length(xind)) stop(xname, ": Dimension of signal matrix and its index do not match.")
  
  # compute design-matrix in s-direction
  Bs <- switch(args$inS, 
               # B-spline basis of specified degree 
               #"smooth" = bsplines(xind, knots=args$knots$s$knots, 
               #                    boundary.knots=args$knots$s$boundary.knots, 
               #                    degree=args$degree),
               "smooth" = mboost_intern(xind, knots = args$knots$s$knots, 
                                        boundary.knots = args$knots$s$boundary.knots, 
                                        degree = args$degree, 
                                        fun = "bsplines"),
               "linear" = matrix(c(rep(1, length(xind)), xind), ncol = 2),
               "constant"=  matrix(c(rep(1, length(xind))), ncol = 1))
  
  colnames(Bs) <- paste(xname, 1:ncol(Bs), sep="")
  
  # integration weights 
  L <- args$intFun(X1=X1, xind=xind)
  # print(L[1,])
  
  ## Weighting with matrix of functional covariate
  #X1 <- L*X1 ## -> do the integration weights more sophisticated!!
  
  #   # set up design matrix for historical model and s<=t with s and t equal to xind
  #   # expand matrix of original observations to lower triangular matrix 
  #   X1des0 <- matrix(0, ncol=ncol(X1), nrow=ncol(X1)*nrow(X1))
  #   for(i in 1:ncol(X1des0)){
  #     #print(nrow(X1)*(i-1)+1)
  #     X1des0[(nrow(X1)*(i-1)+1):nrow(X1des0) ,i] <- X1[,i] # use fun. variable * integration weights
  #   }
  
  ## set up design matrix for historical model according to args$limits()
  # use the argument limits (Code taken of function ff(), package refund)
  limits <- args$limits
  if (!is.null(limits)) {
    if (!is.function(limits)) {
      if (!(limits %in% c("s<t", "s<=t"))) {
        stop("supplied <limits> argument unknown")
      }
      if (limits == "s<t") {
        limits <- function(s, t) {
          s < t
        }
      }
      else {
        if (limits == "s<=t") {
          limits <- function(s, t) {
            (s < t) | (s == t)
          }
        }
      }
    }
  }else{
    stop("<limits> argument cannot be NULL.")
  }
  
  ## save the limits function in the arguments
  args$limits <- limits
  
  ### use function limits to set up design matrix according to function limits 
  ### by setting 0 at the time-points that should not be used
  if(args$format == "wide"){
    ## expand yind by replication to the yind of all observations together
    ind0 <- !t(outer( xind, rep(yind, each=nobs), limits) )
    yindHelp <- rep(yind, each=nobs)
  }else{
    ## yind is over all observations in long format
    ind0 <- !t(outer( xind, yind, limits) )
    yindHelp <- yind
  }  
  
  ### Compute the design matrix as sparse or normal matrix 
  ### depending on dimensions of the final design matrix
  MATRIX <- any(c(nrow(ind0), ncol(Bs)) > c(50, 50)) #MATRIX <- any(dim(X) > c(500, 50))
  MATRIX <- MATRIX && options("mboost_useMatrix")$mboost_useMatrix 

  if(MATRIX){
    #message("use sparse matrix in X_hist")
    diag <- Diagonal
    ###### more efficient construction of X1des directly as sparse matrix 
    #     ### compute the design matrix as sparse matrix
    #     if(args$format == "wide"){
    #       tempIndexDesign <- which(!ind0, arr.ind=TRUE)
    #       tempIndexX1 <- cbind(rep(1:nobs, length.out=nrow(tempIndexDesign)), tempIndexDesign[,2] )
    #       X1des <- sparseMatrix(i=tempIndexDesign[,1], j=tempIndexDesign[,2],
    #                             x=X1[tempIndexX1], dims=dim(ind0))  
    #       # object.size(X1des)
    #       rm(tempIndexX1, tempIndexDesign)
    #     }else{ # long format
    #       tempj <- unlist(apply(!ind0, 1, which)) # in which columns are the values? 
    #       ## i: row numbers: one row number per observation of response, 
    #       #     repeat the row number for each entry
    #       X1des <- sparseMatrix(i=rep(1:length(id), times=rowSums(!ind0)), j=tempj,
    #                             x=X1[cbind(rep(id, t=rowSums(!ind0)), tempj)], dims=dim(ind0))
    #       # object.size(X1des)
    #       rm(tempj)       
    #     }
    
    ###### instead: build the matrix as dense matrix and convert it into a sparse matrix
    if(args$format == "wide"){
      ### expand the design matrix for all observations (yind is equal for all observations!)
      ### the response is a vector (y1(t1), y2(t1), ... , yn(t1), yn(tG))
      X1des <- X1[rep(1:nobs, times=length(yind)), ]
    } else{ # yind is over all observations in long format
      X1des <- X1[id, ] 
    }
    X1des[ind0] <- 0    
    X1des <- Matrix(X1des, sparse=TRUE) # convert into sparse matrix
    
  }else{ # small matrices: do not use Matrix
    if(args$format == "wide"){
      ### expand the design matrix for all observations (yind is equal for all observations!)
      ### the response is a vector (y1(t1), y2(t1), ... , yn(t1), yn(tG))
      X1des <- X1[rep(1:nobs, times=length(yind)), ]
    } else{ # yind is over all observations in long format
      X1des <- X1[id, ] 
    }
    X1des[ind0] <- 0
  }
  
  ## set up matrix with adequate integration and standardization weights
  ## start with a matrix of integration weights
  ## case of "no" standardization
  Lnew <- args$intFun(X1des, xind)
  Lnew[ind0] <- 0
  
  ## Standardize with exact length of integration interval
  ##  (1/t-t0) \int_{t0}^t f(s) ds
  if(args$standard == "length"){
    ## use fundamental theorem of calculus 
    ## \lim t->t0- (1/t-t0) \int_{t0}^t f(s) ds = f(t0)
    ## -> integration weight in s-direction should be 1
    ## integration weights in s-direction always sum exactly to 1, 
    ## good for small number of observations!
    args$vecStand <- rowSums(Lnew)
    args$vecStand[args$vecStand==0] <- 1 ## cannnot divide 0/0, instead divide 0/1
    Lnew <- Lnew * 1/args$vecStand
  } 
  
  ## use time of current observation for standardization
  ##  (1/t) \int_{t0}^t f(s) ds
  if(args$standard=="time"){
    if(any(yindHelp <= 0)) stop("For standardization with time, time must be positive.")
    ## Lnew <- matrix(1, ncol=ncol(X1des), nrow=nrow(X1des))
    ## Lnew[ind0] <- 0  
    ## use fundamental theorem of calculus 
    ## \lim t->0+ (1/t) \int_0^t f(s) ds = f(0), if necessary
    ## (as previously X*L, use now X*(1/L) for cases with one single point)
    yindHelp[yindHelp==0] <- L[1,1] # impossible! 
    # standFact <- 1/yindHelp 
    args$vecStand <- yindHelp
    Lnew <- Lnew * 1/yindHelp 
  }
  ## print(round(Lnew, 2))
  ## print(rowSums(Lnew))
  ## print(args$vecStand)
  
  # multiply design matrix with integration weights and standardization weights
  X1des <- X1des * Lnew
  
  # Design matrix is product of expanded X1 and basis expansion over xind 
  X1des <- X1des %*% Bs
  
  
  ## see Scheipl and Greven (2016): Identifiability in penalized function-on-function regression models  
  if(args$check.ident && args$inS == "smooth"){
    K1 <- diff(diag(ncol(Bs)), differences = args$differences)
    K1 <- crossprod(K1)
    # use the limits function to compute check measures on corresponding subsets of x(s) and B_j
    res_check <- check_ident(X1 = X1, L = L, Bs = Bs, K = K1, xname = xname, 
                             penalty = args$penalty, 
                             limits = args$limits, 
                             yind = yindHelp, id = id, # yind in long format
                             X1des = X1des, ind0 = ind0, xind = xind)
    args$penalty <- res_check$penalty
    args$logCondDs <- res_check$logCondDs
    args$logCondDs_hist <- res_check$logCondDs_hist
    args$overlapKe <- res_check$overlapKe
    args$cumOverlapKe <- res_check$cumOverlapKe
    args$maxK <- res_check$maxK
  }
  
  # wide: design matrix over index of response for one response
  # long: design matrix over index of response (yind has long format!)
  Bt <- switch(args$inTime, 
               # B-spline basis of specified degree 
               #"smooth" = bsplines(yind, knots=args$knots$time$knots, 
               #                    boundary.knots=args$knots$time$boundary.knots, 
               #                    degree=args$degree),
               "smooth" = mboost_intern(yind, knots = args$knots$time$knots, 
                                        boundary.knots = args$knots$time$boundary.knots, 
                                        degree = args$degree, 
                                        fun = "bsplines"),
               "linear" = matrix(c(rep(1, length(yind)), yind), ncol = 2),
               "constant"=  matrix(c(rep(1, length(yind))), ncol = 1))
    
  # stack design-matrix of response nobs times in wide format
  if(args$format == "wide"){
    Bt <- Bt[rep(1:length(yind), each=nobs), ]
  }
  
  if(! mboost_intern(Bt, fun = "isMATRIX") ) Bt <- matrix(Bt, ncol=1)
    
  # calculate row-tensor
  # X <- (X1 %x% t(rep(1, ncol(X2))) ) * ( t(rep(1, ncol(X1))) %x% X2  )
  dimnames(Bt) <- NULL # otherwise warning "dimnames [2] mismatch..."
  X <- X1des[,rep(1:ncol(Bs), each=ncol(Bt))] * Bt[,rep(1:ncol(Bt), times=ncol(Bs))]
  
  if(! mboost_intern(X, fun = "isMATRIX") ) X <- matrix(X, ncol=1)
  
  colnames(X) <- paste0(xname, 1:ncol(X))
  
  ### Penalty matrix: product differences matrix for smooth effect
  if(args$inS == "smooth"){
    K1 <- diff(diag(ncol(Bs)), differences = args$differences)
    K1 <- crossprod(K1)    
    if(args$penalty == "pss"){
      # instead of using 0.1, allow for flexible shrinkage parameter in penalty_pss()? 
      K1 <- penalty_pss(K = K1, difference = args$difference, shrink = 0.1)
    }    
  }else{ # Ridge-penalty
    K1 <- diag(ncol(Bs))
  }
  #K1 <- matrix(0, ncol=ncol(Bs), nrow=ncol(Bs))
  #print(args$penalty)

  if(args$inTime == "smooth"){
    K2 <- diff(diag(ncol(Bt)), differences = args$differences)
    K2 <- crossprod(K2)  
  }else{
    K2 <- diag(ncol(Bt))
  }

  # compute penalty matrix for the whole effect
  suppressMessages(K <- kronecker(K1, diag(ncol(Bt))) +
                      kronecker(diag(ncol(Bs)), K2))
  
  ## compare specified degrees of freedom to dimension of null space
  if (!is.null(args$df)){
    rns <- ncol(K) - qr(as.matrix(K))$rank # compute rank of null space
    if (rns == args$df)
      warning( sQuote("df"), " equal to rank of null space ",
               "(unpenalized part of P-spline);\n  ",
               "Consider larger value for ", sQuote("df"),
               ## " or set ", sQuote("center = TRUE"), 
               ".", immediate.=TRUE)
    if (rns > args$df)
      stop("not possible to specify ", sQuote("df"),
           " smaller than the rank of the null space\n  ",
           "(unpenalized part of P-spline). Use larger value for ",
           sQuote("df"), 
           ## " or set ", sQuote("center = TRUE"), 
           ".")
  }
  
  # save matrices to compute identifiability checks
  args$Bs <- Bs
  args$X1des <- X1des
  args$K1 <- K1
  args$L <- L
  
  # tidy up workspace 
  rm(Bs, Bt, ind0, X1des, X1, L)
  
  return(list(X = X, K = K, args = args))
}



### P-spline base learner for signal matrix with index vector
### for historical model according to function limit, defaults to s<=t
#' @rdname bsignal
#' @export
bhist <- function(x, s, time, index = NULL, #by = NULL, 
                  limits="s<=t", standard=c("no", "time", "length"), ##, "transform"
                  intFun=integrationWeightsLeft, 
                  inS=c("smooth","linear","constant"), inTime=c("smooth","linear","constant"),
                  knots = 10, boundary.knots = NULL, degree = 3, differences = 1, df = 4,
                  lambda = NULL, #center = FALSE, cyclic = FALSE
                  penalty = c("ps", "pss"), check.ident = FALSE
){
  
  if (!is.null(lambda)) df <- NULL
  
  cll <- match.call()
  cll[[1]] <- as.name("bhist")
  #print(cll)
  penalty <- match.arg(penalty)
  #print(penalty)
  
  standard <- match.arg(standard)
  
  inS <- match.arg(inS)
  inTime <- match.arg(inTime)
  #print(inS)
  
  if(! mboost_intern(x, fun = "isMATRIX") ) stop("signal has to be a matrix")
  if(ncol(x) != length(s)) stop("Dimension of x and s do not match.")
  
  varnames <- all.vars(cll)
  
  if(!is.atomic(s)) stop("index of signal has to be a vector")
  if(!is.atomic(time)) stop("index of response has to be a vector")
  
  if(substitute(time)==substitute(s)){
    #warning("Do not use the same variable t as time-variable in y(t) and in x(t).")
    warning("Do not use the same variable for s and time in bhist().")
  }
  
  # compare range of index signal and index response
  # minimal value of the signal-index has to be smaller than the response-index
  if(!is.function(limits)){
    if(limits=="s<=t" & min(s) > min(time) ) stop("Index of response has values before index of signal.")
  }
  
  # Reshape mfL so that it is the dataframe of the signal with 
  # the index of the signal and the index of the response as attributes
  xname <- varnames[1]
  indname <- varnames[2]
  indnameY <- varnames[3]
  if(length(varnames)==2) indnameY <- varnames[2]
  if(is.null(colnames(x))) colnames(x) <- paste(xname, 1:ncol(x), sep="_")
  attr(x, "signalIndex") <- s
  attr(x, "xname") <- xname
  attr(x, "indname") <- indname 
  attr(x, "indexY") <- time
  attr(x, "indnameY") <- indnameY
  attr(x, "id") <- index
  
  mf <- data.frame("z"=I(x))
  names(mf) <- xname 
  
  if(!all( abs(colMeans(x, na.rm = TRUE)) < .Machine$double.eps*10^10)){
    message(xname, " is not centered per column, inducing a non-centered effect.")
  }
  
  #   mf <- mfL
  #   names(mf) <- varnames
  
  vary <- ""
  
  # CC <- all(Complete.cases(mf))
  CC <- all(mboost_intern(mf, fun = "Complete.cases"))
  if (!CC)
    warning("base-learner contains missing values;\n",
            "missing values are excluded per base-learner, ",
            "i.e., base-learners may depend on different",
            " numbers of observations.")
  
  #index <- NULL 
  
  ## call X_hist in oder to compute parameter settings, e.g. 
  ## the transformation matrix Z, shrinkage penalty, identifiability problems...
  if(is.null(index)){
    ### X_hist for data in wide format with regular response
    temp <- X_hist(mf, vary, 
                      args = hyper_hist(mf, vary, knots = knots, boundary.knots = boundary.knots, 
                                        degree = degree, differences = differences,
                                        df = df, lambda = lambda, center = FALSE, cyclic = FALSE,
                                        s = s, time=time, limits = limits, 
                                        standard = standard, intFun = intFun, 
                                        inS = inS, inTime = inTime, 
                                        penalty = penalty, check.ident = check.ident, 
                                        format="wide"))
  }else{
    ### X_hist for data in long format with irregular response
    temp <- X_hist(mf, vary, 
                   args = hyper_hist(mf, vary, knots = knots, boundary.knots = boundary.knots, 
                                     degree = degree, differences = differences,
                                     df = df, lambda = lambda, center = FALSE, cyclic = FALSE,
                                     s = s, time=time, limits = limits, 
                                     standard = standard, intFun = intFun, 
                                     inS = inS, inTime = inTime, 
                                     penalty = penalty, check.ident = check.ident, 
                                     format="long"))
  }
  temp$args$check.ident <- FALSE
  
  ret <- list(model.frame = function() 
    if (is.null(index)) return(mf) else{
      mftemp <- mf
      mf <- mftemp[index,,drop = FALSE] # this is necessary to pass the attributes
      attributes(mftemp[,xname])
      attr(mf[,xname], "signalIndex") <- attr(mftemp[,xname], "signalIndex")
      attr(mf[,xname], "xname") <- attr(mftemp[,xname], "xname")
      attr(mf[,xname], "indname") <- attr(mftemp[,xname], "indname")
      attr(mf[,xname], "indexY") <- attr(mftemp[,xname], "indexY")
      attr(mf[,xname], "indnameY") <- attr(mftemp[,xname], "indnameY")
      attr(mf[,xname], "id") <- attr(mftemp[,xname], "id")
      return(mf)
    } ,
    get_call = function(){
      cll <- deparse(cll, width.cutoff=500L)
      if (length(cll) > 1)
        cll <- paste(cll, collapse="")
      cll
    },
    get_data = function() mf,
    ## set index to NULL, as the index is treated within X_hist()
    ##get_index = function() index, 
    get_index = function() NULL,
    get_vary = function() vary,
    get_names = function(){
      attr(xname, "indname") <- indname 
      attr(xname, "indnameY") <- indnameY
      xname 
    }, #colnames(mf),
    set_names = function(value) {
      #if(length(value) != length(colnames(mf)))
      if(length(value) != names(mf[1]))
        stop(sQuote("value"), " must have same length as ",
             sQuote("names(mf[1])"))
      for (i in 1:length(value)){
        cll[[i+1]] <<- as.name(value[i])
      }
      attr(mf, "names") <<- value
    })
  class(ret) <- "blg"
  
  ### X_hist is for data in wide format with regular response
  # ret$dpp <- bl_lin(ret, Xfun = X_hist, args = temp$args) 
  ret$dpp <- mboost_intern(ret, Xfun = X_hist, args = temp$args, fun = "bl_lin")
  
  #   ## function that comoutes a design matrix such that new_des %*% hat{theta} = beta(s)
  #   ## use ng equally spaced observation points 
  #   ret$get_des <- function(ng = 40){
  #     
  #     ## use a new grid of s-values with ng grid points 
  #     s_grid <- seq(min(s), max(s), l = ng)
  #     time_grid <- seq(min(time), max(time), l = ng)
  #     ## matrix with inverse integraion weights
  #     dummyX <- I( diag(length(s_grid)) /  intFun(diag(length(s_grid)), s_grid ) )
  #     
  #     ## setup for X_signal 
  #     attr(dummyX, "signalIndex") <- s_grid
  #     attr(dummyX, "xname") <- xname
  #     attr(dummyX, "indname") <- indname
  #     attr(x, "indexY") <- time_grid
  #     attr(x, "indnameY") <- indnameY
  #     attr(x, "id") <- index
  #     
  #     new_mf <- data.frame("z" = I(dummyX))
  #     names(new_mf) <- xname
  #     
  #     ## use vary and args like in bl  
  #     new_des <- X_hist(mf = new_mf, vary = vary, args = temp$args)$X
  #     
  #     ## give arguments to new_des for easier use in the plot-function 
  #     attr(new_des, "x") <- s_grid
  #     attr(new_des, "xlab") <- indname
  #     attr(new_des, "y") <- time_grid
  #     attr(new_des, "ylab") <- indnameY
  #     attr(new_des, "varname") <- xname
  #     
  #     return(new_des)
  #   }
  
  return(ret)
}

# testX <- I(matrix(rnorm(40), ncol=5))
# s <- seq(0,1,l=5)
# time <- s
# test <- bhist(testX, s, time, knots=5, df=5)
# test$get_names()
# test$get_data()
# names(test$dpp(rep(1,nrow(testX))))
# extract(test)[1:10, 1:20]


### hyper parameters for signal baselearner with eigenfunctions as bases, FPCA-based
hyper_fpc <- function(mf, vary, df = 4, lambda = NULL, 
                      pve = 0.99, npc = NULL, npc.max = 15, getEigen=TRUE, 
                      s=NULL, penalty = "identity") {
  ## prediction is usually set in/by newX() 
  list(df = df, lambda = lambda, pve = pve, npc = npc, npc.max = npc.max, 
       getEigen = getEigen, s = s, penalty = penalty, prediction = FALSE)
}

### model.matrix for FPCA based functional base-learner
X_fpc <- function(mf, vary, args) {

  stopifnot(is.data.frame(mf))
  xname <- names(mf)
  X1 <- as.matrix(mf)
  xind <- attr(mf[[1]], "signalIndex")
  if(is.null(xind)) xind <- args$s # if the attribute is NULL use the s of the model fit
  #print(xind)
  
  if(ncol(X1) != length(xind)) stop(xname, ": Dimension of signal matrix and its index do not match.")
  
  ## do FPCA on X1 (code of refund::ffpc adapted) using xind as argvals 
  if(is.null(args$klX)){
    
    decomppars <- list(argvals = xind, pve = args$pve, npc = args$npc, useSymm = TRUE)
    decomppars$Y <- X1
    ## functional covariate is per default centered per time-point
    klX <- do.call(refund::fpca.sc, decomppars)
    
    ## add the solution of the decomposition to args
    args$klX <- klX
    args$klX$xind <- xind
    
    ## only use part of the eigen-functions! 
    args$subset <- 1:min(ncol(klX$scores), args$npc.max)
    ## args$a <- max(xind) - min(xind)
    
    ## scores \xi_{ik}: rows i=1,..., N and columns k=1,...,K
    ## are the design matrix
    X <- klX$scores[ , args$subset, drop = FALSE]
    
    ## scores can be computed as \xi_{ik}=\int X1cen_i(s) \phi_k(s) ds 
    ## all(round(klX$scores,6) == round(scale(X1, center=klX$mu, scale=FALSE) %*% klX$efunctions, 6))
    
  }else{
    
    klX <- args$klX 
    ## compute scores on new X1 observations
    if(ncol(X1) == length(klX$mu) && all(args$s == xind)){
      ## is the same as "X <- klX$scores[ , args$subset, drop = FALSE]" if klX is FPCA on X1 
      X <- (scale(X1, center=klX$mu, scale=FALSE) %*% klX$efunctions)[ , args$subset, drop = FALSE]
      ## use integration weights?
      # X <- 1/args$a*(scale(X1, center=klX$mu, scale=FALSE) %*% klX$efunctions)[ , args$subset, drop = FALSE]
    }else{
      ##stop("In bfpc the grid for the functional covariate has to be the same as in the model fit!")
      ## linear interpolation of the basis functions 
      approxEfunctions <- matrix(NA, nrow=length(xind), ncol=length(args$subset))
      for(i in 1:ncol(klX$efunctions[ , args$subset, drop = FALSE])){
        approxEfunctions[,i] <- approx(x=args$klX$xind, y=klX$efunctions[,i], xout=xind)$y
      }
      approxMu <- approx(x=args$klX$xind, y=klX$mu, xout=xind)$y
      X <-(scale(X1, center=approxMu, scale=FALSE) %*% approxEfunctions)
      ## use integration weights?
      # X <- 1/args$a*(scale(X1, center=approxMu, scale=FALSE) %*% approxEfunctions)
    } 
    
  }

  colnames(X) <- paste(xname, ".PC", 1:ncol(X), sep = "")
  
  ## set up the penalty matrix 
  K <- switch(args$penalty, 
              ### use the identity matrix for penalization
              ### all eigenfunctions are penalized with the same strength
              identity = diag(rep(1, length = length(args$subset))),  
              ## Penalty matrix: diagonal matrix of inverse eigen-values
              ## implicit assumption: important eigen-functions of X process 
              ## are more important in shape of beta
              inverse = diag(1 / klX$evalues[args$subset]), 
              ### no penalty at all, as regularization is done by truncating the number of PCs used
              ### gives bad estimates
              no = matrix(0, ncol = length(args$subset), nrow = length(args$subset))
  )
  
  return(list(X = X, K = K, args = args))
}



###############################################################################
### FPCA based base-learner for signal matrix with index vector
### inspired by refund::fpca.sc
#' @rdname bsignal
#' @export
bfpc <- function(x, s, index = NULL, df = 4, 
                 lambda = NULL, penalty = c("identity", "inverse", "no"), 
                 pve = 0.99, npc = NULL, npc.max = 15, getEigen=TRUE
){
  
  if (!is.null(lambda)) df <- NULL
  
  cll <- match.call()
  cll[[1]] <- as.name("bfpc")
  penalty <- match.arg(penalty)
  
  if (!requireNamespace("refund", quietly = TRUE))
    stop("The package refund is needed for the function 'fpca.sc'.\nTo use the bfpc baseleaner, please install the package 'refund'.")
  if(! mboost_intern(x, fun = "isMATRIX") ) 
    stop("signal has to be a matrix")
  
  varnames <- all.vars(cll)
  
  # Reshape mfL so that it is the dataframe of the signal with the index as attribute
  xname <- varnames[1]
  indname <- varnames[2] 
  if(is.null(colnames(x))) colnames(x) <- paste(xname, 1:ncol(x), sep="_")
  attr(x, "signalIndex") <- s
  attr(x, "xname") <- xname
  attr(x, "indname") <- indname
  
  mf <- data.frame("z"=I(x))
  names(mf) <- xname
  
  vary <- ""
  
  ## improvement: for a FPCA based base-learner the X can contain missings!
  # CC <- all(Complete.cases(mf))
  CC <- all(mboost_intern(mf, fun = "Complete.cases"))
  if (!CC)
    warning("base-learner contains missing values;\n",
            "missing values are excluded per base-learner, ",
            "i.e., base-learners may depend on different",
            " numbers of observations.")
  
  #index <- NULL
  
  ## call X_fpc in oder to compute parameter settings, e.g. 
  ## the basis functions, based on FPCA 
  temp <- X_fpc(mf, vary, 
                args = hyper_fpc(mf, vary, df = df, lambda = lambda, 
                                 pve = pve, npc = npc, npc.max = npc.max, 
                                 s = s, penalty = penalty))
  ## save the FPCA in args
  ##str(temp$args)
  
  ret <- list(model.frame = function() 
    if (is.null(index)) return(mf) else{
      mftemp <- mf
      mf <- mftemp[index,,drop = FALSE] # this is necessary to pass the attributes
      attributes(mftemp[,xname])
      attr(mf[,xname], "signalIndex") <- attr(mftemp[,xname], "signalIndex")
      attr(mf[,xname], "xname") <- attr(mftemp[,xname], "xname")
      attr(mf[,xname], "indname") <- attr(mftemp[,xname], "indname")
      return(mf)
    } ,
    get_call = function(){
      cll <- deparse(cll, width.cutoff=500L)
      if (length(cll) > 1)
        cll <- paste(cll, collapse="")
      cll
    },
    get_data = function() mf,
    get_index = function() index,
    get_vary = function() vary,
    get_names = function(){
      attr(xname, "indname") <- indname 
      xname 
    }, #colnames(mf),
    set_names = function(value) {
      #if(length(value) != length(colnames(mf)))
      if(length(value) != names(mf[1]))
        stop(sQuote("value"), " must have same length as ",
             sQuote("names(mf[1])"))
      for (i in 1:length(value)){
        cll[[i+1]] <<- as.name(value[i])
      }
      attr(mf, "names") <<- value
    })
  class(ret) <- "blg"
  
  # ret$dpp <- bl_lin(ret, Xfun = X_fpc, args = temp$args)
  ret$dpp <- mboost_intern(ret, Xfun = X_fpc, args = temp$args, fun = "bl_lin")
  
  rm(temp)
  
  return(ret)
}






#######################################################################################
# Base-learner with constraints for smooth varying scalar covariate

### almost equal to X_bbs() of package mboost
### difference: implements sum-to-zero-constraint over index of response
X_bbsc <- function(mf, vary, args) {

  stopifnot(is.data.frame(mf))
  mm <- lapply(which(colnames(mf) != vary), function(i) {
    # X <- bsplines(mf[[i]],
    #              knots = args$knots[[i]]$knots,
    #              boundary.knots = args$knots[[i]]$boundary.knots,
    #              degree = args$degree, 
    #              Ts_constraint = args$Ts_constraint,
    #              deriv = args$deriv)
    X <- mboost_intern(mf[[i]],
                       knots = args$knots[[i]]$knots,
                       boundary.knots = args$knots[[i]]$boundary.knots,
                       degree = args$degree, 
                       Ts_constraint = args$Ts_constraint,
                       deriv = args$deriv, extrapolation = args$prediction, 
                       fun = "bsplines")
    if (args$cyclic) {
      # X <- cbs(mf[[i]],
      #         knots = args$knots[[i]]$knots,
      #         boundary.knots = args$knots[[i]]$boundary.knots,
      #         degree = args$degree,
      #         deriv = args$deriv)
      X <- mboost_intern(mf[[i]],
                         knots = args$knots[[i]]$knots,
                         boundary.knots = args$knots[[i]]$boundary.knots,
                         degree = args$degree,
                         deriv = args$deriv,
                         fun = "cbs")
    }
    class(X) <- "matrix"
    return(X)
  })  ### options
  MATRIX <- any(sapply(mm, dim) > c(500, 50)) || (length(mm) > 1)
  MATRIX <- MATRIX && options("mboost_useMatrix")$mboost_useMatrix
  if (MATRIX) {
    diag <- Diagonal
    for (i in 1:length(mm)){
      tmp <- attributes(mm[[i]])[c("degree", "knots", "Boundary.knots")]
      mm[[i]] <- Matrix(mm[[i]])
      attributes(mm[[i]])[c("degree", "knots", "Boundary.knots")] <- tmp
    }
  }
  
  if (length(mm) == 1) {
    X <- mm[[1]]
    if (vary != "") {
      by <- model.matrix(as.formula(paste("~", vary, collapse = "")),
                         data = mf)[ , -1, drop = FALSE] # drop intercept
      DM <- lapply(1:ncol(by), function(i) {
        ret <- X * by[, i]
        colnames(ret) <- paste(colnames(ret), colnames(by)[i], sep = ":")
        ret
      })
      X <- do.call("cbind", DM)
    }
    if (args$differences > 0){
      if (!args$cyclic) {
        K <- diff(diag(ncol(mm[[1]])), differences = args$differences)
      } else {
        ## cyclic P-splines
        differences <- args$differences
        K <- diff(diag(ncol(mm[[1]]) + differences),
                  differences = differences)
        tmp <- K[,(1:differences)]   # save first "differences" columns
        K <- K[,-(1:differences)]    # drop first "differences" columns
        indx <- (ncol(mm[[1]]) - differences + 1):(ncol(mm[[1]]))
        K[,indx] <- K[,indx] + tmp   # add first "differences" columns
      }
    } else {
      if (args$differences != 0)
        stop(sQuote("differences"), " must be an non-neative integer")
      K <- diag(ncol(mm[[1]]))
    }
    
    if (vary != "" && ncol(by) > 1){       # build block diagonal penalty
      suppressMessages(K <- kronecker(diag(ncol(by)), K))
    }
    if (args$center) {
      tmp <- attributes(X)[c("degree", "knots", "Boundary.knots")]
      center <- match.arg(as.character(args$center),
                          choices = c("TRUE", "differenceMatrix", "spectralDecomp"))
      if (center == "TRUE") center <- "differenceMatrix"
      X <- switch(center,
                  ### L = t(D) in Section 2.3. of Fahrmeir et al. (2004, Stat Sinica)
                  "differenceMatrix" = tcrossprod(X, K) %*% solve(tcrossprod(K)),
                  ### L = \Gamma \Omega^1/2 in Section 2.3. of
                  ### Fahrmeir et al. (2004, Stat Sinica)
                  "spectralDecomp" = {
                    SVD <- eigen(crossprod(K), symmetric = TRUE)
                    ev <- SVD$vector[, 1:(ncol(X) - args$differences), drop = FALSE]
                    ew <- SVD$values[1:(ncol(X) - args$differences), drop = FALSE]
                    X %*% ev %*% diag(1/sqrt(ew))
                  }
      )
      attributes(X)[c("degree", "knots", "Boundary.knots")] <- tmp
      K <- diag(ncol(X)) 
    } else {
      K <- crossprod(K)
    }
    if (!is.null(attr(X, "Ts_constraint"))) {
      D <- attr(X, "D")
      K <- crossprod(D, K) %*% D
    }
  }
  
  if (length(mm) == 2) {
    suppressMessages(
      X <- kronecker(mm[[1]], matrix(1, ncol = ncol(mm[[2]]))) *
        kronecker(matrix(1, ncol = ncol(mm[[1]])), mm[[2]])
    )
    if (vary != "") {
      by <- model.matrix(as.formula(paste("~", vary, collapse = "")),
                         data = mf)[ , -1, drop = FALSE] # drop intercept
      DM <- X * by[,1]
      if (ncol(by) > 1){
        for (i in 2:ncol(by))
          DM <- cbind(DM, (X * by[,i]))
      }
      X <- DM
      ### <FIXME> Names of X if by is given
    }
    if (args$differences > 0){
      if (!args$cyclic) {
        Kx <- diff(diag(ncol(mm[[1]])), differences = args$differences)
        Ky <- diff(diag(ncol(mm[[2]])), differences = args$differences)
      } else {
        ## cyclic P-splines
        differences <- args$differences
        Kx <- diff(diag(ncol(mm[[1]]) + differences),
                   differences = differences)
        Ky <- diff(diag(ncol(mm[[2]]) + differences),
                   differences = differences)
        
        tmp <- Kx[,(1:differences)]   # save first "differences" columns
        Kx <- Kx[,-(1:differences)]    # drop first "differences" columns
        indx <- (ncol(mm[[1]]) - differences + 1):(ncol(mm[[1]]))
        Kx[,indx] <- Kx[,indx] + tmp   # add first "differences" columns
        
        tmp <- Ky[,(1:differences)]   # save first "differences" columns
        Ky <- Ky[,-(1:differences)]    # drop first "differences" columns
        indx <- (ncol(mm[[2]]) - differences + 1):(ncol(mm[[2]]))
        Ky[,indx] <- Ky[,indx] + tmp   # add first "differences" columns
      }
    } else {
      if (args$differences != 0)
        stop(sQuote("differences"), " must be an non-negative integer")
      Kx <- diag(ncol(mm[[1]]))
      Ky <- diag(ncol(mm[[2]]))
    }
    
    Kx <- crossprod(Kx)
    Ky <- crossprod(Ky)
    suppressMessages(
      K <- kronecker(Kx, diag(ncol(mm[[2]]))) +
        kronecker(diag(ncol(mm[[1]])), Ky)
    )
    if (vary != "" && ncol(by) > 1){       # build block diagonal penalty
      suppressMessages(K <- kronecker(diag(ncol(by)), K))
    }
    if (!identical(args$center, FALSE)) {
      ### L = \Gamma \Omega^1/2 in Section 2.3. of Fahrmeir et al.
      ### (2004, Stat Sinica), always
      L <- eigen(K, symmetric = TRUE)
      L$vectors <- L$vectors[,1:(ncol(X) - args$differences^2), drop = FALSE]
      L$values <- sqrt(L$values[1:(ncol(X) - args$differences^2), drop = FALSE])
      L <- L$vectors %*% (diag(length(L$values)) * (1/L$values))
      X <- as(X %*% L, "matrix")
      K <- as(diag(ncol(X)), "matrix")
    }
  }
  
  if (length(mm) > 2)
    stop("not possible to specify more than two variables in ",
         sQuote("..."), " argument of smooth base-learners")
  
  #----------------------------------
  ### <SB> Calculate constraints

  ## for center = TRUE, design matrix does not contain constant part 
  if(args$center != FALSE){ 
    
    ## center the columns of the design matrix 
    ## Z contains column means
    # If the argument Z is not NULL use the given Z (important for prediction!)
    if(is.null(args$Z)){
      args$Z <- colMeans(X)
    }
    
    ### Transform design and penalty matrix 
    ## use column means of original design matrix
    X <- scale(X, center = args$Z, scale = FALSE)
    
  }else{
    ## sum-to-zero constraint - orthogonal to constant part, 
    ## cf. Web Appendix A of Brockhaus et al. 2015
    
    # If the argument Z is not NULL use the given Z (important for prediction!)
    if(is.null(args$Z)){
      C <- t(X) %*% rep(1, nrow(X))
      Q <- qr.Q(qr(C), complete=TRUE) # orthonormal matrix of QR decomposition
      args$Z <- Q[  , 2:ncol(Q)] # only keep last columns    
    }
    
    ### Transform design and penalty matrix
    X <- X %*% args$Z
    K <- t(args$Z) %*% K %*% args$Z
    #print(args$Z)
  }

  #----------------------------------
  
  ## compare specified degrees of freedom to dimension of null space
  if (!is.null(args$df)){
    rns <- ncol(K) - qr(as.matrix(K))$rank # compute rank of null space
    if (rns == args$df)
      warning( sQuote("df"), " equal to rank of null space ",
               "(unpenalized part of P-spline);\n  ",
               "Consider larger value for ", sQuote("df"),
               " or set ", sQuote("center != FALSE"), ".", immediate.=TRUE)
    if (rns > args$df)
      stop("not possible to specify ", sQuote("df"),
           " smaller than the rank of the null space\n  ",
           "(unpenalized part of P-spline). Use larger value for ",
           sQuote("df"), " or set ", sQuote("center != FALSE"), ".")
  }
  
  return(list(X = X, K = K, args = args))
}

## add the parameter Z to the arguments of hyper_bbs()
hyper_bbsc <- function(Z, ...){
  ret <- c(mboost_intern(..., fun = "hyper_bbs"), list(Z=Z))
  return(ret)
}


###############################################################################

#' Constrained Base-learners for Scalar Covariates
#' 
#' Constrained base-learners for fitting effects of scalar covariates in models 
#' with functional response
#' 
#' @param ... one or more predictor variables or one matrix or data 
#' frame of predictor variables. 
#' @param by an optional variable defining varying coefficients, 
#' either a factor or numeric variable.
#' @param index a vector of integers for expanding the variables in \code{...}.
#' @param knots either the number of knots or a vector of the positions 
#' of the interior knots (for more details see \code{\link[mboost]{bbs}}).
#' @param boundary.knots boundary points at which to anchor the B-spline basis 
#' (default the range of the data). A vector (of length 2) 
#' for the lower and the upper boundary knot can be specified.
#' @param degree degree of the regression spline.
#' @param differences a non-negative integer, typically 1, 2 or 3. 
#' If \code{differences} = \emph{k}, \emph{k}-th-order differences are used as 
#' a penalty (\emph{0}-th order differences specify a ridge penalty).
#' @param df trace of the hat matrix for the base-learner defining the 
#' base-learner complexity. Low values of \code{df} correspond to a 
#' large amount of smoothing and thus to "weaker" base-learners.
#' @param lambda smoothing parameter of the penalty, computed from \code{df} when 
#' \code{df} is specified. 
#' @param K in \code{bolsc} it is possible to specify the penalty matrix K
#' @param weights experiemtnal! weights that are used for the computation of the transformation matrix Z.
#' @param center See \code{\link[mboost]{bbs}}. 
#' @param cyclic  if \code{cyclic = TRUE} the fitted values coincide at 
#' the boundaries (useful for cyclic covariates such as day time etc.).
#' @param contrasts.arg Note that a special \code{contrasts.arg} exists in 
#' package \code{mboost}, namely "contr.dummy". This contrast is used per default 
#' in \code{brandomc}. It leads to a 
#' dummy coding as returned by \code{model.matrix(~ x - 1)} were the 
#' intercept is implicitly included but each factor level gets a 
#' separate effect estimate (for more details see \code{\link[mboost]{brandom}}).
#' @param intercept if \code{intercept = TRUE} an intercept is added to the design matrix 
#' of a linear base-learner. 
#' 
#' @details The base-learners \code{bbsc}, \code{bolsc} and \code{brandomc} are 
#' the base-learners \code{\link[mboost]{bbs}}, \code{\link[mboost]{bols}} and 
#' \code{\link[mboost]{brandom}} with additional identifiability constraints. 
#' The constraints enforce that 
#' \eqn{\sum_{i} \hat h(x_i, t) = 0} for all \eqn{t}, so that 
#' effects varying over \eqn{t} can be interpreted as deviations 
#' from the global functional intercept, see Web Appendix A of 
#' Scheipl et al. (2015). 
#' The constraint is enforced by a basis transformation of the design and penalty matrix. 
#' In particular, it is sufficient to apply the constraint on the covariate-part of the design 
#' and penalty matrix and thus, it is not necessary to change the basis in $t$-direction.  
#' See Appendix A of Brockhaus et al. (2015) for technical details on how to enforce this sum-to-zero constraint.   
#' 
#' Cannot deal with any missing values in the covariates.
#' 
#' @return Equally to the base-learners of package \code{mboost}: 
#' 
#' An object of class \code{blg} (base-learner generator) with a 
#' \code{dpp} function (data pre-processing) and other functions. 
#' 
#' The call to \code{dpp} returns an object of class 
#' \code{bl} (base-learner) with a \code{fit} function. The call to 
#' \code{fit} finally returns an object of class \code{bm} (base-model).
#' 
#' @seealso \code{\link{FDboost}} for the model fit. 
#' \code{\link[mboost]{bbs}}, \code{\link[mboost]{bols}} and \code{\link[mboost]{brandom}} for the 
#' corresponding base-learners in \code{mboost}. 
#' 
#' @references 
#' Brockhaus, S., Scheipl, F., Hothorn, T. and Greven, S. (2015): 
#' The functional linear array model. Statistical Modelling, 15(3), 279-300.
#' 
#' Scheipl, F., Staicu, A.-M. and Greven, S. (2015):  
#' Functional Additive Mixed Models, Journal of Computational and Graphical Statistics, 24(2), 477-501.
#' 
#' @author Sarah Brockhaus, Almond Stoecker 
#' 
#' @examples 
#' #### simulate data with functional response and scalar covariate (functional ANOVA)
#' n <- 60   ## number of cases
#' Gy <- 27  ## number of observation poionts per response curve 
#' dat <- list()
#' dat$t <- (1:Gy-1)^2/(Gy-1)^2
#' set.seed(123)
#' dat$z1 <- rep(c(-1, 1), length = n)
#' dat$z1_fac <- factor(dat$z1, levels = c(-1, 1), labels = c("1", "2"))
#' # dat$z1 <- runif(n)
#' # dat$z1 <- dat$z1 - mean(dat$z1)
#' 
#' # mean and standard deviation for the functional response 
#' mut <- matrix(2*sin(pi*dat$t), ncol = Gy, nrow = n, byrow = TRUE) + 
#'         outer(dat$z1, dat$t, function(z1, t) z1*cos(pi*t) ) # true linear predictor
#' sigma <- 0.1
#' 
#' # draw respone y_i(t) ~ N(mu_i(t), sigma)
#' dat$y <- apply(mut, 2, function(x) rnorm(mean = x, sd = sigma, n = n)) 
#' 
#' ## fit function-on-scalar model with a linear effect of z1
#' m1 <- FDboost(y ~ 1 + bolsc(z1_fac, df = 1), timeformula = ~ bbs(t, df = 6), data = dat)
#' 
#' # look for optimal mSTOP using cvrisk() or validateFDboost()
#'  \dontrun{
#' cvm <- cvrisk(m1, grid = 1:500)
#' m1[mstop(cvm)]
#' }
#' m1[200] # use 200 boosting iterations 
#' 
#' # plot true and estimated coefficients 
#' plot(dat$t, 2*sin(pi*dat$t), col = 2, type = "l", main = "intercept")
#' plot(m1, which = 1, lty = 2, add = TRUE)
#' 
#' plot(dat$t, 1*cos(pi*dat$t), col = 2, type = "l", main = "effect of z1")
#' lines(dat$t, -1*cos(pi*dat$t), col = 2, type = "l")
#' plot(m1, which = 2, lty = 2, col = 1, add = TRUE)
#' 
#' 
#' @keywords models
#' @aliases brandomc bolsc
#' @export
bbsc <- function(..., by = NULL, index = NULL, knots = 10, boundary.knots = NULL,
                 degree = 3, differences = 2, df = 4, lambda = NULL, center = FALSE,
                 cyclic = FALSE) {
  
  #----------------------------------
  ## <SB> arguments constraint and dervi of bbs() are set to their defaults
  ## i.e. no constraints and no derivatives 
  # constraint <- match.arg(constraint)
  constraint <- "none"
  deriv <- 0
  #----------------------------------
  
  if (!is.null(lambda)) df <- NULL
  
  cll <- match.call()
  cll[[1]] <- as.name("bbsc")
  
  mf <- list(...)
  if (length(mf) == 1 && ((is.matrix(mf[[1]]) || is.data.frame(mf[[1]])) &&
                            ncol(mf[[1]]) > 1 )) {
    mf <- as.data.frame(mf[[1]])
  } else {
    mf <- as.data.frame(mf)
    cl <- as.list(match.call(expand.dots = FALSE))[2][[1]]
    colnames(mf) <- sapply(cl, function(x) deparse(x))
  }
  stopifnot(is.data.frame(mf))
  if(!(all(sapply(mf, is.numeric)))) {
    if (ncol(mf) == 1){
      warning("cannot compute ", sQuote("bbsc"),
              " for non-numeric variables; used ",
              sQuote("bols"), " instead.")
      return(bols(mf, by = by, index = index))
    }
    stop("cannot compute bbsc for non-numeric variables")
  }
  vary <- ""
  if (!is.null(by)){
    mf <- cbind(mf, by)
    colnames(mf)[ncol(mf)] <- vary <- deparse(substitute(by))
  }
  
  # CC <- all(Complete.cases(mf))
  CC <- all(mboost_intern(mf, fun = "Complete.cases"))
  if (!CC)
    warning("base-learner contains missing values;\n",
            "missing values are excluded per base-learner, ",
            "i.e., base-learners may depend on different",
            " numbers of observations.")
  
  ################# do not use the index option as then Z is always computed as for balanced data
  ################# make an option available for this? as this means centering per group
  ### option
  DOINDEX <- (nrow(mf) > options("mboost_indexmin")[[1]])
  if (is.null(index)) {
    if (!CC || DOINDEX) {
      # index <- get_index(mf)
      index <- mboost_intern(mf, fun = "get_index")
      mf <- mf[index[[1]],,drop = FALSE]
      index <- index[[2]]
    }
  }
  
  ## call X_bbsc in oder to compute the transformation matrix Z
  if(is.null(index)){
    temp <- X_bbsc(mf, vary, 
                   args = hyper_bbsc(mf, vary, knots = knots, boundary.knots =
                                       boundary.knots, degree = degree, differences = differences,
                                     df = df, lambda = lambda, center = center, cyclic = cyclic, 
                                     constraint = constraint, deriv = deriv, 
                                     Z = NULL))
  }else{
    temp <- X_bbsc(mf, vary, 
                   args = hyper_bbsc(mf[index,,drop = FALSE], vary, knots = knots, boundary.knots =
                                       boundary.knots, degree = degree, differences = differences,
                                     df = df, lambda = lambda, center = center, cyclic = cyclic, 
                                     constraint = constraint, deriv = deriv, 
                                     Z = NULL))
  }

  Z <- temp$args$Z
  
  ret <- list(model.frame = function()
    if (is.null(index)) return(mf) else return(mf[index,,drop = FALSE]),
    get_call = function(){
      cll <- deparse(cll, width.cutoff=500L)
      if (length(cll) > 1)
        cll <- paste(cll, collapse="")
      cll
    },
    get_data = function() mf,
    get_index = function() index,
    get_vary = function() vary,
    get_names = function() colnames(mf),
    set_names = function(value) {
      if(length(value) != length(colnames(mf)))
        stop(sQuote("value"), " must have same length as ",
             sQuote("colnames(mf)"))
      for (i in 1:length(value)){
        cll[[i+1]] <<- as.name(value[i])
      }
      attr(mf, "names") <<- value
    })
  class(ret) <- "blg"
  
  # ret$dpp <- bl_lin(ret, Xfun = X_bbsc, args = temp$args)
  ret$dpp <- mboost_intern(ret, Xfun = X_bbsc, args = temp$args, fun = "bl_lin")
  
  #   ## function that comoutes a design matrix such that new_des %*% hat{theta} = beta(s)
  #   ## use ng equally spaced observation points 
  #   ret$get_des <- function(ng = 40){
  #     
  #     if(ncol(mf) > 1) stop("get_des() in bbsc() is only implemented for one variable.")
  #     
  #     ## use a new grid of x-values with ng grid points 
  #     x_grid <- seq(min(mf[ , 1]), max(mf[ , 1]), l = ng)
  # 
  #     new_mf <- data.frame("z" = I(x_grid))
  #     names(new_mf) <- names(mf)[1]
  #     
  #     ## use vary and args like in bl  
  #     new_des <- X_bbsc(mf = new_mf, vary = vary, args = temp$args)$X
  #     
  #     ## give arguments to new_des for easier use in the plot-function 
  #     attr(new_des, "x") <- new_mf[ , 1]
  #     attr(new_des, "xlab") <- names(new_mf)[1]
  # 
  #     return(new_des)
  #   }
  
  return(ret)
}


# z2 <- rnorm(17)
# blz <- bbsc(z=z2, df=3)
# blz$get_call()
# blz$get_names()
# str(blz$get_data())
# str(blz$dpp(weights=rep(1,17)))

#################
### model.matrix for constrained ols base-learner with penalty matrix K
X_olsc <- function(mf, vary, args) {
  
  if ( mboost_intern(mf, fun = "isMATRIX") ) {
    X <- mf
    contr <- NULL
  } else {
    ### set up model matrix
    fm <- paste("~ ", paste(colnames(mf)[colnames(mf) != vary],
                            collapse = "+"), sep = "")
    fac <- sapply(mf[colnames(mf) != vary], is.factor)
    if (any(fac)){
      if (!is.list(args$contrasts.arg)){
        ## first part needed to prevent warnings from calls such as
        ## contrasts.arg = contr.treatment(4, base = 1):
        if (is.character(args$contrasts.arg) &&
              args$contrasts.arg == "contr.dummy"){
          if (!args$intercept)
            stop('"contr.dummy" can only be used with ',
                 sQuote("intercept = TRUE"))
          fm <- paste(fm, "-1")
        } else {
          txt <- paste("list(", paste(colnames(mf)[colnames(mf) != vary][fac],
                                      "= args$contrasts.arg", collapse = ", "),")")
          args$contrasts.arg <- eval(parse(text=txt))
        }
      } else {
        ## if contrasts are given as list check if "contr.dummy" is specified
        if (any(args$contrasts.arg == "contr.dummy"))
          stop('"contr.dummy"',
               " can only be used for all factors at the same time.\n",
               "Use ", sQuote('contrasts.arg = "contr.dummy"'),
               " to achieve this.")
      }
    } else {
      args$contrasts.arg <- NULL
    }
    X <- model.matrix(as.formula(fm), data = mf, contrasts.arg = args$contrasts.arg)
    if (!is.null(args$contrasts.arg) && args$contrasts.arg == "contr.dummy")
      attr(X, "contrasts") <- lapply(attr(X, "contrasts"),
                                     function(x) x <- "contr.dummy")
    contr <- attr(X, "contrasts")
    if (!args$intercept)
      X <- X[ , -1, drop = FALSE]
    MATRIX <- any(dim(X) > c(500, 50)) && any(fac)
    MATRIX <- MATRIX && options("mboost_useMatrix")$mboost_useMatrix
    if (MATRIX) {
      diag <- Diagonal
      if (!is(X, "Matrix"))
        X <- Matrix(X)
    }
    if (vary != "") {
      by <- model.matrix(as.formula(paste("~", vary, collapse = "")),
                         data = mf)[ , -1, drop = FALSE] # drop intercept
      DM <- lapply(1:ncol(by), function(i) {
        ret <- X * by[, i]
        colnames(ret) <- paste(colnames(ret), colnames(by)[i], sep = ":")
        ret
      })
      X <- do.call("cbind", DM)
    }
  }
  
  #----------------------------------
  ## <SB> use given penalty-matrix K
  if(is.null(args$K)){
    ### <FIXME> penalize intercepts???
    ### set up penalty matrix
    ANOVA <- (!is.null(contr) && (length(contr) == 1)) && (ncol(mf) == 1)
    K <- diag(ncol(X))
    ### for ordered factors use difference penalty
    if (ANOVA && any(sapply(mf[, names(contr), drop = FALSE], is.ordered))) {
      K <- diff(diag(ncol(X) + 1), differences = 1)[, -1, drop = FALSE]
      if (vary != "" && ncol(by) > 1){       # build block diagonal penalty
        suppressMessages(K <- kronecker(diag(ncol(by)), K))
      }
      K <- crossprod(K)
    }
  }else{
    ## <SB> check dimensions of given K
    stopifnot(dim(args$K)==rep(ncol(X), 2))
    K <- args$K
  }
  #----------------------------------
  
  #----------------------------------
  ### <SB> Calculate constraints

  # Only compute Z in model fit, not in model prediction  
  #if(!args$prediction){ ## computes Z 3 times during model fit
  if(is.null(args$Z)){
    C <- t(X) %*% rep(1, nrow(X))
    Q <- qr.Q(qr(C), complete=TRUE) # orthonormal matrix of QR decomposition
    args$Z <- Q[  , 2:ncol(Q)] # only keep last columns
  }
  
  ### Transform design and penalty matrix
  X <- X %*% args$Z
  K <- t(args$Z) %*% K %*% args$Z
  
  #----------------------------------
  
  if (is(X, "Matrix") && !is(K, "Matrix"))
    K <- Matrix(K)
  
  ### <SB> return the transformation matrix Z as well
  list(X = X, K = K, args = args)
}


#' @rdname bbsc
#' @export
### Linear base-learner, potentially Ridge-penalized (but not by default)
### one can specify the penalty matrix K
### with sum-to-zero constraint over index of response
bolsc <- function(..., by = NULL, index = NULL, intercept = TRUE, df = NULL,
                  lambda = 0, K = NULL, weights = NULL, contrasts.arg = "contr.treatment") {
  
  if (!is.null(df)) lambda <- NULL
  
  cll <- match.call()
  cll[[1]] <- as.name("bolsc")
  
  mf <- list(...)

  if(!intercept && length(mf)==1) stop("Intercept has to be TRUE for bolsc with one covariate.")
  
  if (length(mf) == 1 && (( mboost_intern(mf[[1]], fun = "isMATRIX") || 
                            is.data.frame(mf[[1]])) &&
                            ncol(mf[[1]]) > 1 )) {
    mf <- mf[[1]]
    ### spline bases should be matrices
    if ( mboost_intern(mf, fun = "isMATRIX") && !is(mf, "Matrix"))
      class(mf) <- "matrix"
  } else {
    mf <- as.data.frame(mf)
    cl <- as.list(match.call(expand.dots = FALSE))[2][[1]]
    colnames(mf) <- sapply(cl, function(x) as.character(x))
  }
  if(!intercept && !any(sapply(mf, is.factor)) &&
       !any(sapply(mf, function(x){uni <- unique(x);
                                   length(uni[!is.na(uni)])}) == 1)){
    ## if no intercept is used and no covariate is a factor
    ## and if no intercept is specified (i.e. mf[[i]] is constant)
    if (any(sapply(mf, function(x) abs(mean(x, na.rm=TRUE) / sd(x,na.rm=TRUE))) > 0.1))
      ## if covariate mean is not near zero
      warning("covariates should be (mean-) centered if ",
              sQuote("intercept = FALSE"))
  }
  vary <- ""
  if (!is.null(by)){
    stopifnot(is.data.frame(mf))
    mf <- cbind(mf, by)
    colnames(mf)[ncol(mf)] <- vary <- deparse(substitute(by))
  }
  
  # CC <- all(Complete.cases(mf))
  CC <- all(mboost_intern(mf, fun = "Complete.cases"))
  if (!CC)
    warning("base-learner contains missing values;\n",
            "missing values are excluded per base-learner, ",
            "i.e., base-learners may depend on different",
            " numbers of observations.")

  ### option
  DOINDEX <- is.data.frame(mf) &&
    (nrow(mf) > options("mboost_indexmin")[[1]] || is.factor(mf[[1]]))
  if (is.null(index)) {
    ### try to remove duplicated observations or
    ### observations with missings
    if (!CC || DOINDEX) {
      # index <- get_index(mf)
      index <- mboost_intern(mf, fun = "get_index")
      mf <- mf[index[[1]],,drop = FALSE]
      index <- index[[2]]
    }
  }

  ### call X_olsc in order to compute the transformation matrix Z, 
  ## Z is saved in args$Z and is used after the model fit.  
  ### use index, as otherwise Z is computed as for one observation per factor level/ per unique observation
  ## this is equivalent to the same number of observations in each factor level 
  ### use weights, as otherwise the weights are not used for the computation of Z, 
  ## but weights here are NOT the weights in model call as they are an argument to bolsc()
  if(is.null(index)){
    
    if(is.null(weights)){ ## use weights 
      w <- 1:nrow(mf)
    }else{
      w <- rep(1:nrow(mf), weights)
    }
    
    temp <- X_olsc(mf[w, , drop = FALSE], vary, 
                   args = hyper_olsc(df = df, lambda = lambda, 
                                     intercept = intercept, contrasts.arg = contrasts.arg,
                                     K = K, Z = NULL)) 
  }else{
    
    if(is.null(weights)){  ## use weights
      w <- 1:nrow(mf[index, , drop = FALSE])
    }else{
      w <- rep(1:nrow(mf[index, , drop = FALSE]), weights)
    }
    
    temp <- X_olsc(mf = (mf[index, , drop = FALSE])[w, , drop = FALSE], vary = vary, 
                   args = hyper_olsc(df = df, lambda = lambda, 
                                     intercept = intercept, contrasts.arg = contrasts.arg,
                                     K = K, Z = NULL))  
  }

  
  ret <- list(model.frame = function()
    if (is.null(index)) return(mf) else return(mf[index,,drop = FALSE]),
    get_call = function(){
      cll <- deparse(cll, width.cutoff=500L)
      if (length(cll) > 1)
        cll <- paste(cll, collapse="")
      cll
    },
    get_data = function() mf,
    get_index = function() index,
    get_names = function() colnames(mf),
    get_vary = function() vary,
    set_names = function(value) {
      if(length(value) != length(colnames(mf)))
        stop(sQuote("value"), " must have same length as ",
             sQuote("colnames(mf)"))
      for (i in 1:length(value)){
        cll[[i+1]] <<- as.name(value[i])
      }
      attr(mf, "names") <<- value
    })
  class(ret) <- "blg"
  
  # ret$dpp <- bl_lin(ret, Xfun = X_olsc, args = temp$args)
  ret$dpp <- mboost_intern(ret, Xfun = X_olsc, args = temp$args, fun = "bl_lin")
  return(ret)
}


### hyper parameters for olsc base-learner
## add the parameter Z and K to the arguments of hyper_ols()
hyper_olsc <- function(Z = NULL, K = NULL, ...){
  
  ## prediction is usually set in/by newX()
  ret <- c(mboost_intern(..., fun = "hyper_ols"), list(Z = Z, K = K, prediction = FALSE))
  
  return(ret)
}


#' @rdname bbsc
#' @export
# random-effects (Ridge-penalized ANOVA) base-learner
# almost equal to brandom, but with sum-to-zero-constraint over index of t
brandomc <- function (..., contrasts.arg = "contr.dummy", df = 4) {
  cl <- cltmp <- match.call()
  if (is.null(cl$df))
    cl$df <- df
  if (is.null(cl$contrasts.arg))
    cl$contrasts.arg <- contrasts.arg
  cl[[1L]] <- as.name("bolsc")
  ret <- eval(cl, parent.frame())
  cltmp[[1]] <- as.name("brandomc")
  assign("cll", cltmp, envir = environment(ret$get_call))
  ret
}


