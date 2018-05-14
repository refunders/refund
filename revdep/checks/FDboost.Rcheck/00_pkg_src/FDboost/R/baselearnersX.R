#### Base-learner for historic effect of functional covariate
### with integral over specific limits, e.g. s<=t

### hyper parameters for signal baselearner with P-splines
hyper_histx <- function(mf, vary, knots = 10, boundary.knots = NULL, degree = 3,
                       differences = 1, df = 4, lambda = NULL, center = FALSE,
                       cyclic = FALSE, constraint = "none", deriv = 0L, 
                       Z=NULL, limits=NULL, 
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
  
  nm <- colnames(mf)[colnames(mf) != vary]
  if (is.list(knots)) if(!all(names(knots) %in% nm))
    stop("variable names and knot names must be the same")
  if (is.list(boundary.knots)) if(!all(names(boundary.knots) %in% nm))
    stop("variable names and boundary.knot names must be the same")
  if (!identical(center, FALSE) && cyclic)
    stop("centering of cyclic covariates not yet implemented")
  ret <- vector(mode = "list", length = length(nm))
  names(ret) <- nm
  for (n in nm)
    ret[[n]] <- knotf(getTime(mf[[n]]), 
                      knots=if(is.list(knots)) knots[[n]] else knots,
                      boundary.knots = if(is.list(boundary.knots)) boundary.knots[[n]] else boundary.knots)
  if (cyclic & constraint != "none")
    stop("constraints not implemented for cyclic B-splines")
  stopifnot(is.numeric(deriv) & length(deriv) == 1)
  
  ## prediction is usually set in/by newX()
  list(knots = ret, degree = degree, differences = differences,
       df = df, lambda = lambda, center = center, cyclic = cyclic,
       Ts_constraint = constraint, deriv = deriv, 
       Z = Z, limits = limits, 
       standard = standard, intFun = intFun, 
       inS = inS, inTime = inTime, 
       penalty = penalty, check.ident = check.ident, format = format, 
       prediction = FALSE)
}


### model.matrix for P-splines base-learner of signal matrix mf
### for response observed over a common grid, args$format="wide"
### or irregularly observed reponse, args$format="long" 
X_histx <- function(mf, vary, args) {
  
  stopifnot(is.data.frame(mf))
  varnames <- names(mf)
  
  xname <- getXLab(mf[[1]])
  X1 <- getX(mf[[1]])
  xind <- getArgvals(mf[[1]])
  yind <- getTime(mf[[1]])
  
  # force user to supply indices
  #if(is.null(xind)) xind <- args$s # if the attribute is NULL use the s of the model fit
  #if(is.null(yind)) yind <- args$time # if the attribute is NULL use the time of the model fit  
  
  nobs <- nrow(X1)
  
  ## get id-variable
  id <- getId(mf[[1]])  ### id is always there, in long and wide format!!!
  ## id has values 1, 2, 3, ... for response in long format
  
  # compute design-matrix in s-direction
  Bs <- switch(args$inS, 
               # B-spline basis of specified degree 
               #"smooth" = bsplines(xind, knots=args$knots[[varnames]]$knots, 
               #                     boundary.knots=args$knots[[varnames]]$boundary.knots, 
               #                    degree=args$degree),
               "smooth" = mboost_intern(xind, knots = args$knots[[varnames]]$knots, 
                                        boundary.knots = args$knots[[varnames]]$boundary.knots, 
                                        degree = args$degree, 
                                        fun = "bsplines"),
               "linear" = matrix(c(rep(1, length(xind)), xind), ncol = 2),
               "constant"=  matrix(c(rep(1, length(xind))), ncol = 1))
  
  colnames(Bs) <- paste(xname, 1:ncol(Bs), sep="")
  
  # integration weights 
  L <- args$intFun(X1=X1, xind=xind)
  
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
  
  # use function limits to set up design matrix according to function limits 
  # by setting 0 at the time-points that should not be used
  # yind is over all observations in long format
  ind0 <- !t(outer( xind, yind, limits) )
  yindHelp <- yind
  #  }  
  
  ### Compute the design matrix as sparse or normal matrix 
  ### depending on dimensions of the final design matrix
  MATRIX <- any(c(nrow(ind0), ncol(Bs)) > c(500, 50)) #MATRIX <- any(dim(X) > c(500, 50))
  MATRIX <- MATRIX && options("mboost_useMatrix")$mboost_useMatrix 
  
  if(MATRIX){
    # message("use sparse matrix in X_hist")
    diag <- Diagonal
    
    ## <enhancement> it would be nicer to construct the matrix directly as sparse matrix 
    
    X1des <- X1[id, ] 
    X1des[ind0] <- 0    
    X1des <- Matrix(X1des, sparse=TRUE) # convert into sparse matrix
    
  }else{ # small matrices: do not use Matrix
    
    X1des <- X1[id, ] 
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
  
  ## see Scheipl and Greven (2016) Identifiability in penalized function-on-function regression models  
  if(args$check.ident && args$inS == "smooth"){
    K1 <- diff(diag(ncol(Bs)), differences = args$differences)
    K1 <- crossprod(K1)
    # use the limits function to compute check measures on corresponding subsets of x(s) and B_j
    res_check <- check_ident(X1 = X1, L = L, Bs = Bs, K = K1, xname = xname, 
                             penalty = args$penalty, 
                             limits = args$limits, 
                             yind = yind, id = id, # yind is always in long format 
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
               # "smooth" = bsplines(yind, knots=args$knots[[varnames]]$knots, 
               #                   boundary.knots=args$knots[[varnames]]$boundary.knots, 
               #                   degree=args$degree),
               "smooth" = mboost_intern(yind, knots = args$knots[[varnames]]$knots, 
                                        boundary.knots = args$knots[[varnames]]$boundary.knots, 
                                        degree = args$degree, 
                                        fun = "bsplines"),
               "linear" = matrix(c(rep(1, length(yind)), yind), ncol=2),
               "constant"=  matrix(c(rep(1, length(yind))), ncol=1))

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
               " or set ", sQuote("center = TRUE"), ".", immediate.=TRUE)
    if (rns > args$df)
      stop("not possible to specify ", sQuote("df"),
           " smaller than the rank of the null space\n  ",
           "(unpenalized part of P-spline). Use larger value for ",
           sQuote("df"), " or set ", sQuote("center = TRUE"), ".")
  }
  
  # save matrices to compute numbers for identifiability checks
  args$Bs <- Bs
  args$X1des <- X1des
  args$K1 <- K1
  args$L <- L
  
  # tidy up workspace 
  rm(Bs, Bt, ind0, X1des, X1, L)
  
  return(list(X = X, K = K, args = args))
}



###############################################################################

#' Base-learners for Functional Covariates
#' 
#' Base-learners that fit historical functional effects that can be used with the 
#' tensor product, as, e.g., \code{hbistx(...) \%X\% bolsc(...)}, to form interaction 
#' effects (Ruegamer et al., 2018).  
#' For expert use only! May show unexpected behavior  
#' compared to other base-learners for functional data!
#' 
#' @param x object of type \code{hmatrix} containing time, index and functional covariate; 
#' note that \code{timeLab} in the \code{hmatrix}-object must be equal to 
#' the name of the time-variable in \code{timeformula} in the \code{FDboost}-call
#' @param knots either the number of knots or a vector of the positions 
#' of the interior knots (for more details see \code{\link[mboost]{bbs})}.
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
#' @param penalty by default, \code{penalty="ps"}, the difference penalty for P-splines is used, 
#' for \code{penalty="pss"} the penalty matrix is transformed to have full rank, 
#' so called shrinkage approach by Marra and Wood (2011)
#' @param check.ident use checks for identifiability of the effect, based on Scheipl and Greven (2016); 
#' see Brockhaus et al. (2017) for identifiability checks that take into account the integration limits
#' @param standard the historical effect can be standardized with a factor. 
#' "no" means no standardization, "time" standardizes with the current value of time and 
#' "lenght" standardizes with the lenght of the integral 
#' @param intFun specify the function that is used to compute integration weights in \code{s} 
#' over the functional covariate \eqn{x(s)}
#' @param inS historical effect can be smooth, linear or constant in s, 
#' which is the index of the functional covariates x(s). 
#' @param inTime historical effect can be smooth, linear or constant in time, 
#' which is the index of the functional response y(time). 
#' @param limits defaults to \code{"s<=t"} for an historical effect with s<=t;  
#' either one of \code{"s<t"} or \code{"s<=t"} for [l(t), u(t)] = [T1, t]; 
#' otherwise specify limits as a function for integration limits [l(t), u(t)]: 
#' function that takes \eqn{s} as the first and \code{t} as the second argument and returns 
#' \code{TRUE} for combinations of values (s,t) if \eqn{s} falls into the integration range for 
#' the given \eqn{t}.  
#' 
#' @details \code{bhistx} implements a base-learner for functional covariates with 
#' flexible integration limits \code{l(t)}, \code{r(t)} and the possibility to
#' standardize the effect by \code{1/t} or the length of the integration interval. 
#' The effect is \code{stand * int_{l(t)}^{r_{t}} x(s)beta(t,s) ds}. 
#' The base-learner defaults to a historical effect of the form 
#' \eqn{\int_{T1}^{t} x_i(s)beta(t,s) ds}, 
#' where \eqn{T1} is the minimal index of \eqn{t} of the response \eqn{Y(t)}. 
#' \code{bhistx} can only be used if \eqn{Y(t)} and \eqn{x(s)} are observd over
#' the same domain \eqn{s,t \in [T1, T2]}. 
#' The base-learner \code{bhistx} can be used to set up complex interaction effects 
#' like factor-specific historical  effects as discussed in Ruegamer et al. (2018). 
#' 
#' Note that the data has to be supplied as a \code{hmatrix} object for 
#' model fit and predictions. 
#' 
#' @return Equally to the base-learners of package mboost: 
#' 
#' An object of class \code{blg} (base-learner generator) with a 
#' \code{dpp} function (dpp, data pre-processing). 
#' 
#' The call of \code{dpp} returns an object of class 
#' \code{bl} (base-learner) with a \code{fit} function. The call to 
#' \code{fit} finally returns an object of class \code{bm} (base-model).
#' 
#' @seealso \code{\link{FDboost}} for the model fit and \code{\link{bhist}} 
#' for simple hisotorical effects. 
#'  
#' @keywords models
#' 
#' @references 
#' Brockhaus, S., Melcher, M., Leisch, F. and Greven, S. (2017): 
#' Boosting flexible functional regression models with a high number of functional historical effects,  
#' Statistics and Computing, 27(4), 913-926. 
#' 
#' Marra, G. and Wood, S.N. (2011): Practical variable selection for generalized additive models. 
#' Computational Statistics & Data Analysis, 55, 2372-2387.
#' 
#' Ruegamer D., Brockhaus, S., Gentsch K., Scherer, K., Greven, S. (2018). 
#' Boosting factor-specific functional historical models for the detection of synchronization in bioelectrical signals. 
#' Journal of the Royal Statistical Society: Series C (Applied Statistics), 67, 621-642.
#' 
#' Scheipl, F., Staicu, A.-M. and Greven, S. (2015): 
#' Functional Additive Mixed Models, Journal of Computational and Graphical Statistics, 24(2), 477-501.
#' \url{http://arxiv.org/abs/1207.5947} 
#' 
#' Scheipl, F. and Greven, S. (2016): Identifiability in penalized function-on-function regression models. 
#' Electronic Journal of Statistics, 10(1), 495-526. 
#'  
#' @examples 
#' if(require(refund)){
#' ## simulate some data from a historical model
#' ## the interaction effect is in this case not necessary
#' n <- 100
#' nygrid <- 35
#' data1 <- pffrSim(scenario = c("int", "ff"), limits = function(s,t){ s <= t }, 
#'                 n = n, nygrid = nygrid)
#' data1$X1 <- scale(data1$X1, scale = FALSE) ## center functional covariate                  
#' dataList <- as.list(data1)
#' dataList$tvals <- attr(data1, "yindex")
#'
#' ## create the hmatrix-object
#' X1h <- with(dataList, hmatrix(time = rep(tvals, each = n), id = rep(1:n, nygrid), 
#'                              x = X1, argvals = attr(data1, "xindex"), 
#'                              timeLab = "tvals", idLab = "wideIndex", 
#'                              xLab = "myX", argvalsLab = "svals"))
#' dataList$X1h <- I(X1h)   
#' dataList$svals <- attr(data1, "xindex")
#' ## add a factor variable 
#' dataList$zlong <- factor(gl(n = 2, k = n/2, length = n*nygrid), levels = 1:3)  
#' dataList$z <- factor(gl(n = 2, k = n/2, length = n), levels = 1:3)
#'
#' ## do the model fit with main effect of bhistx() and interaction of bhistx() and bolsc()
#' mod <- FDboost(Y ~ 1 + bhistx(x = X1h, df = 5, knots = 5) + 
#'                bhistx(x = X1h, df = 5, knots = 5) %X% bolsc(zlong), 
#'               timeformula = ~ bbs(tvals, knots = 10), data = dataList)
#'               
#' ## alternative parameterization: interaction of bhistx() and bols()
#' mod <- FDboost(Y ~ 1 + bhistx(x = X1h, df = 5, knots = 5) %X% bols(zlong), 
#'               timeformula = ~ bbs(tvals, knots = 10), data = dataList)
#'
#' \dontrun{
#'   # find the optimal mstop over 5-fold bootstrap (small example to reduce run time)
#'   cv <- cvrisk(mod, folds = cv(model.weights(mod), B = 5))
#'   mstop(cv)
#'   mod[mstop(cv)]
#'   
#'   appl1 <- applyFolds(mod, folds = cv(rep(1, length(unique(mod$id))), type = "bootstrap", B = 5))
#'
#'  # plot(mod)
#' }
#'}
#' 
#' @export
### P-spline base-learner for effect with time-specific integration limits using a hmatrix-object
bhistx <- function(x, 
                  limits = "s<=t", standard = c("no", "time", "length"), 
                  intFun = integrationWeightsLeft, 
                  inS = c("smooth","linear","constant"), inTime = c("smooth", "linear", "constant"),
                  knots = 10, boundary.knots = NULL, degree = 3, differences = 1, df = 4,
                  lambda = NULL, #center = FALSE, cyclic = FALSE
                  penalty = c("ps", "pss"), check.ident = FALSE
){
  
  index <- NULL # index is treated within hmatrix
  
  if (!is.null(lambda)) df <- NULL
  
  cll <- match.call()
  cll[[1]] <- as.name("bhistx")
  #print(cll)
  penalty <- match.arg(penalty)
  #print(penalty)
  
  standard <- match.arg(standard)
  
  inS <- match.arg(inS)
  inTime <- match.arg(inTime)
  #print(inS)
  
  if(!is.hmatrix(x)) stop("x has to be a hmatrix-object.")

  varnames <- all.vars(cll)[1] # only the first name, the name of x is a variable name
  
  # Reshape mfL so that it is the dataframe of the signal with 
  # the index of the signal and the index of the response as attributes
  xname <- getXLab(x)
  indname <- getArgvalsLab(x)
  indnameY <- getTimeLab(x)

  ## save the hmatrix-object into a data.frame
  mf <- data.frame(z=I(x))
  names(mf) <- varnames
  
  if(!all( abs(colMeans(getX(x), na.rm = TRUE)) < .Machine$double.eps*10^10)){
    message(xname, " is not centered per column, inducing a non-centered effect.")
  }
  
  vary <- ""
  
  # CC <- all(Complete.cases(mf))
  CC <- all(mboost_intern(mf, fun = "Complete.cases"))
  if (!CC)
    warning("base-learner contains missing values;\n",
            "missing values are excluded per base-learner, ",
            "i.e., base-learners may depend on different",
            " numbers of observations.")
  
  ## call X_histx in order to compute parameter settings, e.g. 
  ## the transformation matrix Z, shrinkage penalty, identifiability problems...
  # X_histx for data in long format, as hmatrix is always in long-format
  temp <- X_histx(mf, vary, 
                 args = hyper_histx(mf, vary, knots = knots, boundary.knots = boundary.knots, 
                                   degree = degree, differences = differences,
                                   df = df, lambda = lambda, center = FALSE, cyclic = FALSE,
                                   limits = limits, 
                                   standard = standard, intFun = intFun, 
                                   inS = inS, inTime = inTime, 
                                   penalty = penalty, check.ident = check.ident, 
                                   format="long"))
  temp$args$check.ident <- FALSE
  
  ret <- list(model.frame = function() mf,
    get_call = function(){
      cll <- deparse(cll, width.cutoff=500L)
      if (length(cll) > 1)
        cll <- paste(cll, collapse="")
      cll
    },
    get_data = function() mf,
    ## set index to NULL, as the index is treated within X_hist()
    ## get_index = function() index, 
    get_index = function() NULL,
    get_vary = function() vary,
    get_names = function(){
      varnames 
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
  
  ### call fitter with X_histx 
  # ret$dpp <- bl_lin(ret, Xfun = X_histx, args = temp$args) 
  ret$dpp <- mboost_intern(ret, Xfun = X_histx, args = temp$args, fun = "bl_lin")
  
  return(ret)
}
