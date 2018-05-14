#' Constrained row tensor product 
#' 
#' Combining single base-learners to form new, more complex base-learners, with
#' an identifiability constraint to center the interaction around the intercept and
#' around the two main effects. Suitable for functional response. 
#' @param bl1 base-learner 1, e.g. \code{bols(x1)}
#' @param bl2 base-learner 2, e.g. \code{bols(x2)}
#' 
#' @details Similar to \code{\%X\%} in package mboost, see \code{\link[mboost]{\%X\%}}, 
#' a row tensor product of linear base-learners is returned by \code{\%Xc\%}. 
#' \code{\%Xc\%} applies a sum-to-zero constraint to the design matrix suitable for
#' functional response if an interaction of two scalar covariates is specified 
#' in the case that the model contains a global intercept and both main effects, 
#' as the interaction is centerd around the intercept and centered around the two main effects. 
#' See Web Appendix A of Brockhaus et al. (2015) for details on how to enforce the constraint 
#' for the functional intercept.   
#' Use, e.g., in a model call to \code{FDboost}, following the scheme, 
#' \code{y ~ 1 + bolsc(x1) + bolsc(x2) + bols(x1) \%Xc\% bols(x2)}, 
#' where \code{1} induces a global intercept and \code{x1}, \code{x2} are factor variables, 
#' see Ruegamer et al. (2018).  
#' 
#' @references 
#' Brockhaus, S., Scheipl, F., Hothorn, T. and Greven, S. (2015): 
#' The functional linear array model. Statistical Modelling, 15(3), 279-300.
#' 
#' Ruegamer D., Brockhaus, S., Gentsch K., Scherer, K., Greven, S. (2018). 
#' Boosting factor-specific functional historical models for the detection of synchronization in bioelectrical signals. 
#' Journal of the Royal Statistical Society: Series C (Applied Statistics), 67, 621-642.
#' 
#' @author Sarah Brockhaus, David Ruegamer
#' 
#' @examples  
#' ######## Example for function-on-scalar-regression with interaction effect of two scalar covariates 
#' data("viscosity", package = "FDboost") 
#' ## set time-interval that should be modeled
#' interval <- "101"
#' 
#' ## model time until "interval" and take log() of viscosity
#' end <- which(viscosity$timeAll == as.numeric(interval))
#' viscosity$vis <- log(viscosity$visAll[,1:end])
#' viscosity$time <- viscosity$timeAll[1:end]
#' # with(viscosity, funplot(time, vis, pch = 16, cex = 0.2))
#' 
#' ## fit model with interaction that is centered around the intercept 
#' ## and the two main effects 
#' mod1 <- FDboost(vis ~ 1 + bolsc(T_C, df=1) + bolsc(T_A, df=1) + 
#'                 bols(T_C, df=1) %Xc% bols(T_A, df=1),
#'                 timeformula = ~bbs(time, df=6),
#'                 numInt = "equal", family = QuantReg(),
#'                 offset = NULL, offset_control = o_control(k_min = 9),
#'                 data = viscosity, control=boost_control(mstop = 100, nu = 0.4))
#'                 
#' ## check centering around intercept
#' colMeans(predict(mod1, which = 4))
#' 
#' ## check centering around main effects
#' colMeans(predict(mod1, which = 4)[viscosity$T_A == "low", ])
#' colMeans(predict(mod1, which = 4)[viscosity$T_A == "high", ])
#' colMeans(predict(mod1, which = 4)[viscosity$T_C == "low", ])
#' colMeans(predict(mod1, which = 4)[viscosity$T_C == "low", ])
#'
#' ## find optimal mstop using cvrsik() or validateFDboost()
#' ## ... 
#' 
#' ## look at interaction effect in one plot
#' # funplot(mod1$yind, predict(mod1, which=4))
#' 
#' @export
"%Xc%" <- function(bl1, bl2) {
  
  if (is.list(bl1) && !inherits(bl1, "blg"))
    return(lapply(bl1, "%Xc%", bl2 = bl2))
  
  if (is.list(bl2) && !inherits(bl2, "blg"))
    return(lapply(bl2, "%Xc%", bl1 = bl1))
  
  cll <- paste(bl1$get_call(), "%Xc%",
               bl2$get_call(), collapse = "")
  stopifnot(inherits(bl1, "blg"))
  stopifnot(inherits(bl2, "blg"))
  
  ## Check that the used base-learners contain an intercept
  used_bl <- c( deparse(match.call()$bl1[[1]]), 
                deparse(match.call()$bl2[[1]]) )
  if(any(used_bl == "bolsc")) stop("Use bols instead of bolsc with %Xc%.")
  if(any(used_bl == "brandomc")) stop("Use brandom instead of brandomc with %Xc%.")
  if(any(used_bl == "bbsc")) stop("Use bbs instead of bbsc with %Xc%.")
  if( (!is.null(match.call()$bl1$intercept) &&  match.call()$bl1$intercept != TRUE) |
      (!is.null(match.call()$bl2$intercept) &&  match.call()$bl2$intercept != TRUE) ){
    stop("Set intercept = TRUE in base-learners used with %Xc%.")
  }
  
  if(any(!used_bl %in% c("bols", "brandom", "bbs")) ){
    warning("%Xc% is intended to combine base-learners bols, brandom and bbs.")
  }
  
  stopifnot(!any(colnames(mboost_intern(bl1, fun = "model.frame.blg")) %in%
                   colnames(mboost_intern(bl2, fun = "model.frame.blg"))))
  mf <- cbind( mboost_intern(bl1, fun = "model.frame.blg"), 
               mboost_intern(bl2, fun = "model.frame.blg") )
  index1 <- bl1$get_index()
  index2 <- bl2$get_index()
  if (is.null(index1)) index1 <- 1:nrow(mf)
  if (is.null(index2)) index2 <- 1:nrow(mf)
  
  mfindex <- cbind(index1, index2)
  index <- NULL
  
  # CC <- all(Complete.cases(mf))
  CC <- all(mboost_intern(mf, fun = "Complete.cases"))
  if (!CC)
    warning("base-learner contains missing values;\n",
            "missing values are excluded per base-learner, ",
            "i.e., base-learners may depend on different",
            " numbers of observations.")
  ### option
  DOINDEX <- (nrow(mf) > options("mboost_indexmin")[[1]])
  if (is.null(index)) {
    if (!CC || DOINDEX) {
      index <- mboost_intern(mfindex, fun = "get_index")
      mf <- mf[index[[1]],,drop = FALSE]
      index <- index[[2]]
    }
  }
  
  vary <- ""
  
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
    ## <FIXME> Is this all we want to change if we set names here?
    set_names = function(value) attr(mf, "names") <<- value)
  ## </FIXME>
  class(ret) <- "blg"
  
  args1 <- environment(bl1$dpp)$args
  args2 <- environment(bl2$dpp)$args
  l1 <- args1$lambda
  l2 <- args2$lambda
  if (!is.null(l1) && !is.null(l2)) {
    args <- list(lambda = 1, df = NULL)
  } else {
    args <- list(lambda = NULL,
                 df = ifelse(is.null(args1$df), 1, args1$df) *
                   ifelse(is.null(args2$df), 1, args2$df))
  }
  
  Xfun <- function(mf, vary, args) {
    
    ## <SB> set prediciton to FALSE, if it is NULL
    if(is.null(args$prediction)) args$prediction <- FALSE
    
    newX1 <- environment(bl1$dpp)$newX
    newX2 <- environment(bl2$dpp)$newX
    
    X1 <- newX1(mf[, bl1$get_names(), drop = FALSE],
                prediction = args$prediction)
    K1 <- X1$K
    X1 <- X1$X
    if (!is.null(l1)) K1 <- l1 * K1
    MATRIX <- options("mboost_useMatrix")$mboost_useMatrix
    if (MATRIX & !is(X1, "Matrix"))
      X1 <- Matrix(X1)
    if (MATRIX & !is(K1, "Matrix"))
      K1 <- Matrix(K1)
    
    X2 <- newX2(mf[, bl2$get_names(), drop = FALSE],
                prediction = args$prediction)
    K2 <- X2$K
    X2 <- X2$X
    if (!is.null(l2)) K2 <- l2 * K2
    if (MATRIX & !is(X2, "Matrix"))
      X2 <- Matrix(X2)
    if (MATRIX & !is(K2, "Matrix"))
      K2 <- Matrix(K2)
    suppressMessages(
      X <- kronecker(X1, Matrix(1, ncol = ncol(X2),
                                dimnames = list("", colnames(X2))),
                     make.dimnames = TRUE) *
        kronecker(Matrix(1, ncol = ncol(X1),
                         dimnames = list("", colnames(X1))),
                  X2, make.dimnames = TRUE)
    )
    suppressMessages(
      K <- kronecker(K1, diag(ncol(X2))) +
        kronecker(diag(ncol(X1)), K2)
    )
    
    #----------------------------------
    ### <SB> Calculate constraints
    
    ## If model is fitted -> compute Z; but if model is predicted use the Z from the model fit
    ## if(!args$prediction){
    ## compute QR-decompotition only once
    if(is.null(args$Z)){
      ## put all effects of the two main effects + intercept into the constraints 
      # C <- t(X) %*% cbind(rep(1, nrow(X)), X1[ , -1], X2[ , -1])
      ## use whole matrices of marginal effects for constraints as Almond suggested  
      C <- t(X) %*% cbind(rep(1, nrow(X)), X1, X2)    
      qr_C <- qr(C)  ## , tol = 1e-10 ## time?
      if( any(class(qr_C) == "sparseQR") ){
        rank_C <- qr_C@Dim[2]
      }else{
        rank_C <- qr_C$rank 
      } 
      Q <- qr.Q(qr_C, complete=TRUE) # orthonormal matrix of QR decomposition
      args$Z <- Q[  , (rank_C + 1) : ncol(Q) ] # only keep last columns    
    } 
    
    ### Transform design and penalty matrix
    X <- X %*% args$Z
    K <- t(args$Z) %*% K %*% args$Z
    ## print(args$Z)
    #----------------------------------
    
    list(X = X, K = K, args = args)  
  }
  
  ## compute the transformation matrix Z
  temp <- Xfun(mf = mf, vary = vary, args = args)
  args$Z <- temp$args$Z
  #print(head(args$Z))
  #image(t(as.matrix(args$Z)))
  rm(temp)
  
  # ret$dpp <- bl_lin(ret, Xfun = Xfun, args = args)
  ret$dpp <- mboost_intern(ret, Xfun = Xfun, args = args, fun = "bl_lin")
  
  return(ret)
}



#############################################################################################
#############################################################################################

### workhorse for fitting matrix-response (Ridge-penalized) baselearners
### Y = kronecker(X2, X1)
### see Currie, Durban, Eilers (2006, JRSS B)
### args has extra argument isotropic TRUE/FALSE
bl_lin_matrix_a <- function(blg, Xfun, args) {
  
  mf <- blg$get_data()
  index <- blg$get_index()
  vary <- blg$get_vary()
  
  newX <- function(newdata = NULL, prediction = FALSE) {
    if (!is.null(newdata)) {
      mf <- mboost_intern(newdata, blg, mf, to.data.frame = FALSE, 
                          fun = "check_newdata")
    }
    ## this argument is currently only used in X_bbs --> bsplines
    args$prediction <- prediction
    return(Xfun(mf, vary, args))
  }
  X <- newX()
  K <- X$K
  X <- X$X
  c1 <- ncol(X$X1)
  c2 <- ncol(X$X2)
  n1 <- nrow(X$X1)
  n2 <- nrow(X$X2)
  
  G <- function(x) {
    one <- matrix(rep(1, ncol(x)), nrow = 1)
    suppressMessages(
      ret <- kronecker(x, one) * kronecker(one, x)
    )
    ret
  }
  
  dpp <- function(weights) {
    
    if (!is.null(attr(X$X1, "deriv")) || !is.null(attr(X$X2, "deriv")))
      stop("fitting of derivatives of B-splines not implemented")
    
    W <- matrix(weights, nrow = n1, ncol = n2)
    
    ### X = kronecker(X2, X1)
    XtX <- crossprod(G(X$X1), W) %*% G(X$X2)
    mymatrix <- matrix
    if (is(XtX, "Matrix")) mymatrix <- Matrix
    XtX <- array(XtX, c(c1, c1, c2, c2))
    XtX <- mymatrix(aperm(XtX, c(1, 3, 2, 4)), nrow = c1 * c2)
    
    ### If lambda was given in both baselearners, we
    ### directly multiply the marginal penalty matrices by lambda
    ### and then compute the total penalty as the kronecker sum.
    ### args$lambda is NA in this case and we don't compute
    ### the corresponding df's (unlike bl_lin)
    if (is.null(args$lambda)) {
      
      if(!is.null(args$isotropic) && !args$isotropic){  ## %A%
        
        ## get the design and penalty matrices of the marginal base-learners  
        # X1 <- X$X1
        # X2 <- X$X2
        # K1 <- args$K1 
        # K2 <- args$K2 
        
        #  ## per default do not expand the marginal design matrices 
        #  expand_index1 <- 1:nrow(X$X1)
        #  expand_index2 <- 1:nrow(X$X2)
        
        ## weights-matrix W: weights are for single observations in the matrix Y
        ## but the marginal bl work either on columns or rows of Y 
        
        ## easy case: all weights are 1
        if(all(W == 1)){ 
          w1 <- w2 <- 1 ## set weights for bl1 and bl2 to 1
          
          ## more difficult case:  different weights in W 
        }else{  

          ## use colMeans and rowMeans of W as weights
          ## weights in mboost are rescaled such that: 
          ## sum(W) == nrow(X$X1)*nrow(X$X2), if weight are not integres
          ## see mboost:::rescale_weights 
          w1 <- rowMeans(W)  ## rowSums(W) 
          w2 <- colMeans(W)  ## colSums(W)
          
          ## w1, w2 are correct, when the matrix W can be computed from them as 
          w1w2 <- w1 %*% t(w2)  ## all( W == w1w2 )

          if( any(abs(W - w1w2) > .Machine$double.eps*10^10) ){
            
            ## check whether w1w2 and W differ just by a factor 
            multFactor <- unique(c(W / w1w2)) 
            multFactor <- multFactor[!is.na(multFactor)] # get rid of NaN -> division by 0 
            ## multFactor is equal up to numerical inaccuracies
            if( all(abs(multFactor - multFactor[1] ) < .Machine$double.eps*10^10) ) multFactor <- multFactor[1]
            
            ## case that W and w1w2 just differ by a factor 
            if( all((W == w1w2)[w1w2 == 0]) & all((W == w1w2)[W == 0]) & ## check positions of zeros 
                length(multFactor) == 1 ){ # check that only 1 multiplicative factor 
              
              ## it is impossible to know whether multFactor is multiplied to w1 or w2! 
              # w2 <- w2 * multFactor
              # all( W == (w1 %*% t(w2 * multFactor)))
              # all( W == ( (w1 * multFactor) %*% t(w2)))
              
              ## just assume that resampling was done on the level of curves 
              w2 <- w2 * multFactor
              
              ## do not warn in the case that all rows of W are equal -> no resampling
              if(nrow(unique(W, MARGIN = 1)) != 1){
                warning("Assume that resampling is such that whole observations of blg1 are used.")
              }
              
            }else{
              warning("Set all weights = 1 for computation of lambda1 and lambda2 in %A%.")
              w1 <- w2 <- 1 
            } 
            
          } # end if that W != w1w2 
          
          ## more detailed warnings are given when necessary 
          # warning("rowMeans and colMeans of W do not multiply back to W, 
          #         thus anisotropic lambdas only roughly correct.")
          
          ### idea: blow up the two marginal design matrices and use weights on them
          ### but: this does not work correctly: problem with factor remains 
          #   ## W cannot be computed from w1 and w2, 
          #   ## -> blow up the marginal design matrices and use W with them, 
          #   expand_index1 <- rep(1:nrow(X$X1), times = nrow(X$X2)) 
          #   expand_index2 <- rep(1:nrow(X$X2), each = nrow(X$X1))  
          #   ## all( c(W) == weights) is TRUE, ordering of weights must match to blown-up marginal design matrices 
          #   ## standardize weights to compensate for the blow-up of the marginal design-matrices  
          #   #w1 <- c(W) / mean(rowSums(W)) ## for some special cases (e.g. BS on rows): mean(rowSums(W)) == nrow(X$X2)
          #   #w2 <- c(W) / mean(colSums(W)) ## for some special cases: mean(colSums(W)) == nrow(X$X1) 
          #   w1 <- c(W) / nrow(X$X2) ## for some special cases (e.g. BS on rows): mean(rowSums(W)) == nrow(X$X2)
          #   w2 <- c(W) / nrow(X$X1) ## for some special cases: mean(colSums(W)) == nrow(X$X1)
          
        } ## end computation of w1, w2
        
        
        ## case that df equals nr columns of design matrix -> no penalty -> lambda = 0
        if( abs(ncol(X$X1) - args$df1) < .Machine$double.eps*10^10 ){ 
          args$lambda1 <- 0 
        }else{
          ## call df2lambda for marginal bl1 
          al1 <- mboost_intern(X = X$X1,  # X$X1[expand_index1 , , drop = FALSE],
                               df = args$df1, lambda = NULL, ## lambda = args$df1, do not allow for lambda in %A%
                               dmat = args$K1, weights = w1, XtX = NULL, 
                               fun = "df2lambda")
          args$lambda1 <- al1["lambda"]
        }
        
        
        ## case that df equals nr columns of design matrix
        if( abs(ncol(X$X2) - args$df2) < .Machine$double.eps*10^10 ){ 
          args$lambda2 <- 0 
        }else{
          ## call df2lambda for marginal bl2 
          al2 <- mboost_intern(X = X$X2,  # X$X2[expand_index2 , , drop = FALSE],
                               df = args$df2, lambda = NULL,
                               dmat = args$K2, weights = w2, XtX = NULL, 
                               fun = "df2lambda")
          args$lambda2 <- al2["lambda"]
        }
        
        ## compute the penalty matrix which includes lambda1 and lambda2 for K1 and K2
        suppressMessages(
          K <- kronecker(args$lambda2 * args$K2, diag(ncol(X$X1))) +
            kronecker(diag(ncol(X$X2)), args$lambda1 * args$K1)
        ) 
        
      } ## end of anisotropic penalty 
      
      
      ### <FIXME>: is there a better way to feed XtX into lambdadf?
      lambdadf <- mboost_intern(X = diag(rankMatrix(X$X1, method = 'qr', warn.t = FALSE) *
                                           rankMatrix(X$X2, method = 'qr', warn.t = FALSE)),
                                df = args$df, lambda = args$lambda,
                                dmat = K, weights = weights, XtX = XtX, 
                                fun = "df2lambda")
      ### </FIXME>
      
      lambda <- lambdadf["lambda"]
      K <- lambda * K
      
      ## save the two marginal lambdas as attributes of lambdadf, as 
      ## lambdadf can be called via object$basemodel[[1]]$df()
      attr(lambdadf, "anisotropic") <- list(df1 = args$df1, lambda1 = args$lambda1, 
                                            df2 = args$df2, lambda2 = args$lambda2)
      
      
    } else {
      lambdadf <- args[c("lambda", "df")]
    }

    ### note: K already contains the lambda penalty parameter(s)
    XtX <- XtX + K
    
    ### nnls
    constr <- (!is.null(attr(X$X1, "constraint"))) +
      (!is.null(attr(X$X2, "constraint")))
    
    if (constr == 2)
      stop("only one dimension may be subject to constraints")
    constr <- constr > 0
    
    ## matrizes of class dgeMatrix are dense generic matrices; they should
    ## be coerced to class matrix and handled in the standard way
    if (is(XtX, "Matrix") && !extends(class(XtX), "dgeMatrix") && !extends(class(XtX), "dsyMatrix")) {
      XtXC <- Cholesky(forceSymmetric(XtX))
      mysolve <- function(y) {
        Y <- matrix(y, nrow = n1) * W
        if (constr)
          return( mboost_intern(X, as(XtXC, "matrix"), Y, fun = "nnls2D") )
        XWY <- as.vector(crossprod(X$X1, Y) %*% X$X2)
        solve(XtXC, XWY)  ## special solve routine from
        ## package Matrix
      }
    } else {
      if (is(XtX, "Matrix")) {
        ## coerce Matrix to matrix
        XtX <- as(XtX, "matrix")
      }
      mysolve <- function(y) {
        Y <- matrix(y, nrow = n1) * W
        if (constr)
          return( mboost_intern(X, as(XtX, "matrix"), Y, fun = "nnls2D") )
        XWY <- crossprod(X$X1, Y) %*% X$X2
        solve(XtX, matrix(as(XWY, "matrix"), ncol = 1),
              LINPACK = FALSE)
      }
    }
    
    cfprod <- function(b) tcrossprod(X$X1 %*% b, X$X2)
    
    fit <- function(y) {
      coef <- as(mysolve(y), "matrix")
      if (nrow(coef) != c1) coef <- matrix(as.vector(coef), nrow = c1)
      f <- cfprod(coef)
      f <- as(f, "matrix")
      if (options("mboost_Xmonotone")$mboost_Xmonotone) {
        md <- apply(f, 1, function(x) min(diff(x)))
        if (any(md < -(.Machine$double.eps)^(1/3))) {
          coef <- matrix(0, nrow = nrow(coef), ncol = ncol(coef))
          f <- matrix(0, nrow = nrow(f), ncol = ncol(f))
        }
      }
      ret <- list(model = coef,
                  fitted = function() as.vector(f))
      class(ret) <- c("bm_lin", "bm")
      ret
    }
    
    ### <FIXME> check for large n, option?
    hatvalues <- function() {
      return(NULL)
    }
    ### </FIXME>
    
    ### actually used degrees of freedom (trace of hat matrix)
    df <- function() lambdadf
    
    ### prepare for computing predictions
    predict <- function(bm, newdata = NULL, aggregate = c("sum", "cumsum", "none")) {
      cf <- lapply(bm, function(x) x$model)
      if(!is.null(newdata)) {
        index <- NULL
        X <- newX(newdata, prediction = TRUE)$X
      }
      ncfprod <- function(b)
        as.vector(as(tcrossprod(X$X1 %*% b, X$X2), "matrix"))
      aggregate <- match.arg(aggregate)
      pr <- switch(aggregate, "sum" = {
        cf2 <- 0
        for (b in cf) cf2 <- cf2 + b
        ncfprod(cf2)
      },
      "cumsum" = {
        cf2 <- 0
        ret <- c()
        for (b in cf) {
          cf2 <- cf2 + b
          ret <- cbind(ret, ncfprod(cf2))
        }
        ret
      },
      "none" = {
        ret <- c()
        for (b in cf) {
          ret <- cbind(ret, ncfprod(b))
        }
        ret
      })
      return(pr)
    }
    
    Xnames <- outer(colnames(X$X1), colnames(X$X2), paste, sep = "_")
    ret <- list(fit = fit, hatvalues = hatvalues,
                predict = predict, df = df,
                Xnames = as.vector(Xnames))
    class(ret) <- c("bl_lin", "bl")
    return(ret)
    
  } ## end dpp()
  
  return(dpp)
}



#' Kronecker product or row tensor product of two base-learners with anisotropic penalty 
#' 
#' Kronecker product or row tensor product of two base-learners allowing for anisotropic penalties. 
#' For the Kronecker product, \code{\%A\%} works in the general case, \code{\%A0\%} for the special case where 
#' the penalty is zero in one direction. 
#' For the row tensor product, \code{\%Xa0\%} works for the special case where 
#' the penalty is zero in one direction.  
#' 
#' @param bl1 base-learner 1, e.g. \code{bbs(x1)}
#' @param bl2 base-learner 2, e.g. \code{bbs(x2)}
#' 
#' @details 
#' When \code{\%O\%} is called with a specification of \code{df} in both base-learners, 
#' e.g. \code{bbs(x1, df = df1) \%O\% bbs(t, df = df2)}, the global \code{df} for the 
#' Kroneckered base-learner is computed as \code{df = df1 * df2}. 
#' And thus the penalty has only one smoothness parameter lambda resulting in an isotropic penalty, 
#' \deqn{P = lambda * [(P1 o I) + (I o P2)],} 
#' with overall penalty \eqn{P}, Kronecker product \eqn{o}, 
#' marginal penalty matrices \eqn{P1, P2} and identity matrices \eqn{I}.  
#' (Currie et al. (2006) introduced the generalized linear array model, which has a design matrix that 
#' is composed of the Kronecker product of two marginal design matrices, which was implemented in mboost 
#' as \code{\%O\%}.  
#' See Brockhaus et al. (2015) for the application of array models to functional data.)  
#' 
#' In contrast, a Kronecker product with anisotropic penalty is obtained by \code{\%A\%}, 
#' which allows for a different amount of smoothness in the two directions. 
#' For example \code{bbs(x1, df = df1) \%A\% bbs(t, df = df2)} results in computing two
#' different values for lambda for the two marginal design matrices and a global value of 
#' lambda to adjust for the global \code{df}, i.e. 
#' \deqn{P = lambda * [(lambda1 * P1 o I) +  (I o lambda2 * P2)],} 
#' with Kronecker product \eqn{o}, 
#' where \eqn{lambda1} is computed individually for \eqn{df1} and \eqn{P1}, 
#' \eqn{lambda2} is computed individually for \eqn{df2}  and \eqn{P2}, 
#' and \eqn{lambda} is computed such that the global \eqn{df} hold \eqn{df = df1 * df2}. 
#' For the computation of \eqn{lambda1} and \eqn{lambda2} weights specified in the model 
#' call can only be used when the weights, are such that they are specified on the level 
#' of rows and columns of the response matrix Y, e.g. resampling weights on the level of 
#' rows of Y and integration weights on the columns of Y are possible. 
#' If this the weights cannot be separated to blg1 and blg2 all
#' weights are set to 1 for the computation of \eqn{lambda1} and \eqn{lambda2} which implies that 
#' \eqn{lambda1} and \eqn{lambda2} are equal over 
#' folds of \code{cvrisk}. The computation of the global \eqn{lambda} considers the 
#' specified \code{weights}, such the global \eqn{df} are correct.    
#' 
#' The operator \code{\%A0\%} treats the important special case where \eqn{lambda1 = 0} or 
#' \eqn{lambda2 = 0}. In this case it suffices to compute the global lambda and computation gets
#' faster and arbitrary weights can be specified. Consider \eqn{lambda1 = 0} then the penalty becomes 
#' \deqn{P = lambda * [(1 * P1 o I) +  (I o lambda2 * P2)] = lambda * lambda2 * (I o P2),}  
#' and only one global \eqn{lambda} is computed which is then  \eqn{lambda * lambda2}.  
#'      
#' If the \code{formula} in \code{FDboost} contains base-learners connected by 
#' \code{\%O\%}, \code{\%A\%} or \code{\%A0\%}, 
#' those effects are not expanded with \code{timeformula}, allowing for model specifications 
#' with different effects in time-direction.  
#' 
#' \code{\%Xa0\%} computes like \code{\%X\%} the row tensor product of two base-learners, 
#' with the difference that it sets the penalty for one direction to zero. 
#' Thus, \code{\%Xa0\%} behaves to \code{\%X\%} analogously like \code{\%A0\%} to \code{\%O\%}. 
#'  
#' 
#' @references 
#' Brockhaus, S., Scheipl, F., Hothorn, T. and Greven, S. (2015): 
#' The functional linear array model. Statistical Modelling, 15(3), 279-300. 
#' 
#' Currie, I.D., Durban, M. and Eilers P.H.C. (2006):  
#' Generalized linear array models with applications to multidimensional smoothing. 
#' Journal of the Royal Statistical Society, Series B-Statistical Methodology, 68(2), 259-280.
#' 
#' @examples  
#' ######## Example for anisotropic penalty  
#' data("viscosity", package = "FDboost") 
#' ## set time-interval that should be modeled
#' interval <- "101"
#' 
#' ## model time until "interval" and take log() of viscosity
#' end <- which(viscosity$timeAll == as.numeric(interval))
#' viscosity$vis <- log(viscosity$visAll[,1:end])
#' viscosity$time <- viscosity$timeAll[1:end]
#' # with(viscosity, funplot(time, vis, pch = 16, cex = 0.2))
#' 
#' ## isotropic penalty, as timeformula is kroneckered to each effect using %O% 
#' ## only for the smooth intercept %A0% is used, as 1-direction should not be penalized 
#' mod1 <- FDboost(vis ~ 1 + 
#'                 bolsc(T_C, df = 1) + 
#'                 bolsc(T_A, df = 1) + 
#'                 bols(T_C, df = 1) %Xc% bols(T_A, df = 1),
#'                 timeformula = ~ bbs(time, df = 3),
#'                 numInt = "equal", family = QuantReg(),
#'                 offset = NULL, offset_control = o_control(k_min = 9),
#'                 data = viscosity, control=boost_control(mstop = 100, nu = 0.4))
#' ## cf. the formula that is passed to mboost
#' mod1$formulaMboost
#' 
#' ## anisotropic effects using %A0%, as lambda1 = 0 for all base-learners
#' ## in this case using %A% gives the same model, but three lambdas are computed explicitly 
#' mod1a <- FDboost(vis ~ 1 + 
#'                 bolsc(T_C, df = 1) %A0% bbs(time, df = 3) + 
#'                 bolsc(T_A, df = 1) %A0% bbs(time, df = 3) + 
#'                 bols(T_C, df = 1) %Xc% bols(T_A, df = 1) %A0% bbs(time, df = 3),
#'                 timeformula = ~ bbs(time, df = 3),
#'                 numInt = "equal", family = QuantReg(),
#'                 offset = NULL, offset_control = o_control(k_min = 9),
#'                 data = viscosity, control=boost_control(mstop = 100, nu = 0.4)) 
#' ## cf. the formula that is passed to mboost
#' mod1a$formulaMboost
#' 
#' ## alternative model specification by using a 0-matrix as penalty 
#' ## only works for bolsc() as in bols() one cannot specify K 
#' ## -> model without interaction term 
#' K0 <- matrix(0, ncol = 2, nrow = 2)
#' mod1k0 <- FDboost(vis ~ 1 + 
#'                  bolsc(T_C, df = 1, K = K0) +
#'                  bolsc(T_A, df = 1, K = K0), 
#'                  timeformula = ~ bbs(time, df = 3), 
#'                  numInt = "equal", family = QuantReg(), 
#'                  offset = NULL, offset_control = o_control(k_min = 9), 
#'                  data = viscosity, control=boost_control(mstop = 100, nu = 0.4))
#' ## cf. the formula that is passed to mboost
#' mod1k0$formulaMboost
#'                 
#' ## optimize mstop for mod1, mod1a and mod1k0
#' ## ...
#'                 
#' ## compare estimated coefficients
#' \dontrun{
#' par(mfrow=c(4, 2))
#' plot(mod1, which = 1)
#' plot(mod1a, which = 1)
#' plot(mod1, which = 2)
#' plot(mod1a, which = 2)
#' plot(mod1, which = 3)
#' plot(mod1a, which = 3)
#' funplot(mod1$yind, predict(mod1, which=4))
#' funplot(mod1$yind, predict(mod1a, which=4))
#' }
#' 
#' @name anisotropic_Kronecker
NULL



#' @rdname anisotropic_Kronecker
#' @export
"%A%" <- function(bl1, bl2) {

  ### code snippet from %O% in package mboost
  #   if (is.list(bl1) && !inherits(bl1, "blg"))
  #     return(lapply(bl1, "%X%", bl2 = bl2))
  #   
  #   if (is.list(bl2) && !inherits(bl2, "blg"))
  #     return(lapply(bl2, "%X%", bl1 = bl1))
  
  cll <- paste(bl1$get_call(), "%A%",
               bl2$get_call(), collapse = "")
  stopifnot(inherits(bl1, "blg"))
  stopifnot(inherits(bl2, "blg"))
  
  mf1 <- mboost_intern(bl1, fun = "model.frame.blg")
  mf2 <- mboost_intern(bl2, fun = "model.frame.blg")
  stopifnot(!any(colnames(mf1) %in%
                   colnames(mf2)))
  mf <- c(mf1, mf2)
  stopifnot(all(complete.cases(mf[[1]])))
  stopifnot(all(complete.cases(mf[[2]])))
  
  index <- NULL
  
  vary <- ""
  
  ret <- list(model.frame = function()
    return(mf),
    get_call = function(){
      cll <- deparse(cll, width.cutoff=500L)
      if (length(cll) > 1)
        cll <- paste(cll, collapse="")
      cll
    },
    get_data = function() mf,
    get_index = function() index,
    get_vary = function() vary,
    get_names = function() names(mf),
    ## <FIXME> Is this all we want to change if we set names here?
    set_names = function(value) attr(mf, "names") <<- value)
  ## </FIXME>
  class(ret) <- "blg"
  
  args1 <- environment(bl1$dpp)$args
  args2 <- environment(bl2$dpp)$args
  l1 <- args1$lambda
  l2 <- args2$lambda
  if (xor(is.null(l1), is.null(l2)))
    stop("you cannot mix lambda and df in ",
         sQuote("%A%"))
  
  if (!is.null(l1) && !is.null(l2)) {
    ### there is no common lambda!
    args <- list(lambda = NA, df = NA)
    
  }else{

    ### <SB> anisotropic penalty matrix
    ### save lambda, df and penalty matrices of the marginal base-learners 
    args <- list(lambda = NULL,
                 df = ifelse(is.null(args1$df), 1, args1$df) *
                   ifelse(is.null(args2$df), 1, args2$df), 
                 lambda1 = args1$lambda,  
                 df1 = args1$df, 
                 K1 = environment(bl1$dpp)$K, 
                 lambda2 = args2$lambda,
                 df2 = args2$df, 
                 K2 = environment(bl2$dpp)$K)
  }
  
  Xfun <- function(mf, vary, args) {

    newX1 <- environment(bl1$dpp)$newX
    newX2 <- environment(bl2$dpp)$newX
    
    X1 <- newX1(as.data.frame(mf[bl1$get_names()]),
                prediction = args$prediction)
    K1 <- X1$K
    X1 <- X1$X
    if (!is.null(l1)) K1 <- l1 * K1
    MATRIX <- options("mboost_useMatrix")$mboost_useMatrix
    if (MATRIX & !is(X1, "Matrix"))
      X1 <- Matrix(X1)
    if (MATRIX & !is(K1, "Matrix"))
      K1 <- Matrix(K1)
    
    X2 <- newX2(as.data.frame(mf[bl2$get_names()]),
                prediction = args$prediction)
    K2 <- X2$K
    X2 <- X2$X
    if (!is.null(l2)) K2 <- l2 * K2
    if (MATRIX & !is(X2, "Matrix"))
      X2 <- Matrix(X2)
    if (MATRIX & !is(K2, "Matrix"))
      K2 <- Matrix(K2)
    suppressMessages(
      K <- kronecker(K2, diag(ncol(X1))) +
        kronecker(diag(ncol(X2)), K1)
    )
    list(X = list(X1 = X1, X2 = X2), K = K)
  }
  
  args$isotropic <- FALSE
  
  # ret$dpp <- mboost:::bl_lin_matrix(ret, Xfun = Xfun, args = args)
  ret$dpp <- bl_lin_matrix_a(ret, Xfun = Xfun, args = args)
  
  return(ret)
}



###################################################################################


# Kronecker product of two base-learners with penalty in one direction  
# Only works for the special case were lambda1 or lambda2 is 0. 
# Computes only one global lambda for the penalty. 
#' @rdname anisotropic_Kronecker
#' @export
"%A0%" <- function(bl1, bl2) {
  
  #   if (is.list(bl1) && !inherits(bl1, "blg"))
  #     return(lapply(bl1, "%X%", bl2 = bl2))
  #   
  #   if (is.list(bl2) && !inherits(bl2, "blg"))
  #     return(lapply(bl2, "%X%", bl1 = bl1))
  
  cll <- paste(bl1$get_call(), "%A0%",
               bl2$get_call(), collapse = "")
  stopifnot(inherits(bl1, "blg"))
  stopifnot(inherits(bl2, "blg"))
  
  mf1 <- mboost_intern(bl1, fun = "model.frame.blg")
  mf2 <- mboost_intern(bl2, fun = "model.frame.blg")
  stopifnot(!any(colnames(mf1) %in%
                   colnames(mf2)))
  mf <- c(mf1, mf2)
  stopifnot(all(complete.cases(mf[[1]])))
  stopifnot(all(complete.cases(mf[[2]])))
  
  index <- NULL
  
  vary <- ""
  
  ret <- list(model.frame = function()
    return(mf),
    get_call = function(){
      cll <- deparse(cll, width.cutoff=500L)
      if (length(cll) > 1)
        cll <- paste(cll, collapse="")
      cll
    },
    get_data = function() mf,
    get_index = function() index,
    get_vary = function() vary,
    get_names = function() names(mf),
    ## <FIXME> Is this all we want to change if we set names here?
    set_names = function(value) attr(mf, "names") <<- value)
  ## </FIXME>
  class(ret) <- "blg"
  
  args1 <- environment(bl1$dpp)$args
  args2 <- environment(bl2$dpp)$args
  l1 <- args1$lambda
  l2 <- args2$lambda
  if (xor(is.null(l1), is.null(l2)))
    stop("you cannot mix lambda and df in ",
         sQuote("%A0%"))
  if (!is.null(l1) && !is.null(l2)) {
    ### there is no common lambda!
    args <- list(lambda = NA, df = NA)
  } else {
    
    ### <SB> anisotropic penalty matrix 
    df1 <- args1$df
    df2 <- args2$df
    
    args <- list(lambda = NULL,
                 df = ifelse(is.null(df1), 1, df1) *
                   ifelse(is.null(df2), 1, df2))

    ## case that df equals nr columns of design matrix -> no penalty -> lambda = 0
    if( ncol(environment(bl1$dpp)$X) - df1 < .Machine$double.eps*10^10){ 
      args$lambda1 <- 0 
      if( ncol(environment(bl1$dpp)$X) - df1 < - .Machine$double.eps*10^10){
        warning("Specified df in ", bl1$get_call(), " are higher than the number of columns, ",  
                "which is ", ncol(environment(bl1$dpp)$X), ".")
      }
    }else{
      args$lambda1 <- 1
    }
    
    ## case that df equals nr columns of design matrix
    if( ncol(environment(bl2$dpp)$X) - df2 < .Machine$double.eps*10^10){ 
      args$lambda2 <- 0 
      if( ncol(environment(bl2$dpp)$X) - df2 < - .Machine$double.eps*10^10){
        warning("Specified df in ", bl2$get_call(), " are higher than the number of columns, ",  
                "which is ", ncol(environment(bl2$dpp)$X), ".")
      }
    }else{
      args$lambda2 <- 1
    }

    if(args$lambda1 != 0 & args$lambda2 != 0) 
      stop("%A0% can only be used when smoothing parameter is zero for one direction.")
    
    l1 <- args$lambda1
    l2 <- args$lambda2 
    
  }
  
  Xfun <- function(mf, vary, args) {
    
    newX1 <- environment(bl1$dpp)$newX
    newX2 <- environment(bl2$dpp)$newX
    
    X1 <- newX1(as.data.frame(mf[bl1$get_names()]),
                prediction = args$prediction)
    K1 <- X1$K
    X1 <- X1$X
    if (!is.null(l1)) K1 <- l1 * K1
    MATRIX <- options("mboost_useMatrix")$mboost_useMatrix
    if (MATRIX & !is(X1, "Matrix"))
      X1 <- Matrix(X1)
    if (MATRIX & !is(K1, "Matrix"))
      K1 <- Matrix(K1)
    
    X2 <- newX2(as.data.frame(mf[bl2$get_names()]),
                prediction = args$prediction)
    K2 <- X2$K
    X2 <- X2$X
    if (!is.null(l2)) K2 <- l2 * K2
    if (MATRIX & !is(X2, "Matrix"))
      X2 <- Matrix(X2)
    if (MATRIX & !is(K2, "Matrix"))
      K2 <- Matrix(K2)
    suppressMessages(
      K <- kronecker(K2, diag(ncol(X1))) +
        kronecker(diag(ncol(X2)), K1)
    )
    list(X = list(X1 = X1, X2 = X2), K = K)
  }
  
  ret$dpp <- mboost_intern(ret, Xfun = Xfun, args = args, 
                           fun = "bl_lin_matrix")
  
  return(ret)
}




# Row tensor product of two base-learners with penalty in one direction  
# Only works for the special case were lambda1 or lambda2 is 0. 
# Computes only one global lambda for the penalty. 
#' @rdname anisotropic_Kronecker
#' @export
"%Xa0%" <- function(bl1, bl2) {

  if (is.list(bl1) && !inherits(bl1, "blg"))
    return(lapply(bl1, "%Xa0%", bl2 = bl2))
  
  if (is.list(bl2) && !inherits(bl2, "blg"))
    return(lapply(bl2, "%Xa0%", bl1 = bl1))
  
  cll <- paste(bl1$get_call(), "%Xa0%",
               bl2$get_call(), collapse = "")
  stopifnot(inherits(bl1, "blg"))
  stopifnot(inherits(bl2, "blg"))
  
  stopifnot(!any(colnames(mboost_intern(bl1, fun = "model.frame.blg")) %in%
                   colnames(mboost_intern(bl2, fun = "model.frame.blg"))))
  mf <- cbind( mboost_intern(bl1, fun = "model.frame.blg"), 
               mboost_intern(bl2, fun = "model.frame.blg") )
  
  index1 <- bl1$get_index()
  index2 <- bl2$get_index()
  if (is.null(index1)) index1 <- 1:nrow(mf)
  if (is.null(index2)) index2 <- 1:nrow(mf)
  
  mfindex <- cbind(index1, index2)
  index <- NULL
  
  # CC <- all(Complete.cases(mf))
  CC <- all(mboost_intern(mf, fun = "Complete.cases"))
  if (!CC)
    warning("base-learner contains missing values;\n",
            "missing values are excluded per base-learner, ",
            "i.e., base-learners may depend on different",
            " numbers of observations.")
  ### option
  DOINDEX <- (nrow(mf) > options("mboost_indexmin")[[1]])
  if (is.null(index)) {
    if (!CC || DOINDEX) {
      index <- mboost_intern(mfindex, fun = "get_index")
      mf <- mf[index[[1]],,drop = FALSE]
      index <- index[[2]]
    }
  }
  
  vary <- ""
  
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
    ## <FIXME> Is this all we want to change if we set names here?
    set_names = function(value) attr(mf, "names") <<- value)
  ## </FIXME>
  class(ret) <- "blg"
  
  args1 <- environment(bl1$dpp)$args
  args2 <- environment(bl2$dpp)$args
  l1 <- args1$lambda
  l2 <- args2$lambda
  
  if (xor(is.null(l1), is.null(l2)))
    stop("you cannot mix lambda and df in ",
         sQuote("%Xa0%"))
  
  if (!is.null(l1) && !is.null(l2)) {
    args <- list(lambda = 1, df = NULL)
    
  } else {

    ### <SB> anisotropic penalty matrix 
    df1 <- args1$df
    df2 <- args2$df
    
    args <- list(lambda = NULL,
                 df = ifelse(is.null(df1), 1, df1) * 
                   ifelse(is.null(df2), 1, df2))
    
    ## case that df equals nr columns of design matrix -> no penalty -> lambda = 0
    if( ncol(environment(bl1$dpp)$X) - df1 < .Machine$double.eps*10^10){ 
      args$lambda1 <- 0 
      if( ncol(environment(bl1$dpp)$X) - df1 < - .Machine$double.eps*10^10){
        warning("Specified df in ", bl1$get_call(), " are higher than the number of columns, ",  
                "which is ", ncol(environment(bl1$dpp)$X), ".")
      }
    }else{
      args$lambda1 <- 1
    }
    
    ## case that df equals nr columns of design matrix
    if( ncol(environment(bl2$dpp)$X) - df2 < .Machine$double.eps*10^10){ 
      args$lambda2 <- 0 
      if( ncol(environment(bl2$dpp)$X) - df2 < - .Machine$double.eps*10^10){
        warning("Specified df in ", bl2$get_call(), " are higher than the number of columns, ",  
                "which is ", ncol(environment(bl2$dpp)$X), ".")
      }
    }else{
      args$lambda2 <- 1
    }
    
    if(args$lambda1 != 0 & args$lambda2 != 0) 
      stop("%Xa0% can only be used when smoothing parameter is zero for one direction.")
    
    l1 <- args$lambda1
    l2 <- args$lambda2 

  }
  
  Xfun <- function(mf, vary, args) {

    newX1 <- environment(bl1$dpp)$newX
    newX2 <- environment(bl2$dpp)$newX
    
    X1 <- newX1(mf[, bl1$get_names(), drop = FALSE],
                prediction = args$prediction)
    K1 <- X1$K
    X1 <- X1$X
    if (!is.null(l1)) K1 <- l1 * K1
    MATRIX <- options("mboost_useMatrix")$mboost_useMatrix
    if (MATRIX & !is(X1, "Matrix"))
      X1 <- Matrix(X1)
    if (MATRIX & !is(K1, "Matrix"))
      K1 <- Matrix(K1)
    
    X2 <- newX2(mf[, bl2$get_names(), drop = FALSE],
                prediction = args$prediction)
    K2 <- X2$K
    X2 <- X2$X
    if (!is.null(l2)) K2 <- l2 * K2
    if (MATRIX & !is(X2, "Matrix"))
      X2 <- Matrix(X2)
    if (MATRIX & !is(K2, "Matrix"))
      K2 <- Matrix(K2)
    suppressMessages(
      X <- kronecker(X1, Matrix(1, ncol = ncol(X2),
                                dimnames = list("", colnames(X2))),
                     make.dimnames = TRUE) *
        kronecker(Matrix(1, ncol = ncol(X1),
                         dimnames = list("", colnames(X1))),
                  X2, make.dimnames = TRUE)
    )
    suppressMessages(
      K <- kronecker(K1, diag(ncol(X2))) +
        kronecker(diag(ncol(X1)), K2)
    )
    list(X = X, K = K)
  }
  
  ret$dpp <- mboost_intern(ret, Xfun = Xfun, args = args, fun = "bl_lin")
  
  return(ret)
}
