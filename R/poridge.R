#' @export
Predict.matrix.pco.smooth <- function(object, data){

  # somehow data gets passed as a list here :-\
  dat <- matrix(data[[object$term]], ncol=object$dim)

  return(dat)
}
#' Faster multi-dimensional scaling
#'
#' This is a modified version of \code{\link{cmdscale}} that uses the Lanczos
#' procedure (\code{\link[mgcv]{slanczos}}) instead of \code{eigen}. Called by
#' \code{\link{smooth.construct.pco.smooth.spec}}.
#'
#' @param d a distance structure as returned by \code{\link{dist}}, or a full
#'   symmetric matrix of distances or dissimilarities.
#' @param k the maximum dimension of the space which the data are to be
#'   represented in; must be in \code{\{1, 2, ..., n-1\}}.
#' @param eig logical indicating whether eigenvalues should be returned.
#' @param add logical indicating if the additive constant of Cailliez (1983)
#'   should be computed, and added to the non-diagonal dissimilarities such that
#'   the modified dissimilarities are Euclidean.
#' @param x.ret indicates whether the doubly centred symmetric distance matrix
#'   should be returned.
#' @return as \code{\link{cmdscale}}
#'
#' @author David L Miller, based on code by R Core.
#' @seealso \code{\link{smooth.construct.pco.smooth.spec}}
#' @importFrom mgcv slanczos
#' @export
#' @references
#' Cailliez, F. (1983). The analytical solution of the additive constant problem.
#' \emph{Psychometrika}, 48, 343-349.
cmdscale_lanczos <- function(d, k = 2, eig = FALSE, add = FALSE, x.ret = FALSE){

  if (anyNA(d))
    stop("NA values not allowed in 'd'")
  if (is.null(n <- attr(d, "Size"))) {
    if(add) d <- as.matrix(d)
    x <- as.matrix(d^2)
    storage.mode(x) <- "double"
    if ((n <- nrow(x)) != ncol(x))
      stop("distances must be result of 'dist' or a square matrix")
    rn <- rownames(x)
  } else {
    rn <- attr(d, "Labels")
    x <- matrix(0, n, n) # must be double
    if (add) d0 <- x
    x[row(x) > col(x)] <- d^2
    x <- x + t(x)
    if (add) {
      d0[row(x) > col(x)] <- d
      d <- d0 + t(d0)
    }
  }
  n <- as.integer(n)

  ## we need to handle nxn internally in dblcen
  if(is.na(n) || n > 46340) stop("invalid value of 'n'")
  if((k <- as.integer(k)) > n - 1 || k < 1)
    stop("'k' must be in {1, 2, .. n - 1}")

  ## NB: this alters argument x, which is OK as it is re-assigned.
  #x <- .Call(stats::C_DoubleCentre, x)
  x <- scale(t(scale(t(x), scale=FALSE)),scale=FALSE)

  if(add) { ## solve the additive constant problem
    ## it is c* = largest eigenvalue of 2 x 2 (n x n) block matrix Z:
    i2 <- n + (i <- 1L:n)
    Z <- matrix(0, 2L*n, 2L*n)
    Z[cbind(i2,i)] <- -1
    Z[ i, i2] <- -x
#    Z[i2, i2] <- .Call(stats::C_DoubleCentre, 2*d)
    Z[i2, i2] <- scale(t(scale(t(2*d), scale=FALSE)),scale=FALSE)

    ###### this is where Dave modified things
    add.c <- max(slanczos(Z, k=1, kl=1)$values)
    #e <- eigen(Z, symmetric = FALSE, only.values = TRUE)$values
    #add.c <- max(Re(e))
    ## and construct a new x[,] matrix:
    x <- matrix(double(n*n), n, n)
    non.diag <- row(d) != col(d)
    x[non.diag] <- (d[non.diag] + add.c)^2
    #x <- .Call(stats::C_DoubleCentre, x)
    x <- scale(t(scale(t(x), scale=FALSE)),scale=FALSE)
  }

  ###### this is where Dave modified things
  e <- slanczos(-x/2, k=k)
  ev <- e$values#[seq_len(k)]
  evec <- e$vectors#[, seq_len(k), drop = FALSE]
  k1 <- sum(ev > 0)

  if(k1 < k) {
    warning(gettextf("only %d of the first %d eigenvalues are > 0", k1, k),
            domain = NA)
    evec <- evec[, ev > 0, drop = FALSE]
    ev <- ev[ev > 0]
  }

  points <- evec * rep(sqrt(ev), each=n)
  dimnames(points) <- list(rn, NULL)
  if (eig || x.ret || add) {
    evalus <- e$values # Cox & Cox have sum up to n-1, though
    list(points = points, eig = if(eig) evalus, x = if(x.ret) x,
         ac = if(add) add.c else 0,
         GOF = sum(ev)/c(sum(abs(evalus)), sum(pmax(evalus, 0))) )
  } else points
}
#' Make predictions using pco basis terms
#'
#' This function performs the necessary preprocessing for making predictions
#' with \code{\link[mgcv]{gam}} models that include \code{\link{pco}} basis
#' terms. The function \code{pco_predict_preprocess} builds a \code{data.frame}
#' (or augments an existing one) to be used with the usual \code{predict}
#' function.
#'
#' Models with \code{\link{pco}} basis terms are fitted by inputting distances
#' among the observations and then regressing (with a ridge penalty) on leading
#' principal coordinates arising from these distances. To perform prediction, we
#' must input the distances from the new data points to the original points, and
#' then "insert" the former into the principal coordinate space by the
#' interpolation method of Gower (1968) (see also Miller, 2012).
#'
#' An example of how to use this function in practice is shown in
#' \code{\link{smooth.construct.pco.smooth.spec}}.
#'
#' @param model a fitted \code{\link[mgcv]{gam}} model with at least one term of
#'   class "\code{pco.smooth}".
#' @param newdata data frame including the new values for any
#'   non-\code{\link{pco}} terms in the original fit. If there were none, this
#'   can be left as \code{NULL}.
#' @param dist_list a list of \code{n} \eqn{\times} \code{n*} matrices, one per
#'   \code{\link{pco}} term in the model, giving the distances from the
#'   \code{n*} prediction points to the \code{n} design points (original
#'   observations). List entry names should correspond to the names of the terms
#'   in the model (e.g., if the model includes a \code{s(x)} term,
#'   \code{dist_list} must include an element named "\code{x}").
#'
#' @return a \code{\link{data.frame}} with the coordinates for the new data
#'   inserted into principal coordinate space, in addition to the supplied
#'   \code{newdata} if this was non-\code{NULL}. This can be used as the
#'   \code{newdata} argument in a call to \code{\link[mgcv]{predict.gam}}.
#' @author David L Miller
#' @export
#' @references Gower, J. C. (1968). Adding a point to vector diagrams in
#' multivariate analysis. Biometrika, 55(3), 582-585.
#' \url{http://doi.org/10.2307/2334268}
#'
#' Miller, D. L. (2012). On smooth models for complex domains and distances. PhD
#' dissertation, Department of Mathematical Sciences, University of Bath.
#' @seealso \code{\link{smooth.construct.pco.smooth.spec}}
pco_predict_preprocess <- function(model, newdata=NULL, dist_list){

  # populate newdata
  destroy_col <- FALSE
  if(is.null(newdata)){
    newdata <- data.frame(newdata_dummy_data = rep(NA,ncol(dist_list[[1]])))
    # reminder to destroy that extra column
    destroy_col <- TRUE
  }

  # which terms in the model are pco terms?
  which_pco <- which(unlist(lapply(model$smooth,
                             function(x) any(class(x)=="pco.smooth"))))

  # die if there are no pco terms
  if(length(which_pco)==0){
    stop("There are no pco smooths in the model")
  }

  # loop over smooths of the right type
  for(i in which_pco){

    # this term
    term_lab <- model$smooth[[i]]$term

    # get this distance matrix
    distmat <- dist_list[[term_lab]]

    # goofy - build the data to put into the equation later on
    mds.obj <- model$smooth[[i]]$xt$mds.obj
    X <- mds.obj$x
    S <- -1/2*X
    d <- distmat
    non.diag <- row(d) != col(d)
    d[non.diag] <- (d[non.diag] + mds.obj$ac)
    d <- -(d^2-diag(S))

    # lambda^-1, 1/eigenvalues in a diagonal matrix
    ev <- mds.obj$eig[1:ncol(mds.obj$points)]
    if(length(ev)>1){
      lambda.inverse <- diag(1/ev)
    }else{
      lambda.inverse <- as.matrix(1/ev)
    }

    # conversion from maths into code of equation 10 in Gower 1968
    mds.points <- t(1/2*(lambda.inverse %*% t(mds.obj$points) %*% d))

    # make the prediction frame for this term
    preddat <- as.data.frame(mds.points)

    # give the columns names
    names(preddat) <- paste0("xdummy", "_" ,1:ncol(mds.points))

    # splodge it into the prediction frame
    newdata[[term_lab]] <- mds.points
  }

  # get rid of the column we needed to setup the data
  if(destroy_col){
    newdata[["newdata_dummy_data"]] <- NULL
  }

  return(newdata)
}

#'Principal coordinate ridge regression
#'
#'Smooth constructor function for principal coordinate ridge regression fitted
#'by \code{\link[mgcv]{gam}}. When the principal coordinates are defined by a
#'relevant distance among functional predictors, this is a form of nonparametric
#'scalar-on-function regression. Reiss et al. (2016) describe the approach and
#'apply it to dynamic time warping distances among functional predictors.
#'
#'@aliases pco smooth.construct.pco.smooth.spec Predict.matrix.pco.smooth
#'  poridge
#'@export
#'@importFrom stats cmdscale
#'
#'@param object a smooth specification object, usually generated by a term of
#'  the form \code{s(dummy, bs="pco", k, xt)}; see Details.
#'@param data a list containing just the data.
#'@param knots IGNORED!
#'
#'@return An object of class \code{pco.smooth}. The resulting object has an
#'  \code{xt} element which contains details of the multidimensional scaling,
#'  which may be interesting.
#'
#'@section Details: The constructor is not normally called directly, but is
#'  rather used internally by \code{\link{gam}}.
#'
#'  In a \code{\link[mgcv]{gam}} term of the above form \code{s(dummy, bs="pco",
#'  k, xt)}, \itemize{ \item \code{dummy} is an arbitrary vector (or name of a
#'  column in \code{data}) whose length is the number of observations. This is
#'  not actually used, but is required as part of the input to
#'  \code{\link[mgcv]{s}}. Note that if multiple \code{pco} terms are used in
#'  the model, there must be multiple unique term names (e.g., "\code{dummy1}",
#'  "\code{dummy2}", etc). \item \code{k} is the number of principal coordinates
#'  (e.g., \code{k=9} will give a 9-dimensional projection of the data). \item
#'  \code{xt} is a list supplying the distance information, in one of two ways.
#'  (i) A matrix \code{Dmat} of distances can be supplied directly via
#'  \code{xt=list(D=Dmat,\dots)}. (ii) Alternatively, one can use
#'  \code{xt=list(realdata=\dots, dist_fn=\dots, \dots)} to specify a data
#'  matrix \code{realdata} and distance function \code{dist_fn}, whereupon a
#'  distance matrix \code{dist_fn(realdata)} is created. } The list \code{xt}
#'  also has the following optional elements: \itemize{ \item \code{add}: Passed
#'  to \code{\link{cmdscale}} when performing multidimensional scaling; for
#'  details, see the help for that function. (Default \code{FALSE}.)\cr \item
#'  \code{fastcmd}: if \code{TRUE}, multidimensional scaling is performed by
#'  \code{\link{cmdscale_lanczos}}, which uses Lanczos iteration to
#'  eigendecompose the distance matrix; if \code{FALSE}, MDS is carried out by
#'  \code{\link{cmdscale}}. Default is \code{FALSE}, to use \code{cmdscale}. }
#'
#'@author David L Miller, based on code from Lan Huo and Phil Reiss
#'
#'@references Reiss, P. T., Miller, D. L., Wu, P.-S., and Wen-Yu Hua, W.-Y.
#'Penalized nonparametric scalar-on-function regression via principal
#'coordinates. Under revision. Available at
#'\url{https://works.bepress.com/phil_reiss/42/}.
#'
#' @examples
#' \dontrun{
#' # a simulated example
#' library(refund)
#' library(mgcv)
#' require(dtw)
#'
#' ## First generate the data
#' Xnl <- matrix(0, 30, 101)
#' set.seed(813)
#' tt <- sort(sample(1:90, 30))
#' for(i in 1:30){
#'   Xnl[i, tt[i]:(tt[i]+4)] <- -1
#'   Xnl[i, (tt[i]+5):(tt[i]+9)] <- 1
#' }
#' X.toy <- Xnl + matrix(rnorm(30*101, ,0.05), 30)
#' y.toy <- tt + rnorm(30, 0.05)
#' y.rainbow <- rainbow(30, end=0.9)[(y.toy-min(y.toy))/
#'                                    diff(range(y.toy))*29+1]
#'
#' ## Display the toy data
#' par(mfrow=c(2, 2))
#' matplot((0:100)/100, t(Xnl[c(4, 25), ]), type="l", xlab="t", ylab="",
#'         ylim=range(X.toy), main="Noiseless functions")
#' matplot((0:100)/100, t(X.toy[c(4, 25), ]), type="l", xlab="t", ylab="",
#'         ylim=range(X.toy), main="Observed functions")
#' matplot((0:100)/100, t(X.toy), type="l", lty=1, col=y.rainbow, xlab="t",
#'         ylab="", main="Rainbow plot")
#'
#' ## Obtain DTW distances
#' D.dtw <- dist(X.toy, method="dtw", window.type="sakoechiba", window.size=5)
#'
#' ## Compare PC vs. PCo ridge regression
#'
#' # matrix to store results
#' GCVmat <- matrix(NA, 15, 2)
#' # dummy response variable
#' dummy <- rep(1,30)
#'
#' # loop over possible projection dimensions
#' for (k. in 1:15){
#'   # fit PC (m1) and PCo (m2) ridge regression
#'   m1 <- gam(y.toy ~ s(dummy, bs="pco", k=k.,
#'             xt=list(realdata=X.toy, dist_fn=dist)), method="REML")
#'   m2 <- gam(y.toy ~ s(dummy, bs="pco", k=k., xt=list(D=D.dtw)), method="REML")
#'   # calculate and store GCV scores
#'   GCVmat[k., ] <- length(y.toy) * c(sum(m1$residuals^2)/m1$df.residual^2,
#'                    sum(m2$residuals^2)/m2$df.residual^2)
#' }
#'
#' ## plot the GCV scores per dimension for each model
#' matplot(GCVmat, lty=1:2, col=1, pch=16:17, type="o", ylab="GCV",
#'         xlab="Number of principal components / coordinates",
#'         main="GCV score")
#' legend("right", c("PC ridge regression", "DTW-based PCoRR"), lty=1:2, pch=16:17)
#'
#' ## example of making a prediction
#'
#' # fit a model to the toy data
#' m <- gam(y.toy ~ s(dummy, bs="pco", k=2, xt=list(D=D.dtw)), method="REML")
#'
#' # first build the distance matrix
#' # in this case we just subsample the original matrix
#' # see ?pco_predict_preprocess for more information on formatting this data
#' dist_list <- list(dummy = as.matrix(D.dtw)[, c(1:5,10:15)])
#'
#' # preprocess the prediction data
#' pred_data <- pco_predict_preprocess(m, newdata=NULL, dist_list)
#'
#' # make the prediction
#' p <- predict(m, pred_data)
#'
#' # check that these are the same as the corresponding fitted values
#' print(cbind(fitted(m)[ c(1:5,10:15)],p))
#'
#' }
smooth.construct.pco.smooth.spec <- function(object, data, knots){

  ## test what we got given
  # all the extra stuff gets put in object$xt
  xt <- object$xt
  if(all(c("D","realdata","dist_fn") %in% names(xt)) |
     all(c("D","dist_fn") %in% names(xt)) |
     all(c("D","realdata") %in% names(xt)) |
     !any(c("D", "realdata", "dist_fn") %in% names(xt))){
    stop("Please supply either a distance matrix or data and distance function!")
  }

  # distance matrix
  if(all(c("realdata","dist_fn") %in% names(xt))){
    D <- xt$dist_fn(xt$realdata)
  }else{
    D <- xt$D
  }

  # projection dimension
  pdim <- object$bs.dim

  ## do some input checking
  # either K or D must be supplied
  if(is.null(D)) {
    stop("No distance matrix!")
  }
  # default to use regular cmdscale
  if(is.null(xt$fastcmd)){
    xt$fastcmd <- FALSE
  }

  # use the additive constant?
  # default to FALSE
  if(!is.null(xt$add)){
    add <- xt$add
  }else{
    add <- FALSE
  }

  # what do cmdscale options mean?!
  # k     - dimension of MDS projection
  # eig   - return eigenvalues
  # x.ret - return (double centered distance matrix)
  # add   - add a constant so D is Euclidean (cov() gives -ve
  #         values, which is not a property of a distance.
  if(xt$fastcmd){
    # use lanczos for the eigendecomposition
    mds.obj <- cmdscale_lanczos(D, k=pdim, eig=TRUE, x.ret=TRUE, add=add)
  }else{
    mds.obj <- cmdscale(D, k=pdim, eig=TRUE, x.ret=TRUE, add=add)
  }

  if(sum(mds.obj$eig>0) < pdim){
    stop("Only the first ",sum(mds.obj$eig>0)," eigenvalues are positive, this is the maximum projection dimension without setting 'add=TRUE', see ?smooth.construct.pco.smooth.spec for further information.")
  }

  ## four required additions to the return object:
  # model matrix
  object$X <- mds.obj$points
  colnames(object$X) <- paste0("pco_", 1:ncol(object$X))
  # penalty matrix
  object$S <- list(diag(nrow = pdim))
  # penalty rank
  object$rank <- array(pdim, dim=c(1,1))
  # null space dimension
  object$null.space.dim <- 0

  # set the label (for plots and summary)
  object$label <- object$term

  # store dimension
  object$dim <- pdim

  ## extra options
  # don't allow tensor products
  object$te.ok <- 0
  # probably best not to plot this with plot.gam
  # (what would the x axis be?)
  object$plot.me <- FALSE

  # see ?smooth.construct for what these mean!
  object$C <- array(dim=c(0, 1))
  object$side.constrain <- FALSE

  object$no.rescale <- TRUE

  # save mds object returned by cmdscale
  object$xt$mds.obj <- mds.obj

  # give it some class
  class(object) <- "pco.smooth"
  return(object)
}
