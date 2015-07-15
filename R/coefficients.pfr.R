#' Extract coefficient functions from a fitted pfr-object
#' 
#' This function is used to extract a coefficient from a fitted `pfr` model, in
#' particular smooth functions resulting from including functional terms specified
#' with \code{lf}, \code{af}, etc. It can also be used to extract smooths
#' genereated using \code{mgcv}'s \code{s}, \code{te}, or \code{t2}.
#' 
#' @param object return object from \code{\link{pfr}}
#' @param select integer indicating the index of the desired smooth term
#'   in \code{object$smooth}. Enter 0 to request the raw coefficients
#'   (i.e., \code{object$coefficients}) and standard errors (if \code{se==TRUE}).
#' @param coords named list indicating the desired coordinates where the
#'   coefficient function is to be evaluated. Names must match the argument names
#'   in \code{object$smooth[[select]]$term}. If \code{NULL}, uses \code{n}
#'   to generate equally-spaced coordinates.
#' @param n integer vector indicating the number of equally spaced coordinates
#'   for each argument. If length 1, the same number is used for each argument.
#'   Otherwise, the length must match \code{object$smooth[[select]]$dim}.
#' @param se if \code{TRUE}, returns pointwise standard error estimates. Defaults
#'   to \code{FALSE} if raw coefficients are being returned; otherwise \code{TRUE}.
#' @param seWithMean if \code{TRUE} the standard errors include uncertainty about
#'   the overall mean; if \code{FALSE}, they relate purely to the centered
#'   smooth itself. Marra and Wood (2012) suggests that \code{TRUE} results in
#'   better coverage performance for GAMs.
#' @param useVc if \code{TRUE}, standard errors are calculated using a covariance
#'   matrix that has been corrected for smoothing parameter uncertainty. This
#'   matrix will only be available under ML or REML smoothing.
# @param limit if \code{TRUE}, checks if a \code{limits} function was
#   used to generate the term, and if so, applies the function to the
#   output to only produce estimates over regions that were within "limits".
#' @param ... these arguments are ignored
#'   
#' @details
#' Should add a description here about entering coordinates.
#'   
#' @return a data frame containing the evaluation points, 
#'    coefficient function values and optionally the SE's for the term indicated
#'    by \code{select}.
#' 
#' @author Fabian Scheipl and Jonathan Gellar
#' 
#' @references
#' Marra, G and S.N. Wood (2012) Coverage Properties of Confidence Intervals for
#' Generalized Additive Model Components. Scandinavian Journal of Statistics.
#' 
#' @importFrom stats coefficients
#' @export

coefficients.pfr <- function(object, select=1, coords=NULL, n=NULL,
                             se=ifelse(length(object$smooth) & select, TRUE, FALSE),
                             seWithMean=FALSE, useVc=TRUE, ...) {
  
  V <- if (useVc & !is.null(object$Vc)) object$Vc else object$Vp
  
  if (!length(object$smooth) | select==0) {
    ret <- if (se) list(coefficients=object$coefficients, se = sqrt(diag(V)))
    else object$coefficients
    return(ret)
  }
  
  object$coefficients[is.na(object$coefficients)] <- 0 #Zero out singular coefs
  smooth.i <- object$smooth[[select]]
  
  coef.i <- if ("random.effect" %in% class(smooth.i)) {
    stop("Random effects not yet implemented for coef.pfr")
    
  } else {
    if (is.null(coords)) {
      # Create new coords
      if (is.null(n))
        n <- sapply(object$model[smooth.i$term], ndefault)
      else if (length(n)==1)
        n <- rep(n, smooth.i$dim)
      else if (length(n)!=smooth.i$dim)
        stop("length of n must match the number of terms of the smooth")
      coords <- mapply(function(x,y) {
        seq(min(x), max(x), length=y)
      }, object$model[smooth.i$term], n, SIMPLIFY=FALSE)
    } else {
      # Check coords to make sure they match, otherwise throw a warning
      cnms <- names(coords)
      if (length(smooth.i$term)==1) {
        if (!is.null(cnms))
          if (!any(is.na(cnms), smooth.i$term==cnms,
                   modify_nm(smooth.i$term)==cnms))
            warning("Re-naming coordinates to those used to generate the term")
        names(coords) <- smooth.i$term
      } else if (is.null(names(coords))) {
        if (length(coords)==1) {
          warning("Using same coordinates for all arguments (dimensions)")
          coords <- rep(coords, length(smooth.i$term))
        } else if (length(coords)==length(smooth.i$term)) {
          warning("Assigning coordinate names from smooth.i$term")
        } else
          stop("coords must be a named list with names in smooth[[select]]$term")
        names(coords) <- smooth.i$term
      } else {
        # Check if user used the modified coordinate name
        for (i in 1:length(coords))
          for (j in 1:length(smooth.i$term))
            if (cnms[i] == modify_nm(smooth.i$term[j]))
              # Rename coordinate with unmodified name
              names(coords)[i] <- smooth.i$term[j]
        if (!all(smooth.i$term %in% names(coords))) {
          newtrms <- smooth.i$term[!(smooth.i$term %in% names(coords))]
          if (is.null(n))
            n <- sapply(object$model[newtrms], ndefault)
          else if (length(n)==1)
            n <- rep(n, length(newtrms))
          else
            stop("length(n) must match the number of arguments not included in coords")
          newcoords <- mapply(function(x,y) {
            seq(min(x), max(x), length=y)
          }, object$model[newtrms], n, SIMPLIFY=FALSE)
          coords <- c(coords, newcoords)
        }
      }
    }
    
    # Create grid and get prediction matrix
    coef.i <- do.call(expand.grid, coords)[smooth.i$term]
    if (smooth.i$by!="NA")
      coef.i[smooth.i$by] <- 1
    pmat <- mgcv::PredictMat(smooth.i, coef.i)
    
    # Get coefficient function
    coef.i[smooth.i$by] <- NULL
    first <- smooth.i$first.para
    last  <- smooth.i$last.para
    coef.i$value <- pmat %*% object$coefficients[first:last]
    
    #### GET SE: COPIED FROM plot.mgcv.smooth
    if (se) { ## get standard errors for fit
      ## test whether mean variability to be added to variability (only for centred terms)
      if (seWithMean && attr(smooth.i, "nCons")>0) {
        if (length(object$cmX) < ncol(V))
          object$cmX <- c(object$cmX, rep(0,ncol(V)-length(object$cmX)))
        
        X1 <- matrix(object$cmX, nrow(coef.i), ncol(V), byrow=TRUE)
        meanL1 <- smooth.i$meanL1
        if (!is.null(meanL1)) {
          X1 <- X1 / meanL1
        } 
        X1[,first:last] <- pmat
        coef.i$se  <- sqrt(pmax(0, rowSums((X1 %*% V) * X1)))
      } else {
        ## se in centred (or anyway unconstained) space only
        coef.i$se  <-
          sqrt(pmax(0, rowSums((pmat %*% V[first:last, first:last, drop=FALSE])
                               * pmat)))
      }
      
    }
    
    ###### LIMITS NOT YET IMPLEMENTED
    ## Apply limits if requested
    #if (limit) {
    #  # Check if limits function exists
    #  tf.env <- environment(object$pfr$t.funcs[[
    #    which(sapply(object$pfr$smoothmap, function(x) select %in% x))
    #    ]])
    #  limits.i <- tf.env$limits
    #  if (is.function(limits.i)) {
    #    args <- coef.i[1:length(formals(limits.i))]
    #    names(args) <- names(formals(limits.i))
    #    coef.i <- coef.i[do.call(limits.i, args), ]
    #  }
    #}
    coef.i
  }
  
  # Modify coordinate names
  names(coef.i)[1:smooth.i$dim] <- sapply(smooth.i$term, modify_nm)
  coef.i
}


#' @rdname coefficients.pfr
#' @export
coef.pfr <- coefficients.pfr

modify_nm <- function(x) {
  x <- gsub("\\.omat", "", x)
  x <- gsub("tmat", "argvals", x)
  x
}

ndefault <- function(x) {
  x <- unique(as.vector(x))
  x <- x[order(x)]
  
  diffs <- unique(round(diff(x), 10))
  n <- if (length(diffs)<5) {
    length(seq(min(x), max(x), by=min(diffs)))
  } else 101
  while(n>101)
    n <- round(n/2)
  n
}
