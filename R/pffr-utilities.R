# Utility functions for pffr:

safeDeparse <- function(expr){
	# turn an expression into a _single_ string, regardless of the expression's length
	ret <- paste(deparse(expr), collapse="")
	#rm whitespace
	gsub("[[:space:]][[:space:]]+", " ", ret)
}

#' Return call with all possible arguments
#'
#' Return a call in which all of the arguments which were supplied or have presets are specified by their full names and their supplied or default values.
#'
#' @param definition a function. See \code{\link[base]{match.call}}.
#' @param call an unevaluated call to the function specified by definition. See \code{\link[base]{match.call}}.
#' @param expand.dots logical. Should arguments matching ... in the call be included or left as a ... argument? See \code{\link[base]{match.call}}.
#' @return An object of mode "\code{\link[base]{call}}".
#' @author Fabian Scheipl
#' @seealso \code{\link[base]{match.call}}
expand.call <- function(definition=NULL, call=sys.call(sys.parent(1)),
  expand.dots = TRUE)
  {
    call <- match.call(definition, call, expand.dots)
    #given args:
    ans <- as.list(call)

    this_function <- safeDeparse(ans[[1]])
    # rm namespace operators so that function finding works (mlr-org/mlr/#1559)
    if (grepl(":", this_function)) {
      this_function <- gsub("(^[^:]+:+)(.+$)", "\\2", this_function)
      if (!exists(this_function)) {
        warning("expand.call couldn't find ", this_function,
          " and returned unchanged call:\n <", call, ">")
        return(call)
      }
    }
    frmls <- formals(this_function)
    #remove formal args with no presets:
    frmls <- frmls[!sapply(frmls, is.symbol)]
    add <- which(!(names(frmls) %in% names(ans)))
    return(as.call(c(ans, frmls[add])))
  }



list2df <- function(l){
# make a list into a dataframe -- matrices are left as matrices!
    nrows <- sapply(l, function(x) nrow(as.matrix(x)))
    stopifnot(length(unique(nrows)) == 1)
    ret <- data.frame(rep(NA, nrows[1]))
    for(i in 1:length(l)) ret[[i]] <- l[[i]]
    names(ret) <- names(l)
    return(ret)
}

## TODO: this does not always yield unique labels, e.g. if you have
##   s(g, bs="re") + s(g, bs="mrf", xt=somepenalty)
getShrtlbls <- function(object){
# make short labels for display/coef-list, etc...

    labelmap <- object$pffr$labelmap

    ret <- sapply(names(unlist(labelmap)),
            function(x){
                #make a parseable expression for ffpc terms
                if(grepl("^ffpc", x)){
                    ffpcnumber <- gsub("(^.+))([0-9]+$)","\\2", x)
                    x <- gsub(")[0-9]+",")",x)
                   }
                exprx <- parse(text=x)

                #remove argument names
                x <- gsub("((?:[A-Za-z]+))(\\s)(=)(\\s)", "", x, perl=TRUE)

                #remove whitespace
                x <- gsub("\\s", "", x)

                #remove everything after and including first quoted argument
                if(any(regexpr(",\".*",x)[[1]]>0)) {
                    x <- gsub("([^\"]*)(,[c\\(]*\".*)(\\)$)", "\\1\\3", x, perl=TRUE)
                }

                #remove everything after last variable:
                lstvrbl <- tail(all.vars(exprx),1)
                x <- gsub(paste("\\(.*(^.*?(?=",lstvrbl,")",lstvrbl,")(.*$)",sep=""), "\\1", x, perl=TRUE)

                #match braces
                openbr <- sum(grepl("\\(", strsplit(x, "")[[1]]))
                closebr <- sum(grepl("\\)", strsplit(x, "")[[1]]))
                if(openbr>closebr) x <- paste(c(x, rep(")", openbr-closebr)), sep="",collapse="")

                #add number of PC for ffpc terms
                if(grepl("^ffpc", x)){
                    x <- paste(x, ffpcnumber, sep="")
                }
                return(x)
            })
    # correct labels for factor variables:
    if(any(sapply(labelmap, length) > 1 & !sapply(names(labelmap), function(x) grepl("^ffpc", x)))){
        which <- which(sapply(labelmap, length) > 1 & !sapply(names(labelmap), function(x) grepl("^ffpc", x)))
        inds <- c(0, cumsum(sapply(labelmap, length)))
        for(w in which){
            ret[(inds[w]+1):inds[w+1]] <- {
                lbls <- labelmap[[w]]
                bylevels <- sapply(object$smooth[lbls], function(x) x$by.level)
                by <- object$smooth[[lbls[1]]]$by
                paste(by, bylevels, "(", object$pffr$yindname, ")", sep="")
            }
        }
    }
    #append labels for varying coefficient terms
    if(any(!grepl("\\(", ret))){
        which <- which(!grepl("\\(", ret))
        ret[which] <- paste(ret[which],"(", object$pffr$yindname, ")", sep="")
    }

    return(ret)
}

#' Simulate example data for pffr
#'
#' Simulates example data for \code{\link{pffr}} from a variety of terms.
#' Scenario "all" generates data from a complex multivariate model \deqn{Y_i(t)
#' = \mu(t) + \int X_{1i}(s)\beta_1(s,t)ds + xlin \beta_3(t) + f(xte1, xte2) +
#' f(xsmoo, t) + \beta_4 xconst + f(xfactor, t) +  \epsilon_i(t)}. Scenarios "int", "ff", "lin",
#' "te", "smoo", "const", "factor", generate data from simpler models containing only the
#' respective term(s)  in the model equation given above. Specifying a
#' vector-valued scenario will generate data from a combination of the
#' respective terms. Sparse/irregular response trajectories can be generated by
#' setting \code{propmissing} to something greater than 0 (and smaller than 1).
#' The return object then also includes a \code{ydata}-item with the sparsified
#' data.
#'
#' See source code for details.\cr
#'
#' @param scenario see Description
#' @param n number of observations
#' @param nxgrid number of evaluation points of functional covariates
#' @param nygrid number of evaluation points of the functional response
#' @param SNR the signal-to-noise ratio for the generated data: empirical
#'   variance of the additive predictor divided by variance of the errors.
#' @param propmissing proportion of missing data in the response, default = 0.
#'   See Details.
#' @param limits a function that defines an integration range, see
#'   \code{\link{ff}}
#' @export
#' @importFrom splines spline.des
#' @importFrom stats dnorm rnorm
#' @return a named list with the simulated data, and the true components of the
#'   predictor etc as attributes.
pffrSim <- function(
        scenario="all",
        n = 100,
        nxgrid = 40,
        nygrid = 60,
        SNR = 10,
        propmissing=0,
        limits = NULL){

   mc <- match.call()
   for(i in 2:length(mc)) if(is.symbol(mc[[i]]))
     mc[[i]] <- get(deparse(mc[[i]]), envir=parent.frame())

    ## generates random functions...
    rf <- function(x=seq(0,1,length=100), bs.dim=7, center=FALSE) {
      nk <- bs.dim - 2
      xu <- max(x)
      xl <- min(x)
      xr <- xu - xl
      xl <- xl - xr * 0.001
      xu <- xu + xr * 0.001
      dx <- (xu - xl)/(nk - 1)
      kn <- seq(xl - dx * 3, xu + dx * 3,
          length = nk + 4 + 2)
      X <- splines::spline.des(kn, x, 4, x * 0)$design

      drop(X %*% rnorm(bs.dim))
    }

    test1 <- function(s, t){
      s*cos(pi*abs(s-t)) - .19
    }
#     test2 <- function(s, t, ss=0.3, st=0.4)
#     {
#       cos(pi*s)*sin(pi*t) + (s*t)^2 - 0.11
#     }
     s <- seq(0, 1, length=nxgrid)
     t <- seq(0, 1, length=nygrid)


    mu.t <- matrix(1 + dbeta(t, 2,7), nrow=n, ncol=nygrid, byrow=TRUE)

    data <- list()
    etaTerms <- list()
    etaTerms$int <- mu.t

    #functional covariates
    data$X1 <- I(t(replicate(n, rf(s))))

    L <- matrix(1/nxgrid, ncol=nxgrid, nrow=n)
    LX1 <- L*data$X1
    beta1.st <- outer(s, t, test1)
    if(!is.null(limits)){
       range <- outer(s, t, limits)
       beta1.st  <- beta1.st * range
    }
    etaTerms$X1 <- LX1 %*% beta1.st

#     data$X2 <- I(t(replicate(n, rf(s))))
#     LX2 <- L*data$X2
#     beta2.st <- outer(s, t, test2)
#     etaTerms$X2 <- LX2%*%beta2.st

    #scalar covariates
    data$xlin <- I(rnorm(n))
    beta.t <- matrix(scale(-dnorm(4 * (t - .2))), nrow = n, ncol = nygrid,
      byrow = T)
    etaTerms$xlin <- data$xlin * beta.t

    data$xsmoo <- I(rnorm(n))
    etaTerms$xsmoo <- outer(drop(scale(cos(data$xsmoo))), (t-.5), "*")

    data$xfactor <- sample(gl(3, n/3), replace = TRUE)
    etaTerms$xfactor <-  2* as.numeric(data$xfactor) +
      sin(2 * outer(as.numeric(data$xfactor), t))
    if ("2factor" %in% scenario) {
      data$x2factor <- sample(gl(3, n/3), replace = TRUE)
      etaTerms$x2factor <-  2 * as.numeric(data$x2factor) +
        cos(2 * outer(as.numeric(data$x2factor), t))
    }

    data$xte1 <- I(rnorm(n))
    data$xte2 <- I(rnorm(n))
    etaTerms$xte <- matrix(drop(scale(-data$xte1*data$xte2^2)),
                           ncol=nygrid,
                           nrow=n)

    data$xconst <- I(rnorm(n))
    etaTerms$xconst <- matrix(2*data$xconst,
                              ncol=nygrid,
                              nrow=n)

    if(length(scenario)==1){
      eta <- mu.t + switch(scenario,
                           "int" = 0,
                           "all" = Reduce("+", etaTerms),
                           "ff" =  etaTerms$X1, #+ etaTerms$X2,
                           "lin" = etaTerms$xlin,
                           "smoo" = etaTerms$xsmoo,
                           "te" = etaTerms$xte,
                           "const" = etaTerms$xconst,
                           "factor" = etaTerms$xfactor)
    } else {
      stopifnot(all(scenario %in% c("int" ,"ff", "lin", "smoo", "te", "const",
        "factor", "2factor")))
      eta <- 0*mu.t
      if("int" %in% scenario) eta <- eta + mu.t
      if("ff" %in% scenario) eta <- eta + etaTerms$X1 #+ etaTerms$X2
      if("lin" %in% scenario) eta <- eta + etaTerms$xlin
      if("smoo" %in% scenario) eta <- eta + etaTerms$xsmoo
      if("te" %in% scenario) eta <- eta + etaTerms$xte
      if("const" %in% scenario) eta <- eta + etaTerms$xconst
      if("factor" %in% scenario) eta <- eta + etaTerms$xfactor
      if("2factor" %in% scenario) eta <- eta + etaTerms$x2factor
    }


    eps <-  sd(as.vector(eta))/sqrt(SNR) * matrix(scale(rnorm(n*nygrid)),
                                                  nrow=n)
    data$Y <- I(eta + eps)



    if(propmissing == 0){
      return(structure(as.data.frame(data, rownames=1:n), xindex=s, yindex=t,
                       truth=list(eta=eta, etaTerms=etaTerms), call=mc))
    } else {
      missing <- sample(c(rep(T, propmissing*n*nygrid),
                          rep(F, n*nygrid-propmissing*n*nygrid)))
      data <- as.data.frame(data, rownames=1:n)

      ydata <- data.frame(.obs = rep(1:n, each=nygrid)[!missing],
                          .index = rep(t, times=n)[!missing],
                          .value = as.vector(t(data$Y))[!missing])

      return(structure(list(data=data, ydata=ydata), xindex=s, yindex=t,
                       truth=list(eta=eta, etaTerms=etaTerms),
                       call=mc))
    }
}



#' P-spline constructor with modified 'shrinkage' penalty
#'
#' Construct a B-spline basis with a modified difference penalty
#' of full rank (i.e., that also penalizes low-order polynomials).
#'
#' This penalty-basis combination is useful to avoid non-identifiability issues for \code{\link{ff}} terms.
#' See 'ts' or 'cs' in \code{\link[mgcv]{smooth.terms}}
#' for similar "shrinkage penalties" for thin plate and cubic regression splines.
#' The basic idea is to replace the k-th zero eigenvalue of the original penalty by
#' \eqn{s^k \nu_m}, where \eqn{s} is the shrinkage factor (defaults to 0.1)
#' and \eqn{\nu_m} is the smallest non-zero eigenvalue. See reference for the
#' original idea, implementation follows that in the 'ts' and 'cs' constructors
#' (see \code{\link[mgcv]{smooth.terms}}).
#'
#' @param object see \code{\link[mgcv]{smooth.construct}}. The shrinkage factor can be speficied via \code{object$xt$shrink}
#' @param data see \code{\link[mgcv]{smooth.construct}}.
#' @param knots see \code{\link[mgcv]{smooth.construct}}.
#'
#' @author Fabian Scheipl;  adapted from 'ts' and 'cs' constructors by S.N. Wood.
#'
#' @references Marra, G., & Wood, S. N. (2011). Practical variable selection for generalized additive models.
#' \emph{Computational Statistics & Data Analysis}, 55(7), 2372-2387.
#' @method smooth.construct pss.smooth.spec
#' @export
#' @importFrom mgcv smooth.construct.ps.smooth.spec
smooth.construct.pss.smooth.spec<-function(object,data,knots)
{

  shrink <- ifelse(is.null(object$xt$shrink), 0.1, object$xt$shrink)
  stopifnot(shrink>0, shrink<1)

  object <- smooth.construct.ps.smooth.spec(object,data,knots)
  nk <- object$bs.dim
  difforder <- object$m[2]
  ## add shrinkage term to penalty:
  ## Modify the penalty by increasing the penalty on the
  ## unpenalized space from zero...
  es <- eigen(object$S[[1]],symmetric=TRUE)
  ## now add a penalty on the penalty null space
  es$values[(nk-difforder+1):nk] <- es$values[nk-difforder]*shrink
  ## ... so penalty on null space is still less than that on range space.
  object$S[[1]] <- es$vectors%*%(as.numeric(es$values)*t(es$vectors))
  object$rank <- nk
  object$null.space.dim <- 0
  class(object) <- "pss.smooth"
  object
}

#' @importFrom mgcv Predict.matrix.pspline.smooth
#' @export
Predict.matrix.pss.smooth<-function(object,data)
{
  Predict.matrix.pspline.smooth(object,data)
}

#' @importFrom mgcv smooth.construct.ps.smooth.spec
#' @export
smooth.construct.ps_c.smooth.spec<-function(object,data,knots) {
  object <- smooth.construct.ps.smooth.spec(object,data,knots)
  object$C <- object$xt$C1
  object
}

#' @importFrom mgcv smooth.construct.ps.smooth.spec
#' @export
smooth.construct.fame.smooth.spec<-function(object,data,knots){
  object$bs <- "ps"
  object <- smooth.construct.ps.smooth.spec(object,data,knots)
  # anti-kernel penalty:
  object$S[[1]] <- object$xt$penK

  object$rank <- qr(object$S[[1]])$rank
  object$null.space.dim <-  object$bs.dim - object$rank

  return(object)
}

# compute dim of overlap of the span of two orthonormal matrices
trace_lv <- function(A, B, tol=1e-10){
  ## A, B orthnormal!!
  #Rolf Larsson, Mattias Villani (2001)
  #"A distance measure between cointegration spaces"

  if(NCOL(A)==0 | NCOL(B)==0){
    return(0)
  }

  if(NROW(A) != NROW(B) | NCOL(A) > NROW(A) | NCOL(B) > NROW(B)){
    return(NA)
  }

  trace <- if(NCOL(B)<=NCOL(A)){
    sum(diag(t(B) %*% A %*% t(A) %*% B))
  } else {
    sum(diag(t(A) %*% B %*% t(B) %*% A))
  }
  trace
}


