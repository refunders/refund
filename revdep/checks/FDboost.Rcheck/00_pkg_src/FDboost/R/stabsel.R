
#' Stability Selection 
#' 
#' Function for stability selection with functional response. Per default the sampling is done
#' on the level of curves and if the model contains a smooth functional intercept, this intercept 
#' is refittedn in each sampling fold. 
#' 
#' @param x fitted FDboost-object 
#' @param refitSmoothOffset logical, should the offset be refitted in each learning sample? 
#' Defaults to \code{TRUE}. 
#' @param cutoff cutoff between 0.5 and 1. Preferably a value between 0.6 and 0.9 should be used.
#' @param q number of (unique) selected variables (or groups of variables depending on the model) 
#' that are selected on each subsample.
#' @param PFER upper bound for the per-family error rate. This specifies the amount of falsely 
#' selected base-learners, which is tolerated. See details of \code{\link[mboost]{stabsel}}.
#' @param folds a weight matrix with number of rows equal to the number of observations, 
#' see \code{{cvLong}}. Usually one should not change the default here as subsampling 
#' with a fraction of 1/2 is needed for the error bounds to hold. One usage scenario where 
#' specifying the folds by hand might be the case when one has dependent data (e.g. clusters) and 
#' thus wants to draw clusters (i.e., multiple rows together) not individuals.
#' @param assumption Defines the type of assumptions on the distributions of the selection probabilities 
#' and simultaneous selection probabilities. Only applicable for \code{sampling.type = "SS"}. 
#' For \code{sampling.type = "MB"} we always use \code{"none"}.
#' @param sampling.type use sampling scheme of of Shah & Samworth (2013), i.e., with complementary pairs 
#' (\code{sampling.type = "SS"}), or the original sampling scheme of Meinshausen & Buehlmann (2010).
#' @param B number of subsampling replicates. Per default, we use 50 complementary pairs for the error 
#' bounds of Shah & Samworth (2013) and 100 for the error bound derived in Meinshausen & Buehlmann (2010). 
#' As we use \code{B} complementary pairs in the former case this leads to \code{2B} subsamples.
#' @param papply (parallel) apply function, defaults to mclapply. Alternatively, parLapply can be used. 
#' In the latter case, usually more setup is needed (see example of cvrisk for some details).
#' @param verbose logical (default: TRUE) that determines wether warnings should be issued.
#' @param eval logical. Determines whether stability selection is evaluated (\code{eval = TRUE}; default) 
#' or if only the parameter combination is returned.
#' @param ... additional arguments to \code{\link[mboost]{cvrisk}} or \code{\link{validateFDboost}}.
#' 
#' @details The number of boosting iterations is an important hyper-parameter of the boosting algorithms 
#' and can be chosen using the functions \code{cvrisk.FDboost} and \code{validateFDboost} as they compute
#' honest, i.e. out-of-bag, estimates of the empirical risk for different numbers of boosting iterations. 
#' The weights (zero weights correspond to test cases) are defined via the folds matrix, 
#' see \code{\link[mboost]{cvrisk}} in package mboost. 
#' See Hofner et al. (2015) for the combination of stability selection and component-wise boosting. 
#' 
#' @seealso \code{\link[mboost]{stabsel}} to perform stability selection for a mboost-object.
#' 
#' @references 
#' B. Hofner, L. Boccuto and M. Goeker (2015), Controlling false discoveries in 
#' high-dimensional situations: boosting with stability selection. 
#' BMC Bioinformatics, 16, 1-17.
#'  
#' N. Meinshausen and P. Buehlmann (2010), Stability selection. 
#' Journal of the Royal Statistical Society, Series B, 72, 417-473.
#' 
#' R.D. Shah and R.J. Samworth (2013), Variable selection with error control: 
#' another look at stability selection. Journal of the Royal Statistical Society, Series B, 75, 55-80. 
#' 
#' @return An object of class \code{stabsel} with a special print method. 
#' For the elements of the object, see \code{\link[mboost]{stabsel}}
#' 
#' @examples
#' ######## Example for function-on-scalar-regression
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
#' ## fit a model cotaining all main effects 
#' modAll <- FDboost(vis ~ 1 
#'           + bolsc(T_C, df=1) %A0% bbs(time, df=5) 
#'           + bolsc(T_A, df=1) %A0% bbs(time, df=5)
#'           + bolsc(T_B, df=1) %A0% bbs(time, df=5)
#'           + bolsc(rspeed, df=1) %A0% bbs(time, df=5)
#'           + bolsc(mflow, df=1) %A0% bbs(time, df=5), 
#'        timeformula = ~bbs(time, df=5), 
#'        numInt = "Riemann", family = QuantReg(), 
#'        offset = NULL, offset_control = o_control(k_min = 10),
#'        data = viscosity, 
#'        control = boost_control(mstop = 100, nu = 0.2))
#' 
#' 
#' ## create folds for stability selection  
#' ## only 5 folds for a fast example, usually use 50 folds 
#' set.seed(1911)
#' folds <- cvLong(modAll$id, weights = rep(1, l = length(modAll$id)), 
#'                 type = "subsampling", B = 5) 
#'     
#' \dontrun{        
#' ## stability selection with refit of the smooth intercept 
#' stabsel_parameters(q = 3, PFER = 1, p = 6, sampling.type = "SS")
#' sel1 <- stabsel(modAll, q = 3, PFER = 1, folds = folds, grid = 1:200, sampling.type = "SS")
#' sel1
#' 
#' ## stability selection without refit of the smooth intercept 
#' sel2 <- stabsel(modAll, refitSmoothOffset = FALSE, q = 3, PFER = 1, 
#'                 folds = folds, grid = 1:200, sampling.type = "SS")
#' sel2
#' }
#' 
#' @export
## code for stabsel.mboost taken from mboost 2.6-0
## stabsel method for FDboost; requires stabs
stabsel.FDboost <- function(x, refitSmoothOffset = TRUE, 
                            cutoff, q, PFER,
                            # folds = subsample(model.weights(x), B = B),
                            folds = cvLong(x$id, weights = rep(1, l = length(x$id)), type = "subsampling", B = B), 
                            B = ifelse(sampling.type == "MB", 100, 50),
                            assumption = c("unimodal", "r-concave", "none"),
                            sampling.type = c("SS", "MB"),
                            papply = mclapply, verbose = TRUE, eval = TRUE, ...) {

  cll <- match.call()
  p <- length(variable.names(x))
  ibase <- 1:p
  
  sampling.type <- match.arg(sampling.type)
  if (sampling.type == "MB")
    assumption <- "none"
  else
    assumption <- match.arg(assumption)
  
  B <- ncol(folds)
  
  pars <- stabsel_parameters(p = p, cutoff = cutoff, q = q,
                             PFER = PFER, B = B,
                             verbose = verbose, sampling.type = sampling.type,
                             assumption = assumption)
  ## return parameter combination only if eval == FALSE
  if (!eval)
    return(pars)
  
  cutoff <- pars$cutoff
  q <- pars$q
  PFER <- pars$PFER
  
  fun <- function(model) {
    xs <- selected(model)
    qq <- sapply(1:length(xs), function(x) length(unique(xs[1:x])))
    xs[qq > q] <- xs[1]
    xs
  }
  if (sampling.type == "SS") {
    ## use complementary pairs
    folds <- cbind(folds, model.weights(x) - folds)
  }
  
  ## for scalar response and/or scalar offset, use the more efficient cvrisk()
  if( any(class(x) == "FDboostScalar" )  ) refitSmoothOffset <- FALSE
  if( !is.null(x$call$offset) && x$call$offset == "scalar" ) refitSmoothOffset <- FALSE

  if(refitSmoothOffset){
    message("Use validateFDboost() to recompute the smooth offset in each fold.")
    ## folds are on level of single observations 
    ## but validateFDboost() expects folds on the level of curves 
    folds <- folds[! duplicated(x$id), ]
    ss <- validateFDboost(x, fun = fun,
                 folds = folds, getCoefCV = FALSE, ...)$fun_ret #, papply = papply 
                  
  }else{
    ss <- cvrisk(x, fun = fun,
                 folds = folds,
                 papply = papply, ...)
  }

  
  if (verbose){
    qq <- sapply(ss, function(x) length(unique(x)))
    sum_of_violations <- sum(qq < q)
    if (sum_of_violations > 0)
      warning(sQuote("mstop"), " too small in ",
              sum_of_violations, " of the ", ncol(folds),
              " subsampling replicates to select ", sQuote("q"),
              " base-learners; Increase ", sQuote("mstop"),
              " bevor applying ", sQuote("stabsel"))
  }
  
  
  ## if grid specified in '...'
  if (length(list(...)) >= 1 && "grid" %in% names(list(...))) {
    m <- max(list(...)$grid)
  } else {
    m <- mstop(x)
  }
  ret <- matrix(0, nrow = length(ibase), ncol = m)
  for (i in 1:length(ss)) {
    tmp <- sapply(ibase, function(x)
      ifelse(x %in% ss[[i]], which(ss[[i]] == x)[1], m + 1))
    ret <- ret + t(sapply(tmp, function(x) c(rep(0, x - 1), rep(1, m - x + 1))))
  }
  
  phat <- ret / length(ss)
  rownames(phat) <- names(variable.names(x))
  if (extends(class(x), "glmboost"))
    rownames(phat) <- variable.names(x)
  ret <- list(phat = phat, selected = which((mm <- apply(phat, 1, max)) >= cutoff),
              max = mm, cutoff = cutoff, q = q, PFER = PFER, p = p, B = B, 
              sampling.type = sampling.type, assumption = assumption,
              call = cll)
  ret$call[[1]] <- as.name("stabsel")
  class(ret) <- c("stabsel_FDboost", "stabsel_mboost", "stabsel")
  ret
}


