#' Cross-validate a pfr model
#' 
#' This function performs cross-validation for a pfr model. This would most
#' commonly be used to optimize a user-specified tuning parameter, such as
#' the \code{k} value for a basis, or the \code{ncomp} argument of \code{fpc}.
#' It can also be used to try a list of different options identified by a
#' character string, such as \code{bs} or \code{integration}.
#' 
#' @param object Either a fitted \code{pfr} model, or the formula to fit a
#'   new \code{pfr} model
#' @param xval a list that defines which parameters to cross-validate; see
#'   Details below.
#' @param K the number of folds for K-fold cross-validation; enter
#'   \code{nrow(data)} for N-fold cross-validations
#' @param ... Other arguments to pass on to other functions
#' 
#' @details
#' These are the details that need to be filled in
#' 
#' @author Jonathan Gellar
#' 
#' @examples
#' xval <- list(list(term="X1", nbasis=5:10), list(term="X2", bs=c("tp", "bs")))
#' xval <- list(list(term="cca",  bs=c("ps", "tp", "cr"), k=10:12),
#'              list(term="rcst", bs=c("ps", "tp", "cr"), k=13:15))
#' 




pfr.xval <- function(formula, xval, data, K=10, ...) {
  
  #form <- if (class(object)=="formula") {
  #  object
  #} else if ("pfr" %in% class(object)) {
  #  object$formula
  #} else stop("Unrecognized type of object")
  
  tf <- terms.formula(formula, specials = c("s", "te", "t2", "lf", "af",
                                            "lf.vd", "re", "peer", "fpc"))
  responsename <- attr(tf, "variables")[2][[1]]
  response <- data[[responsename]]
  
  # Turn xval into a more convenient list: one element per parameter
  nms <- do.call("c", lapply(xval, function(xval.i) {
    paste(xval.i$term, names(xval.i)[-1], sep=".")
  }))
  xval <- do.call("c", lapply(xval, function(xval.i) xval.i[-1]))
  names(xval) <- nms
  
  xval.grid <- expand.grid(xval)
  folds <- rep(1:K, ceiling(N/K))[1:N]
  xval.res <- apply(xval.grid[1:3,], 1, function(xval.i) {
    print(xval.i)
    newfrml <- modifyForm(form, xval.i)
    folds <- sample(folds)
    preds <- rep(NA, N)
    for (k in 1:K) {
      cat(".")
      fit.k <- pfr(newfrml, data=data[folds!=k,], ...)
      preds[folds==k] <- predict(fit.k, newdata = data[folds==k,,drop=FALSE])
    }
    mean((preds-response)^2)
    cat("\n")
  })
  
}


