#' Prediction from a fitted bayes_fosr model
#'
#' Takes a fitted \code{fosr}-object produced by \code{\link{bayes_fosr}} and produces predictions given a
#' new set of values for the model covariates or the original values used for the model fit.
#' 
#' @param object a fitted \code{fosr} object as produced by \code{\link{bayes_fosr}}
#' @param newdata a named list containing the values of the model covariates at which predictions
#' are required. If this is not provided then predictions corresponding to the original data are
#' returned. All variables provided to newdata should be in the format supplied to the model fitting 
#' function.
#' @param ... additional (unused) arguments
#' 
#' @return ...
#' 
#' @author Jeff Goldsmith \email{jeff.goldsmith@@columbia.edu}
#' @seealso \code{\link{bayes_fosr}}
#' @export
#' @examples
#' \dontrun{
#' library(reshape2)
#' library(dplyr)
#' library(ggplot2)
#' 
#' ##### Cross-sectional real-data example #####
#' 
#' ## organize data
#' data(DTI)
#' DTI = subset(DTI, select = c(cca, case, pasat))
#' DTI = DTI[complete.cases(DTI),]
#' DTI$gender = factor(sample(c("male","female"), dim(DTI)[1], replace = TRUE))
#' DTI$status = factor(sample(c("RRMS", "SPMS", "PPMS"), dim(DTI)[1], replace = TRUE))
#' 
#' ## fit models
#' VB = bayes_fosr(cca ~ pasat, data = DTI, Kp = 4, Kt = 10)
#' 
#' ## obtain predictions
#' pred = predict(VB, sample_n(DTI, 10))
#' }
#' 
predict.fosr <- function (object, newdata, ...) {
  
  if (!missing(newdata)) {
    X.design = model.matrix(object$terms, newdata)
    y = X.design %*% object$beta.hat
  }
  
  else {
    y = object$Yhat
  }
  
  return(y)
  
}

