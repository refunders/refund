#' Summary for a pfr fit
#'
#' Take a fitted \code{pfr}-object and produce summaries from it.
#' See \code{\link[mgcv]{summary.gam}()} for details.
#'
#' @param object a fitted \code{pfr}-object
#' @param ... see \code{\link[mgcv]{summary.gam}()} for options.
#'
#' @return A list with summary information, see \code{\link[mgcv]{summary.gam}()}
#' @method summary pfr
#'
#' @details
#' This function currently simply strips the \code{"pfr"} class label and
#' calls \code{\link[mgcv]{summary.gam}}.
#'
#' @author Jonathan Gellar \email{JGellar@@mathematica-mpr.com}
#' @export
summary.pfr <- function(object, ...) {
  class(object) <- class(object)[-1]
  summary(object)
}
