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
#' @author Jonathan Gellar \email{JGellar@@mathematica-mpr.com}, Fabian Scheipl
#' @export
summary.pfr <- function (object, ...) {
  call <- match.call()
  call[[1]] <- mgcv::summary.gam
  ## drop "pfr" class and replace <object> with changed value s.t. method dispatch
  ## works without glitches
  ## if we don't do this, summary.gam will call coef.pfr on the object and that
  ## doesn't work
  class(object) <- class(object)[!(class(object) %in% "pfr")]
  call$object <- as.name("object")
  eval(call)
  #TODO: modify term names to correspond to pfr formula notation, see summary.pffr
}
