#' Random effects constructor for fgam
#'
#' Sets up a random effect for the levels of \code{x}.
#' Use the \code{by}-argument to request random slopes.
#'
#' See \code{\link[mgcv]{random.effects}} in \pkg{mgcv}.
#'
#' @param x a grouping variable: must be a \code{factor}
#' @param ... further arguments handed over to \code{\link[mgcv]{s}},
#' see \code{\link[mgcv]{random.effects}}
#' @return A list with components \code{call} (a rewritten call to
#'   \code{\link[mgcv]{s}} with \code{bs = "re"}) and \code{data} (a named
#'   list containing the grouping factor).
#' @seealso \code{\link[mgcv]{random.effects}}
#' @export
re <- function(x, ...) {
  # TODO: add `cov`-arg, then call bs="mrf" to allow for correlated effects.
  data <- list(x)
  xsymbol <- substitute(x)
  names(data) <- deparse(xsymbol)
  call <- match.call()
  call[[1]] <- quote(s)
  call$bs <- "re"
  list(call = call, data = data)
}
