#' Compute quadrature weights
#'
#' Utility function for numerical integration.
#' @param argvals function arguments.
#' @param method quadrature method. Can be either \code{trapedoidal} or \code{midpoint}.
#' @return a vector of quadrature weights for the points supplied in \code{argvals}.
#' @author Clara Happ, with modifications by Philip Reiss 

# TODO: check about 'midpoint'
# TODO: have this function called by lf, af, etc.
# TODO: harmonize with quadrature implemented elsewhere: 'simpson' and 'riemann' options?

quadWeights<- function(argvals, method = "trapezoidal")
{
  ret <- switch(method,
                trapezoidal = {D <- length(argvals)
                               1/2*c(argvals[2] - argvals[1], argvals[3:D] -argvals[1:(D-2)], argvals[D] - argvals[D-1])},
                midpoint = c(0,diff(argvals)),  # why is this called 'midpoint'???
                stop("function quadWeights: choose either trapezoidal or midpoint quadrature rule"))
  
  return(ret)  
}
