#' Random effect parser
#' 
#' Internal function used to parse random effects in the formula interface
#' to other functions.
#' 
#' @param variable id used to identify clusters
#' 
#' @author Julia Wrobel \email{jw3134@@cumc.columbia.edu}
#' 
re.fosr <-function(variable){
  ret <- as.list(rep(NA, 2))
  ret[[1]] <- deparse(substitute(variable))
  ret[[2]] <- (paste("(1|",deparse(substitute(variable)), ")", sep=""))
  return(ret)
}
