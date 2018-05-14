predict.pspline.face <- function(object,argvals.new,...){
  
  #require(splines)
  stopifnot(class(object)=="pspline.face")
  
  knots <- object$knots
  p <- object$p
  B = spline.des(knots=knots, x = argvals.new, ord = p+1,outer.ok = TRUE,sparse=TRUE)$design
  return(as.vector(B%*%object$theta))
  
}