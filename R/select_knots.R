
select_knots <- function(t,knots=10,p=3,option="equally-spaced"){
  
  qs <- seq(0,1,length=knots+1)
  
  if(option=="equally-spaced"){
    knots <- (max(t)-min(t))*qs + min(t)
  }
  if(option=="quantile"){
    knots <- as.vector(quantile(t,qs))
  }
  
  K <- length(knots)
  knots_left <- 2*knots[1]-knots[p:1+1]
  knots_right <- 2*knots[K] - knots[K-(1:p)]
  
  return(c(knots_left,knots,knots_right))
}