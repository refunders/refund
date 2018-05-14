select.knots <- function(t,knots=10,p=3,option="equally-spaced"){
  
  qs <- seq(0,1,length=knots+1)
  max_t <- max(t) 
  min_t <- min(t) 
  range_t <- max_t - min_t 
  
  min_t <- min_t - range_t * 0.001 # follow Simon Wood's implementation
  max_t <- max_t + range_t * 0.001 #
  
  if(option=="equally-spaced"){
    knots <- (max_t-min_t)*qs + min_t 
  }
  if(option=="quantile"){
    loc_max <- which.max(t)
    loc_min <- which.min(t)
    s <- t+ rnorm(length(t))*range_t/100
    s[loc_max] <- t[loc_max]
    s[loc_min] <- t[loc_min]
    t <- s
    knots <- as.vector(quantile(t,qs))
  }
  
  K <- length(knots)
  knots_left <- 2*knots[1]-knots[p:1+1]
  knots_right <- 2*knots[K] - knots[K-(1:p)]
  
  if(p>0) return(c(knots_left,knots,knots_right))
  if(p==0) return(knots)
}
