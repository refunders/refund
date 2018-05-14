construct.knots <- function(argvals,knots,knots.option,p){
  
if(length(knots)==1){
  allknots <- select.knots(argvals,knots,p=p,option=knots.option)
}

if(length(knots)>1){
  K = length(knots)-1 
  knots_left <- 2*knots[1]-knots[p:1+1]
  knots_right <- 2*knots[K] - knots[K-(1:p)]
  if(p>0) allknots <- c(knots_left,knots,knots_right)
  if(p==0) allknots <- knots
}

return(allknots)

}
