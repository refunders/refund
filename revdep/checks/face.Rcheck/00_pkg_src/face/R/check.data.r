check.data <- function(data,type="fit"){
  if(type=="fit") {
    if(!is.data.frame(data)|is.null(data$y)|is.null(data$subj)|is.null(data$argvals))
      stop("'data' should be a data frame with three variables:argvals,subj and y")
    
    if(sum(is.na(data))>0) stop("No NA values are allowed in the data")
  
  }
  
  if(type=="predict"){
    if(!is.data.frame(data)|is.null(data$y)|is.null(data$subj)|is.null(data$argvals))
      stop("'newdata' should be a data frame with three variables:argvals,subj and y") 
  }
  return(0)
}