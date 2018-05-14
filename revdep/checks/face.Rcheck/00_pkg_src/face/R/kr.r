kr <- function (A, B, w, byrow = TRUE) 
{
    if (byrow) {
    if (nrow(A) != nrow(B)) 
      stop("Dimensions of the matrices do not match.")
    if (missing(w)) 
      w <- rep(1, nrow(A))
    if (nrow(A) != length(w)) 
      stop("Length of the weight does not match with the dimension of the matrices.")
    cola <- ncol(A)
    colb <- ncol(B)
    colab <- cola * colb
    
    expr <- paste("rbind(", paste(rep("A", colb), collapse = ","), 
                             ")", sep = "")
    A <- eval(parse(text = expr))
    A <- matrix(c(A), nrow(B), ncol = colab)
    A <- w * A
    expr2 <- paste("cbind(", paste(rep("B", cola), collapse = ","), 
                              ")", sep = "")
    B <- eval(parse(text = expr2))
  }
  else {
    if (ncol(A) != ncol(B)) 
      stop("Dimensions of the matrices do not match.")
    if (missing(w)) 
      w <- rep(1, ncol(A))
    if (ncol(A) != length(w)) 
      stop("Length of the weight does not match with the dimension of the matrices.")
    rowa <- nrow(A)
    rowb <- nrow(B)
    rowab <- rowa * rowb
    A <- matrix(rep(A, each = rowb), rowab, )
    A <- A * matrix(w, nrow(A), ncol(A), byrow = TRUE)
    expr <- paste("rbind(", paste(rep("B", rowa), collapse = ","), 
                             ")", sep = "")
    B <- eval(parse(text = expr))
  }
  return(A * B)
}
