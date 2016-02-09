pspline.setting <- function(x,knots=select_knots(x,35),p=3,m=2,periodicity=FALSE,weight=NULL){
  
# x: the marginal data points
# knots: the list of interior knots or the numbers of interior knots
# p: degrees for B-splines, with defaults values 3
# m: orders of difference penalty, with default values 2
#require(splines)
#require(Matrix)

### design matrix 
K = length(knots)-2*p-1
B = spline.des(knots=knots, x=x, ord = p+1,outer.ok = TRUE)$design
if(periodicity){
  Bint = B[,-c(1:p,K+1:p)]
  Bleft = B[,1:p]
  Bright = B[,K+1:p]
  B = cbind(Bint,Bleft+Bright)
}


difference.penalty <-function(m,p,K,periodicity=FALSE){
  
  # parameter  m: difference order
  # parameter  p: degree of B-splines
  # parameter  K: number of interior knots
  c = rep(0,m+1)
  
  for(i in 0:m)
    c[i+1] = (-1)^(i+1)*factorial(m)/(factorial(i)*factorial(m-i))
  
  if(!periodicity){
  
    M = matrix(0,nrow=K+p-m,ncol=K+p)
    for(i in 1:(K+p-m)) M[i,i:(i+m)] = c
  }
  if(periodicity){
    
    M = matrix(0,nrow=K,ncol=K)
    for(i in 1:(K-m)) M[i,i:(i+m)] = c
    for(i in (K-m+1):K) M[i,c(i:K,1:(m-K+i))] = c
  }
  
  return(M)
}


P = difference.penalty(m,p,K,periodicity)
P1 = Matrix(P)
P2 = Matrix(t(P))
P = P2%*%P1

MM <- function(A,s,option=1){
  if(option==2)
    return(A*(s%*%t(rep(1,dim(A)[2]))))
  if(option==1)
    return(A*(rep(1,dim(A)[1])%*%t(s)))
}

if(is.null(weight)) weight <- rep(1,length(x))


B1 = Matrix(MM(t(B),weight))
B = Matrix(B) 
Sig = B1%*%B
eSig = eigen(Sig)
V = eSig$vectors
E = eSig$values
if(min(E)<=0.0000001) {#cat("Warning! t(B)%*%B is singular!\n");
               #cat("A small identity matrix is added!\n");
               E <- E + 0.000001;
               
}
Sigi_sqrt = MM(V,1/sqrt(E))%*%t(V)

#Sigi = V%*%diag(1/E)%*%t(V)
tUPU = Sigi_sqrt%*%(P%*%Sigi_sqrt)
Esig = eigen(tUPU,symmetric=TRUE)
U = Esig$vectors
s = Esig$values
if(!periodicity) s[(K+p-m+1):(K+p)]=0
if(periodicity) s[K] = 0
A = B%*%(Sigi_sqrt%*%U)

List = list(
        "A" = A,
        "B" = B,
        "s" = s,
        "Sigi.sqrt"=Sigi_sqrt,
        "U" = U,
        "P" = P)

return(List)
}
