pspline.setting <-
function(x,knots=35,p=3,m=2){
  
# x: the marginal data points
# K: the list of knots or the numbers of knots
# p: degrees for B-splines, with defaults values 3
# m: orders of difference penalty, with default values 2
#require(Matrix)
if(length(knots)==1)
{
  K = knots
  knots=seq(-p,K+p,length=K+1+2*p)/K
  knots = knots*(max(x)-min(x)) + min(x)
  }

if(length(knots)>1) 
  { knots = knots
    K = length(knots)-2*p-1
}


difference.penalty <-function(m,K){
  
  # parameter  m: difference order
  # parameter  K: size  
  
  M = matrix(0,nrow=K-m,ncol=K)
  c = rep(0,m+1)
  
  for(i in 0:m)
    c[i+1] = (-1)^(i+1)*factorial(m)/(factorial(i)*factorial(m-i))
  
  for(i in 1:(K-m)) M[i,i:(i+m)] = c
  
  return(M)
}


P = difference.penalty(m,K+p)
P1 = Matrix(P)
P2 = Matrix(t(P))
P = P2%*%P1

MM <- function(A,s,option=1){
  if(option==2)
    return(A*(s%*%t(rep(1,dim(A)[2]))))
  if(option==1)
    return(A*(rep(1,dim(A)[1])%*%t(s)))
}

### design matrix and some pre-calculation 
### for the penalty without interaction
### The idea about pre-calculation, when lambda
## is changed, only eigenvalues change.

B = spline.des(knots=knots, x=x, ord = p+1,outer.ok = TRUE)$design
B1 = Matrix(t(B))
B = Matrix(B) 
Sig = B1%*%B
eSig = eigen(Sig)
V = eSig$vectors
E = eSig$values
Sigi_sqrt = MM(V,1/sqrt(E))%*%t(V)

#Sigi = V%*%diag(1/E)%*%t(V)
tUPU = Sigi_sqrt%*%(P%*%Sigi_sqrt)
Esig = eigen(tUPU)
U = Esig$vectors
s = Esig$values
s[(K+p-m+1):(K+p)]=0
A = B%*%(Sigi_sqrt%*%U)

List = list(
        "A" = A,
        "B" = B,
        "s" = s,
        "Sigi.sqrt"=Sigi_sqrt,
        "U" = U)

return(List)
}
