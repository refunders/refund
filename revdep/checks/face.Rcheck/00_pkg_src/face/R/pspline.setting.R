pspline.setting <- function(x,knots=select.knots(x,35),
                            p=3,m=2,weight=NULL,type="full",
                            knots.option="equally-spaced"){
  
# x: the marginal data points
# knots: the list of interior knots or the numbers of interior knots
# p: degrees for B-splines, with defaults values 3
# m: orders of difference penalty, with default values 2
# knots.option: type of knots placement, with default values "equally-spaced"

#require(splines)
#require(Matrix)

### design matrix 
K = length(knots)-2*p-1
B = spline.des(knots=knots, x=x, ord = p+1,outer.ok = TRUE,sparse=TRUE)$design

bs = "ps"

if(knots.option == "quantile"){
  bs = "bs"
}

s.object = s(x=x, bs=bs, k=K+p,m=c(p-1,2), sp=NULL)
object  = smooth.construct(s.object,data = data.frame(x=x),knots=list(x=knots))
P =  object$S[[1]]
if(knots.option == "quantile") P = P / max(abs(P))*10 # rescaling
  
if(is.null(weight)) weight <- rep(1,length(x))

if(type=="full"){

Sig = crossprod(matrix.multiply(B,weight,option=2),B)
eSig = eigen(Sig)
V = eSig$vectors
E = eSig$values
if(min(E)<=0.0000001) {#cat("Warning! t(B)%*%B is singular!\n");
               #cat("A small identity matrix is added!\n");
               E <- E + 0.000001;
               
}
Sigi_sqrt = matrix.multiply(V,1/sqrt(E))%*%t(V)

tUPU = Sigi_sqrt%*%(P%*%Sigi_sqrt)
Esig = eigen(tUPU,symmetric=TRUE)
U = Esig$vectors
s = Esig$values
s[(K+p-m+1):(K+p)]=0
A = B%*%(Sigi_sqrt%*%U)
}

if(type=="simple"){
  
  A = NULL
  s = NULL
  Sigi_sqrt = NULL
  U = NULL
}
List = list(
        "A" = A,
        "B" = B,
        "s" = s,
        "Sigi.sqrt" = Sigi_sqrt,
        "U" = U,
        "P" = P)

return(List)
}
