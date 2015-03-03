##' Sandwich smoother for matrix data
##'
##' A fast bivariate \emph{P}-spline method for smoothing matrix data.
##'
##' The smoothing parameter can be user-specified; otherwise, the function uses
##' grid searching method or \command{optim} for selecting the smoothing
##' parameter.
##'
##' @param data n1 by n2 data matrix without missing data
##' @param covariates list of two vectors of covariates of lengthes n1 and n2;
##' if NULL, then generates equidistant covariates
##' @param knots list of two vectors of knots or number of equidistant knots
##' for all dimensions; defaults to 35
##' @param p degrees of B-splines; defaults to 3
##' @param m order of differencing penalty; defaults to 2
##' @param lambda user-specified smoothing parameters; defaults to NULL
##' @param alpha the coefficient of the penalty for degrees of freedom in the
##' GCV criterion; defaults to 1
##' @param search.grid logical; defaults to TRUE, if FALSE, uses
##' \command{optim}
##' @param search.length number of equidistant (log scale) smoothing parameter;
##' defaults to 100
##' @param method see \command{optim}; defaults to \command{L-BFGS-B}
##' @param lower,upper bounds for log smoothing parameter, passed to
##' \command{optim}; defaults are -20 and 20.
##' @param control see \command{optim}
##' @return A list with components \item{lambda}{vector of length 2 of selected
##' smoothing parameters} \item{Yhat}{fitted data} \item{trace}{trace of the
##' overall smoothing matrix} \item{gcv}{value of generalized cross validation}
##' \item{Theta}{matrix of estimated coefficients}
##' @author Luo Xiao \email{lxiao@@jhsph.edu}
##' @export
##' @importFrom Matrix kronecker as.matrix
##' @references Xiao, L., Li, Y., and Ruppert, D. (2013). Fast bivariate
##' \emph{P}-splines: the sandwich smoother. \emph{Journal of the Royal
##' Statistical Society: Series B}, 75(3), 577--599.
##' @examples
##'
##' ##########################
##' #### True function   #####
##' ##########################
##' n1 <- 60
##' n2 <- 80
##' x <- (1:n1)/n1-1/2/n1
##' z <- (1:n2)/n2-1/2/n2
##' MY <- array(0,c(length(x),length(z)))
##'
##' sigx <- .3
##' sigz <- .4
##' for(i in 1:length(x))
##' for(j in 1:length(z))
##' {
##' #MY[i,j] <- .75/(pi*sigx*sigz) *exp(-(x[i]-.2)^2/sigx^2-(z[j]-.3)^2/sigz^2)
##' #MY[i,j] <- MY[i,j] + .45/(pi*sigx*sigz) *exp(-(x[i]-.7)^2/sigx^2-(z[j]-.8)^2/sigz^2)
##' MY[i,j] = sin(2*pi*(x[i]-.5)^3)*cos(4*pi*z[j])
##' }
##' ##########################
##' #### Observed data   #####
##' ##########################
##' sigma <- 1
##' Y <- MY + sigma*rnorm(n1*n2,0,1)
##' ##########################
##' ####   Estimation    #####
##' ##########################
##'
##' est <- fbps(Y,list(x=x,z=z))
##' mse <- mean((est$Yhat-MY)^2)
##' cat("mse of fbps is",mse,"\n")
##' cat("The smoothing parameters are:",est$lambda,"\n")
##' ########################################################################
##' ########## Compare the estimated surface with the true surface #########
##' ########################################################################
##'
##' par(mfrow=c(1,2))
##' persp(x,z,MY,zlab="f(x,z)",zlim=c(-1,2.5), phi=30,theta=45,expand=0.8,r=4,
##'       col="blue",main="True surface")
##' persp(x,z,est$Yhat,zlab="f(x,z)",zlim=c(-1,2.5),phi=30,theta=45,
##'       expand=0.8,r=4,col="red",main="Estimated surface")
##' @importFrom stats optim
fbps <-
function(data, covariates = NULL, knots=35, p=3,m=2,lambda=NULL,alpha=1,
         search.grid = T, search.length = 100, method="L-BFGS-B",
         lower= -20, upper=20, control=NULL){

# return a smoothed matrix using fbps

# data: a matrix
# covariates: the list of data points for each dimension
# knots: to specify either the number/numbers of  knots  or the vector/vectors of knots for each dimension; defaults to 35
# p: the degrees of B-splines
# m: the order of difference penalty
# lambda: the user-selected smoothing parameters
# method: see "optim"
# lower, upper, control: see "optim"

#library(fBasics)
#source("pspline.setting.R")

## data dimension
data_dim = dim(data)
n1 = data_dim[1]
n2 = data_dim[2]

## covariates for the two axis
if(!is.list(covariates)) {

  x=(1:n1)/n1-1/2/n1; ## if NULL, assume equally distributed
  z = (1:n2)/n2-1/2/n2
  }
if(is.list(covariates)){

  x = covariates[[1]]
  z = covariates[[2]]
 }

## B-spline basis setting
p1 = rep(p,2)[1]
p2 = rep(p,2)[2]
m1 = rep(m,2)[1]
m2 = rep(m,2)[2]

## knots
if(!is.list(knots)){

  K1 = rep(knots,2)[1]
  xknots=seq(-p1,K1+p1,length=K1+1+2*p1)/K1
  xknots = xknots*(max(x)-min(x)) + min(x)
  K2 = rep(knots,2)[2]
  zknots=seq(-p2,K2+p2,length=K2+1+2*p2)/K2
  zknots = xknots*(max(z)-min(z)) + min(z)

}

if(is.list(knots)){

  xknots = knots[[1]]
  K1 = length(xknots)-2*p1-1
  zknots= knots[[2]]
  K2 = length(zknots)-2*p2-1
}
#######################################################################################
Y = data

###################   precalculation for fbps smoothing  ##########################################66

List = pspline.setting(x,xknots,p1,m1)
A1 = List$A
s1 = List$s
Sigi1_sqrt = List$Sigi.sqrt
U1 = List$U

List = pspline.setting(z,zknots,p2,m2)
A2 = List$A
s2 = List$s
Sigi2_sqrt = List$Sigi.sqrt
U2 = List$U

##B1 and B2 are the B-spline design matrices
B1 = spline.des(knots=xknots, x=x, ord = p1+1,outer.ok = TRUE,sparse=FALSE)$design
B2 = spline.des(knots=zknots, x=z, ord = p2+1,outer.ok = TRUE,sparse=FALSE)$design
#################select optimal penalty ####################################

tr <-function(A){ return(sum(diag(A)))} ## the trace of a square matrix

Ytilde = as.matrix(t(A1)%*%Y%*%A2)
Y_sum = sum(Y^2)
ytilde = as.vector(Ytilde)

fbps_gcv =function(x){

lambda=exp(x)
## two lambda's are the same
if(length(lambda)==1)
{
lambda1 = lambda
lambda2 = lambda
}
## two lambda's are different
if(length(lambda)==2){
 lambda1=lambda[1]
 lambda2=lambda[2]
}


sigma = kronecker(1/(1+lambda2*s2),1/(1+lambda1*s1))
sigma.2 = sqrt(sigma)

gcv = Y_sum + sum((ytilde*sigma)^2) - 2*sum((ytilde*sigma.2)^2)
trace = sum(1/(1+lambda1*s1))*sum(1/(1+lambda2*s2))
gcv = gcv/(1-alpha*trace/(n1*n2))^2
return(gcv)
}

fbps_est =function(x){

lambda=exp(x)
## two lambda's are the same
if(length(lambda)==1)
{
lambda1 = lambda
lambda2 = lambda
}
## two lambda's are different
if(length(lambda)==2){
 lambda1=lambda[1]
 lambda2=lambda[2]
}

sigma = kronecker(1/(1+lambda2*s2),1/(1+lambda1*s1))
sigma.2 = sqrt(sigma)

gcv = Y_sum + sum((ytilde*sigma)^2) - 2*sum((ytilde*sigma.2)^2)
trace = sum(1/(1+lambda1*s1))*sum(1/(1+lambda2*s2))
gcv = gcv/(1-alpha*trace/(n1*n2))^2

Theta = Sigi1_sqrt%*%U1%*%diag(1/(1+lambda1*s1))%*%Ytilde
Theta = as.matrix(Theta%*%diag(1/(1+lambda2*s2))%*%t(U2)%*%Sigi2_sqrt)
Yhat = as.matrix(B1%*%Theta%*%t(B2))
result=list(lambda=c(lambda1,lambda2),Yhat=Yhat,trace=trace,gcv=gcv,Theta=Theta)
return(result)
}

if(is.null(lambda)){

if(search.grid ==T){

  Lambda = seq(lower,upper,length = search.length)
  lambda.length = length(Lambda)
  GCV = matrix(0,lambda.length,lambda.length)
  for(j in 1:lambda.length)
    for(k in 1:lambda.length){
      GCV[j,k] = fbps_gcv(c(Lambda[j],Lambda[k]))
    }
  location = which(GCV==min(GCV))
  j0 = location%%lambda.length
  if(j0==0) j0 = lambda.length
  k0 = (location-j0)/lambda.length+1
  lambda = exp(c(Lambda[j0],Lambda[k0]))
} ## end of search.grid

if(search.grid == F){
fit = optim(0,fbps_gcv,method=method,control=control,
  lower=rep(lower,2)[1],upper=rep(upper,2)[1])

fit = optim(c(fit$par,fit$par),fbps_gcv,method=method,control=control,
  lower=rep(lower,2)[1:2],upper=rep(upper,2)[1:2])
if(fit$convergence>0) {
  expression = paste("Smoothing failed! The code is:",fit$convergence)
  print(expression)
}
lambda = exp(fit$par)
} ## end of optim

} ## end of finding smoothing parameters
lambda = rep(lambda,2)[1:2]

return(fbps_est(log(lambda)))
}
