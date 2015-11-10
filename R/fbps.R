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
fbps <- function(data,subj=NULL,covariates = NULL, knots=35, knots.option="equally-spaced",
                 periodicity = c(FALSE,FALSE), p=3,m=2,lambda=NULL,
                 selection = "GCV",
                 search.grid = T, search.length = 100, method="L-BFGS-B",
                 lower= -20, upper=20, control=NULL){
    
    # return a smoothed matrix using fbps
    
    # data: a matrix 
    # covariates: the list of data points for each dimension
    # knots: to specify either the number/numbers of  knots  or the vector/vectors of knots for each dimension; defaults to 35
    # p: the degrees of B-splines
    # m: the order of difference penalty
    # lambda: the user-selected smoothing parameters  
    # lscv: for leave-one-subject-out cross validation, the columns are subjects
    # method: see "optim"
    # lower, upper, control: see "optim"
    
    #require(splines)
    #require(Matrix)
    #source("pspline.setting.R")
    
    ## data dimension
    data_dim = dim(data)
    n1 = data_dim[1]
    n2 = data_dim[2]
    
    ## subject ID
    if(is.null(subj)) subj = 1:n2
    subj_unique = unique(subj)
    I = length(subj_unique)
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
      xknots = select_knots(x,knots=K1,option=knots.option)
      K2 = rep(knots,2)[2]
      zknots = select_knots(z,knots=K2,option=knots.option)
    }
    
    if(is.list(knots)){
      
      xknots = knots[[1]]
      K1 = length(xknots)-1 
      knots_left <- 2*xknots[1]-xknots[p1:1+1]
      knots_right <- 2*xknots[K1] - xknots[K1-(1:p1)]
      xknots <- c(knots_left,xknots,knots_right)
      
      zknots= knots[[2]]
      K2 = length(zknots)-1
      knots_left <- 2*zknots[1]- zknots[p2:1+1]
      knots_right <- 2*zknots[K2] - zknots[K2-(1:p2)]
      zknots <- c(knots_left,zknots,knots_right)
    }
  #######################################################################################
    Y = data 
    
    ###################   precalculation for fbps smoothing  ##########################################66
    
    List = pspline.setting(x,xknots,p1,m1,periodicity[1])
    A1 = List$A
    B1 = List$B
    Bt1 = Matrix(t(as.matrix(B1)))
    s1 = List$s
    Sigi1_sqrt = List$Sigi.sqrt
    U1 = List$U
    A01 = Sigi1_sqrt%*%U1
    c1 = length(s1)
    
    List = pspline.setting(z,zknots,p2,m2,periodicity[2])
    A2 = List$A
    B2 = List$B
    Bt2 = Matrix(t(as.matrix(B2)))
    s2 = List$s
    Sigi2_sqrt = List$Sigi.sqrt
    U2 = List$U
    A02 = Sigi2_sqrt%*%U2
    c2 = length(s2)
    #################select optimal penalty ####################################
    
    tr <-function(A){ return(sum(diag(A)))} ## the trace of a square matrix
    
    Ytilde = Bt1%*%(Y%*%B2)
    Ytilde = t(A01)%*%Ytilde%*%A02
    Y_sum = sum(Y^2)
    ytilde = as.vector(Ytilde)
    if(selection=="iGCV"){
      
      KH = function(A,B){
        C = matrix(0,dim(A)[1],dim(A)[2]*dim(B)[2])
        for(i in 1:dim(A)[1])
          C[i,] = kronecker(A[i,],B[i,])
        return(C)
      }
      G = rep(0,I)
      Ybar = matrix(0,c1,I)
      C = matrix(0,c2,I)
      
      Y2 = Bt1%*%Y
      Y2 = matrix(t(A01)%*%Y2,c1,n2)
      for(i in 1:I){
        sel = (1:n2)[subj==subj_unique[i]]
        len = length(sel)
        G[i] = len
        Ybar[,i] = as.vector(matrix(Y2[,sel],ncol=len)%*%rep(1,len))
        C[,i] = as.vector(t(matrix(A2[sel,],nrow=len))%*%rep(1,len))
      }
      g1 = diag(Ybar%*%diag(G)%*%t(Ybar))
      g2 = diag(Ybar%*%t(Ybar))
      g3 = ytilde*as.vector(Ybar%*%diag(G)%*%t(C))
      g4 = ytilde*as.vector(Ybar%*%t(C))
      #g5 = as.vector(C%*%diag(G)%*%t(C))
      #g6 = as.vector(C%*%t(C))
      #g7 = as.vector(KH(Ytilde,Ytilde))
      g5 = diag(C%*%diag(G)%*%t(C))
      g6 = diag(C%*%t(C))
      #cat("Processing completed\n") 
    }
    
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
      
      sigma2 = 1/(1+lambda2*s2)
      sigma1 = 1/(1+lambda1*s1)
      sigma2_sum = sum(sigma2)
      sigma1_sum = sum(sigma1)
      sigma = kronecker(sigma2,sigma1)
      sigma.2 = kronecker(sqrt(sigma2),sqrt(sigma1))
      
      if(selection=="iGCV"){
                            d = 1/(1-(1+lambda1*s1)/sigma2_sum*n2)
                            gcv = sum((ytilde*sigma)^2) - 2*sum((ytilde*sigma.2)^2)
                            gcv = gcv + sum(d^2*g1) - 2*sum(d*g2)
                            gcv = gcv - 2*sum(g3*kronecker(sigma2,sigma1*d^2))
                            gcv = gcv + 4*sum(g4*kronecker(sigma2,sigma1*d))
                            #sigma22 = kronecker(sigma2,sigma2)
                            #gcv = gcv + sum(g7*kronecker(sigma22*g5,sigma1^2*d^2))
                            #gcv = gcv - 2*sum(g7*kronecker(sigma22*g6,sigma1^2*d))
                            gcv = gcv + sum(ytilde^2*kronecker(sigma2^2*g5,sigma1^2*d^2))
                            gcv = gcv - 2*sum(ytilde^2*kronecker(sigma2^2*g6,sigma1^2*d))
                            
      }
      
      if(selection=="GCV") {
                           
                            gcv = sum((ytilde*sigma)^2) - 2*sum((ytilde*sigma.2)^2)
                            gcv = Y_sum + gcv
                            trace = sigma2_sum*sigma1_sum
                            gcv = gcv/(1-trace/(n1*n2))^2
      }
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
      
      sigma2 = 1/(1+lambda2*s2)
      sigma1 = 1/(1+lambda1*s1)
      sigma2_sum = sum(sigma2)
      sigma1_sum = sum(sigma1)
      sigma = kronecker(sigma2,sigma1)
      sigma.2 = kronecker(sqrt(sigma2),sqrt(sigma1))
      
      Theta = A01%*%diag(sigma1)%*%Ytilde
      Theta = as.matrix(Theta%*%diag(sigma2)%*%t(A02))
      Yhat = as.matrix(as.matrix(B1%*%Theta)%*%Bt2)
      #result=list(lambda=c(lambda1,lambda2),Yhat=Yhat,Theta=Theta)
      result=list(lambda=c(lambda1,lambda2),Yhat=Yhat,Theta=Theta,
                  setting = list(x = list(knots = xknots, p = p1, m = m1), 
                                 z = list(knots = zknots, p = p2, m = m2)))
      class(result) ="fbps"
      return(result)
    }
    
    if(is.null(lambda)){
      
      if(search.grid ==T){
        lower2 <- lower1 <- lower[1]
        upper2 <- upper1 <- upper[1]
        search.length2 <- search.length1 <- search.length[1]
        if(length(lower)==2) lower2 <- lower[2]
        if(length(upper)==2) upper2 <- upper[2]
        if(length(search.length)==2) search.length2 <- search.length[2]
        
        Lambda1 = seq(lower1,upper1,length = search.length1)
        Lambda2 = seq(lower2,upper2,length = search.length2)
        
        lambda.length1 = length(Lambda1)
        lambda.length2 = length(Lambda2)
        
        GCV = matrix(0,lambda.length1,lambda.length2)
        for(j in 1:lambda.length1)
          for(k in 1:lambda.length2){
            GCV[j,k] = fbps_gcv(c(Lambda1[j],Lambda2[k]))
          }
        location = which(GCV==min(GCV))[1]
        j0 = location%%lambda.length1
        if(j0==0) j0 = lambda.length1
        k0 = (location-j0)/lambda.length1+1
        lambda = exp(c(Lambda1[j0],Lambda2[k0]))
      } ## end of search.grid
      
      if(search.grid == F){
        fit = optim(0,fbps_gcv,method=method,control=control,
                    lower=lower[1],upper=upper[1])
        
        fit = optim(c(fit$par,fit$par),fbps_gcv,method=method,control=control,
                    lower=lower[1],upper=upper[1])
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
