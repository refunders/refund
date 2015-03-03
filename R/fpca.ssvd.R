#'    Smoothed FPCA via iterative penalized rank one SVDs.
#'
#'    Implements the algorithm of Huang, Shen, Buja (2008) for finding smooth
#'    right singular vectors of a matrix \code{X} containing (contaminated) evaluations
#'    of functional random variables on a regular, equidistant grid. If the number of smooth SVs to
#'    extract is not specified, the function hazards a guess for the appropriate number based
#'    on the asymptotically optimal truncation threshold under the assumption of a low rank matrix
#'    contaminated with i.i.d. Gaussian noise with unknown variance derived in Donoho, Gavish (2013).
#'  Please note that Donoho, Gavish (2013) should be regarded as experimental for functional PCA, and will typically not work well if you have more observations than grid points.
#'
#'    @param Y data matrix (rows: observations; columns: grid of eval. points)
#'    @param argvals the argument values where functions are evaluated. It is implemented yet.
#'    @param npc how many smooth SVs to try to extract, if \code{NA} (the default) the hard thresholding
#'     rule of Donoho, Gavish (2013) is used (see Details, References).
#'    @param center center \code{Y} so that its column-means are 0? Defaults to \code{TRUE}
#'    @param maxiter how many iterations of the power algorithm to perform at most (defaults to 15)
#'    @param tol convergence tolerance for power algorithm (defaults to 1e-4)
#'    @param diffpen difference penalty order controlling the desired smoothness of the right singular vectors,
#'        defaults to 3 (i.e., deviations from local quadratic polynomials).
#'    @param gridsearch use \code{\link[stats]{optimize}} or a grid search to find GCV-optimal
#'        smoothing parameters? defaults to \code{TRUE}.
#'    @param alphagrid  grid of smoothing parameter values for grid search
#'    @param lower.alpha lower limit for for smoothing parameter if \code{!gridsearch}
#'    @param upper.alpha upper limit for smoothing parameter if \code{!gridsearch}
#'    @param verbose generate graphical summary of progress and diagnostic messages?
#'        defaults to \code{FALSE}
#'    @return a list like the returned object from \code{\link{fpca.sc}},  with entries
#'    \code{Yhat}, the smoothed trajectories, \code{scores}, the estimated FPC loadings, \code{mu},
#'    the column means of \code{Y} (or a vector of zeroes if \code{!center}),  \code{efunctions},
#'    the estimated smooth FPCs (note that these are orthonormal vectors, not evaluations of orthonormal functions...), \code{evalues}, their associated eigenvalues, and \code{npc}, the
#'     number of smooth components that were extracted.
#'    @seealso  \code{\link{fpca.sc}} and \code{\link{fpca.face}} for FPCA based on smoothing a covariance estimate;
#'      \code{\link{fpca2s}} for a faster SVD-based approach.
#'    @author Fabian Scheipl
#'    @export
#'    @importFrom stats optimize
#'    @references Huang, J. Z., Shen, H., and Buja, A. (2008).
#'     Functional principal components analysis via penalized rank one approximation.
#'     \emph{Electronic Journal of Statistics}, 2, 678-695
#'
#'     Donoho, D.L., and Gavish, M. (2013). The Optimal Hard Threshold for Singular Values is 4/sqrt(3).
#'     eprint arXiv:1305.5870. Available from \url{http://arxiv.org/abs/1305.5870}.
#'    @examples
#'  ## as in Sec. 6.2 of Huang, Shen, Buja (2008):
#'  set.seed(2678695)
#'  n <- 101
#'  m <- 101
#'  s1 <- 20
#'  s2 <- 10
#'  s <- 4
#'  t <- seq(-1, 1, l=m)
#'  v1 <- t + sin(pi*t)
#'  v2 <- cos(3*pi*t)
#'  V <- cbind(v1/sqrt(sum(v1^2)), v2/sqrt(sum(v2^2)))
#'  U <- matrix(rnorm(n*2), n, 2)
#'  D <- diag(c(s1^2, s2^2))
#'  eps <- matrix(rnorm(m*n, sd=s), n, m)
#'  Y <- U%*%D%*%t(V) + eps
#'
#'  smoothSV <- fpca.ssvd(Y, verbose=TRUE)
#'
#'  layout(t(matrix(1:4, nr=2)))
#'  clrs <- sapply(rainbow(n), function(c)
#'            do.call(rgb, as.list(c(col2rgb(c)/255, .1))))
#'  matplot(V, type="l", lty=1, col=1:2, xlab="",
#'          main="FPCs: true", bty="n")
#'  matplot(smoothSV$efunctions, type="l", lty=1, col=1:5, xlab="",
#'          main="FPCs: estimate", bty="n")
#'  matplot(1:m, t(U%*%D%*%t(V)), type="l", lty=1, col=clrs, xlab="", ylab="",
#'          main="true smooth Y", bty="n")
#'  matplot(1:m, t(smoothSV$Yhat), xlab="", ylab="",
#'          type="l", lty=1,col=clrs, main="estimated smooth Y", bty="n")
fpca.ssvd <- function(Y, argvals=NULL, npc=NA, center=TRUE,
        maxiter=15, tol=1e-4,
        diffpen=3,
        gridsearch=TRUE,
        alphagrid=1.5^(-20:40),
        lower.alpha=1e-5, upper.alpha=1e7,
        verbose=FALSE){

    if(!is.null(argvals)) stop("argvals is not supported.")

    #GCV criterion from eq. (10), App. C:
    gcv <- function(alpha, w, m, lambda){
        sqrt(sum( (w * (alpha*lambda)/(1 + alpha*lambda))^2 ))/
                (1- 1/m * sum( 1/(1 + alpha*lambda) ))
    }
    # Recursion for difference operator matrix
    makeDiffOp <- function(degree, dim){
        if(degree==0){
            return(diag(dim))
        } else {
            return(diff(makeDiffOp(degree-1, dim)))
        }
    }


    if(any(is.na(Y))) stop("No missing values in <Y> allowed.")
    m <- ncol(Y)
    n <- nrow(Y)

    if(is.na(npc)){
      npc <- getNPC.DonohoGavish(Y)
    }


    if(!is.numeric(npc)) stop("Invalid <npc>.")
    if(npc<1 | npc>min(m,n)) stop("Invalid <npc>.")
    if(verbose){
      cat("Using ", npc , "smooth components based on Donoho/Gavish (2013).\n")
    }

    if(!is.numeric(alphagrid)) stop("Invalid <alphagrid>.")
    if(any(is.na(alphagrid)) | any(alphagrid<.Machine$double.eps)) stop("Invalid <alphagrid>.")
    uhoh <- numeric(0)
    if(verbose & interactive()){
        par(ask=TRUE)
        on.exit(par(ask=FALSE))
    }
    if(verbose){
       cat("Singular values of smooth and non-smooth ('noise') parts:")
    }



    Omega  <- crossprod(makeDiffOp(degree=diffpen, dim=m))
    eOmega <- eigen(Omega, symmetric=TRUE)
    lambda <- eOmega$values
    Gamma  <- eOmega$vectors

    U <- matrix(NA, nrow=n, ncol=npc)
    V <- matrix(NA, nrow=m, ncol=npc)
    d <- rep(NA, npc)

    if(center){
        meanY <- predict(smooth.spline(x=1:m, y=colMeans(Y)), x=1:m)$y
        Y <- t(t(Y) - meanY)
    } else {
        meanY <- rep(0, ncol(Y))
    }

    Ynow <- Y

    for(k in 1:npc){
        ## 'Power algorithm', Section 3 / Appendix C :
        YnowGamma <- Ynow%*%Gamma
        vold <- svd(Ynow, nu=0, nv=1)$v[,1]
        u <- Ynow %*% vold
        iter <- 1; reldiff <- tol+1
        while(all(c(reldiff > tol, iter<=maxiter))){
            w <- t(YnowGamma)%*%u
            if(gridsearch){
                gridmin <- which.min(sapply(alphagrid, gcv, w=w, m=m, lambda=lambda))
                minalpha <- alphagrid[gridmin]
            } else {
                minalpha <- optimize(f=gcv, interval=c(lower.alpha, upper.alpha),
                        w=w, m=m, lambda=lambda)$minimum
            }
            # v = S(alpha)^-1 Y u = Gamma%*%(I + minalpha*Lambda)^-1 Gamma' Y u)
            vnew <- Gamma%*%(w/(1+minalpha*lambda))

            vnew <- vnew/sqrt(sum(vnew^2))
            reldiff <- sum((vold-vnew)^2)/sum(vold^2)
            iter <- iter +1
            vold <- vnew
            u <- Ynow %*% vold
        } # end while(reldiff)
        if(reldiff > tol){
            warning("Not converged for SV ", k,
                    "; relative difference was ", format(reldiff),".")
        }

        U[,k] <- u/sqrt(sum(u^2))
        V[,k] <- vnew
        d[k]  <- sqrt(sum(u^2))

        if(verbose){

            layout(t(matrix(1:6, ncol=2)))
            matlplot <- function(...) matplot(..., type="l", lty=1, col=1, lwd=.1)

            matlplot(t(Ynow), ylim=range(Ynow), main=bquote(Ynow[.(k)]), xlab="", ylab="", bty="n")
            matlplot(t(U[,k, drop=FALSE]%*%t(V[,k, drop=FALSE]*d[k])),
                    ylim=range(Ynow), main=bquote((UDV^T)[.(k)]), xlab="", ylab="", bty="n")
            matlplot(t(Ynow - U[,k, drop=FALSE]%*%t(V[,k, drop=FALSE]*d[k])),
                    ylim=range(Ynow), main=bquote(Ynow[.(k)] - (UDV^T)[.(k)]), xlab="", ylab="", bty="n")

            matlplot(t(Y), ylim=range(Y), main=bquote(Y), xlab="", ylab="", bty="n")
            matlplot(t(U[,1:k, drop=FALSE]%*%(t(V[,1:k, drop=FALSE])*d[1:k])),
                    ylim=range(Y), main=bquote(Y[.(k)]), xlab="", ylab="", bty="n")
            matlplot(t(Ynow - U[,k, drop=FALSE]%*%(t(V[,k, drop=FALSE])*d[k])),
                    ylim=range(Y), main=bquote(Y - Y[.(k)]), xlab="", ylab="", bty="n")
        }

        Ynow <- Ynow - U[, k, drop=FALSE] %*% (t(V[, k, drop=FALSE]) * d[k])
        noisesv <- svd(Ynow, nu=0, nv=0)$d[1]
        if(verbose){
            cat("k:",k, "-- smooth:", d[k], "-- 'noise':", noisesv, "-- alpha:", minalpha, "\n")
        }
        if(noisesv > 1.1 * d[k]){
            uhoh <- c(uhoh, k)
        }
    }# end for(k)
    if(length(uhoh)) warning("First SV for remaining un-smooth signal larger than ",
                "SV found for smooth signal for component(s) ", paste(uhoh, collapse=","))

#     return(list(smooth=list(d=d, u=U, v=V),
#                     noise=svd(Ynow, nu=min(dim(Y))-npc, nv=min(dim(Y))-npc),
#                     mean=meanY))
    scores <- U%*%(d*diag(length(d)))

    return(list(
      Yhat= t(meanY + t(scores%*%t(V))),
      scores=scores,
      mu=meanY,
      efunctions=V,
      evalues=d^2,
      npc=npc))

}


getNPC.DonohoGavish <- function(X){
  # use Donoho, Gavish (2013) for estimating suitable number of sv's to extract:

  n <- nrow(X)
  m <- ncol(X)

  beta <- n/m

  if(beta > 1 | beta < 1e-3){
    warning("Approximation for \\beta(\\omega) may be invalid.")
  }
  #            ## approx for omega.beta below eq. (25):
  #            betaplus  <- 1 + sqrt(beta)^2
  #            betaminus <- 1 - sqrt(beta)^2
  #            marcenkopastur <- function(x){
  #                abs(integrate(function(t){
  #                          sqrt((betaplus - t) * (t - betaminus))/(2*pi*t)
  #                        }, lower=betaminus, upper=x, subdivisions=1e5,
  #                        stop.on.error = FALSE)$value - 0.5)
  #            }
  #            mu.beta <- optimize(marcenkopastur,
  #                    interval=c(betaminus, betaplus))$minimum
  #            # eq. (10)
  #            lambda.beta <- sqrt(2*(beta+1) + (8*beta)/(beta + 1 + sqrt(beta^2 + 14*beta +1)))
  #            omega.beta <- lambda.beta/sqrt(mu.beta)

  omega.beta <- .56*beta^3 - 0.95*beta^2 + 1.82*beta + 1.43
  y <- svd(X, nu=0, nv=0)$d
  rankY <- min(which(cumsum(y[y>0])/sum(y[y>0]) > .995))
  y.med <- median(y)

  npc <- min(max(1, sum(y > omega.beta * y.med)),  rankY)
  return(npc)
}
