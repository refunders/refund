##' Two-step function-on-scalar regression
##'
##' This function performs linear regression with functional responses and
##' scalar predictors by (1) fitting a separate linear model at each point
##' along the function, and then (2) smoothing the resulting coefficients to
##' obtain coefficient functions.
##'
##' Unlike \code{\link{fosr}} and \code{\link{pffr}}, which obtain smooth
##' coefficient functions by minimizing a penalized criterion, this function
##' introduces smoothing only as a second step. The idea was proposed by Fan
##' and Zhang (2000), who employed local polynomials rather than roughness
##' penalization for the smoothing step.
##'
##' @param Y the functional responses, given as an \eqn{n\times d} matrix.
##' @param X \eqn{n\times p} model matrix, whose columns represent scalar
##' predictors. Should ordinarily include a column of 1s.
##' @param argvals the \eqn{d} argument values at which the functional
##' responses are evaluated, and at which the coefficient functions will be
##' evaluated.
##' @param nbasis number of basis functions used to represent the coefficient
##' functions.
##' @param norder norder of the spline basis, when \code{basistype="bspline"}
##' (the default, 4, gives cubic splines).
##' @param pen.order order of derivative penalty.
##' @param basistype type of basis used. The basis is created by an appropriate
##' constructor function from the \pkg{fda} package; see
##' \code{\link[fda]{basisfd}}. Only \code{"bspline"} and \code{"fourier"} are
##' supported.
##' @return An object of class \code{fosr}, which is a list with the following
##' elements: \item{fd}{object of class \code{"\link{fd}"} representing the
##' estimated coefficient functions. Its main components are a basis and a
##' matrix of coefficients with respect to that basis. }
##' \item{raw.coef}{\eqn{d\times p} matrix of coefficient estimates from
##' regressing on \code{X} separately at each point along the function. }
##' \item{raw.se}{\eqn{d\times p} matrix of standard errors of the raw
##' coefficient estimates. } \item{yhat}{\eqn{n\times d} matrix of fitted
##' values. } \item{est.func}{\eqn{d\times p} matrix of coefficient function
##' estimates, obtained by smoothing the columns of \code{raw.coef}. }
##' \item{se.func}{\eqn{d\times p} matrix of coefficient function standard
##' errors. } \item{argvals}{points at which the coefficient functions are
##' evaluated. } \item{lambda}{smoothing parameters (chosen by REML) used to
##' smooth the \eqn{p} coefficient functions with respect to the supplied
##' basis. }
##' @author Philip Reiss \email{phil.reiss@@nyumc.org} and Lan Huo
##' @seealso \code{\link{fosr}}, \code{\link{pffr}}
##' @references Fan, J., and Zhang, J.-T. (2000). Two-step estimation of
##' functional linear models with applications to longitudinal data.
##' \emph{Journal of the Royal Statistical Society, Series B}, 62(2), 303--322.
##' @examples
##'
##' require(fda)
##'
##' # Effect of latitude on daily mean temperatures
##' tempmat = t(CanadianWeather$dailyAv[,,1])
##' latmat = cbind(1, scale(CanadianWeather$coord[ , 1], TRUE, FALSE))  # centred!
##' fzmod <- fosr2s(tempmat, latmat, argvals=day.5, basistype="fourier", nbasis=25)
##'
##' par(mfrow=1:2)
##' ylabs = c("Intercept", "Latitude effect")
##' for (k in 1:2) {
##' 	with(fzmod,matplot(day.5, cbind(raw.coef[,k],raw.coef[,k]-2*raw.se[,k],
##' 	     raw.coef[,k]+2*raw.se[,k],est.func[,k],est.func[,k]-2*se.func[,k],
##' 	     est.func[,k]+2*se.func[,k]), type=c("p","l","l","l","l","l"),pch=16,
##' 	     lty=c(1,2,2,1,2,2),col=c(1,1,1,2,2,2), cex=.5,axes=FALSE,xlab="",ylab=ylabs[k]))
##'     axesIntervals()
##'     box()
##'     if (k==1) legend("topleft", legend=c("Raw","Smoothed"), col=1:2, lty=2)
##' }
##'
##' @export
##' @importFrom fda create.bspline.basis create.fourier.basis eval.basis getbasispenalty
##' @importFrom mgcv gam
fosr2s <- function (Y, X, argvals = seq(0,1,,ncol(Y)), nbasis = 15, norder = 4,
                       pen.order=norder-2, basistype="bspline")
{
    # Stage 1: raw estimates
    n = dim(X)[1]
    p = dim(X)[2]
    XtX.inv = solve(crossprod(X))
    raw.coef = t(XtX.inv %*% crossprod(X, Y))
    resmat = Y - X %*% t(raw.coef)
    covmat = cov(resmat)  # TODO: add more sophisticated covariance methods?
    sigma2 = apply(resmat, 2, crossprod) / (n-p)
    raw.se = t(sqrt(diag(XtX.inv) %o% sigma2))

    # Stage 2: smooth raw estimates to get coefficient functions
    if (basistype=="bspline") bss = create.bspline.basis(range(argvals), nbasis = nbasis, norder = norder)
    else if (basistype=="fourier") bss = create.fourier.basis(range(argvals), nbasis = nbasis)

    Bmat = eval.basis(argvals, bss)
    P = getbasispenalty(bss, pen.order)
    Bt.Sig.B <- crossprod(Bmat, covmat %*% Bmat)

    lambda=c()
    coefmat <- matrix(NA, nbasis, p)
    est.func <- se.func <- matrix(NA,length(argvals),p)
    for (j in 1:p) {
    	swmod <- gam(raw.coef[ , j] ~ Bmat-1, paraPen=list(Bmat=list(P)), method="REML")
    	lambda[j] <- swmod$sp
    	coefmat[ , j] <- swmod$coef
    	est.func[ , j] <- fitted(swmod)
    	m1 <- Bmat %*% solve(crossprod(Bmat) + swmod$sp * P)
    	se.func[ , j] <- sqrt(XtX.inv[j,j] * rowSums(m1 * (m1 %*% Bt.Sig.B)))
    }

    list(fd=fd(coef=coefmat, basisobj=bss), raw.coef=raw.coef, raw.se=raw.se,
         yhat=tcrossprod(X, est.func),
         est.func=est.func, se.func=se.func, argvals=argvals, lambda=lambda)
}
