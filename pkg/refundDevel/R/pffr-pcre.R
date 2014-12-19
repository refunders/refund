#' mgcv-style constructor for PC-basis functional random effects
#' 
#' Sets up design matrix for functional random effects based on the PC loadings
#' of the covariance operator of the random effect process. 
#' See \code{\link[mgcv]{smooth.construct.re.smooth.spec}} for more details on \code{mgcv}-style smoother specification
#' and \code{\link{pcre}} for the corresponding \code{pffr()}-formula wrapper. 
#' 
#' @param object a smooth specification object, see \code{\link[mgcv]{smooth.construct}}
#' @param data  see \code{\link[mgcv]{smooth.construct}}
#' @param knots see \code{\link[mgcv]{smooth.construct}}
#' @method smooth.construct pcre.smooth.spec
#' @return An object of class \code{"random.effect"}. See \code{\link[mgcv]{smooth.construct}} 
#'  for the elements that this object will contain. 
#' @author Fabian Scheipl;  adapted from 're' constructor by S.N. Wood.
smooth.construct.pcre.smooth.spec <- function(object, data, knots){
    if (!is.null(object$id)) 
        stop("random effects don't work with ids.")
    form <- as.formula(paste("~", paste(object$term[1], ":", 
                                        paste("(",paste(object$term[-1],
                                                        collapse="+"), 
                                              ")")), 
                    "-1"))
    # structure of X: [eigenf.1 for all ids | eigenf.2 for all ids | ...]
    object$X <- model.matrix(form, data)#Matrix:::sparse.
    object$bs.dim <- ncol(object$X)
    object$S <- list(diag(object$bs.dim))
    object$rank <- object$bs.dim
    object$null.space.dim <- 0
    object$C <- #diag(length(object$term)-1) %x% 
        #t(rep(1, nlevels(data[[object$term[1]]]))) 
        matrix(0, 0, ncol(object$X))
    object$form <- form
    object$side.constrain <- FALSE
    object$plot.me <- TRUE
    object$te.ok <- 2
    class(object) <- c("pcre.random.effect", "random.effect")
    object
}


#' pffr-constructor for functional principal component-based functional random intercepts.
#' 
#' @section Details: Fits functional random intercepts \eqn{B_i(t)} for a grouping variable \code{id} 
#' using as a basis the functions \eqn{\phi_m(t)} in \code{efunctions} with variances \eqn{\lambda_m} in \code{evalues}:
#' \eqn{B_i(t) \approx \sum_m^M \phi_m(t)\delta_{im}} with  
#' independent \eqn{\delta_{im} \sim N(0, \sigma^2\lambda_m)}, where \eqn{\sigma^2}
#' is (usually) estimated and controls the overall contribution of the \eqn{B_i(t)} while the relative importance
#' of the \eqn{M} basisfunctions is controlled by the supplied variances \code{lambda_m}.
#' Can be used to model smooth residuals if \code{id} is simply an index of observations.
#' 
#' \code{efunctions} and \code{evalues} are typically eigenfunctions and eigenvalues of an estimated 
#' covariance operator for the functional process to be modeled, i.e., they are
#' a functional principal components basis. 
#'  
#' @param id grouping variable a factor
#' @param efunctions matrix of eigenfunction evaluations on gridpoints \code{yind} (<length of \code{yind}> x <no. of used eigenfunctions>)
#' @param evalues eigenvalues associated with \code{efunctions}
#' @param yind vector of gridpoints on which \code{efunctions} are evaluated.
#' @param ... not used
#' @return a list used internally for constructing an appropriate call to \code{mgcv::gam}
#' @author Fabian Scheipl 
#' @examples \dontrun{
#' residualfunction <- function(t){
#' #generate quintic polynomial error functions 
#'     drop(poly(t, 5)%*%rnorm(5, sd=sqrt(2:6)))
#' }
#' # generate data Y(t) = mu(t) + E(t) + white noise
#' set.seed(1122)
#' n <- 50
#' T <- 30
#' t <- seq(0,1, l=T)
#' # E(t): smooth residual functions
#' E <- t(replicate(n, residualfunction(t)))
#' int <- matrix(scale(3*dnorm(t, m=.5, sd=.5) - dbeta(t, 5, 2)), byrow=T, n, T)
#' Y <- int + E + matrix(.2*rnorm(n*T), n, T)
#' data <- data.frame(Y=I(Y))
#' # fit model under independence assumption:
#' summary(m0 <- pffr(Y ~ 1, yind=t, data=data))
#' # get first 5 eigenfunctions of residual covariance
#' # (i.e. first 5 functional PCs of empirical residual process)
#' Ehat <- resid(m0)
#' fpcE <- fpca.sc(Ehat, npc=5)
#' efunctions <- fpcE$efunctions
#' evalues <- fpcE$evalues
#' data$id <- factor(1:nrow(data))
#' # refit model with fpc-based residuals
#' m1 <- pffr(Y ~ 1 + pcre(id=id, efunctions=efunctions, evalues=evalues, yind=t), yind=t, data=data)
#' t1 <- predict(m1, type="terms")
#' summary(m1)
#' #compare squared errors
#' mean((int-fitted(m0))^2)
#' mean((int-t1[[1]])^2)
#' mean((E-t1[[2]])^2)
#' # compare fitted & true smooth residuals and fitted intercept functions:
#' layout(t(matrix(1:4,2,2)))
#' matplot(t(E), lty=1, type="l", ylim=range(E, t1[[2]]))
#' matplot(t(t1[[2]]), lty=1, type="l", ylim=range(E, t1[[2]]))
#' plot(m1, select=1, main="m1", ylim=range(Y))
#' lines(t, int[1,], col=rgb(1,0,0,.5))
#' plot(m0, select=1, main="m0", ylim=range(Y))
#' lines(t, int[1,], col=rgb(1,0,0,.5))
#' }
pcre <- function(id, 
        efunctions,
        evalues,
        yind,
        ...
){
    # check args
    stopifnot(is.factor(id), nrow(efunctions)==length(yind), 
              ncol(efunctions)==length(evalues), all(evalues>0))
    
    phiname <- deparse(substitute(efunctions))
    idname <- paste(deparse(substitute(id)),".vec",sep="")
    
    #scale eigenfunctions by their eigenvalues:
    efunctions <- t(t(efunctions)*sqrt(evalues))
    
    #assign unique names based on the given args
    colnames(efunctions) <- paste(phiname,".PC", 1:ncol(efunctions), sep="")
    
    call <- as.call(c(as.symbol("s"),
                    as.symbol(substitute(idname)),
                    sapply(colnames(efunctions), function(x) as.symbol(x)),
                    bs=c("pcre")))
    
    return(list(efunctions=efunctions, yind=yind, idname=idname, 
                id=id, call=call, ...))
}#end pcre()
