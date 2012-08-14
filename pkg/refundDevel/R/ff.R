#' Construct a function-on-function regression term
#' 
#' Defines a term \eqn{\int^{s_{hi, i}}_{s_{lo, i}} X_i(s)\beta(t,s)ds} 
#' for inclusion in an \code{mgcv::gam}-formula (or \code{bam} or \code{gamm} or \code{gamm4:::gamm}) as constructed
#' by \code{\link{pffr}}.
#' Defaults to a cubic tensor product B-spline with marginal first difference penalties for \eqn{\beta(t,s)} and
#' numerical integration over the entire range \eqn{[s_{lo, i}, s_{hi, i}] = [\min(s_i), \max(s_i)]} by using Simpson weights. 
#' Can't deal with any missing \eqn{X(s)}, unequal lengths of \eqn{X_i(s)} not (yet?) possible.
#' Unequal ranges for different \eqn{X_i(s)} should work. \eqn{X_i(s)} is assumed to be numeric.\cr
#' If \code{check.ident==TRUE} (the default), the routine tries to determine the effective rank of
#' the covariance operator of the \eqn{X_i(s)} and, if necessary, adjusts the number of marginal
#' basis functions for \code{xind} downwards from the default or supplied value in \code{splinepars} 
#' to ensure identifiability of \eqn{\beta(s,t)}. Specifically, the number of basis functions is limited
#' to be at most the number of eigenvalues accounting for at least .99 of the total variance in \eqn{X_i(s)}. 
#' 
#' @param X an n by \code{ncol(xind)} matrix of function evaluations \eqn{X_i(s_{i1}),\dots, X_i(s_{iS})}; \eqn{i=1,\dots,n}.
#' @param yind matrix (or vector) of indices of evaluations of \eqn{Y_i(t)}
#' @param xind matrix (or vector) of indices of evaluations of \eqn{X_i(s)}; i.e. matrix with rows \eqn{(s_{i1},\dots,s_{iS})}
#' @param basistype defaults to "\code{\link[mgcv]{te}}", i.e. a tensor product spline to represent \eqn{\beta(t,s)}. 
#' 		Alternatively, use \code{"s"} for bivariate basis functions (see \code{mgcv}'s \code{\link[mgcv]{s}}) 
#'      or \code{"t2"} for an alternative parameterization of tensor product splines (see \code{mgcv}'s \code{\link[mgcv]{t2}}).
#' @param integration method used for numerical integration. Defaults to \code{"simpson"}'s rule for calculating entries in \code{L}. 
#'  Alternatively and for non-equidistant grids, \code{"trapezoidal"} or \code{"riemann"}. \code{"riemann"} integration is always used 
#'  if \code{limits} is specified
#' @param L optional: an n by \code{ncol(xind)} matrix giving the weights for the numerical integration over \eqn{s}. 
#' @param limits : defaults to NULL for integration across the entire range of \eqn{X(s)}, otherwise 
#' specifies the integration limits \eqn{s_{hi, i}, s_{lo, i}}:  
#' either one of \code{"s<t"} or \code{"s<=t"} for \eqn{(s_{hi, i}, s_{lo, i}) = (0, t)} or 
#' a function that takes \code{s} as the first and \code{t} as the second argument and returns TRUE for combinations
#' of values \code{(s,t)} if \code{s} falls into the integration range for the given \code{t}. This is an experimental feature and not
#' well tested yet, use at your own risk.   
#' @param splinepars optional arguments supplied to the \code{basistype}-term. Defaults to a cubic tensor product 
#' 	B-spline with marginal first difference penalties, i.e. \code{list(bs="ps", m=c(2, 1))} See \code{\link[mgcv]{te}} or \code{\link[mgcv]{s}} in \code{mgcv} for details
#' @param check.ident check rank of \eqn{\operatorname{Cov}(X_i(s),X_i(s'))} and adjust number of basis functions 
#'  if necessary. See Details. Defaults to TRUE.
#' 
#' 
#' @seealso \code{mgcv}'s \code{\link[mgcv]{linear.functional.terms}}
#' @return a list containing \itemize{
#'  \item \code{call} a "call" to \code{\link[mgcv]{te}} (or \code{\link[mgcv]{s}}, \code{\link[mgcv]{t2}}) using the appropriately constructed covariate
#' 	 and weight matrices 
#'  \item \code{data} a list containing the necessary covariate and weight matrices
#' } 
#'  
#' @author Fabian Scheipl, Sonja Greven
# FIXME: weights for simpson's rule on non-equidistant grids 
# TODO: allow X to be of class fd (?)
# TODO: allow X to be a factor -- would result in one beta(s,t) surface for each level? (?)
# TODO: by variables
ff <- function(X,
		yind,
		xind=seq(0, 1, l=ncol(X)),
		basistype= c("te", "t2", "s"),
		integration=c("simpson", "trapezoidal", "riemann"), 
		L=NULL,
		limits=NULL, 
		splinepars=list(bs="ps", m=list(c(2, 1), c(2,1))),
        check.ident=TRUE
){
	n <- nrow(X)
	nxgrid <- ncol(X)
	stopifnot(all(!is.na(X)))
	
	
# check & format index for Y 
	if(!missing(yind))
	if(is.null(dim(yind))){
		yind <- t(t(yind))
	} 
	nygrid <- nrow(yind)
	
# check & format index for X
	if(is.null(dim(xind))){
		xind <- t(xind)
		stopifnot(ncol(xind) == nxgrid)
		if(nrow(xind)== 1){
			xind <- matrix(as.vector(xind), nrow=n, ncol=nxgrid, byrow=T) 
		} 
		stopifnot(nrow(xind) == n)  
	}	
	stopifnot(all.equal(order(xind[1,]), 1:nxgrid), all.equal(order(yind), 1:nygrid))
	
	basistype <- match.arg(basistype)
	integration <- match.arg(integration)
	
# scale xind to [0, 1] and check for reasonably equidistant gridpoints
	xind.sc <- xind - min(xind)
	xind.sc <- xind.sc/max(xind.sc)
	diffXind <- t(round(apply(xind.sc, 1, diff), 3))
	if(is.null(L) & any(apply(diffXind, 1, function(x) length(unique(x))) != 1) && # gridpoints for any  X_i(s) not equidistant?
			integration=="simpson"){
		warning("Non-equidistant grid detected for ", deparse(substitute(X)), ".\n Changing to trapezoidal rule for integration.")
		integration <- "trapezoidal"
	}
    if(!is.null(limits) && integration != "riemann"){
        integration <- "riemann" 
        warning("<limits>-argument detected. Changing to Riemann sums for numerical integration.")
    }
# FIXME: figure out weights for simpson's rule on non-equidistant grids instead of all this...
	
#make weight matrix for by-term
	if(!is.null(L)){
		stopifnot(nrow(L) == n, ncol(L) == nxgrid)
		#TODO: check whether supplied L is compatibel with limits argument
	} else {
		
			L <- switch(integration,
					"simpson" = {
						# \int^b_a f(t) dt = (b-a)/gridlength/3 * [f(a) + 4*f(t_1) + 2*f(t_2) + 4*f(t_3) + 2*f(t_3) +...+ f(b)]
						((xind[,nxgrid]-xind[,1])/nxgrid)/3 * 
								matrix(c(1, rep(c(4, 2), length=nxgrid-2), 1), nrow=n, ncol=nxgrid, byrow=T)
					}, 
					"trapezoidal" = {
						# \int^b_a f(t) dt = .5* sum_i (t_i - t_{i-1}) f(t_i) + f(t_{i-1}) = 
						#	(t_2 - t_1)/2 * f(a=t_1) + sum^{nx-1}_{i=2} ((t_i - t_i-1)/2 + (t_i+1 - t_i)/2) * f(t_i) + ... + 
						#			+ (t_nx - t_{nx-1})/2 * f(b=t_n)
						diffs <- t(apply(xind, 1, diff))
						.5 * cbind(diffs[,1], t(apply(diffs, 1, filter, filter=c(1,1)))[,-(nxgrid-1)], diffs[,(nxgrid-1)])
					},
                    "riemann" = {
                        # simple quadrature rule:
                        # \int^b_a f(t) dt = sum_i (t_i-t_{i-1})*(f(t_i)) 
                        diffs <- t(apply(xind, 1, diff))
                        #assume delta(t_0=a, t_1) = avg. delta
                        cbind(rep(mean(diffs),n), diffs)
                    } 
            )
		
	}
	LX <- L*X
    LX.stacked <- LX[rep(1:n, each=nygrid),]
    
    if(!is.null(limits)){
        if(!is.function(limits)){
            if(!(limits %in% c("s<t","s<=t"))){
                stop("supplied <limits> argument unknown")  
            }
            if(limits=="s<t"){
                limits <- function(s, t){
                    s < t
                }    
            } else {
                if(limits=="s<=t"){
                    limits <- function(s, t){
                        (s < t) | (s == t)
                    }    
                }   
            }
        }
        ind0 <- !t(outer(xind[1,], rep(yind, times=n), limits))
        LX.stacked[ind0] <- 0
    }    
    
    
	#expand for stacked Y-observations and assign unique names based on the given args
	xindname <- paste(deparse(substitute(X)), ".smat", sep="")
	yindname <- paste(deparse(substitute(X)), ".tmat", sep="")
	LXname <- paste("L.", deparse(substitute(X)), sep="")
	
	data <- list(
			xind[rep(1:n, times=nygrid), ], #stack xind nygrid-times
			matrix(rep(yind, times=n), nrow=n*nygrid, ncol=nxgrid), #repeat each entry of yind n times s.t. rows are constant  
            LX.stacked)
	names(data)  <- c(xindname, yindname, LXname)
	
# make call
	splinefun <- as.symbol(basistype) # if(basistype=="te") quote(te) else quote(s)
	frmls <- formals(getFromNamespace(deparse(splinefun), ns="mgcv"))
	frmls <- modifyList(frmls[names(frmls) %in% names(splinepars)], splinepars)
	call <- as.call(c(
					list(splinefun,
							x = as.symbol(substitute(xindname)), 
							z = as.symbol(substitute(yindname)),
							by =as.symbol(substitute(LXname))),
					frmls))
	
    
    if(check.ident){
        ## check whether (number of basis functions) < (number of relevant eigenfunctions of X)
        evls <- svd(X, nu=0, nv=0)$d^2
        evls[evls<0] <- 0
        maxK <- max(1, min(which((cumsum(evls)/sum(evls)) >= .99)))
        if(maxK <= 4) 
            warning("Very low effective rank of ", deparse(match.call()$X), " detected.",maxK," largest eigenvalues alone account for >99% of variability.")
        if(basistype!="s"){
            bsdim <- eval(call)$margin[[1]]$bs.dim
            if(maxK < bsdim){
                warning("<k> larger than effective rank. Adjusting marginal number of basis functions for ",deparse(match.call()$X)," downwards from ", bsdim, " to ", maxK,".")
                if(is.null(call$k)){
                    call$k <- eval(maxK) 
                } else {
                    call$k[1] <- eval(maxK)  
                }
            } 
        }
    }
    
    
	return(list(call=call, data = data, yind=yind, xind=xind[1, ], L=L,
                    xindname=xindname, yindname=yindname, LXname=LXname, limits=limits))
}#end ff()
