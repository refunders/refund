#' Construct a PC-based function-on-function regression term
#' 
#' Defines a term \eqn{\int X_i(s)\beta(t,s)ds} 
#' for inclusion in an \code{mgcv::gam}-formula (or \code{bam} or \code{gamm} or \code{gamm4:::gamm4}) as constructed
#' by \code{\link{pffr}}. 
#' 
#' In contrast to \code{\link{ff}}, \code{ffpc}  
#' does an FPCA decomposition \eqn{X(s) \approx \sum^K_{k=1} \xi_{ik} \Phi_k(s)} using \code{\link{fpca.sc}} and
#' represents \eqn{\beta(t,s)} in the function space spanned by these \eqn{\Phi_k(s)}.
#' That is, since 
#' \deqn{\int X_i(s)\beta(t,s)ds = \sum^K_{k=1} \xi_{ik} \int \Phi_k(s) \beta(s,t) ds = \sum^K_{k=1} \xi_{ik} \tilde \beta_k(t),}  
#' the function-on-function term can be represented as a sum of \eqn{K} univariate functions \eqn{\tilde \beta_k(t)} in \eqn{t} each multiplied by the FPC 
#' loadings \eqn{\xi_{ik}}. The truncation parameter \eqn{K} is chosen as described in \code{\link{fpca.sc}}.
#' Using this instead of \code{ff()} can be beneficial if the covariance operator of the \eqn{X_i(s)}
#' has low effective rank (i.e., if \eqn{K} is small). If the covariance operator of the \eqn{X_i(s)} 
#' is of (very) high rank, i.e., if \eqn{K} is large, \code{ffpc()} will not be very efficient.
#' 
#' To reduce model complexity, the \eqn{\tilde \beta_k(t)} all have a single joint smoothing parameter 
#' (in \code{mgcv}, they get the same \code{id}, see \code{\link[mgcv]{s}}).\cr 
#' 
#' Please see \code{\link[refund]{pffr}} for details on model specification and 
#' implementation. 
#' 
#' @param X an n by \code{ncol(xind)} matrix of function evaluations \eqn{X_i(s_{i1}),\dots, X_i(s_{iS})}; \eqn{i=1,\dots,n}.
#' @param yind matrix (or vector) of indices of evaluations of \eqn{Y_i(t)}
#' @param xind matrix (or vector) of indices of evaluations of \eqn{X_i(t)}, defaults to \code{seq(0, 1, length=ncol(X))}.
#' @param splinepars optional arguments supplied to the \code{basistype}-term. Defaults to a cubic 
#' 	B-spline with first difference penalties and 8 basis functions for each \eqn{\tilde \beta_k(t)}.
#' @param decomppars  parameters for the FPCA performed with \code{\link{fpca.sc}}.
#' @param npc.max maximal number \eqn{K} of FPCs to use, regardless of \code{decomppars}; defaults to 15
#' @return a list containing the necessary information to construct a term to be included in a \code{mgcv::gam}-formula.
#'  
#' @author Fabian Scheipl
#' @examples \dontrun{
#' set.seed(1122)
#' n <- 55
#' S <- 60
#' T <- 50
#' s <- seq(0,1, l=S)
#' t <- seq(0,1, l=T)
#' 
#' #generate X from a polynomial FPC-basis:
#' rankX <- 5
#' Phi <- cbind(1/sqrt(S), poly(s, degree=rankX-1))
#' lambda <- rankX:1
#' Xi <- sapply(lambda, function(l) 
#'             scale(rnorm(n, sd=sqrt(l)), scale=FALSE))
#' X <- Xi %*% t(Phi)
#' 
#' beta.st <- outer(s, t, function(s, t) cos(2 * pi * s * t))
#' 
#' y <- (1/S*X) %*% beta.st + 0.1 * matrix(rnorm(n * T), nrow=n, ncol=T)
#' 
#' data <- list(y=y, X=X)
#' # set number of FPCs to true rank of process for this example:
#' m.pc <- pffr(y ~ c(1) + 0 + ffpc(X, yind=t, decomppars=list(npc=rankX)), 
#'         data=data, yind=t)
#' summary(m.pc)
#' m.ff <- pffr(y ~ c(1) + 0 + ff(X, yind=t), data=data, yind=t)
#' summary(m.ff)
#' 
#' # fits are very similar:
#' all.equal(fitted(m.pc), fitted(m.ff))
#' 
#' # plot implied coefficient surfaces:
#' layout(t(1:3))
#' persp(t, s, t(beta.st), theta=50, phi=40, main="Truth", 
#'     ticktype="detailed")
#' plot(m.ff, select=1, pers=TRUE, zlim=range(beta.st), theta=50, phi=40,
#'     ticktype="detailed") 
#' title(main="ff()")
#' ffpcplot(m.pc, type="surf", auto.layout=FALSE, theta = 50, phi = 40)
#' title(main="ffpc()")
#' 
#' # show default ffpcplot:
#' ffpcplot(m.pc)
#' }
ffpc <- function(X,
        yind,
        xind=seq(0, 1, length=ncol(X)),
        splinepars=list(bs="ps", m=c(2, 1), k=8),
        decomppars=list(pve = .99, useSymm = TRUE),
        npc.max=15
){
    # check & format index for Y 
    if(!missing(yind))
        if(is.null(dim(yind))){
            yind <- t(t(yind))
        } 
    nygrid <- length(yind)
    nxgrid <- length(xind)
    
    # check & format index for X
    stopifnot(length(xind) == ncol(X))
    stopifnot(all.equal(order(xind), 1:nxgrid))
    
    decomppars$Y <- X
    klX <- do.call(fpca.sc, decomppars)
    xiMat <- klX$scores[,1:min(ncol(klX$scores),npc.max), drop=FALSE]
   
    #expand for stacked Y-observations and assign unique names based on the given args
    colnames(xiMat) <- paste(make.names(deparse(substitute(X))),".PC", 1:ncol(xiMat), sep="")
    xiMat <- xiMat[rep(1:nrow(xiMat), each=nygrid), ,drop=FALSE]
    id <- paste(make.names(deparse(substitute(X))),".ffpc", sep="")
    return(list(data=xiMat, PCMat=klX$efunctions[,1:min(ncol(klX$scores),npc.max), drop=FALSE],
                    meanX=klX$mu,
                    eigenvalues=klX$evalues, 
                    xind=xind, id=id, splinepars=splinepars))
}#end ffpc()

#' Plot PC-based function-on-function regression terms
#' 
#' Convenience function for graphical summaries of \code{ffpc}-terms from a 
#' \code{pffr} fit.  
#' 
#' @param object a fitted \code{pffr}-model
#' @param type one of "fpc+surf", "surf" or "fpc": "surf" shows a perspective plot of the coefficient surface implied 
#'          by the estimated effect functions of the FPC loadings, "fpc" shows three plots: 
#'          1) a scree-type plot of the estimated eigenvalues of the functional covariate, 2) the estimated eigenfunctions,
#'          and 3) the estimated coefficient functions associated with the FPC loadings. Defaults to showing both.
#' @param se.mult display estimated coefficient functions associated with the FPC loadings with plus/minus this number time the estimated standard error.
#' Defaults to 2. 
#' @param pages  the number of pages over which to spread the output. Defaults to 1. (Irrelevant if \code{auto.layout=FALSE}.)
#' @param ticktype see \code{\link[graphics]{persp}}.
#' @param theta see \code{\link[graphics]{persp}}.
#' @param phi see \code{\link[graphics]{persp}}.
#' @param plot Produce plots or only return plotting data? Defaults to \code{TRUE}. 
#' @param auto.layout should the the function set a suitable layout automatically? Defaults to TRUE
#' @return primarily produces plots, invisibly returns a list containing 
#' the data used for the plots.
#'  
#' @author Fabian Scheipl
#' @examples \dontrun{
#'  #see ?ffpc
#' }


ffpcplot <- function(object, type=c("fpc+surf", "surf", "fpc"), pages=1, 
        se.mult=2,  ticktype="detailed", theta = 30, phi = 30, plot=TRUE, 
        auto.layout=TRUE){
   
    type <- match.arg(type)
    T <- object$pffr$nyindex
    nterms <- length(object$pffr$ffpc)
    ffpcnames <- names(object$pffr$ffpc)
    
    
    betadata <- object$model[1:T,]
    betadata[, grep(".PC[[:digit:]]+$", colnames(betadata))] <- 1
    termsffpc <- predict.gam(object, newdata=betadata, type="iterms", se.fit=TRUE)
    
    betatilde <- termsffpc$fit[,grep(".PC[[:digit:]]+$", colnames(termsffpc$fit)), drop=FALSE]
    betatilde.se <- termsffpc$se.fit[,grep(".PC[[:digit:]]+$", colnames(termsffpc$se.fit)), drop=FALSE]
    betatilde.up <-  betatilde + se.mult*betatilde.se
    betatilde.lo <-  betatilde - se.mult*betatilde.se
    betatildemap <- lapply(ffpcnames, function(n) which(colnames(betatilde) %in% object$pffr$labelmap[[n]]))
    
    phibeta <- vector(length=nterms, mode="list")
    names(phibeta) <- sapply(object$pffr$ffpc, "[[", "id")
    for(i in 1:nterms){
        ### betatilde_k (t) = \int phi_k(s) beta(s,t) ds \approx 1/S t(Phi_k) %*% beta(s,t)
        ### --> beta(s,t) = S * Phi %*% [betatilde_1(t) ... betatilde_K(t)]
        phibeta[[i]] <- length(object$pffr$ffpc[[i]]$xind) * object$pffr$ffpc[[i]]$PCMat %*% t(betatilde[,betatildemap[[i]], drop=FALSE])
    }   
    
    
   
    if(plot){
        if(auto.layout){
            nplots <- switch(type,
                    "surf"=nterms,
                    "fpc"=3*nterms,
                    "fpc+surf"=4*nterms)
            
            #define layout 
            plotsperpage <- ceiling(nplots/pages)
            columns <- switch(type,
                    "surf"=plotsperpage,
                    "fpc"=3,
                    "fpc+surf"=4)
            layout(matrix(1:plotsperpage, ncol=columns, nrow=ceiling(plotsperpage/columns), byrow=TRUE))
        }
        
        for(i in 1:nterms){
            trm <- object$pffr$ffpc[[i]]
            if(type=="fpc+surf" | type=="fpc"){#
                npc <- ncol(trm$PCMat)
                plot(1:npc, trm$eigenvalues[1:npc], col=1:npc, xlab="FPC", type="b", pch=19,
                        ylab=paste0("Estimated eigenvalues: ", trm$id), bty="n")
                matplot(trm$xind, trm$PCMat[,npc:1], type="l", lty=1, col=npc:1, xlab = "", 
                        ylab = paste0("Estimated FPCs: ", trm$id), bty="n")
                matplot(object$pffr$yind, betatilde[,rev(betatildemap[[i]])], type="l", lty=1,
                        ylim = range(betatilde.up[,betatildemap[[i]]], betatilde.lo[,betatildemap[[i]]]), col=npc:1,
                        xlab = "", ylab = paste0("Effects of FPC loadings"), bty="n")
                abline(h=0, col="grey", lwd=.5)
                secol <- length(betatildemap[[i]])
                for(j in rev(betatildemap[[i]])){
                    polygon(cbind(x=c(object$pffr$yind, rev(object$pffr$yind))), 
                            y=c(betatilde.up[,j], rev(betatilde.lo[,j])), 
                            col=do.call(rgb, as.list(c(col2rgb(secol)/255,.1))), border=NA)
                    secol <- secol-1
                }
            }
            if(type=="fpc+surf" | type=="surf"){
              
                persp(trm$xind, object$pffr$yind,  z=phibeta[[i]],
                        theta=theta, phi=phi,
                        ticktype=ticktype, xlab="x.index", ylab="y.index", 
                        zlim=range(as.vector(phibeta[[i]])),
                        zlab=trm$id)
            }
        } #end for(i)
    } 
    invisible(list(betatilde=betatilde, betatilde.se=betatilde.se, phibeta=phibeta))
}#end ffpcplot()