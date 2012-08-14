#' pffr-constructor for PC-basis functional random effects.
#' 
#' Fits functional random effects \eqn{B_i(t)} for a grouping variable \code{id} 
#' using as a basis the eigenfunctions in \code{Phi} with eigenvalues \code{lambda}:
#' \eqn{B_i(t) \approx \sum_m^M \phi_m(t)\delta_{im}} with  
#' independent \eqn{\delta_{im} \sim N(0, \sigma^2\lambda_m)}, where \eqn{\sigma^2}
#' is estimated and controls the overall contribution of the \eqn{B_i(t)} while the relative importance
#' of the \eqn{M} eigenfunctions is controlled by the supplied eigenvalues \code{lambda_m}.
#' Can be used to model smooth residuals if \code{id} is simply an index of observations. This is an experimental feature and not
#' well tested yet, use at your own risk.   
#' 
#' @param id grouping variable a factor
#' @param Phi matrix of eigenfunction evaluations on gripoints \code{yind} (<length of \code{yind}> x <no. of used eigenfunctions>)
#' @param lambda eigenvalues associated with \code{Phi}
#' @param yind vector of gridpoints on which responses \eqn{Y(t)} are evaluated.
#' @return a list used internally for constructing an appropriate call to \code{mgcv::gam} 
pcre <- function(id, 
        efunctions,
        evalues,
        yind
){
    
    # check args
    stopifnot(is.factor(id), nrow(efunctions)==length(yind), ncol(efunctions)==length(evalues), all(evalues>0))
    
    nygrid <- length(yind)
    
    efunctionsname <- deparse(substitute(efunctions))
    idname <- paste(deparse(substitute(id)),".vec",sep="")
    
    #scale eigenfunctions by their eigenvalues:
    efunctions <- t(t(efunctions)*sqrt(evalues))
    
    #expand for stacked Y-observations and assign unique names based on the given args
    colnames(efunctions) <- paste(efunctionsname,".PC", 1:ncol(efunctions), sep="")
    efunctionsmat <- efunctions[rep(1:nrow(efunctions), times=length(id)), ]
    
    idvec <- id[rep(1:length(id), each=nygrid)]
    
    
    data <- data.frame(id=idvec, efunctions=efunctionsmat)
    names(data) <- c(idname, colnames(efunctions)) 
    
    call <- as.call(c(as.symbol("s"),
                    as.symbol(substitute(idname)),
                    sapply(colnames(efunctions), function(x) as.symbol(x)),
                    bs=c("pcre")))
    
    return(list(data=data, efunctions=efunctions, yind=yind, id=id, call=call))
}#end pcre()
