# Utility functions for pffr:

safeDeparse <- function(expr){
	# turn an expression into a _single_ string, regardless of the expression's length 	
	ret <- paste(deparse(expr), collapse="")
	#rm whitespace
	gsub("[[:space:]][[:space:]]+", " ", ret)
}

#' Return call with all possible arguments
#' 
#' Return a call in which all of the arguments which were supplied or have presets are specified by their full names and their supplied or default values. 
#' 
#' @param definition a function. See \code{\link[base]{match.call}}.
#' @param call an unevaluated call to the function specified by definition. See \code{\link[base]{match.call}}.
#' @param expand.dots logical. Should arguments matching ... in the call be included or left as a ... argument? See \code{\link[base]{match.call}}.
#' @return An object of mode "\code{\link[base]{call}}". 
#' @author Fabian Scheipl
#' @export
#' @seealso \code{\link[base]{match.call}}
expand.call <-
        function(definition=NULL, call=sys.call(sys.parent(1)), expand.dots = TRUE)
{
    call <- match.call(definition, call, expand.dots)
    #given args:
    ans <- as.list(call)
    
    #possible args:
    frmls <- formals(safeDeparse(ans[[1]]))
    #remove formal args with no presets:
    frmls <- frmls[!sapply(frmls, is.symbol)]
    
    add <- which(!(names(frmls) %in% names(ans)))
    return(as.call(c(ans, frmls[add])))
}



list2df <- function(l){
# make a list into a dataframe -- matrices are left as matrices!
    nrows <- sapply(l, function(x) nrow(as.matrix(x)))
    stopifnot(length(unique(nrows)) == 1)
    ret <- data.frame(rep(NA, nrows[1]))
    for(i in 1:length(l)) ret[[i]] <- l[[i]]
    names(ret) <- names(l)
    return(ret)
}

## TODO: this does not always yield unique labels, e.g. if you have
##   s(g, bs="re") + s(g, bs="mrf", xt=somepenalty) 
getShrtlbls <- function(object){
# make short labels for display/coef-list, etc...
    
    labelmap <- object$pffr$labelmap
    
    ret <- sapply(names(unlist(labelmap)),
            function(x){
                #make a parseable expression for ffpc terms
                if(grepl("^ffpc", x)){
                    ffpcnumber <- gsub("(^.+))([0-9]+$)","\\2", x)
                    x <- gsub(")[0-9]+",")",x)
                   }
                exprx <- parse(text=x)
                
                #remove argument names
                x <- gsub("((?:[A-Za-z]+))(\\s)(=)(\\s)", "", x, perl=TRUE)
                
                #remove whitespace
                x <- gsub("\\s", "", x)
                
                #remove everything after and including first quoted argument
                if(any(regexpr(",\".*",x)[[1]]>0)) {
                    x <- gsub("([^\"]*)(,[c\\(]*\".*)(\\)$)", "\\1\\3", x, perl=TRUE) 
                }
                
                #remove everything after last variable:
                lstvrbl <- tail(all.vars(exprx),1)                
                x <- gsub(paste("(^.*?(?=",lstvrbl,")",lstvrbl,")(.*$)",sep=""), "\\1", x, perl=TRUE)
                
                #match braces
                openbr <- sum(grepl("\\(", strsplit(x, "")[[1]]))
                closebr <- sum(grepl("\\)", strsplit(x, "")[[1]]))
                if(openbr>closebr) x <- paste(c(x, rep(")", openbr-closebr)), sep="",collapse="") 
                
                #add number of PC for ffpc terms
                if(grepl("^ffpc", x)){
                    x <- paste(x, ffpcnumber, sep="")
                }
                return(x)
            })
    # correct labels for factor variables:
    if(any(sapply(labelmap, length) > 1 & !sapply(names(labelmap), function(x) grepl("^ffpc", x)))){
        which <- which(sapply(labelmap, length) > 1 & !sapply(names(labelmap), function(x) grepl("^ffpc", x)))
        inds <- c(0, cumsum(sapply(labelmap, length)))
        for(w in which){
            ret[(inds[w]+1):inds[w+1]] <- {
                lbls <- labelmap[[w]]
                bylevels <- sapply(object$smooth[lbls], function(x) x$by.level)
                by <- object$smooth[[lbls[1]]]$by
                paste(by, bylevels, "(", object$pffr$yindname, ")", sep="")
            }
        }
    }        
    #append labels for varying coefficient terms
    if(any(!grepl("\\(", ret))){
        which <- which(!grepl("\\(", ret))
        ret[which] <- paste(ret[which],"(", object$pffr$yindname, ")", sep="")
    }
    
    return(ret)
}

#' Simulate example data for pffr
#' 
#' Simulates example data for \code{\link{pffr}}.
#' Scenario 1 generates data from a simple univariate model 
#' \eqn{Y_i(t) = \mu(t) +\int X_i(s)\beta(s,t)ds + \epsilon_i(t)}.
#' Scenario 2 generates data from a complex multivariate model 
#' \eqn{Y_i(t) = \mu(t) + \int X_{1i}(s)\beta_1(s,t)ds + 
#' \int X_{2i}(s)\beta_2(s,t)ds + x_3 \beta_3(t) + 
#' f(x_4, x_5) + f(x6, t) + \beta_4 x_7 + \epsilon_i(t)}.
#' 
#' See source code for details.
#' 
#' @param scenario see description
#' @param n number of observations
#' @param nxgrid number of evaluation points of functional covariates  
#' @param nygrid number of evaluation points of the functional response
#' @param SNR the signal-to-noise ratio for the generated data:  
#'  empirical variance of the additive predictor divided by variance of the errors.
#' @return a named list with the simulated data, and the true components of the predictor etc as attributes.  
pffrSim <- function(
        scenario=1,
        n = 100,
        nxgrid = 40,
        nygrid = 60,
        SNR = 10
){
    ## generates random functions...
    rf <- function(x=seq(0,1,length=100), bs.dim=7, center=FALSE) {
      nk <- bs.dim - 2
      xu <- max(x)
      xl <- min(x)
      xr <- xu - xl
      xl <- xl - xr * 0.001
      xu <- xu + xr * 0.001
      dx <- (xu - xl)/(nk - 1)
      kn <- seq(xl - dx * 3, xu + dx * 3, 
          length = nk + 4 + 2) 
      X <- splines::spline.des(kn, x, 4, x * 0)$design
      
      drop(X %*% rnorm(bs.dim))
    } 
    
    test1 <- function(s, t){
      s*cos(pi*abs(s-t)) - .19
    }
    test2 <- function(s, t, ss=0.3, st=0.4)
    { 
      cos(pi*s)*sin(pi*t) + (s*t)^2 - 0.11
    }
    s <- seq(0, 1, length=nxgrid)
    t <- seq(0, 1, length=nygrid)
    
    
    mu.t <- matrix(dbeta(t, 2,7), nrow=n, ncol=nygrid, byrow=TRUE)
    
    data <- list()
    etaTerms <- list()
    etaTerms$int <- mu.t
    
    #functional covariates
    data$X1 <- I(t(replicate(n, rf(s))))
    
    L <- matrix(1/nxgrid, ncol=nxgrid, nrow=n)
    LX1 <- L*data$X1 
    beta1.st <- outer(s, t, test1)
    etaTerms$X1 <- LX1%*%beta1.st
    
    data$X2 <- I(t(replicate(n, rf(s))))
    LX2 <- L*data$X2 
    beta2.st <- outer(s, t, test2)
    etaTerms$X2 <- LX2%*%beta2.st
    
    #scalar covariates
    data$xlin <- I(rnorm(n))
    beta.t <- matrix(scale(-dnorm(4*(t-.2))), nrow=n, ncol=nygrid, byrow=T)
    etaTerms$xlin <- data$xlin*beta.t
    
    data$xsmoo <- I(rnorm(n))
    etaTerms$xsmoo <- outer(drop(scale(cos(data$xsmoo))), (t-.5), "*")
    
    data$xte1 <- I(rnorm(n))
    data$xte2 <- I(rnorm(n))
    etaTerms$xte <- drop(scale(-data$xte1*data$xte2^2))
    
    data$xconst <- I(rnorm(n))
    etaTerms$xconst <- 2*data$xconst
    
    eta <- mu.t + switch(scenario,
            "1" = etaTerms$X1,
            "2" = Reduce("+", etaTerms))
    
    eps <-  sd(as.vector(eta))/sqrt(SNR) * matrix(scale(rnorm(n*nygrid)), nrow=n)
    data$Y <- I(eta + eps)
    
    return(structure(as.data.frame(data), xindex=s, yindex=t, 
                    truth=list(eta=eta, etaTerms=etaTerms), call=match.call()))
}