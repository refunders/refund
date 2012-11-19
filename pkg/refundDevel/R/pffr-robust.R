#' Penalized function-on-function regression with non-i.i.d. residuals
#' 
#' THIS IS AN EXPERIMENTAL VERSION AND NOT WELL TESTED YET -- USE AT YOUR OWN RISK.\cr
#' Implements additive regression for functional and scalar covariates and functional responses.
#' This function is a wrapper for \code{mgcv}'s \code{\link[mgcv]{gam}} and its siblings to fit models of the general form \cr
#' \eqn{Y_i(t) = \mu(t) + \int X_i(s)\beta(s,t)ds + f(z_{1i}, t) + f(z_{2i}) + z_{3i} \beta_3(t) + \dots  + E_i(t))}\cr
#' with a functional (but not necessarily continuous) response \eqn{Y(t)}, 
#' (optional) smooth intercept \eqn{\mu(t)}, (multiple) functional covariates \eqn{X(t)} and scalar covariates
#' \eqn{z_1}, \eqn{z_2}, etc. The residual functions \eqn{E_i(t) \sim GP(0, K(t,t'))} are assumed to be i.i.d. 
#' realizations of a Gaussian process. An estimate of the covariance operator \eqn{K(t,t')} evaluated on \code{yind} 
#' has to be supplied in the \code{hatSigma}-argument.
#' Note that this has to be positive definite. If \code{hatSigma} is close to positive \emph{semi-}definite or badly conditioned, 
#' estimated standard errors become unstable (typically much too small). \code{pffrGLS} will try to diagnose this and issue a warning.
#' The danger is especially big if the number of functional observations is smaller than the number of gridpoints 
#' (i.e, \code{length(yind)}), since the raw covariance estimate will not have full rank.     
#' 
#' @section Details: 
#' Please see \code{\link[refund]{pffr}} for details on model specification and 
#' implementation. \cr THIS IS AN EXPERIMENTAL VERSION AND NOT WELL TESTED YET -- USE AT YOUR OWN RISK.
#'  
#' @param formula a formula with special terms as for \code{\link[mgcv]{gam}}, with additional special terms \code{\link{ff}()} and \code{c()}. See \code{\link[refund]{pffr}}.
#' @param yind a vector with length equal to the number of columns of the matrix of functional responses giving the vector of evaluation points \eqn{(t_1, \dots ,t_{G})}.
#'   see \code{\link[refund]{pffr}}
#' @param algorithm the name of the function used to estimate the model. Defaults to \code{\link[mgcv]{gam}} if the matrix of functional responses has less than \code{2e5} data points
#' 	 and to \code{\link[mgcv]{bam}} if not. "gamm" (see \code{\link[mgcv]{gamm}}) and "gamm4" (see \code{\link[gamm4]{gamm4}}) are valid options as well.  
#' @param hatSigma (an estimate of) the within-observation covariance (along the responses' index), evaluated at \code{yind}. See Details.
#' @param method See \code{\link[refund]{pffr}}
#' @param bs.yindex See \code{\link[refund]{pffr}}
#' @param bs.int See \code{\link[refund]{pffr}}
#' @param tensortype See \code{\link[refund]{pffr}}
#' @param ... additional arguments that are valid for \code{\link[mgcv]{gam}} or \code{\link[mgcv]{bam}}. See \code{\link[refund]{pffr}}.
#' 
#' @return a fitted \code{pffr}-object, see \code{\link[refund]{pffr}}.
#' @seealso \code{\link[refund]{pffr}}, \code{\link[refund]{fpca.sc}}
#' @author Fabian Scheipl
pffrGLS <- function(
        formula,
        yind,
        hatSigma,
        algorithm = NA, 
        method="REML",
        tensortype = c("te", "t2"),
        bs.yindex = list(bs="ps", k=5, m=c(2, 1)), # only bs, k, m are propagated...
        bs.int = list(bs="ps", k=20, m=c(2, 1)), # only bs, k, m are propagated...
        ...
){
    # TODO: need check whether yind supplied in ff-terms and pffr are identical!  
    # TODO: allow term-specific overrides of bs.yindex    
    # TODO: weights, subset, offset args!
    # TODO: write missing locations into pffr so fitted etc will work...    
    
    call <- match.call()
    tensortype <- as.symbol(match.arg(tensortype))
    
    ## warn if any entries in ... are not arguments for gam/gam.fit or gamm4/lmer 
    dots <- list(...)
    if(length(dots)){
        validDots <- if(!is.na(algorithm) && algorithm=="gamm4"){
                    c(names(formals(gamm4)), names(formals(lmer)))         
                } else {
                    c(names(formals(gam)), names(formals(gam.fit)))
                }
        notUsed <- names(dots)[!(names(dots) %in% validDots)]
        if(length(notUsed))
            warning("Arguments <", paste(notUsed, collapse=", "), "> supplied but not used." )
    }
    
    
    tf <- terms.formula(formula, specials=c("s", "te", "t2", "ff", "c", "sff", "ffpc", "pcre"))
    trmstrings <- attr(tf, "term.labels")
    terms <- sapply(trmstrings, function(trm) as.call(parse(text=trm))[[1]], simplify=FALSE) 
    #ugly, but getTerms(formula)[-1] does not work for terms like I(x1:x2) 
    frmlenv <- environment(formula)
    
    
    where.c <- attr(tf, "specials")$c - 1    # indices of scalar offset terms
    where.ff <- attr(tf, "specials")$ff - 1  # function-on-function terms
    where.sff <- attr(tf, "specials")$sff - 1  #smooth function-on-function terms
    where.s <- attr(tf, "specials")$s - 1    # smooth terms 
    where.te <- attr(tf, "specials")$te - 1  # tensor product terms
    where.t2 <- attr(tf, "specials")$t2 - 1  # type 2 tensor product terms
    where.ffpc <- attr(tf, "specials")$ffpc - 1  # PC-based function-on-function terms
    where.pcre <- attr(tf, "specials")$pcre - 1  # functional random effects/residuals with PC-basis 
    if(length(trmstrings)) {
        where.par <- which(!(1:length(trmstrings) %in%
                            c(where.c, where.ff, where.sff, where.ffpc, where.pcre, where.s, where.te, where.t2))) # indices of linear/factor terms with varying coefficients over yind.
    } else where.par <- numeric(0)
    
    responsename <- attr(tf,"variables")[2][[1]]
    
    #start new formula
    newfrml <- paste(responsename, "~", sep="")
    newfrmlenv <- new.env()
    evalenv <- if("data" %in% names(call)) eval.parent(call$data) else NULL
    
    nobs <- nrow(eval(responsename,  envir=evalenv, enclos=frmlenv))
    nyindex <- ncol(eval(responsename,  envir=evalenv, enclos=frmlenv))
    
    
    if(missing(algorithm)||is.na(algorithm)){
        algorithm <- ifelse(nobs > 1e5, "bam", "gam")
    } 
    algorithm <- as.symbol(algorithm)
    if(as.character(algorithm)=="bam" && !("chunk.size" %in% names(call))){
        call$chunk.size <- max(nobs/5, 10000) 
        #same default as in bam
    }
    ## 
    ###<GLS
    if(as.character(algorithm)=="gamm4") stop("pffrGLS not implemented for gamm4")
    ###GLS>
    
    #if missing, define y-index or get it from first ff/sff-term, then assign expanded versions to newfrmlenv
    if(missing(yind)){
        if(length(c(where.ff, where.sff))){
            if(length(where.ff)){
                ffcall <- expand.call(ff, as.call(terms[where.ff][1])[[1]])  
            }  else ffcall <- expand.call(sff, as.call(terms[where.sff][1])[[1]]) 
            if(!is.null(ffcall$yind)){
                yind <- eval(ffcall$yind, envir=evalenv, enclos=frmlenv)
                yindname <- deparse(ffcall$yind)
            } else {
                yind <- 1:nyindex
                yindname <- "yindex"	
            }			
        } else {
            yind <- 1:nyindex
            yindname <- "yindex"	
        }
    } else {
        stopifnot(is.vector(yind), is.numeric(yind), 
                length(yind) == nyindex)
        yindname <- deparse(substitute(yind))
    }
    #make sure it's a valid name 
    if(length(yindname)>1) yindname <- "yindex"	
    # make sure yind is sorted
    stopifnot(all.equal(order(yind), 1:nyindex))
    
    
    yindvec <- rep(yind, times = nobs)
    yindvecname <- as.symbol(paste(yindname,".vec",sep=""))
    assign(x=deparse(yindvecname), value=yindvec, envir=newfrmlenv)
    
    
    ###<GLS
    sqrtSigmaInv <- {
        eSigma <- eigen(hatSigma, symmetric=TRUE)
        cond <- max(eSigma$values)/min(eSigma$values)
        if(cond > 1e5){
            diag(hatSigma) <- 1.05*diag(hatSigma)
            eSigma <- eigen(hatSigma, symmetric = TRUE) 
            condnew <- max(eSigma$values)/min(eSigma$values)
            warning("Supplied hatSigma badly conditioned with condition number ", round(cond),
                    "\n   -- inflated diagonal for new working covariance matrix (new condition number: ", round(condnew),").")
        }
        with(eSigma, vectors%*%diag(1/sqrt(values))%*%t(vectors))
    } 
    
    hatSigmaname <- deparse(substitute(hatSigmaname))
    assign(x=deparse(hatSigmaname), value=hatSigma, envir=frmlenv)
    
   
    originalresponsename <- paste(deparse(responsename),".Original", sep="")
    assign(x=originalresponsename, 
            value=as.vector(t(eval(responsename, envir=evalenv, enclos=frmlenv))), 
            envir=newfrmlenv)
    ## 'decorrelate' Y
    ytilde <- sqrtSigmaInv%*%t(eval(responsename, envir=evalenv, enclos=frmlenv))
    #assign (decorrelated) response in _long_ format to newfrmlenv
    assign(x=deparse(responsename), value=as.vector(ytilde), 
            envir=newfrmlenv)
    
    
    if(any(is.na(get(as.character(responsename), newfrmlenv)))){
        stop("no missings in y for GLS version of pffr allowed")
    }
    ####GLS>
    
    ##################################################################################
    #modify formula terms.... 
    newtrmstrings <- attr(tf, "term.labels")
    
    #if intercept, add \mu(yindex)
    if(attr(tf, "intercept")){ 
        # have to jump thru some hoops to get bs.yindex handed over properly 
        # without having yind evaluated within the call    
        arglist <- c(name="s", x = as.symbol(yindvecname), bs.int)
        intcall <- NULL
        assign(x= "intcall", value= do.call("call", arglist, envir=newfrmlenv), envir=newfrmlenv)
        newfrmlenv$intcall$x <- as.symbol(yindvecname)  
        
        intstring <- deparse(newfrmlenv$intcall)
        rm(intcall, envir=newfrmlenv)
        
        newfrml <- paste(newfrml, intstring, sep=" ")
        addFint <- TRUE
        names(intstring) <- paste("Intercept(",yindname,")",sep="") 
    } else{
        newfrml <-paste(newfrml, "0", sep="")
        addFint <- FALSE
    } 
    
    #transform: c(foo) --> foo
    if(length(where.c)){ 
        newtrmstrings[where.c] <- sapply(trmstrings[where.c], function(x){
                    sub("\\)$", "", sub("^c\\(", "", x)) #c(BLA) --> BLA
                })
    }
    
    #prep function-on-function-terms
    if(length(c(where.ff, where.sff))){ 
        ffterms <- lapply(terms[c(where.ff, where.sff)], function(x){
                    eval(x, envir=evalenv, enclos=frmlenv)
                })
        
        newtrmstrings[c(where.ff, where.sff)] <- sapply(ffterms, function(x) {
                    safeDeparse(x$call)
                })
        #assign newly created data to newfrmlenv
        lapply(ffterms, function(x){
                    lapply(names(x$data), function(nm){
                                assign(x=nm, value=x$data[[nm]], envir=newfrmlenv)
                                invisible(NULL)
                            })
                    invisible(NULL)
                })  
        ffterms <- lapply(ffterms, function(x) x[names(x)!="data"])
    } else ffterms <- NULL
    if(length(where.ffpc)){
        ffpcterms <- lapply(terms[where.ffpc], function(x){
                    eval(x, envir=evalenv, enclos=frmlenv)
                })
        #assign newly created data to newfrmlenv
        lapply(ffpcterms, function(trm){
                    lapply(colnames(trm$data), function(nm){
                                assign(x=nm, value=trm$data[,nm], envir=newfrmlenv)
                                invisible(NULL)
                            })
                    invisible(NULL)
                })  
        
        
        getFfpcFormula <- function(trm) {
            frmls <- lapply(colnames(trm$data), function(pc) {
                        arglist <- c(name="s", x = as.symbol(yindvecname), by= as.symbol(pc), id=trm$id, trm$splinepars)
                        call <- do.call("call", arglist, envir=newfrmlenv)
                        call$x <- as.symbol(yindvecname)
                        call$by <- as.symbol(pc)
                        safeDeparse(call)
                    })
            return(paste(unlist(frmls), collapse=" + "))
        }
        newtrmstrings[where.ffpc] <- sapply(ffpcterms, getFfpcFormula)
        
        ffpcterms <- lapply(ffpcterms, function(x) x[names(x)!="data"])
    } else ffpcterms <- NULL
    
    #prep PC-based random effects
    if(length(where.pcre)){
        pcreterms <- lapply(terms[where.pcre], function(x){
                    eval(x, envir=evalenv, enclos=frmlenv)
                })
        #assign newly created data to newfrmlenv
        lapply(pcreterms, function(trm){
                    lapply(colnames(trm$data), function(nm){
                                assign(x=nm, value=trm$data[,nm], envir=newfrmlenv)
                                invisible(NULL)
                            })
                    invisible(NULL)
                })  
        
        newtrmstrings[where.pcre] <- sapply(pcreterms, function(x) {
                    safeDeparse(x$call)
                })
        
        pcereterms <- lapply(pcreterms, function(x) x[names(x)!="data"])
    }else pcreterms <- NULL
    
    #transform: s(x, ...), te(x, z,...), t2(x, z, ...) --> te(x, (z,)  yindex, ..., <bs.yindex>)
    makeSTeT2 <- function(x){
        
        xnew <- x
        if(deparse(x[[1]]) == "te" && as.character(algorithm) == "gamm4") xnew[[1]] <- quote(t2)
        if(deparse(x[[1]]) == "s"){
            xnew[[1]] <- if(as.character(algorithm) != "gamm4") {
                        tensortype 
                    } else quote(t2)
            #accomodate multivariate s()-terms
            xnew$d <- if(!is.null(names(xnew))){
                        c(length(all.vars(xnew[names(xnew)==""])), 1)
                    } else c(length(all.vars(xnew)), 1)  
        } else {
            if("d" %in% names(x)){ #either expand given d...
                xnew$d <- c(eval(x$d), 1)
            } else {#.. or default to univariate marginal bases
                xnew$d <- rep(1, length(all.vars(x))+1)
            }
        }   
        xnew[[length(xnew)+1]] <- yindvecname
        this.bs.yindex <- if("bs.yindex" %in% names(x)){
                    x$bs.yindex
                } else bs.yindex
        xnew <- xnew[names(xnew) != "bs.yindex"]
        
        xnew$bs <- if("bs" %in% names(x)){
                    if("bs" %in% names(this.bs.yindex)){
                        c(eval(x$bs), this.bs.yindex$bs)
                    } else {
                        c(xnew$bs, "tp")
                    }
                } else {
                    if("bs" %in% names(this.bs.yindex)){
                        c(rep("tp", length(xnew$d)-1), this.bs.yindex$bs)
                    } else {
                        rep("tp", length(all.vars(x))+1)
                    }   
                }
        xnew$m <- if("m" %in% names(x)){
                    if("m" %in% names(this.bs.yindex)){
                        warning("overriding bs.yindex for m in ", deparse(x))
                    }
                    #TODO: adjust length if necessary, m can be a list for bs="ps","cp","ds"!
                    x$m
                } else {
                    if("m" %in% names(this.bs.yindex)){
                        this.bs.yindex$m
                    } else {
                        NA
                    }   
                }
        #defaults to 8 basis functions 
        xnew$k <- if("k" %in% names(x)){
                    if("k" %in% names(this.bs.yindex)){
                        c(xnew$k, this.bs.yindex$k)
                    } else {
                        c(xnew$k, 8)
                    } 
                } else {
                    if("k" %in% names(this.bs.yindex)){
                        c(pmax(8, 5^head(xnew$d, -1)), this.bs.yindex$k)
                    } else {
                        pmax(8, 5^xnew$d)
                    }   
                }
        
        xnew$xt <- if("xt" %in% names(x)){
                    
                    add.impose <- function(lst){
                        # use this in order to propagate xt-args to gam WITHOUT evaluating them,
                        # because this can break (stupid parse-cutoff!) and 
                        # shows up very ugly in summary etc for stuff like large penalty matrices    
                        for(i in 2:length(lst)){
                            llst <- length(lst[[i]])
                            if(llst){
                                lst[[i]][[llst+1]] <- TRUE
                                names(lst[[i]])[length(names(lst[[i]]))] <- "impose.ffregC"
                                lst[[i]][[llst+2]] <- nobs
                                names(lst[[i]])[length(names(lst[[i]]))] <- "nobs"
                            } else {
                                lst[[i]] <- bquote(list(impose.ffregC = .(TRUE), 
                                                nobs=.(nobs)))
                            }
                        }    
                        return(lst)
                    }
                    # xt has to be supplied as a list, with length(x$d) entries,
                    # each of which is a list or NULL:
                    stopifnot(x$xt[[1]]==as.symbol("list") && 
                                    length(x$xt)==length(xnew$d) &&  # =length(x$d)+1, since first element in parse tree is ``list''
                                    all(sapply(2:length(x$xt), function(i) 
                                                        x$xt[[i]][[1]] == as.symbol("list") ||
                                                                is.null(eval(x$xt[[i]][[1]])))))
                    xtra <- add.impose(x$xt)
                    xtra[[length(xnew$d)+1]] <- bquote(list(impose.ffregC = .(TRUE), 
                                    nobs=.(nobs)))
                    xtra
                } else {
                    imposearg <- bquote(list(impose.ffregC = .(TRUE), 
                                    nobs=.(nobs)))
                    
                    xtra <- bquote(list(.(imposearg)))
                    for (i in 3:(length(xnew$d)+1)) 
                        xtra[[i]] <- imposearg 
                    xtra
                }
        
        ret <- safeDeparse(xnew)
        return(ret)
    }
    
    if(length(c(where.s, where.te, where.t2))){ 
        newtrmstrings[c(where.s, where.te, where.t2)] <- 
                sapply(terms[c(where.s, where.te, where.t2)], makeSTeT2)
    }
    
    #transform: x --> s(YINDEX, by=x)
    if(length(where.par)){ 
        newtrmstrings[where.par] <- sapply(terms[where.par], function(x){
                    xnew <- bs.yindex
                    xnew <- as.call(c(quote(s), yindvecname, by=x, xnew))
                    safeDeparse(xnew)
                })
        
    }
    
    #... & assign expanded/additional variables to newfrmlenv
    where.notff <- c(where.c, where.par, where.s, where.te, where.t2)
    if(length(where.notff)){
        if("data" %in% names(call)) frmlenv <- list2env(eval.parent(call$data), frmlenv) 
        lapply(terms[where.notff], 
                function(x){
                    #nms <- all.vars(x)
                    if(any(unlist(lapply(terms[where.c], function(s) 
                                                if(length(s)==1){
                                                    s==x
                                                } else {
                                                    s[[1]]==x
                                                })))){
                        # drop c()
                        # FIXME: FUGLY! 
                        x <- formula(paste("~", gsub("\\)$", "",
                                                gsub("^c\\(", "", deparse(x)))))[[2]]
                    } 
                    ## remove names in xt, k, bs,  information (such as variable names for MRF penalties etc)
                    nms <- if(!is.null(names(x))){
                                all.vars(x[names(x) == ""]) 
                            }  else all.vars(x)
                    
                    
                    sapply(nms, function(nm){
                                stopifnot(length(get(nm, envir=frmlenv)) == nobs)
                                assign(x=nm, 
                                        value=rep(get(nm, envir=frmlenv), each=nyindex),
                                        envir=newfrmlenv)	
                                invisible(NULL)
                            })
                    invisible(NULL)
                })
    }
    newfrml <- formula(paste(c(newfrml, newtrmstrings), collapse="+"))
    environment(newfrml) <- newfrmlenv
    
    pffrdata <- list2df(as.list(newfrmlenv))
    
    newcall <- expand.call(pffr, call)
    newcall$yind <- newcall$tensortype <- newcall$bs.int <-
            newcall$bs.yindex <- newcall$algorithm <- NULL
    newcall$formula <- newfrml
    newcall$data <- quote(pffrdata)
    newcall[[1]] <- algorithm
    
    
    # add appropriate centering constraints for smooth effects
    # (not possible for gamm4, as the constraint destroys the simple 
    # diagonal structure of the t2-penalties)
    if(!(as.character(algorithm) %in% c("gamm4"))){
        suppressMessages(
                trace(mgcv:::smooth.construct.tensor.smooth.spec, 
                        at = max(which(sapply(as.list(body(mgcv:::smooth.construct.tensor.smooth.spec)), function(x) any(grepl(x, pattern="object$C", fixed=TRUE))))) + 1, 
                        print=FALSE,
                        tracer = quote({
                                    #browser()
                                    if(!is.null(object$margin[[length(object$margin)]]$xt$impose.ffregC)){
                                        ## constraint: sum_i f(z_i, t) = 0 \forall t
                                        cat("imposing constraints..\n")
                                        nygrid <- length(unique(data[[object$margin[[length(object$margin)]]$term]]))
                                        
                                        
                                        ## C = ((1,...,1) \otimes I_G) * B
                                        Ctmp <- kronecker(t(rep(1, object$margin[[length(object$margin)]]$xt$nobs)), diag(nygrid))
                                        if(ncol(Ctmp) > nrow(object$X)){
                                            #drop rows for missing obs.
                                            Ctmp <- Ctmp[, - get(".PFFRmissingResponses", .GlobalEnv)]
                                        }
                                        C <- Ctmp %*% object$X
                                        
                                        ## we need the number of effective constraints <nC> to correspond to the
                                        ## rank of the constraint matrix C, which is the min. of number of basis
                                        ## functions for t (object$margin[[length(object$margin)]]$bs.dim) and
                                        ## timepoints (<nygrid>)
                                        nC <- min(nygrid, object$margin[[length(object$margin)]]$bs.dim)
                                        object$C <- object$Cp <- C[seq(1, nygrid, length.out = nC), ]
                                    }}))
        )
        
        suppressMessages(
                trace(mgcv:::smooth.construct.t2.smooth.spec, 
                        at = max(which(sapply(as.list(body(mgcv:::smooth.construct.t2.smooth.spec)), function(x) any(grepl(x, pattern="object$Cp", fixed=TRUE))))) + 1, 
                        print=FALSE,
                        tracer = quote({
                                    if(!is.null(object$margin[[length(object$margin)]]$xt$impose.ffregC) &&
                                            object$margin[[length(object$margin)]]$xt$impose.ffregC){
                                        ## constraint: sum_i f(z_i, t) = 0 \forall t
                                        cat("imposing constraints..\n")
                                        nygrid <- length(unique(data[[object$margin[[length(object$margin)]]$term]]))
                                        
                                        ## C = ((1,...,1) \otimes I_G) * B
                                        Ctmp <- kronecker(t(rep(1, object$margin[[length(object$margin)]]$xt$nobs)), diag(nygrid))
                                        if(ncol(Ctmp) > nrow(object$X)){
                                            #drop rows for missing obs.
                                            Ctmp <- Ctmp[, - get(".PFFRmissingResponses", .GlobalEnv)]
                                        }
                                        C <- Ctmp %*% object$X
                                        
                                        ## we need the number of effective constraints <nC> to correspond to the
                                        ## rank of the constraint matrix C, which is the min. of number of basis
                                        ## functions for t (object$margin[[length(object$margin)]]$bs.dim) and
                                        ## timepoints (<nygrid>)
                                        nC <- min(nygrid, object$margin[[length(object$margin)]]$bs.dim)
                                        object$C <- object$Cp <- C[seq(1, nygrid, length.out = nC), ]
                                    }}))
        )     
        
        
        on.exit({
                    suppressMessages(try(untrace(mgcv:::smooth.construct.tensor.smooth.spec), silent = TRUE))
                    suppressMessages(try(untrace(mgcv:::smooth.construct.t2.smooth.spec), silent = TRUE))
                })
    }
    
    
    ###<GLS
    newcall$GLSinfo <- list(yind=yind, nobs=nobs, sqrtSigmaInv=sqrtSigmaInv)
    if(as.character(algorithm) == "gam"){
        trace(gam, at=4, quote({
                            #cat("premultiply Design...\n")
                            GLSinfo <- list(...)$GLSinfo
                            # modify design matrix: premultiply submatrix for each obs by sqrt(hatSigma)^-1
                            obsvec <- seq(1, GLSinfo$nobs*length(GLSinfo$yind), by=length(GLSinfo$yind))
                            from <- 1
                            to <- length(GLSinfo$yind)
                            for(i in 1:GLSinfo$nobs){
                                G$X[from:to, ] <- GLSinfo$sqrtSigmaInv%*%G$X[from:to, ]
                                from <- from + length(GLSinfo$yind)
                                to <- to + length(GLSinfo$yind)
                            }
                        }))
        on.exit(untrace(gam), add=TRUE)
    }
    
    if(as.character(algorithm) == "bam"){
        trace(bam, at=40, quote({
                            #browser()
                            cat("premultiply Design...\n")
                            GLSinfo <- list(...)$GLSinfo
                             # modify design matrix: premultiply submatrix for each obs by sqrt(hatSigma)^-1
                            obsvec <- seq(1, GLSinfo$nobs*length(GLSinfo$yind), by=length(GLSinfo$yind))
                            from <- 1
                            to <- length(GLSinfo$yind)
                            for(i in 1:GLSinfo$nobs){
                                G$X[from:to, ] <- GLSinfo$sqrtSigmaInv%*%G$X[from:to, ]
                                from <- from + length(GLSinfo$yind)
                                to <- to + length(GLSinfo$yind)
                            }  
                        }))
        on.exit(untrace(bam), add=TRUE)
    }
    #browser()
    
    m <- eval(newcall)
    # browser()
     
    # overwrite decorrelated response, fitted values etc s.t. summary etc are (more) correct
    m$y <- newfrmlenv[[originalresponsename]]
    # check that we really fitted the object 
    # (need this for coefboot.pffr), if yes overwrite
    if(is.null(newcall$fit) | (!is.null(newcall$fit) && eval(newcall$fit))){
        #browser()
        
        if(!is.null(m$model)){
            m$model[, deparse(responsename)] <- newfrmlenv[[originalresponsename]]
        }
        m$fitted <- m$fitted.values <- predict(m)
        m$residuals <- m$fitted - m$y
        # m$scale <- m$sig2 <- var(m$residuals)
        # TODO: fix Vp (drop Ve): 
        # Vp= (XWX + S^-1)^-1; Ve=(Vp^-1) XWX (Vp^-1)
#        X <- predict(m, type="lpmatrix")
#        XWX <- {
#            sqrtw <- if(is.null(m$prior.weights)){
#                    rep(1, length(m$y))             
#				}  else {
#                    sqrt(m$prior.weights)
#                }
#            as.matrix(crossprod((Diagonal(length(sqrtw), sqrtw)%*%
#                                kronecker(Diagonal(nobs), sqrtSigmaInv))%*%X))
#        }
#        Vp <- solve(XWX)
    
    
#        XtildeWXtilde <- {
#           xcall <- newcall
#           xcall$fit <- FALSE
#           sqrtw <- if(is.null(m$weights)){
#                    rep(1, length(m$y))             
#				}  else {
#                    sqrt(m$weights)
#                }
#           crossprod(diag(sqrtw)%*%eval(newcall)$X)
#        }
#        X <- predict(m, type="lpmatrix")
#        all.equal(as.matrix(crossprod(kronecker(Diagonal(nobs), solve(sqrtSigma))%*%X)), XtildeWXtilde)
    }
    
    
    #browser()
    
    ###GLS>
    
    
    m.smooth <- if(as.character(algorithm) %in% c("gamm4","gamm")){
                m$gam$smooth
            } else m$smooth
    
    #return some more info s.t. custom predict/plot/summary will work
    trmmap <- newtrmstrings
    names(trmmap) <- names(terms)
    if(addFint) trmmap <- c(trmmap, intstring)
    
    # map labels to terms -- 
    # ffpc are associated with multiple smooths
    # parametric are associated with multiple smooths if covariate is a factor
    labelmap <- as.list(trmmap)
    lbls <- sapply(m.smooth, function(x) x$label)
    if(length(c(where.par, where.ffpc))){
        if(length(where.par)){
            for(w in where.par) 
                labelmap[[w]] <- {
                    #covariates for parametric terms become by-variables:   
                    where <- sapply(m.smooth, function(x) x$by) == names(labelmap)[w]
                    sapply(m.smooth[where], function(x) x$label)
                }
        }
        if(length(where.ffpc)){
            ind <- 1
            for(w in where.ffpc){
                labelmap[[w]] <- {
                    #PCs for X become by-variables:   
                    where <- sapply(m.smooth, function(x) x$id) == ffpcterms[[ind]]$id  
                    sapply(m.smooth[where], function(x) x$label)
                }
                ind <- ind+1
            } 
        }
        labelmap[-c(where.par, where.ffpc)] <- lbls[pmatch(
                        sapply(labelmap[-c(where.par, where.ffpc)], function(x){
                                    tmp <- eval(parse(text=x))
                                    if(is.list(tmp)){
                                        return(tmp$label)  
                                    } else {
                                        return(x)   
                                    } 
                                }), lbls)]
    } else{
        labelmap[1:length(labelmap)] <-  lbls[pmatch(
                        sapply(labelmap, function(x){
                                    tmp <- eval(parse(text=x))
                                    if(is.list(tmp)){
                                        return(tmp$label)  
                                    } else {
                                        return(x)   
                                    }
                                }), lbls)]
    } 
    # check whether any parametric terms were left out & add them
    if(any(nalbls <- sapply(labelmap, function(x) any(is.na(x))))){
        labelmap[nalbls] <- trmmap[nalbls]
    }
    
    names(m.smooth) <- lbls
    if(as.character(algorithm) %in% c("gamm4","gamm")){
        m$gam$smooth <- m.smooth 
    } else{
        m$smooth  <- m.smooth
    } 
    
    ret <-  list(
            call=match.call(),
            formula=formula, 
            termmap=trmmap, 
            labelmap=labelmap, 
            responsename = responsename,
            nobs=nobs,
            nyindex=nyindex,
            yindname = yindname,
            yind=yind,
            where=list(
                    where.c=where.c,
                    where.ff=where.ff,
                    where.ffpc=where.ffpc,
                    where.sff=where.sff,
                    where.s=where.s,
                    where.te= where.te,
                    where.t2=where.t2,
                    where.par=where.par
            ),
            ff=ffterms,
            ffpc=ffpcterms)
    
    if(as.character(algorithm) %in% c("gamm4","gamm")){
        m$gam$pffr <- ret
        class(m$gam) <- c("pffr", class(m$gam))
    } else {
        m$pffr <- ret
        class(m) <- c("pffr", class(m))
    }
    return(m)
}# end pffr()

#' Simple bootstrap CIs for pffr
#' 
#' This function resamples observations in the data set to obtain approximate CIs for different
#' terms and coefficient functions that correct for the effects of dependency and heteroskedasticity
#' of the residuals along the index of the functional response, i.e., it aims for correct inference
#' if the residuals along the index of the functional response are not i.i.d.
#' 
#' @param object a fitted \code{\link{pffr}}-model 
#' @param n1 see \code{\link{coef.pffr}}
#' @param n2 see \code{\link{coef.pffr}}
#' @param n3 see \code{\link{coef.pffr}}
#' @param B  number of bootstrap replicates, defaults to (a measly) 100
#' @param parallel see \code{\link[boot]{boot}}
#' @param cl see \code{\link[boot]{boot}}
#' @param ncpus see \code{\link[boot]{boot}}. Defaults to \code{getOption("boot.ncpus", 1L)} (like \code{boot}).
#' @param conf desired levels of bootstrap CIs, defaults to 0.90 and 0.95
#' @param type type of bootstrap interval, see \code{\link[boot]{boot.ci}}. Defaults to "percent" for percentile-based CIs. 
#' @param showProgress TRUE/FALSE 
#' @param ... not used
#' @return a list with similar structure as the return value of \code{\link{coef.pffr}}, containing the 
#'         original point estimates of the various terms along with their bootstrap CIs.
#' @author Fabian Scheipl
coefboot.pffr <- function(object,  
        n1=100, n2=40, n3=20, 
        B = 100, ncpus = getOption("boot.ncpus", 1),
        parallel = c("no", "multicore", "snow"), cl=NULL,
        conf=c(.9,.95), type="percent", showProgress=TRUE, ...){
    
    if(is.null(ncpus)) ncpus <- getOption("boot.ncpus", 1L) 
    ## 
    modcall <- object$pffr$call
    modcall$formula <- object$pffr$formula
    modcall$data <- environment(object$pffr$formula)$data
    #modcall$fit <- FALSE
    if(is.symbol(modcall$yind)){
        modcall$yind <- eval(modcall$yind, environment(object$pffr$formula))   
    }
    if(is.symbol(modcall$hatSigma)){
        modcall$hatSigma <- eval(modcall$hatSigma, environment(object$pffr$formula))   
    }

    if(!is.null(modcall$algorithm)){
        if(modcall$algorithm != "gam"){
            stop("bootstrap implemented only for gam-fits.\n")
        }
    }
    #modcall$start <- object$coefficients
    #modcall$G <- eval(modcall)
    modcall$fit <- TRUE
    yind <- modcall$yind 
    
    ## refit models on bootstrap data sets & save fits
    cat("starting bootstrap ")
    # check whether this is already a child process, and
    # switch of parallelization if yes
    if(parallel:::isChild()) parallel <- "no"
    bootfct <- function(modcall, data, indices) {
        modcall$weights <- rep(table(factor(indices, levels=1:nrow(data))), each = length(yind))
        #modcall$G$w <- rep(table(factor(indices, levels=1:nrow(data))), each = length(yind))        
        mb <- eval(as.call(modcall))
        coefs <- coef(mb, se = FALSE, n1=n1, n2=n2, n3=n3)
        ##save results as one long vecotr as boot can't deal with lists
        coefvec <- c(coefs$pterms[,"value"], 
                sapply(coefs$smterms, function(x) x$value))
        if(showProgress) cat(".")
        return(unlist(coefvec))
    }
    bootreps <- boot(data=modcall$data, statistic=bootfct, R=B, 
            modcall=modcall, parallel=parallel, ncpus=ncpus)
    if(showProgress) cat("done.\n")
    
    ## disentangle bootreps
    coefboot <- coef(object, se = FALSE, n1=n1, n2=n2, n3=n3)
    
    ptrms <- rownames(coefboot$pterms)
    ptrmsboot <- matrix(bootreps$t[,1:length(ptrms)], )
    ptrmsboot <- t(ptrmsboot)
    
    smtrms <- names(coefboot$smterms)
    smtrmvallengths <- sapply(coefboot$smterms, function(x) length(x$value))
    
    #drop NULL-entries
    NULLentries <- smtrmvallengths==0
    if(any(NULLentries)){
        smtrms <- smtrms[!NULLentries]
        smtrmvallengths <- smtrmvallengths[!NULLentries]
        coefboot$smterms <- coefboot$smterms[!NULLentries]
    }
        
    start <- c(1, cumsum(head(smtrmvallengths,-1))+1) 
    stop <- cumsum(smtrmvallengths)
    smtrmsboot <- vector(length(smtrms), mode="list")
    for(i in 1:length(start)){
        smtrmsboot[[i]] <- t((bootreps$t[,-(1:length(ptrms))])[,start[i]:stop[i]])
    }
    names(smtrmsboot) <- smtrms
    
    
    ## get bootstrap CIs
    cat("calculating bootstrap CIs....")
    mylapply <- if(parallel=="multicore"){
                parallel::mclapply  
            } else {
                if(parallel=="snow"){
                    if(is.null(cl)){
                        cl <- parallel::makePSOCKcluster(rep("localhost", 
                                        ncpus))  
                    } 
                    if (RNGkind()[1L] == "L'Ecuyer-CMRG") {
                        parallel::clusterSetRNGStream(cl)
                    }
                    function(X, FUN) {
                        res <- parallel::parLapply(cl, X, FUN)
                        parallel::stopCluster(cl)
                        return(res)
                    }    
                } else lapply
            }    
    
    getCIs <- function(which=1:length(bootreps$t0), conf=.95, type="bca"){
        ret <- mylapply(which, function(i){
                    ci <- boot.ci(bootreps, conf=conf, index=i, 
                            #fix stupidity in boot.ci: option name has to be "perc", 
                            # return value is named "percent"
                            type=ifelse(type=="percent", "perc", type))
                    ret <- ci[[type]][,4:5]
                    return(as.vector(ret))          
                })
        ret <- do.call(cbind, ret)
        rownames(ret) <- paste(c(((1-conf)/2), 1-(1-conf)/2)*100,"%", sep="")
        return(t(ret))
    }
    coefboot$pterms <- cbind(coefboot$pterms, getCIs(which=1:length(ptrms), conf=conf, type=type))
    for(i in 1:length(start)){
        coefboot$smterms[[i]] <- cbind(coefboot$smterms[[i]]$coef,
                getCIs(which=length(ptrms)+(start[i]:stop[i]), conf=conf, type=type))
    }
    names(coefboot$smterms) <- smtrms
    
    return(structure(coefboot, call=match.call()))
}
