# wrapper functions: do penalized function-on-function via mgcv::gam(m)
# 
# Author: fabians
# 13.06.2011, 10:54:19: initial commit
#
# TODO: add functionality for concurrent model
# TODO: supplying k as a variable name breaks down because it is propagated without evaluation
#       and so k is a list(<expression>, default.k)  where a vector is expected.
###############################################################################
#' Penalized function-on-function regression
#' 
#' Implements additive regression for functional and scalar covariates and functional responses.
#' This function is a wrapper for \code{mgcv}'s \code{\link[mgcv]{gam}} and its siblings to fit models of the general form \cr
#' \eqn{E(Y_i(t)) = g(\mu(t) + \int X_i(s)\beta(s,t)ds + f(z_{1i}, t) + f(z_{2i}) + z_{3i} \beta_3(t) + \dots )}\cr
#' with a functional (but not necessarily continuous) response \eqn{Y(t)}, response function \eqn{g},
#' (optional) smooth intercept \eqn{\mu(t)}, (multiple) functional covariates \eqn{X(t)} and scalar covariates
#' \eqn{z_1}, \eqn{z_2}, etc. 
#' 
#' @section Details: The routine can estimate
#' \enumerate{
#'     \item (nonlinear, and possibly multivariate) effects of (one or multiple) scalar covariates \eqn{z} that vary smoothly over the index
#' \eqn{t} of \eqn{Y(t)} (e.g. \eqn{f(z_{1i}, t)}, specified in the \code{formula} simply as \code{~s(z1)}),
#' \item (nonlinear) effects of scalar covariates that are constant over \eqn{t} (e.g. \eqn{f(z_{2i})}, specified as \code{~c(s(z2))}, or
#' \eqn{\beta_2 z_{2i}}, specified as \code{~c(z2)}),
#'  \item linear functional effects of scalar (numeric or factor) covariates that vary smoothly over \eqn{t} (e.g. \eqn{z_{3i} \beta_3(t)}, specified as \code{~z3}),
#'  \item function-on-function regression terms (e.g. \eqn{\int X_i(s)\beta(s,t)ds}, specified as \code{~ff(X, yindex=t, xindex=s)}, see \code{\link{ff}}). 
#' }
#' Use the \code{c()}-notation to denote model terms that are constant over the index of the functional response.\cr
#' 
#' Internally, univariate smooth terms without a \code{c()}-wrapper are expanded into bivariate smooth terms in the original covariate and the index
#' of the functional response. Bivariate smooth terms (\code{s(), te()} or \code{t2()}) without a \code{c()}-wrapper are expanded into trivariate smooth
#' terms in the original covariates and the index of the functional response. Linear terms for scalar covariates or categorical covariates are expanded
#' into varying coefficient terms, varying smoothly over the index of the functional response. For factor variables, a separate smooth function 
#' with its own smoothing parameter is estimated for each level of the factor.\cr
#' Functional random intercepts \eqn{B_{0g(i)}(t)} for a grouping variable \code{g} can be specified via \code{~s(g, bs="re")}),
#' functional random slopes \eqn{u_i B_{1g(i)}(t)} in a numeric variable \code{u} via \code{~s(g, u, bs="re")}).\cr 
#' The marginal spline basis used for the index of the 
#' the functional response is specified via the \emph{global} argument \code{bs.yindex}. If necessary, this can be overriden for any specific term by 
#' supplying a \code{bs.yindex}-argument to that term in the formula, e.g. \code{~s(x, bs.yindex=list(bs="tp", k=7))} would yield a tensor
#' product spline for which
#' the marginal basis for the index of the response are 7 cubic thin-plate spline functions overriding the global default for the basis and penalty on
#' the index of the response given by the \emph{global} \code{bs.yindex}-argument .\cr 
#' Use \code{~-1 + c(1) + ...} to specify a model with only a constant and no functional intercept. \cr
#' 
#' The functional covariates have to be supplied as a \eqn{n} by <no. of evaluations> matrices, 
#' i.e. each row is one functional observation. 
#' For data on a regular grid, the functional response
#' is supplied in the same format, i.e. as a matrix-valued entry in \code{data},  which can contain missing values.
#' If the functional responses are sparse or irregular (i.e., not evaluated on the same evaluation points
#' across all observations), the \code{ydata}-argument can be used to specify the responses: 
#' \code{ydata} must be a \code{data.frame} with 3 columns called \code{'.obs', '.index', '.value'}
#'  which specify which curve the point belongs to (\code{'.obs'}=\eqn{i}), at which \eqn{t} it was
#'  observed (\code{'.index'}=\eqn{t}), and the observed value (\code{'.value'}=\eqn{Y_i(t)}). Note that  
#'  the vector of unique sorted entries in \code{ydata$.obs} must be equal to 
#'  \eqn{(1, \ldots, <no. of observations>)} to ensure the correct association of entries in \code{ydata} to the 
#'  corresponding rows of \code{data}.
#' For both regular and irregular functional responses, the model is then fitted with the data in long 
#' format, i.e., for data on a grid the rows of the
#' matrix of the functional response evaluations \eqn{Y_i(t)} are stacked
#' into one long vector and the covariates are expanded/repeated correspondingly. This means 
#' the models get quite big fairly fast,
#' since the effective number of rows in the design matrix is number of observations times 
#' number of evaluations of \eqn{Y(t)} per observation.\cr 
#' 
#' Note that \code{pffr} overrides the default identifiability constraints 
#' (\eqn{\sum_{i,t} \hat f(z_i, x_i, t) = 0}) implemented in \code{mgcv} 
#' for tensor product terms whose marginals include the index \eqn{t} of the functional 
#' response.  Instead, \eqn{\sum_i \hat f(z_i, x_i, t) = 0} for all \eqn{t}
#' is enforced, so that effects varying over \eqn{t} can be interpreted 
#' as local deviations from the global functional intercept. 
#' We recommend using centered scalar covariates for terms like
#' \eqn{z \beta(t)} (\code{~z}) and centered functional covariates with \eqn{\sum_i X_i(t) = 0}
#' for all \eqn{t} in \code{ff}-terms
#' so that the global functional intercept can be interpreted as the global mean function.      
#' For irregular or sparse \eqn{Y_i(t)} supplied via \code{ydata} it is not trivial to specify these 
#' custom constraints. They are enforced on the synthetic grid of \eqn{t}-values returned in 
#' the \code{pffr$yind} entry of the return object. 
#'   
#' @param formula a formula with special terms as for \code{\link[mgcv]{gam}}, with additional special terms \code{\link{ff}()} and \code{c()}
#' @param yind a vector with length equal to the number of columns of the matrix of functional responses giving the vector of evaluation points \eqn{(t_1, \dots ,t_{G})}.
#'   If not supplied, \code{yind} is \code{1:ncol(<response>)}.   
#' @param algorithm the name of the function used to estimate the model. Defaults to \code{\link[mgcv]{gam}} if
#' the matrix of functional responses has less than \code{2e5} data points
#' 	 and to \code{\link[mgcv]{bam}} if not. "gamm" (see \code{\link[mgcv]{gamm}}) and "gamm4" (see \code{\link[gamm4]{gamm4}}) are valid options as well.  
#' @param data an (optional) \code{data.frame} or a named list containing the data. 
#' The functional response and functional covariates have to be supplied as n by <no. of evaluations> matrices, i.e. each row is one functional observation.
#' The model is then fitted with the data in long format, i.e., the rows of the matrix of the functional response evaluations \eqn{Y_i(t)} are stacked
#' into one long vector and the covariates are expanded/repeated correspondingly. This means the models get quite big fairly fast,
#' since the effective number of rows in the design matrix is number of observations times number of evaluations of \eqn{Y(t)} per observation.\cr 
#' @param ydata an (optional) \code{data.frame} supplying functional responses that are not observed on a regular grid. See Details. 
#' @param method Defaults to \code{"REML"}-estimation, including of unknown scale. See \code{\link[mgcv]{gam}} for details.  
#' @param bs.yindex a named (!) list giving the parameters for spline bases on the index of the functional response. 
#'  Defaults to \code{list(bs="ps", k=5, m=c(2, 1))}, i.e. 5 cubic B-splines bases with first order difference penalty. 
#' @param bs.int a named (!) list giving the parameters for the spline basis for the global functional intercept.  
#'  Defaults to \code{list(bs="ps", k=20, m=c(2, 1))}, i.e. 20 cubic B-splines bases with first order difference penalty.  
#' @param tensortype which typ of tenor product splines to use. One of "\code{\link[mgcv]{te}}" or "\code{\link[mgcv]{t2}}", defaults to \code{te}
#' @param ... additional arguments that are valid for \code{\link[mgcv]{gam}} or \code{\link[mgcv]{bam}}. \code{weights, subset, offset} are not yet implemented!
#' @return a fitted \code{pffr}-object, which is a \code{\link[mgcv]{gam}}-object with some additional information in an \code{pffr}-entry. If \code{algorithm} is
#'   \code{"gamm"} or \code{"gamm4"}, only the \code{$gam} part of the returned list is modified in this way.
#' @author Fabian Scheipl, Sonja Greven
#' @seealso \code{\link[mgcv]{smooth.terms}} for details of \code{mgcv} syntax and available spline bases and penalties.
#' @references Ivanescu, A., Staicu, A.-M., Scheipl, F. and Greven, S. (2012). 
#'  Penalized function-on-function regression. (under revision) 
#' \url{http://biostats.bepress.com/jhubiostat/paper240/}
#' 
#'  Scheipl, F., Staicu, A.-M. and Greven, S. (2012). Functional Additive Mixed Models.
#'  (submitted) \url{http://arxiv.org/abs/1207.5947}
#' @examples  
#' ###############################################################################
#' # univariate model: 
#' # Y(t) = f(t)  + \int X1(s)\beta(s,t)ds + eps  
#' data1 <- pffrSim(scenario="ff")
#' t <- attr(data1, "yindex")
#' s <- attr(data1, "xindex")
#' m1 <- pffr(Y ~ ff(X1, xind=s), yind=t, data=data1)
#' summary(m1)
#' plot(m1, pers=TRUE)
#' 
#' \dontrun{
#' ###############################################################################
#' # multivariate model: 
#' # Y(t) = f0(t)  + \int X1(s)\beta1(s,t)ds + \int X2(s)\beta2(s,t)ds +
#' #		xlin \beta3(t) + f1(xte1, xte2) + f2(xsmoo, t) + beta4 xconst + eps  
#' data2 <- pffrSim(scenario="all", n=200)
#' t <- attr(data2, "yindex")
#' s <- attr(data2, "xindex")
#' m2 <- pffr(Y ~  ff(X1, xind=s) + #linear function-on-function		
#'                 ff(X2, xind=s) + #linear function-on-function
#'                 xlin  +  #varying coefficient term
#'                 c(te(xte1, xte2)) + #bivariate smooth term in xte1 & xte2, const. over Y-index
#'                 s(xsmoo) + #smooth effect of xsmoo varying over Y-index
#'                 c(xconst), # linear effect of xconst constant over Y-index
#'         yind=t,
#'         data=data2)
#' summary(m2)
#' plot(m2, pers=TRUE)
#' 
#' #############################################################################
#' # sparse data (keep only 20% of function evaluations): 
#' set.seed(121456)
#' ydata1 <- data.frame(cbind(.obs=as.vector(row(data1$Y)),
#'                            .index=rep(attr(data1, "yindex"), e=nrow(data1)),
#'                            .value=as.vector(data1$Y)))
#' ydata1 <- ydata1[sample(1:nrow(ydata1), 0.2*nrow(ydata1)), ]
#' 
#' s <- attr(data1, "xindex")
#' t <- attr(data1, "yindex")
#' m1.sparse <- pffr(Y ~ ff(X1, xind=s), data=data1, ydata=ydata1, yind=t)
#' summary(m1.sparse)
#' plot(m1.sparse, pers=TRUE)
#' } 
#FIXME: te-constraints for bam broken cause bam.fit now samples data to set up bases
pffr <- function(
  formula,
  yind,
  data=NULL,
  ydata=NULL,
  algorithm = NA, 
  method="REML",
  tensortype = c("te", "t2"),
  bs.yindex = list(bs="ps", k=5, m=c(2, 1)), # only bs, k, m are propagated...
  bs.int = list(bs="ps", k=20, m=c(2, 1)), # only bs, k, m are propagated...
  ...
){
  # TODO: weights, subset, offset args!
  
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
  
  sparseOrNongrid <- !is.null(ydata)
  if(sparseOrNongrid){
    stopifnot(ncol(ydata)==3)
    stopifnot(c(".obs", ".index", ".value") == colnames(ydata))
  }
  
  tf <- terms.formula(formula, specials=c("s", "te", "t2", "ff", "c", "sff", "ffpc", "pcre"))
  trmstrings <- attr(tf, "term.labels")
  terms <- sapply(trmstrings, function(trm) as.call(parse(text=trm))[[1]], simplify=FALSE) 
  #ugly, but getTerms(formula)[-1] does not work for terms like I(x1:x2) 
  frmlenv <- environment(formula)
  
  
  where.c <- attr(tf, "specials")$c - 1    # indices of scalar (not index-varying) terms
  where.ff <- attr(tf, "specials")$ff - 1  # function-on-function terms
  where.sff <- attr(tf, "specials")$sff - 1  #smooth function-on-function terms
  where.s <- attr(tf, "specials")$s - 1    # smooth terms 
  where.te <- attr(tf, "specials")$te - 1  # tensor product terms
  where.t2 <- attr(tf, "specials")$t2 - 1  # type 2 tensor product terms
  where.ffpc <- attr(tf, "specials")$ffpc - 1  # PC-based function-on-function terms 
  where.pcre <- attr(tf, "specials")$pcre - 1  # functional random effects/residuals with PC-basis 
  if(length(trmstrings)) {
    where.par <- which(!(1:length(trmstrings) %in%
                           c(where.c, where.ff, where.sff, where.ffpc, where.pcre, where.s, where.te, where.t2))) 
    # indices of linear/factor terms with varying coefficients over yind
  } else where.par <- numeric(0)
  
  responsename <- attr(tf,"variables")[2][[1]]
  
  #start new formula
  newfrml <- paste(responsename, "~", sep="")
  newfrmlenv <- new.env()
  evalenv <- if("data" %in% names(call)) eval.parent(call$data) else NULL
  
  if(sparseOrNongrid){
    nobs <- length(unique(ydata$.obs))
    stopifnot(all(sort(unique(ydata$.obs)) == 1:nobs))
    
    #works for data-lists or matrix-valued covariates as well:
    nobs.data <- nrow(as.matrix(data[[1]])) 
    stopifnot(nobs == nobs.data)
    ntotal <- nrow(ydata)
    
    #generate yind for estimates/predictions etc
    yind <- if(length(unique(ydata$.index))>100){
      seq(min(ydata$.index), max(ydata$.index), l=100) 
    } else {
      sort(unique(ydata$.index))
    }
    nyindex <- length(yind)
  } else {
    nobs <- nrow(eval(responsename,  envir=evalenv, enclos=frmlenv))
    nyindex <- ncol(eval(responsename,  envir=evalenv, enclos=frmlenv))    
    ntotal <- nobs*nyindex
  }
  
  if(missing(algorithm)||is.na(algorithm)){
    algorithm <- ifelse(ntotal > 1e5, "bam", "gam")
  } 
  algorithm <- as.symbol(algorithm)
  if(as.character(algorithm)=="bam" && !("chunk.size" %in% names(call))){
    call$chunk.size <- max(ntotal/5, 10000) 
    #same default as in bam
  }
  ## no te-terms possible in gamm4:
  if(as.character(algorithm)=="gamm4") stopifnot(length(where.te)<1)
  
  
  if(!sparseOrNongrid){
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
      yindname <- deparse(substitute(yind))
      if(!is.null(data) && !is.null(data[[deparse(substitute(yind))]])){
        yind <- data[[deparse(substitute(yind))]]
      } 
      stopifnot(is.vector(yind), is.numeric(yind), 
                length(yind) == nyindex)
    }
    #make sure it's a valid name 
    if(length(yindname)>1) yindname <- "yindex"	
    # make sure yind is sorted
    stopifnot(all.equal(order(yind), 1:nyindex))
    
    
    yindvec <- rep(yind, times = nobs)
    yindvecname <- as.symbol(paste(yindname,".vec",sep=""))
    assign(x=deparse(yindvecname), value=yindvec, envir=newfrmlenv)
    
    #assign response in _long_ format to newfrmlenv
    assign(x=deparse(responsename), value=as.vector(t(eval(responsename, envir=evalenv, 
                                                           enclos=frmlenv))), 
           envir=newfrmlenv)
    
    
    if(any(is.na(get(as.character(responsename), newfrmlenv)))){
      missingindname <- paste0(".PFFRmissingResponses.", 
                               paste(sample(letters, 6), collapse=""))
      missingind <- which(is.na(get(as.character(responsename), newfrmlenv)))
    } else missingind <- NULL 
    
    # repeat which row in <data> how many times 
    stackpattern <- rep(1:nobs, each=nyindex)
    
  } else {
    
    yindname <- "yindex"
    yindvec <- ydata$.index
    yindvecname <- as.symbol(paste(yindname,".vec",sep=""))
    assign(x=deparse(yindvecname), value=ydata$.index, envir=newfrmlenv)
    
    #assign response in _long_ format to newfrmlenv
    assign(x=deparse(responsename), value=ydata$.value, envir=newfrmlenv)
    
    missingind <- NULL
    
    # repeat which row in <data> how many times: 
    stackpattern <- ydata$.obs
  }
  
  
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
    
    #apply limits function and assign stacked data to newfrmlenv
    lapply(ffterms, function(x){
      assign(x=x$yindname, 
             value=matrix(yindvec, nrow=length(yindvec), ncol=length(x$xind)), 
             envir=newfrmlenv)
      assign(x=x$xindname, 
             value=matrix(x$xind, nrow=length(yindvec), ncol=length(x$xind), byrow=TRUE), 
             envir=newfrmlenv) 
      if(!is.null(x$LX)){
        # for ff: stack weights * covariate
        LStacked <- x$LX[stackpattern,]
      } else {
        # for sff: stack weights only
        LStacked <- x$L[stackpattern,]
      } 
      if(!is.null(x$limits)){
        ind0 <- !t(outer(get(x$xindname, newfrmlenv)[1,], 
                         get(x$yindname, newfrmlenv)[,1], x$limits))
        LStacked[ind0] <- 0
      }    
      assign(x=x$LXname, 
             value=LStacked, 
             envir=newfrmlenv)
      #assign additional X-term for sff:
      if(!is.null(x$X)){
        assign(x=x$xname, 
               value=x$X[stackpattern, ], 
               envir=newfrmlenv)
      }
      invisible(NULL)
    })  
  } else ffterms <- NULL
  
  if(length(where.ffpc)){ ##TODO for sparse
    ffpcterms <- lapply(terms[where.ffpc], function(x){
      eval(x, envir=evalenv, enclos=frmlenv)
    })
    
    lapply(ffpcterms, function(trm){
      lapply(colnames(trm$data), function(nm){
        assign(x=nm, value=trm$data[stackpattern, nm], envir=newfrmlenv)
        invisible(NULL)
      })
      invisible(NULL)
    })  
    
    
    getFfpcFormula <- function(trm) {
      frmls <- lapply(colnames(trm$data), function(pc) {
        arglist <- c(name="s", x = as.symbol(yindvecname), by= as.symbol(pc), 
                     id=trm$id, trm$splinepars)
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
      if(!sparseOrNongrid && all(trm$yind==yind)){
        lapply(colnames(trm$efunctions), function(nm){
          assign(x=nm, value=trm$efunctions[rep(1:nyindex, times=nobs), nm], 
                 envir=newfrmlenv)
          invisible(NULL)
        })    
      } else {
        # don't ever extrapolate eigenfunctions: 
        stopifnot(min(trm$yind)<=min(yind))
        stopifnot(max(trm$yind)>=max(yind))
        
        # interpolate given eigenfunctions to observed index values:
        lapply(colnames(trm$efunctions), function(nm){
          tmp <- approx(x=trm$yind, 
                        y=trm$efunctions[, nm], 
                        xout=yindvec,
                        method = "linear")$y
          assign(x=nm, value=tmp, 
                 envir=newfrmlenv)
          invisible(NULL)
        }) 
      }
      assign(x=trm$idname, value=trm$id[stackpattern], envir=newfrmlenv)
      invisible(NULL)
    })  
    
    newtrmstrings[where.pcre] <- sapply(pcreterms, function(x) {
      safeDeparse(x$call)
    })
  }else pcreterms <- NULL
  
  
  #transform: s(x, ...), te(x, z,...), t2(x, z, ...) --> <te|t2>(x, <z,> yindex, ..., <bs.yindex>)
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
            if(exists("missingindname")){
              lst[[i]][[llst+3]] <- missingindname
              names(lst[[i]])[length(names(lst[[i]]))] <- "missingindname"
            }
          } else {
            if(exists("missingindname")){
              lst[[i]] <-   bquote(list(impose.ffregC = .(TRUE), 
                                        nobs=.(nobs),
                                        missingindname=.(missingindname)))
            } else {
              lst[[i]] <- bquote(list(impose.ffregC = .(TRUE), 
                                      nobs=.(nobs)))
            }
          }
        }    
        return(lst)
      }
      # xt has to be supplied as a list, with length(x$d) entries,
      # each of which is a list or NULL:
      stopifnot(x$xt[[1]]==as.symbol("list") && 
                  length(x$xt)==length(xnew$d) &&  
                  # =length(x$d)+1, since first element in parse tree is ``list''
                  all(sapply(2:length(x$xt), function(i) 
                    x$xt[[i]][[1]] == as.symbol("list") ||
                               is.null(eval(x$xt[[i]][[1]])))))
      xtra <- add.impose(x$xt)
      
      append <- if(exists("missingindname")) {
        bquote(list(impose.ffregC = .(TRUE), 
                    nobs=.(nobs),
                    missingindname=.(missingindname)))
      } else {
        bquote(list(impose.ffregC = .(TRUE), 
                    nobs=.(nobs)))
      }
      xtra[[length(xnew$d)+1]] <- append
      xtra
    } else {
      imposearg <- if(exists("missingindname")) {
        bquote(list(impose.ffregC = .(TRUE), 
                    nobs=.(nobs),
                    missingindname=.(missingindname)))
      } else {
        bquote(list(impose.ffregC = .(TRUE), 
                    nobs=.(nobs)))
      }
      
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
    # evalenv below used to be list2env(eval.parent(call$data)), frmlenv),
    # but that assigned everything in <data> to the global workspace if frmlenv was the global
    # workspace.  
    evalenv <- if("data" %in% names(call))  {
      list2env(eval.parent(call$data)) 
      } else frmlenv 
    lapply(terms[where.notff], 
           function(x){
             #nms <- all.vars(x)
             isC <- deparse(x) %in% sapply(terms[where.c], deparse) 
             if(isC){
               # drop c()
               # FIXME: FUGLY! 
               x <- formula(paste("~", gsub("\\)$", "",
                                            gsub("^c\\(", "", deparse(x)))))[[2]]
             } 
             ## remove names in xt, k, bs,  information (such as variable names for MRF penalties etc)
             nms <- if(!is.null(names(x))){
               all.vars(x[names(x) %in% c("", "by")]) 
             }  else all.vars(x)
             
             
             sapply(nms, function(nm){
               stopifnot(length(get(nm, envir=evalenv)) == nobs)
               assign(x=nm, 
                      value=get(nm, envir=evalenv)[stackpattern],
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
  # (not possible for sparse/nongrid as constraint definition does
  # not carry over to that case...)
  if(!(as.character(algorithm) %in% c("gamm4"))){
    suppressMessages(
      trace(mgcv::smooth.construct.tensor.smooth.spec, 
            at = max(which(sapply(as.list(body(mgcv::smooth.construct.tensor.smooth.spec)), 
                                  function(x) any(grepl(x, pattern="object$C", fixed=TRUE))))) + 1, 
            print=FALSE,
            tracer = quote({
              M <- length(object$margin) 
              if(!is.null(object$margin[[M]]$xt$impose.ffregC)){
                ## constraint: sum_i f(z_i, t) = 0 \forall t
                #cat("imposing constraints..\n")
                #some really ugly shit right here: look for "sparseOrNongrid" etc 
                #  anywhere in the call stack
                pffrframe <- which(sapply(lapply(sys.frames(), ls), 
                                    function(x) any(grepl(pattern="sparseOrNongrid", x=x))))[1]
                irreg <- get("sparseOrNongrid", sys.frame(pffrframe))
                
                if(irreg){
                  unique.obs <- !duplicated(get("stackpattern", envir=sys.frame(pffrframe)))
                  unique.yind <- !duplicated(get("ydata", envir=sys.frame(pffrframe))$.index)
                  yind <- get("ydata", envir=sys.frame(pffrframe))$.index[unique.yind]
                  yind.grid <- get("yind", envir=sys.frame(pffrframe))
                  
                  ## eval tensor basis on gridded data, use this to define constraint:
                  
                  # 1) extract design matrix for covariates for each observed curve 
                  # and stack'em observationwise:
                  Xlist <- list()
                  for(m in 1:(M-1)){
                    Xlist[[m]] <- Xm[[m]][unique.obs,]
                  }
                  Xx <- if(length(Xlist)>1) {
                    tensor.prod.model.matrix(Xlist)
                  } else Xlist[[1]]
                  Xx <- Xx[rep(1:object$margin[[M]]$xt$nobs, e=length(yind.grid)),]
                  
                  # 2) evaluate yindex-basis on a regular grid by interpolating the basis 
                  #  given for the irregular yindex-values and stack'em observationwise:
                  Xt <- apply(Xm[[M]][unique.yind,], 2, 
                              function(basis){
                                approx(x=yind, y=basis, xout=yind.grid)$y
                              })[rep(1:length(yind.grid), t=object$margin[[M]]$xt$nobs),]
                  # 
                  Xgrid <- tensor.prod.model.matrix(list(Xx,Xt))
                }
                
                nygrid <- ifelse(irreg, length(yind.grid), 
                                 length(unique(data[[object$margin[[M]]$term]])))
                
                ## Ctmp =  ((1,...,1) \otimes I_G) sums over all observations 
                ##  for each value of yindex
                ## C = ((1,...,1) \otimes I_G) * B
                Ctmp <- kronecker(t(rep(1, object$margin[[M]]$xt$nobs)), diag(nygrid))
                if(!irreg & !is.null(object$Cp) & (ncol(Ctmp) > nrow(object$X))){
                  #drop rows for missing obs.
                  Ctmp <- Ctmp[, -get("missingind", envir=sys.frame(pffrframe))]
                }
                Xmat <- if(irreg){
                  Xgrid
                } else {
                  object$X
                }                               
                C <- Ctmp %*% Xmat
                
                ## we need the number of effective constraints <nC> to correspond to the
                ## rank of the constraint matrix C, which is the min. of number of basis
                ## functions for t (object$margin[[length(object$margin)]]$bs.dim) and
                ## timepoints (<nygrid>)
                nC <- min(nygrid, object$margin[[M]]$bs.dim)
                object$C <- object$Cp <- C[seq(1, nygrid, length.out = nC), ]
              }}))
    )
    
    suppressMessages(
      trace(mgcv::smooth.construct.t2.smooth.spec, 
            at = max(which(sapply(as.list(body(mgcv::smooth.construct.t2.smooth.spec)), function(x) any(grepl(x, pattern="object$Cp", fixed=TRUE))))) + 1, 
            print=FALSE,
            tracer = quote({
              if(!is.null(object$margin[[length(object$margin)]]$xt$impose.ffregC) &&
                   object$margin[[length(object$margin)]]$xt$impose.ffregC){
                ## constraint: sum_i f(z_i, t) = 0 \forall t
                #cat("imposing constraints..\n")
                nygrid <- length(unique(data[[object$margin[[length(object$margin)]]$term]]))
                
                ## C = ((1,...,1) \otimes I_G) * B
                Ctmp <- kronecker(t(rep(1, object$margin[[length(object$margin)]]$xt$nobs)), diag(nygrid))
                if(ncol(Ctmp) > nrow(object$X)){
                  #drop rows for missing obs.
                  #Ctmp <- Ctmp[, - get(object$margin[[length(object$margin)]]$xt$missingindname, .GlobalEnv)]
                  Ctmp <- Ctmp[, -get("missingind", envir=sys.frame(3))]
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
      suppressMessages(try(untrace(mgcv::smooth.construct.tensor.smooth.spec), silent = TRUE))
      suppressMessages(try(untrace(mgcv::smooth.construct.t2.smooth.spec), silent = TRUE))
    })
  }
  
  
  
  # call algorithm to estimate model
  m <- eval(newcall)
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
      for(w in where.par){
        # only combine if <by>-variable is a factor!
        if(is.factor(get(names(labelmap)[w], envir=newfrmlenv))){
          labelmap[[w]] <- {
            #covariates for parametric terms become by-variables:   
            where <- sapply(m.smooth, function(x) x$by) == names(labelmap)[w]   
            sapply(m.smooth[where], function(x) x$label)
          }    
        } else {
          labelmap[[w]] <- paste0("s(",yindvecname,"):",names(labelmap)[w])
        }
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
        ## FUGLY: check whether x is a function call of some sort
        ##  or simply a variable name.
        if(length(parse(text=x)[[1]]) != 1){
          tmp <- eval(parse(text=x))
          return(tmp$label)  
        } else {
          return(x)   
        } 
      }), lbls)]
  } else{
    labelmap[1:length(labelmap)] <-  lbls[pmatch(
      sapply(labelmap[1:length(labelmap)], function(x){
        ## FUGLY: check whether x is a function call of some sort
        ##  or simply a variable name.
        if(length(parse(text=x)[[1]]) != 1){
          tmp <- eval(parse(text=x))
          return(tmp$label)  
        } else {
          return(x)   
        } 
      }), lbls)]
  } 

  
  # check whether any parametric terms were left out & add them
  if(any(nalbls <- sapply(labelmap, function(x) any(is.null(x) | (!is.null(x)&&is.na(x)))))){
    labelmap[nalbls] <- trmmap[nalbls]
  }
  
  names(m.smooth) <- lbls
  if(as.character(algorithm) %in% c("gamm4","gamm")){
    m$gam$smooth <- m.smooth 
  } else{
    m$smooth  <- m.smooth
  } 
  
  ret <-  list(
    call=call,
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
      where.pcre=where.pcre,
      where.par=where.par
    ),
    ff=ffterms,
    ffpc=ffpcterms,
    pcreterms=pcreterms,
    missingind = missingind,
    sparseOrNongrid=sparseOrNongrid, 
    ydata=ydata)
  
  if(as.character(algorithm) %in% c("gamm4","gamm")){
    m$gam$pffr <- ret
    class(m$gam) <- c("pffr", class(m$gam))
  } else {
    m$pffr <- ret
    class(m) <- c("pffr", class(m))
  }
  return(m)
}# end pffr()	