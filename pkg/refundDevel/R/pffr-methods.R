# methods for pffr-objects
# 
# 
# Author: fabians
# 16.08.2011, 13:01:24
###############################################################################

#' Prediction for penalized function-on-function regression 
#' 
#'  Takes a fitted \code{pffr}-object produced by \code{\link{pffr}()} and produces 
#'  predictions given a new set of values for the model covariates or the original 
#'  values used for the model fit. Predictions can be accompanied by standard errors,
#'  based on the posterior distribution of the model coefficients. This is a wrapper
#'  function for \code{\link[mgcv]{predict.gam}()}
#' 
#' @param object a fitted \code{pffr}-object
#' @param newdata  A named list containing the values of the model covariates at which predictions are required.
#'   If this is not provided then predictions corresponding to the original data are returned. 
#'  If \code{newdata} is provided then it should contain all the variables needed for prediction, 
#'  in the format supplied to \code{pffr}, i.e., 
#'  functional predictors must be supplied as matrices with each row corresponding to one observed function.
#'  Index variables for the functional covariates are reused from the fitted model object 
#'  and cannot be supplied with newdata.  
#'  Prediction is always for the entire range of \eqn{Y(t)} as defined in the original fit. 
#' @param reformat logical, defaults to TRUE. Should predictions be returned in matrix form (default) or 
#' in the long vector shape returned by \code{predict.gam()}? 
#' @param type see \code{\link[mgcv]{predict.gam}()} for details. 
#'  Note that \code{type == "lpmatrix"} will force \code{reformat} to FALSE.
#' @param se.fit see \code{\link[mgcv]{predict.gam}()}
#' @param ...  additional arguments passed on to \code{\link[mgcv]{predict.gam}()}
#' @seealso \code{\link[mgcv]{predict.gam}()}
#' @return If \code{type == "lpmatrix"}, the design matrix for the supplied covariate values in long format.
#'  If \code{se == TRUE}, a list with entries \code{fit} and \code{se.fit} containing fits and standard errors, respectively.
#'  If \code{type == "terms"} or \code{"iterms"} each of these lists is a list of matrices of the same dimension as the response for \code{newdata}
#'  containing the linear predictor and its se for each term. 
#' @method predict pffr
#' @author Fabian Scheipl
predict.pffr <- function(object,
        newdata,
        reformat=TRUE,
        type = "link",
        se.fit = FALSE,
        ...){
    #browser()
    
    call <- match.call()
    nyindex <- object$pffr$nyindex

    ## warn if any entries in ... are not arguments for predict.gam 
    dots <- list(...)
    if(length(dots)){
        validDots <- names(formals(predict.gam))
        notUsed <- names(dots)[!(names(dots) %in% validDots)]
        if(length(notUsed))
            warning("Arguments <", paste(notUsed, collapse=", "), "> supplied but not used." )
    }
    
    
    if(!missing(newdata)){
        nobs <- nrow(as.matrix(newdata[[1]]))
        
        # check if the supplied data already has the shape expected by predict.gam
        # and dispatch immediately if so (need this so summary works as expected!)
        if(!(all(names(newdata) %in% names(object$model))) | 
          !(paste0(object$pffr$yindname,".vec") %in% names(newdata))){
            # check lengths
            stopifnot(length(unique(sapply(newdata, function(x) ifelse(is.matrix(x), nrow(x), length(x))))) ==1)
            #        #FIXME: better leave this check to predict.gam....  
            #        covnames <- mapply(gsub, 
            #                pattern=c(".[st]mat$"), 
            #                replacement="", x=unique(unlist(sapply(object$smooth, function(x) x$term))))
            #        covnames <- unique(covnames[covnames != paste(object$pffr$yindname, ".vec", sep="")])
            #        stopifnot(all(covnames %in% names(newdata)))
            
            
            #get newdata into the shape expected by predict gam:    
            gamdata <- list()
            #y-index
            gamdata[[paste(object$pffr$yindname, ".vec", sep="")]] <- rep(object$pffr$yind, times=nobs)
            
            # which covariates occur in which terms?
            varmap <- sapply(names(object$pffr$labelmap), function(x) all.vars(formula(paste("~", x))))
            
            # don't include response
            covnames <- names(newdata)[names(newdata)!=deparse(object$formula[[2]])]
            for(cov in covnames){
                #find the term(s) <cov> is associated with
                trms <- which(sapply(varmap, function(x) any(grep(paste("^",cov,"$",sep=""), x)))) 
                if(!is.null(dots$terms)) trms <- trms[names(trms) %in% dots$terms]
                if(length(trms)!=0){
                    for(trm in trms){
                        is.ff <- trm %in% object$pffr$where$where.ff
                        is.sff <- trm %in% object$pffr$where$where.sff
                        is.ffpc <- trm %in% object$pffr$where$where.ffpc
                        is.pcre <- trm %in% object$pffr$where$where.pcre
                        #if ff(X) or sff(X), generate (X.mat), X.tmat, X.smat, L.X ...
                        if(is.ff){
                            ff <- object$pffr$ff[[grep(paste(cov,"[,\\)]",sep=""), names(object$pffr$ff))]]
                            #... but don't generate new data unless <cov> is the functional covariate.
                            if(grepl(paste(cov,"\\.[st]mat",sep=""), deparse(ff$call$x))){
                                # make L-matrix for new obs:
                                L <- ff$L 
                                if(any(apply(L, 2, function(x) length(unique(x)))!=1)){
                                    stop("Error for ", names(varmap)[trm],
                                            "-- Prediction for ff-terms with varying rows in integration operator L not implememented yet.")
                                }
                                if(!is.null(ff$limits)){
                                    #TODO implement prediction with limits
                                    stop("Error for ", names(varmap)[trm],
                                            "-- Prediction for ff-terms with <limits> not implememented yet.")
                                }
                                
                                predL <- matrix(L[1,], byrow=TRUE, nrow=nrow(newdata[[cov]]), ncol=ncol(L))
                                
                                
                                gamdata[[paste(cov, ".smat", sep="")]] <- 
                                        matrix(ff$xind, byrow=TRUE, ncol=length(ff$xind), nrow=nobs*nyindex)
                                gamdata[[paste(cov, ".tmat", sep="")]] <- 
                                        matrix(rep(ff$yind, times=nobs), ncol=length(ff$xind), nrow=nobs*nyindex)
                                gamdata[[paste("L.", cov, sep="")]] <-  
                                        (predL*newdata[[cov]])[rep(1:nobs, each=nyindex),]
                            }
                        }
                        if(is.sff){
                            sff <- object$pffr$ff[[grep(paste(cov,"[,\\)]",sep=""), names(object$pffr$ff))]]
                            #... but don't generate new data unless <cov> is the functional covariate.
                            if(grepl(paste(cov,"\\.[st]mat",sep=""), deparse(sff$call$x))){
                                # make L-matrix for new obs:
                                L <- sff$L 
                                if(any(apply(L, 2, function(x) length(unique(x)))!=1)){
                                    stop("Error for ", names(varmap)[trm],
                                            "-- Prediction for sff-terms with varying rows in integration operator L not implememented yet.")
                                } 
                                predL <- matrix(L[1,], byrow=TRUE, nrow=nrow(newdata[[cov]]), ncol=ncol(L)) 
                                
                                gamdata[[paste(cov, ".mat", sep="")]] <- newdata[[cov]][rep(1:nobs, e=nyindex),]
                                gamdata[[paste(cov, ".smat", sep="")]] <- 
                                        matrix(sff$xind, byrow=TRUE, ncol=length(sff$xind), nrow=nobs*nyindex)
                                gamdata[[paste(cov, ".tmat", sep="")]] <- 
                                        matrix(rep(sff$yind, times=nobs), ncol=length(sff$xind), nrow=nobs*nyindex)
                                gamdata[[paste("L.", cov, sep="")]] <-  predL[rep(1:nobs, e=nyindex),]
                            }
                        }
                        if(is.pcre){
                            pcre <- object$pffr$pcre[[grep(cov, names(object$pffr$pcre))]]
                            gamdata[[paste(cov, ".vec", sep="")]] <- rep(pcre$id, each=nyindex)
                            for(nm in colnames(pcre$efunctions)){
                                gamdata[[nm]] <- pcre$efunctions[rep(1:nyindex, times=nobs), nm]
                            }
                        }
                        if(is.ffpc){
                            ffpc <- object$pffr$ffpc[[grep(paste(cov,"[,\\)]",sep=""), names(object$pffr$ffpc))]]
                            
                            # Xc' = Phi xi' + error --> get loadings for new data:
                            Xct <- t(newdata[[cov]]) - as.vector(ffpc$meanX)
                            xiMat <- t(qr.coef(qr(ffpc$PCMat), Xct))
                            colnames(xiMat) <- paste(make.names(cov),".PC", 1:ncol(xiMat), sep="")
                            xiMat <- xiMat[rep(1:nobs, each=nyindex), ]
                            for(nm in colnames(xiMat)){
                                gamdata[[nm]] <- xiMat[,nm] 
                            }
                        }
                        if(!(is.ff | is.sff | is.ffpc | is.pcre)) {
                            #just repeat each entry nyindex-times to correspond to vec(<Response>)  
                            gamdata[[cov]] <- drop(newdata[[cov]])[rep(1:nobs, each=nyindex)]
                        }
                        
                    } 
                }
            }
            gamdata <- list2df(gamdata)
            call[["newdata"]] <- gamdata
        } 
    } else {
        call$newdata <- eval(call$newdata)
        nobs <- object$pffr$nobs
    }
    
 
    
    
    
    #call predict.gam
    call[[1]] <- mgcv::predict.gam
    call$object <- as.name("object")
    ret <- eval(call)
    
    if(type=="lpmatrix" && reformat){
        reformat <- FALSE
        warning("Setting reformat to FALSE for type=\"lpmatrix\".")
    }
    
    #reformat into matrices with same shape as <Response>

    if(reformat){
        
        if(missing(newdata) && !is.null(object$pffr$missingind)){
            #pad with NAs at the appropriate locations so that fits are nobs x nyindex:
            insertNA <- function(x){
                if(length(x) != nobs*object$pffr$nyindex){
                    tmp <- rep(NA, nobs*object$pffr$nyindex)
                    tmp[-object$pffr$missingind] <- x
                    return(tmp)
                } else {
                    return(x)
                }
            }
        } else insertNA <- function(x) return(x)
        
        
        if(se.fit){
            if(type %in% c("terms", "iterms")){
                ret <- lapply(ret, function(x)
                            do.call(list,
                                    sapply(1:ncol(x), function(i){
                                                #browser()
                                                d <- list(I(matrix(insertNA(x[,i]), nrow=nobs, ncol=object$pffr$nyindex, byrow=TRUE)))
                                                names(d)  <- colnames(x)[i]
                                                return(d)
                                            })))
                                               
            } else {
                ret <- lapply(ret, function(x) matrix(insertNA(x), nrow=nobs, ncol=object$pffr$nyindex, byrow=TRUE))  
            } 
        } else {
            if(type %in% c("terms", "iterms")){
                ret <- do.call(list, sapply(1:ncol(ret), function(i){
                                    #browser()
                                    d <- list(I(matrix(insertNA(ret[,i]), nrow=nobs, ncol=object$pffr$nyindex, byrow=TRUE)))
                                    names(d)  <- colnames(ret)[i]
                                    return(d)
                                }))
            } else ret <- matrix(insertNA(ret), nrow=nobs, ncol=object$pffr$nyindex, byrow=TRUE)
        }
    }
    return(ret)
}

#' Obtain model matrix for a pffr fit
#' 
#' @param object a fitted \code{pffr}-object 
#' @param ... other arguments, passed to \code{\link[mgcv]{predict.gam}}.
#' 
#' @return A model matrix
#' @method model.matrix pffr
#' @author Fabian Scheipl
model.matrix.pffr <- function (object, ...) 
{
    if (!inherits(object, "pffr")) 
        stop("`object' is not of class \"pffr\"")
    predict(object, type = "lpmatrix", reformat=FALSE, ...)
}

#' Obtain residuals for a pffr fit
#' 
#' @param object a fitted \code{pffr}-object 
#' @param reformat logical, defaults to TRUE. Should residuals be returned in \code{n x yindex} matrix form (default) or 
#' in the long vector shape returned by \code{resid.gam()}? 
#' @param ... other arguments, passed to \code{\link[mgcv]{residuals.gam}}.
#' 
#' @return A matrix or vector of residuals
#' @method residuals pffr
#' @author Fabian Scheipl
residuals.pffr <- function (object, reformat=TRUE, ...) 
{
    if (!inherits(object, "pffr")) 
        stop("`object' is not of class \"pffr\"")
    ret <- mgcv:::residuals.gam(object, ...)
    if(reformat){
        if(!(length(ret)==object$pffr$nobs*object$pffr$nyindex)){
            tmp <- rep(NA, object$pffr$nobs*object$pffr$nyindex)
            tmp[-object$pffr$missingind] <- ret
            ret <- tmp
        }
        ret <- matrix(ret, nrow=object$pffr$nobs, ncol=object$pffr$nyindex, byrow=TRUE)  
    } 
    return(ret)
}

#' Obtain fitted values for a pffr fit
#' 
#' Returns the linear predictor for the model data. See \code{\link{predict.pffr}}
#' for alternative options.
#' 
#' @param object a fitted \code{pffr}-object 
#' @param reformat logical, defaults to TRUE. Should residuals be returned in \code{n x yindex} matrix form (default) or 
#' in the long vector shape returned by \code{fitted.default()}? 
#' @param ... not used
#' 
#' @return A matrix or vector of fitted values
#' @method fitted pffr
#' @author Fabian Scheipl
fitted.pffr <- function (object, reformat=TRUE, ...) 
{
    if (!inherits(object, "pffr")) 
        stop("`object' is not of class \"pffr\"")
    ret <- object$fitted.values
    if(reformat){
        if(!(length(ret)==object$pffr$nobs*object$pffr$nyindex)){
            tmp <- rep(NA, object$pffr$nobs*object$pffr$nyindex)
            tmp[-object$pffr$missingind] <- ret
            ret <- tmp
        }
        ret <- matrix(ret, nrow=object$pffr$nobs, ncol=object$pffr$nyindex, byrow=TRUE)  
    } 
    return(ret)
}

#' Plot a pffr fit
#' 
#' Plot a fitted pffr-object. Simply dispatches to \code{\link[mgcv]{plot.gam}}.
#' 
#' @param x a fitted \code{pffr}-object 
#' @param ... arguments handed over to \code{\link[mgcv]{plot.gam}}
#' 
#' @return This function only generates plots.
#' @method plot pffr
#' @author Fabian Scheipl
plot.pffr <- function (x, ...) 
{
    call <- match.call()
    call[[1]] <- mgcv::plot.gam
    #drop "pffr" class and replace <object> with changed value s.t. method dispatch works without glitches
    class(x) <- class(x)[-1]
    invisible(eval(call))
}


#' Get estimated coefficients from a pffr fit
#' 
#' Returns estimated coefficient functions/surfaces \eqn{\beta(t), \beta(s,t)}
#' and estimated smooth effects \eqn{f(z), f(x,z)} or \eqn{f(x, z, t)} and their point-wise estimated standard errors. 
#' Not implemented for smooths in more than 3 dimensions.
#' 
#' The \code{seWithMean}-option corresponds to the \code{"iterms"}-option in \code{\link[mgcv]{predict.gam}}.
#' The \code{sandwich}-options works as follows: Assuming that the residual vectors \eqn{\epsilon_i(t), i=1,\dots,n} are i.i.d.
#' realizations of a mean zero Gaussian process with covariance \eqn{K(t,t')}, we can construct an estimator for  
#' \eqn{K(t,t')} from the \eqn{n} replicates of the observed residual vectors. The covariance matrix of the stacked observations
#' vec\eqn{(Y_i(t))} is then given by a block-diagonal matrix with \eqn{n} copies of the estimated \eqn{K(t,t')} on the diagonal.
#' This block-diagonal matrix is used to construct the "meat" of a sandwich covariance estimator, similar to Chen et al. (2012), 
#' see reference below.  
#' 
#' 
#' @param object a fitted \code{pffr}-object
#' @param raw logical, defaults to FALSE. If TRUE, the function simply returns \code{object$coefficients}
#' @param se logical, defaults to TRUE. Return estimated standard error of the estimates?
#' @param freq logical, defaults to FALSE. If FALSE, use posterior variance \code{object$Vp} for variability estimates, 
#'  else use \code{object$Ve}. See \code{\link[mgcv]{gamObject}}
#' @param sandwich logical, defaults to FALSE. Use a Sandwich-estimator for approximate variances? See Details. 
#' 	THIS IS AN EXPERIMENTAL FEATURE, USE A YOUR OWN RISK.
#' @param seWithMean logical, defaults to TRUE. Include uncertainty about the intercept/overall mean in  standard errors returned for smooth components?
#' @param n1 see below
#' @param n2 see below
#' @param n3 \code{n1, n2, n3} give the number of gridpoints for 1-/2-/3-dimensional smooth terms 
#' used in the marginal equidistant grids over the range of the covariates at which the estimated effects are evaluated. 
#' @param Ktt (optional) an estimate of the covariance operator of the residual process \eqn{\epsilon_i(t) \sim N(0, K(t,t'))}, 
#' evaluated on \code{yind} of \code{object}. If not supplied, this is estimated from the crossproduct matrices of the
#' observed residual vectors. Only relevant for sandwich CIs.
#' @param ... other arguments, not used.
#' 
#' @return If \code{raw==FALSE}, a list containing \itemize{
#'  \item \code{pterms} a matrix containing the parametric / non-functional coefficients (and, optionally, their se's)
#'  \item \code{smterms} a named list with one entry for each smooth term in the model. Each entry contains
#'     \itemize{
#'          \item \code{coef} a matrix giving the grid values over the covariates, the estimated effect (and, optionally, the se's). 
#'                          The first covariate varies the fastest.
#'          \item \code{x, y, z} the unique gridpoints used to evaluate the smooth/coefficient function/coefficient surface
#'          \item \code{xlim, ylim, zlim} the extent of the x/y/z-axes
#'          \item \code{xlab, ylab, zlab} the names of the covariates for the x/y/z-axes
#'          \item \code{dim} the dimensionality of the effect
#'          \item \code{main} the label of the smooth term (a short label, same as the one used in \code{summary.pffr})
#' }} 
#' @references Chen H., Wang Y., Paik C.M., Choi A. (2012). 
#' A marginal approach to reduced-rank penalized spline smoothing for multilevel data. 
#' \emph{Journal of the American Statistical Association}, under revision. 
#' \url{http://www.columbia.edu/~yw2016/Marginal Spline6.pdf}
#' @method coef pffr
#' @seealso \code{\link[mgcv]{plot.gam}}, \code{\link[mgcv]{predict.gam}} which this routine is
#'   based on.
#' @author Fabian Scheipl
coef.pffr <- function(object, raw=FALSE, se=TRUE, freq=FALSE, sandwich=FALSE, 
        seWithMean=TRUE, n1=100, n2=40, n3=20, Ktt=NULL, ...){
    if(raw){
        return(object$coefficients)  
    } else {
        getCoefs <- function(i){
            ## this constructs a grid over the range of the covariates
            ## and returns estimated values on this grid, with 
            ## by-variables set to 1
            ## cf. mgcv:::plots.R (plot.mgcv.smooth etc..) for original code
            safeRange <- function(x){
                if(is.factor(x)) return(c(NA, NA))
                return(range(x, na.rm=TRUE))
            }
            
            makeDataGrid <- function(trm){
                #generate grid of values in range of original data
                if(trm$dim==1){
                    x <- get.var(trm$term, object$model)
                    xg <- if(is.factor(x)) {
                                unique(x)
                            } else seq(min(x), max(x), length=n1)
                    d <- data.frame(xg)
                    colnames(d) <- trm$term
                    attr(d, "xm") <- xg
                }
                if(trm$dim > 1){
                    
                    ng <- ifelse(trm$dim==2, n2, n3) 
                    
                    varnms <- trm$term
                    
                    x <- get.var(trm$term[1], object$model)
                    xg <- if(is.factor(x)) {
                                unique(x)
                            } else seq(min(x), max(x),length=ng)
                    y <- get.var(trm$term[2], object$model)
                    yg <- if(is.factor(y)) {
                                unique(y)
                            } else seq(min(y), max(y),length=ng)
                    if(length(varnms)==2){
                        d <- expand.grid(xg, yg)
                        attr(d, "xm") <- xg
                        attr(d, "ym") <- yg    
                    } else {
                        z <- get.var(trm$term[3], object$model)
                        zg <- if(is.factor(z)) {
                                    unique(z)
                                } else seq(min(z), max(z), length=ng)
                        d <- expand.grid(xg, yg, zg)
                        attr(d, "xm") <- xg
                        attr(d, "ym") <- yg
                        attr(d, "zm") <- zg
                    }
                    colnames(d) <- varnms 
                    
                }
                if(trm$by!="NA"){
                    d$by <- 1
                    colnames(d) <- c(head(colnames(d),-1), trm$by)
                } 
                return(d)
            }
            
            getP <- function(trm, d){
                #return an object similar to what plot.mgcv.smooth etc. return 
                X <- PredictMat(trm, d)
                P <- if(trm$dim==1){
                            list(x=attr(d, "xm"), xlab=trm$term, xlim=safeRange(attr(d, "xm")))
                        } else {
                            varnms <-  trm$term
                            if(length(varnms) == 2){
                                list(x=attr(d, "xm"), y=attr(d, "ym"), xlab=varnms[1], ylab=varnms[2],
                                        ylim=safeRange(attr(d, "ym")), xlim=safeRange(attr(d, "xm")))
                            } else {
                                if(trm$dim==3){
                                    list(x=attr(d, "xm"), y=attr(d, "ym"), z=attr(d, "zm"), 
                                            xlab=varnms[1], ylab=varnms[2], zlab=varnms[3],
                                            ylim=safeRange(attr(d, "ym")), xlim=safeRange(attr(d, "xm")), zlim=safeRange(attr(d, "zm")))
                                }
                            }
                        }
                trmind <- trm$first.para:trm$last.para
                P$value <- X%*%object$coefficients[trmind]
                P$coef <- cbind(d, "value"=P$value)
                if(se){
                    # use seWithMean if possible: 
                    if(seWithMean & attr(trm,"nCons")>0){
                        cat("using seWithMean for ", trm$label,".\n")
                        X1 <- matrix(object$cmX,nrow(X),ncol(object$Vp),byrow=TRUE)
                        meanL1 <- trm$meanL1
                        if (!is.null(meanL1)) X1 <- X1 / meanL1
                        X1[,trmind] <- X
                        P$se <- sqrt(rowSums((X1%*%covmat)*X1))
                    } else {
                        P$se <- sqrt(rowSums((X%*%covmat[trmind, trmind])*X))
                    }
                    P$coef <- cbind(P$coef, se=P$se)
                }
                
                P$dim <- trm$dim
                return(P)
            }
            
            trm <- object$smooth[[i]]
            if(trm$dim > 3){
                warning("can't deal with smooths with more than 3 dimensions, returning NULL for ", 
                        shrtlbls[names(object$smooth)[i] == unlist(object$pffr$labelmap)])
                return(NULL)
            }
            
            d <- makeDataGrid(trm)
            P <- getP(trm, d)
            
            #browser()
            # get proper labeling
            P$main <- shrtlbls[names(object$smooth)[i] == unlist(object$pffr$labelmap)]  
            which <- match(names(object$smooth)[i], object$pffr$labelmap)
            if(which %in% object$pffr$where$where.ff){
                which.ff <- which(object$pffr$where$where.ff == which)
                P$ylab <- object$pffr$yindname
                xlab <- deparse(as.call(formula(paste("~",names(object$pffr$ff)[which.ff]))[[2]])$xind)
                if(xlab=="NULL") xlab <- "xindex"
                P$xlab <- xlab
            }
            if(which %in% object$pffr$where$where.sff){
                which.sff <- which(object$pffr$where$where.sff == which)
                P$ylab <- object$pffr$yindname
                xlab <- deparse(as.call(formula(paste("~",names(object$pffr$ff)[which.sff]))[[2]])$xind)
                if(xlab=="NULL") xlab <- "xindex"
                P$xlab <- xlab
                P$zlab <- gsub(".mat$", "", object$pffr$ff[[which.sff]]$xname) 
            }
            
            return(P)
        }
        
        bread <- if(freq){
                    object$Ve
                } else {
                    object$Vp
                }
        if(sandwich){
            X <- predict(object, type = "lpmatrix", reformat=FALSE)
            bread <- bread/object$sig2
            res <- residuals(object)
            if(is.null(Ktt)){
                # get estimate of Cov(eps_i(t)) = K(t,t')
                stopifnot(require(Matrix))
                Ktt <- Reduce("+",  lapply(1:nrow(res), function(i) tcrossprod(res[i,])))/nrow(res)
            }
            #Chen/Wang, Sec. 2.1: M = X' V^-1 (Y-eta)(Y-eta)' V^-1 X  with V ^-1 = diag(sigma^-2) 
            # since the estimate is under working independence.... 
            meat <- (t(X)%*%kronecker(Diagonal(nrow(res)), Ktt))%*%X / (object$scale^2)
            covmat <- as.matrix(bread %*% meat %*% bread)
        } else {
            covmat <- bread
        }  
        
        ret <- list()
        smind <- unlist(sapply(object$smooth, function(x){
                            seq(x$first.para, x$last.para)
                        }))
        ret$pterms <- cbind(value=object$coefficients[-smind])
        if(se) ret$pterms <- cbind(ret$pterms, se=sqrt(diag(covmat)[-smind]))
        
        shrtlbls <- getShrtlbls(object)
        
        ret$smterms <- lapply(1:length(object$smooth), getCoefs)         
        names(ret$smterms) <- sapply(seq_along(ret$smterms), function(i){
                    ret$smterms[[i]]$main
                })
        return(ret)    
    }
}

#' Summary for a pffr fit
#' 
#' Take a fitted \code{pffr}-object and produce summaries from it.
#' See \code{\link[mgcv]{summary.gam}()} for details.
#'  
#' @param object a fitted \code{pffr}-object 
#' @param ... see \code{\link[mgcv]{summary.gam}()} for options.
#' 
#' @return A list with summary information, see \code{\link[mgcv]{summary.gam}()}
#' @method summary pffr
#' @author Fabian Scheipl, adapted from \code{\link[mgcv]{summary.gam}()} by Simon Wood, Henric Nilsson
summary.pffr <- function (object, ...) {
    call <- match.call()
    call[[1]] <- mgcv::summary.gam
    ## drop "pffr" class and replace <object> with changed value s.t. method dispatch works without glitches
    ## if we don't do this, summary.gam will call predict on the object if n>3000 & freq==TRUE 
    ## and this predict-call gets dispatched to predict.pffr which dispatches back
    ## to predict.gam. Somewhere along the way an index variable get's lost and
    ## shit breaks down. 
    class(object) <- class(object)[!(class(object) %in% "pffr")]
    call$object <- as.name("object")
    ret <- eval(call)
    
    ret$formula <- object$pffr$formula
    
    # make short labels for display
    shrtlbls <- getShrtlbls(object)
    
    rownames(ret$s.table) <- sapply(rownames(ret$s.table), 
            function(x){
                shrtlbls[pmatch(x, unlist(object$pffr$labelmap))]     
            })
    class(ret) <- c("summary.pffr", class(ret))
    ret$n  <- paste(ret$n, "(", object$pffr$nobs," x ", object$pffr$nyindex, ")", sep="") 
    return(ret)
}

#' Print method for summary of a pffr fit
#' 
#' Pretty printing for a \code{summary.pffr}-object.
#' See \code{\link[mgcv]{print.summary.gam}()} for details.
#'  
#' @param x a fitted \code{pffr}-object 
#' @param digits controls number of digits printed in output. 
#' @param signif.stars Should significance stars be printed alongside output?
#' @param ... not used 
#' 
#' @return A \code{\link{summary.pffr}} object
#' @method print summary.pffr
#' @author Fabian Scheipl, adapted from \code{\link[mgcv]{print.summary.gam}()} by Simon Wood, Henric Nilsson
print.summary.pffr <- function(x, digits = max(3, getOption("digits") - 3), 
        signif.stars = getOption("show.signif.stars"), ...){
# mostly identical to print.summary.gam
    print(x$family)
    cat("Formula:\n")
    print(x$formula)
    if (length(x$p.coeff)>0)
    { cat("\nConstant coefficients:\n")
        printCoefmat(x$p.table, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)
    }
    cat("\n")
    if(x$m>0)
    { cat("Smooth terms & functional coefficients:\n")
        printCoefmat(x$s.table, digits = digits, signif.stars = signif.stars, has.Pvalue = TRUE, na.print = "NA",cs.ind=1, ...)
    }
    cat("\nR-sq.(adj) = ",formatC(x$r.sq,digits=3,width=5))
    if (length(x$dev.expl)>0) cat("   Deviance explained = ",formatC(x$dev.expl*100,digits=3,width=4),"%\n",sep="")
    
    if (!is.null(x$method)&&!(x$method%in%c("PQL","lme.ML","lme.REML")))  
        cat(x$method," score = ",formatC(x$sp.criterion,digits=5),sep="")
    
    cat("  Scale est. = ",formatC(x$scale,digits=5,width=8,flag="-"),"  n = ",x$n,"\n",sep="")
    invisible(x)
}



