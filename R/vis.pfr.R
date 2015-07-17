#' Visualization of PFR objects
#'
#' Produces perspective or contour plot views of an estimated surface corresponding
#' smooths over two or more dimensions. Alternatively plots \dQuote{slices} of the
#' estimated surface or estimated second derivative surface with one of its arguments fixed.
#' Corresponding twice-standard error \dQuote{Bayesian} confidence bands are
#' constructed using the method in Marra and Wood (2012). See the details.
#' 
#' @param object an \code{pfr} object, produced by \code{\link{pfr}}
#' @param select index for the smooth term to be plotted, according to its position
#'   in the model formula (and in \code{object$smooth}). Not needed if only one
#'   multivariate term is present.
#   Can be entered as an integer index, or as a
#   character string containing the name of the functional predictor to be plotted.
#' @param xval a number in the range of functional predictor to be plotted.  The surface will be plotted
#'   with the first argument of the estimated surface fixed at this value
#' @param tval a number in the domain of the functional predictor to be plotted.  The surface will be
#'   plotted with the second argument of the estimated surface fixed at this value. Ignored if \code{xval}
#'   is specified.
#' @param deriv2 logical; if \code{TRUE}, plot the estimated second derivative surface along with
#'   Bayesian confidence bands.  Only implemented for the "slices" plot from either \code{xval} or
#'   \code{tval} being specified
#' @param theta numeric; viewing angle; see \code{\link{persp}}
#' @param plot.type one of \code{"contour"} (to use \code{\link{levelplot}}) or \code{"persp"}
#'   (to use \code{\link{persp}}).  Ignored if either \code{xval} or \code{tval} is specified
#' @param ticktype how to draw the tick marks if \code{plot.type="persp"}.  Defaults to \code{"detailed"}
#' @param ... other options to be passed to \code{\link{persp}}, \code{\link{levelplot}}, or
#'   \code{\link{plot}}
#' 
#' @details The confidence bands used when plotting slices of the estimated surface or second derivative
#' surface are the ones proposed in Marra and Wood (2012).  These are a generalization of the "Bayesian"
#' intevals of Wahba (1983) with an adjustment for the uncertainty about the model intercept. The
#' estimated covariance matrix of the model parameters is obtained from assuming a particular Bayesian
#' model on the parameters.
#' 
#' @return Simply produces a plot
#' @references
#' McLean, M. W., Hooker, G., Staicu, A.-M., Scheipl, F., and Ruppert, D. (2014). Functional
#' generalized additive models. \emph{Journal of Computational and Graphical Statistics}, \bold{23(1)},
#' pp. 249-269.  Available at \url{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3982924}.
#'
#' Marra, G., and Wood, S. N. (2012) Coverage properties of confidence intervals for generalized
#' additive model components. \emph{Scandinavian Journal of Statistics}, \bold{39(1)}, pp. 53--74.
#'
#' Wabha, G. (1983) "Confidence intervals" for the cross-validated smoothing spline. \emph{Journal of the
#' Royal Statistical Society, Series B}, \bold{45(1)}, pp. 133--150.
#' 
#' @author Mathew W. McLean \email{mathew.w.mclean@@gmail.com}
#' @seealso \code{\link{vis.gam}}, \code{\link{plot.gam}}, \code{\link{pfr}}, \code{\link{persp}},
#'   \code{\link{levelplot}}
#' @importFrom mgcv vis.gam predict.gam
#' @importFrom lattice levelplot
#' @importFrom graphics persp
#' @export
#' 
#' @examples
#' ################# DTI Example #####################
#' data(DTI)
#' 
#' ## only consider first visit and cases (since no PASAT scores for controls),
#' ## and remove missing data
#' DTI <- DTI[DTI$visit==1 & DTI$case==1 & complete.cases(DTI$cca),]
#' 
#' ## Fit the PFR using FA measurements along corpus
#' ## callosum as functional predictor with PASAT as response
#' ## using 8 cubic B-splines for each marginal bases with
#' ## third order marginal difference penalties.
#' ## Specifying gamma>1 enforces more smoothing when using GCV
#' ## to choose smoothing parameters
#' fit <- pfr(pasat ~ af(cca, basistype="te", k=c(8,8), m=list(c(2,3),c(2,3)), bs="ps"),
#'            method="GCV.Cp", gamma=1.2, data=DTI)
#'
#' ## contour plot of the fitted surface
#' vis.pfr(fit, plot.type='contour')
#'
#' ## similar to Figure 5 from McLean et al.
#' ## Bands seem too conservative in some cases
#' xval <- runif(1, min(fit$pfr$ft[[1]]$Xrange), max(fit$pfr$ft[[1]]$Xrange))
#' tval <- runif(1, min(fit$pfr$ft[[1]]$xind), max(fit$pfr$ft[[1]]$xind))
#' par(mfrow=c(2, 2))
#' vis.pfr(fit, deriv2=FALSE, xval=xval)
#' vis.pfr(fit, deriv2=FALSE, tval=tval)
#' vis.pfr(fit, deriv2=TRUE, xval=xval)
#' vis.pfr(fit, deriv2=TRUE, tval=tval)

vis.pfr=function(object, select=1, xval = NULL, tval = NULL, deriv2 = FALSE, theta = 50,
                  plot.type = "persp", ticktype = "detailed", ...){
  
  trm.dim <- sapply(object$smooth, function(x) x$dim)
  n.mv <- sum(trm.dim>1)
  if (!all(length(object$pfr), length(object$smooth), n.mv>0))
    stop('Model contains no smooth terms for visualization')
  
  if (!(object$smooth[[select]]$dim > 1)) {
    if (select==1 & n.mv==1) {
      select <- which(trm.dim>1)
    } else {
      stop('Term indicated by \"select\" is not a multivariate smooth')
    }
  }
  
  
  # ttypes are the term types for object$smooth; fttypes for object$pfr$ft
  #ftlist  <- c("lf", "af", "lf.vd", "peer") # This vector identifies term types that end up in ft
  ttypes  <- object$pfr$termtype[object$pfr$termtype != "par"]
  #fttypes <- ttypes[ttypes %in% ftlist]
  #ftnames <- sapply(strsplit(sapply(object$pfr$ft, function(x) x$tindname),
  #                           ".", fixed=TRUE), function(x) x[1])
  #stopifnot(length(ttypes)  == length(object$smooth))
  #stopifnot(length(fttypes) == length(object$pfr$ft))
  #
  #if (is.character(select)) {
  #  af.ind <- which(ftnames == select & fttypes=="af")
  #  if (!length(af.ind))  stop("Unrecognized \"select\" index")
  #  if (length(af.ind)>1) stop("Multiple \"select\" indicess matched?")
  #} else if (is.numeric(select)) {
  af.ind <- which(ttypes%in%c("lf","af","peer","fpc","lf.vd"))[select]
  #  select <- ftnames[af.ind]
  #  if (is.na(af.ind)) stop("select entry exceeds the number of af terms")
  #} else {
  #  stop("Uncrecognized select type")
  #}
  
  # Corresponding coefficient indices
  smooth.i <- object$smooth[[select]]
  tind <- 1:length(object$coefficients) %in%
    c(smooth.i$first.para:smooth.i$last.para)
  tvar <- modify_nm(smooth.i$term[2])
  mtitle <- if (!is.null(smooth.i$QT))
    bquote(paste(hat(F),'(p,t),   p=',hat(G)[t],'(',.(tvar),')',sep=''))
  else
    bquote(paste(hat(F),'(',.(tvar),',t)',sep=''))
  
  if(!(length(tval)+length(xval))) {
    # Plot entire surface (not slices)
    temp <- list()
    temp[[smooth.i$by]] <- 1
    if(plot.type=='persp' ) {
      vis.gam(object, view=paste0(tvar, c(".omat", ".tmat")), cond=temp,
              ticktype=ticktype, theta=theta, contour.col=rev(heat.colors(100)),
              xlab=tvar, ylab='\nt', zlab='', main=mtitle, ...)
      
    } else if (plot.type=='contour'){
      # not making use of vis.gam because want colour key/legend
      est <- coef(object, select=select)
      ylab=ifelse(is.null(smooth.i$QT), tvar, 'p')
      names(est)[1:2] <- c("t", "x")
      
      lattice::levelplot(value~t*x, data=est, contour=TRUE, labels=TRUE, 
                         pretty=TRUE, xlab='t', ylab=ylab,
                         col.regions=rev(heat.colors(100)),
                         main=as.expression(mtitle), ...)
    }
  } else {
    # Plot slice of surface
    
    if(length(xval)+length(tval)>1){
      warning('Specify only one value for either xval or tval.  Only first specified value will be used')
      if(length(xval)){
        xval=xval[1]
        tval=NULL
      } else{
        tval=tval[1]
      }
    }
    
    if (length(xval)) {
      # x fixed
      nfine <- ndefault(object$pfr$ft[[af.ind]]$xind)
      trange <- range(object$pfr$ft[[af.ind]]$xind)
      tvals <- seq(trange[1],trange[2],l=nfine)
      xvals <- rep(xval,l=nfine)
      xlab='t'
      x=tvals
      if(!deriv2) {
        main=bquote(paste(hat(F)(.(round(xval,3)),t),' by t',sep=''))
        ylab=bquote(paste(hat(F)(.(round(xval,3)),t),sep=''))
      } else {
        main=bquote(paste(frac(partialdiff^2,partialdiff*.(tvar)^2)*hat('F')(.(tvar),t),
                          '|',phantom()[.(tvar)==.(round(xval,3))],' by t',sep=''))
        ylab=bquote(paste(frac(partialdiff^2,partialdiff*.(tvar)^2)*hat('F')(.(tvar),t),
                          '|',phantom()[.(tvar)==.(round(xval,3))],sep=''))
      }
    } else{
      # t fixed
      nfine <- ndefault(object$model[[smooth.i$term[2]]])
      Xrange <- object$pfr$ft[[af.ind]]$Xrange
      tvals <- rep(tval,nfine)
      xvals <- seq(Xrange[1],Xrange[2],l=nfine)
      xlab='x'
      x=xvals
      if (!deriv2) {
        main=bquote(paste(hat(F)(.(tvar),.(round(tval,3))),' by ',.(tvar),sep=''))
        ylab=bquote(paste(hat(F)(.(tvar),.(round(tval,3))),sep=''))
      } else {
        main=bquote(paste(frac(partialdiff^2,partialdiff*.(tvar)^2)*hat('F')(.(tvar),.(round(tval,3))),' by ',.(tvar),sep=''))
        ylab=bquote(paste(frac(partialdiff^2,partialdiff*.(tvar)^2)*hat('F')(.(tvar),.(round(tval,3))),sep=''))
      }
      
    }
    newdata <- list()
    newdata[[paste(tvar,'.omat',sep='')]] <- xvals
    newdata[[paste(tvar,'.tmat',sep='')]] <- tvals
    newdata[[paste('L.',tvar,sep='')]] <- rep(1,l=length(xvals))
    varnames <- all.vars(terms(object))
    varnames <- varnames[!(varnames %in% c(paste('L.',tvar,sep=''),paste(tvar,'.tmat',sep=''),paste(tvar,'.omat',sep='')))]
    varnames <- varnames[-1]
    if(length(varnames)){
      for(i in 1:length(varnames)){
        newdata[[varnames[i]]] <- rep(object$pfr$datameans[names(object$pfr$datameans)==varnames[i]],l=length(xvals))#matrix(object$pfr$datameans[names(object$pfr$datameans)==varnames[i]],nr=1,nc=length(newX))
      }
    }
    
    lpmat <- predict.gam(object,newdata=newdata,type='lpmatrix')[,tind]
    if(deriv2){
      eps <- 1e-7 ## finite difference interval
      newdata[[paste(tvar,'.tmat',sep='')]] <- tvals + eps
      X1 <- predict.gam(object,newdata=newdata,type='lpmatrix')[,tind]
      
      newdata[[paste(tvar,'.tmat',sep='')]] <- tvals + 2*eps
      X2 <- predict.gam(object,newdata=newdata,type='lpmatrix')[,tind]
      
      lpmat <- (X2-2*X1+lpmat)/eps^2
    }
    
    preds <- sum(object$cmX)+lpmat%*%object$coef[tind]
    se.preds <- rowSums(lpmat%*%object$Vp[tind,tind]*lpmat)^.5
    ucl <- preds+2*se.preds
    lcl <- preds-2*se.preds
    ylim <- range(ucl,lcl,preds)
    
    par(mar=c(5.1,5.7,4.1,0.5))
    plot(x=x,preds,ylim=ylim,type='l',main=main,ylab=ylab,xlab=xlab,...)
    lines(x,ucl,col=2,lty=2)
    lines(x,lcl,col=2,lty=2)
    if(deriv2)
      abline(h=sum(object$cmX),lty=2)
  }
}

modify_nm <- function(x) {
  x <- gsub("\\.omat", "", x)
  x <- gsub("tmat", "argvals", x)
  x
}