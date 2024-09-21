#' Visualization of FGAM objects
#'
#' Produces perspective or contour plot views of an estimated surface corresponding to \code{{af}}
#' terms fit using \code{{fgam}} or plots \dQuote{slices} of the estimated surface or estimated
#' second derivative surface with one of its arguments fixed and corresponding twice-standard error
#' \dQuote{Bayesian} confidence bands constructed using the method in Marra and Wood (2012).  See the details.
#' @param object an \code{fgam} object, produced by \code{{fgam}}
#' @param af.term character; the name of the functional predictor to be plotted.  Only important
#' if multiple \code{af} terms are fit.  Defaults to the first \code{af} term in \code{object$call}
#' @param xval a number in the range of functional predictor to be plotted.  The surface will be plotted
#' with the first argument of the estimated surface fixed at this value
#' @param tval a number in the domain of the functional predictor to be plotted.  The surface will be
#' plotted with the second argument of the estimated surface fixed at this value. Ignored if \code{xval}
#' is specified
#' @param deriv2 logical; if \code{TRUE}, plot the estimated second derivative surface along with
#' Bayesian confidence bands.  Only implemented for the "slices" plot from either \code{xval} or
#' \code{tval} being specified
#' @param theta numeric; viewing angle; see \code{{persp}}
#' @param plot.type one of \code{"contour"} (to use \code{{levelplot}}) or \code{"persp"}
#' (to use \code{{persp}}).  Ignored if either \code{xval} or \code{tval} is specified
#' @param ticktype how to draw the tick marks if \code{plot.type="persp"}.  Defaults to \code{"detailed"}
#' @param ... other options to be passed to \code{{persp}}, \code{{levelplot}}, or
#' \code{{plot}}
#' @details The confidence bands used when plotting slices of the estimated surface or second derivative
#' surface are the ones proposed in Marra and Wood (2012).  These are a generalization of the "Bayesian"
#' intervals of Wahba (1983) with an adjustment for the uncertainty about the model intercept. The
#' estimated covariance matrix of the model parameters is obtained from assuming a particular Bayesian
#' model on the parameters.
#' @return Simply produces a plot
#' @references McLean, M. W., Hooker, G., Staicu, A.-M., Scheipl, F., and Ruppert, D. (2014). Functional
#' generalized additive models. \emph{Journal of Computational and Graphical Statistics}, \bold{23(1)},
#' pp. 249-269.
#'
#' Marra, G., and Wood, S. N. (2012) Coverage properties of confidence intervals for generalized
#' additive model components. \emph{Scandinavian Journal of Statistics}, \bold{39(1)}, pp. 53--74.
#'
#' Wabha, G. (1983) "Confidence intervals" for the cross-validated smoothing spline. \emph{Journal of the
#' Royal Statistical Society, Series B}, \bold{45(1)}, pp. 133--150.
#' @author Mathew W. McLean \email{mathew.w.mclean@@gmail.com}
#' @seealso \code{{vis.gam}}, \code{{plot.gam}}, \code{{fgam}}, \code{{persp}},
#' \code{{levelplot}}
#' @importFrom mgcv vis.gam predict.gam
#' @importFrom lattice levelplot
#' @importFrom graphics persp
#' @export
#' @examples
#' ################# DTI Example #####################
#' data(DTI)
#'
#' ## only consider first visit and cases (since no PASAT scores for controls)
#' y <- DTI$pasat[DTI$visit==1 & DTI$case==1]
#' X <- DTI$cca[DTI$visit==1 & DTI$case==1,]
#'
#' ## remove samples containing missing data
#' ind <- rowSums(is.na(X))>0
#'
#' y <- y[!ind]
#' X <- X[!ind,]
#'
#' ## fit the fgam using FA measurements along corpus
#' ## callosum as functional predictor with PASAT as response
#' ## using 8 cubic B-splines for each marginal bases with
#' ## third order marginal difference penalties
#' ## specifying gamma>1 enforces more smoothing when using GCV
#' ## to choose smoothing parameters
#' #fit <- fgam(y~af(X,splinepars=list(k=c(8,8),m=list(c(2,3),c(2,3)))),gamma=1.2)
#'
#' ## contour plot of the fitted surface
#' #vis.fgam(fit,plot.type='contour')
#'
#' ## similar to Figure 5 from McLean et al.
#' ## Bands seem too conservative in some cases
#' #xval <- runif(1, min(fit$fgam$ft[[1]]$Xrange), max(fit$fgam$ft[[1]]$Xrange))
#' #tval <- runif(1, min(fit$fgam$ft[[1]]$xind), max(fit$fgam$ft[[1]]$xind))
#' #par(mfrow=c(4, 1))
#' #vis.fgam(fit, af.term='X', deriv2=FALSE, xval=xval)
#' #vis.fgam(fit, af.term='X', deriv2=FALSE, tval=tval)
#' #vis.fgam(fit, af.term='X', deriv2=TRUE, xval=xval)
#' #vis.fgam(fit, af.term='X', deriv2=TRUE, tval=tval)
vis.fgam=function(object, af.term, xval = NULL, tval = NULL, deriv2 = FALSE, theta = 50,
                  plot.type = "persp", ticktype = "detailed", ...){

  if(!length(object$fgam) | !length(object$fgam$where$where.af))
    stop('Model contains no af terms for plotting')

  if(missing(af.term)){
    tnames <- names(object$model)
    af.term <- tnames[grep(paste('','.omat',sep=''),tnames)[1]]
    af.term <- strsplit(af.term,'.omat')
  }

  tnames <- object$fgam$labelmap
  ## index of desired af among all predictors
  afind <- grep(paste('te[(]',af.term,sep = ""), tnames)
  ## index of desired af among all func. predictors
  af.ind <- grep(paste('af[(]',af.term,sep = ""), names(object$fgam$ft))
  if(!length(afind))
    stop('The specified af.term is not valid')
  tname <- tnames[[afind]]
  basistype <- strsplit(tname,'[(]')[[1]][1]
  sstring <- paste(basistype,'[(]',af.term,'\\.omat,',af.term,'\\.tmat','[)]:L\\.',af.term,sep='')
  tind <- grepl(sstring,names(object$coef))

    if(!(length(tval)+length(xval))){
    temp <- list()
    temp[[paste('L.',af.term,sep='')]] <- 1
    if(plot.type=='persp' ){

      tvar <- tolower(af.term)

      if(object$fgam$ft[[af.ind]]$Qtransform){
        mtitle <- bquote(paste(hat(F),'(p,t),   p=',hat(G)[t],'(',.(tvar),')',sep=''))
      }else{
        mtitle <- bquote(paste(hat(F),'(',.(tvar),',t)',sep=''))
      }
      tvar <- paste('\n',tvar,sep='')
      vis.gam(object,view=c(paste(af.term,'.omat',sep=''),paste(af.term,'.tmat',sep='')),cond=temp,
            ticktype=ticktype,theta=theta,contour.col=rev(heat.colors(100)),
            xlab=tvar,ylab='\nt',zlab='',main=mtitle,...)

    }else if(plot.type=='contour'){
      # not making use of vis.gam because want colour key/legend
      nfine <- 51
      trange <- range(object$fgam$ft[[af.ind]]$xind)
      tvals <- seq(trange[1],trange[2],l=nfine)
      Xrange <- object$fgam$ft[[af.ind]]$Xrange
      xvals <- seq(Xrange[1],Xrange[2],l=nfine)
      newdata <- expand.grid(tvals,xvals)
      newot <- newdata[[1]]
      newX <- newdata[[2]]
      newdata <- list()
      newdata[[paste(af.term,'.omat',sep='')]] <- matrix(newX,ncol=1)
      newdata[[paste(af.term,'.tmat',sep='')]] <- matrix(newot,ncol=1)
      newdata[[paste('L.',af.term,sep='')]] <- matrix(1,ncol=1,nrow=length(newX))
      #attr(newdata[[paste(af.term,sep='')]],'L') <- matrix(1,nr=1,nc=length(newX))
      #attr(newdata[[paste(af.term,sep='')]],'tmat') <- matrix(newot,nr=1)

      varnames <- all.vars(terms(object))
      varnames <- varnames[!(varnames %in% c(paste('L.',af.term,sep=''),paste(af.term,'.tmat',sep=''),paste(af.term,'.omat',sep='')))]
      varnames <- varnames[-1]
      if(length(varnames)){
        for(i in 1:length(varnames)){
          newdata[[varnames[i]]] <- rep(object$fgam$datameans[names(object$fgam$datameans)==varnames[i]],l=length(newX))#matrix(object$fgam$datameans[names(object$fgam$datameans)==varnames[i]],nr=1,nc=length(newX))
        }
      }
      #newdata <- list2df(newdata)
      dmat <- predict.gam(object,newdata=newdata,type='lpmatrix',terms=af.term,na.action=na.exclude)[,tind]

      preds <- dmat%*%object$coef[tind]
      tvar <- tolower(af.term)
      if(object$fgam$ft[[af.ind]]$Qtransform){
        mtitle <- bquote(paste(hat(F),'(p,t),   p=',hat(G)[t],'(',.(tvar),')',sep=''))
      }else{
        mtitle <- bquote(paste(hat(F),'(',.(tvar),',t)',sep=''))
      }

      xlab=ifelse(object$fgam$ft[[af.ind]]$Qtransform,tvar,'p')
      lattice::levelplot(preds~newX*newot,contour=TRUE,labels=TRUE,pretty=TRUE,ylab='t',xlab=xlab,
                col.regions=rev(heat.colors(100)),main=as.expression(mtitle),...)

    }
  }else{
    if(length(xval)+length(tval)>1){
      warning('Specify only one value for either xval or tval.  Only first specified value will be used')
      if(length(xval)){
        xval=xval[1]
        tval=NULL
      }else{
        tval=tval[1]
      }
    }
    nfine <- 200
    parnames <- names(object$fgam$labelmap[sapply(object$fgam$labelmap,length)==0])
    parind=names(object$coef) %in% c('(Intercept)',parnames)
    ind=as.logical(parind+tind)
    afterm <- tolower(af.term)

    if(length(xval)){
      trange <- range(object$fgam$ft[[af.ind]]$xind)
      tvals <- seq(trange[1],trange[2],l=nfine)
      xvals <- rep(xval,l=nfine)
      xlab='t'
      x=tvals
      if(!deriv2){
        main=bquote(paste(hat(F)(.(round(xval,3)),t),' by t',sep=''))
        ylab=bquote(paste(hat(F)(.(round(xval,3)),t),sep=''))
      }else{
        main=bquote(paste(frac(partialdiff^2,partialdiff*.(afterm)^2)*hat('F')(.(afterm),t),
                          '|',phantom()[.(afterm)==.(round(xval,3))],' by t',sep=''))
        ylab=bquote(paste(frac(partialdiff^2,partialdiff*.(afterm)^2)*hat('F')(.(afterm),t),
                                       '|',phantom()[.(afterm)==.(round(xval,3))],sep=''))
      }
    }else{ #t fixed
      Xrange <- object$fgam$ft[[af.ind]]$Xrange
      tvals <- rep(tval,nfine)
      xvals <- seq(Xrange[1],Xrange[2],l=nfine)
      xlab='x'
      x=xvals
      if(!deriv2){
        main=bquote(paste(hat(F)(.(afterm),.(round(tval,3))),' by ',.(afterm),sep=''))
        ylab=bquote(paste(hat(F)(.(afterm),.(round(tval,3))),sep=''))
      }else{
        main=bquote(paste(frac(partialdiff^2,partialdiff*.(afterm)^2)*hat('F')(.(afterm),.(round(tval,3))),' by ',.(afterm),sep=''))
        ylab=bquote(paste(frac(partialdiff^2,partialdiff*.(afterm)^2)*hat('F')(.(afterm),.(round(tval,3))),sep=''))
      }

    }
    newdata <- list()
    newdata[[paste(af.term,'.omat',sep='')]] <- xvals
    newdata[[paste(af.term,'.tmat',sep='')]] <- tvals
    newdata[[paste('L.',af.term,sep='')]] <- rep(1,l=length(xvals))
    varnames <- all.vars(terms(object))
    varnames <- varnames[!(varnames %in% c(paste('L.',af.term,sep=''),paste(af.term,'.tmat',sep=''),paste(af.term,'.omat',sep='')))]
    varnames <- varnames[-1]
    if(length(varnames)){
      for(i in 1:length(varnames)){
        newdata[[varnames[i]]] <- rep(object$fgam$datameans[names(object$fgam$datameans)==varnames[i]],l=length(xvals))#matrix(object$fgam$datameans[names(object$fgam$datameans)==varnames[i]],nr=1,nc=length(newX))
      }
    }

    lpmat <- predict.gam(object,newdata=newdata,type='lpmatrix')[,ind]
    if(deriv2){
      eps <- 1e-7 ## finite difference interval
      newdata[[paste(af.term,'.omat',sep='')]] <- xvals + eps
      X1 <- predict.gam(object,newdata=newdata,type='lpmatrix')[,ind]

      newdata[[paste(af.term,'.omat',sep='')]] <- xvals + 2*eps
      X2 <- predict.gam(object,newdata=newdata,type='lpmatrix')[,ind]

      lpmat <- (X2-2*X1+lpmat)/eps^2
    }

    preds <- sum(object$cmX)+lpmat%*%object$coef[ind]
    se.preds <- rowSums(lpmat%*%object$Vp[ind,ind]*lpmat)^.5
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
