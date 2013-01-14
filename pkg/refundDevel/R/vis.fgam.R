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
  afind <- grep(paste('te[(]',af.term,sep = ""), tnames) #index of desired af among all predictors
  af.ind <- grep(paste('af[(]',af.term,sep = ""), names(object$fgam$ft)) #index of desired af among all func. predictors
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
      levelplot(preds~newX*newot,contour=TRUE,labels=TRUE,pretty=TRUE,ylab='t',xlab=xlab,
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