irreg2mat <- function(ydata, binning=FALSE, maxbins=1e3){
    # TODO: arg checks
    
    ## drop any missings:
    ydata <- ydata[complete.cases(ydata), ]
    
    ## turn into row/column indices for new matrix
    nobs <- length(unique(ydata$.id))
    # make sure newid takes values 1:nobs
    newid <- as.numeric(as.factor(ydata$.id))
    
    ## bin y-index, if necessary
    bins <- sort(unique(ydata$.index))
    if(binning && (length(bins) > maxbins)){
        # linear binning;
        # bin-borders go from just below min to just above max
        # TODO: quantile-based binning?
        binvalues <- seq((1-.001*sign(bins[1]))*bins[1], 
                         (1+.001*sign(bins[length(bins)]))*bins[length(bins)],
                         l=maxbins+1)
        bins <- binvalues
        binvalues <- head(filter(binvalues, c(.5, .5)), -1)
     } else {
        binvalues <- bins
        bins <- c((1-.001*sign(bins[1]))*bins[1], 
                  bins[-length(bins)],
                  (1+.001*sign(bins[length(bins)]))*bins[length(bins)])
        # take care of edge cases:
        if(bins[1] == 0) bins[1] <- -.001
        if(bins[length(bins)] == 0) bins[length(bins)] <- .001
    }
    newindex <- cut(ydata$.index, breaks=bins, include.lowest=TRUE)
    
    Y <- matrix(NA, nrow=nobs, ncol=nlevels(newindex))
    colnames(Y) <- binvalues
    attr(Y, "index") <- binvalues
    Y[cbind(newid, as.numeric(newindex))] <- ydata$.value
    return(Y)
}