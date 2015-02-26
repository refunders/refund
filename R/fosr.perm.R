fosr.perm <-
function(Y=NULL, fdobj=NULL, X, con=NULL, X0=NULL, con0=NULL, argvals = NULL,   
lambda=NULL, lambda0=NULL, multi.sp=FALSE, nperm, level=.05, plot=TRUE, xlabel="", title=NULL, prelim=15, ...) {
    fpobj1 = fosr.perm.fit(Y=Y, fdobj=fdobj, X=X, con=con, X0=X0, con0=con0, argvals=argvals, lambda=lambda, lambda0=lambda0, multi.sp=multi.sp, nperm=nperm, prelim=prelim, ...)
    fpobj2 = fosr.perm.test(fpobj1, level=level)
    if (plot) plot(fpobj2, xlabel=xlabel, title=title)
    fpobj2
}

