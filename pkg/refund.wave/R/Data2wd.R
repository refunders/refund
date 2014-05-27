Data2wd <-
function(y, xfuncs, covt = NULL, min.scale = 0, nfeatures = NULL, filter.number = 10,
                    wavelet.family = 'DaubLeAsymm') {
    if (!is.null(covt)) covt <- as.matrix(covt) 
    n <- length(y)
    dim.sig <- length(dim(xfuncs)) - 1
    d <- dim(xfuncs)[2]
    p <- prod(dim(xfuncs)[-1]) + ifelse(dim.sig == 2, 1, 0)
    if (is.null(nfeatures)) nfeatures <- p
    setup <- waveletSetup(xfuncs = xfuncs, filter.number = filter.number, 
                          family = wavelet.family, min.scale = min.scale)
    criteria <- apply(setup$coef[[1]], 2, var)
    names(criteria) <- 1:p
	sorted <- sort(criteria, decreasing = TRUE, na.last = TRUE)[1:max(nfeatures)]
	subset <- setup$coef[[1]][, as.numeric(names(sorted))]                      
    coef.red <- subset[,1:nfeatures]
    X <- cbind(covt, coef.red)
    setup$ncovt = ifelse(is.null(covt), 0, ncol(covt))
    setup$dim.sig = dim.sig
    setup$d = dim(xfuncs)[2]
    return(list(X = X, info = setup))	
}
