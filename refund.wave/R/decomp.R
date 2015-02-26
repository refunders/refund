decomp <- function(wdobj, min.scale = 0) {
    if (class(wdobj) != "wd") {
        stop("wdobj must be of class 'wd'")
    }    
    rowNum <- 2 ^ wdobj$nlevels
    level <- min.scale + 1
    first.last.d <- wdobj$fl.dbase$first.last.d
    first.level <- first.last.d[level, 1]
    last.level <- first.last.d[level, 2]
    offset.level <- first.last.d[level, 3]
    n <- last.level - first.level + 1
    coef <- c(wdobj$D[seq(offset.level + n)], accessC(wdobj, level = min.scale, boundary = T))
    ncoef <- length(coef)
    l <- list(coef = coef, rowNum = rowNum, ncoef = ncoef, min.scale = min.scale, callInfo = list(filter = wdobj$filter, type = wdobj$type, bc = wdobj$bc, index = list(offset.level = offset.level, n = n)), date = wdobj$date) 
    class(l) <- "decomp"   
    return(l)
}
