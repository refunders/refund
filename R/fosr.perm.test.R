#' @export
#' @rdname fosr.perm
fosr.perm.test <-
function(x, level=.05) {
	if (!inherits(x, "fosr.perm")) stop("First argument must be an object of class 'fosr.perm")
    x$level = level
    critval = quantile(apply(x$F.perm, 1, max), 1-level, type=1)
    x$critval=critval
    if (length(level)==1) {
        x$signif = (x$F > critval)
        x$n2s = which(diff(c(0, x$signif))==1)
        x$s2n = which(diff(c(x$signif, 0))==-1)
    }
    x
}

