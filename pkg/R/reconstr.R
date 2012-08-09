reconstr <- function(decompobj) {  
	ctmp <- class(decompobj)
	if (is.null(ctmp)) {
	    stop("decompobj has no class")
	}
	if (ctmp != "decomp") {
	    stop("decompobj is not of class decompobj")
	}
	coeff <- decompobj$coef
	rowNum <- decompobj$rowNum
	min.scale <- decompobj$min.scale
	bc <- decompobj$callInfo$bc
	filter.number <- decompobj$callInfo$filter$filter.number
	family <- decompobj$callInfo$filter$family
	type <- decompobj$callInfo$type
	index <- decompobj$callInfo$index
    wds <- wd(rep(0, rowNum), filter.number = filter.number, family = family, type = type, bc = bc, min.scale = min.scale)
    offset.level <- index$offset.level
    n <- index$n  
    beta <- NULL
    wds$D[seq(offset.level + n)] <- coeff[seq(offset.level + n)]
    wds <- putC(wds, level = min.scale, boundary = TRUE, coeff[-seq(offset.level + n)])
    beta <- wr(wds, start.level = min.scale)
    return(beta)
}
