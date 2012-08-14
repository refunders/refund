reconstr2d <- function(decomp2dobj){
	ctmp <- class(decomp2dobj)
	if (is.null(ctmp)){
	    stop('decomp2dobj has no class')
	} else if (ctmp != 'decomp2d'){
	    stop('decomp2dobj is not of class decomp2d')
	}
	require("wavethresh")
    rowNum <- decomp2dobj$rowNum
    min.scale <- decomp2dobj$callInfo$min.scale
    filter.number <- decomp2dobj$callInfo$filter$filter.number
    family <- decomp2dobj$callInfo$filter$family
    type <- decomp2dobj$callInfo$type
    bc <- decomp2dobj$callInfo$bc
    wdobj <- imwd_test(array(0, dim = c(rowNum, rowNum)), filter.number = filter.number, family = family, type = type, bc = bc)

    # coefficient D
    # listIdx - iterator through the coefficients list of w
    listIdx <- 1
    # arrIdx - iterator through the columns of coeffArr
    arrIdx <- 1
    # replace the element in the 'wavelet' coefficients list with the
    # corresponding coefficients in the coeffArr
    while(listIdx <= (wdobj$nlevels - min.scale) * 4){
       if (listIdx %% 4 != 1){
           wdobj[[6 + listIdx]] <-
                decomp2dobj$coef[arrIdx : (arrIdx - 1 + length(wdobj[[6 + listIdx]]))]
           arrIdx <- arrIdx + length(wdobj[[6 + listIdx]])
       }
       listIdx <- listIdx + 1
    }
    # coefficient C
    wdobj[[7 + (wdobj$nlevels - min.scale - 1) * 4]] <-
    decomp2dobj$coef[arrIdx : (arrIdx - 1 + length(wdobj[[7 + (wdobj$nlevels-min.scale-1)*4]]))]
    arrIdx <- arrIdx + length(wdobj[[7 + (wdobj$nlevels-min.scale-1)*4]])
    # w0Lconstant when min.scale = 0
    if (min.scale == 0){
        wdobj[[6+listIdx]] <- decomp2dobj$coef[arrIdx]
    } else{
        smooth <- matrix(wdobj[[7+(wdobj$nlevels-min.scale-1)*4]],nrow=2^ min.scale)
        smoothW <- imwd_test(smooth, family = family, filter.number = filter.number,
                             type = type, bc = bc)
        tail <- length(wdobj)
        for (smoothIdx in length(smoothW) : 7){
            wdobj[[tail]] <- smoothW[[smoothIdx]]
            tail <- tail - 1
        }
    }
    arr <- tryCatch(imwr(wdobj), error = function(e) print(e$message))
    return (arr)
}
