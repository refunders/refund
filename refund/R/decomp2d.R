decomp2d <- function(imwdobj, min.scale = 0){
    if(class(imwdobj) != "imwd"){
        stop("imwdobj must be of class 'imwd'")
    }    
    if (min.scale >= imwdobj$nlevels){
        stop(cat("exceed the maximum decomposition level: ", imwdobj$nlevels - 1))
    }
    rowNum <- 2 ^ imwdobj$nlevels

    coef <- rep(0, rowNum^2 + 1)

    listIdx <- arrIdx <- 1

    while(listIdx <= (log2(rowNum) - min.scale) * 4){
        if (listIdx %% 4 != 1){
            lD <- length(imwdobj[[listIdx + 6]])
            coef[arrIdx : (arrIdx + lD - 1)] <- imwdobj[[listIdx + 6]]
            arrIdx <- arrIdx + lD
        }
        listIdx <- listIdx + 1
    }
    lC <- length(imwdobj[[7 + (log2(rowNum) - min.scale - 1) * 4]])
    coef[arrIdx : (arrIdx +lC - 1)] <- imwdobj[[7 + (log2(rowNum) - min.scale - 1) * 4]]

    if (min.scale == 0){
         coef[length(coef)] <- imwdobj[[listIdx + 6]]
    }
    l <- list(coef = coef, rowNum = rowNum, callInfo = list(min.scale = min.scale, 
              filter = imwdobj$filter, type = imwdobj$type, bc = imwdobj$bc))
    class(l) <- 'decomp2d'
    return (l)
}
