.onAttach <- function(libname, pkgname) {

    ## get package version
    vers <- packageDescription("FDboost")[["Version"]]

    packageStartupMessage("This is FDboost ", vers, ". ", 
                          appendLF = TRUE)
    return(TRUE)
}

.onLoad <- function(libname, pkgname) {
    options("mboost_indexmin" = +Inf) ### in FDboost, do not use ties in the data  
}

.onUnload <- function(libpath) {
    if (is.infinite(options("mboost_indexmin")[[1]]) )
      options("mboost_indexmin" = 10000)  ### use mboost-default to use ties again
}
