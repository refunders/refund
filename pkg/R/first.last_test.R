first.last_test <-
function (LengthH, DataLength, type = "wavelet", bc = "periodic",
    current.scale = 0)
{
    if (type == "station" && bc != "periodic")
        stop("Can only do periodic boundary conditions with station")
    if (type != "station" && type != "wavelet")
        stop("Type can only be wavelet or station")
    levels <- log(DataLength)/log(2)
    first.last.c <- matrix(0, nrow = levels + 1, ncol = 3, dimnames =
                           list(NULL, c("First", "Last", "Offset")))
    first.last.d <- matrix(0, nrow = levels - current.scale,
        ncol = 3, dimnames = list(NULL, c("First", "Last", "Offset")))
    if (bc == "periodic") {
        if (type == "wavelet") {
            first.last.c[, 1] <- rep(0, levels + 1)
            first.last.c[, 2] <- 2^(0:levels) - 1
            first.last.c[, 3] <- rev(c(0, cumsum(rev(1 + first.last.c[,
                2]))[1:levels]))
            first.last.d[, 1] <- rep(0, levels)
            first.last.d[, 2] <- 2^(0:(levels - 1)) - 1

          # print(first.last.d[1,])

          # print(length(first.last.d[, 3]))
          # print(length(rev(c(0, cumsum(rev(1 + first.last.d[, 2]))[1:(levels - 1)]))))


          ######################################################################
          # When the input is a 2*2 image, levels == 1, boundary case error in
          # [1:(levels - 1)]
          # Original code:
          # first.last.d[, 3] <- rev(c(0, cumsum(rev(1 + first.last.d[,
          #      2]))[1:(levels - 1)]))
          #
          # Code modified from here:
            if (levels == 1){
                first.last.d[, 3] <- 0
            } else{
            first.last.d[, 3] <- rev(c(0, cumsum(rev(1 + first.last.d[,
                2]))[1:(levels - 1)]))
            }
          ######################################################################

            ntotal <- 2 * DataLength - 1
            ntotal.d <- DataLength - 1
        }
        else if (type == "station") {
            first.last.c[, 1] <- rep(0, levels + 1)
            first.last.c[, 2] <- 2^levels - 1
            first.last.c[, 3] <- rev(c(0, cumsum(rev(1 + first.last.c[,
                2]))[1:levels]))
            first.last.d[, 1] <- rep(0, levels)
            first.last.d[, 2] <- 2^levels - 1
            first.last.d[, 3] <- rev(c(0, cumsum(rev(1 + first.last.d[,
                2]))[1:(levels - 1)]))
            ntotal <- (levels + 1) * 2^levels
            ntotal.d <- levels * 2^levels
        }
    }
    else if (bc == "symmetric") {
        first.last.c[levels + 1, 1] <- 0
        first.last.c[levels + 1, 2] <- DataLength - 1
        first.last.c[levels + 1, 3] <- 0
        ntotal <- first.last.c[levels + 1, 2] - first.last.c[levels +
            1, 1] + 1
        ntotal.d <- 0
        for (i in levels:1) {
            first.last.c[i, 1] <- trunc(0.5 * (1 - LengthH +
                first.last.c[i + 1, 1]))
            first.last.c[i, 2] <- trunc(0.5 * first.last.c[i +
                1, 2])
            first.last.c[i, 3] <- first.last.c[i + 1, 3] + first.last.c[i +
                1, 2] - first.last.c[i + 1, 1] + 1
            first.last.d[i, 1] <- trunc(0.5 * (first.last.c[i +
                1, 1] - 1))
            first.last.d[i, 2] <- trunc(0.5 * (first.last.c[i +
                1, 2] + LengthH - 2))
            if (i != levels) {
                first.last.d[i, 3] <- first.last.d[i + 1, 3] +
                  first.last.d[i + 1, 2] - first.last.d[i + 1,
                  1] + 1
            }
            ntotal <- ntotal + first.last.c[i, 2] - first.last.c[i,
                1] + 1
            ntotal.d <- ntotal.d + first.last.d[i, 2] - first.last.d[i,
                1] + 1
        }
    }
    else if (bc == "interval") {
        first.last.d[, 1] <- rep(0, levels - current.scale)
        first.last.d[, 3] <- 2^(current.scale:(levels - 1))
        first.last.d[, 2] <- first.last.d[, 3] - 1
        first.last.c <- c(0, 2^current.scale - 1, 0)
        return(list(first.last.c = first.last.c, first.last.d = first.last.d))
    }
    else {
        stop("Unknown boundary correction method")
    }
    names(ntotal) <- NULL
    names(ntotal.d) <- NULL
    list(first.last.c = first.last.c, ntotal = ntotal, first.last.d = first.last.d,
         ntotal.d = ntotal.d)
}
