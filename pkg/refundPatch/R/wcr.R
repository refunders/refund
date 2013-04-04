wcr <- function(y, xfuncs, min.scale, nfeatures, ncomp, method = c("pcr", "pls"), 
                mean.signal.term = FALSE, covt = NULL, filter.number = 10, 
                wavelet.family = "DaubLeAsymm", family = "gaussian", 
                cv1 = FALSE, nfold = 5, compare.fits = FALSE, store.cv = FALSE, 
                store.glm = TRUE) {
    method <- match.arg(method)
    # Error handling
    n <- length(y)
    if (!is.array(xfuncs) || !length(dim(xfuncs)) %in% 2:3)
        stop("Argument 'xfuncs' is invalid: must be either a 2D or 3D array.")
    dim.sig <- length(dim(xfuncs)) - 1
    if (dim(xfuncs)[1] != n)
        stop("Arguments y and xfuncs have invalid lengths: ", length(y), " and ", dim(xfuncs)[1], ".")
    if (dim(xfuncs)[2] != dim(xfuncs)[1 + dim.sig])
        stop("Number of rows and columns in image must be identical.")
    d <- dim(xfuncs)[2]
    if (as.integer(log2(d)) != log2(d))
        stop("Argument 'xfuncs' is invalid: length must be a power of 2.")
    if (sum(!min.scale %in% 0:(log2(d) - 1)) != 0)
        stop("Argument 'min.scale' is invalid: must be integer(s) between 0 and ", log2(d) - 1, ".")
    if (sum(!nfeatures %in% 1:(d^dim.sig)) != 0)
        stop("Argument 'nfeatures' is invalid: must be integer(s) between 1 and ", d^dim.sig, ".")
    if (sum(!ncomp %in% 1:max(nfeatures)) != 0)
        stop("Argument 'ncomp' is invalid: must be integers() between 1 and ", max(nfeatures), ".")
    if (!family %in% c("gaussian", "binomial")) 
        stop("Only 'gaussian' and 'binomial' models are implemented in the current version.")
    if (!nfold %in% 1:n) 
        stop("Argument 'nfold' is invalid: must be an integer between 1 and ", n, ".") 
    # Determine CV or not 
    if (length(min.scale) + length(nfeatures) + length(ncomp) > 3 || cv1) {
        do.cv <- TRUE
        # Set up CV
        groups <- split(sample(1:n), rep(1:nfold, length=n))
    }
    else {
        do.cv <- FALSE
        store.cv <- FALSE
        compare.fits <- FALSE
        nfold <- 1
  	}
    
    # Wavelet decomposition
    if (dim.sig == 1){
        wave.decomp <- wd
        dec <- decomp
        rec <- reconstr
    } 
    else {
        wave.decomp <- imwd
        dec <- decomp2d
        rec <- reconstr2d	
    }
    wdobjs <- apply(xfuncs, 1, wave.decomp, filter.number = filter.number, family = wavelet.family)
    temp <- dec(wdobjs[[1]]) 
    p <- length(temp$coef)
    n.unpen.cols <- 1 + mean.signal.term + ifelse(is.null(covt), 0, ncol(as.matrix(covt)))
    fhat.eigen <- array(0, dim = c(max(ncomp), d^dim.sig))
    if (dim.sig == 2) dim(xfuncs) <- c(n, d^dim.sig)
  
    # Begin CV
    cv.table <- array(0, dim = c(length(min.scale), length(nfeatures), length(ncomp)))
    dimnames(cv.table) <- list(paste("ms =", min.scale), paste("nfeatures =", nfeatures), 
                               paste("ncomp =", ncomp))
    if (compare.fits) {
    	fhat.table <- array(0, dim = c(d^dim.sig, nfold, length(min.scale), length(nfeatures), length(ncomp)))
    	dimnames(fhat.table) <- list(NULL, paste("nfold =", 1 : nfold), paste("ms =", min.scale), 
    	                             paste("nfeatures =", nfeatures), paste("ncomp =", ncomp))
    }
    for (ims in 1 : length(min.scale)){
        coef <- t(array(unlist(sapply(wdobjs, dec, min.scale = min.scale[ims])[1,]), dim = c(p, n)))
        dimnames(coef) <- list(NULL, 1:p)  
        for (ifold in 1 : nfold) {
            if (do.cv) {
                idxTest <- groups[[ifold]]
                idxTrain <- (1:n)[-idxTest]
            }
            else idxTest <- idxTrain <- 1:n
            if (method == "pcr") criteria <- apply(coef[idxTrain, ], 2, var)
            else criteria <- abs(as.vector(cov(coef[idxTrain,], y[idxTrain])))
            names(criteria) <- 1:p
            sorted <- sort(criteria, decreasing = TRUE, na.last = TRUE)[1:max(nfeatures)]
            subset <- coef[, as.numeric(names(sorted))]
            temp$callInfo$min.scale <- min.scale[ims]
            for (infeatures in 1 : length(nfeatures)) {
                cat("min.scale:", min.scale[ims], "\tfold:", ifold, "\tnfeatures:",nfeatures[infeatures], "\t")
                X0 <- if (mean.signal.term) as.matrix(rowMeans(subset[, 1:nfeatures[infeatures]])) else NULL
                if (!is.null(covt)) X0 <- cbind(X0, as.matrix(covt))
                sigs.decor <- if (is.null(X0)) scale(subset[idxTrain, 1:nfeatures[infeatures]], scale = FALSE)
                              else lm(subset[idxTrain, 1:nfeatures[infeatures]] ~ X0[idxTrain,] - 1)$resid
                if (method == "pcr") V <- svd(sigs.decor)$v
                else {
                    e <- sigs.decor
                    f <- y[idxTrain]
                    W <- P <- array(0, dim = c(nfeatures[infeatures], min(max(ncomp), nfeatures[infeatures])))
                    for (i in 1 : min(max(ncomp), nfeatures[infeatures])) {
                        svdEF <- svd(crossprod(e, f), nu = 1, nv = 1)
                        W[, i] <- svdEF$u
                        scoret <- e %*% W[, i]
                        normt <- scoret / drop(sqrt(crossprod(scoret)))
                        P[, i] <- t(e) %*% normt
                        e <- e - tcrossprod(normt) %*% e
                        f <- f - tcrossprod(normt) %*% f
                    }
                    V <- W %*% ginv(t(P) %*% W)
                }
                for (i in 1 : min(max(ncomp), nfeatures[infeatures])){
                    temp$coef <- rep(0, p)
                    temp$coef[as.numeric(colnames(subset[, 1:nfeatures[infeatures]]))] <- V[,i]
                    fhat.eigen[i,] <- matrix(rec(temp), nrow = 1)
                }
                cat("Number of components: ")
                for (icomp in 1 : length(ncomp)) {
                    if (ncomp[icomp] <= nfeatures[infeatures]) {
                        cat(ncomp[icomp], "\t")
                        X <- cbind(X0[idxTrain,], sigs.decor %*% V[,1:ncomp[icomp]])
                        obje <- glm(y[idxTrain] ~ X, family = get(family))                      
                        fhat <- t(matrix(fhat.eigen[1:ncomp[icomp], ], ncol = d^dim.sig)) %*% 
                                obje$coef[-(1:n.unpen.cols)]
                        if (compare.fits) fhat.table[,ifold, ims, infeatures, icomp] <- fhat
                        undecor.coef <- obje$coef[1:n.unpen.cols] - 
                                        ginv(cbind(rep(1, length(idxTrain)), X0[idxTrain,])) %*% 
                                        xfuncs[idxTrain, ] %*% fhat
                        X0.tst <- cbind(matrix(1,n,1), X0)[idxTest, ]
                        yhat <- X0.tst %*% undecor.coef + xfuncs[idxTest, ] %*% fhat
                        if (family == "gaussian")
                            cv.table[ims, infeatures, icomp] <- cv.table[ims, infeatures, icomp] + 
                                                                mean((yhat - y[idxTest]) ^ 2)                                                                                                             
                        else if (family == "binomial") {
                            phat <- exp(yhat) / (1 + exp(yhat))
                            phat <- replace(phat, exp(yhat) == Inf, 1)
                            cv.table[ims, infeatures, icomp] <- cv.table[ims, infeatures, icomp] + 
                                                                mean((phat > mean(y[idxTrain])) != y[idxTest])
                        }
                    }
                }
                cat("\n") 
            }   # nfeatures loop
        }       # fold loop
    }           # min.scale loop
    if (do.cv) {
        idxmin <- which(cv.table == min(cv.table[cv.table != 0], na.rm = TRUE), arr.ind = TRUE)
        if (nrow(idxmin) > 1) idxmin <- idxmin[1,]
        dim(xfuncs) <- c(n, rep(d, dim.sig))
        min.scale <- min.scale[idxmin[1]]
        nfeatures <- nfeatures[idxmin[2]]
        ncomp <- ncomp[idxmin[3]]
        obje <- wcr(y=y, xfuncs=xfuncs, min.scale=min.scale, nfeatures=nfeatures, ncomp=ncomp, 
                    method=method, mean.signal.term=mean.signal.term, covt = covt,
                    filter.number=filter.number, wavelet.family=wavelet.family,
                    family=family, cv1 = FALSE, store.glm = store.glm)
        if (store.cv) obje$cv.table <- cv.table
        else obje$cv.table <- min(cv.table[cv.table!=0], na.rm = TRUE)
        if (compare.fits) {
            obje$stability <- apply(apply(fhat.table, 3:5, apply, 1, scale, center=TRUE, scale=FALSE)^2, 2:4, sum)/(d^dim.sig*(nfold-1))
        }
    } 
    else {
        undecor.coef <- obje$coef[1:n.unpen.cols] - ginv(cbind(rep(1, n), X0)) %*% xfuncs %*% fhat
        if (!store.glm) {
    		yhat <- obje$fitted.values
    		obje <- list()
    		obje$fitted.values <- yhat
    	}
    	obje$undecor.coef <- matrix(undecor.coef, nrow = 1)
    	colnames(obje$undecor.coef) <- if (is.null(dimnames(covt)) || is.null(dimnames(covt)[[2]])) 
                                           paste("X", 0:(n.unpen.cols-1), sep="")
    	                               else c("Intercept", dimnames(covt)[[2]])
        if (dim.sig == 2) dim(fhat) <- c(d, d)
        obje$fhat <- fhat	
        obje$min.scale <- min.scale
        obje$nfeatures <- nfeatures
        obje$ncomp <- ncomp
    }
    return(obje)
}

