wcr <- function(y, xfuncs, min.scale, nfeatures, ncomp, method = c("pcr", "pls"), 
                mean.signal.term = FALSE, covt = NULL, filter.number = 10, 
                wavelet.family = "DaubLeAsymm", family = "gaussian", cv1 = FALSE, 
                nfold = 5, nsplit = 1, store.cv = FALSE, store.glm = FALSE, seed = NULL) {
    checkError(y = y, xfuncs = xfuncs, covt = covt, 
               vars = list(min.scale = min.scale, nfeatures = nfeatures, ncomp = ncomp))    
    method <- match.arg(method)
    n <- length(y)
    dim.sig <- length(dim(xfuncs)) - 1
    d <- dim(xfuncs)[2]
    multiValFlag <- c(length(min.scale) > 1, length(nfeatures) > 1, length(ncomp) > 1)   
    names(multiValFlag) <- c("min.scale", "nfeatures", "ncomp")
    if (any(multiValFlag) || cv1) {
        do.cv <- TRUE
        if (!is.null(seed)) set.seed(seed)
        splits <- replicate(nsplit, split(sample(1:n), rep(1:nfold, length=n)))
    } else {
        do.cv <- FALSE
        store.cv <- FALSE
        nfold <- 1
        nsplit <- 1
  	}   
  	inm <- ifelse(length(which(multiValFlag == TRUE)) == 0, ifelse(nfold > 1, 'nfold', ''), 
                  switch(max(which(multiValFlag == TRUE)), 
                         ifelse(nfold > 1, 'nfold', 'ms'), 'nf', 'ncomp', '')) 	
    p <- prod(dim(xfuncs)[-1]) + ifelse(dim.sig == 2, 1, 0)
    n.unpen.cols <- 1 + mean.signal.term + ifelse(is.null(covt), 0, ncol(as.matrix(covt)))
    fhat.eigen <- array(0, dim = c(max(ncomp), d^dim.sig))  
    cv.table <- array(0, dim = c(nsplit, nfold, length(min.scale), length(nfeatures), length(ncomp)))
    dimnames(cv.table) <- list(NULL, NULL, paste("ms =", min.scale), paste("nfeatures =", nfeatures), 
                               paste("ncomp =", ncomp))
  	setup <- waveletSetup(xfuncs = xfuncs, filter.number = filter.number, family = wavelet.family, 
                          min.scale = min.scale)
    if(do.cv && any(multiValFlag)){
        if (!multiValFlag["min.scale"]) cat("min.scale =", min.scale, " ")
        if (!multiValFlag["nfeatures"]) cat("nfeatures =", nfeatures, " ")
        if (!multiValFlag["ncomp"]) cat("ncomp =", ncomp, " ")
        cat("\n")
    }
    dim(xfuncs) <- c(n, d^dim.sig)
    for (isplit in 1 : nsplit){
    	if (do.cv){
    	    if (nsplit != 1) cat("split", isplit, "\n")
            if (nfold > 1) groups <- splits[,isplit] else groups <- splits[isplit]
    	}
    	if (inm == 'ms') cat("min.scale = ")
	    for (ims in 1 : length(min.scale)){
	    	if (multiValFlag["min.scale"]) if (inm == 'ms') cat(min.scale[ims], " ")
	    	setup$temp$callInfo$min.scale <- min.scale[ims]
	    	if (inm == 'nfold'){
	    		if (multiValFlag["min.scale"]) cat("min.scale =", min.scale[ims], " ")
	    		cat("fold ")
	    	}
	        for (ifold in 1 : nfold) {
	        	if (nfold > 1) if (inm == 'nfold') cat(ifold, " ")
	            if (do.cv) {
	                idxTest <- groups[[ifold]]
	                idxTrain <- (1:n)[-idxTest]
	            } else {
	            	idxTest <- idxTrain <- 1:n
	            }
	            if (method == "pcr") criteria <- apply(setup$coef[[ims]][idxTrain, ], 2, var)
	            else criteria <- abs(as.vector(cov(setup$coef[[ims]][idxTrain, ], y[idxTrain])))
	            names(criteria) <- 1:p
	            sorted <- sort(criteria, decreasing = TRUE, na.last = TRUE)[1:max(nfeatures)]
	            subset <- setup$coef[[ims]][, as.numeric(names(sorted))]
	            if (inm == 'nf'){
	            	if (multiValFlag["min.scale"]) cat("min.scale =", min.scale[ims], " ")
	            	if (nfold > 1) cat("fold", ifold, " ")
	            	cat("nfeatures = ")
	            }
	            for (infeatures in 1 : length(nfeatures)) {
	            	if (multiValFlag["nfeatures"]) if (inm == 'nf') cat(nfeatures[infeatures], " ")
	                X0 <- if (mean.signal.term) as.matrix(rowMeans(subset[, 1:nfeatures[infeatures]])) else NULL
	                if (!is.null(covt)) X0 <- cbind(X0, as.matrix(covt))
	                sigs.decor <- if (is.null(X0)) scale(subset[idxTrain, 1:nfeatures[infeatures]], scale = FALSE)
	                              else lm(subset[idxTrain, 1:nfeatures[infeatures]] ~ X0[idxTrain,] - 1)$resid
	                if (method == "pcr") {
	                	V <- svd(sigs.decor)$v
	                } else {
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
	                    setup$temp$coef <- rep(0, p)
	                    setup$temp$coef[as.numeric(colnames(subset[, 1:nfeatures[infeatures]]))] <- V[,i]
	                    fhat.eigen[i,] <- matrix(setup$rec(setup$temp), nrow = 1)
	                }
	                if (multiValFlag["ncomp"]){
	                	if (multiValFlag["min.scale"]) cat("min.scale =", min.scale[ims], " ")
	                	if (nfold > 1) cat("fold", ifold, " ")
	                	if (multiValFlag["nfeatures"]) cat("nfeatures =", nfeatures[infeatures], " ")
	                	cat("ncomp = ")
	                }
	                for (icomp in 1 : length(ncomp)) {
	                    if (ncomp[icomp] <= nfeatures[infeatures]) {
	                    	if (multiValFlag["ncomp"]) cat(ncomp[icomp], " ")
	                        X <- cbind(X0[idxTrain,], sigs.decor %*% V[,1:ncomp[icomp]])
	                        obje <- glm(y[idxTrain] ~ X, family = get(family))                      
	                        fhat <- t(matrix(fhat.eigen[1:ncomp[icomp], ], ncol = d^dim.sig)) %*% 
	                                obje$coef[-(1:n.unpen.cols)]
	                        undecor.coef <- obje$coef[1:n.unpen.cols] - 
	                                        ginv(cbind(rep(1, length(idxTrain)), X0[idxTrain,])) %*% 
	                                        xfuncs[idxTrain, ] %*% fhat
	                        X0.tst <- cbind(matrix(1,n,1), X0)[idxTest, ]
	                        yhat <- X0.tst %*% undecor.coef + xfuncs[idxTest, ] %*% fhat
	                        if (family == "gaussian"){
	                            cv.table[isplit, ifold, ims, infeatures, icomp] <- mean((yhat - y[idxTest]) ^ 2)                                                                                                             
	                        } else if (family == "binomial") {
	                        	# negative log likelihood
	                            phat <- exp(yhat) / (1 + exp(yhat))
	                            phat <- replace(phat, exp(yhat) == Inf, 1)
	                            cv.table[isplit, ifold, ims, infeatures, icomp] <- 
	                                 - sum(log(phat)[y[idxTest] == 1]) - sum(log((1-phat))[y[idxTest] == 0])
	                        }
	                    }
	                } # ncomp
	                if (inm == 'ncomp') cat('\n')
	            }   # nfeatures loop
	            if (inm == 'nf') cat('\n')
	        }       # fold loop
	        if (inm == 'nfold') cat('\n')
	    }           # min.scale loop
	    if (inm == 'ms') cat('\n')
    }               # nsplit loop
    if (do.cv) {
    	res <- waveletGetCV(cv.table = cv.table, family = family)
    	param <- array(0, dim = c(2, 3))
    	param[1, ] <- res$idxmin
    	param[2, ] <- c(min.scale[param[1,1]], nfeatures[param[1,2]], ncomp[param[1,3]])
    	dimnames(param) = list(c("index", "value"), c("min.scale", "nfeatures", "ncomp"))
        dim(xfuncs) <- c(n, rep(d, dim.sig))
        l <- wcr(y = y, xfuncs = xfuncs, min.scale = param[2,1], nfeatures = param[2,2], 
                    ncomp = param[2,3], method=method, mean.signal.term=mean.signal.term, 
                    covt = covt, filter.number=filter.number, wavelet.family=wavelet.family,
                    family=family, cv1 = FALSE, store.glm = store.glm)
        if (store.cv) {
        	l$cv.table <- res$cv.table
        	l$se.cv <- res$se.cv
        } else {
        	l$cv.table <- min(res$cv.table[res$cv.table !=0], na.rm=TRUE)
        }
        l$tuning.params <- param
    } else {
    	l <- list()       
    	if (store.glm){
    		l$glm <- obje
    	} 
    	l$fitted.values <- obje$fitted.values
    	l$param.coef <- obje$coef[1:n.unpen.cols]
    	l$undecor.coef <- matrix(obje$coef[1:n.unpen.cols] - ginv(cbind(rep(1, n), X0)) %*% xfuncs %*% fhat, nrow=1)
        l$fhat <- array(fhat, dim = rep(d, dim.sig))	
        l$Rsq <- getRsq(y = y, yhat = obje$fitted.values, family = family)
    }
    l$family <- family
    class(l) <- 'wcr'
    return(l)
}