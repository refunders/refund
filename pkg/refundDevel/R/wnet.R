wnet <- function(y, xfuncs, covt = NULL, min.scale, nfeatures, alpha, lambda = NULL, 
                 standardize = FALSE, pen.covt = FALSE, filter.number = 10, 
                 wavelet.family = 'DaubLeAsymm', family = 'gaussian', nfold = 5, 
                 nsplit = 1, compare.fits = FALSE, store.cv = FALSE, seed = NULL, 
                 ...){
    n <- length(y)
    if (!is.array(xfuncs) || !length(dim(xfuncs)) %in% 2:4)
        stop("Argument xfuncs is invalid: must be a 2D, 3D or 4D array.")
    dim.sig <- length(dim(xfuncs)) - 1
    if (dim(xfuncs)[1] != n)
        stop("Arguments y and xfuncs has invalid lengths: ", length(y), " and ",
             dim(xfuncs)[1], ".")
    if (dim(xfuncs)[2] != dim(xfuncs)[1 + dim.sig])
        stop("Number of rows and columns in image are not identical: ", 
             dim(xfuncs)[2], " and ", dim(xfuncs)[1 + dim.sig])
    d <- dim(xfuncs)[2]
    if (as.integer(log2(d)) != log2(d))
        stop("Argument xfuncs is invalid: the length of xfuncs must be of 
             power of 2.")
    if (sum(!min.scale %in% 0:(log2(d) - 1)) != 0)
        stop("Argument min.scale is invalid: must be integer(s) between 0 and ",
             log2(d) - 1, ".")
    if (alpha < 0 || alpha > 1)
        stop("Argument alpha s invalid: must in [0,1].")
    if (!nfold %in% 1:n)
        stop("Argument nfold is invalid: must be an integer between 1 and ", n, ".") 
    if (is.null(lambda) || 
        (length(min.scale) + length(nfeatures)+length(alpha)+length(lambda)) > 4){
    	do.cv <- TRUE
    	if (!is.null(seed)) set.seed(seed)
    	splits <- replicate(nsplit, split(sample(1:n), rep(1:nfold, length=n)))
    } else{
    	do.cv <- FALSE; nfold <- nsplit <- 1
    }
    
    wave.decomp <- switch(dim.sig, wd, imwd, wd3D)
    dec <- switch(dim.sig, decomp, decomp2d, decomp3d)
    rec <- switch(dim.sig, reconstr, reconstr2d, reconstr3d)
    
    wdobjs <- apply(xfuncs, 1, wave.decomp, filter.number = filter.number, 
                    family = wavelet.family)
    temp <- dec(wdobjs[[1]]) 
    p <- length(temp$coef)
    type.gaussian <- ifelse(p > n, "naive", "covariance")
    n.covt <- if(is.null(covt)) 0 else ncol(as.matrix(covt))
    penalty.factor <- if(pen.covt) rep(1,n.covt+p) else c(rep(0,n.covt), rep(1,p))
    cv.table <- lambda.table <- 
                array(0, dim = c(length(min.scale), length(nfeatures),
                                 length(alpha), 
                                 ifelse(is.null(lambda), 100, length(lambda))))
    dimnames(cv.table) <- list(paste("ms=", min.scale, sep=""), 
                               paste('nfeatures=', nfeatures, sep=''), 
                               paste("alpha=", alpha, sep=""), 
                               if(is.null(lambda)) NULL 
                               else paste("lambda=", lambda, sep = ""))
    if (compare.fits | (!do.cv)) {
    	fhat.table <- array(0, dim=c(d^dim.sig, nsplit*nfold, length(min.scale),
                                     length(nfeatures), length(alpha),
                                     ifelse(is.null(lambda), 100, length(lambda))))
    	dimnames(fhat.table) <- list(NULL, NULL, paste("ms =", min.scale),
                                     paste('nfeatures=', nfeatures), 
                                     paste("alpha =", alpha),
    	                             paste("lambda",if(is.null(lambda))1:100 
    	                                            else lambda))
    }
    
    for (isplit in 1 : nsplit){    	
    	if (do.cv){
    		cat('nsplit:', isplit, '\n')
    		groups <- splits[,isplit]
    	}
    	for (ims in 1 : length(min.scale)){
    		coef <- t(array(unlist(sapply(wdobjs, dec, min.scale = min.scale[ims])[1,]), dim = c(p, n)))
    		dimnames(coef) <- list(NULL, 1:p)
    		temp$callInfo$min.scale <- min.scale[ims]
    		for (ifold in 1 : nfold){
    			if (do.cv){
    				idxTest <- groups[[ifold]]
    				idxTrain <- (1:n)[-idxTest]
    			} else{
    				idxTest <- idxTrain <- 1:n
    			}
    			criteria <- apply(coef[idxTrain, ], 2, var)
    			names(criteria) <- 1:p
	            sorted <- sort(criteria, decreasing = TRUE, na.last = TRUE)[1:max(nfeatures)]
	            subset <- coef[, as.numeric(names(sorted))]
	            temp$callInfo$min.scale <- min.scale[ims]	
	            for (infeatures in 1:length(nfeatures)){
	            	coef.red <- subset[,1:nfeatures[infeatures]]
	            	cat("min.scale:", min.scale[ims], "fold:", ifold, "\tnfeatures:", nfeatures[infeatures], 
	            		    "\talpha: ")
	            	for (ialpha in 1 : length(alpha)){
	            		cat(alpha[ialpha], '\t')
	            		if(!is.null(lambda)){
	            			lambda.table[ims, infeatures, ialpha, ] <- lambda
	            		} else if(sum(lambda.table[ims, infeatures, ialpha, ] !=0)==0){
	            			obje <- glmnet(x = as.matrix(cbind(covt, coef.red)), y = y, family = family, 
	            			               alpha = alpha[ialpha], standardize = standardize, 
	            			               type.gaussian = type.gaussian, penalty.factor = penalty.factor, 
	            			               ...)
	            			templam <- range(obje$lambda)
	            			lambda.table[ims,infeatures,ialpha,] <- seq(templam[1],templam[2],length=100)
	            		}
	            		obje <- glmnet(x = as.matrix(cbind(covt, coef.red)[idxTrain,]), y = y[idxTrain], 
	            		               lambda = lambda.table[ims, infeatures, ialpha,], family = family,
	            		               alpha = alpha[ialpha], standardize = standardize, 
	            		               type.gaussian = type.gaussian, penalty.factor = penalty.factor, 
	            		               ...)
	            		if(compare.fits | (!do.cv)){
	            			theta.w <- matrix(predict(obje, s = lambda.table[ims, infeatures, ialpha, ], 
	            			                          type = 'coefficients'), ncol = dim(lambda.table)[4])	            			
	            			for (ilambda in 1 : length(lambda.table[ims, infeatures, ialpha,])){
	            				temp$coef <- rep(0, p)
	            				temp$coef[as.numeric(colnames(coef.red))] <- 
	            				    theta.w[-(1:ncol(cbind(1, covt))), ilambda]
	            				fhat.table[, (isplit-1)*nfold+ifold, ims, infeatures, ialpha, ilambda] <- 
	            				    as.vector(rec(temp))
	            			}
	            		}
	            		yhat <- predict(obje, newx = as.matrix(cbind(covt, coef.red)[idxTest, ]),
	            		                s = lambda.table[ims, infeatures, ialpha,], type = 'response')
	            		if (family == 'gaussian'){
	            			cv.table[ims, infeatures, ialpha,] <- cv.table[ims, infeatures, ialpha,] + 
	            			                                      colMeans((y[idxTest] - yhat)^2)
	            		} else if (family == 'binomial'){
	            			cv.table[ims,infeatures, ialpha, ] <- cv.table[ims, infeatures, ialpha, ] - 
	                                                              apply(as.matrix(log(yhat)[y[idxTest] == 1,]), 2, sum) - 
	                                                              apply(as.matrix(log((1-yhat))[y[idxTest] == 0, ]), 2, sum)
	            		}         
	            	} # alpha
	                cat('\n')
	            } # nfeatures
    		} # fold
    	} # min.scale
    } # split
    if (do.cv){
    	idxmin <- which(cv.table == min(cv.table[cv.table != 0], na.rm = TRUE), arr.ind = TRUE)
    	if (nrow(idxmin) > 1) idxmin = idxmin[1,]
    	min.scale <- min.scale[idxmin[1]]
    	nfeatures <- nfeatures[idxmin[2]]
    	alpha <- alpha[idxmin[3]]
    	lambda <- lambda.table[idxmin[1], idxmin[2], idxmin[3], idxmin[4]]
    	obje <- wnet(y = y, xfuncs = xfuncs, covt = covt, min.scale = min.scale, nfeatures = nfeatures, 
    	             alpha = alpha, lambda = lambda, standardize = standardize, pen.covt = pen.covt, 
    	             filter.number = filter.number, wavelet.family = wavelet.family, family = family, 
    	             nfold = 1, nsplit = 1, compare.fits = TRUE)
        if (store.cv) {
        	obje$cv.table <- 2 * cv.table / n / nfold / nsplit
        } else{
        	obje$cv.table <- 2 * min(cv.table[cv.table!=0]) / n / nfold / nsplit
        }
        if (compare.fits){
        	obje$stability <- apply(apply(fhat.table, 3:6, apply, 1, scale, center=TRUE, scale=FALSE)^2, 2:5, sum)/
        	                  (d^dim.sig*(nsplit*nfold-1))
        }
    } else{
    	obje$fhat <- array(fhat.table[,1,1,1,1,1], dim=dim(xfuncs)[-1])
    	obje$const <- theta.w[1:ncol(cbind(1, covt)), 1]
    	obje$min.scale <- min.scale
    	obje$nfeatures <- nfeatures
    	obje$alpha <- alpha
    	obje$lambda <- lambda
    	if (family == 'binomial'){
            pp = mean(y)
            logL0 <- n * (pp * log(pp) + (1-pp) * log(1-pp))
            func <- function(fhat, const){
                dim(fhat) = c(d^dim.sig, 1)
                dim(xfuncs) = c(n, d^dim.sig)
                dim(const) = c(length(const), 1)
                yhat = cbind(rep(1, n), covt) %*% const + xfuncs %*% fhat
                phat = 1 / (1 + exp(-yhat))
                sum(log(phat)[y == 1]) + sum(log(1-phat)[y == 0])
            }
            obje$Rsq = 1 - func(fhat.table[,1,1,1,1,1], obje$const)/logL0
    	}
    }
    obje$family = family
    class(obje) <- 'wnet'
    return(obje)
}