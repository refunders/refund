wnet <- function(y, xfuncs, covt = NULL, min.scale = 0, nfeatures = NULL, alpha = 1,
                 lambda = NULL, standardize = FALSE, pen.covt = FALSE, filter.number = 10, 
                 wavelet.family = 'DaubLeAsymm', family = 'gaussian', nfold = 5, 
                 nsplit = 1, store.cv = FALSE, store.glmnet = FALSE, seed = NULL, ...){
    if (!is.null(covt)){
    	covt <- as.matrix(covt)
    }
    checkError(y = y, xfuncs = xfuncs, covt = covt, 
               vars = list(min.scale = min.scale, nfeatures = nfeatures, alpha = alpha))    
    n <- length(y)
    dim.sig <- length(dim(xfuncs)) - 1
    d <- dim(xfuncs)[2]
    multiValFlag <- c(length(min.scale) > 1, length(nfeatures) > 1, length(alpha) > 1, 
                   is.null(lambda) || (length(lambda) > 1))       
    names(multiValFlag) <- c("min.scale", "nfeatures", "alpha", "lambda")
    if (any(multiValFlag)){
    	do.cv <- TRUE
    	if (!is.null(seed)) set.seed(seed)
    	splits <- replicate(nsplit, split(sample(1:n), rep(1:nfold, length=n)))
    } else{
    	do.cv <- FALSE 
    	store.cv <- FALSE
    	nfold <- 1
    	nsplit <- 1
    }
    inm <- ifelse(length(which(multiValFlag[-4] == TRUE)) == 0, ifelse(nfold > 1, 'nfold', ''), 
                  switch(max(which(multiValFlag[-4] == TRUE)), 
                         ifelse(nfold > 1, 'nfold', 'ms'), 'nf', 'alpha', '')) 
    p <- prod(dim(xfuncs)[-1]) + ifelse(dim.sig == 2, 1, 0)
    if (is.null(nfeatures)) nfeatures <- p
    nlambda <- ifelse(do.cv, 100, 1)
    type.gaussian <- ifelse(p > n, "naive", "covariance")
    n.covt <- if(is.null(covt)) 0 else ncol(as.matrix(covt))
    penalty.factor <- if(pen.covt) rep(1,n.covt+p) else c(rep(0,n.covt), rep(1,p))
    lambda.table <- array(0, dim = c(length(min.scale), length(nfeatures), length(alpha), nlambda))
    dimnames(lambda.table) <- list(paste("min.scale=", min.scale, sep=""), paste("nfeatures=", nfeatures, sep=""),
                                   paste("alpha=", alpha, sep=""), paste("lambda", 1:nlambda, sep=""))
    cv.table <- array(0, dim = c(nsplit, nfold, length(min.scale), length(nfeatures), length(alpha), nlambda))
    dimnames(cv.table) <- list(NULL, NULL, paste("min.scale=", min.scale, sep=""), 
                               paste("nfeatures=", nfeatures, sep=""), paste("alpha=", alpha, sep=""), NULL) 
    setup <- waveletSetup(xfuncs = xfuncs, filter.number = filter.number, family = wavelet.family, 
                          min.scale = min.scale)
    if(do.cv && any(multiValFlag)){
        if (!multiValFlag["min.scale"]) cat("min.scale =", min.scale, " ")
        if (!multiValFlag["nfeatures"]) cat("nfeatures =", nfeatures, " ")
        if (!multiValFlag["alpha"]) cat("alpha =", alpha, " ")
        if (!multiValFlag["lambda"]) cat("lambda =", lambda, " ")
        cat("\n")
    }
    for (isplit in 1 : nsplit){    	
    	if (do.cv){
    		if (nsplit != 1) cat('split', isplit, '\n')
    		if (nfold > 1) groups <- splits[,isplit] else groups <- splits[isplit]
    	}
    	if (inm == 'ms') cat("min.scale = ")
    	for (ims in 1 : length(min.scale)){
    		if (multiValFlag["min.scale"]) if (inm == 'ms') cat(min.scale[ims], " ") 
    		setup$temp$callInfo$min.scale <- min.scale[ims]
    		if (inm == 'nfold') {
    			if (multiValFlag["min.scale"]) cat("min.scale =", min.scale[ims], " ")
    			cat("fold ")
    		}
    		for (ifold in 1 : nfold){
    			if (nfold > 1) if (inm == 'nfold') cat(ifold, " ") 
    			if (do.cv){
    				idxTest <- groups[[ifold]]
    				idxTrain <- (1:n)[-idxTest]
    			} else{
    				idxTest <- idxTrain <- 1:n
    			}
    			criteria <- apply(setup$coef[[ims]][idxTrain, ], 2, var)
    			names(criteria) <- 1:p
	            sorted <- sort(criteria, decreasing = TRUE, na.last = TRUE)[1:max(nfeatures)]
	            subset <- setup$coef[[ims]][, as.numeric(names(sorted))]
	            if (inm == 'nf') {
	            	if (multiValFlag["min.scale"]) cat("min.scale =", min.scale[ims], " ")
	            	if (nfold > 1) cat("fold", ifold, " ") 
	            	cat("nfeatures = ")
	            }
	            for (infeatures in 1:length(nfeatures)){
	            	if (multiValFlag["nfeatures"]) if (inm == 'nf') cat(nfeatures[infeatures], " ")
	            	coef.red <- subset[,1:nfeatures[infeatures]]
	            	if (multiValFlag["alpha"]) {
	            		if (multiValFlag["min.scale"]) cat("min.scale =", min.scale[ims], " ")
	            		if (nfold > 1) cat("fold", ifold, " ")
	            		if (multiValFlag["nfeatures"]) cat("nfeatures =", nfeatures[infeatures], " ")
	            		cat("alpha = ")
	            	}
	            	for (ialpha in 1 : length(alpha)){
	            		if (multiValFlag["alpha"]) cat(alpha[ialpha], " ")
	            		if(!is.null(lambda)){
	            	        lambda.table[ims, infeatures, ialpha, ] <- lambda          			
	            		} else if(sum(lambda.table[ims, infeatures, ialpha, ] !=0)==0){
	            			obje <- glmnet(x = as.matrix(cbind(covt, coef.red)), y = y, family = family, 
	            			               alpha = alpha[ialpha], standardize = standardize, 
	            			               type.gaussian = type.gaussian, penalty.factor = penalty.factor, 
	            			               ...)
	            			templam <- range(obje$lambda)
	            			lambda.table[ims,infeatures,ialpha,] <- seq(templam[1],templam[2],length=nlambda)
	            		}
	            		obje <- glmnet(x = as.matrix(cbind(covt, coef.red)[idxTrain,]), y = y[idxTrain], 
	            		               lambda = lambda.table[ims, infeatures, ialpha,], family = family,
	            		               alpha = alpha[ialpha], standardize = standardize, 
	            		               type.gaussian = type.gaussian, penalty.factor = penalty.factor, 
	            		               ...)
	            		yhat <- predict(obje, newx = as.matrix(cbind(covt, coef.red)[idxTest, ]),
	            		                s = lambda.table[ims, infeatures, ialpha,], type = 'response')
	            		if (do.cv){
	            			if (family == 'gaussian'){
	            			    cv.table[isplit, ifold, ims, infeatures, ialpha,] <- colMeans((y[idxTest] - yhat)^2)
	            		    } else if (family == 'binomial'){
	                            cv.table[isplit, ifold, ims,infeatures, ialpha, ] <- 
	                                            (- apply(as.matrix(log(yhat)[y[idxTest] == 1,]), 2, sum) - 
	                                              apply(as.matrix(log((1-yhat))[y[idxTest] == 0, ]), 2, sum)) / length(y[idxTest])
	            		    }  
	            		} else {
	            			theta.w <- predict(obje, s = lambda.table[ims, infeatures, ialpha, ], type = 'coefficients')
	            			setup$temp$coef <- rep(0, p)
	            			setup$temp$coef[as.numeric(colnames(coef.red))] <- theta.w[-(1:ncol(cbind(1, covt)))]
	            		    fhat <- as.vector(setup$rec(setup$temp))	            			
	            		}
	            	} # alpha
	                if (inm == 'alpha') cat('\n')
	            } # nfeatures
	            if (inm == 'nf') cat('\n')
    		} # fold
    		if (inm == 'nfold') cat('\n')
    	} # min.scale
    	if (inm == 'ms') cat('\n')
    } # split
    if (do.cv){
    	res <- waveletGetCV(cv.table = cv.table, family = family)
    	param <- array(0, dim = c(2, 4))
    	param[1, ] <- res$idxmin
    	param[2, ] <- c(min.scale[param[1,1]], nfeatures[param[1,2]], alpha[param[1,3]], 
    	                lambda.table[param[1,1], param[1,2], param[1,3], param[1,4]])
    	dimnames(param) = list(c("index", "value"), c("min.scale", "nfeatures", "alpha", "lambda"))
    	l <- wnet(y = y, xfuncs = xfuncs, covt = covt, min.scale = param[2,1], nfeatures = param[2,2], 
    	             alpha = param[2,3], lambda = param[2,4], standardize = standardize, pen.covt = pen.covt, 
    	             filter.number = filter.number, wavelet.family = wavelet.family, family = family, 
    	             nfold = 1, nsplit = 1, store.glmnet = store.glmnet, ...)
        if (store.cv) {
            l$cv.table <- res$cv.table
            l$se.cv <- res$se.cv
        }
        else {
            l$cv.table <- min(res$cv.table[res$cv.table !=0], na.rm=TRUE)
        }
        l$lambda.table <- lambda.table
        l$tuning.param <- param
    } else{
    	l <- list()
    	if (store.glmnet){
    		l$glmnet <- obje
    	}
        l$fitted.value <- yhat    
        l$coef.params <- theta.w[1:(n.covt+1),]  
    	l$fhat <- array(fhat, dim = rep(d, dim.sig))  
    	l$Rsq <- getRsq(y = y, yhat = yhat, family = family) 	
    }
    l$family <- family
    class(l) <- 'wnet'
    return(l)
}