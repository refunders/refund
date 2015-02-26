lpfr <-
function(Y, subj, covariates=NULL, funcs, kz=30, kb=30, smooth.cov=FALSE, family="gaussian", method = "REML", ...) {
	
	kb = min(kz, kb)
	n = length(Y)
	p = ifelse(is.null(covariates), 0, dim(covariates)[2])
	N_subj=length(unique(subj))

	if(is.matrix(funcs)){
		Funcs = list(length=1)
		Funcs[[1]] = funcs
	} else {
		Funcs = funcs
	}
	
	# functional predictors
	N.Pred = length(Funcs)
	
	t = phi = psi = CJ = list(length=N.Pred)
	for(i in 1:N.Pred){
		t[[i]] = seq(0, 1, length = dim(Funcs[[i]])[2])
		N_obs = length(t[[i]])
		
		# de-mean the functions
		meanFunc=apply(Funcs[[i]], 2, mean, na.rm=TRUE)
		resd=sapply(1:length(t[[i]]), function(u) Funcs[[i]][,u]-meanFunc[u])
		Funcs[[i]]=resd

		# construct and smooth covariance matrices
		G.sum <- matrix(0, N_obs, N_obs)
		G.count <- matrix(0, N_obs, N_obs)

		for(j in 1:dim(resd)[1]){
		    row.ind=j
	    	temp=resd[row.ind, ] %*% t( resd[row.ind, ])
			G.sum <- G.sum + replace(temp, which(is.na(temp)), 0)
			G.count <- G.count + as.numeric(!is.na(temp))
		}
		G <- ifelse(G.count==0, NA,  G.sum/G.count)   
		
		## get the eigen decomposition of the smoothed variance matrix
		if(smooth.cov){
			G2 <- G
			M <- length(t[[i]])
			diag(G2)= rep(NA, M)
			g2 <- as.vector(G2)
			## define a N*N knots for bivariate smoothing
			N <- 10

			## bivariate smoothing using the gamm function
			t1 <- rep(t[[i]], each=M)
			t2 <- rep(t[[i]], M)
			newdata <- data.frame(t1 = t1, t2 = t2)
			K.0  <- matrix(predict(gam(as.vector(g2) ~ te(t1, t2, k = N)), newdata), M, M) # smooth K.0
			K.0 <- (K.0 + t(K.0)) / 2    
		
			eigenDecomp <- eigen(K.0)
		} else {
			eigenDecomp <- eigen(G)
		}
	
		psi[[i]] = eigenDecomp$vectors[,1:kz]

		# set the basis to be used for beta(t)
		num=kb-2
		qtiles <- seq(0, 1, length = num + 2)[-c(1, num + 2)]
		knots <- quantile(t[[i]], qtiles)
		phi[[i]] = cbind(1, t[[i]], sapply(knots, function(k) ((t[[i]] - k > 0) * (t[[i]] - k))))

		C=matrix(0, nrow=dim(resd)[1], ncol=kz)

		for(j in 1:dim(resd)[1]){
			C[j,] <-  replace(resd[j,], which(is.na(resd[j,])), 0) %*% psi[[i]][ ,1:kz ] 
		}

		J = t(psi[[i]]) %*% phi[[i]]
		CJ[[i]] = C %*% J
	}
	
	Z1=matrix(0, nrow=length(Y), ncol=length(unique(subj)))
	for(i in 1:length(unique(subj))){
		Z1[which(subj==unique(subj)[i]),i]=1
	}

	X = cbind(1, covariates, Z1)
	for(i in 1:N.Pred){
		X = cbind(X, CJ[[i]])
	}
	
	D = list(length=N.Pred + 1)
	D[[1]] = diag(c(rep(0, 1+p), rep(1, N_subj), rep(0, length = N.Pred * (kb))))
	for(i in 2:(N.Pred+1)){
		D[[i]] = diag(c(rep(0, 1+p+N_subj), rep(0, kb*(i-2)), c(rep(0, 2), rep(1, kb-2)), rep(0, kb*(N.Pred+1-i))))
	}
	
	## fit the model
	fit = gam(Y~X-1, paraPen=list(X=D), family=family, method = method, ...)

	## get the coefficient and betaHat estimates
	coefs = fit$coef
	fitted.vals <- as.matrix(X[,1:length(coefs)]) %*% coefs
	
	beta.covariates = coefs[1:(p+1)]
	BetaHat = varBeta = varBetaHat = Bounds = list(length(N.Pred))
	for(i in 1:N.Pred){
		BetaHat[[i]] = phi[[i]] %*% coefs[(2+p+N_subj+kb*(i-1)):(1+p+N_subj+kb*(i))]
		
		## get the covariance matrix of the estimated functional coefficient
		varBeta[[i]]=fit$Vp[(2+p+N_subj+kb*(i-1)):(1+p+N_subj+kb*(i)),(2+p+N_subj+kb*(i-1)):(1+p+N_subj+kb*(i))]
		varBetaHat[[i]]=phi[[i]]%*%varBeta[[i]]%*%t(phi[[i]])
	
		## construct upper and lower bounds for betahat
		Bounds[[i]] = cbind(BetaHat[[i]] + 1.96*(sqrt(diag(varBetaHat[[i]]))), 
			BetaHat[[i]] - 1.96*(sqrt(diag(varBetaHat[[i]]))))
	}

	ranef = coefs[(2+p):(1+p+N_subj)]

		
	ret <- list(fit, fitted.vals, BetaHat, beta.covariates, ranef, X, phi, psi, varBetaHat, Bounds)
	names(ret) <- c("fit", "fitted.vals", "BetaHat", "beta.covariates", "ranef", "X", "phi", 
		"psi", "varBetaHat", "Bounds")
	ret
}

