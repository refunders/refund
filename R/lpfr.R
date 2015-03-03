##' Longitudinal penalized functional regression
##'
##' Implements longitudinal penalized functional regression (Goldsmith et al.,
##' 2010) for generalized linear functional models with scalar outcomes and
##' subject-specific random intercepts.
##'
##' Functional predictors are entered as a matrix or, in the case of multiple
##' functional predictors, as a list of matrices using the \code{funcs}
##' argument. Missing values are allowed in the functional predictors, but it
##' is assumed that they are observed over the same grid. Functional
##' coefficients and confidence bounds are returned as lists in the same order
##' as provided in the \code{funcs} argument, as are principal component and
##' spline bases.
##'
##' @param Y vector of all outcomes over all visits
##' @param subj vector containing the subject number for each observation
##' @param covariates matrix of scalar covariates
##' @param funcs matrix or list of matrices containing observed functional
##' predictors as rows. NA values are allowed.
##' @param kz dimension of principal components basis for the observed
##' functional predictors
##' @param kb dimension of the truncated power series spline basis for the
##' coefficient function
##' @param smooth.cov logical; do you wish to smooth the covariance matrix of
##' observed functions? Increases computation time, but results in smooth
##' principal components
##' @param family generalized linear model family
##' @param method method for estimating the smoothing parameters; defaults to
##' REML
##' @param ... additional arguments passed to \code{\link[mgcv]{gam}} to fit
##' the regression model.
##' @return \item{fit }{result of the call to \code{gam}} \item{fitted.vals
##' }{predicted outcomes} \item{betaHat }{list of estimated coefficient
##' functions} \item{beta.covariates }{parameter estimates for scalar
##' covariates} \item{ranef }{vector of subject-specific random intercepts}
##' \item{X }{design matrix used in the model fit} \item{phi }{list of
##' truncated power series spline bases for the coefficient functions}
##' \item{psi }{list of principal components basis for the functional
##' predictors} \item{varBetaHat }{list containing covariance matrices for the
##' estimated coefficient functions} \item{Bounds }{list of bounds of a 95\%
##' confidence interval for the estimated coefficient functions}
##' @author Jeff Goldsmith <jeff.goldsmith@@columbia.edu>
##' @references Goldsmith, J., Crainiceanu, C., Caffo, B., and Reich, D.
##' (2012). Longitudinal penalized functional regression for cognitive outcomes
##' on neuronal tract measurements. \emph{Journal of the Royal Statistical
##' Society: Series C}, 61(3), 453--469.
##' @examples
##'
##' \dontrun{
##' ##################################################################
##' # use longitudinal data to regress continuous outcomes on
##' # functional predictors (continuous outcomes only recorded for
##' # case == 1)
##' ##################################################################
##'
##' data(DTI)
##'
##' # subset data as needed for this example
##' cca = DTI$cca[which(DTI$case == 1),]
##' rcst = DTI$rcst[which(DTI$case == 1),]
##' DTI = DTI[which(DTI$case == 1),]
##'
##'
##' # note there is missingness in the functional predictors
##' apply(is.na(cca), 2, mean)
##' apply(is.na(rcst), 2, mean)
##'
##'
##' # fit two models with single functional predictors and plot the results
##' fit.cca = lpfr(Y=DTI$pasat, subj=DTI$ID, funcs = cca, smooth.cov=FALSE)
##' fit.rcst = lpfr(Y=DTI$pasat, subj=DTI$ID, funcs = rcst, smooth.cov=FALSE)
##'
##' par(mfrow = c(1,2))
##' matplot(cbind(fit.cca$BetaHat[[1]], fit.cca$Bounds[[1]]),
##'   type = 'l', lty = c(1,2,2), col = c(1,2,2), ylab = "BetaHat",
##'   main = "CCA")
##' matplot(cbind(fit.rcst$BetaHat[[1]], fit.rcst$Bounds[[1]]),
##'   type = 'l', lty = c(1,2,2), col = c(1,2,2), ylab = "BetaHat",
##'   main = "RCST")
##'
##'
##' # fit a model with two functional predictors and plot the results
##' fit.cca.rcst = lpfr(Y=DTI$pasat, subj=DTI$ID, funcs = list(cca,rcst),
##'   smooth.cov=FALSE)
##'
##' par(mfrow = c(1,2))
##' matplot(cbind(fit.cca.rcst$BetaHat[[1]], fit.cca.rcst$Bounds[[1]]),
##'   type = 'l', lty = c(1,2,2), col = c(1,2,2), ylab = "BetaHat",
##'   main = "CCA")
##' matplot(cbind(fit.cca.rcst$BetaHat[[2]], fit.cca.rcst$Bounds[[2]]),
##'   type = 'l', lty = c(1,2,2), col = c(1,2,2), ylab = "BetaHat",
##'   main = "RCST")
##' }
##' @export
##' @importFrom mgcv gam predict.gam
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

