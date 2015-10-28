postprocess.pfr <- function(fit=NULL, X=NULL, p=NULL, N_subj=NULL, phi=NULL, subj=NULL, N.Pred=NULL, kb=NULL){ 
  coefs = fit$coef
  fitted.vals <- as.matrix(X[, 1:length(coefs)]) %*% coefs
  beta.covariates = coefs[1:(p + 1)]
  ## assign level 0 (population) and level 1 (subject) fittings
  ## assigning depends on subj=NULL and equivalently N_subj =0;
  if(is.null(subj) & N_subj == 0){fitted.vals.level.0 <- fitted.vals
                                  fitted.vals.level.1 <- NULL
                                }else{ fitted.vals.level.1 <- fitted.vals
                                       ## for population level (level 0), remove subject specific columns and coefficients
                                       fitted.vals.level.0 <- as.matrix(X[,-1*c((p+2):(N_subj+p+1))])%*%fit$coef[-1*c((p+2):(N_subj+p+1))]
                                     }
  BetaHat <- varBeta <- varBetaHat <- Bounds <- list()
  for(i in 1:N.Pred){
    BetaHat[[i]] = phi[[i]] %*% coefs[-1*(1:(N_subj+p+1))][((i-1)*kb+1):(kb*i)]
    varBeta[[i]] = fit$Vp[-1*(1:(N_subj+p+1)),-1*(1:(N_subj+p+1))][((i-1)*kb+1):(kb*i),((i-1)*kb+1):(kb*i)]
    varBetaHat[[i]] = phi[[i]] %*% varBeta[[i]] %*% t(phi[[i]])
    Bounds[[i]] = cbind(BetaHat[[i]] + 1.96 * (sqrt(diag(varBetaHat[[i]]))),
            BetaHat[[i]] - 1.96 * (sqrt(diag(varBetaHat[[i]]))))
  }

  ## old pfr (v. XX.XX.Y) would not return C, J, or CJ.  We do that here:
  ret <- list(fit, fitted.vals, fitted.vals.level.0, fitted.vals.level.1,
              BetaHat, beta.covariates,
              varBetaHat, Bounds)
  names(ret) <- c("fit", "fitted.vals", "fitted.vals.level.0", "fitted.vals.level.1",
                  "BetaHat", "beta.covariates",
                  "varBetaHat", "Bounds")
  ret
}
