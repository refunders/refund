pfr  <- function (Y, subj=NULL, covariates = NULL, funcs, argvals = NULL, kz = 10, kb = 30, nbasis=10,
                            family = "gaussian", method="REML", smooth.option="fpca.sc", pve=0.99,...)
{
  if (!is.null(argvals))
    stop("argvals is not supported in the current version of refund.")
  ## Step 1:
  ## parse formulae, etc.
  ## parse.pfr() In progress.

  ## Step 2:
  ## Preprocess in prep for gam() fit.
  pre <- preprocess.pfr(subj=subj, covariates=covariates, funcs=funcs, kz=kz, kb=kb, nbasis=nbasis, smooth.option=smooth.option, pve=pve)

  ## Step 3:
  ## gam() fit.
  fit = with(pre, gam(Y ~ X - 1, paraPen = list(X = D), method = method, family = family, ...))

  ## Step 4:
  ## Postprocess objects within "fit" to be of use to user and rlrt.pfr(), predict.pfr(), plot.pfr()
  pos <- postprocess.pfr(fit=fit, X=pre$X, p=pre$p, N_subj=pre$N_subj, phi=pre$phi, subj=subj,
                         N.Pred=pre$N.Pred, kb=kb)  

  ## Step 5:  return everything.
  ret <- c(pos, pre, list(Y=Y))
  ret
}
