parse.predict.pfr <- function(pfr.obj, new.data){
  ## parse out new.data
  subj.new       <- new.data$subj
  covariates.new <- new.data$covariates
  funcs.new      <- new.data$funcs
  ## parse out pfr.obj
  funcs.old      <- pfr.obj$funcs
  kb.old         <- pfr.obj$kb 
  kz.old         <- pfr.obj$kz 
  nbasis.old     <- pfr.obj$nbasis
  alpha.old      <- pfr.obj$beta.covariates[1]
  beta.old       <- pfr.obj$beta.covariates[-1] ## what happens if covariates = NULL
  p.old          <- length(beta.old)
  N_subj.old     <- ifelse(is.null(subj.new), 0, ncol(pfr.obj$Z1))
  if(is.null(subj.new)){rand.int.old <- Inf
  }else                 rand.int.old <-  matrix(pfr.obj$fit$coef[c((p.old+2):(N_subj.old+p.old+1))], ncol=1)
  subj.old       <- pfr.obj$subj
  W              <- pfr.obj$BetaHat
  

  ## need to manage old and new subjects for level 1 predictions
  ## this code chunk determines which subjects were in original
  ## fit, and which ones are new.  This is important and we follow
  ## what lme() does:  everyone can have a fixed effect level 0 predicted
  ## value that utilizes no random effects; however, for level 1 individual level
  ## predictions we can only utilize random effects if the individual that is in the
  ## prediction was first in the fitted set
  subj.old.key <- unique(subj.old)
  subj.new.key <- unique(subj.new)
  in.both <- subj.old.key %in% subj.new.key
  subj.ext.key <- subj.old.key[in.both]
  rand.int.ext <- rand.int.old[in.both]
  extracted <- subj.new.key %in% subj.ext.key
  rand.int.new <- matrix(rep(NA, length(subj.new.key)), ncol=1)
  rand.int.new[extracted] <- rand.int.ext


  
  ret <- list(subj.new, covariates.new, funcs.new,
              kb.old, kz.old, nbasis.old, alpha.old, beta.old, p.old,
              N_subj.old, rand.int.old, subj.old,
              W,
              rand.int.new,
              funcs.old)
  names(ret) <- c("subj.new", "covariates.new", "funcs.new",
              "kb.old", "kz.old", "nbasis.old", "alpha.old", "beta.old", "p.old",
              "N_subj.old", "rand.int.old", "subj.old",
              "W",
              "rand.int.new",
              "funcs.old")
  ret
}
