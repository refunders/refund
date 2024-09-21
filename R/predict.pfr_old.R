#' @export
predict.pfr_old <- function(object, new.data=NULL, levels=NULL, ...){
  ## predict.pfr() will provide fitted values at the population and subject level for pfr objects
  ## predict.pfr() will also provide predictions for new.data at each level where applicable.
  ## if new.data is null, then return fitted.vals.level.0 and fitted.vals.level.1 from object
  if(is.null(new.data)){
    fitted.vals.level.0 <- object$fitted.vals.level.0
    fitted.vals.level.1 <- object$fitted.vals.level.1
  }else{
    ## if new.data is not null, then we need to do a few parsing steps and calculate predictions.
    par <- parse.predict.pfr(object, new.data)
    ## prep new.data for matrix multiplication; calculate subject specific scores and loadings for funcs.new
    ## pre      <- with(par,preprocess.pfr(subj=subj.new,
    ##                                     covariates=covariates.new, funcs=funcs.old, kz=kz.old, kb=kb.old,
    ##                                     nbasis=nbasis.old,
    ##                                     funcs.new=funcs.new))
    pre      <- preprocess.pfr(subj=par$subj.new,
                               covariates=par$covariates.new, funcs=par$funcs.old, kz=par$kz.old, kb=par$kb.old,
                               nbasis=par$nbasis.old,
                               funcs.new=par$funcs.new,smooth.option=par$smooth.option.old)
    #    psi            <- data.calc$psi
    #    C              <- data.calc$C
    #    Z1             <- data.calc$Z1
    ## calculate all the functional pieces into one.sum; need loop.
    one.sum        <- rep(0, nrow(pre$C[[1]]))
    for(i in 1:length(pre$C)){one.sum <- one.sum + pre$C[[i]]%*%t(pre$psi[[i]])%*%(unlist(par$W[[i]]))}
    ## calc level 0 and level 1
    ## fitted.vals.level.0 <- with(par,
    ##                             if(is.null(covariates.new)){ alpha.old + one.sum
    ##                             }else alpha.old + covariates.new%*%beta.old + one.sum)
    fitted.vals.level.0 <- if(is.null(par$covariates.new)){ par$alpha.old + one.sum
    }else{ par$alpha.old + par$covariates.new%*%par$beta.old + one.sum}
    ## for Z1, need for loop to assign NAs into Z1 and replace NA intercepts with 0
    ## this is so that the fitted.vals.level.1 will be NA if the subj.new was not in subj.old
    if(!is.null(par$subj.new)){
      for(i in 1:ncol(pre$Z1)){ pre$Z1[pre$Z1[,i]==1, i] <- par$rand.int.new[i]}
      ##par$rand.int.new[is.na(par$rand.int.new)] <- 0
      ##fitted.vals.level.1 <- fitted.vals.level.0 + pre$Z1%*%par$rand.int.new
      rs <- rowSums(pre$Z1)
      fitted.vals.level.1 <- fitted.vals.level.0 + rs
    }else{ fitted.vals.level.1 <- rep(NA, length(fitted.vals.level.0))}

  }
  ret <- list(fitted.vals.level.0, fitted.vals.level.1)
  names(ret) <- c("fitted.vals.level.0", "fitted.vals.level.1")
  ret
}

