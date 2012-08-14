ccb.fpc <-
function(Y, nbasis = 10, pve = .99, n.boot = 100, simul = FALSE, sim.alpha = .95){
  require(mgcv)
  require(MASS)
  
  set.seed(10)
  
  D = dim(Y)[2]   # size of grid
  I = dim(Y)[1]   # number of curves

  ## create lists to store curve estimates and variances
  Yhat.boot = MODEL.VAR = list(length = I)
  for(i in 1:I){
    Yhat.boot[[i]] = matrix(NA, n.boot, D)
    MODEL.VAR[[i]] = matrix(0, D, D)
  }
  
  ## begin bootstrap sampling
  n.succ = 0
  for(i.boot in 1:n.boot){
    set.seed(i.boot)
    cat("Iteration:", i.boot, "\n")    # print iteration number
    
    ## draw bootstrap sample
    boot.samp = sample(1:I, I, replace = TRUE)
    Y.Boot = Y[boot.samp,]

    ## do decomposition for this sample; predict curves for full data
    Fit.Iter = try(fpca.sc(Y = Y.Boot, Y.pred = Y, nbasis = nbasis, pve = pve, var = TRUE, simul = FALSE))

    ## save estimates and variances
    if(class(Fit.Iter) != "try-error"){
      n.succ = n.succ + 1
      for(i.subj in 1:I){
        Yhat.boot[[i.subj]][i.boot,] = Fit.Iter$Yhat[i.subj,]
        MODEL.VAR[[i.subj]] = MODEL.VAR[[i.subj]] + Fit.Iter$VarMats[[i.subj]]
      }    
    }
  }

  ## create matrices to store results
  VarMats = list(length = I)
  for(i in 1:I){
    VarMats[[i]] = matrix(NA, D, D)
  }
  Yhat = matrix(NA, I, D)
  diag.var = matrix(NA, I, D)
  crit.val = rep(0, I)
  
  ## for all subjects, estimate total variance using iterated variance formula
  for(i.subj in 1:I){
    Exp.Model.Var = MODEL.VAR[[i.subj]]/n.succ           # E(Var | FPC decomp)
    Var.Model.Exp = var(Yhat.boot[[i.subj]], na.rm=TRUE)     # Var(E | FPC decomp)
    VarMats[[i.subj]] = Exp.Model.Var + Var.Model.Exp    # E(Var | FPC decomp) + Var(E | FPC decomp)
      
    diag.var[i.subj, ] = diag(VarMats[[i.subj]])
    Yhat[i.subj,] = apply(Yhat.boot[[i.subj]], 2, mean, na.rm = TRUE) # E(E | FPC decomp)

    ## estimate critical values for simultaneous intervals
    if(simul){
      norm.samp = mvrnorm(2500, mu = rep(0, D), Sigma = VarMats[[i.subj]])/
        matrix(sqrt(diag(VarMats[[i.subj]])), nrow = 2500, ncol = D, byrow = TRUE)
      crit.val[i.subj] = quantile(apply(abs(norm.samp), 1, max), sim.alpha)      
    }
  }

  ## return results
  if(simul){
    ret = list(Yhat, Yhat.boot, diag.var, VarMats, crit.val)
    names(ret)= c("Yhat", "Yhat.boot", "diag.var", "VarMats", "crit.val")
  } else if(!simul){
    ret = list(Yhat, Yhat.boot, diag.var, VarMats)
    names(ret)= c("Yhat", "Yhat.boot", "diag.var", "VarMats")
  }

  return(ret)
}
