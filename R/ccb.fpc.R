##' Corrected confidence bands using functional principal components
##'
##' Uses iterated expectation and variances to obtain corrected estimates and
##' inference for functional expansions.
##'
##' To obtain corrected curve estimates and variances, this function accounts
##' for uncertainty in FPC decomposition objects. Observed curves are
##' resampled, and a FPC decomposition for each sample is constructed. A
##' mixed-model framework is used to estimate curves and variances conditional
##' on each decomposition, and iterated expectation and variances combines both
##' model-based and decomposition-based uncertainty.
##'
##' @param Y matrix of observed functions for which estimates and covariance
##' matrices are desired.
##' @param argvals numeric; function argument.
##' @param nbasis number of splines used in the estimation of the mean function
##' and the bivariate smoothing of the covariance matrix
##' @param pve proportion of variance explained used to choose the number of
##' principal components to be included in the expansion.
##' @param n.boot number of bootstrap iterations used to estimate the
##' distribution of FPC decomposition objects.
##' @param simul TRUE or FALSE, indicating whether critical values for
##' simultaneous confidence intervals should be estimated
##' @param sim.alpha alpha level of the simultaneous intervals.
##' @return \item{Yhat }{a matrix whose rows are the estimates of the curves in
##' \code{Y}.} \item{Yhat.boot }{a list containing the estimated curves within
##' each bootstrap iteration.} \item{diag.var }{diagonal elements of the
##' covariance matrices for each estimated curve.} \item{VarMats }{a list
##' containing the estimated covariance matrices for each curve in \code{Y}.}
##' \item{crit.val }{estimated critical values for constructing simultaneous
##' confidence intervals.}
##' @author Jeff Goldsmith \email{jeff.goldsmith@@columbia.edu}
##' @references Goldsmith, J., Greven, S., and Crainiceanu, C. (2013).
##' Corrected confidence bands for functional data using principal components.
##' \emph{Biometrics}, 69(1), 41--51.
##' @examples
##'
##' \dontrun{
##' data(cd4)
##'
##' # obtain a subsample of the data with 25 subjects
##' set.seed(1236)
##' sample = sample(1:dim(cd4)[1], 25)
##' Y.sub = cd4[sample,]
##'
##' # obtain a mixed-model based FPCA decomposition
##' Fit.MM = fpca.sc(Y.sub, var = TRUE, simul = TRUE)
##'
##' # use iterated variance to obtain curve estimates and variances
##' Fit.IV = ccb.fpc(Y.sub, n.boot = 25, simul = TRUE)
##'
##' # for one subject, examine curve estimates, pointwise and simultaneous itervals
##' EX = 2
##' EX.IV =  cbind(Fit.IV$Yhat[EX,],
##'       Fit.IV$Yhat[EX,] + 1.96 * sqrt(Fit.IV$diag.var[EX,]),
##'       Fit.IV$Yhat[EX,] - 1.96 * sqrt(Fit.IV$diag.var[EX,]),
##'       Fit.IV$Yhat[EX,] + Fit.IV$crit.val[EX] * sqrt(Fit.IV$diag.var[EX,]),
##'       Fit.IV$Yhat[EX,] - Fit.IV$crit.val[EX] * sqrt(Fit.IV$diag.var[EX,]))
##'
##' EX.MM =  cbind(Fit.MM$Yhat[EX,],
##'       Fit.MM$Yhat[EX,] + 1.96 * sqrt(Fit.MM$diag.var[EX,]),
##'       Fit.MM$Yhat[EX,] - 1.96 * sqrt(Fit.MM$diag.var[EX,]),
##'       Fit.MM$Yhat[EX,] + Fit.MM$crit.val[EX] * sqrt(Fit.MM$diag.var[EX,]),
##'       Fit.MM$Yhat[EX,] - Fit.MM$crit.val[EX] * sqrt(Fit.MM$diag.var[EX,]))
##'
##' # plot data for one subject, with curve and interval estimates
##' d = as.numeric(colnames(cd4))
##' plot(d[which(!is.na(Y.sub[EX,]))], Y.sub[EX,which(!is.na(Y.sub[EX,]))], type = 'o',
##'   pch = 19, cex=.75, ylim = range(0, 3400), xlim = range(d),
##'     xlab = "Months since seroconversion", lwd = 1.2, ylab = "Total CD4 Cell Count",
##'       main = "Est. & CI - Sampled Data")
##'
##' matpoints(d, EX.IV, col = 2, type = 'l', lwd = c(2, 1, 1, 1, 1), lty = c(1,1,1,2,2))
##' matpoints(d, EX.MM, col = 4, type = 'l', lwd = c(2, 1, 1, 1, 1), lty = c(1,1,1,2,2))
##'
##' legend("topright", c("IV Est", "IV PW Int", "IV Simul Int",
##'     expression(paste("MM - ", hat(theta), " Est", sep = "")),
##'     expression(paste("MM - ", hat(theta), " PW Int", sep = "")),
##'     expression(paste("MM - ", hat(theta), " Simul Int", sep = ""))),
##'     lty=c(1,1,2,1,1,2), lwd = c(2.5,.75,.75,2.5,.75,.75),
##'     col = c("red","red","red","blue","blue","blue"))
##' }
##' @importFrom MASS mvrnorm
##' @importFrom stats quantile
##' @export
ccb.fpc <-
function(Y, argvals=NULL, nbasis = 10, pve = .99, n.boot = 100, simul = FALSE, sim.alpha = .95){
  
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
    message("Iteration:", i.boot, "\n")    # print iteration number

    ## draw bootstrap sample
    boot.samp = sample(1:I, I, replace = TRUE)
    Y.Boot = Y[boot.samp,]

    ## do decomposition for this sample; predict curves for full data
    Fit.Iter = try(fpca.sc(Y = Y.Boot, Y.pred = Y, argvals = argvals, nbasis = nbasis, pve = pve, var = TRUE, simul = FALSE))

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