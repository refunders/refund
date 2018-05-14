## ----ops, echo=FALSE-----------------------------------------------------
knitr::opts_chunk$set(comment=NA, warning=FALSE, message=FALSE)

## ----pkg-attach, echo = FALSE, warning=FALSE, message=FALSE--------------
library(FDboost)

## ----load-data, echo = TRUE----------------------------------------------
data(fuelSubset)
fuel <- fuelSubset
str(fuel)

# # normalize the wavelength to 0-1
# fuel$nir.lambda0 <- (fuel$nir.lambda - min(fuel$nir.lambda)) / 
#   (max(fuel$nir.lambda) - min(fuel$nir.lambda)) 
# fuel$uvvis.lambda0 <- (fuel$uvvis.lambda - min(fuel$uvvis.lambda)) / 
#   (max(fuel$uvvis.lambda) - min(fuel$uvvis.lambda))

# compute first derivatives as first order differences
fuel$dUVVIS <- t(apply(fuel$UVVIS, 1, diff))
fuel$dNIR <- t(apply(fuel$NIR, 1, diff)) 

# get the wavelength for the derivatives
fuel$duvvis.lambda <- fuel$uvvis.lambda[-1]
fuel$dnir.lambda <- fuel$nir.lambda[-1]
# fuel$duvvis.lambda0 <- fuel$uvvis.lambda0[-1]
# fuel$dnir.lambda0 <- fuel$nir.lambda0[-1]

## ----humidity-model, echo = TRUE, eval=FALSE-----------------------------
#  modH2O <- FDboost(h2o ~ bsignal(UVVIS, uvvis.lambda, knots=40, df=4)
#                      + bsignal(NIR, nir.lambda, knots=40, df=4)
#                      + bsignal(dUVVIS, duvvis.lambda, knots=40, df=4)
#                      + bsignal(dNIR, dnir.lambda, knots=40, df=4),
#                      timeformula=~bols(1), data=fuel)
#  
#  set.seed(212)
#  cvmH2O <- suppressWarnings(cvrisk(modH2O, grid=seq(100, 5000, by=100),
#                                folds=cv( model.weights(modH2O),
#                                type = "bootstrap", B = 10), mc.cores=10))
#  
#  par(mfrow=c(1,2))
#  plot(cvmH2O)
#  
#  modH2O[mstop(cvmH2O)]
#  #modH2O[2400]
#  
#  #### create new variable of predicted h2o
#  h2o.fit <- modH2O$fitted()
#  
#  plot(fuel$h2o, h2o.fit)
#  abline(0,1)

## ----model-spec, echo=TRUE-----------------------------------------------
formula <- formula(heatan ~ bsignal(UVVIS, uvvis.lambda, knots=40, df=4.41) 
                   + bsignal(NIR, nir.lambda, knots=40, df=4.41))

## do a model fit:
mod <- FDboost(formula, timeformula=~bols(1), data=fuel)
mod <- mod[198]

## ----cv-model-spec, echo=TRUE, eval=FALSE--------------------------------
#  ## get optimal mstop and do bootstrapping for coefficient estimates
#  set.seed(2703)
#  val <- validateFDboost(mod,
#                         folds=cv(model.weights(mod), type = "bootstrap", B = 50),
#                         grid = 10:500, mc.cores=10)
#  
#  mopt <- val$grid[which.min(colMeans(val$oobrisk))]
#  print(mopt)
#  
#  ## use optimal mstop
#  mod <- mod[mopt] # 198

## ----model-spec-plot, echo=TRUE, fig.cap='Coefficient estimates for the effects of the two spectra.', fig.height=4, fig.width=10----
par(mfrow=c(1,2))
plot(mod, which=1, lwd=2, lty=5, rug=FALSE,
     ylab="", xlab="wavelength [nm]")

plot(mod, which=2, lwd=2, lty=5, rug=FALSE, 
     ylab="", xlab="wavelength [nm]")

