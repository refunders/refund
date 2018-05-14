## ----ops, echo=FALSE-----------------------------------------------------
knitr::opts_chunk$set(comment=NA, warning=FALSE, message=FALSE)

## ----pkg-attach, echo = FALSE--------------------------------------------
library(FDboost)

## function based on funplot() to plot values on logscale with nice labels
## for viscosity data
funplotLogscale <- function(x, y, ylim = NULL, col=1, add=FALSE, ...){
  
  dots <- list(...)
  if(!add){
    if(is.null(ylim)) ylim <- range(y, na.rm=TRUE)
    
    plot(x, rep(1, length(x)), col = "white", ylim = ylim, yaxt = "n",
         ylab = "", xlab = "",...)
    
    mtext("time [s]", 1, line=2)#, cex=1.5)
    mtext("viscosity [mPas]", 2, line=2)#, cex=1.5)
    
    abline(h = log(1*10^(0:9)), col="gray")
    axis(2, at = log(1*10^(0:9)), labels=1*10^(0:9))
    
    #axis(4, at=log(0.001*10^(0:9)), labels=FALSE)
    
    axis(2, at=log(2*10^(0:9)), labels=FALSE)
    axis(2, at=log(3*10^(0:9)), labels=FALSE)
    axis(2, at=log(4*10^(0:9)), labels=FALSE)
    axis(2, at=log(5*10^(0:9)), labels=FALSE)
    axis(2, at=log(6*10^(0:9)), labels=FALSE)
    axis(2, at=log(7*10^(0:9)), labels=FALSE)
    axis(2, at=log(8*10^(0:9)), labels=FALSE)
    axis(2, at=log(9*10^(0:9)), labels=FALSE)
    
    #axis(4, at=seq(0, 30, by=5))
    if(diff(ylim) > 5) axis(4, at=seq(-20, 20, by=2))
    if(diff(ylim) > 2 & diff(ylim) < 5) axis(4, at=seq(-10, 10, by=1))
    if(diff(ylim) > 1 & diff(ylim) < 2) axis(4, at=seq(-10, 10, by=0.5))  
    if(diff(ylim) < 1) axis(4, at=seq(-10, 10, by=0.25))
  }
  
  funplot(x, y, add=TRUE, col=col, type="l", ...) 
  
}

# function to color-code according to a factor
getCol2 <- function(x, cols = rainbow(18)[1:length(table(x))]){
  ret <- c()
  for(i in 1:length(cols)){
    ret[x==names(table(x))[i]] <- cols[i] 
  }
  return(ret)
}

## ----load-data, echo = TRUE----------------------------------------------
# load("viscosity.RData")
data(viscosity)
str(viscosity)

## set time-interval that should be modeled
interval <- "509"

## model time until "interval"
end <- which(viscosity$timeAll==as.numeric(interval))
viscosity$vis <- log(viscosity$visAll[,1:end])
viscosity$time <- viscosity$timeAll[1:end]

## set up interactions by hand
vars <- c("T_C", "T_A", "T_B", "rspeed", "mflow")
for(v in 1:length(vars)){
  for(w in v:length(vars))
  viscosity[[paste(vars[v], vars[w], sep="_")]] <- factor(
    (viscosity[[vars[v]]]:viscosity[[vars[w]]]=="high:high")*1)
}

#str(viscosity)
names(viscosity)

## ----plot-data, echo=TRUE, fig.width=4, fig.height=3, fig.cap='Viscostiy over time with temperature of tools ($T_C$) and temerature of resin ($T_A$) color coded.', fig.align='center'----
par(mfrow=c(1,1), mar=c(3, 3, 1, 2))#, cex=1.5)
mycol <- gray(seq(0, 0.8, l=4), alpha=0.8)[c(1,3,2,4)]
int_T_CA <- with(viscosity, paste(T_C,"-", T_A, sep=""))
with(viscosity, funplotLogscale(time, vis, 
                                col=getCol2(int_T_CA, cols=mycol[4:1])))
legend("bottomright", fill=mycol, 
       legend=c("T_C low, T_A low", "T_C low, T_A high", 
                "T_C high, T_A low", "T_C high, T_A high"), cex = 0.8)

## ----complete-model, echo = TRUE, eval=FALSE-----------------------------
#  set.seed(1911)
#  modAll <- FDboost(vis ~ 1
#                    + bols(T_C) # main effects
#                    + bols(T_A)
#                    + bols(T_B)
#                    + bols(rspeed)
#                    + bols(mflow)
#                    + bols(T_C_T_A) # interactions T_WZ
#                    + bols(T_C_T_B)
#                    + bols(T_C_rspeed)
#                    + bols(T_C_mflow)
#                    + bols(T_A_T_B) # interactions T_A
#                    + bols(T_A_rspeed)
#                    + bols(T_A_mflow)
#                    + bols(T_B_rspeed) # interactions T_B
#                    + bols(T_B_mflow)
#                    + bols(rspeed_mflow), # interactions rspeed
#                    timeformula=~bbs(time, lambda=100),
#                    numInt="Riemann", family=QuantReg(),
#                    offset=NULL, offset_control = o_control(k_min = 10),
#                    data=viscosity, check0=FALSE,
#                    control=boost_control(mstop = 100, nu = 0.2))

## ----cv-complete, echo = TRUE, eval=FALSE--------------------------------
#  set.seed(1911)
#  folds <- cv(weights=rep(1, modAll$ydim[1]), type="bootstrap", B=10)
#  cvmAll <- suppressWarnings(validateFDboost(modAll, folds = folds,
#                                    getCoefCV=FALSE,
#                                    grid=seq(10, 500, by=10), mc.cores=10))
#  mstop(cvmAll) # 180
#  # modAll <- modAll[mstop(cvmAll)]
#  # summary(modAll)
#  # cvmAll

## ----stabsel1, echo = TRUE, eval=FALSE-----------------------------------
#  set.seed(1911)
#  folds <- cvMa(ydim=modAll$ydim, weights=model.weights(modAll),
#                type = "subsampling", B = 50)
#  
#  stabsel_parameters(q=5, PFER=2, p=16, sampling.type = "SS")
#  sel1 <- stabsel(modAll, q=5, PFER=2, folds=folds, grid=1:100,
#                  sampling.type="SS", mc.cores=10)
#  sel1
#  # selects effects T_C, T_A, T_C_T_A

## ----selected-model, echo = TRUE-----------------------------------------
set.seed(1911)
mod1 <- FDboost(vis ~ 1 + bols(T_C) + bols(T_A) + bols(T_C_T_A),
                timeformula = ~bbs(time, lambda = 100),  
                numInt = "Riemann", family = QuantReg(), check0 = FALSE, 
                offset = NULL, offset_control = o_control(k_min = 10),
                data = viscosity, control = boost_control(mstop = 200, nu = 0.2))

## ----cv-selected-model0, echo = TRUE-------------------------------------
mod1 <- mod1[430]

## ----cv-selected-model, echo = TRUE, eval=FALSE--------------------------
#  set.seed(1911)
#  folds <- cv(weights = rep(1, mod1$ydim[1]), type = "bootstrap", B = 10)
#  cvm1 <- validateFDboost(mod1, folds = folds, getCoefCV = FALSE,
#                          grid = seq(10, 500, by = 10), mc.cores = 10)
#  mstop(cvm1) # 430
#  mod1 <- mod1[mstop(cvm1)]
#  # summary(mod1)

## ----coefs-selected-model, echo = TRUE-----------------------------------
# set up dataframe containing systematically all variable combinations
newdata <- list(T_C=factor(c(1,1,2,2), levels=1:2, labels=c("low","high")) ,  
             T_A=factor(c(1, 2, 1, 2), levels=1:2, labels=c("low","high")), 
             T_C_T_A=factor(c(1, 1, 1, 2)), time=mod1$yind)
intercept <- 0

## effect of T_C
pred2 <- predict(mod1, which=2, newdata=newdata)
intercept <- intercept + colMeans(pred2)
pred2 <- t(t(pred2)-intercept)

## effect of T_A
pred3 <- predict(mod1, which=3, newdata=newdata)
intercept <- intercept + colMeans(pred3)
pred3 <- t(t(pred3)-colMeans(pred3))

## interaction effect T_C_T_A
pred4 <- predict(mod1, which=4, newdata=newdata)
intercept <- intercept + colMeans(pred4[3:4,])
pred4 <- t(t(pred4)-colMeans(pred4[3:4,]))

# offset+intercept 
smoothIntercept <- mod1$predictOffset(newdata$time) + intercept 

## ----plot-selected-model, echo = 11:16, fig.cap='Viscosity over time and estimated coefficient functions. On the left hand side the viscosity measures are plotted over time with temperature of tools ($T_C$) and temperature of resin ($T_A$) color-coded. On the right hand side the estimated coefficient functions are plotted.', fig.height=4, fig.width=10----
par(mfrow=c(1,2), mar=c(3, 3, 1, 2))#, cex=1.5)

mycol <- gray(seq(0, 0.8, l=4), alpha=0.8)[c(1,3,2,4)]
int_T_CA <- with(viscosity, paste(T_C,"-", T_A, sep=""))
with(viscosity, funplotLogscale(time, vis, 
                                col=getCol2(int_T_CA, cols=mycol[4:1])))
legend("bottomright", fill=mycol, 
       legend=c("T_C low, T_A low", "T_C low, T_A high", 
                "T_C high, T_A low", "T_C high, T_A high"))

mycol <- gray(seq(0, 0.5, l=3), alpha=0.8)
funplotLogscale(mod1$yind, pred2[3:4,], col=mycol[1], ylim=c(-0.5,6), lty=2, lwd=2)
lines(mod1$yind, pred3[2,], col=mycol[2], lty=3, lwd=2)
lines(mod1$yind, pred4[4,], col=mycol[3], lty=4, lwd=2)
legend("topright", lty=2:4, lwd=2, col=mycol,
       legend=c("effect T_C high", "effect T_A high", "effect T_C, T_A high"))


