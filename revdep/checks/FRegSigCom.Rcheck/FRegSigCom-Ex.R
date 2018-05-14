pkgname <- "FRegSigCom"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "FRegSigCom-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('FRegSigCom')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("air")
### * air

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: air
### Title: Air quality data
### Aliases: air
### Keywords: datasets

### ** Examples

 data(air)
 str(air)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("air", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("cv.ff.interaction")
### * cv.ff.interaction

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: cv.ff.interaction
### Title: Cross-validation for function-on-function regression models with
###   specified main effects and two-way interaction terms
### Aliases: cv.ff.interaction

### ** Examples


library(FRegSigCom)
data(ocean)

Y=ocean$Salinity
X=list()
X[[1]]=ocean$Potential.density
X[[2]]=ocean$Temperature
X[[3]]=ocean$Oxygen
X[[4]]=ocean$Chloropigment
n.curves=length(X)
ntot=dim(Y)[1]
ntrain=50
ntest=ntot-ntrain
X.uncent=X
for(i in 1:n.curves){
  X[[i]]=scale(X.uncent[[i]],center=TRUE, scale=FALSE)
}
lengthX=dim(X[[1]])[2]
lengthY=dim(Y)[2]
t.x=seq(0,1,length=lengthX)
t.y=seq(0,1,length=lengthY)
I.train=sample(1:ntot, ntrain)
X.train=list()
X.test=list()
t.x.all=list()
for(i in 1:n.curves)
{
  X.train[[i]]=X[[i]][I.train,]
  X.test[[i]]=X[[i]][-I.train,]
  t.x.all[[i]]=t.x
}
Y.train=Y[I.train, ]
Y.test=Y[-I.train, ]


#########################################################################
# an interaction model with given main effects and two-way interactions
#########################################################################
main.effect=c(1,2,3)
inter.effect=rbind(c(1,1),c(1,2))
fit.inter=cv.ff.interaction(X.train, Y.train, t.x.all, t.y,  main.effect, inter.effect)
Y.pred=pred.ff.interaction(fit.inter,  X.test)
error.inter=mean((Y.pred-Y.test)^2)
print(c("error.inter=", error.inter))
#coef.obj=getcoef.ff.interaction(fit.inter)
#str(coef.obj)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("cv.ff.interaction", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("cv.hd")
### * cv.hd

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: cv.hd
### Title: Cross-validation for sparse linear function-on-function
###   regression with a large number of functional predictors
### Aliases: cv.hd

### ** Examples

#########################################
#toy example using the air quality data with p=1
#########################################
data(air)
t.x=seq(0,1,length=24)
t.y=seq(0,1,length=24)
air.cv=cv.hd(X=list(air[[2]][1:20,]), Y=air[[1]][1:20,], list(t.x), t.y,
             K.cv=2, s.n.basis = 8, t.n.basis = 8)
air.pred=pred.hd(air.cv, list(air[[2]][1:2,]))
predict.error=mean((air.pred-air[[1]][1:2,])^2)
print(c("predict error=", predict.error))





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("cv.hd", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("cv.nonlinear")
### * cv.nonlinear

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: cv.nonlinear
### Title: Cross-validation for nonlinear function-on-function regression
### Aliases: cv.nonlinear

### ** Examples

#########################################
# toy example
# fit a nonlinear regression model with p=1 in the air quality data
#########################################

ptm <- proc.time()
data(air)
t.x=seq(0,1,length=24)
t.y=seq(0,1,length=24)
air.nonlinear.cv=cv.nonlinear(X=list(air[[2]][1:60,]), Y=air[[1]][1:60,], list(t.x),
           t.y, K.cv=2, s.n.basis = 8, x.n.basis = 10, t.n.basis = 8)
air.pred=pred.nonlinear(air.nonlinear.cv, list(air[[2]][1:2,]))
est.coefficient=getcoef.nonlinear(air.nonlinear.cv)
print(proc.time()-ptm)






base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("cv.nonlinear", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("cv.sigcom")
### * cv.sigcom

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: cv.sigcom
### Title: Cross-validation for linear function-on-function regression
### Aliases: cv.sigcom

### ** Examples

#################################################################
# Example 1: function-on-function model without scalar predictors
#########################################################

ptm <- proc.time()
library(MASS)

n.curves <- 2 #number of predictor curves
nnew <- 50 # number of observations in new data
ntrain <- 10 # number of observations in training data
ntot <- nnew+ntrain
t.x <- seq(0,2,length=50) # all the four predictor curves are observed
              # at 50 equally spaced points in [0,2].
t.y <- seq(0,1,length=50) # the response curve is observed at 50
              # equally spaced points in [0,1].

# functions for mu(t) and beta_i(s,t)
mu.fun <- function(t){
    2*exp(-(t-1)^2)*sin(2*pi*t)
}
beta.fun.1 <- function(s,t)
{
    sin(1.5*pi*s)*sin( pi*t)
}
beta.fun.2 <- function(s,t)
{
    exp(-5*(s-0.5)^2-3*(t-0.5))+2*exp(-5*(s-1.5)^2-3*(t-0.5))
}

#generate the predictor curves using sin() and cos() basis functions.
 X=list()
 Sigma=matrix(0.5, 2, 2)
 diag(Sigma)=1
 sins=lapply(1:5, function(i){t(sin(i*pi*t.x))/i})
 coss=lapply(1:5, function(i){t(cos(i*pi*t.x))/i})
 coefs1=lapply(1:5, function(i){mvrnorm(ntot,rep(0,2),Sigma)})
 coefs2=lapply(1:5, function(i){mvrnorm(ntot,rep(0,2),Sigma)})
 X[[1]]=Reduce("+",lapply(1:5, function(i){coefs1[[i]][,1] %*% sins[[i]]
             +coefs2[[i]][,1] %*% coss[[i]]}))
 X[[2]]=Reduce("+",lapply(1:5, function(i){coefs1[[i]][,2] %*% sins[[i]]
             +coefs2[[i]][,2] %*% coss[[i]]}))
 mu.val=mu.fun(t.y)
 beta.val=list()
 beta.val[[1]]<- outer(t.x,t.y,beta.fun.1)
 beta.val[[2]]<- outer(t.x,t.y,beta.fun.2)
# generate sample curves for the noise and the respose function.
 E<-matrix(rnorm(ntot*length(t.y)),ntot, length(t.y))
 delta<-(max(t.x)-min(t.x))/length(t.x)
 Y<-t(sapply(1:ntot, function(i){mu.val+(X[[1]][i,]%*% beta.val[[1]]
   +X[[2]][i,]%*% beta.val[[2]])*(max(t.x)-min(t.x))/length(t.x)+ E[i,]}))

# gereate the training data and perform CV
 X.train=lapply(1:n.curves, function(j){X[[j]][1:ntrain,]})
 t.x.list=lapply(1:n.curves,function(j){t.x})
 Y.train <- Y[1:ntrain, ]
 fit.cv=cv.sigcom(X.train, Y.train, t.x.list, t.y, s.n.basis=20, t.n.basis=20)

# prediction and estimation error
 X.new=lapply(1:n.curves, function(j){X[[j]][-(1:ntrain),]})
 Y.new <- Y[-(1:ntrain), ]
 E.new=E[-(1:ntrain), ]
 Y.pred=pred.sigcom(fit.cv, X.new)
 error <- mean((Y.pred-Y.new)^2)
 print(c(" prediction error=", error))
 est.error <- mean((Y.pred-Y.new+E.new)^2)
# extract the estimated intercept and coefficient functions
 print(c("estimation error for regression function (or signal function)=", est.error))
 coef.est.all=getcoef.sigcom(fit.cv)
 mu.est=coef.est.all[[1]]
 beta.est=coef.est.all[[2]]
print(proc.time()-ptm)

#################################################################
# Example 2: function-on-function model with scalar predictors
#########################################################


#################################################################
#Example 3: application to the DTI data in 'refund' package
#########################################################




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("cv.sigcom", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("getcoef.ff.interaction")
### * getcoef.ff.interaction

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: getcoef.ff.interaction
### Title: Get the estimated coefficient functions for function-on-function
###   interaction model
### Aliases: getcoef.ff.interaction

### ** Examples
 #See the examples in cv.ff.interaction() and step.ff.interaction().


base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("getcoef.ff.interaction", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("getcoef.hd")
### * getcoef.hd

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: getcoef.hd
### Title: Get the estimated intercept and coefficient functions for sparse
###   linear FOF models
### Aliases: getcoef.hd

### ** Examples
 #See the examples in cv.hd().


base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("getcoef.hd", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("getcoef.nonlinear")
### * getcoef.nonlinear

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: getcoef.nonlinear
### Title: Get the estimated intercept and nonlinear functions
### Aliases: getcoef.nonlinear

### ** Examples
#See the examples in cv.nonlinear().



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("getcoef.nonlinear", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("getcoef.sigcom")
### * getcoef.sigcom

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: getcoef.sigcom
### Title: Get the estimated intercept and coefficient functions for linear
###   FOF models
### Aliases: getcoef.sigcom

### ** Examples
#See the examples in cv.sigcom().



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("getcoef.sigcom", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ocean")
### * ocean

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ocean
### Title: Hawaii ocean data
### Aliases: ocean

### ** Examples

 data(ocean)
 str(ocean)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ocean", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("pred.ff.interaction")
### * pred.ff.interaction

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: pred.ff.interaction
### Title: Prediction for a linear FOF regression model with two-way
###   interactions
### Aliases: pred.ff.interaction

### ** Examples
 #See the examples in cv.ff.interaction() and step.ff.interaction().


base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("pred.ff.interaction", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("pred.hd")
### * pred.hd

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: pred.hd
### Title: Prediction for sparse linear function-on-function regression
### Aliases: pred.hd

### ** Examples
 #See the examples in cv.hd().


base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("pred.hd", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("pred.nonlinear")
### * pred.nonlinear

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: pred.nonlinear
### Title: Prediction for nonlinear function-on-function regression
### Aliases: pred.nonlinear

### ** Examples
#See the examples in cv.nonlinear().


base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("pred.nonlinear", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("pred.sigcom")
### * pred.sigcom

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: pred.sigcom
### Title: Prediction for linear function-on-function regression using
###   signal compression
### Aliases: pred.sigcom

### ** Examples
#See the examples in cv.sigcom().



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("pred.sigcom", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("step.ff.interaction")
### * step.ff.interaction

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: step.ff.interaction
### Title: Stepwise variable selection procedure for FOF regression models
###   with two-way interactions
### Aliases: step.ff.interaction

### ** Examples


library(FRegSigCom)
data(ocean)

Y=ocean$Salinity
X=list()
X[[1]]=ocean$Potential.density
X[[2]]=ocean$Temperature
X[[3]]=ocean$Oxygen
n.curves=length(X)
ntot=dim(Y)[1]
ntrain=50
ntest=ntot-ntrain
X.uncent=X
for(i in 1:n.curves){
  X[[i]]=scale(X.uncent[[i]],center=TRUE, scale=FALSE)
}
lengthX=dim(X[[1]])[2]
lengthY=dim(Y)[2]
t.x=seq(0,1,length=lengthX)
t.y=seq(0,1,length=lengthY)
I.train=sample(1:ntot, ntrain)
X.train=list()
X.test=list()
t.x.all=list()
for(j in 1:n.curves){
  X.train[[j]]=X[[j]][I.train,]
  X.test[[j]]=X[[j]][-I.train,]
  t.x.all[[j]]=t.x
}
Y.train=Y[I.train, ]
Y.test=Y[-I.train, ]


###############################
#model selection
###############################

fit.step=step.ff.interaction(X.train, Y.train, t.x.all, t.y)
Y.pred=pred.ff.interaction(fit.step,  X.test)
error.selected=mean((Y.pred-Y.test)^2)
print(c("error.selected=", error.selected))
#coef.obj=getcoef.ff.interaction(fit.step)
#str(coef.obj)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("step.ff.interaction", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
