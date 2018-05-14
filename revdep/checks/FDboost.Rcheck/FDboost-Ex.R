pkgname <- "FDboost"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "FDboost-Ex.timings", pos = 'CheckExEnv')
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
library('FDboost')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("FDboost")
### * FDboost

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: FDboost
### Title: Model-based Gradient Boosting for Functional Response
### Aliases: FDboost
### Keywords: models, nonlinear

### ** Examples

######## Example for function-on-scalar-regression 
data("viscosity", package = "FDboost") 
## set time-interval that should be modeled
interval <- "101"

## model time until "interval" and take log() of viscosity
end <- which(viscosity$timeAll == as.numeric(interval))
viscosity$vis <- log(viscosity$visAll[,1:end])
viscosity$time <- viscosity$timeAll[1:end]
# with(viscosity, funplot(time, vis, pch = 16, cex = 0.2))

## fit median regression model with 100 boosting iterations,
## step-length 0.4 and smooth time-specific offset
## the factors are coded such that the effects are zero for each timepoint t
## no integration weights are used!
mod1 <- FDboost(vis ~ 1 + bolsc(T_C, df = 2) + bolsc(T_A, df = 2),
               timeformula = ~ bbs(time, df = 4),
               numInt = "equal", family = QuantReg(),
               offset = NULL, offset_control = o_control(k_min = 9),
               data = viscosity, control=boost_control(mstop = 100, nu = 0.4))

## Not run: 
##D   #### find optimal mstop over 5-fold bootstrap, small number of folds for example
##D   #### do the resampling on the level of curves
##D   
##D   ## possibility 1: smooth offset and transformation matrices are refitted 
##D   set.seed(123)
##D   appl1 <- applyFolds(mod1, folds = cv(rep(1, length(unique(mod1$id))), B = 5), 
##D                       grid = 1:500)
##D   ## plot(appl1)
##D   mstop(appl1)
##D   mod1[mstop(appl1)]
##D   
##D   ## possibility 2: smooth offset is refitted, 
##D   ## computes oob-risk and the estimated coefficients on the folds
##D   set.seed(123)
##D   val1 <- validateFDboost(mod1, folds = cv(rep(1, length(unique(mod1$id))), B = 5), 
##D                         grid = 1:500)
##D   ## plot(val1)
##D   mstop(val1)
##D   mod1[mstop(val1)]
##D 
##D   ## possibility 3: very efficient 
##D   ## using the function cvrisk; be careful to do the resampling on the level of curves
##D   folds1 <- cvLong(id = mod1$id, weights = model.weights(mod1), B = 5)
##D   cvm1 <- cvrisk(mod1, folds = folds1, grid = 1:500)
##D   ## plot(cvm1)
##D   mstop(cvm1)
##D   
##D ## look at the model
##D summary(mod1)
##D coef(mod1)
##D plot(mod1)
##D plotPredicted(mod1, lwdPred = 2)
## End(Not run)

######## Example for scalar-on-function-regression 
data("fuelSubset", package = "FDboost")

## center the functional covariates per observed wavelength
fuelSubset$UVVIS <- scale(fuelSubset$UVVIS, scale = FALSE)
fuelSubset$NIR <- scale(fuelSubset$NIR, scale = FALSE)

## to make mboost:::df2lambda() happy (all design matrix entries < 10)
## reduce range of argvals to [0,1] to get smaller integration weights
fuelSubset$uvvis.lambda <- with(fuelSubset, (uvvis.lambda - min(uvvis.lambda)) / 
                                          (max(uvvis.lambda) - min(uvvis.lambda) ))
fuelSubset$nir.lambda <- with(fuelSubset, (nir.lambda - min(nir.lambda)) / 
                                          (max(nir.lambda) - min(nir.lambda) )) 

## model fit with scalar response 
## include no intercept as all base-learners are centered around 0
mod2 <- FDboost(heatan ~ bsignal(UVVIS, uvvis.lambda, knots = 40, df = 4, check.ident = FALSE) 
               + bsignal(NIR, nir.lambda, knots = 40, df = 4, check.ident = FALSE), 
               timeformula = NULL, data = fuelSubset, control = boost_control(mstop = 200)) 
               
## additionally include a non-linear effect of the scalar variable h2o 
mod2s <- FDboost(heatan ~ bsignal(UVVIS, uvvis.lambda, knots = 40, df = 4, check.ident = FALSE) 
               + bsignal(NIR, nir.lambda, knots = 40, df = 4, check.ident = FALSE) 
               + bbs(h2o, df = 4), 
               timeformula = NULL, data = fuelSubset, control = boost_control(mstop = 200)) 
               
## alternative model fit as FLAM model with scalar response; as timeformula = ~ bols(1)  
## adds a penalty over the index of the response, i.e., here a ridge penalty
## thus, mod2f and mod2 have different penalties 
mod2f <- FDboost(heatan ~ bsignal(UVVIS, uvvis.lambda, knots = 40, df = 4, check.ident = FALSE) 
               + bsignal(NIR, nir.lambda, knots = 40, df = 4, check.ident = FALSE), 
               timeformula = ~ bols(1), data = fuelSubset, control = boost_control(mstop = 200))
               
## Not run: 
##D    
##D   ## bootstrap to find optimal mstop takes some time
##D   set.seed(123)      
##D   folds2 <- cv(weights = model.weights(mod2), B = 10)     
##D   cvm2 <- cvrisk(mod2, folds = folds2, grid = 1:1000)
##D   mstop(cvm2) ## mod2[327]
##D   summary(mod2) 
##D   ## plot(mod2)
## End(Not run)

## Example for function-on-function-regression 
if(require(fda)){

  data("CanadianWeather", package = "fda")
  CanadianWeather$l10precip <- t(log(CanadianWeather$monthlyPrecip))
  CanadianWeather$temp <- t(CanadianWeather$monthlyTemp)
  CanadianWeather$region <- factor(CanadianWeather$region)
  CanadianWeather$month.s <- CanadianWeather$month.t <- 1:12
  
  ## center the temperature curves per time-point
  CanadianWeather$temp <- scale(CanadianWeather$temp, scale = FALSE)
  rownames(CanadianWeather$temp) <- NULL ## delete row-names
  
  ## fit model with cyclic splines over the year
  mod3 <- FDboost(l10precip ~ bols(region, df = 2.5, contrasts.arg = "contr.dummy") 
                   + bsignal(temp, month.s, knots = 11, cyclic = TRUE, 
                             df = 2.5, boundary.knots = c(0.5,12.5), check.ident = FALSE), 
                  timeformula = ~ bbs(month.t, knots = 11, cyclic = TRUE, 
                                      df = 3, boundary.knots = c(0.5, 12.5)), 
                  offset = "scalar", offset_control = o_control(k_min = 5), 
                  control = boost_control(mstop = 60), 
                  data = CanadianWeather) 
 
 ## Not run: 
##D                   
##D    #### find the optimal mstop over 5-fold bootstrap 
##D    ## using the function applyFolds 
##D    set.seed(123)
##D    folds3 <- cv(rep(1, length(unique(mod3$id))), B = 5)
##D    appl3 <- applyFolds(mod3, folds = folds3, grid = 1:200)
##D  
##D    ## use function cvrisk; be careful to do the resampling on the level of curves
##D    set.seed(123)
##D    folds3long <- cvLong(id = mod3$id, weights = model.weights(mod3), B = 5)
##D    cvm3 <- cvrisk(mod3, folds = folds3long, grid = 1:200)
##D    mstop(cvm3) ## mod3[64]
##D    
##D    summary(mod3)
##D    ## plot(mod3, pers = TRUE)
##D  
## End(Not run)
}

######## Example for functional response observed on irregular grid
######## Delete part of observations in viscosity data-set
data("viscosity", package = "FDboost")
## set time-interval that should be modeled
interval <- "101"

## model time until "interval" and take log() of viscosity
end <- which(viscosity$timeAll == as.numeric(interval))
viscosity$vis <- log(viscosity$visAll[,1:end])
viscosity$time <- viscosity$timeAll[1:end]
# with(viscosity, funplot(time, vis, pch = 16, cex = 0.2))

## only keep one eighth of the observation points
set.seed(123)
selectObs <- sort(sample(x = 1:(64*46), size = 64*46/4, replace = FALSE))
dataIrregular <- with(viscosity, list(vis = c(vis)[selectObs], 
                                      T_A = T_A, T_C = T_C,  
                                      time = rep(time, each = 64)[selectObs], 
                                      id = rep(1:64, 46)[selectObs]))

## fit median regression model with 50 boosting iterations,
## step-length 0.4 and smooth time-specific offset
## the factors are in effect coding -1, 1 for the levels
## no integration weights are used!
mod4 <- FDboost(vis ~ 1 + bols(T_C, contrasts.arg = "contr.sum", intercept = FALSE)
                + bols(T_A, contrasts.arg = "contr.sum", intercept=FALSE),
                timeformula = ~ bbs(time, lambda = 100), id = ~id, 
                numInt = "Riemann", family = QuantReg(),
                offset = NULL, offset_control = o_control(k_min = 9),
                data = dataIrregular, control = boost_control(mstop = 50, nu = 0.4))
## summary(mod4)
## plot(mod4)
## plotPredicted(mod4, lwdPred = 2)

## Not run: 
##D   ## Find optimal mstop, small grid/low B for a fast example
##D   set.seed(123)
##D   folds4 <- cv(rep(1, length(unique(mod4$id))), B = 3)
##D   appl4 <- applyFolds(mod4, folds = folds4, grid = 1:50)
##D   ## val4 <- validateFDboost(mod4, folds = folds4, grid = 1:50)
##D 
##D   set.seed(123)
##D   folds4long <- cvLong(id = mod4$id, weights = model.weights(mod4), B = 3)
##D   cvm4 <- cvrisk(mod4, folds = folds4long, grid = 1:50)
##D   mstop(cvm4)
## End(Not run)

## Be careful if you want to predict newdata with irregular response,  
## as the argument index is not considered in the prediction of newdata. 
## Thus, all covariates have to be repeated according to the number of observations 
## in each response trajectroy. 
## Predict four response curves with full time-observations 
## for the four combinations of T_A and T_C. 
newd <- list(T_A = factor(c(1,1,2,2), levels = 1:2, 
                        labels = c("low", "high"))[rep(1:4, length(viscosity$time))], 
             T_C = factor(c(1,2,1,2), levels = 1:2, 
                        labels = c("low", "high"))[rep(1:4, length(viscosity$time))], 
             time = rep(viscosity$time, 4))
             
pred <- predict(mod4, newdata = newd)
## funplot(x = rep(viscosity$time, 4), y = pred, id = rep(1:4, length(viscosity$time)))
                  
                



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("FDboost", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("FDboostLSS")
### * FDboostLSS

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: FDboostLSS
### Title: Model-based Gradient Boosting for Functional GAMLSS
### Aliases: FDboostLSS
### Keywords: models, nonlinear

### ** Examples

########### simulate Gaussian scalar-on-function data
n <- 500 ## number of observations
G <- 120 ## number of observations per functional covariate
set.seed(123) ## ensure reproducibility
z <- runif(n) ## scalar covariate
z <- z - mean(z)
s <- seq(0, 1, l=G) ## index of functional covariate
## generate functional covariate
if(require(splines)){
   x <- t(replicate(n, drop(bs(s, df = 5, int = TRUE) %*% runif(5, min = -1, max = 1))))
}else{
  x <- matrix(rnorm(n*G), ncol = G, nrow = n)
}
x <- scale(x, center = TRUE, scale = FALSE) ## center x per observation point

mu <- 2 + 0.5*z + (1/G*x) %*% sin(s*pi)*5 ## true functions for expectation
sigma <- exp(0.5*z - (1/G*x) %*% cos(s*pi)*2) ## for standard deviation

y <- rnorm(mean = mu, sd = sigma, n = n) ## draw respone y_i ~ N(mu_i, sigma_i)

## save data as list containing s as well 
dat_list <- list(y = y, z = z, x = I(x), s = s)

## model fit with noncyclic algorithm assuming Gaussian location scale model 
m_boost <- FDboostLSS(list(mu = y ~ bols(z, df = 2) + bsignal(x, s, df = 2, knots = 16), 
                           sigma = y ~ bols(z, df = 2) + bsignal(x, s, df = 2, knots = 16)), 
                           timeformula = NULL, data = dat_list, method = "noncyclic")
summary(m_boost)

## Not run: 
##D  if(require(gamboostLSS)){
##D   ## find optimal number of boosting iterations on a grid in 1:1000
##D   ## using 5-fold bootstrap
##D   ## takes some time, easy to parallelize on Linux
##D   set.seed(123) 
##D   cvr <- cvrisk(m_boost, folds = cv(model.weights(m_boost[[1]]), B = 5),
##D                 grid = 1:1000, trace = FALSE)
##D   ## use model at optimal stopping iterations 
##D   m_boost <- m_boost[mstop(cvr)] ## 832
##D    
##D   ## plot smooth effects of functional covariates for mu and sigma
##D   par(mfrow = c(1,2))
##D   plot(m_boost$mu, which = 2, ylim = c(0,5))
##D   lines(s, sin(s*pi)*5, col = 3, lwd = 2)
##D   plot(m_boost$sigma, which = 2, ylim = c(-2.5,2.5))
##D   lines(s, -cos(s*pi)*2, col = 3, lwd = 2)
##D  }
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("FDboostLSS", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("anisotropic_Kronecker")
### * anisotropic_Kronecker

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: anisotropic_Kronecker
### Title: Kronecker product or row tensor product of two base-learners
###   with anisotropic penalty
### Aliases: anisotropic_Kronecker %A% %A0% %Xa0%

### ** Examples

 
######## Example for anisotropic penalty  
data("viscosity", package = "FDboost") 
## set time-interval that should be modeled
interval <- "101"

## model time until "interval" and take log() of viscosity
end <- which(viscosity$timeAll == as.numeric(interval))
viscosity$vis <- log(viscosity$visAll[,1:end])
viscosity$time <- viscosity$timeAll[1:end]
# with(viscosity, funplot(time, vis, pch = 16, cex = 0.2))

## isotropic penalty, as timeformula is kroneckered to each effect using %O% 
## only for the smooth intercept %A0% is used, as 1-direction should not be penalized 
mod1 <- FDboost(vis ~ 1 + 
                bolsc(T_C, df = 1) + 
                bolsc(T_A, df = 1) + 
                bols(T_C, df = 1) %Xc% bols(T_A, df = 1),
                timeformula = ~ bbs(time, df = 3),
                numInt = "equal", family = QuantReg(),
                offset = NULL, offset_control = o_control(k_min = 9),
                data = viscosity, control=boost_control(mstop = 100, nu = 0.4))
## cf. the formula that is passed to mboost
mod1$formulaMboost

## anisotropic effects using %A0%, as lambda1 = 0 for all base-learners
## in this case using %A% gives the same model, but three lambdas are computed explicitly 
mod1a <- FDboost(vis ~ 1 + 
                bolsc(T_C, df = 1) %A0% bbs(time, df = 3) + 
                bolsc(T_A, df = 1) %A0% bbs(time, df = 3) + 
                bols(T_C, df = 1) %Xc% bols(T_A, df = 1) %A0% bbs(time, df = 3),
                timeformula = ~ bbs(time, df = 3),
                numInt = "equal", family = QuantReg(),
                offset = NULL, offset_control = o_control(k_min = 9),
                data = viscosity, control=boost_control(mstop = 100, nu = 0.4)) 
## cf. the formula that is passed to mboost
mod1a$formulaMboost

## alternative model specification by using a 0-matrix as penalty 
## only works for bolsc() as in bols() one cannot specify K 
## -> model without interaction term 
K0 <- matrix(0, ncol = 2, nrow = 2)
mod1k0 <- FDboost(vis ~ 1 + 
                 bolsc(T_C, df = 1, K = K0) +
                 bolsc(T_A, df = 1, K = K0), 
                 timeformula = ~ bbs(time, df = 3), 
                 numInt = "equal", family = QuantReg(), 
                 offset = NULL, offset_control = o_control(k_min = 9), 
                 data = viscosity, control=boost_control(mstop = 100, nu = 0.4))
## cf. the formula that is passed to mboost
mod1k0$formulaMboost
                
## optimize mstop for mod1, mod1a and mod1k0
## ...
                
## compare estimated coefficients
## Not run: 
##D par(mfrow=c(4, 2))
##D plot(mod1, which = 1)
##D plot(mod1a, which = 1)
##D plot(mod1, which = 2)
##D plot(mod1a, which = 2)
##D plot(mod1, which = 3)
##D plot(mod1a, which = 3)
##D funplot(mod1$yind, predict(mod1, which=4))
##D funplot(mod1$yind, predict(mod1a, which=4))
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("anisotropic_Kronecker", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("applyFolds")
### * applyFolds

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: applyFolds
### Title: Cross-Validation and Bootstrapping over Curves
### Aliases: applyFolds cvMa cvLong cvrisk.FDboost cvrisk.FDboost cvLong
###   cvMa

### ** Examples

Ytest <- matrix(rnorm(15), ncol = 3) # 5 trajectories, each with 3 observations 
Ylong <- as.vector(Ytest)
## 4-folds for bootstrap for the response in long format without integration weights
cvMa(ydim = c(5,3), type = "bootstrap", B = 4)  
cvLong(id = rep(1:5, times = 3), type = "bootstrap", B = 4)

if(require(fda)){
 ## load the data
 data("CanadianWeather", package = "fda")
 
 ## use data on a daily basis 
 canada <- with(CanadianWeather, 
                list(temp = t(dailyAv[ , , "Temperature.C"]),
                     l10precip = t(dailyAv[ , , "log10precip"]),
                     l10precip_mean = log(colMeans(dailyAv[ , , "Precipitation.mm"]), base = 10),
                     lat = coordinates[ , "N.latitude"],
                     lon = coordinates[ , "W.longitude"],
                     region = factor(region),
                     place = factor(place),
                     day = 1:365,  ## corresponds to t: evaluation points of the fun. response 
                     day_s = 1:365))  ## corresponds to s: evaluation points of the fun. covariate
 
## center temperature curves per day 
canada$tempRaw <- canada$temp
canada$temp <- scale(canada$temp, scale = FALSE) 
rownames(canada$temp) <- NULL ## delete row-names 
  
## fit the model  
mod <- FDboost(l10precip ~ 1 + bolsc(region, df = 4) + 
                 bsignal(temp, s = day_s, cyclic = TRUE, boundary.knots = c(0.5, 365.5)), 
               timeformula = ~ bbs(day, cyclic = TRUE, boundary.knots = c(0.5, 365.5)), 
               data = canada)
mod <- mod[75]

## Not run: 
##D   #### create folds for 3-fold bootstrap: one weight for each curve
##D   set.seed(123)
##D   folds_bs <- cv(weights = rep(1, mod$ydim[1]), type = "bootstrap", B = 3)
##D 
##D   ## compute out-of-bag risk on the 3 folds for 1 to 75 boosting iterations  
##D   cvr <- applyFolds(mod, folds = folds_bs, grid = 1:75)
##D 
##D   ## weights per observation point  
##D   folds_bs_long <- folds_bs[rep(1:nrow(folds_bs), times = mod$ydim[2]), ]
##D   attr(folds_bs_long, "type") <- "3-fold bootstrap"
##D   ## compute out-of-bag risk on the 3 folds for 1 to 75 boosting iterations  
##D   cvr3 <- cvrisk(mod, folds = folds_bs_long, grid = 1:75)
## End(Not run)

## Not run: 
##D   ## plot the out-of-bag risk
##D   par(mfrow = c(1,3))
##D   plot(cvr); legend("topright", lty=2, paste(mstop(cvr)))
##D   plot(cvr3); legend("topright", lty=2, paste(mstop(cvr3)))
## End(Not run)

}




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("applyFolds", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("bbsc")
### * bbsc

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: bbsc
### Title: Constrained Base-learners for Scalar Covariates
### Aliases: bbsc brandomc bolsc bolsc brandomc
### Keywords: models

### ** Examples

#### simulate data with functional response and scalar covariate (functional ANOVA)
n <- 60   ## number of cases
Gy <- 27  ## number of observation poionts per response curve 
dat <- list()
dat$t <- (1:Gy-1)^2/(Gy-1)^2
set.seed(123)
dat$z1 <- rep(c(-1, 1), length = n)
dat$z1_fac <- factor(dat$z1, levels = c(-1, 1), labels = c("1", "2"))
# dat$z1 <- runif(n)
# dat$z1 <- dat$z1 - mean(dat$z1)

# mean and standard deviation for the functional response 
mut <- matrix(2*sin(pi*dat$t), ncol = Gy, nrow = n, byrow = TRUE) + 
        outer(dat$z1, dat$t, function(z1, t) z1*cos(pi*t) ) # true linear predictor
sigma <- 0.1

# draw respone y_i(t) ~ N(mu_i(t), sigma)
dat$y <- apply(mut, 2, function(x) rnorm(mean = x, sd = sigma, n = n)) 

## fit function-on-scalar model with a linear effect of z1
m1 <- FDboost(y ~ 1 + bolsc(z1_fac, df = 1), timeformula = ~ bbs(t, df = 6), data = dat)

# look for optimal mSTOP using cvrisk() or validateFDboost()
 ## Not run: 
##D cvm <- cvrisk(m1, grid = 1:500)
##D m1[mstop(cvm)]
## End(Not run)
m1[200] # use 200 boosting iterations 

# plot true and estimated coefficients 
plot(dat$t, 2*sin(pi*dat$t), col = 2, type = "l", main = "intercept")
plot(m1, which = 1, lty = 2, add = TRUE)

plot(dat$t, 1*cos(pi*dat$t), col = 2, type = "l", main = "effect of z1")
lines(dat$t, -1*cos(pi*dat$t), col = 2, type = "l")
plot(m1, which = 2, lty = 2, col = 1, add = TRUE)





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("bbsc", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("bhistx")
### * bhistx

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: bhistx
### Title: Base-learners for Functional Covariates
### Aliases: bhistx
### Keywords: models

### ** Examples

if(require(refund)){
## simulate some data from a historical model
## the interaction effect is in this case not necessary
n <- 100
nygrid <- 35
data1 <- pffrSim(scenario = c("int", "ff"), limits = function(s,t){ s <= t }, 
                n = n, nygrid = nygrid)
data1$X1 <- scale(data1$X1, scale = FALSE) ## center functional covariate                  
dataList <- as.list(data1)
dataList$tvals <- attr(data1, "yindex")

## create the hmatrix-object
X1h <- with(dataList, hmatrix(time = rep(tvals, each = n), id = rep(1:n, nygrid), 
                             x = X1, argvals = attr(data1, "xindex"), 
                             timeLab = "tvals", idLab = "wideIndex", 
                             xLab = "myX", argvalsLab = "svals"))
dataList$X1h <- I(X1h)   
dataList$svals <- attr(data1, "xindex")
## add a factor variable 
dataList$zlong <- factor(gl(n = 2, k = n/2, length = n*nygrid), levels = 1:3)  
dataList$z <- factor(gl(n = 2, k = n/2, length = n), levels = 1:3)

## do the model fit with main effect of bhistx() and interaction of bhistx() and bolsc()
mod <- FDboost(Y ~ 1 + bhistx(x = X1h, df = 5, knots = 5) + 
               bhistx(x = X1h, df = 5, knots = 5) %X% bolsc(zlong), 
              timeformula = ~ bbs(tvals, knots = 10), data = dataList)
              
## alternative parameterization: interaction of bhistx() and bols()
mod <- FDboost(Y ~ 1 + bhistx(x = X1h, df = 5, knots = 5) %X% bols(zlong), 
              timeformula = ~ bbs(tvals, knots = 10), data = dataList)

## Not run: 
##D   # find the optimal mstop over 5-fold bootstrap (small example to reduce run time)
##D   cv <- cvrisk(mod, folds = cv(model.weights(mod), B = 5))
##D   mstop(cv)
##D   mod[mstop(cv)]
##D   
##D   appl1 <- applyFolds(mod, folds = cv(rep(1, length(unique(mod$id))), type = "bootstrap", B = 5))
##D 
##D  # plot(mod)
## End(Not run)
}




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("bhistx", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("bootstrapCI")
### * bootstrapCI

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: bootstrapCI
### Title: Function to compute bootstrap confidence intervals
### Aliases: bootstrapCI

### ** Examples

if(require(refund)){
#########
# model with linear functional effect, use bsignal()
# Y(t) = f(t) + \int X1(s)\beta(s,t)ds + eps
set.seed(2121)
data1 <- pffrSim(scenario = "ff", n = 40)
data1$X1 <- scale(data1$X1, scale = FALSE)
dat_list <- as.list(data1)
dat_list$t <- attr(data1, "yindex")
dat_list$s <- attr(data1, "xindex")

## model fit by FDboost 
m1 <- FDboost(Y ~ 1 + bsignal(x = X1, s = s, knots = 8, df = 3), 
              timeformula = ~ bbs(t, knots = 8), data = dat_list)

}
              
## Not run: 
##D              
##D # a short toy example with to few folds  
##D # and up to 200 boosting iterations 
##D bootCIs <- bootstrapCI(m1[200], B_inner = 2, B_outer = 5) 
##D 
##D # look at stopping iterations
##D bootCIs$mstops
##D 
##D # plot bootstrapped coefficient estimates
##D plot(bootCIs, ask = FALSE)
## End(Not run)

## now speed things up by defining the inner resampling
## function with parallelization based on mclapply (does not work on Windows)

my_inner_fun <- function(object){ 
cvrisk(object, folds = cvLong(id = object$id, weights = 
model.weights(object), 
B = 10 # 10-fold for inner resampling
), mc.cores = 10) # use ten cores
}

## Not run: 
##D bootCIs <- bootstrapCI(m1, resampling_fun_inner = my_inner_fun)
## End(Not run)

## We can also use the ... argument to parallelize the applyFolds
## function in the outer resampling 

## Not run: 
##D bootCIs <- bootstrapCI(m1, mc.cores = 30)
## End(Not run)

## Now let's parallelize the outer resampling and use 
## crossvalidation instead of bootstrap for the inner resampling

my_inner_fun <- function(object){ 
cvrisk(object, folds = cvLong(id = object$id, weights = 
model.weights(object), type = "kfold", # use CV
B = 10, # 10-fold for inner resampling
),
mc.cores = 10) # use ten cores
}

# use applyFolds for outer function to avoid messing up weights
my_outer_fun <- function(object, fun){
applyFolds(object = object,
folds = cv(rep(1, length(unique(object$id))), 
type = "bootstrap", B = 100), fun = fun,
mc.cores = 10) # parallelize on 10 cores
}

######## Example for scalar-on-function-regression with bsignal() 
data("fuelSubset", package = "FDboost")

## center the functional covariates per observed wavelength
fuelSubset$UVVIS <- scale(fuelSubset$UVVIS, scale = FALSE)
fuelSubset$NIR <- scale(fuelSubset$NIR, scale = FALSE)

## to make mboost:::df2lambda() happy (all design matrix entries < 10)
## reduce range of argvals to [0,1] to get smaller integration weights
fuelSubset$uvvis.lambda <- with(fuelSubset, (uvvis.lambda - min(uvvis.lambda)) /
(max(uvvis.lambda) - min(uvvis.lambda) ))
fuelSubset$nir.lambda <- with(fuelSubset, (nir.lambda - min(nir.lambda)) /
(max(nir.lambda) - min(nir.lambda) ))

## model fit with scalar response and two functional linear effects 
## include no intercept as all base-learners are centered around 0    

mod2 <- FDboost(heatan ~ bsignal(UVVIS, uvvis.lambda, knots = 40, df = 4, check.ident = FALSE) 
               + bsignal(NIR, nir.lambda, knots = 40, df=4, check.ident = FALSE), 
               timeformula = NULL, data = fuelSubset) 


## Not run: 
##D # takes some time, because of defaults: B_outer = 100, B_inner = 25
##D bootCIs <- bootstrapCI(mod2)
## End(Not run)

## run with a larger number of outer bootstrap samples
## and only 10-fold for validation of each outer fold
## WARNING: This may take very long!
## Not run: 
##D bootCIs <- bootstrapCI(mod2, B_outer = 1000, B_inner = 10)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("bootstrapCI", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("bsignal")
### * bsignal

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: bsignal
### Title: Base-learners for Functional Covariates
### Aliases: bsignal bconcurrent bhist bfpc bconcurrent bhist bfpc
### Keywords: models

### ** Examples

######## Example for scalar-on-function-regression with bsignal()  
data("fuelSubset", package = "FDboost")

## center the functional covariates per observed wavelength
fuelSubset$UVVIS <- scale(fuelSubset$UVVIS, scale = FALSE)
fuelSubset$NIR <- scale(fuelSubset$NIR, scale = FALSE)

## to make mboost:::df2lambda() happy (all design matrix entries < 10)
## reduce range of argvals to [0,1] to get smaller integration weights
fuelSubset$uvvis.lambda <- with(fuelSubset, (uvvis.lambda - min(uvvis.lambda)) /
                                  (max(uvvis.lambda) - min(uvvis.lambda) ))
fuelSubset$nir.lambda <- with(fuelSubset, (nir.lambda - min(nir.lambda)) /
                                (max(nir.lambda) - min(nir.lambda) ))

## model fit with scalar response and two functional linear effects 
## include no intercept 
## as all base-learners are centered around 0 
mod2 <- FDboost(heatan ~ bsignal(UVVIS, uvvis.lambda, knots = 40, df = 4, check.ident = FALSE) 
               + bsignal(NIR, nir.lambda, knots = 40, df=4, check.ident = FALSE), 
               timeformula = NULL, data = fuelSubset) 
summary(mod2) 
## plot(mod2)


###############################################
### data simulation like in manual of pffr::ff

if(require(refund)){

#########
# model with linear functional effect, use bsignal()
# Y(t) = f(t) + \int X1(s)\beta(s,t)ds + eps
set.seed(2121)
data1 <- pffrSim(scenario = "ff", n = 40)
data1$X1 <- scale(data1$X1, scale = FALSE)
dat_list <- as.list(data1)
dat_list$t <- attr(data1, "yindex")
dat_list$s <- attr(data1, "xindex")

## model fit by FDboost 
m1 <- FDboost(Y ~ 1 + bsignal(x = X1, s = s, knots = 5), 
              timeformula = ~ bbs(t, knots = 5), data = dat_list, 
              control = boost_control(mstop = 21))

## search optimal mSTOP
## Not run: 
##D   set.seed(123)
##D   cv <- validateFDboost(m1, grid = 1:100) # 21 iterations
## End(Not run)

## model fit by pffr
t <- attr(data1, "yindex")
s <- attr(data1, "xindex")
m1_pffr <- pffr(Y ~ ff(X1, xind = s), yind = t, data = data1)

## Not run: 
##D   par(mfrow = c(2, 2))
##D   plot(m1, which = 1); plot(m1, which = 2) 
##D   plot(m1_pffr, select = 1, shift = m1_pffr$coefficients["(Intercept)"]) 
##D   plot(m1_pffr, select = 2)
## End(Not run)


############################################
# model with functional historical effect, use bhist() 
# Y(t) = f(t)  + \int_0^t X1(s)\beta(s,t)ds + eps
set.seed(2121)
mylimits <- function(s, t){
  (s < t) | (s == t)
}
data2 <- pffrSim(scenario = "ff", n = 40, limits = mylimits)
data2$X1 <- scale(data2$X1, scale = FALSE)
dat2_list <- as.list(data2)
dat2_list$t <- attr(data2, "yindex")
dat2_list$s <- attr(data2, "xindex")

## model fit by FDboost 
m2 <- FDboost(Y ~ 1 + bhist(x = X1, s = s, time = t, knots = 5), 
              timeformula = ~ bbs(t, knots = 5), data = dat2_list, 
              control = boost_control(mstop = 40))
              
## search optimal mSTOP
## Not run: 
##D   set.seed(123)
##D   cv2 <- validateFDboost(m2, grid = 1:100) # 40 iterations
## End(Not run)               

## model fit by pffr
t <- attr(data2, "yindex")
s <- attr(data2, "xindex")
m2_pffr <- pffr(Y ~ ff(X1, xind = s, limits = "s<=t"), yind = t, data = data2)

## Not run: 
##D par(mfrow = c(2, 2))
##D plot(m2, which = 1); plot(m2, which = 2)
##D ## plot of smooth intercept does not contain m1_pffr$coefficients["(Intercept)"]
##D plot(m2_pffr, select = 1, shift = m2_pffr$coefficients["(Intercept)"]) 
##D plot(m2_pffr, select = 2) 
##D 
## End(Not run)


}





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("bsignal", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("emotion")
### * emotion

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: emotion
### Title: EEG and EMG recordings in a computerised gambling study
### Aliases: emotion
### Keywords: datasets

### ** Examples

data("emotion", package = "FDboost")

# fit function-on-scalar model with random effect and power effect
fos_random_power <- FDboost(EMG ~ 1 + brandomc(subject, df = 2)
                            + bolsc(power, df = 2),
                            timeformula = ~ bbs(t, df = 3),
                            data = emotion)
## Not run: 
##D                             
##D # fit function-on-function model with intercept and historical EEG effect
##D # where limits specifies the used lag between EMG and EEG signal
##D fof_historical <- FDboost(EMG ~ 1 + bhist(EEG, s = s, time = t,
##D                           limits = function(s,t) s < t - 3),
##D                           timeformula = ~ bbs(t, df = 3), data = emotion,
##D                           control = boost_control(mstop = 200))                            
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("emotion", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fuelSubset")
### * fuelSubset

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fuelSubset
### Title: Spectral data of fossil fuels
### Aliases: fuelSubset
### Keywords: datasets

### ** Examples


    data("fuelSubset", package = "FDboost")
    
    ## center the functional covariates per observed wavelength
    fuelSubset$UVVIS <- scale(fuelSubset$UVVIS, scale = FALSE)
    fuelSubset$NIR <- scale(fuelSubset$NIR, scale = FALSE)

    ## to make mboost::df2lambda() happy (all design matrix entries < 10)
    ## reduce range of argvals to [0,1] to get smaller integration weights
    fuelSubset$uvvis.lambda <- with(fuelSubset, (uvvis.lambda - min(uvvis.lambda)) /
                                          (max(uvvis.lambda) - min(uvvis.lambda) ))
    fuelSubset$nir.lambda <- with(fuelSubset, (nir.lambda - min(nir.lambda)) /
                                          (max(nir.lambda) - min(nir.lambda) ))


    ### fit mean regression model with 100 boosting iterations,
    ### step-length 0.1 and
    mod <- FDboost(heatan ~ bsignal(UVVIS, uvvis.lambda, knots=40, df=4, check.ident=FALSE)
                   + bsignal(NIR, nir.lambda, knots=40, df=4, check.ident=FALSE),
                   timeformula = NULL, data = fuelSubset)
    summary(mod)
    ## plot(mod)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fuelSubset", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("funplot")
### * funplot

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: funplot
### Title: Plot functional data with linear interpolation of missing values
### Aliases: funplot

### ** Examples

## Not run: 
##D ### examples for regular data in wide format
##D data(viscosity)
##D with(viscosity, funplot(timeAll, visAll, pch=20))
##D if(require(fda)){
##D   with(fda::growth, funplot(age, t(hgtm)))
##D }
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("funplot", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("grapes-Xc-grapes")
### * grapes-Xc-grapes

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: %Xc%
### Title: Constrained row tensor product
### Aliases: %Xc%

### ** Examples

 
######## Example for function-on-scalar-regression with interaction effect of two scalar covariates 
data("viscosity", package = "FDboost") 
## set time-interval that should be modeled
interval <- "101"

## model time until "interval" and take log() of viscosity
end <- which(viscosity$timeAll == as.numeric(interval))
viscosity$vis <- log(viscosity$visAll[,1:end])
viscosity$time <- viscosity$timeAll[1:end]
# with(viscosity, funplot(time, vis, pch = 16, cex = 0.2))

## fit model with interaction that is centered around the intercept 
## and the two main effects 
mod1 <- FDboost(vis ~ 1 + bolsc(T_C, df=1) + bolsc(T_A, df=1) + 
                bols(T_C, df=1) %Xc% bols(T_A, df=1),
                timeformula = ~bbs(time, df=6),
                numInt = "equal", family = QuantReg(),
                offset = NULL, offset_control = o_control(k_min = 9),
                data = viscosity, control=boost_control(mstop = 100, nu = 0.4))
                
## check centering around intercept
colMeans(predict(mod1, which = 4))

## check centering around main effects
colMeans(predict(mod1, which = 4)[viscosity$T_A == "low", ])
colMeans(predict(mod1, which = 4)[viscosity$T_A == "high", ])
colMeans(predict(mod1, which = 4)[viscosity$T_C == "low", ])
colMeans(predict(mod1, which = 4)[viscosity$T_C == "low", ])

## find optimal mstop using cvrsik() or validateFDboost()
## ... 

## look at interaction effect in one plot
# funplot(mod1$yind, predict(mod1, which=4))




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("grapes-Xc-grapes", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("hmatrix")
### * hmatrix

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: hmatrix
### Title: A S3 class for univariate functional data on a common grid
### Aliases: hmatrix

### ** Examples

## Example for a hmatrix object
t1 <- rep((1:5)/2, each = 3)
id1 <- rep(1:3, 5)
x1 <- matrix(1:15, ncol = 5) 
s1 <- (1:5)/2 
myhmatrix <- hmatrix(time = t1, id = id1, x = x1, argvals = s1, 
                     timeLab = "t1", argvalsLab = "s1", xLab = "test")

# extract with [ keeps attributes 
# select observations of subjects 2 and 3
myhmatrixSub <- myhmatrix[id1 %in% c(2, 3), ]  
str(myhmatrixSub)
getX(myhmatrixSub)
getX(myhmatrix)

# get time
myhmatrix[ , 1] # as column matrix as drop = FALSE
getTime(myhmatrix) # as vector

# get id
myhmatrix[ , 2] # as column matrix as drop = FALSE
getId(myhmatrix) # as vector

# subset hmatrix on the basis of an index, which is defined on the curve level
reweightData(data = list(hmat = myhmatrix), vars = "hmat", index = c(1, 1, 2))
# this keeps only the unique x values in attr(,'x') but multiplies the corresponding
# ids and times in the time id matrix 
# for bhistx baselearner, there may be an additional id variable for the tensor product
newdat <- reweightData(data = list(hmat = myhmatrix, 
  repIDx = rep(1:nrow(attr(myhmatrix,'x')), length(attr(myhmatrix,"argvals")))), 
  vars = "hmat", index = c(1,1,2), idvars="repIDx")
length(newdat$repIDx) 

## use hmatrix within a data.frame
mydat <- data.frame(I(myhmatrix), z=rnorm(3)[id1])
str(mydat)
str(mydat[id1 %in% c(2, 3), ])
str(myhmatrix[id1 %in% c(2, 3), ])
 



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("hmatrix", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("integrationWeights")
### * integrationWeights

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: integrationWeights
### Title: Functions to compute integration weights
### Aliases: integrationWeights integrationWeightsLeft
###   integrationWeightsLeft

### ** Examples

## Example for trapezoidal integration weights
xind0 <- seq(0,1,l = 5)
xind <- c(0, 0.1, 0.3, 0.7, 1)
X1 <- matrix(xind^2, ncol = length(xind0), nrow = 2)

# Regualar observation points
integrationWeights(X1, xind0)
# Irregular observation points
integrationWeights(X1, xind)

# with missing value
X1[1,2] <- NA
integrationWeights(X1, xind0)
integrationWeights(X1, xind)

## Example for left integration weights
xind0 <- seq(0,1,l = 5)
xind <- c(0, 0.1, 0.3, 0.7, 1)
X1 <- matrix(xind^2, ncol = length(xind0), nrow = 2)

# Regular observation points
integrationWeightsLeft(X1, xind0, leftWeight = "mean") 
integrationWeightsLeft(X1, xind0, leftWeight = "first") 
integrationWeightsLeft(X1, xind0, leftWeight = "zero")

# Irregular observation points
integrationWeightsLeft(X1, xind, leftWeight = "mean") 
integrationWeightsLeft(X1, xind, leftWeight = "first") 
integrationWeightsLeft(X1, xind, leftWeight = "zero")

# obervation points that do not start with 0
xind2 <- xind + 0.5
integrationWeightsLeft(X1, xind2, leftWeight = "zero")
 



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("integrationWeights", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("reweightData")
### * reweightData

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: reweightData
### Title: Function to Reweight Data
### Aliases: reweightData

### ** Examples

## load data
data("viscosity", package = "FDboost")
interval <- "101"
end <- which(viscosity$timeAll == as.numeric(interval))
viscosity$vis <- log(viscosity$visAll[ , 1:end])
viscosity$time <- viscosity$timeAll[1:end]

## what does data look like
str(viscosity)

## do some reweighting
# correct weights
str(reweightData(viscosity, vars=c("vis", "T_C", "T_A", "rspeed", "mflow"), 
    argvals = "time", weights = c(0, 32, 32, rep(0, 61))))

str(visNew <- reweightData(viscosity, vars=c("vis", "T_C", "T_A", "rspeed", "mflow"), 
    argvals = "time", weights = c(0, 32, 32, rep(0, 61))))
# check the result
# visNew$vis[1:5, 1:5] ## image(visNew$vis)

# incorrect weights
str(reweightData(viscosity, vars=c("vis", "T_C", "T_A", "rspeed", "mflow"), 
    argvals = "time", weights = sample(1:64, replace = TRUE)), 1)

# supply meaningful index
str(visNew <- reweightData(viscosity, vars = c("vis", "T_C", "T_A", "rspeed", "mflow"), 
              argvals = "time", index = rep(1:32, each = 2)))
# check the result
# visNew$vis[1:5, 1:5]

# errors
if(FALSE){
   reweightData(viscosity, argvals = "")
   reweightData(viscosity, argvals = "covThatDoesntExist", index = rep(1,64))
   }
   



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("reweightData", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("stabsel.FDboost")
### * stabsel.FDboost

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: stabsel.FDboost
### Title: Stability Selection
### Aliases: stabsel.FDboost

### ** Examples

######## Example for function-on-scalar-regression
data("viscosity", package = "FDboost")
## set time-interval that should be modeled
interval <- "101"

## model time until "interval" and take log() of viscosity
end <- which(viscosity$timeAll == as.numeric(interval))
viscosity$vis <- log(viscosity$visAll[,1:end])
viscosity$time <- viscosity$timeAll[1:end]
# with(viscosity, funplot(time, vis, pch = 16, cex = 0.2))

## fit a model cotaining all main effects 
modAll <- FDboost(vis ~ 1 
          + bolsc(T_C, df=1) %A0% bbs(time, df=5) 
          + bolsc(T_A, df=1) %A0% bbs(time, df=5)
          + bolsc(T_B, df=1) %A0% bbs(time, df=5)
          + bolsc(rspeed, df=1) %A0% bbs(time, df=5)
          + bolsc(mflow, df=1) %A0% bbs(time, df=5), 
       timeformula = ~bbs(time, df=5), 
       numInt = "Riemann", family = QuantReg(), 
       offset = NULL, offset_control = o_control(k_min = 10),
       data = viscosity, 
       control = boost_control(mstop = 100, nu = 0.2))


## create folds for stability selection  
## only 5 folds for a fast example, usually use 50 folds 
set.seed(1911)
folds <- cvLong(modAll$id, weights = rep(1, l = length(modAll$id)), 
                type = "subsampling", B = 5) 
    
## Not run: 
##D         
##D ## stability selection with refit of the smooth intercept 
##D stabsel_parameters(q = 3, PFER = 1, p = 6, sampling.type = "SS")
##D sel1 <- stabsel(modAll, q = 3, PFER = 1, folds = folds, grid = 1:200, sampling.type = "SS")
##D sel1
##D 
##D ## stability selection without refit of the smooth intercept 
##D sel2 <- stabsel(modAll, refitSmoothOffset = FALSE, q = 3, PFER = 1, 
##D                 folds = folds, grid = 1:200, sampling.type = "SS")
##D sel2
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("stabsel.FDboost", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("subset_hmatrix")
### * subset_hmatrix

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: subset_hmatrix
### Title: Subsets hmatrix according to an index
### Aliases: subset_hmatrix

### ** Examples

t1 <- rep((1:5)/2, each = 3)
id1 <- rep(1:3, 5)
x1 <- matrix(1:15, ncol = 5) 
s1 <- (1:5)/2 
hmat <- hmatrix(time = t1, id = id1, x = x1, argvals = s1, timeLab = "t1", 
                argvalsLab = "s1", xLab = "test")

index1 <- c(1, 1, 3)
index2 <- c(2, 3, 3)
resMat <- subset_hmatrix(hmat, index = index1)
try(resMat2 <- subset_hmatrix(resMat, index = index2))
resMat <- subset_hmatrix(hmat, index = index1, compress = FALSE)
try(resMat2 <- subset_hmatrix(resMat, index = index2))




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("subset_hmatrix", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("truncateTime")
### * truncateTime

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: truncateTime
### Title: Function to truncate time in functional data
### Aliases: truncateTime

### ** Examples

if(require(fda)){
  dat <- fda::growth
  dat$hgtm <- t(dat$hgtm[,1:10])
  dat$hgtf <- t(dat$hgtf[,1:10])
  
  ## only use time-points 1:16 of variable age
  datTr <- truncateTime(funVar=c("hgtm","hgtf"), time="age", newtime=1:16, data=dat)
  
  ## Not run: 
##D   par(mfrow=c(1,2))
##D   with(dat, funplot(age, hgtm, main="Original data"))
##D   with(datTr, funplot(age, hgtm, main="Yearly data"))
##D   par(mfrow=c(1,1))   
##D   
## End(Not run)
}



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("truncateTime", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("update.FDboost")
### * update.FDboost

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: update.FDboost
### Title: Function to update FDboost objects
### Aliases: update.FDboost

### ** Examples

######## Example from \code{?FDboost}
data("viscosity", package = "FDboost") 
## set time-interval that should be modeled
interval <- "101"

## model time until "interval" and take log() of viscosity
end <- which(viscosity$timeAll == as.numeric(interval))
viscosity$vis <- log(viscosity$visAll[,1:end])
viscosity$time <- viscosity$timeAll[1:end]
# with(viscosity, funplot(time, vis, pch = 16, cex = 0.2))

mod1 <- FDboost(vis ~ 1 + bolsc(T_C, df = 2) + bolsc(T_A, df = 2),
               timeformula = ~ bbs(time, df = 4),
               numInt = "equal", family = QuantReg(),
               offset = NULL, offset_control = o_control(k_min = 9),
               data = viscosity, control=boost_control(mstop = 10, nu = 0.4))
               
# update nu
mod2 <- update(mod1, control=boost_control(nu = 1)) # mstop will stay the same
# update mstop
mod3 <- update(mod2, control=boost_control(mstop = 100)) # nu=1 does not get changed
mod4 <- update(mod1, formula = vis ~ 1 + bolsc(T_C, df = 2)) # drop one term



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("update.FDboost", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("validateFDboost")
### * validateFDboost

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: validateFDboost
### Title: Cross-Validation and Bootstrapping over Curves
### Aliases: validateFDboost

### ** Examples

## Not run: 
##D if(require(fda)){
##D  ## load the data
##D  data("CanadianWeather", package = "fda")
##D  
##D  ## use data on a daily basis 
##D  canada <- with(CanadianWeather, 
##D                 list(temp = t(dailyAv[ , , "Temperature.C"]),
##D                      l10precip = t(dailyAv[ , , "log10precip"]),
##D                      l10precip_mean = log(colMeans(dailyAv[ , , "Precipitation.mm"]), base = 10),
##D                      lat = coordinates[ , "N.latitude"],
##D                      lon = coordinates[ , "W.longitude"],
##D                      region = factor(region),
##D                      place = factor(place),
##D                      day = 1:365,  ## corresponds to t: evaluation points of the fun. response 
##D                      day_s = 1:365))  ## corresponds to s: evaluation points of the fun. covariate
##D  
##D ## center temperature curves per day 
##D canada$tempRaw <- canada$temp
##D canada$temp <- scale(canada$temp, scale = FALSE) 
##D rownames(canada$temp) <- NULL ## delete row-names 
##D   
##D ## fit the model  
##D mod <- FDboost(l10precip ~ 1 + bolsc(region, df = 4) + 
##D                  bsignal(temp, s = day_s, cyclic = TRUE, boundary.knots = c(0.5, 365.5)), 
##D                timeformula = ~ bbs(day, cyclic = TRUE, boundary.knots = c(0.5, 365.5)), 
##D                data = canada)
##D mod <- mod[75]
##D 
##D   #### create folds for 3-fold bootstrap: one weight for each curve
##D   set.seed(123)
##D   folds_bs <- cv(weights = rep(1, mod$ydim[1]), type = "bootstrap", B = 3)
##D 
##D   ## compute out-of-bag risk on the 3 folds for 1 to 75 boosting iterations  
##D   cvr <- applyFolds(mod, folds = folds_bs, grid = 1:75)
##D 
##D   ## compute out-of-bag risk and coefficient estimates on folds  
##D   cvr2 <- validateFDboost(mod, folds = folds_bs, grid = 1:75)
##D 
##D   ## weights per observation point  
##D   folds_bs_long <- folds_bs[rep(1:nrow(folds_bs), times = mod$ydim[2]), ]
##D   attr(folds_bs_long, "type") <- "3-fold bootstrap"
##D   ## compute out-of-bag risk on the 3 folds for 1 to 75 boosting iterations  
##D   cvr3 <- cvrisk(mod, folds = folds_bs_long, grid = 1:75)
##D 
##D   ## plot the out-of-bag risk
##D   par(mfrow = c(1,3))
##D   plot(cvr); legend("topright", lty=2, paste(mstop(cvr)))
##D   plot(cvr2)
##D   plot(cvr3); legend("topright", lty=2, paste(mstop(cvr3)))
##D 
##D   ## plot the estimated coefficients per fold
##D   ## more meaningful for higher number of folds, e.g., B = 100 
##D   par(mfrow = c(2,2))
##D   plotPredCoef(cvr2, terms = FALSE, which = 2)
##D   plotPredCoef(cvr2, terms = FALSE, which = 3)
##D   
##D   ## compute out-of-bag risk and predictions for leaving-one-curve-out cross-validation
##D   cvr_jackknife <- validateFDboost(mod, folds = cvLong(unique(mod$id), 
##D                                    type = "curves"), grid = 1:75)
##D   plot(cvr_jackknife)
##D   ## plot oob predictions per fold for 3rd effect 
##D   plotPredCoef(cvr_jackknife, which = 3) 
##D   ## plot coefficients per fold for 2nd effect
##D   plotPredCoef(cvr_jackknife, which = 2, terms = FALSE)
##D 
##D }
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("validateFDboost", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("viscosity")
### * viscosity

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: viscosity
### Title: Viscosity of resin over time
### Aliases: viscosity
### Keywords: datasets

### ** Examples


 data("viscosity", package = "FDboost")
 ## set time-interval that should be modeled
 interval <- "101"

 ## model time until "interval" and take log() of viscosity
 end <- which(viscosity$timeAll==as.numeric(interval))
 viscosity$vis <- log(viscosity$visAll[,1:end])
 viscosity$time <- viscosity$timeAll[1:end]
 # with(viscosity, funplot(time, vis, pch=16, cex=0.2))

 ## fit median regression model with 100 boosting iterations,
 ## step-length 0.4 and smooth time-specific offset
 ## the factors are in effect coding -1, 1 for the levels
 mod <- FDboost(vis ~ 1 + bols(T_C, contrasts.arg = "contr.sum", intercept=FALSE)
                + bols(T_A, contrasts.arg = "contr.sum", intercept=FALSE),
                timeformula=~bbs(time, lambda=100),
                numInt="equal", family=QuantReg(),
                offset=NULL, offset_control = o_control(k_min = 9),
                data=viscosity, control=boost_control(mstop = 100, nu = 0.4))
 summary(mod)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("viscosity", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
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
