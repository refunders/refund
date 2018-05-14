

library(FDboost)

################################################################
######### simulate some data 

if(require(refund)){
  
  ## simulate a small data set 
  set.seed(230)
  pffr_data <- pffrSim(n = 25, nxgrid = 21, nygrid = 19)
  pffr_data$X1 <- scale(pffr_data$X1, scale = FALSE)
  
  dat <- as.list(pffr_data)
  dat$tvals <- attr(pffr_data, "yindex")
  dat$svals <- attr(pffr_data, "xindex")
  
  dat$Y_scalar <- dat$Y[ , 10]
  
  dat$Y_long <- c(dat$Y)
  dat$tvals_long <- rep(dat$tvals, each = nrow(dat$Y))
  dat$id_long <- rep(1:nrow(dat$Y), ncol(dat$Y))
  
  
  ################################################################
  ######### model fit 
  
  ## response matrix for response observed on one common grid 
  m <- FDboost(Y ~ 1 + bhist(X1, svals, tvals, knots = 6, df = 12) 
               + bsignal(X1, svals, knots = 6, df = 4)
               + bbsc(xsmoo, knots = 6, df = 4) 
               + bolsc(xte1, df = 4)
               + brandomc(xte2, df = 4), 
               timeformula = ~ bbs(tvals, knots = 9, df = 3, differences = 1), 
               control = boost_control(mstop = 10), data = dat)
  
  ## response in long format
  ml <- FDboost(Y_long ~ 1 + bhist(X1, svals, tvals_long, knots = 6, df = 12) 
                + bsignal(X1, svals, knots = 6, df = 4)
                + bbsc(xsmoo, knots = 6, df = 4) 
                + bolsc(xte1, df = 4)
                + brandomc(xte2, df = 4), 
                timeformula = ~ bbs(tvals_long, knots = 8, df = 3, differences = 1), 
                id = ~ id_long, 
                offset_control = o_control(k_min = 10), 
                control = boost_control(mstop = 10), data = dat)
  
  ## scalar response 
  ms <- FDboost(Y_scalar ~ 1 + bsignal(X1, svals, knots = 6, df = 2)
                + bbs(xsmoo, knots = 6, df = 2, differences = 1) 
                + bols(xte1, df = 2) 
                + bols(xte2, df = 2), 
                timeformula = NULL, 
                control = boost_control(mstop = 50), data = dat)
  
  ## GAMLSS with functional response 
  mlss <- FDboostLSS(Y ~ 1 + bsignal(X1, svals, knots = 6, df = 4), 
                     timeformula = ~ bbs(tvals, knots = 9, df = 3, differences = 1), 
                     control = boost_control(mstop = 10), data = dat)
  
  
  ################################################################
  ######### test some methods and utility functions 
  
  ## test plot()
  par(mfrow = c(4,4))
  plot(m, ask = FALSE)
  plot(ml, ask = FALSE)
  plot(ms, ask = FALSE)
  plot(mlss$mu, ask = FALSE); plot(mlss$sigma, ask = FALSE)
  
  ## test applyFolds()
  set.seed(123)
  applyFolds(m, folds = cv(rep(1, length(unique(m$id))), B = 2), grid = 0:5)
  applyFolds(ml, folds = cv(rep(1, length(unique(ml$id))), B = 2), grid = 0:5)
  applyFolds(ms, folds = cv(rep(1, length(unique(ms$id))), B = 2), grid = 0:5)
  
  
}

