context("Testing rlrt.pfr")
library(refundDevel)

test_that("basis rlrt tests are working", {
  skip_on_cran()

  data(DTI2)
  O  <- DTI2$pasat ## PASAT outcome
  id <- DTI2$id    ## subject id
  W1 <- DTI2$cca   ## Corpus Callosum
  W2 <- DTI2$rcst  ## Right corticospinal
  V  <- DTI2$visit ## visit

 ## prep scalar covariate
 visit.1.rest <- matrix(as.numeric(V > 1), ncol=1)
 covar.in <- visit.1.rest


  ## fit two univariate models, then one model with both functional predictors
  suppressWarnings({
    pfr.obj.t1 <- pfr(Y = O, covariates=covar.in, funcs = list(W1),     subj = id, kz = 10, kb = 50)
    pfr.obj.t2 <- pfr(Y = O, covariates=covar.in, funcs = list(W2),     subj = id, kz = 10, kb = 50)
    pfr.obj.t3 <- pfr(Y = O, covariates=covar.in, funcs = list(W1, W2), subj = id, kz = 10, kb = 50)
  })

 ## do some testing
 t1 <- rlrt.pfr(pfr.obj.t1, "constancy")
 t2 <- rlrt.pfr(pfr.obj.t2, "constancy")
 t3 <- rlrt.pfr(pfr.obj.t3, "inclusion")

  expect_is(t1, "list")
  expect_is(t2, "list")
  expect_is(t3, "list")
})

test_that("tests work with subj = NULL", {
  skip_on_cran()

  suppressWarnings({
    pfr.obj.t1 <- pfr(Y = O, covariates=covar.in, funcs = list(W1),     subj = NULL, kz = 10, kb = 50)
    pfr.obj.t2 <- pfr(Y = O, covariates=covar.in, funcs = list(W2),     subj = NULL, kz = 10, kb = 50)
    pfr.obj.t3 <- pfr(Y = O, covariates=covar.in, funcs = list(W1, W2), subj = NULL, kz = 10, kb = 50)
  })
  t1 <- rlrt.pfr(pfr.obj.t1, "constancy")
  t2 <- rlrt.pfr(pfr.obj.t2, "constancy")
  t3 <- rlrt.pfr(pfr.obj.t3, "inclusion")

  expect_is(t1, "list")
  expect_is(t2, "list")
  expect_is(t3, "list")
})

test_that("lpeer with pentype='DECOMP' works", {
 skip_on_cran()

 ##Load Data
 data(PEER.Sim)

 ## Extract values for arguments for lpeer() from given data
 K<- 100
 W<- PEER.Sim[,c(3:(K+2))]
 Y<- PEER.Sim[,K+3]
 t<- PEER.Sim[,2]
 id<- PEER.Sim[,1]

 ##Load Q matrix containing structural information
 data(Q)

 ##2.1 Fit the model with two component function
 ##    gamma(t,s)=gamm0(s) + t*gamma1(s)
 expect_message(lpeer(Y=Y, subj=id, t=t, covariates=cbind(t), funcs=W,
              pentype='DECOMP', f_t=cbind(1,t), Q=Q, se=TRUE),
     "The fit is successful.")

 ## Fit1$Beta
 ## plot(Fit1)

 ##2.2 Fit the model with three component function
 ##    gamma(t,s)=gamm0(s) + t*gamma1(s) + t^2*gamma1(s)
 Fit2<- lpeer(Y=Y, subj=id, t=t, covariates=cbind(t), funcs=W,
 		     pentype='DECOMP', f_t=cbind(1,t, t^2), Q=Q, se=TRUE)

 Fit2$Beta
 plot(Fit2)

 ##2.3 Fit the model with two component function with different penalties
 ##    gamma(t,s)=gamm0(s) + t*gamma1(s)
 Q1<- cbind(Q, Q)
 expect_message(lpeer(Y=Y, subj=id, t=t, covariates=cbind(t), comm.pen=FALSE, funcs=W,
 		     pentype='DECOMP', f_t=cbind(1,t), Q=Q1, se=TRUE), "The fit is successful.")
})

test_that("lpeer with pentype='USER' works", {
 skip_on_cran()

 ##Load Data
 data(PEER.Sim)

 ## Extract values for arguments for lpeer() from given data
 K<- 100
 W<- PEER.Sim[,c(3:(K+2))]
 Y<- PEER.Sim[,K+3]
 t<- PEER.Sim[,2]
 id<- PEER.Sim[,1]

 ##Load Q matrix containing structural information
 data(Q)

 ##2.4 Fit the model with two component function with user defined penalties
 ##    gamma(t,s)=gamm0(s) + t*gamma1(s)
 phia<- 10^3
 P_Q <- t(Q)%*%solve(Q%*%t(Q))%*% Q
 L<- phia*(diag(K)- P_Q) + 1*P_Q
 expect_message(lpeer(Y=Y, subj=id, t=t, covariates=cbind(t), funcs=W,
 		     pentype='USER', f_t=cbind(1,t), L=L, se=TRUE), "The fit is successful.")

 L1<- adiag(L, L)
 expect_message(lpeer(Y=Y, subj=id, t=t, covariates=cbind(t), comm.pen=FALSE, funcs=W,
              pentype='USER', f_t=cbind(1,t), L=L1, se=TRUE), "The fit is successful.")
})
