context("Testing peer")
library(refundDevel)

test_that("lpeer with ridge penalty works", {
  skip_on_cran()

 ##Load Data
 data(DTI)

 ## Extract values for arguments for lpeer() from given data
 cca = DTI$cca[which(DTI$case == 1),]
 DTI = DTI[which(DTI$case == 1),]

 ##1.1 Fit the model with single component function
 ##    gamma(t,s)=gamm0(s)
 t<- DTI$visit
 expect_message(lpeer(Y=DTI$pasat, t=t, subj=DTI$ID, funcs = cca),
                 "The fit is successful.")

 ## plot(fit.cca.lpeer1)

 ##1.2 Fit the model with two component function
 ##    gamma(t,s)=gamm0(s) + t*gamma1(s)
 fit.cca.lpeer2 = lpeer(Y=DTI$pasat, t=t, subj=DTI$ID, funcs = cca,
                       f_t=t, se=TRUE)
 plot(fit.cca.lpeer2)
})

test_that("lpeer with pentype='DECOMP' works", {
 skip_on_cran()

 ##Load Data
 data(PEER.Sim)

 ## Extract values for arguments for lpeer() from given data
 K<- 100
 W <- PEER.Sim[, "W"]
 Y <- PEER.Sim[, "Y"]
 t <- PEER.Sim[, "t"]
 id <- PEER.Sim[, "id"]

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
 K <- 100
 W <- PEER.Sim[, "W"]
 Y <- PEER.Sim[, "Y"]
 t <- PEER.Sim[, "t"]
 id <- PEER.Sim[, "id"]

 ##Load Q matrix containing structural information
 data(Q)

 ##2.4 Fit the model with two component function with user defined penalties
 ##    gamma(t,s)=gamm0(s) + t*gamma1(s)
 phia<- 10^3
 P_Q <- t(Q)%*%solve(Q%*%t(Q))%*% Q
 L<- phia*(diag(K)- P_Q) + 1*P_Q
 expect_message(lpeer(Y=Y, subj=id, t=t, covariates=cbind(t), funcs=W,
 		     pentype='USER', f_t=cbind(1,t), L=L, se=TRUE), "The fit is successful.")

 L1<- magic::adiag(L, L)
 expect_message(lpeer(Y=Y, subj=id, t=t, covariates=cbind(t), comm.pen=FALSE, funcs=W,
              pentype='USER', f_t=cbind(1,t), L=L1, se=TRUE), "The fit is successful.")
})
