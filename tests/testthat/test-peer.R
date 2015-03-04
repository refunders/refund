context("Testing peer")
library(refundDevel)

test_that("peer with D2 penalty", {
  skip_on_cran()

  data(DTI)

  ## Extract values for arguments for peer() from given data
  cca = DTI$cca[which(DTI$case == 1),]
  DTI = DTI[which(DTI$case == 1),]

  ##1.1 Fit the model
  expect_message(peer(Y=DTI$pasat, funcs = cca, pentype='D2', se=TRUE), "The fit is successful.")
})

test_that("peer with structured penalty works", {
 skip_on_cran()

 data(PEER.Sim)

 ## Extract values for arguments for peer() from given data
 PEER.Sim1<- subset(PEER.Sim, t==0)
 K<- 100
 W<- PEER.Sim1[,c(3:(K+2))]
 Y<- PEER.Sim1[,K+3]

 ##Load Q matrix containing structural information
 data(Q)

 ##2.1 Fit the model
 Fit1<- peer(Y=Y, funcs=W, pentype='Decomp', Q=Q, se=TRUE)
 expect_is(Fit1, "peer")
 plot(Fit1)
})
