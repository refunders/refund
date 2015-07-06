context("Testing peer")
library(refundDevel)

test_that("peer with D2 penalty", {
  skip_on_cran()

  data(DTI)

  ## Extract values for arguments for peer() from given data
  cca = DTI$cca[which(DTI$case == 1),]
  DTI = DTI[which(DTI$case == 1),]

  ##1.1 Fit the model
  fit.D2 <- pfr(pasat ~ peer(cca, pentype="D"), data=DTI)
  expect_is(fit.D2, "pfr")
})

test_that("peer with structured penalty works", {
 skip_on_cran()
 
 data(PEER.Sim)
 data(Q)
 PEER.Sim1<- subset(PEER.Sim, t==0)
 
 # Setting k to max possible value
 fit.decomp <- pfr(Y ~ peer(W, pentype="Decomp", Q=Q, k=99), data=PEER.Sim1)
 expect_is(fit.decomp, "pfr")
})
