context("Testing pfr's peer()")

test_that("peer with D2 penalty", {
  skip_on_cran()

  data(DTI)
  DTI = DTI[which(DTI$case == 1),]
  fit.D2 <- pfr(pasat ~ peer(cca, pentype="D"), data=DTI)
  expect_is(fit.D2, "pfr")
})

test_that("peer with structured penalty works", {
 skip_on_cran()

 data(PEER.Sim, Q)

 # Setting k to max possible value
 fit.decomp <- pfr(Y ~ peer(W, pentype="Decomp", Q=Q, k=99),
   data=subset(PEER.Sim, t==0))
 expect_is(fit.decomp, "pfr")
})
