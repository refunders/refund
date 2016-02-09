context("Testing fpca.xxx")
library(refund)

set.seed(12212)
n <- 100
ngrid <- 40
t <- seq(0, 1, l=ngrid)
efcts <- poly(t, 2)
Y <- outer(2 * drop(scale(rnorm(n))), efcts[, 1]) +
  outer(drop(scale(rnorm(n))), efcts[, 2])

test_that("all fpca functions agree on toy example", {
  sc <- fpca.sc(Y)
  face <- fpca.face(Y)
  ssvd <- fpca.ssvd(Y)
  twos <- fpca2s(Y)

  expect_equal(sc$Yhat, unname(face$Yhat), tolerance=.01)
  expect_equal(sc$Yhat, ssvd$Yhat, tolerance=.01)
  expect_equal(sc$Yhat, twos$Yhat, tolerance=.01)

  if(FALSE){
    ##TODO: - fix quadrature weights first
    ##      - flip sign of efunctions if necessary
    expect_equal(sc$efunctions, face$efunctions, tolerance=.01)
    expect_equal(sc$efunctions, ssvd$efunctions, tolerance=.01)
    expect_equal(sc$efunctions, twos$efunctions, tolerance=.01)

    expect_equal(sc$evalues, face$evalues, tolerance=.01)
    expect_equal(sc$evalues, ssvd$evalues, tolerance=.01)
    expect_equal(sc$evalues, twos$evalues, tolerance=.01)
  }
})

test_that("fpca.sc options  wrk", {
  sc <- fpca.sc(Y)
  sc_cov1 <- fpca.sc(Y, cov.est.method = 1)
  sc_sym <- fpca.sc(Y, useSymm = TRUE)
  sc_int <- fpca.sc(Y, random.int = TRUE)

  expect_equal(sc$Yhat, sc_cov1$Yhat, tolerance=.02)
  expect_equal(sc$Yhat, sc_sym$Yhat, tolerance=.01)
  expect_equal(sc$Yhat, sc_int$Yhat, tolerance=.01)
})
