context("Testing fpca.xxx")

set.seed(12212)
n <- 100
ngrid <- 40
t <- seq(0, 1, l=ngrid)
efcts <- poly(t, 2)
Y <- outer(2 * drop(scale(rnorm(n))), efcts[, 1]) +
  outer(drop(scale(rnorm(n))), efcts[, 2])

flip_efunctions <- function(ef1, ef2) {
  squared_diff <- function(x,y) crossprod(x - y)
  for(i in 1:ncol(ef1)){
    if(squared_diff(ef1[,i], ef2[,i]) > squared_diff(ef1[,i], - ef2[,i])){
      ef2[,i] <- -ef2[,i]
    }
  }
  ef2
}

test_that("all fpca functions agree on toy example", {
  skip_on_cran()

  sc <- fpca.sc(Y)
  face <- fpca.face(Y)
  #ssvd <- fpca.ssvd(Y)
  #twos <- fpca2s(Y)

  expect_equal(sc$Yhat, unname(face$Yhat), tolerance=.01)
  #expect_equal(sc$Yhat, ssvd$Yhat, tolerance=.01)
  #expect_equal(sc$Yhat, twos$Yhat, tolerance=.01)

  #ssvd$efunctions <- flip_efunctions(sc$efunctions, ssvd$efunctions)
  #expect_equal(sc$efunctions, ssvd$efunctions, tolerance=.1)
  #expect_equal(sc$evalues, ssvd$evalues, tolerance=.1)

  #twos$efunctions <- flip_efunctions(sc$efunctions, twos$efunctions)
  #expect_equal(sc$efunctions, twos$efunctions, tolerance=.1)
  #expect_equal(sc$evalues, twos$evalues, tolerance=.1)

  if(FALSE){
    ##TODO: - fix quadrature weights first
    ##      - flip sign of efunctions if necessary
    expect_equal(sc$efunctions, face$efunctions, tolerance=.01)
    #expect_equal(sc$efunctions, twos$efunctions, tolerance=.01)
    expect_equal(sc$evalues, face$evalues, tolerance=.01)
    #expect_equal(sc$evalues, ssvd$evalues, tolerance=.01)
    #expect_equal(sc$evalues, twos$evalues, tolerance=.01)
  }
})

test_that("fpca.sc options work", {
  skip_on_cran()

  sc <- fpca.sc(Y)
  sc_cov1 <- fpca.sc(Y, cov.est.method = 1)
  #sc_sym <- fpca.sc(Y, useSymm = TRUE)
  sc_int <- fpca.sc(Y, random.int = TRUE)

  expect_equal(sc$Yhat, sc_cov1$Yhat, tolerance=.02)
  #expect_equal(sc$Yhat, sc_sym$Yhat, tolerance=.01)
  expect_equal(sc$Yhat, sc_int$Yhat, tolerance=.01)
})


# test_that("fpca.ssvd options work", {
#   skip_on_cran()
# 
#   expect_error(fpca.ssvd(Y = 1:10, ydata=data.frame()), "irregular data")
#   expect_warning(fpca.ssvd(Y = Y, argvals=sqrt(t)), "non-equidistant")
#   ssvd <- fpca.ssvd(Y)
#   ssvd_npc1 <- fpca.ssvd(Y, npc=1)
#   ssvd_d2 <- fpca.ssvd(Y, diffpen = 2)
#   expect_equal(ssvd_npc1$efunctions[,1], ssvd$efunctions[,1])
#   expect_true(ncol(ssvd_npc1$efunctions) == 1)
#   expect_equal(ssvd_d2$efunctions, ssvd$efunctions, tol=.01)
# })
# 
# test_that("fpca2s options work", {
#   skip_on_cran()
# 
#   expect_error(fpca2s(Y = 1:10, ydata=data.frame()), "irregular data")
#   expect_warning(fpca2s(Y = Y, argvals=sqrt(t)), "non-equidistant")
#   twos <- fpca2s(Y)
#   twos_npc1 <- fpca2s(Y, npc=1)
#   expect_equal(twos_npc1$efunctions[,1], twos$efunctions[,1])
#   expect_true(ncol(twos_npc1$efunctions) == 1)
# })
