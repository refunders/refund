context("Testing ccb.fpc")
library(refundDevel)

test_that("ccb.fpc works as expected", {
    skip_on_cran()

  data(cd4)

  ## obtain a subsample of the data with 25 subjects
  set.seed(1236)
  sample = sample(1:dim(cd4)[1], 25)
  Y.sub = cd4[sample,]

  # obtain a mixed-model based FPCA decomposition
  Fit.MM = try(fpca.sc(Y.sub, var = TRUE, simul = TRUE), TRUE)
  if (inherits(Fit.MM, "try-error"))
    skip("fpca.sc seems to not be working")
  expect_equal_to_reference(Fit.MM, file = "ccb.fpca.obj.rds")

  # use iterated variance to obtain curve estimates and variances
  Fit.IV = suppressMessages(ccb.fpc(Y.sub, n.boot = 25, simul = TRUE))
  expect_equal_to_reference(Fit.IV, "ccb.obj.rds")
})
