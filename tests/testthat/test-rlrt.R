context("Testing rlrt.pfr")
library(refundDevel)

test_that("basis rlrt tests are working", {
  skip_on_cran()

  set.seed(97867)
  data(DTI2)

  ## fit two univariate models, then one model with both functional predictors
  suppressWarnings({
    pfr.obj.t1 <- pfr(Y = DTI2$pasat,
      covariates= matrix(as.numeric(DTI2$visit > 1), ncol=1),
      funcs = list(DTI2$cca),
      subj = DTI2$id, kz = 10, kb = 50)
    pfr.obj.t2 <- pfr(Y = DTI2$pasat,
      covariates=matrix(as.numeric(DTI2$visit > 1), ncol=1),
      funcs = list(DTI2$rcst),
      subj = DTI2$id, kz = 10, kb = 50)
    pfr.obj.t3 <- pfr(Y = DTI2$pasat,
      covariates=matrix(as.numeric(DTI2$visit > 1), ncol=1),
      funcs = list(DTI2$cca, DTI2$rcst),
      subj = DTI2$id, kz = 10, kb = 50)
  })

 ## do some testing
 t1 <- rlrt.pfr(pfr.obj.t1, "constancy")
 t2 <- rlrt.pfr(pfr.obj.t2, "constancy")
 suppressWarnings({
   t3 <- rlrt.pfr(pfr.obj.t3, "inclusion")
 })

  expect_is(t1, "list")
  expect_is(t2, "list")
  expect_is(t3, "list")
})

test_that("tests work with subj = NULL", {
  skip_on_cran()

  set.seed(67887)
  data(DTI2)
  O  <- DTI2$pasat ## PASAT outcome
  id <- DTI2$id    ## subject id
  W1 <- DTI2$cca   ## Corpus Callosum
  W2 <- DTI2$rcst  ## Right corticospinal
  V  <- DTI2$visit ## visit

  ## prep scalar covariate
  covar.in <- matrix(as.numeric(V > 1), ncol=1)

  suppressWarnings({
    pfr.obj.t1 <- pfr(Y = DTI2$pasat,
      covariates=matrix(as.numeric(DTI2$visit > 1), ncol=1),
      funcs = list(DTI2$cca),
      subj = NULL, kz = 10, kb = 50)
    pfr.obj.t2 <- pfr(Y = DTI2$pasat,
      covariates=matrix(as.numeric(DTI2$visit > 1), ncol=1),
      funcs = list(DTI2$rcst),
      subj = NULL, kz = 10, kb = 50)
    pfr.obj.t3 <- pfr(Y = DTI2$pasat,
      covariates=matrix(as.numeric(DTI2$visit > 1), ncol=1),
      funcs = list(DTI2$cca, DTI2$rcst),
      subj = NULL, kz = 10, kb = 50)
  })
  t1 <- rlrt.pfr(pfr.obj.t1, "constancy")
  t2 <- rlrt.pfr(pfr.obj.t2, "constancy")
  suppressWarnings(t3 <- rlrt.pfr(pfr.obj.t3, "inclusion"))

  expect_is(t1, "list")
  expect_is(t2, "list")
  expect_is(t3, "list")
})
