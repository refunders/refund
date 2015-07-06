context("Testing rlrt.pfr")
library(refundDevel)

test_that("basis rlrt tests are working", {
  skip_on_cran()

  set.seed(97867)
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
 suppressWarnings({
   t3 <- rlrt.pfr(pfr.obj.t3, "inclusion")
 })

  expect_is(t1, "list")
  expect_is(t2, "list")
  expect_is(t3, "list")
})

test_that("tests work with subj = NULL", {
  skip_on_cran()
  rm(list=ls())
  set.seed(67887)
  data(DTI2)
  O  <- DTI2$pasat ## PASAT outcome
  id <- DTI2$id    ## subject id
  W1 <- DTI2$cca   ## Corpus Callosum
  W2 <- DTI2$rcst  ## Right corticospinal
  V  <- DTI2$visit ## visit

  ## prep scalar covariate
  visit.1.rest <- matrix(as.numeric(V > 1), ncol=1)
  covar.in <- visit.1.rest

  suppressWarnings({
    pfr.obj.t1 <- pfr(Y = O, covariates=covar.in, funcs = list(W1),     subj = NULL, kz = 10, kb = 50)
    pfr.obj.t2 <- pfr(Y = O, covariates=covar.in, funcs = list(W2),     subj = NULL, kz = 10, kb = 50)
    ## pfr.obj.t3 <- pfr(Y = O, covariates=covar.in, funcs = list(W1, W2), subj = NULL, kz = 10, kb = 50)
  })
  t1 <- rlrt.pfr(pfr.obj.t1, "constancy")
  t2 <- rlrt.pfr(pfr.obj.t2, "constancy")
  # t3 <- rlrt.pfr(pfr.obj.t3, "inclusion")

  expect_is(t1, "list")
  expect_is(t2, "list")
  # expect_is(t3, "list")
})
