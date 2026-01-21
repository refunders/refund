###############################################################################
# Tests for pffr with AR(1) correlated errors
###############################################################################

test_that("pffr builds AR.start for dense data when rho is supplied", {
  skip_on_cran()

  sim <- get_ar_data()
  tgrid <- attr(sim, "yindex")
  fit <- pffr(
    Y ~ c(1),
    data = sim,
    yind = tgrid,
    algorithm = "bam",
    rho = 0.35,
    bs.int = list(bs = "ps", k = length(tgrid), m = c(2, 1))
  )

  expect_s3_class(fit, "pffr")
  expect_equal(fit$AR1.rho, 0.35, tolerance = 1e-8)

  smry_lines <- capture.output(summary(fit))
  expect_true(any(grepl("AR\\(1\\) residual correlation", smry_lines)))

  preds <- predict(fit)
  expect_true(is.matrix(preds))
  expect_equal(dim(preds), dim(sim$Y))

  ar_col <- fit$model[["(AR.start)"]]
  expect_true(is.logical(ar_col))
  expect_equal(sum(ar_col), nrow(sim$Y))
})

test_that("unsupported AR settings throw informative errors", {
  skip_on_cran()

  sim <- get_ar_data()
  tgrid <- attr(sim, "yindex")

  expect_error(
    pffr(Y ~ c(1), data = sim, yind = tgrid, rho = 1.2),
    "`rho` must have absolute value strictly less than 1"
  )
  expect_error(
    pffr(Y ~ c(1), data = sim, yind = tgrid, rho = "bad"),
    "`rho` must be a single numeric value"
  )
  expect_error(
    pffr(Y ~ c(1), data = sim, yind = tgrid, algorithm = "gam", rho = 0.2),
    "Autocorrelated errors via `rho` are currently supported only when `algorithm = \"bam\"`\\."
  )
  expect_error(
    pffr(Y ~ c(1), data = sim, yind = tgrid, rho = 0.2, family = binomial()),
    "Autocorrelated errors \\(via `rho`\\) require either a Gaussian identity model or setting `discrete = TRUE`"
  )
  expect_error(
    pffr(
      Y ~ c(1),
      data = sim,
      yind = tgrid,
      rho = 0.2,
      AR.start = rep(TRUE, length = nrow(sim$Y) * ncol(sim$Y))
    ),
    "Please do not supply `AR.start` directly"
  )
})

test_that("pffr AR fits match mgcv::bam on stacked data", {
  skip_on_cran()

  sim <- get_ar_data()
  tgrid <- attr(sim, "yindex")
  fit <- pffr(
    Y ~ c(1),
    data = sim,
    yind = tgrid,
    algorithm = "bam",
    rho = 0.25,
    bs.int = list(bs = "ps", k = length(tgrid), m = c(2, 1))
  )

  mgcv_data <- fit$model
  mgcv_ar <- mgcv_data[["(AR.start)"]]
  mgcv_formula <- fit$formula
  mgcv_fit <- do.call(
    mgcv::bam,
    list(
      formula = mgcv_formula,
      data = mgcv_data,
      rho = 0.25,
      AR.start = mgcv_ar,
      method = "fREML"
    )
  )

  fit_stripped <- fit
  class(fit_stripped) <- setdiff(class(fit_stripped), "pffr")
  expect_equal(fitted(fit_stripped), fitted(mgcv_fit))
  expect_equal(fit_stripped$AR1.rho, mgcv_fit$AR1.rho)
})

test_that("pffr builds AR.start for sparse responses", {
  skip_on_cran()
  set.seed(2504)
  sim_sparse <- pffr_simulate(
    scenario = "smoo",
    n = 12,
    nygrid = 18,
    propmissing = 0.35
  )
  tgrid <- attr(sim_sparse, "yindex")
  ydata <- sim_sparse$ydata[
    order(sim_sparse$ydata$.obs, sim_sparse$ydata$.index),
  ]

  fit_sparse <- pffr(
    Y ~ s(xsmoo),
    data = sim_sparse$data,
    ydata = ydata,
    yind = tgrid,
    algorithm = "bam",
    rho = 0.4,
    bs.int = list(bs = "ps", k = length(tgrid), m = c(2, 1))
  )

  smry_lines <- capture.output(summary(fit_sparse))
  expect_true(any(grepl("AR\\(1\\) residual correlation", smry_lines)))

  ar_sparse <- fit_sparse$model[["(AR.start)"]]
  expect_true(is.logical(ar_sparse))
  expect_equal(sum(ar_sparse), nrow(sim_sparse$data))
})

test_that("binomial models can use rho when discrete sampling is enabled", {
  skip_on_cran()
  set.seed(2505)
  n <- 10
  ny <- 8
  binary_Y <- matrix(rbinom(n * ny, size = 1, prob = 0.4), nrow = n, ncol = ny)
  df <- data.frame(Y = I(binary_Y))
  tgrid <- seq(0, 1, length.out = ny)

  fit_binom <- pffr(
    Y ~ c(1),
    data = df,
    yind = tgrid,
    algorithm = "bam",
    rho = 0.3,
    family = binomial(),
    discrete = TRUE,
    bs.int = list(bs = "ps", k = length(tgrid), m = c(2, 1))
  )
  expect_equal(fit_binom$family$family, "binomial")
  expect_equal(fit_binom$AR1.rho, 0.3, tolerance = 1e-8)
})
