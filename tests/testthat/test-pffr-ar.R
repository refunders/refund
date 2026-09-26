###############################################################################
# Tests for pffr with AR(1) correlated errors
###############################################################################

test_that("pffr builds AR.start for dense data when rho is supplied", {
  skip_on_cran()

  sim <- get_ar_data()
  tgrid <- attr(sim, "yindex")
  fit <- quiet_pffr(
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
  sim_bin <- sim
  sim_bin$Y <- I(1L * (sim$Y > 0))
  fit <- quiet_pffr(
    Y ~ c(1),
    data = sim_bin,
    yind = tgrid,
    rho = 0.2,
    family = binomial(),
    bs.int = list(bs = "ps", k = length(tgrid), m = c(2, 1))
  )
  expect_true(isTRUE(fit$call$discrete))
  expect_error(
    pffr(
      Y ~ c(1),
      data = sim,
      yind = tgrid,
      rho = 0.2,
      family = binomial(),
      discrete = FALSE
    ),
    "Autocorrelated errors \\(via `rho`\\) require either a Gaussian identity model or `discrete = TRUE`"
  )
  expect_error(
    pffr(Y ~ c(1), data = sim, yind = tgrid, rho = 0.2, method = "REML"),
    "Autocorrelated errors via `rho` require `method = \"fREML\"`"
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
  fit <- quiet_pffr(
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

  fit_sparse <- quiet_pffr(
    Y ~ s(xsmoo),
    data = sim_sparse$data,
    ydata = ydata,
    yind = tgrid,
    algorithm = "bam",
    rho = 0.4,
    # sandwich estimators are unavailable with an AR(1) working correlation
    sandwich = "none",
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

  fit_binom <- quiet_pffr(
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

test_that("explicit sandwich requests with rho error before fitting", {
  skip_on_cran()

  sim <- get_ar_data()
  tgrid <- attr(sim, "yindex")
  fit_ar <- function(sandwich) {
    pffr(
      Y ~ c(1),
      data = sim,
      yind = tgrid,
      rho = 0.3,
      sandwich = sandwich,
      bs.int = list(bs = "ps", k = length(tgrid), m = c(2, 1))
    )
  }
  for (type in c("cluster", "cl2", "hc")) {
    expect_error(
      fit_ar(type),
      paste0(
        "sandwich = \"",
        type,
        "\" is not supported for fits with an AR\\(1\\) working correlation"
      )
    )
  }
  expect_error(fit_ar(TRUE), "sandwich = \"cluster\" is not supported")
  expect_error(fit_ar("cluster"), "Use sandwich = \"none\"")
})

test_that("default sandwich resolves to none for AR(1) fits", {
  skip_on_cran()

  sim <- get_ar_data()
  tgrid <- attr(sim, "yindex")
  fit_default <- function(...) {
    pffr(
      Y ~ c(1),
      data = sim,
      yind = tgrid,
      rho = 0.3,
      bs.int = list(bs = "ps", k = length(tgrid), m = c(2, 1)),
      ...
    )
  }

  msgs <- character()
  fit <- withCallingHandlers(
    expect_no_warning(fit_default()),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  expect_length(grep("resolved to \"none\"", msgs), 1)
  expect_false(any(grepl("now defaults to sandwich", msgs)))
  expect_identical(fit$pffr$sandwich, "none")
  expect_null(fit$pffr[["Vsandwich"]])
  expect_equal(fit$AR1.rho, 0.3)

  # explicit "auto" is a request for the policy, not a specific sandwich
  fit_auto <- suppressMessages(fit_default(sandwich = "auto"))
  expect_identical(fit_auto$pffr$sandwich, "none")

  # model-based inference works; post-hoc sandwich requests error
  cf <- coef(fit)
  expect_true(all(is.finite(cf$smterms[[1]]$coef$se)))
  expect_error(coef(fit, sandwich = "cluster"), "AR\\(1\\) working correlation")
  expect_error(coef(fit, sandwich = "cl2"), "AR\\(1\\) working correlation")
  expect_error(vcov(fit, sandwich = TRUE), "AR\\(1\\) working correlation")
  expect_equal(unname(vcov(fit)), unname(fit$Vp))
})

test_that("the shared sandwich guard is a no-op without rho", {
  expect_true(pffr_check_sandwich_ar1(list(AR1.rho = NULL), "cluster"))
  expect_true(pffr_check_sandwich_ar1(list(AR1.rho = 0), "cl2"))
  expect_true(pffr_check_sandwich_ar1(list(AR1.rho = 0.5), "none"))
  expect_error(
    pffr_check_sandwich_ar1(list(AR1.rho = 0.5), "hc"),
    "sandwich = \"hc\""
  )
})
