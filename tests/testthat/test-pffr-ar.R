context("pffr autocorrelated errors")

test_that("pffr builds AR.start for dense data when rho is supplied", {
  set.seed(2501)
  sim <- pffrSim(scenario = "int", n = 8, nygrid = 12)
  tgrid <- attr(sim, "yindex")
  fit <- pffr(
    Y ~ c(1),
    data = sim,
    yind = tgrid,
    rho = 0.35,
    bs.int = list(bs = "ps", k = length(tgrid), m = c(2, 1))
  )

  expect_is(fit, "pffr")
  expect_equal(fit$AR1.rho, 0.35, tolerance = 1e-8)

  smry_lines <- capture.output(summary(fit))
  expect_true(any(grepl("AR\\(1\\) residual correlation", smry_lines)))

  preds <- predict(fit)
  expect_is(preds, "matrix")
  expect_equal(dim(preds), dim(sim$Y))

  ar_col <- fit$model[["(AR.start)"]]
  expect_true(is.logical(ar_col))
  expect_equal(sum(ar_col), nrow(sim$Y))
})

test_that("unsupported AR settings throw informative errors", {
  set.seed(2502)
  sim <- pffrSim(scenario = "int", n = 4, nygrid = 6)
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
  set.seed(2503)
  sim <- pffrSim(scenario = "int", n = 9, nygrid = 10)
  tgrid <- attr(sim, "yindex")
  fit <- pffr(
    Y ~ c(1),
    data = sim,
    yind = tgrid,
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
  set.seed(2504)
  sim_sparse <- pffrSim(
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
    rho = 0.3,
    family = binomial(),
    discrete = TRUE,
    bs.int = list(bs = "ps", k = length(tgrid), m = c(2, 1))
  )
  expect_equal(fit_binom$family$family, "binomial")
  expect_equal(fit_binom$AR1.rho, 0.3, tolerance = 1e-8)
})

test_that("AR modelling maintains nominal size under AR(1) noise", {
  skip_on_cran()
  set.seed(2506)
  n <- 25
  ny <- 15
  nrep <- 500
  rho_true <- 0.6
  alpha <- 0.05
  tgrid <- seq(0, 1, length.out = ny)
  mu <- sin(2 * pi * tgrid)
  rej_ar <- logical(nrep)
  rej_iid <- logical(nrep)

  for (rep in seq_len(nrep)) {
    z <- rnorm(n)
    Y <- matrix(0, nrow = n, ncol = ny)
    for (i in seq_len(n)) {
      eps <- as.numeric(arima.sim(list(ar = rho_true), ny))
      Y[i, ] <- mu + eps
    }
    df <- data.frame(z = z, Y = I(Y))
    fit_ar <- pffr(
      Y ~ c(z),
      data = df,
      yind = tgrid,
      rho = rho_true,
      algorithm = "bam",
      bs.int = list(bs = "ps", k = ny, m = c(2, 1))
    )
    fit_iid <- pffr(
      Y ~ c(z),
      data = df,
      yind = tgrid,
      algorithm = "bam",
      bs.int = list(bs = "ps", k = ny, m = c(2, 1))
    )
    pv_ar <- summary(fit_ar)$p.table["z", "Pr(>|t|)"]
    pv_iid <- summary(fit_iid)$p.table["z", "Pr(>|t|)"]
    rej_ar[rep] <- pv_ar < alpha
    rej_iid[rep] <- pv_iid < alpha
  }

  rate_ar <- mean(rej_ar)
  rate_iid <- mean(rej_iid)
  expect_lt(rate_ar, 0.1)
  expect_gt(rate_iid, rate_ar + 0.1)
})
