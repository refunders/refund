context("Testing pffrSim")

###############################################################################
# Structure and Grid Tests
###############################################################################

test_that("pffrSim handles sparse data with propmissing", {
  skip_on_cran()
  set.seed(53)

  n <- 40
  nygrid <- 40
  propmissing <- 0.3

  dat <- pffr_simulate(
    Y ~ s(xsmoo),
    n = n,
    nygrid = nygrid,
    effects = list(xsmoo = "sine"),
    SNR = 50,
    propmissing = propmissing
  )

  expect_true(is.list(dat))
  expect_true(all(c("data", "ydata") %in% names(dat)))
  expect_equal(length(attr(dat, "yindex")), nygrid)

  n_total <- n * nygrid
  n_missing <- round(propmissing * n_total)
  expect_equal(nrow(dat$ydata), n_total - n_missing)
  expect_true(is.matrix(dat$data$Y))
  expect_equal(dim(dat$data$Y), c(n, nygrid))
})

test_that("pffrSim respects xind/yind grid sizes", {
  skip_on_cran()
  set.seed(54)

  n <- 30

  nxgrid1 <- 20
  nygrid1 <- 25
  s1 <- seq(0, 1, length.out = nxgrid1)

  dat1 <- pffr_simulate(
    Y ~ ff(X1, xind = s1),
    n = n,
    nxgrid = nxgrid1,
    nygrid = nygrid1,
    effects = list(X1 = "cosine"),
    SNR = 50
  )

  expect_equal(length(attr(dat1, "xindex")), nxgrid1)
  expect_equal(length(attr(dat1, "yindex")), nygrid1)
  expect_equal(dim(dat1$Y), c(n, nygrid1))

  nxgrid2 <- 50
  nygrid2 <- 60
  s2 <- seq(0, 1, length.out = nxgrid2)

  dat2 <- pffr_simulate(
    Y ~ ff(X1, xind = s2),
    n = n,
    nxgrid = nxgrid2,
    nygrid = nygrid2,
    effects = list(X1 = "cosine"),
    SNR = 50
  )

  expect_equal(length(attr(dat2, "xindex")), nxgrid2)
  expect_equal(length(attr(dat2, "yindex")), nygrid2)
  expect_equal(dim(dat2$Y), c(n, nygrid2))
})

test_that("pffrSim formula interface returns valid truth structure", {
  set.seed(100)

  n <- 30
  nygrid <- 40
  nxgrid <- 25

  dat <- pffr_simulate(
    Y ~ ff(X1, xind = s) + xlin,
    n = n,
    nxgrid = nxgrid,
    nygrid = nygrid,
    effects = list(X1 = "cosine", xlin = "dnorm"),
    SNR = 50
  )

  expect_valid_truth(dat, n, nygrid)

  truth <- attr(dat, "truth")
  expect_true(is.list(truth$etaTerms))
  expect_true(length(truth$etaTerms) >= 2)
  expect_lt(abs(mean(truth$epsilon)), 0.5)

  term_sum <- Reduce(`+`, truth$etaTerms)
  if ("intercept" %in% names(truth$etaTerms)) {
    expect_equal(truth$eta, term_sum, tolerance = 1e-10)
  }
})


###############################################################################
# Non-Gaussian Simulation Tests
###############################################################################

test_that("pffrSim with family=binomial() produces binary responses", {
  skip_on_cran()
  set.seed(60)

  n <- 50
  nygrid <- 30

  dat <- pffr_simulate(
    Y ~ xlin,
    n = n,
    nygrid = nygrid,
    effects = list(xlin = "dnorm"),
    SNR = 5,
    family = binomial()
  )

  Y <- dat$Y
  expect_true(all(Y %in% c(0, 1)))
  expect_equal(dim(Y), c(n, nygrid))

  truth <- attr(dat, "truth")
  expect_valid_truth(dat, n, nygrid)
  expect_gt(sd(truth$eta), 0.1)
})

test_that("pffrSim with family=poisson() produces count responses", {
  skip_on_cran()
  set.seed(61)

  n <- 50
  nygrid <- 30

  dat <- pffr_simulate(
    Y ~ xlin,
    n = n,
    nygrid = nygrid,
    effects = list(xlin = "dnorm"),
    SNR = 5,
    family = poisson()
  )

  Y <- dat$Y
  expect_true(all(Y >= 0))
  expect_true(all(Y == floor(Y)))
  expect_equal(dim(Y), c(n, nygrid))
  expect_gt(sd(Y), 0)
})

test_that("pffrSim rejects unsupported families with clear error", {
  skip_on_cran()
  set.seed(62)

  expect_error(
    pffr_simulate(Y ~ xlin, n = 10, nygrid = 10, SNR = 5, family = mgcv::nb()),
    "does not support family"
  )
})


###############################################################################
# Effect Selection and Presets
###############################################################################

test_that("pffrSim effect selection via %||% works", {
  skip_on_cran()
  set.seed(70)

  n <- 30
  nygrid <- 30

  custom_effect <- function(t) sin(2 * pi * t)

  dat <- pffr_simulate(
    Y ~ xlin,
    n = n,
    nygrid = nygrid,
    effects = list(xlin = custom_effect),
    SNR = 50
  )

  expect_valid_truth(dat, n, nygrid)
  truth <- attr(dat, "truth")
  expect_true("xlin" %in% names(truth$etaTerms))
})

test_that("intercept presets work correctly", {
  skip_on_cran()
  set.seed(73)

  n <- 40
  nygrid <- 30

  for (intercept_type in c("constant", "beta", "sine", "zero")) {
    dat <- pffr_simulate(
      Y ~ 1,
      n = n,
      nygrid = nygrid,
      intercept = intercept_type,
      SNR = 100
    )

    expect_valid_truth(dat, n, nygrid)

    truth <- attr(dat, "truth")

    if (intercept_type == "zero") {
      expect_equal(max(abs(truth$eta)), 0)
    } else if (intercept_type == "constant") {
      expect_lt(sd(as.vector(truth$eta)), 0.01)
    } else {
      eta_without_noise <- truth$eta - truth$epsilon
      expect_gt(sd(eta_without_noise[1, ]), 0.1)
    }
  }
})

test_that("factor terms with custom 2-arg function work", {
  skip_on_cran()
  set.seed(82)

  n <- 40
  nygrid <- 30

  custom_factor_effect <- function(fac_numeric, t) {
    outer(fac_numeric, sin(2 * pi * t))
  }

  dat <- pffr_simulate(
    Y ~ xfactor,
    n = n,
    nygrid = nygrid,
    effects = list(xfactor = custom_factor_effect),
    SNR = 50
  )

  expect_true(is.factor(dat$xfactor))
  expect_valid_truth(dat, n, nygrid)
})

test_that("pffrSim respects explicit intercept removal (Y ~ 0 + ...)", {
  skip_on_cran()
  set.seed(84)

  n <- 30
  nygrid <- 25

  dat_with <- pffr_simulate(
    Y ~ xlin,
    n = n,
    nygrid = nygrid,
    effects = list(xlin = "dnorm"),
    SNR = 50
  )
  truth_with <- attr(dat_with, "truth")
  expect_true("intercept" %in% names(truth_with$etaTerms))

  dat_without <- pffr_simulate(
    Y ~ 0 + xlin,
    n = n,
    nygrid = nygrid,
    effects = list(xlin = "dnorm"),
    SNR = 50
  )
  truth_without <- attr(dat_without, "truth")
  expect_false("intercept" %in% names(truth_without$etaTerms))

  dat_minus1 <- pffr_simulate(
    Y ~ -1 + xlin,
    n = n,
    nygrid = nygrid,
    effects = list(xlin = "dnorm"),
    SNR = 50
  )
  truth_minus1 <- attr(dat_minus1, "truth")
  expect_false("intercept" %in% names(truth_minus1$etaTerms))

  expect_false(all(truth_with$eta == truth_without$eta))
})

test_that("pffrSim handles empty term sets (Y ~ 0) without error", {
  skip_on_cran()
  set.seed(87)

  n <- 20
  nygrid <- 25

  dat <- pffr_simulate(Y ~ 0, n = n, nygrid = nygrid, SNR = 10)

  Y <- dat$Y
  truth <- attr(dat, "truth")

  expect_equal(dim(Y), c(n, nygrid))
  expect_equal(dim(truth$eta), c(n, nygrid))
  expect_equal(max(abs(truth$eta)), 0)
  expect_gt(sd(Y), 0)

  dat2 <- pffr_simulate(Y ~ -1, n = n, nygrid = nygrid, SNR = 10)
  expect_equal(dim(dat2$Y), c(n, nygrid))
})

test_that("numeric smooth effects produce correct dimensions", {
  skip_on_cran()
  set.seed(85)

  n <- 25
  nygrid <- 30

  dat <- pffr_simulate(
    Y ~ s(xsmoo),
    n = n,
    nygrid = nygrid,
    effects = list(xsmoo = 2),
    SNR = 50
  )

  truth <- attr(dat, "truth")
  expect_equal(dim(truth$eta), c(n, nygrid))
  for (term_name in names(truth$etaTerms)) {
    expect_equal(
      dim(truth$etaTerms[[term_name]]),
      c(n, nygrid),
      info = paste("Term", term_name, "should have n x nygrid dims")
    )
  }
})

test_that("unwrapped te/ti/t2 uses all covariates", {
  skip_on_cran()
  set.seed(86)

  n <- 30
  nygrid <- 25

  dat <- pffr_simulate(
    Y ~ te(xte1, xte2),
    n = n,
    nygrid = nygrid,
    SNR = 50
  )

  expect_true("xte1" %in% names(dat))
  expect_true("xte2" %in% names(dat))

  truth <- attr(dat, "truth")
  expect_equal(dim(truth$eta), c(n, nygrid))

  custom_te_effect <- function(x1, x2, t) {
    outer(x1 * x2, sin(2 * pi * t))
  }

  dat2 <- pffr_simulate(
    Y ~ te(xte1, xte2),
    n = n,
    nygrid = nygrid,
    effects = list("te(xte1,xte2)" = custom_te_effect),
    SNR = 50
  )

  truth2 <- attr(dat2, "truth")
  expect_equal(dim(truth2$eta), c(n, nygrid))
})
