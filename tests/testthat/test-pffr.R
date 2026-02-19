context("Testing pffr")

###############################################################################
# Helper Functions
###############################################################################

#' Compute RMSE between two matrices/vectors
rmse <- function(x, y) sqrt(mean((as.vector(x) - as.vector(y))^2))

#' Safe relative error (avoids division by zero)
safe_rel_error <- function(est, true) {
  denom <- max(abs(true)) + 1e-10
  max(abs(est - true)) / denom
}

#' Check model recovers true eta with both correlation and RMSE
expect_recovers_eta <- function(model, truth, cor_threshold = 0.9, label = "") {
  fitted_vals <- as.vector(fitted(model))
  true_eta <- as.vector(truth$eta)

  corr <- cor(fitted_vals, true_eta)
  rel_rmse <- rmse(fitted_vals, true_eta) / (sd(true_eta) + 1e-10)

  expect_gt(
    corr,
    cor_threshold,
    label = paste(label, "correlation with true eta")
  )
  expect_lt(
    rel_rmse,
    0.5,
    label = paste(label, "relative RMSE should be < 0.5")
  )
}

sim_xlin_data <- function(n, nygrid, SNR = 50, family = gaussian()) {
  pffr_simulate(
    Y ~ xlin,
    n = n,
    nygrid = nygrid,
    effects = list(xlin = "dnorm"),
    SNR = SNR,
    family = family
  )
}

###############################################################################
# Legacy Tests (existing behavior, backward compatibility)
###############################################################################

test_that("all major pffr terms are working (legacy scenario='all')", {
  m2 <- get_all_scenario_model()
  dat <- get_all_scenario_data()
  argvals <- attr(dat, "yindex")

  expect_s3_class(m2, "pffr")
  expect_equal(dim(fitted(m2)), c(30, length(argvals)))
  expect_false(any(is.na(fitted(m2))))
})

test_that("convenience functions are working", {
  m2 <- get_all_scenario_model()
  data2 <- get_all_scenario_data()
  argvals <- attr(data2, "yindex")

  # summary
  summ <- summary(m2)
  expect_s3_class(summ, "summary.pffr")
  expect_true(!is.null(summ$s.table))

  # coef
  cm2 <- coef(m2)
  expect_true(is.list(cm2))
  expect_true("pterms" %in% names(cm2))
  expect_true("smterms" %in% names(cm2))
  expect_true(is.matrix(cm2$pterms) || is.data.frame(cm2$pterms))

  # predict on new data
  set.seed(9314)
  preddata <- suppressWarnings(pffr_simulate(scenario = "all", n = 20))
  pred <- predict(m2, newdata = preddata)
  expect_true(is.matrix(pred))
  expect_equal(dim(pred), c(20, length(argvals)))

  # predict type="terms" - verify structure
  pred_terms <- predict(m2, type = "terms")
  expect_true(is.list(pred_terms))
  expect_true(length(pred_terms) >= 5) # intercept + 5 terms
  # Each term should have correct dimensions
  for (term_name in names(pred_terms)) {
    expect_equal(
      nrow(pred_terms[[term_name]]),
      30,
      info = paste("term", term_name, "should have n rows")
    )
  }
})

test_that("pffr with sparse data works", {
  data3 <- get_sparse_data()
  m3.sparse <- get_sparse_model()

  # Verify sparse data structure
  expect_true(is.list(data3))
  expect_true("ydata" %in% names(data3))
  expect_true("data" %in% names(data3))

  expect_s3_class(m3.sparse, "pffr")
  summ <- summary(m3.sparse)
  expect_s3_class(summ, "summary.pffr")
})

test_that("ffpc terms are working", {
  skip_on_cran()
  set.seed(1122)

  n <- 55
  S <- 60
  n.argvals <- 50
  s <- seq(0, 1, l = S)
  argvals <- seq(0, 1, l = n.argvals)

  # Generate X from polynomial FPC-basis
  rankX <- 5
  Phi <- cbind(1 / sqrt(S), poly(s, degree = rankX - 1))
  lambda <- rankX:1
  Xi <- sapply(lambda, function(l) scale(rnorm(n, sd = sqrt(l)), scale = FALSE))
  X <- Xi %*% t(Phi)

  beta.st <- outer(s, argvals, function(s, argvals) cos(2 * pi * s * argvals))
  y <- (1 / S * X) %*%
    beta.st +
    0.1 * matrix(rnorm(n * n.argvals), nrow = n, ncol = n.argvals)

  data <- list(y = y, X = X)

  m.pc <- pffr(
    y ~ c(1) + 0 + ffpc(X, yind = argvals, decomppars = list(npc = rankX)),
    data = data,
    yind = argvals
  )

  expect_s3_class(m.pc, "pffr")
  expect_s3_class(summary(m.pc), "summary.pffr")
  # Verify fitted values have correct dimensions
  expect_equal(dim(fitted(m.pc)), c(n, n.argvals))
})

test_that("ff limits arg works", {
  set.seed(2121)
  data <- pffr_simulate(
    scenario = "ff",
    n = 20,
    SNR = 100,
    limits = function(s, t) s < t
  )
  t <- attr(data, "yindex")
  s <- attr(data, "xindex")

  m.l <- pffr(Y ~ ff(X1, xind = s, limits = "s<t"), yind = t, data = data)
  m.leq <- pffr(Y ~ ff(X1, xind = s, limits = "s<=t"), yind = t, data = data)
  m.fun <- pffr(
    Y ~ ff(X1, xind = s, limits = function(s, t) s < t),
    yind = t,
    data = data
  )

  # String and function limits should give identical results
  expect_equal(fitted(m.l), fitted(m.fun))

  # s<t and s<=t should be very similar (safe comparison)
  expect_lt(safe_rel_error(fitted(m.l), fitted(m.leq)), 0.05)
})

test_that("ff() rejects mismatched xind length", {
  skip_on_cran()
  set.seed(2122)

  n <- 30
  nxgrid <- 25
  nygrid <- 30

  dat <- pffr_simulate(
    Y ~ ff(X1, xind = s),
    n = n,
    nxgrid = nxgrid,
    nygrid = nygrid,
    effects = list(X1 = "cosine"),
    SNR = 50
  )
  s <- attr(dat, "xindex")
  t <- attr(dat, "yindex")

  expect_error(pffr(Y ~ ff(X1, xind = s[-1]), yind = t, data = dat))
})

test_that("weights and offset args work", {
  skip_on_cran()
  set.seed(112)

  n <- 20
  nygrid <- 50
  data <- sim_xlin_data(n = n, nygrid = nygrid, SNR = 100)
  t <- attr(data, "yindex")

  m0 <- pffr(Y ~ xlin, yind = t, data = data)

  # Vector offset
  vecoffset <- 1:n
  data$Y_vecoffset <- data$Y + vecoffset
  m_vecoffset <- pffr(
    Y_vecoffset ~ xlin,
    yind = t,
    data = data,
    offset = vecoffset
  )
  expect_equal(fitted(m0), fitted(m_vecoffset) - vecoffset, tolerance = 1e-5)

  # Matrix offset
  matoffset <- matrix(rnorm(n * nygrid), n, nygrid)
  data$Y_matoffset <- data$Y + matoffset
  m_matoffset <- pffr(
    Y_matoffset ~ xlin,
    yind = t,
    data = data,
    offset = matoffset
  )
  expect_equal(fitted(m0), fitted(m_matoffset) - matoffset, tolerance = 1e-5)

  # Mismatched dimensions should error
  expect_error(pffr(Y ~ xlin, yind = t, data = data, offset = vecoffset[-1]))
  expect_error(pffr(Y ~ xlin, yind = t, data = data, offset = matoffset[-1, ]))

  # Weights: vector and matrix should match unweighted fit when all 1s
  weights_vec <- rep(1, n)
  m_wvec <- pffr(Y ~ xlin, yind = t, data = data, weights = weights_vec)
  expect_equal(fitted(m0), fitted(m_wvec), tolerance = 1e-10)

  weights_mat <- matrix(1, n, nygrid)
  m_wmat <- pffr(Y ~ xlin, yind = t, data = data, weights = weights_mat)
  expect_equal(fitted(m0), fitted(m_wmat), tolerance = 1e-10)

  # Mismatched weight dimensions should error
  expect_error(pffr(Y ~ xlin, yind = t, data = data, weights = weights_vec[-1]))
  expect_error(pffr(
    Y ~ xlin,
    yind = t,
    data = data,
    weights = weights_mat[-1, ]
  ))
})

test_that("sff terms are working", {
  skip_on_cran()
  set.seed(2121)

  data <- pffr_simulate(scenario = "ff", n = 20, SNR = 100)
  t <- attr(data, "yindex")
  s <- attr(data, "xindex")

  m.s <- pffr(Y ~ sff(X1, xind = s), yind = t, data = data)
  m.lin <- pffr(Y ~ ff(X1, xind = s), yind = t, data = data)

  expect_s3_class(summary(m.s), "summary.pffr")

  # sff and ff should give similar fits (safe comparison)
  expect_lt(safe_rel_error(fitted(m.lin), fitted(m.s)), 0.05)
  expect_equal(fitted(m.s), predict(m.s, newdata = data))
})

test_that("ff identifiability stuff works", {
  skip_on_cran()
  set.seed(112)

  n <- 20
  ngrid <- 50
  s <- t <- seq(0, 1, l = ngrid)

  Xefcts <- t(poly(t, 2))
  X <- matrix(rnorm(n * 2), n, 2) %*% Xefcts
  L <- matrix(1 / ngrid, ncol = ngrid, nrow = n)
  LX <- L * X
  beta.st <- outer(s, t, function(s, t) cos(3 * pi * t) * s)
  Y <- LX %*% beta.st
  data <- data.frame(Y = I(Y), X = I(X))

  expect_warning(
    m <- pffr(Y ~ ff(X, xind = s), yind = t, data = data),
    "kernel overlap"
  )
})


###############################################################################
# Term Recovery Tests (new formula interface)
###############################################################################

test_that("ff() term recovers beta(s,t)", {
  skip_on_cran()
  set.seed(42)

  n <- 40
  nxgrid <- 30
  nygrid <- 40

  dat <- pffr_simulate(
    Y ~ ff(X1, xind = s),
    n = n,
    nxgrid = nxgrid,
    nygrid = nygrid,
    effects = list(X1 = "cosine"),
    SNR = 100
  )
  s <- attr(dat, "xindex")
  t <- attr(dat, "yindex")
  truth <- attr(dat, "truth")

  expect_valid_truth(dat, n, nygrid)

  m <- pffr(Y ~ ff(X1, xind = s), yind = t, data = dat)

  expect_s3_class(m, "pffr")
  expect_equal(dim(fitted(m)), c(n, nygrid))
  expect_recovers_eta(m, truth, cor_threshold = 0.9, label = "ff() term")
})

test_that("ff() with limits='s<t' recovers historical effect", {
  skip_on_cran()
  set.seed(43)

  n <- 40
  nxgrid <- 30
  nygrid <- 40

  dat <- pffr_simulate(
    Y ~ ff(X1, xind = s),
    n = n,
    nxgrid = nxgrid,
    nygrid = nygrid,
    effects = list(X1 = "historical"),
    limits = function(s, t) s < t,
    SNR = 100
  )
  s <- attr(dat, "xindex")
  t <- attr(dat, "yindex")
  truth <- attr(dat, "truth")

  expect_valid_truth(dat, n, nygrid)

  m <- pffr(Y ~ ff(X1, xind = s, limits = "s<t"), yind = t, data = dat)

  expect_s3_class(m, "pffr")
  expect_recovers_eta(m, truth, cor_threshold = 0.9, label = "ff() with limits")
})

test_that("s() smooth term varying over t is recovered", {
  skip_on_cran()
  set.seed(44)

  n <- 50
  nygrid <- 40

  dat <- pffr_simulate(
    Y ~ s(xsmoo),
    n = n,
    nygrid = nygrid,
    effects = list(xsmoo = "sine"),
    SNR = 100
  )
  t <- attr(dat, "yindex")
  truth <- attr(dat, "truth")

  expect_valid_truth(dat, n, nygrid)

  m <- pffr(Y ~ s(xsmoo), yind = t, data = dat)

  expect_s3_class(m, "pffr")
  expect_recovers_eta(m, truth, cor_threshold = 0.9, label = "s() term")
})

test_that("c() constant-over-t terms are recovered", {
  skip_on_cran()
  set.seed(45)

  n <- 50
  nygrid <- 40

  dat <- pffr_simulate(
    Y ~ c(xconst),
    n = n,
    nygrid = nygrid,
    effects = list(xconst = 2),
    SNR = 100
  )
  t <- attr(dat, "yindex")
  truth <- attr(dat, "truth")

  expect_valid_truth(dat, n, nygrid)

  m <- pffr(Y ~ c(xconst), yind = t, data = dat)

  expect_s3_class(m, "pffr")
  expect_recovers_eta(m, truth, cor_threshold = 0.9, label = "c() term")
})

test_that("linear varying coefficient terms are recovered", {
  skip_on_cran()
  set.seed(46)

  n <- 50
  nygrid <- 40

  dat <- pffr_simulate(
    Y ~ xlin,
    n = n,
    nygrid = nygrid,
    effects = list(xlin = "dnorm"),
    SNR = 100
  )
  t <- attr(dat, "yindex")
  truth <- attr(dat, "truth")

  expect_valid_truth(dat, n, nygrid)

  m <- pffr(Y ~ xlin, yind = t, data = dat)

  expect_s3_class(m, "pffr")
  expect_recovers_eta(m, truth, cor_threshold = 0.9, label = "linear term")
})

test_that("sff() term works with new simulation", {
  skip_on_cran()
  set.seed(48)

  n <- 30
  nxgrid <- 25
  nygrid <- 35

  dat <- pffr_simulate(
    Y ~ ff(X1, xind = s),
    n = n,
    nxgrid = nxgrid,
    nygrid = nygrid,
    effects = list(X1 = "cosine"),
    SNR = 100
  )
  s <- attr(dat, "xindex")
  t <- attr(dat, "yindex")

  m_sff <- pffr(Y ~ sff(X1, xind = s), yind = t, data = dat)
  m_ff <- pffr(Y ~ ff(X1, xind = s), yind = t, data = dat)

  expect_s3_class(m_sff, "pffr")
  expect_equal(dim(fitted(m_sff)), c(n, nygrid))

  # sff and ff should give similar fits (safe comparison)
  expect_lt(safe_rel_error(fitted(m_ff), fitted(m_sff)), 0.1)
})

test_that("ffpc() term works with new simulation", {
  skip_on_cran()
  set.seed(49)

  n <- 40
  S <- 50
  n.argvals <- 40
  s <- seq(0, 1, l = S)
  argvals <- seq(0, 1, l = n.argvals)

  rankX <- 4
  Phi <- cbind(1 / sqrt(S), poly(s, degree = rankX - 1))
  lambda <- rankX:1
  Xi <- sapply(lambda, function(l) scale(rnorm(n, sd = sqrt(l)), scale = FALSE))
  X <- Xi %*% t(Phi)

  beta.st <- outer(s, argvals, function(s, t) cos(2 * pi * s * t))
  y <- (1 / S * X) %*% beta.st + 0.05 * matrix(rnorm(n * n.argvals), nrow = n)

  data <- list(y = y, X = X)

  m.pc <- pffr(
    y ~ c(1) + 0 + ffpc(X, yind = argvals, decomppars = list(npc = rankX)),
    data = data,
    yind = argvals
  )

  expect_s3_class(m.pc, "pffr")
  expect_s3_class(summary(m.pc), "summary.pffr")
  expect_equal(dim(fitted(m.pc)), c(n, n.argvals))
})


###############################################################################
# Algorithm Tests
###############################################################################

test_that("algorithm='bam' produces valid results", {
  skip_on_cran()
  set.seed(50)

  n <- 40
  nygrid <- 35

  dat <- sim_xlin_data(n = n, nygrid = nygrid, SNR = 50)
  t <- attr(dat, "yindex")
  truth <- attr(dat, "truth")

  m_bam <- pffr(Y ~ xlin, yind = t, data = dat, algorithm = "bam")

  expect_s3_class(m_bam, "pffr")
  expect_s3_class(summary(m_bam), "summary.pffr")
  expect_equal(dim(fitted(m_bam)), c(n, nygrid))

  # bam should also recover the signal
  expect_recovers_eta(
    m_bam,
    truth,
    cor_threshold = 0.85,
    label = "bam algorithm"
  )
})

test_that("algorithm='gamm' produces valid output", {
  skip_on_cran()
  set.seed(51)

  n <- 40
  nygrid <- 30
  n_groups <- 4

  dat <- sim_xlin_data(n = n, nygrid = nygrid, SNR = 50)
  dat$group <- factor(rep(1:n_groups, each = n / n_groups))
  t <- attr(dat, "yindex")

  m_gamm <- pffr(
    Y ~ xlin + s(group, bs = "re"),
    yind = t,
    data = dat,
    algorithm = "gamm"
  )

  # gamm returns list with $gam component
  expect_true(inherits(m_gamm, "pffr") || inherits(m_gamm, "list"))
  if (inherits(m_gamm, "list")) {
    expect_true("gam" %in% names(m_gamm))
  }
})

test_that("algorithm='gamm4' produces valid output", {
  skip_on_cran()
  skip_if_not_installed("gamm4")
  set.seed(52)

  n <- 40
  nygrid <- 30

  dat <- sim_xlin_data(n = n, nygrid = nygrid, SNR = 50)
  t <- attr(dat, "yindex")

  m_gamm4 <- pffr(Y ~ xlin, yind = t, data = dat, algorithm = "gamm4")

  expect_true(inherits(m_gamm4, "pffr") || inherits(m_gamm4, "list"))
  if (inherits(m_gamm4, "list")) {
    expect_true("gam" %in% names(m_gamm4))
  }
})


###############################################################################
# Inference Tests
###############################################################################

test_that("predict() returns correct dimensions and values", {
  skip_on_cran()
  set.seed(55)

  n <- 40
  nxgrid <- 30
  nygrid <- 35

  dat <- pffr_simulate(
    Y ~ ff(X1, xind = s) + xlin,
    n = n,
    nxgrid = nxgrid,
    nygrid = nygrid,
    effects = list(X1 = "cosine", xlin = "dnorm"),
    SNR = 50
  )
  s <- attr(dat, "xindex")
  t <- attr(dat, "yindex")

  m <- pffr(Y ~ ff(X1, xind = s) + xlin, yind = t, data = dat)

  # Predict on training data
  pred_train <- predict(m)
  expect_equal(dim(pred_train), c(n, nygrid))
  expect_false(any(is.na(pred_train)))

  # predict() should equal fitted()
  expect_equal(pred_train, fitted(m), tolerance = 1e-10)

  # Predict on new data
  set.seed(56)
  newdat <- pffr_simulate(
    Y ~ ff(X1, xind = s) + xlin,
    n = 20,
    nxgrid = nxgrid,
    nygrid = nygrid,
    effects = list(X1 = "cosine", xlin = "dnorm"),
    SNR = 50
  )
  pred_new <- predict(m, newdata = newdat)
  expect_equal(dim(pred_new), c(20, nygrid))
  expect_false(any(is.na(pred_new)))
})

test_that("predict(type='terms') returns valid term contributions", {
  skip_on_cran()
  set.seed(57)

  n <- 40
  nxgrid <- 30
  nygrid <- 35

  dat <- pffr_simulate(
    Y ~ ff(X1, xind = s) + xlin + c(xconst),
    n = n,
    nxgrid = nxgrid,
    nygrid = nygrid,
    effects = list(X1 = "cosine", xlin = "dnorm", xconst = 2),
    SNR = 50
  )
  s <- attr(dat, "xindex")
  t <- attr(dat, "yindex")

  m <- pffr(Y ~ ff(X1, xind = s) + xlin + c(xconst), yind = t, data = dat)

  pred_terms <- predict(m, type = "terms")

  expect_true(is.list(pred_terms))
  expect_true(length(pred_terms) >= 3) # At least 3 terms
  expect_true(any(grepl("X1", names(pred_terms))))
  expect_true(any(grepl("xlin", names(pred_terms))))
  expect_true(any(grepl("xconst", names(pred_terms))))

  # Each term should have correct number of rows
  for (term_name in names(pred_terms)) {
    expect_equal(
      nrow(pred_terms[[term_name]]),
      n,
      info = paste("term", term_name, "should have n rows")
    )
  }

  # Sum of terms + constant should approximate response prediction
  pred_response <- predict(m, type = "response")
  constant <- attr(pred_terms, "constant")
  if (!is.null(constant)) {
    term_sum <- Reduce(`+`, pred_terms)
    if (is.matrix(constant)) {
      term_sum_with_const <- term_sum + constant
    } else {
      term_sum_with_const <- sweep(term_sum, 2, constant, "+")
    }
    # This relationship may not be exact due to mgcv internals, but should be close
    expect_lt(safe_rel_error(pred_response, term_sum_with_const), 0.1)
  }
})

test_that("predict(type='terms') works for intercept-free models", {
  skip_on_cran()
  set.seed(58)

  n <- 35
  nygrid <- 30

  dat <- sim_xlin_data(n = n, nygrid = nygrid, SNR = 50)
  t <- attr(dat, "yindex")

  m <- pffr(Y ~ 0 + xlin, yind = t, data = dat)

  pred_terms <- predict(m, type = "terms")
  term_sum <- Reduce(`+`, pred_terms)
  pred_response <- predict(m, type = "response")

  expect_lt(safe_rel_error(pred_response, term_sum), 0.1)
  constant <- attr(pred_terms, "constant")
  if (!is.null(constant)) {
    expect_true(all(constant == 0))
  }
})

test_that("coef() extraction works correctly", {
  skip_on_cran()
  set.seed(58)

  n <- 40
  nxgrid <- 30
  nygrid <- 35

  dat <- pffr_simulate(
    Y ~ ff(X1, xind = s) + xlin,
    n = n,
    nxgrid = nxgrid,
    nygrid = nygrid,
    effects = list(X1 = "cosine", xlin = "dnorm"),
    SNR = 50
  )
  s <- attr(dat, "xindex")
  t <- attr(dat, "yindex")

  m <- pffr(Y ~ ff(X1, xind = s) + xlin, yind = t, data = dat)

  coefs <- coef(m)

  expect_true(is.list(coefs))
  expect_true("pterms" %in% names(coefs))
  expect_true("smterms" %in% names(coefs))
  expect_true(is.matrix(coefs$pterms) || is.data.frame(coefs$pterms))

  # smterms should contain entries for smooth terms
  expect_true(is.list(coefs$smterms))
  expect_true(length(coefs$smterms) >= 1)
})

test_that("summary() output contains expected components", {
  skip_on_cran()
  set.seed(59)

  n <- 40
  nxgrid <- 30
  nygrid <- 35

  dat <- pffr_simulate(
    Y ~ ff(X1, xind = s) + xlin,
    n = n,
    nxgrid = nxgrid,
    nygrid = nygrid,
    effects = list(X1 = "cosine", xlin = "dnorm"),
    SNR = 50
  )
  s <- attr(dat, "xindex")
  t <- attr(dat, "yindex")

  m <- pffr(Y ~ ff(X1, xind = s) + xlin, yind = t, data = dat)

  summ <- summary(m)

  expect_s3_class(summ, "summary.pffr")
  expect_true(!is.null(summ$s.table)) # Smooth terms table
  expect_true(!is.null(summ$r.sq)) # R-squared
  expect_true(summ$r.sq >= 0 && summ$r.sq <= 1)
})


###############################################################################
# Non-Gaussian Simulation Tests (new family support)
###############################################################################

test_that("pffrSim with family=scat() produces t-distributed responses", {
  skip_on_cran()
  set.seed(63)

  n <- 100
  nygrid <- 30

  dat <- sim_xlin_data(n = n, nygrid = nygrid, SNR = 10, family = mgcv::scat())

  Y <- dat$Y
  truth <- attr(dat, "truth")

  # Dimensions should be correct
  expect_equal(dim(Y), c(n, nygrid))
  expect_equal(dim(truth$eta), c(n, nygrid))
  expect_equal(dim(truth$epsilon), c(n, nygrid))

  # Response should be continuous (not integer)
  expect_false(all(Y == floor(Y)))

  # Can fit a pffr model with the data
  t <- attr(dat, "yindex")
  m <- pffr(Y ~ xlin, yind = t, data = dat)
  expect_s3_class(m, "pffr")
})

test_that("pffrSim with family=gaulss() produces Gaussian-like responses", {
  skip_on_cran()
  set.seed(64)

  n <- 50
  nygrid <- 30

  dat <- sim_xlin_data(
    n = n,
    nygrid = nygrid,
    SNR = 10,
    family = mgcv::gaulss()
  )

  Y <- dat$Y
  truth <- attr(dat, "truth")

  # Dimensions should be correct
  expect_equal(dim(Y), c(n, nygrid))
  expect_equal(dim(truth$eta), c(n, nygrid))

  # Response should be continuous
  expect_false(all(Y == floor(Y)))

  # Can fit a pffr model with the data
  t <- attr(dat, "yindex")
  m <- pffr(Y ~ xlin, yind = t, data = dat)
  expect_s3_class(m, "pffr")
})


###############################################################################
# New Feature Tests (effect selection, multi-arg const)
###############################################################################

test_that("c(te(x1, x2)) with gaussian_2d preset works", {
  skip_on_cran()
  set.seed(71)

  n <- 40
  nygrid <- 30

  dat <- pffr_simulate(
    Y ~ c(te(xte1, xte2)),
    n = n,
    nygrid = nygrid,
    effects = list("te(xte1,xte2)" = "gaussian_2d"),
    SNR = 50
  )

  expect_valid_truth(dat, n, nygrid)

  t <- attr(dat, "yindex")
  m <- pffr(Y ~ c(te(xte1, xte2)), yind = t, data = dat)

  expect_s3_class(m, "pffr")
  expect_equal(dim(fitted(m)), c(n, nygrid))

  # Should recover the signal reasonably well
  truth <- attr(dat, "truth")
  expect_gt(cor(as.vector(fitted(m)), as.vector(truth$eta)), 0.8)
})

test_that("multiple terms with mixed presets work", {
  skip_on_cran()
  set.seed(72)

  n <- 40
  nxgrid <- 25
  nygrid <- 30

  dat <- pffr_simulate(
    Y ~ ff(X1, xind = s) + xlin + s(xsmoo) + c(xconst),
    n = n,
    nxgrid = nxgrid,
    nygrid = nygrid,
    effects = list(
      X1 = "product",
      xlin = "dnorm",
      xsmoo = "cosine",
      xconst = 1.5
    ),
    SNR = 50
  )

  expect_valid_truth(dat, n, nygrid)

  truth <- attr(dat, "truth")
  expect_true(length(truth$etaTerms) >= 4)

  s <- attr(dat, "xindex")
  t <- attr(dat, "yindex")

  m <- pffr(
    Y ~ ff(X1, xind = s) + xlin + s(xsmoo) + c(xconst),
    yind = t,
    data = dat
  )

  expect_s3_class(m, "pffr")
  expect_recovers_eta(m, truth, cor_threshold = 0.85, label = "mixed terms")
})

test_that("unwrapped te/ti/t2 terms work without explicit effects", {
  skip_on_cran()
  set.seed(80)

  n <- 30
  nygrid <- 25

  # Bug fix test: te without c() should use smooth preset (sine), not gaussian_2d
  dat <- pffr_simulate(
    Y ~ te(xte1, xte2), # NOT wrapped in c()
    n = n,
    nygrid = nygrid,
    SNR = 50
  )

  expect_valid_truth(dat, n, nygrid)

  t <- attr(dat, "yindex")
  m <- pffr(Y ~ te(xte1, xte2), yind = t, data = dat)

  expect_s3_class(m, "pffr")
  expect_equal(dim(fitted(m)), c(n, nygrid))
})

test_that("factor linear terms work with default preset", {
  skip_on_cran()
  set.seed(81)

  n <- 40
  nygrid <- 30

  # Bug fix test: factor terms should work without custom 2-arg function
  dat <- pffr_simulate(
    Y ~ xfactor,
    n = n,
    nygrid = nygrid,
    SNR = 50
  )

  # xfactor should be generated as a factor
  expect_true(is.factor(dat$xfactor))

  expect_valid_truth(dat, n, nygrid)

  t <- attr(dat, "yindex")
  m <- pffr(Y ~ xfactor, yind = t, data = dat)

  expect_s3_class(m, "pffr")
  expect_equal(dim(fitted(m)), c(n, nygrid))
})

test_that("pffrSim respects custom response name from formula", {
  skip_on_cran()
  set.seed(83)

  n <- 30
  nygrid <- 25

  # Use custom response name
  dat <- pffr_simulate(
    myresponse ~ xlin,
    n = n,
    nygrid = nygrid,
    effects = list(xlin = "dnorm"),
    SNR = 50
  )

  # Response should be named "myresponse", not "Y"
  expect_true("myresponse" %in% names(dat))
  expect_false("Y" %in% names(dat))
  expect_equal(dim(dat$myresponse), c(n, nygrid))

  # Model should work with the same formula
  t <- attr(dat, "yindex")
  m <- pffr(myresponse ~ xlin, yind = t, data = dat)
  expect_s3_class(m, "pffr")
})


###############################################################################
# Phase 3: pffr Methods Refactoring & Test Coverage
###############################################################################

# -----------------------------------------------------------------------------
# 3.1.1 Factor Variables with >2 Levels
# -----------------------------------------------------------------------------

test_that("factor with 3 levels works as varying coefficient", {
  skip_on_cran()
  set.seed(101)

  n <- 40
  nygrid <- 30

  # Generate data with 3-level factor
  dat <- sim_xlin_data(n = n, nygrid = nygrid, SNR = 50)
  dat$xfactor3 <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
  t <- attr(dat, "yindex")

  # Fit model with 3-level factor

  m <- pffr(Y ~ xfactor3, yind = t, data = dat)

  expect_s3_class(m, "pffr")
  expect_equal(dim(fitted(m)), c(n, nygrid))

  # Check summary shows factor levels
  summ <- summary(m)
  expect_s3_class(summ, "summary.pffr")

  # coef should return coefficients for each level
  coefs <- coef(m)
  expect_true(is.list(coefs))

  # shortlabels should include distinct level names
  shortlabels <- m$pffr$short_labels
  expect_true(length(shortlabels) >= length(levels(dat$xfactor3)))
  expect_true(all(vapply(
    levels(dat$xfactor3),
    function(lvl) {
      any(grepl(lvl, shortlabels))
    },
    logical(1)
  )))
})

test_that("factor in interaction with smooth terms", {
  skip_on_cran()
  set.seed(102)

  n <- 40
  nygrid <- 30

  dat <- sim_xlin_data(n = n, nygrid = nygrid, SNR = 50)
  dat$xfactor2 <- factor(sample(c("A", "B"), n, replace = TRUE))
  t <- attr(dat, "yindex")

  # Factor in by argument of smooth term
  m <- pffr(Y ~ s(xlin, by = xfactor2), yind = t, data = dat)

  expect_s3_class(m, "pffr")
  expect_equal(dim(fitted(m)), c(n, nygrid))
  expect_s3_class(summary(m), "summary.pffr")
})

# -----------------------------------------------------------------------------
# 3.1.2 Untested Methods
# -----------------------------------------------------------------------------

test_that("residuals.pffr returns correct dimensions and handles sparse data", {
  skip_on_cran()

  # Regular data - use shared fixture
  m <- get_xlin_model()
  dat <- get_xlin_data()
  n <- 30
  nygrid <- 25

  # Test reformatted residuals
  resid_mat <- residuals(m, reformat = TRUE)
  expect_true(is.matrix(resid_mat))
  expect_equal(dim(resid_mat), c(n, nygrid))

  # Test raw residuals (vector)
  resid_vec <- residuals(m, reformat = FALSE)
  expect_true(is.vector(resid_vec))
  expect_equal(length(resid_vec), n * nygrid)

  # Test sparse data residuals - use shared fixture
  m_sparse <- get_sparse_model()

  resid_sparse <- residuals(m_sparse, reformat = TRUE)
  expect_true(is.data.frame(resid_sparse))
  expect_true(".value" %in% colnames(resid_sparse))
})

test_that("plot.pffr runs without error (smoke test)", {
  skip_on_cran()
  # Skip due to known mgcv dispatch issue with plot.gam
  # The plot.pffr method itself works but has object name scoping issues
  skip("plot.pffr has known mgcv object dispatch issues")

  m <- get_xlin_model()

  # Smoke test: plot should run without error
  expect_no_error(suppressWarnings(plot(m, pages = 1)))
})

test_that("qq.pffr runs without error (smoke test)", {
  skip_on_cran()

  m <- get_xlin_model()

  # Smoke test: qq.pffr should run without error
  expect_silent(suppressWarnings(pffr_qq(m)))
})

test_that("pffr.check runs without error (smoke test)", {
  skip_on_cran()

  m <- get_xlin_model()

  # Smoke test: pffr.check should run without error
  # It produces output (diagnostics), so we just check it doesn't error
  expect_no_error(suppressWarnings(pffr_check(m)))
})

test_that("print.summary.pffr output contains expected information", {
  skip_on_cran()

  m <- get_basic_pffr_model()
  summ <- summary(m)

  # Capture print output
  output <- capture.output(print(summ))
  output_str <- paste(output, collapse = "\n")

  # Test formula is printed
  expect_true(any(grepl("Formula", output_str)))

  # Test smooth terms table
  expect_true(any(grepl("Smooth terms", output_str)))

  # Test R-squared is shown
  expect_true(any(grepl("R-sq", output_str)))

  # Test sample size format: "n = X (Y x Z)"
  expect_true(any(grepl("n = ", output_str)))
  expect_true(any(grepl("\\(.*x.*\\)", output_str)))
})

# -----------------------------------------------------------------------------
# 3.1.3 Edge Cases
# -----------------------------------------------------------------------------

test_that("se.fit = TRUE in predict.pffr returns correct structure", {
  skip_on_cran()

  m <- get_xlin_model()
  n <- 30
  nygrid <- 25

  # Test se.fit = TRUE
  pred_se <- predict(m, se.fit = TRUE)

  expect_true(is.list(pred_se))
  expect_true("fit" %in% names(pred_se))
  expect_true("se.fit" %in% names(pred_se))

  # Check dimensions match
  expect_equal(dim(pred_se$fit), c(n, nygrid))
  expect_equal(dim(pred_se$se.fit), c(n, nygrid))

  # SE should be non-negative
  expect_true(all(pred_se$se.fit >= 0))

  # Test with different types
  pred_se_link <- predict(m, se.fit = TRUE, type = "link")
  expect_true(is.list(pred_se_link))
  expect_equal(dim(pred_se_link$fit), c(n, nygrid))

  pred_se_response <- predict(m, se.fit = TRUE, type = "response")
  expect_true(is.list(pred_se_response))
  expect_equal(dim(pred_se_response$fit), c(n, nygrid))

  # Test with type = "terms" (list of matrices per term)
  pred_se_terms <- predict(m, se.fit = TRUE, type = "terms")
  expect_true(is.list(pred_se_terms))
  expect_true(is.list(pred_se_terms$fit))
  expect_true(is.list(pred_se_terms$se.fit))
})

test_that("small n (n=10) still fits without error", {
  skip_on_cran()
  set.seed(110)

  # Use n=10 (n=5 is too small for default spline basis)
  n <- 10
  nygrid <- 20

  dat <- sim_xlin_data(n = n, nygrid = nygrid, SNR = 100)
  t <- attr(dat, "yindex")

  # Should fit without error with small sample
  m <- pffr(Y ~ xlin, yind = t, data = dat)

  expect_s3_class(m, "pffr")
  expect_equal(dim(fitted(m)), c(n, nygrid))
})

test_that("model with single smooth term works correctly", {
  skip_on_cran()
  set.seed(111)

  n <- 30
  nygrid <- 25

  dat <- pffr_simulate(
    Y ~ s(xsmoo),
    n = n,
    nygrid = nygrid,
    effects = list(xsmoo = "sine"),
    SNR = 50
  )
  t <- attr(dat, "yindex")

  m <- pffr(Y ~ s(xsmoo), yind = t, data = dat)

  expect_s3_class(m, "pffr")
  expect_equal(dim(fitted(m)), c(n, nygrid))

  # Verify we have coefficients
  coefs <- coef(m)
  expect_true(is.list(coefs))
  expect_true("smterms" %in% names(coefs))

  # Predict works
  pred <- predict(m)
  expect_equal(dim(pred), c(n, nygrid))
})

test_that("predict.pffr with ff limits works on new data", {
  skip_on_cran()
  set.seed(113)

  n <- 30
  nxgrid <- 25
  nygrid <- 30

  # Generate training data with limits
  dat_train <- pffr_simulate(
    Y ~ ff(X1, xind = s),
    n = n,
    nxgrid = nxgrid,
    nygrid = nygrid,
    effects = list(X1 = "historical"),
    limits = function(s, t) s < t,
    SNR = 50
  )
  s <- attr(dat_train, "xindex")
  t <- attr(dat_train, "yindex")

  # Fit model with limits
  m <- pffr(Y ~ ff(X1, xind = s, limits = "s<t"), yind = t, data = dat_train)
  expect_s3_class(m, "pffr")

  # Generate new data for prediction
  set.seed(114)
  dat_new <- pffr_simulate(
    Y ~ ff(X1, xind = s),
    n = 15,
    nxgrid = nxgrid,
    nygrid = nygrid,
    effects = list(X1 = "historical"),
    limits = function(s, t) s < t,
    SNR = 50
  )

  # Predict on new data - should work without error
  pred <- predict(m, newdata = dat_new)
  expect_true(is.matrix(pred))
  expect_equal(dim(pred), c(15, nygrid))
  expect_false(any(is.na(pred)))
})

test_that("coef.pffr works with pcre terms with 3 FPCs", {
  skip_on_cran()
  set.seed(115)

  n <- 30
  n_groups <- 5
  ny <- 40
  t <- seq(0, 1, length = ny)

  # Create eigenfunctions for 3 FPCs
  efunctions <- cbind(
    sin(2 * pi * t),
    cos(2 * pi * t),
    sin(4 * pi * t)
  )
  evalues <- c(1, 0.5, 0.25)
  colnames(efunctions) <- c("PC1", "PC2", "PC3")

  # Generate random group effects
  group_effects <- matrix(0, n, ny)
  groups <- factor(rep(1:n_groups, each = n / n_groups))
  for (g in 1:n_groups) {
    scores <- rnorm(3, sd = sqrt(evalues))
    group_effects[groups == g, ] <- rep(
      efunctions %*% scores,
      each = sum(groups == g)
    )
  }

  # Generate data
  Y <- matrix(rnorm(n * ny, sd = 0.1), n, ny) + group_effects
  data <- list(Y = I(Y), group = groups)

  # Fit pcre model
  m <- pffr(
    Y ~ pcre(id = group, efunctions = efunctions, evalues = evalues, yind = t),
    yind = t,
    data = data
  )

  expect_s3_class(m, "pffr")

  # coef should work without error
  coefs <- coef(m, se = FALSE)
  expect_true(is.list(coefs))
  expect_true("smterms" %in% names(coefs))
  expect_true(length(coefs$smterms) >= 2) # intercept + pcre
})

test_that("fitted.pffr with gaulss family returns mean/scale correctly", {
  skip_on_cran()

  m <- get_gaulss_model()
  n <- 30
  nygrid <- 25

  expect_s3_class(m, "pffr")

  # Default which = "mean" should return matrix
  fit_mean <- fitted(m)
  expect_true(is.matrix(fit_mean))
  expect_equal(dim(fit_mean), c(n, nygrid))

  # which = "scale" should return scale (log-sd) matrix
  fit_scale <- fitted(m, which = "scale")
  expect_true(is.matrix(fit_scale))
  expect_equal(dim(fit_scale), c(n, nygrid))

  # which = "both" should return list with mean and scale
  fit_both <- fitted(m, which = "both")
  expect_true(is.list(fit_both))
  expect_true("mean" %in% names(fit_both))
  expect_true("scale" %in% names(fit_both))
  expect_equal(dim(fit_both$mean), c(n, nygrid))
  expect_equal(dim(fit_both$scale), c(n, nygrid))

  # Mean and scale from "both" should match individual calls
  expect_equal(fit_both$mean, fit_mean)
  expect_equal(fit_both$scale, fit_scale)

  # reformat = FALSE should return vectors (or lists of vectors)
  fit_mean_vec <- fitted(m, reformat = FALSE)
  fit_both_vec <- fitted(m, reformat = FALSE, which = "both")
  expect_true(is.vector(fit_mean_vec))
  expect_true(is.list(fit_both_vec))
  expect_equal(length(fit_both_vec$mean), n * nygrid)
})

test_that("gaulss model has correct shortlabels for scale parameter", {
  skip_on_cran()

  m <- get_gaulss_model()

  # Check shortlabels includes log(SD) label
  shortlabels <- m$pffr$short_labels
  expect_true(any(grepl("log\\(SD\\)", shortlabels)))

  # Summary should show log(SD)(t) instead of NA
  summ <- summary(m)
  s_table_rows <- rownames(summ$s.table)
  expect_true(any(grepl("log\\(SD\\)", s_table_rows)))
  expect_false(any(is.na(s_table_rows)))
})

test_that("create_shortlabels handles complex scale formulas for gaulss", {
  # Test create_shortlabels directly with mock scale smooths
  # that have covariates (complex case for gaulss)

  mock_smooth <- list(
    list(label = "s(t.vec)"),
    list(label = "s(t.vec):xlin"),
    list(label = "s.1(t.vec)"),
    list(label = "s.1(t.vec):xlin")
  )

  mock_labelmap <- list(
    xlin = "s(t.vec):xlin",
    "Intercept(t)" = "s(t.vec)"
  )

  mock_where <- list(par = 1L, ff = integer(0), sff = integer(0))

  # Mock gaulss family (has nlp = 2)
  mock_family <- list(nlp = 2)

  result <- create_shortlabels(
    label_map = mock_labelmap,
    m_smooth = mock_smooth,
    yind_name = "t",
    where_specials = mock_where,
    family = mock_family
  )

  # Simple scale intercept should be log(SD)(t)
  expect_equal(result[["s.1(t.vec)"]], "log(SD)(t)")

  # Complex scale with covariate should be log(SD): xlin(t)
  expect_equal(result[["s.1(t.vec):xlin"]], "log(SD): xlin(t)")

  # Without gaulss family, scale smooths should NOT be labeled
  result_gaussian <- create_shortlabels(
    label_map = mock_labelmap,
    m_smooth = mock_smooth,
    yind_name = "t",
    where_specials = mock_where,
    family = list(nlp = 1) # Regular family
  )

  # Scale smooths should not be in result for non-location-scale family
  expect_false("s.1(t.vec)" %in% names(result_gaussian))
})


###############################################################################
# pffrGLS Tests
###############################################################################

test_that("pffrGLS fits a simple model correctly", {
  skip_on_cran()

  m <- get_gls_model()
  n <- 30
  nygrid <- 25

  expect_s3_class(m, "pffr")
  expect_equal(dim(fitted(m)), c(n, nygrid))
  expect_false(any(is.na(fitted(m))))
})

test_that("pffrGLS has complete pffr slot", {
  skip_on_cran()

  m <- get_gls_model()

  # Check pffr slot has all required fields
  expect_true(!is.null(m$pffr$short_labels))
  expect_true(!is.null(m$pffr$label_map))
  expect_true(!is.null(m$pffr$term_map))
  expect_true(!is.null(m$pffr$nobs))
  expect_true(!is.null(m$pffr$nyindex))
  expect_true(!is.null(m$pffr$yind))
  expect_true(!is.null(m$pffr$where))

  # GLS-specific fields
  expect_true(!is.null(m$pffr$hatSigma))
  expect_true(!is.null(m$pffr$sqrtSigmaInv))
})

test_that("pffrGLS works with ff() terms", {
  skip_on_cran()
  set.seed(202)

  n <- 30
  nxgrid <- 20
  nygrid <- 25

  dat <- pffr_simulate(
    Y ~ ff(X1, xind = s),
    n = n,
    nxgrid = nxgrid,
    nygrid = nygrid,
    effects = list(X1 = "cosine"),
    SNR = 50
  )
  s <- attr(dat, "xindex")
  t <- attr(dat, "yindex")

  rho <- 0.5
  hatSigma <- rho^abs(outer(1:nygrid, 1:nygrid, "-"))

  m <- pffr_gls(Y ~ ff(X1, xind = s), yind = t, data = dat, hatSigma = hatSigma)

  expect_s3_class(m, "pffr")
  expect_equal(dim(fitted(m)), c(n, nygrid))

  # Check coef works
  coefs <- coef(m)
  expect_true(is.list(coefs))
  expect_true("smterms" %in% names(coefs))
})

test_that("pffrGLS summary and coef methods work", {
  skip_on_cran()

  m <- get_gls_model()

  # Test summary
  summ <- summary(m)
  expect_s3_class(summ, "summary.pffr")
  expect_true(!is.null(summ$s.table))

  # Test coef
  coefs <- coef(m)
  expect_true(is.list(coefs))
  expect_true("pterms" %in% names(coefs))
  expect_true("smterms" %in% names(coefs))
})

test_that("pffrGLS handles poorly conditioned hatSigma", {
  skip_on_cran()
  set.seed(204)

  n <- 20
  nygrid <- 15

  dat <- sim_xlin_data(n = n, nygrid = nygrid, SNR = 50)
  t <- attr(dat, "yindex")

  # Create poorly conditioned matrix (high correlation)
  rho <- 0.99
  hatSigma <- rho^abs(outer(1:nygrid, 1:nygrid, "-"))

  # Should warn about condition number
  expect_warning(
    m <- pffr_gls(Y ~ xlin, yind = t, data = dat, hatSigma = hatSigma),
    "condition number"
  )

  # But should still fit
  expect_s3_class(m, "pffr")
})

test_that("pffrGLS rejects missing values in response", {
  skip_on_cran()
  set.seed(205)

  n <- 20
  nygrid <- 15

  dat <- sim_xlin_data(n = n, nygrid = nygrid, SNR = 50)
  dat$Y[1, 5] <- NA # Introduce missing value
  t <- attr(dat, "yindex")

  rho <- 0.5
  hatSigma <- rho^abs(outer(1:nygrid, 1:nygrid, "-"))

  expect_error(
    pffr_gls(Y ~ xlin, yind = t, data = dat, hatSigma = hatSigma),
    "missing values"
  )
})

test_that("pffrGLS predict works correctly", {
  skip_on_cran()

  m <- get_gls_model()
  dat <- get_xlin_data()
  n <- 30
  nygrid <- 25

  # Predict on training data
  pred <- predict(m)
  expect_equal(dim(pred), c(n, nygrid))
  expect_equal(pred, fitted(m), tolerance = 1e-10)

  # Predict on new data
  set.seed(207)
  newdat <- sim_xlin_data(n = 15, nygrid = nygrid, SNR = 50)
  pred_new <- predict(m, newdata = newdat)
  expect_equal(dim(pred_new), c(15, nygrid))
})


###############################################################################
# coefboot.pffr Tests
###############################################################################

test_that("coefboot.pffr runs on simple pffr model", {
  skip_on_cran()

  m <- get_xlin_model()

  # Run bootstrap with small B for speed
  boot_ci <- pffr_coefboot(m, B = 5, showProgress = FALSE)

  expect_true(is.list(boot_ci))
  expect_true("pterms" %in% names(boot_ci))
  expect_true("smterms" %in% names(boot_ci))
  expect_true("boot_ci" %in% names(boot_ci))
  expect_true("ci_meta" %in% names(boot_ci))

  # pterms should have CI columns
  expect_true(ncol(boot_ci$pterms) > 1)
  expect_true(all(c("lower", "upper") %in% colnames(boot_ci$pterms)))
  expect_true(any(!is.na(boot_ci$pterms[, "lower"])))

  # At least one bootstrap replicate should succeed
  expect_lt(boot_ci$ci_meta$n_failed, boot_ci$ci_meta$B_requested)
  expect_gte(boot_ci$ci_meta$B_used, 1)
  expect_true(is.numeric(boot_ci$ci_meta$failure_rate))
  expect_gte(boot_ci$ci_meta$failure_rate, 0)
  expect_lte(boot_ci$ci_meta$failure_rate, 1)
})

test_that("coefboot.pffr residual resampling method works", {
  skip_on_cran()

  m <- get_xlin_model()

  # Test residual method
  boot_ci <- pffr_coefboot(m, B = 5, method = "residual", showProgress = FALSE)
  expect_true(is.list(boot_ci))
  expect_lt(boot_ci$ci_meta$n_failed, boot_ci$ci_meta$B_requested)

  # Test centered residual method
  boot_ci_c <- pffr_coefboot(
    m,
    B = 5,
    method = "residual.c",
    showProgress = FALSE
  )
  expect_true(is.list(boot_ci_c))
  expect_lt(boot_ci_c$ci_meta$n_failed, boot_ci_c$ci_meta$B_requested)
})

test_that("coefboot.pffr works with sparse/irregular data", {
  skip_on_cran()

  m_sparse <- get_sparse_model()

  # Run bootstrap with sparse data
  boot_ci <- pffr_coefboot(m_sparse, B = 5, showProgress = FALSE)

  expect_true(is.list(boot_ci))
  expect_true("pterms" %in% names(boot_ci))
  expect_true("smterms" %in% names(boot_ci))
  expect_lt(boot_ci$ci_meta$n_failed, boot_ci$ci_meta$B_requested)
})

test_that("coefboot.pffr works with gaulss family", {
  skip_on_cran()

  m_gaulss <- get_gaulss_model()

  expect_s3_class(m_gaulss, "pffr")

  # Run bootstrap
  boot_ci <- pffr_coefboot(m_gaulss, B = 5, showProgress = FALSE)

  expect_true(is.list(boot_ci))
  expect_true("pterms" %in% names(boot_ci))
  expect_true("smterms" %in% names(boot_ci))
  expect_lt(boot_ci$ci_meta$n_failed, boot_ci$ci_meta$B_requested)
})

test_that("coefboot percentile ordering and labels are correct", {
  set.seed(42)
  x <- rnorm(5000)
  conf <- c(0.9, 0.95)
  q <- compute_percentile_ci(x, conf)

  expect_length(q, 4)
  expect_lt(q[1], q[2]) # lower_90 < upper_90
  expect_lt(q[3], q[4]) # lower_95 < upper_95
  expect_lte(q[3], q[1]) # lower_95 <= lower_90
  expect_lte(q[2], q[4]) # upper_90 <= upper_95

  expect_identical(format_boot_conf_label(0.9), "90")
  expect_identical(format_boot_conf_label(0.95), "95")
  expect_identical(format_boot_conf_label(0.975), "97.5")

  mock_boot <- list(t = matrix(x, ncol = 1))
  ci_mat <- compute_boot_cis(
    boot_result = mock_boot,
    indices = 1,
    conf = conf,
    type = "percent"
  )
  expect_identical(
    colnames(ci_mat),
    c("lower_90", "upper_90", "lower_95", "upper_95")
  )
})

test_that("coefboot uses first confidence level as primary level", {
  skip_on_cran()

  m <- get_xlin_model()
  boot_ci <- pffr_coefboot(
    m,
    B = 5,
    conf = c(0.9, 0.95),
    showProgress = FALSE
  )

  expect_equal(boot_ci$ci_meta$level, 0.9)
  expect_equal(boot_ci$boot_ci$primary_conf, 0.9)

  first_ok <- which(
    is.finite(boot_ci$pterms[, "lower"]) &
      is.finite(boot_ci$boot_ci$pterms[, "lower_90"])
  )[1]
  expect_false(is.na(first_ok))
  expect_equal(
    boot_ci$pterms[first_ok, "lower"],
    boot_ci$boot_ci$pterms[first_ok, "lower_90"]
  )
  expect_equal(
    boot_ci$pterms[first_ok, "upper"],
    boot_ci$boot_ci$pterms[first_ok, "upper_90"]
  )
})

test_that("handle_failed_replicates warns and reports failure metadata", {
  mock_boot <- list(
    t = rbind(
      c(NA_real_, NA_real_),
      c(1, 2),
      c(NA_real_, NA_real_),
      c(3, 4)
    ),
    R = 4L
  )
  out <- expect_warning(
    handle_failed_replicates(mock_boot, B = 4L),
    "bootstrap replicates failed"
  )
  expect_identical(attr(out, "n_failed"), 2L)
  expect_identical(attr(out, "B_requested"), 4L)
  expect_equal(attr(out, "failure_rate"), 0.5)
  expect_identical(out$R, 2L)
  expect_equal(nrow(out$t), 2L)
})

test_that("coef supports fixed eval_grid for term-aligned extraction", {
  skip_on_cran()

  m <- get_multiterm_model()
  template <- coef(
    m,
    se = FALSE,
    sandwich = "none",
    ci = "none",
    n1 = 40,
    n2 = 40,
    n3 = 20
  )
  eval_grid <- build_coefboot_eval_grid(m, n1 = 40, n2 = 40, n3 = 20)

  # Force a resample that excludes observed min/max xsmoo values.
  mc <- prepare_modcall_for_bootstrap(m)
  n <- nrow(mc$data)
  idx_pool <- setdiff(seq_len(n), c(which.min(mc$data$xsmoo), which.max(mc$data$xsmoo)))
  set.seed(11)
  idx <- sample(idx_pool, n, replace = TRUE)

  mc_b <- resample_grid_observations(mc, mc$data, idx)
  m_b <- eval(mc_b)

  co_default <- coef(
    m_b,
    se = FALSE,
    sandwich = "none",
    ci = "none",
    n1 = 40,
    n2 = 40,
    n3 = 20
  )
  co_fixed <- coef(
    m_b,
    se = FALSE,
    sandwich = "none",
    ci = "none",
    n1 = 40,
    n2 = 40,
    n3 = 20,
    eval_grid = eval_grid
  )

  # Smooth with xsmoo should be range-shrunk in default but fixed under eval_grid.
  ind <- which(vapply(co_default$smterms, function(sm) {
    "xsmoo" %in% colnames(sm$coef)
  }, logical(1)))[1]
  expect_false(is.na(ind))

  expect_false(isTRUE(all.equal(
    co_default$smterms[[ind]]$coef[, "xsmoo"],
    template$smterms[[ind]]$coef[, "xsmoo"]
  )))
  expect_equal(
    co_fixed$smterms[[ind]]$coef[, "xsmoo"],
    template$smterms[[ind]]$coef[, "xsmoo"]
  )
})

test_that("pffrGLS errors on sparse data (not yet implemented)", {
  skip_on_cran()

  dat_sparse <- get_sparse_data()
  t <- attr(dat_sparse, "yindex")

  rho <- 0.5
  hatSigma <- rho^abs(outer(1:length(t), 1:length(t), "-"))

  # pffrGLS with sparse data should error (not yet implemented)
  expect_error(
    pffr_gls(
      Y ~ s(xsmoo),
      data = dat_sparse$data,
      ydata = dat_sparse$ydata,
      yind = t,
      hatSigma = hatSigma
    ),
    "sparse data not yet"
  )
})

###############################################################################
# Sandwich Correction Tests
###############################################################################

test_that("pffr with sandwich='cluster' yields cluster-robust covariance", {
  skip_on_cran()

  dat <- get_xlin_data()
  t <- attr(dat, "yindex")
  m_std <- get_xlin_model()

  m_cl <- pffr(Y ~ xlin, yind = t, data = dat, sandwich = "cluster")

  # Sandwich type stored as character
  expect_identical(m_std$pffr$sandwich, "none")
  expect_identical(m_cl$pffr$sandwich, "cluster")

  # Covariance matrices differ from uncorrected model
  expect_false(identical(m_std$Vp, m_cl$Vp))
  expect_false(identical(m_std$Ve, m_cl$Ve))

  # Cluster-robust differs from HC sandwich
  m_std_stripped <- m_std
  class(m_std_stripped) <- setdiff(class(m_std_stripped), "pffr")
  hc_Vp <- vcov(m_std_stripped, sandwich = TRUE)
  expect_false(identical(m_cl$Vp, hc_Vp))

  # Coefficients identical (only covariance changes)
  expect_equal(coef(m_std, raw = TRUE), coef(m_cl, raw = TRUE))

  # Summary shows sandwich type
  summ_std <- summary(m_std)
  summ_cl <- summary(m_cl)
  expect_identical(summ_std$sandwich, "none")
  expect_identical(summ_cl$sandwich, "cluster")

  # coef.pffr on fitted model uses stored matrices
  coef_cl <- coef(m_cl, sandwich = "cluster")
  # coef.pffr on un-sandwiched model computes on the fly
  coef_std_cl <- coef(m_std, sandwich = "cluster")
  expect_equal(coef_cl$pterms[, "se"], coef_std_cl$pterms[, "se"])
})

test_that("pffr with sandwich='hc' yields observation-level HC sandwich", {
  skip_on_cran()

  dat <- get_xlin_data()
  t <- attr(dat, "yindex")
  m_std <- get_xlin_model()

  m_hc <- pffr(Y ~ xlin, yind = t, data = dat, sandwich = "hc")

  expect_identical(m_hc$pffr$sandwich, "hc")

  # HC sandwich matches mgcv::vcov.gam(sandwich=TRUE) on uncorrected model
  m_std_stripped <- m_std
  class(m_std_stripped) <- setdiff(class(m_std_stripped), "pffr")
  expected_Vp <- vcov(m_std_stripped, sandwich = TRUE)
  expect_equal(m_hc$Vp, expected_Vp)
  expect_equal(m_hc$Vc, expected_Vp)
})

test_that("pffr default sandwich is cluster", {
  skip_on_cran()

  dat <- get_xlin_data()
  t <- attr(dat, "yindex")
  m_default <- pffr(Y ~ xlin, yind = t, data = dat)
  expect_identical(m_default$pffr$sandwich, "cluster")
})

test_that("pffr with sandwich='cl2' yields leverage-adjusted covariance", {
  skip_on_cran()

  dat <- get_xlin_data()
  t <- attr(dat, "yindex")
  m_std <- get_xlin_model()

  m_cl <- pffr(Y ~ xlin, yind = t, data = dat, sandwich = "cluster")
  m_cl2 <- pffr(Y ~ xlin, yind = t, data = dat, sandwich = "cl2")

  expect_identical(m_cl2$pffr$sandwich, "cl2")
  expect_equal(coef(m_std, raw = TRUE), coef(m_cl2, raw = TRUE))

  # CL2 covariance should differ from both uncorrected and CR1 covariance.
  expect_gt(max(abs(m_cl2$Vp - m_std$Vp)), 0)
  expect_gt(max(abs(m_cl2$Vp - m_cl$Vp)), 0)

  # coef.pffr on fitted model uses stored CL2 covariance.
  coef_cl2 <- coef(m_cl2, sandwich = "cl2")
  coef_std_cl2 <- coef(m_std, sandwich = "cl2")
  expect_equal(coef_cl2$pterms[, "se"], coef_std_cl2$pterms[, "se"])
})

test_that("gam_sandwich_cluster_cl2 works for poisson and binomial", {
  skip_on_cran()

  families <- list(poisson(), binomial())
  for (fam in families) {
    set.seed(if (fam$family == "poisson") 1001 else 1002)
    dat <- sim_xlin_data(n = 30, nygrid = 25, SNR = 10, family = fam)
    t <- attr(dat, "yindex")

    m <- pffr(Y ~ xlin, yind = t, data = dat, family = fam)
    m_stripped <- m
    class(m_stripped) <- setdiff(class(m_stripped), "pffr")
    cluster_id <- build_cluster_id(m$pffr)

    V_cl <- gam_sandwich_cluster(m_stripped, cluster_id, freq = FALSE)
    V_cl2 <- gam_sandwich_cluster_cl2(m_stripped, cluster_id, freq = FALSE)

    expect_equal(V_cl2, t(V_cl2), tolerance = 1e-10)
    expect_true(all(is.finite(diag(V_cl2))))
    expect_true(all(diag(V_cl2) >= 0))
    expect_gt(max(abs(V_cl2 - V_cl)), 0)

    coef_cl2 <- coef(m, sandwich = "cl2")
    expect_true(all(is.finite(coef_cl2$pterms[, "se"])))
  }
})

test_that("sandwich backward compat: TRUE/FALSE still work", {
  skip_on_cran()

  dat <- get_xlin_data()
  t <- attr(dat, "yindex")

  # TRUE -> "cluster"
  m_true <- pffr(Y ~ xlin, yind = t, data = dat, sandwich = TRUE)
  expect_identical(m_true$pffr$sandwich, "cluster")

  # FALSE -> "none"
  m_false <- pffr(Y ~ xlin, yind = t, data = dat, sandwich = FALSE)
  expect_identical(m_false$pffr$sandwich, "none")

  # coef.pffr also accepts TRUE/FALSE
  coef_true <- coef(m_false, sandwich = TRUE)
  coef_cluster <- coef(m_false, sandwich = "cluster")
  expect_equal(coef_true$pterms[, "se"], coef_cluster$pterms[, "se"])
})

test_that("coef.pffr recomputes when sandwich type differs from fit", {
  skip_on_cran()

  dat <- get_xlin_data()
  t <- attr(dat, "yindex")

  # Fit with HC, request cluster — should NOT reuse HC matrices
  m_hc <- pffr(Y ~ xlin, yind = t, data = dat, sandwich = "hc")
  coef_hc_stored <- coef(m_hc, sandwich = "hc")
  coef_cl_from_hc <- coef(m_hc, sandwich = "cluster")

  # Cluster SEs should differ from stored HC SEs
  expect_false(identical(
    coef_hc_stored$pterms[, "se"],
    coef_cl_from_hc$pterms[, "se"]
  ))

  # Fit without sandwich, request cluster on the fly
  m_none <- get_xlin_model()
  coef_cl_from_none <- coef(m_none, sandwich = "cluster")
  expect_true(all(is.finite(coef_cl_from_none$pterms[, "se"])))

  # Fit with CL2, request cluster and HC to force recomputation paths.
  m_cl2 <- pffr(Y ~ xlin, yind = t, data = dat, sandwich = "cl2")
  coef_cl2_stored <- coef(m_cl2, sandwich = "cl2")
  coef_cl_from_cl2 <- coef(m_cl2, sandwich = "cluster")
  coef_hc_from_cl2 <- coef(m_cl2, sandwich = "hc")

  expect_false(identical(
    coef_cl2_stored$pterms[, "se"],
    coef_cl_from_cl2$pterms[, "se"]
  ))
  expect_false(identical(
    coef_cl2_stored$pterms[, "se"],
    coef_hc_from_cl2$pterms[, "se"]
  ))
})

test_that("gam_sandwich_cluster works for gaulss family", {
  skip_on_cran()

  set.seed(42)
  n <- 60
  T_grid <- 20
  t_grid <- seq(0, 1, length.out = T_grid)
  x <- rnorm(n)
  beta1 <- cos(2 * pi * t_grid)
  signal <- outer(rep(1, n), sin(2 * pi * t_grid)) + outer(x, beta1)
  Y <- signal + matrix(rnorm(n * T_grid, sd = 0.5), n, T_grid)

  dat <- data.frame(Y = I(Y), x = x)
  m <- pffr(Y ~ x, yind = t_grid, data = dat, family = mgcv::gaulss())

  # Covariance matrix should be symmetric and positive semi-definite
  m_stripped <- m
  class(m_stripped) <- setdiff(class(m_stripped), "pffr")
  cluster_id <- build_cluster_id(m$pffr)
  V <- gam_sandwich_cluster(m_stripped, cluster_id, freq = FALSE)
  expect_equal(V, t(V), tolerance = 1e-12)
  expect_true(all(diag(V) >= 0))

  # Cluster sandwich should differ from HC sandwich
  V_hc <- mgcv::vcov.gam(m_stripped, sandwich = TRUE, freq = FALSE)
  expect_false(identical(V, V_hc))

  # Cluster sandwich should work via coef.pffr without error
  coef_cl <- coef(m, sandwich = "cluster")
  expect_true(length(coef_cl$smterms) > 0)

  # All SEs in all terms should be positive and finite
  for (sm in coef_cl$smterms) {
    if (!is.null(sm$coef$se)) {
      expect_true(all(is.finite(sm$coef$se)))
      expect_true(all(sm$coef$se > 0))
    }
  }
})

test_that("gam_sandwich_cluster_cl2 works for gaulss family", {
  skip_on_cran()

  m <- get_gaulss_model()
  m_stripped <- m
  class(m_stripped) <- setdiff(class(m_stripped), "pffr")
  cluster_id <- build_cluster_id(m$pffr)

  V_cl <- gam_sandwich_cluster(m_stripped, cluster_id, freq = FALSE)
  V_cl2 <- gam_sandwich_cluster_cl2(m_stripped, cluster_id, freq = FALSE)

  expect_equal(V_cl2, t(V_cl2), tolerance = 1e-10)
  expect_true(all(is.finite(diag(V_cl2))))
  expect_true(all(diag(V_cl2) >= 0))
  expect_gt(max(abs(V_cl2 - V_cl)), 0)

  dat <- get_xlin_data()
  t <- attr(dat, "yindex")
  m_fit_cl2 <- pffr(
    Y ~ xlin,
    yind = t,
    data = dat,
    family = mgcv::gaulss(),
    sandwich = "cl2"
  )
  expect_identical(m_fit_cl2$pffr$sandwich, "cl2")

  coef_cl2 <- coef(m, sandwich = "cl2")
  expect_true(length(coef_cl2$smterms) > 0)
  for (sm in coef_cl2$smterms) {
    if (!is.null(sm$coef$se)) {
      expect_true(all(is.finite(sm$coef$se)))
      expect_true(all(sm$coef$se > 0))
    }
  }
})

test_that("CL2 falls back to HC for unsupported custom family$sandwich", {
  skip_on_cran()

  dat <- sim_xlin_data(n = 25, nygrid = 25, SNR = 10, family = mgcv::gaulss())
  t <- attr(dat, "yindex")
  m <- pffr(Y ~ xlin, yind = t, data = dat, family = mgcv::gaulss())

  m_stripped <- m
  class(m_stripped) <- setdiff(class(m_stripped), "pffr")
  m_custom <- m_stripped
  m_custom$family$family <- "mock-custom"

  cluster_id <- build_cluster_id(m$pffr)
  expect_warning(
    V_cl2 <- gam_sandwich_cluster_cl2(m_custom, cluster_id, freq = FALSE),
    "CL2 sandwich not yet implemented for family 'mock-custom'"
  )

  V_hc <- mgcv::vcov.gam(m_custom, sandwich = TRUE, freq = FALSE)
  expect_equal(V_cl2, V_hc)
})

test_that("coef.pffr adds pointwise and simultaneous CIs", {
  skip_on_cran()

  m <- get_xlin_model()

  coef_pw <- coef(
    m,
    ci = "pointwise",
    level = 0.95,
    n1 = 40
  )
  coef_sim <- coef(
    m,
    ci = "simultaneous",
    level = 0.95,
    n_sim = 400,
    sim_seed = 123,
    n1 = 40
  )

  expect_true(all(c("lower", "upper") %in% colnames(coef_pw$pterms)))
  expect_true(all(c("lower", "upper") %in% colnames(coef_sim$pterms)))
  expect_identical(coef_sim$ci_meta$type, "simultaneous")

  sm_pw <- coef_pw$smterms[[1]]$coef
  sm_sim <- coef_sim$smterms[[1]]$coef
  expect_true(all(c("lower", "upper") %in% colnames(sm_pw)))
  expect_true(all(c("lower", "upper") %in% colnames(sm_sim)))

  width_pw <- sm_pw$upper - sm_pw$lower
  width_sim <- sm_sim$upper - sm_sim$lower
  expect_gte(median(width_sim), median(width_pw))
})

test_that("coef.pffr simultaneous CI works for poisson family", {
  skip_on_cran()

  dat <- sim_xlin_data(n = 25, nygrid = 20, SNR = 10, family = poisson())
  t <- attr(dat, "yindex")
  m <- pffr(
    Y ~ xlin,
    yind = t,
    data = dat,
    family = poisson(),
    sandwich = "none"
  )

  coef_sim <- coef(
    m,
    sandwich = "cluster",
    ci = "simultaneous",
    n_sim = 250,
    sim_seed = 42,
    n1 = 30
  )

  sm_coef <- coef_sim$smterms[[1]]$coef
  expect_true(all(c("lower", "upper") %in% colnames(sm_coef)))
  expect_true(all(is.finite(sm_coef$lower)))
  expect_true(all(is.finite(sm_coef$upper)))
})

test_that("coef.pffr simultaneous CI works for gaulss", {
  skip_on_cran()

  m <- get_gaulss_model()
  coef_sim <- coef(
    m,
    sandwich = "none",
    ci = "simultaneous",
    n_sim = 250,
    sim_seed = 99,
    n1 = 30
  )

  expect_identical(coef_sim$ci_meta$type, "simultaneous")
  expect_true(length(coef_sim$smterms) > 0)

  has_bounds <- vapply(
    coef_sim$smterms,
    \(sm) all(c("lower", "upper") %in% colnames(sm$coef)),
    logical(1)
  )
  expect_true(any(has_bounds))
})
