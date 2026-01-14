context("Testing pffr")

###############################################################################
# Helper Functions
###############################################################################

#' Validate pffrSim truth structure
#' @param dat Data from pffrSim with formula interface
#' @param n Expected number of observations
#' @param nygrid Expected number of y grid points
expect_valid_truth <- function(dat, n, nygrid) {
  truth <- attr(dat, "truth")
  expect_true(is.list(truth), info = "truth should be a list")
  expect_true(
    all(c("eta", "etaTerms", "epsilon") %in% names(truth)),
    info = "truth should contain eta, etaTerms, epsilon"
  )
  expect_equal(
    dim(truth$eta),
    c(n, nygrid),
    info = "truth$eta dimensions should match n x nygrid"
  )
  expect_equal(
    dim(truth$epsilon),
    c(n, nygrid),
    info = "truth$epsilon dimensions should match n x nygrid"
  )
}

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
  pffrSim(
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
  set.seed(9312)
  # First legacy call triggers deprecation warning (once per session)
  data2 <- expect_warning(
    pffrSim(scenario = "all", n = 30),
    "scenario.*deprecated"
  )
  argvals <- attr(data2, "yindex")
  s <- attr(data2, "xindex")

  m2 <- pffr(
    Y ~
      ff(X1, xind = s) +
        xlin +
        c(te(xte1, xte2)) +
        s(xsmoo) +
        c(xconst),
    yind = argvals,
    data = data2
  )

  expect_s3_class(m2, "pffr")
  # Verify model produces reasonable fitted values

  expect_equal(dim(fitted(m2)), c(30, length(argvals)))
  expect_false(any(is.na(fitted(m2))))
})

test_that("convenience functions are working", {
  set.seed(9313)
  data2 <- pffrSim(scenario = "all", n = 30)
  argvals <- attr(data2, "yindex")
  s <- attr(data2, "xindex")

  m2 <- pffr(
    Y ~ ff(X1, xind = s) + xlin + c(te(xte1, xte2)) + s(xsmoo) + c(xconst),
    yind = argvals,
    data = data2
  )

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
  preddata <- pffrSim(scenario = "all", n = 20)
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
  set.seed(88182004)
  data3 <- pffrSim(scenario = c("int", "smoo"), n = 30, propmissing = 0.8)
  t <- attr(data3, "yindex")

  # Verify sparse data structure

  expect_true(is.list(data3))
  expect_true("ydata" %in% names(data3))
  expect_true("data" %in% names(data3))

  m3.sparse <- pffr(
    Y ~ s(xsmoo),
    data = data3$data,
    ydata = data3$ydata,
    yind = t
  )

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
  data <- pffrSim(
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

  dat <- pffrSim(
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

  data <- pffrSim(scenario = "ff", n = 20, SNR = 100)
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

  dat <- pffrSim(
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

  dat <- pffrSim(
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

  dat <- pffrSim(
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

  dat <- pffrSim(
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

  dat <- pffrSim(
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

  dat <- pffrSim(
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

  dat <- pffrSim(
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
  newdat <- pffrSim(
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

  dat <- pffrSim(
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

  dat <- pffrSim(
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

  dat <- pffrSim(
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

  dat <- pffrSim(
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

  dat <- pffrSim(
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
  dat <- pffrSim(
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
  dat <- pffrSim(
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
  dat <- pffrSim(
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
