#--------------------------------------
# Numerical integration weights for ff()/sff()
#--------------------------------------
#
# Regression tests for the Simpson weights used to build the integration
# operator L in ff()/sff(). Up to refund 0.1-40 the [1, 4, 2, ..., 4, 1]
# pattern was scaled by (b - a) / (3 * n) instead of (b - a) / (3 * (n - 1)),
# and for even n the alternation ended in 2 before the closing 1. The weights
# therefore summed to less than b - a (a constant on [0, 1] integrated to
# 0.956 at n = 30 and to 0.989 at n = 93), which rescaled estimated ff()
# surfaces by the reciprocal of that factor.

grid_sizes <- c(2:12, 30, 31, 60, 61, 93)

# f(x) = 1 + 2x - 3x^2 + 4x^3, with int_0^1 f = 1 + 1 - 1 + 1 = 2
cubic <- function(x) 1 + 2 * x - 3 * x^2 + 4 * x^3

test_that("simpson weights integrate a constant exactly on [0, 1]", {
  for (n in grid_sizes) {
    w <- refund:::simpson_weights(seq(0, 1, length.out = n))
    expect_equal(sum(w), 1, tolerance = 1e-12, info = paste("n =", n))
  }
})

test_that("simpson weights integrate a constant exactly on a general interval", {
  for (n in grid_sizes) {
    xind <- seq(-2, 5, length.out = n)
    w <- refund:::simpson_weights(xind)
    expect_equal(sum(w), 7, tolerance = 1e-10, info = paste("n =", n))
  }
})

test_that("simpson weights integrate a cubic exactly for odd n", {
  for (n in grid_sizes[grid_sizes %% 2 == 1]) {
    xind <- seq(0, 1, length.out = n)
    w <- refund:::simpson_weights(xind)
    expect_equal(
      sum(w * cubic(xind)),
      2,
      tolerance = 1e-12,
      info = paste("n =", n)
    )
  }
})

test_that("the even-n closing rule (Simpson 3/8) is also exact for cubics", {
  # even n uses composite Simpson on the first n - 3 points plus Simpson's
  # 3/8 rule on the last three intervals, which stays exact for cubics
  for (n in grid_sizes[grid_sizes %% 2 == 0 & grid_sizes >= 4]) {
    xind <- seq(0, 1, length.out = n)
    w <- refund:::simpson_weights(xind)
    expect_equal(
      sum(w * cubic(xind)),
      2,
      tolerance = 1e-12,
      info = paste("n =", n)
    )
  }
  # n = 2 falls back to the trapezoidal rule: exact for linear, not cubic
  w2 <- refund:::simpson_weights(c(0, 1))
  expect_equal(unname(w2), c(0.5, 0.5))
})

test_that("weight patterns match the textbook composite rules", {
  expect_equal(refund:::simpson_pattern(3), c(1, 4, 1) / 3)
  expect_equal(refund:::simpson_pattern(5), c(1, 4, 2, 4, 1) / 3)
  expect_equal(refund:::simpson_pattern(4), c(1, 3, 3, 1) * 3 / 8)
  expect_equal(
    refund:::simpson_pattern(6),
    c(1, 4, 1, 0, 0, 0) / 3 + c(0, 0, 1, 3, 3, 1) * 3 / 8
  )
  # pattern sums to n - 1 so that h * pattern sums to b - a
  for (n in grid_sizes) {
    expect_equal(
      sum(refund:::simpson_pattern(n)),
      n - 1,
      tolerance = 1e-12,
      info = paste("n =", n)
    )
  }
})

test_that("compute_integration_weights() handles matrix input row-wise", {
  xind <- rbind(
    seq(0, 1, length.out = 31),
    seq(2, 5, length.out = 31)
  )
  L <- refund:::compute_integration_weights(xind, "simpson")
  expect_equal(dim(L), c(2L, 31L))
  expect_equal(rowSums(L), c(1, 3), tolerance = 1e-12)
  expect_equal(L[1, ], refund:::simpson_weights(xind[1, ]))
  expect_equal(L[2, ], refund:::simpson_weights(xind[2, ]))
})

test_that("trapezoidal weights also sum to the length of the domain", {
  equi <- rbind(seq(0, 1, length.out = 30), seq(0, 1, length.out = 30))
  expect_equal(
    rowSums(refund:::compute_integration_weights(equi, "trapezoidal")),
    c(1, 1),
    tolerance = 1e-12
  )
  noneq <- rbind(sort(c(0, 1, runif(28))), sort(c(0, 1, runif(28))))
  expect_equal(
    rowSums(refund:::compute_integration_weights(noneq, "trapezoidal")),
    c(1, 1),
    tolerance = 1e-12
  )
})

test_that("simpson_legacy reproduces the pre-0.1-41 weights", {
  legacy_sum <- function(n) {
    sum(refund:::simpson_weights(seq(0, 1, length.out = n), "simpson_legacy"))
  }
  # documented effect sizes for a constant on [0, 1]
  expect_equal(legacy_sum(30), 0.9556, tolerance = 1e-4)
  expect_equal(legacy_sum(60), 0.9778, tolerance = 1e-4)
  expect_equal(legacy_sum(31), 0.9677, tolerance = 1e-4)
  expect_equal(legacy_sum(61), 0.9836, tolerance = 1e-4)
  expect_equal(legacy_sum(93), 0.9892, tolerance = 1e-4)

  # verbatim old formula: ((b - a) / n) / 3 * [1, 4, 2, ..., 4/2, 1]
  for (n in grid_sizes) {
    xind <- seq(0, 1, length.out = n)
    old <- ((xind[n] - xind[1]) / n / 3) *
      c(1, rep(c(4, 2), length.out = n - 2), 1)
    expect_equal(
      refund:::simpson_weights(xind, "simpson_legacy"),
      old,
      info = paste("n =", n)
    )
  }
  # for odd n the fix is a pure rescaling by n / (n - 1)
  for (n in c(31, 61, 93)) {
    xind <- seq(0, 1, length.out = n)
    expect_equal(
      refund:::simpson_weights(xind),
      refund:::simpson_weights(xind, "simpson_legacy") * n / (n - 1)
    )
  }
})

test_that("compute_integration_weights() rejects unknown methods", {
  expect_error(
    refund:::compute_integration_weights(
      matrix(seq(0, 1, length.out = 5), nrow = 1),
      "midpoint"
    ),
    "Unknown `integration` method"
  )
})

test_that("ff() exposes integration = 'simpson_legacy'", {
  set.seed(20260918)
  xind <- seq(0, 1, length.out = 31)
  X <- matrix(rnorm(20 * 31), 20, 31)

  trm_new <- ff(X, xind = xind, check.ident = FALSE)
  trm_old <- ff(
    X,
    xind = xind,
    integration = "simpson_legacy",
    check.ident = FALSE
  )

  expect_equal(rowSums(trm_new$L), rep(1, 20), tolerance = 1e-12)
  expect_equal(trm_old$L[1, ], refund:::simpson_weights(xind, "simpson_legacy"))
  expect_equal(trm_old$L / trm_new$L, matrix(30 / 31, 20, 31))
})

test_that("the ff() surface is rescaled by the old/new weight-sum ratio", {
  skip_on_cran()
  set.seed(20260918)
  nxgrid <- 31
  dat <- pffr_simulate(
    Y ~ ff(X1, xind = s),
    n = 40,
    nxgrid = nxgrid,
    nygrid = 20,
    SNR = 20,
    effects = list(X1 = "cosine")
  )
  yind <- attr(dat, "yindex")
  xind <- attr(dat, "xindex")

  fit_new <- pffr(
    Y ~ ff(X1, xind = xind),
    yind = yind,
    data = dat,
    sandwich = "none"
  )
  fit_old <- pffr(
    Y ~ ff(X1, xind = xind, integration = "simpson_legacy"),
    yind = yind,
    data = dat,
    sandwich = "none"
  )

  # identical fits, only the parameterisation of beta(s, t) changes
  expect_equal(fitted(fit_new), fitted(fit_old), tolerance = 1e-8)

  beta_new <- coef(fit_new, n1 = 12, n2 = 12, n3 = 12)$smterms[["ff(X1)"]]$value
  beta_old <- coef(fit_old, n1 = 12, n2 = 12, n3 = 12)$smterms[["ff(X1)"]]$value

  # sum(w_legacy) / sum(w_simpson) = (n - 1) / n for odd n
  ratio <- (nxgrid - 1) / nxgrid
  expect_gt(max(abs(beta_old)), 0.1)
  expect_equal(beta_new, beta_old * ratio, tolerance = 1e-5)
})
