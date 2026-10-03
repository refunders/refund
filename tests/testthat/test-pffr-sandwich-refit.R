# Regression test: recomputing the sandwich covariance on a fit that was itself
# created with the sandwich must use the model-based Vp/Ve as the bread, not an
# already-robustified matrix. Under the current storage contract $Vp/$Vc/$Ve
# are model-based ALWAYS (the robust covariance lives in $pffr$Vsandwich), so
# this holds by construction; these tests lock the public surface.

test_that("sandwich recomputation on a corrected fit uses model-based bread", {
  skip_on_cran()

  dat <- pffr_simulate(
    Y ~ ff(X1),
    n = 25,
    nxgrid = 20,
    nygrid = 20,
    SNR = 5,
    effects = list(X1 = "random"),
    intercept = "random",
    seed = 1234
  )
  yind <- attr(dat, "yindex")

  fit_none <- pffr(Y ~ ff(X1), data = dat, yind = yind, sandwich = FALSE)
  fit_cl2 <- quiet_pffr(Y ~ ff(X1), data = dat, yind = yind)

  # $Vp/$Vc/$Ve are the model-based matrices on BOTH fits; the robust
  # covariance is stored separately with its metadata.
  expect_identical(fit_cl2$Vp, fit_none$Vp)
  expect_identical(fit_cl2$Ve, fit_none$Ve)
  expect_false(is.null(fit_cl2$pffr[["Vsandwich"]]))
  expect_identical(fit_cl2$pffr$sandwich_info$type, "cl2")
  expect_null(fit_cl2$pffr$model_cov)

  # small evaluation grids: only SE equality matters, not grid resolution
  ff_se <- function(fit, ...) {
    coef(fit, n1 = 20, n2 = 8, ...)$smterms[[2]]$coef$se
  }
  for (sw in c(TRUE, FALSE)) {
    expect_equal(
      ff_se(fit_cl2, sandwich = sw),
      ff_se(fit_none, sandwich = sw),
      tolerance = 1e-8
    )
  }

  # explicit cluster= forces recomputation and reproduces the stored matrix
  expect_equal(
    ff_se(fit_cl2, cluster = seq_len(25)),
    ff_se(fit_cl2),
    tolerance = 1e-8
  )

  # a single cluster is rejected instead of dividing by G - 1 = 0
  expect_error(
    ff_se(fit_none, sandwich = TRUE, cluster = rep(1, 25)),
    "at least two clusters"
  )

  # vcov(): the fit's covariance; sandwich = TRUE/FALSE switches
  expect_equal(vcov(fit_cl2), fit_cl2$pffr$Vsandwich)
  expect_identical(vcov(fit_none), fit_none$Vc)
  expect_identical(vcov(fit_cl2, sandwich = FALSE), fit_cl2$Vc)
  expect_equal(
    unname(matrix(vcov(fit_none, sandwich = TRUE), ncol(fit_none$Vp))),
    unname(matrix(vcov(fit_cl2), ncol(fit_none$Vp))),
    tolerance = 1e-10
  )

  # re-applying apply_sandwich_correction is idempotent: $Vp untouched,
  # identical robust matrices
  twice <- withCallingHandlers(
    refund:::apply_sandwich_correction(fit_cl2),
    pffr_small_G_warning = function(w) invokeRestart("muffleWarning")
  )
  expect_identical(twice$Vp, fit_cl2$Vp)
  expect_identical(twice$Ve, fit_cl2$Ve)
  expect_equal(twice$pffr$Vsandwich, fit_cl2$pffr$Vsandwich, tolerance = 1e-10)
})

test_that("CL2 scores include prior weights (dense reference, unpenalized)", {
  set.seed(42)
  n <- 120
  G <- 30
  cl <- rep(seq_len(G), each = n / G)
  x <- runif(n)
  w <- runif(n, 0.5, 2)
  y <- 1 + 2 * x + rnorm(n, sd = 0.5)
  b <- mgcv::gam(y ~ x, weights = w)

  V <- refund:::gam_sandwich_cluster_cl2(b, cl)
  X <- model.matrix(b)
  s <- sqrt(w / b$sig2)
  Xw <- X * s
  z <- (y - fitted(b)) * s
  M <- diag(n) - Xw %*% b$Vp %*% t(Xw)
  meat <- matrix(0, 2, 2)
  for (g in seq_len(G)) {
    ii <- which(cl == g)
    Bg <- tcrossprod(M[ii, , drop = FALSE])
    ee <- eigen(Bg, symmetric = TRUE)
    A <- ee$vectors %*% diag(1 / sqrt(ee$values)) %*% t(ee$vectors)
    u <- crossprod(Xw[ii, , drop = FALSE], A %*% z[ii])
    meat <- meat + tcrossprod(u)
  }
  V_ref <- G / (G - 1) * b$Vp %*% meat %*% b$Vp + (b$Vp - b$Ve)
  expect_equal(unname(matrix(V, 2)), unname(V_ref), tolerance = 1e-8)
})
