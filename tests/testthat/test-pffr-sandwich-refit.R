# Regression test: recomputing a sandwich covariance on a fit that was itself
# created with a sandwich option (e.g. the default sandwich = "cluster") must
# use the model-based Vp/Ve as the bread, not an already-robustified matrix.
# Under the current storage contract $Vp/$Vc/$Ve are model-based ALWAYS (the
# robust covariance lives in $pffr$Vsandwich), so this holds by construction;
# these tests lock the public surface. Historically, coef(fit, sandwich="cl2")
# on a default-fitted model double-applied the correction (SEs inflated by
# ~1.5-2x on this example) and sandwich = "none" returned robust instead of
# model-based SEs.

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

  fit_none <- pffr(Y ~ ff(X1), data = dat, yind = yind, sandwich = "none")
  fit_cluster <- pffr(Y ~ ff(X1), data = dat, yind = yind, sandwich = "cluster")

  # $Vp/$Vc/$Ve are the model-based matrices on BOTH fits; the robust
  # covariance is stored separately with its metadata.
  expect_identical(fit_cluster$Vp, fit_none$Vp)
  expect_identical(fit_cluster$Ve, fit_none$Ve)
  expect_false(is.null(fit_cluster$pffr$Vsandwich))
  expect_false(is.null(fit_cluster$pffr$Vsandwich_freq))
  expect_identical(fit_cluster$pffr$sandwich_info$type, "cluster")
  expect_null(fit_cluster$pffr$model_cov)

  # small evaluation grids: only SE equality matters, not grid resolution
  ff_se <- function(fit, ...) {
    coef(fit, n1 = 20, n2 = 8, ...)$smterms[[2]]$coef$se
  }

  for (type in c("cl2", "hc", "none")) {
    se_ref <- ff_se(fit_none, sandwich = type)
    se_re <- ff_se(fit_cluster, sandwich = type)
    expect_equal(se_re, se_ref, tolerance = 1e-8, label = paste0("se_", type))
  }

  # frequentist bread must be model-based too
  se_ref <- ff_se(fit_none, sandwich = "cl2", freq = TRUE)
  se_re <- ff_se(fit_cluster, sandwich = "cl2", freq = TRUE)
  expect_equal(se_re, se_ref, tolerance = 1e-8)

  # explicit cluster= forces recomputation even for the fitted type
  se_ref <- ff_se(fit_none, sandwich = "cluster")
  se_re <- ff_se(fit_cluster, sandwich = "cluster", cluster = seq_len(25))
  expect_equal(se_re, se_ref, tolerance = 1e-8)

  # a single cluster is rejected instead of dividing by G - 1 = 0
  expect_error(
    ff_se(fit_none, sandwich = "cluster", cluster = rep(1, 25)),
    "at least two clusters"
  )

  # vcov(): returns the model-based covariance on every fit (mgcv's default is
  # Vp); sandwich = TRUE recomputes the HC sandwich from the model-based bread
  expect_equal(
    unname(vcov(fit_cluster)),
    unname(vcov(fit_none)),
    tolerance = 1e-10
  )
  expect_equal(
    unname(vcov(fit_cluster)),
    unname(fit_cluster$Vp),
    tolerance = 1e-10
  )
  expect_equal(
    vcov(fit_cluster, sandwich = TRUE),
    vcov(fit_none, sandwich = TRUE),
    tolerance = 1e-8
  )

  # re-applying apply_sandwich_correction is idempotent: $Vp untouched,
  # identical robust matrices
  twice <- refund:::apply_sandwich_correction(
    fit_cluster,
    "gam",
    type = "cluster"
  )
  expect_identical(twice$Vp, fit_cluster$Vp)
  expect_identical(twice$Ve, fit_cluster$Ve)
  expect_equal(
    twice$pffr$Vsandwich,
    fit_cluster$pffr$Vsandwich,
    tolerance = 1e-10
  )
})

test_that("CR1 scores include prior weights (matches sandwich::vcovCL)", {
  skip_on_cran()
  skip_if_not_installed("sandwich")

  set.seed(42)
  n <- 200
  G <- 40
  cl <- rep(seq_len(G), each = n / G)
  x <- runif(n)
  w <- runif(n, 0.5, 2)
  y <- 1 + 2 * x + rnorm(n, sd = 0.5)

  b <- mgcv::gam(y ~ x, weights = w)
  fit_glm <- glm(y ~ x, weights = w)

  V_refund <- refund:::gam_sandwich_cluster(b, cl, freq = TRUE)
  V_ref <- sandwich::vcovCL(fit_glm, cluster = cl, type = "HC0", cadjust = TRUE)
  expect_equal(unname(V_refund), unname(V_ref), tolerance = 1e-6)
})
