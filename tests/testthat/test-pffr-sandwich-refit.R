# Regression test: recomputing a sandwich covariance on a fit that was itself
# created with a sandwich option (e.g. the default sandwich = "cluster") must
# use the model-based Vp/Ve as the bread, not the already-robustified matrices.
# Before the model_cov stash + restore_model_cov(), coef(fit, sandwich = "cl2")
# on a default-fitted model double-applied the correction (SEs inflated by
# ~1.5-2x on this example); sandwich = "none" returned robust instead of
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

  expect_false(is.null(fit_cluster$pffr$model_cov))
  # stash holds the model-based matrices, not the robust ones
  expect_equal(fit_cluster$pffr$model_cov$Vp, fit_none$Vp, tolerance = 1e-8)
  expect_equal(fit_cluster$pffr$model_cov$Ve, fit_none$Ve, tolerance = 1e-8)

  # small evaluation grids: only SE equality matters, not grid resolution
  ff_se <- function(fit, ...) {
    coef(fit, n1 = 20, n2 = 8, ...)$smterms[[2]]$coef$se
  }

  for (type in c("cl2", "hc", "none")) {
    se_ref <- ff_se(fit_none, sandwich = type)
    se_re <- ff_se(fit_cluster, sandwich = type)
    expect_equal(se_re, se_ref, tolerance = 1e-8, label = paste0("se_", type))
  }

  # frequentist bread must be restored too
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

  # vcov(): default returns the stored (robust) matrix; sandwich = TRUE
  # recomputes the HC sandwich from the restored model-based bread
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

  # re-applying apply_sandwich_correction is idempotent
  twice <- refund:::apply_sandwich_correction(
    fit_cluster,
    "gam",
    type = "cluster"
  )
  expect_equal(twice$Vp, fit_cluster$Vp, tolerance = 1e-10)
  expect_equal(twice$Ve, fit_cluster$Ve, tolerance = 1e-10)
  expect_equal(twice$pffr$model_cov, fit_cluster$pffr$model_cov)

  # pre-stash objects (older refund versions) warn instead of silently
  # double-applying / mislabeling robust matrices as model-based
  fit_old <- fit_cluster
  fit_old$pffr$model_cov <- NULL
  w_coef <- capture_warnings(ff_se(fit_old, sandwich = "cl2"))
  expect_true(any(grepl("ON TOP", w_coef)))
  w_apply <- capture_warnings(
    refund:::apply_sandwich_correction(fit_old, "gam", type = "cluster")
  )
  expect_true(any(grepl("ON TOP", w_apply)))
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
