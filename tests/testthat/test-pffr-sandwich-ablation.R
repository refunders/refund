# X5/X6 internal sandwich ablation switches: `b2 = FALSE` drops the additive
# Bayesian smoothing-bias term B2 = Vp - Ve; `center_scores = TRUE` centers the
# per-cluster score sums (U_g - (sum_g U_g)/G) before forming the meat. Both
# default to the shipped behavior and are reachable through pffr_vcov() and
# coef(..., b2 =, center_scores =) via `...`. These tests assert the exact
# algebraic identities and that the defaults are byte-unchanged.

make_ablation_fit <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      dat <- pffr_simulate(
        Y ~ ff(X1),
        n = 25,
        nxgrid = 20,
        nygrid = 20,
        SNR = 5,
        effects = list(X1 = "random"),
        intercept = "random",
        seed = 5150
      )
      yind <- attr(dat, "yindex")
      cache <<- list(
        fit_none = pffr(Y ~ ff(X1), data = dat, yind = yind, sandwich = "none"),
        fit_cluster = pffr(
          Y ~ ff(X1),
          data = dat,
          yind = yind,
          sandwich = "cluster"
        )
      )
    }
    cache
  }
})

test_that("(i) b2 = FALSE drops exactly the additive B2 term (cluster & cl2)", {
  skip_on_cran()
  fit <- make_ablation_fit()$fit_none
  B2 <- fit$Vp - fit$Ve

  for (sw in c("cluster", "cl2")) {
    V_def <- refund:::pffr_vcov(fit, sandwich = sw)
    V_nob2 <- refund:::pffr_vcov(fit, sandwich = sw, b2 = FALSE)
    # V_default - V_nob2 == B2 = Vp - Ve, exactly (compare values; CL2 carries
    # extra leverage attributes)
    expect_equal(as.numeric(V_def - V_nob2), as.numeric(B2), tolerance = 1e-10)
    # SEs <= default everywhere (B2 adds a PSD term to the diagonal)
    expect_true(all(diag(V_nob2) <= diag(V_def) + 1e-10))
  }

  # same through the public coef() path (seWithMean = FALSE => term SE^2 =
  # diag(Xp V[trm] Xp')): the nob2 SEs never exceed the default SEs
  co_def <- coef(fit, sandwich = "cluster", seWithMean = FALSE, n1 = 20, n2 = 8)
  co_nob2 <- coef(
    fit,
    sandwich = "cluster",
    b2 = FALSE,
    seWithMean = FALSE,
    n1 = 20,
    n2 = 8
  )
  se_def <- co_def$smterms[[2]]$coef$se
  se_nob2 <- co_nob2$smterms[[2]]$coef$se
  expect_true(all(se_nob2 <= se_def + 1e-8))
  # and the SE^2 gap equals the B2 quadratic form (must be strictly positive
  # somewhere -- B2 is load-bearing)
  expect_gt(max(se_def^2 - se_nob2^2), 0)
})

test_that("(ii) center_scores = TRUE subtracts exactly the rank-one penalty term", {
  skip_on_cran()
  fit <- make_ablation_fit()$fit_none
  sh <- refund:::pffr_sandwich_shares(fit)

  # CR1 (cluster) meat: centering removes hc1 * Vp (Sbeta Sbeta'/G) Vp exactly,
  # where Sbeta = sum_g U_g is the total score. B2 cancels in the difference.
  V_unc <- refund:::pffr_vcov(fit, sandwich = "cluster")
  V_cen <- refund:::pffr_vcov(fit, sandwich = "cluster", center_scores = TRUE)
  expected <- sh$hc1 * sh$Vp %*% (tcrossprod(sh$Sbeta) / sh$G) %*% sh$Vp
  expect_equal(
    as.numeric(V_unc - V_cen),
    as.numeric(expected),
    tolerance = 1e-10
  )
  # centering can only shrink the penalty direction: overall Frobenius norm drops
  expect_lt(norm(V_cen, "F"), norm(V_unc, "F"))

  # CL2 path: centering changes the meat and stays symmetric (its U_g are the
  # leverage-adjusted contributions, so the shift is rank-one in those, not in
  # Sbeta -- checked structurally here)
  V2_unc <- refund:::pffr_vcov(fit, sandwich = "cl2")
  V2_cen <- refund:::pffr_vcov(fit, sandwich = "cl2", center_scores = TRUE)
  expect_gt(max(abs(V2_unc - V2_cen)), 0)
  expect_equal(V2_cen, t(V2_cen), tolerance = 1e-10)
})

test_that("(iii) defaults are unchanged (b2 = TRUE, center_scores = FALSE)", {
  skip_on_cran()
  fits <- make_ablation_fit()
  fit <- fits$fit_none

  for (sw in c("cluster", "cl2")) {
    V_flagged <- refund:::pffr_vcov(
      fit,
      sandwich = sw,
      b2 = TRUE,
      center_scores = FALSE
    )
    V_plain <- refund:::pffr_vcov(fit, sandwich = sw)
    expect_equal(V_flagged, V_plain, tolerance = 1e-12)
  }

  # the recomputed cluster default equals the fit-time stored robust matrix
  V_plain_cluster <- refund:::pffr_vcov(fit, sandwich = "cluster")
  expect_equal(
    fits$fit_cluster$pffr$Vsandwich,
    V_plain_cluster,
    tolerance = 1e-12
  )

  # coef() with the default flags is identical to coef() without them
  co_a <- coef(fit, sandwich = "cluster", n1 = 20, n2 = 8)
  co_b <- coef(
    fit,
    sandwich = "cluster",
    b2 = TRUE,
    center_scores = FALSE,
    n1 = 20,
    n2 = 8
  )
  expect_equal(
    co_a$smterms[[2]]$coef$se,
    co_b$smterms[[2]]$coef$se,
    tolerance = 1e-12
  )
})

test_that("pffr_sandwich_shares returns consistent pen_share and fro_ratio", {
  skip_on_cran()
  fit <- make_ablation_fit()$fit_none
  sh <- refund:::pffr_sandwich_shares(fit)

  expect_true(is.finite(sh$pen_share) && sh$pen_share >= 0)
  expect_true(is.finite(sh$fro_ratio) && sh$fro_ratio >= 0)
  expect_identical(sh$G, 25L)

  # fro_ratio ties to the accessor: ||B2||_F / ||hc1 Vp meat Vp||_F, where the
  # cluster core is exactly the b2 = FALSE covariance and B2 = default - nob2
  V_def <- refund:::pffr_vcov(fit, sandwich = "cluster")
  V_nob2 <- refund:::pffr_vcov(fit, sandwich = "cluster", b2 = FALSE)
  fro_from_cov <- norm(V_def - V_nob2, "F") / norm(V_nob2, "F")
  expect_equal(sh$fro_ratio, fro_from_cov, tolerance = 1e-8)

  # pen_share = ||Sbeta||^2 / sum_g ||U_g||^2; sum_g ||U_g||^2 = tr(meat)
  # (trace of crossprod(U) = sum of all squared entries of U)
  expect_equal(
    sh$pen_share,
    sum(sh$Sbeta^2) / sum(diag(sh$meat)),
    tolerance = 1e-10
  )
})
