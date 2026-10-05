# A18 exact-CL2 hardening gate. The standalone script writes the same table
# for release records; this fast test keeps its family/design coverage in CI.
# Fixtures live in helper-exactcl2.R.

test_that("exact CL2 is robust on high-leverage family and basis fixtures", {
  results <- run_exactcl2_hardening()

  expect_equal(nrow(results), 16L)
  expect_true(all(is.finite(results$n_adjusted)))
  expect_gt(sum(results$n_adjusted), 0)
  expect_true(all(is.finite(results$min_block_eig)))
  expect_true(all(is.finite(results$max_block_kappa)))
  expect_true(all(is.finite(results$max_abs_exact_nocap)))
  expect_true(all(results$se_finite_positive))
  expect_true(all(results$se_ratio_sane))
})

test_that("the fit-time CL2 adjustment is exact and recorded", {
  fixture <- make_exactcl2_fixture("gaussian", 4L, "influential")
  fit <- suppressWarnings(suppressMessages(pffr(
    Y ~ xlin,
    data = fixture$data,
    yind = fixture$yind,
    bs.yindex = list(bs = "ps", k = fixture$k, m = c(2, 1))
  )))
  expect_identical(fit$pffr$sandwich, "cl2")
  expect_identical(fit$pffr$sandwich_info$cl2_adjustment, "exact")
  expect_true(is.numeric(fit$pffr$sandwich_info$n_adjusted))
  expect_equal(
    fit$pffr$Vsandwich,
    suppressWarnings(refund:::gam_sandwich_cluster_cl2(
      refund:::pffr_model_based_gam(fit),
      refund:::build_cluster_id(fit$pffr)
    ))
  )
})

# ---------------------------------------------------------------------------
# Study-LB P-LB5 follow-up (2026-08-04): numerically exploded interval widths
# on degenerate Poisson fits must be DETECTED, not returned silently.
# ---------------------------------------------------------------------------

test_that("pffr_hat_invariant_violation() flags only impossible hat values", {
  viol <- refund:::pffr_hat_invariant_violation

  # Every invariant holds -> NULL.
  expect_null(viol(max_obs_leverage = 0.99))
  expect_null(viol(max_obs_leverage = NA_real_))
  # Round-off-sized excess is tolerated.
  expect_null(viol(max_obs_leverage = 1 + 1e-9))
  # Real violations are described, with the offending number in the message.
  expect_match(
    viol(max_obs_leverage = 22.6),
    "per-observation leverage is 22.6"
  )
  expect_match(viol(min_hat_eig = -0.5), "hat eigenvalue")
  expect_match(viol(min_block_eig_rel = -0.3), "positive semi-definite")
  # A non-finite monitor must not trip the check.
  expect_null(viol(max_obs_leverage = Inf, min_hat_eig = NaN))
})

test_that("the exact leverage floor engages silently", {
  skip_on_cran()
  # The "influential" design saturates one cluster's leverage; the exact path
  # floors the residual-block eigenvalues at (1 - cap)^2, a routine numerical
  # safeguard that does not warn.
  fixture <- make_exactcl2_fixture("poisson", 4L, "influential")
  fit <- suppressWarnings(suppressMessages(pffr(
    Y ~ xlin,
    data = fixture$data,
    yind = fixture$yind,
    family = fixture$family,
    bs.yindex = list(bs = "ps", k = fixture$k, m = c(2, 1)),
    sandwich = FALSE
  )))
  b <- refund:::pffr_model_based_gam(fit)
  cid <- refund:::build_cluster_id(fit$pffr, cluster = fixture$cluster)
  V_exact <- expect_no_warning(refund:::gam_sandwich_cluster_cl2(b, cid))
  expect_gt(attr(V_exact, "n_adjusted"), 0)
  expect_null(attr(V_exact, "hat_invariant_violation"))
})

test_that("a degenerate fit warns exactly once instead of a silent 1e30", {
  skip_on_cran()
  # amp = 10 drives individual fitted means many orders of magnitude above the
  # nominal marginal mean, the study-LB mechanism; the penalized hat then
  # breaks its own bounds by dozens of orders of magnitude.
  fit <- fit_lb5_fixture(make_lb5_fixture(amp = 10, n_grid = 30L, k = 12L, 21L))
  w <- capture_warnings(V <- refund:::pffr_vcov(fit, sandwich = TRUE))
  expect_length(w, 1L)
  expect_match(w, "NOT trustworthy")
  expect_gt(attr(V, "max_obs_leverage"), 1)
  expect_type(attr(V, "hat_invariant_violation"), "character")
  expect_true(is.finite(sqrt(max(abs(diag(V))))))
})

test_that("a benign Poisson CL2 fit reports leverages inside the bounds", {
  skip_on_cran()
  fit <- fit_lb5_fixture(make_lb5_fixture(amp = 1, n_grid = 20L, k = 8L, 21L))
  V <- expect_no_warning(refund:::pffr_vcov(fit, sandwich = TRUE))
  expect_lte(attr(V, "max_obs_leverage"), 1)
  expect_gte(attr(V, "min_obs_leverage"), 0)
  expect_null(attr(V, "hat_invariant_violation"))
})
