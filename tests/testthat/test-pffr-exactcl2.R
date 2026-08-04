# A18 exact-CL2 hardening gate. The standalone script writes the same table
# for release records; this fast test keeps its family/design coverage in CI.
# Fixtures live in helper-exactcl2.R (they used to be sourced from a top-level
# script that was removed from the package tree, which left this file unable to
# run at all).

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

test_that("CL2 adjustment selection can be forced and is recorded", {
  fixture <- make_exactcl2_fixture("gaussian", 4L, "influential")
  fit_auto <- suppressWarnings(suppressMessages(pffr(
    Y ~ xlin,
    data = fixture$data,
    yind = fixture$yind,
    bs.yindex = list(bs = "ps", k = fixture$k, m = c(2, 1)),
    sandwich = "cl2"
  )))
  fit_exact <- suppressWarnings(suppressMessages(pffr(
    Y ~ xlin,
    data = fixture$data,
    yind = fixture$yind,
    bs.yindex = list(bs = "ps", k = fixture$k, m = c(2, 1)),
    sandwich = "cl2",
    cl2_adjustment = "exact"
  )))
  fit_shortcut <- suppressWarnings(suppressMessages(pffr(
    Y ~ xlin,
    data = fixture$data,
    yind = fixture$yind,
    bs.yindex = list(bs = "ps", k = fixture$k, m = c(2, 1)),
    sandwich = "cl2",
    cl2_adjustment = "shortcut"
  )))

  expect_identical(fit_exact$pffr$sandwich_info$cl2_adjustment, "exact")
  expect_identical(fit_shortcut$pffr$sandwich_info$cl2_adjustment, "shortcut")
  expect_identical(fit_auto$pffr$sandwich_info$cl2_adjustment, "exact")
  expect_equal(fit_auto$pffr$Vsandwich, fit_exact$pffr$Vsandwich)
  expect_true(is.numeric(fit_exact$pffr$sandwich_info$n_adjusted))
  expect_true(is.numeric(fit_shortcut$pffr$sandwich_info$n_adjusted))
})

test_that("automatic exact CL2 uses the documented relevance and cost rule", {
  resolve <- refund:::resolve_cl2_adjustment
  expect_identical(resolve("auto", G = 100, maxDg = 50, p = 20), "exact")
  expect_identical(resolve("auto", G = 101, maxDg = 50, p = 20), "shortcut")
  # 9e8 ops: below the 5e9 cost cap (raised from 5e8, see notes 2026-07-21)
  expect_identical(resolve("auto", G = 20, maxDg = 500, p = 300), "exact")
  # 6.4e10 ops: above the cap even at eligible G
  expect_identical(resolve("auto", G = 100, maxDg = 1000, p = 800), "shortcut")
  expect_identical(resolve("exact", G = 1000, maxDg = 5000, p = 500), "exact")
  expect_identical(resolve("shortcut", G = 2, maxDg = 2, p = 2), "shortcut")
})

# ---------------------------------------------------------------------------
# Study-LB P-LB5 follow-up (2026-08-04): numerically exploded interval widths
# on degenerate Poisson fits must be DETECTED, not returned silently.
# ---------------------------------------------------------------------------

test_that("pffr_hat_invariant_violation() flags only impossible hat values", {
  viol <- refund:::pffr_hat_invariant_violation

  # Every invariant holds -> NULL.
  expect_null(viol(max_obs_leverage = 0.99, max_leverage = 0.999))
  expect_null(viol(max_obs_leverage = NA_real_, max_leverage = NA_real_))
  # Round-off-sized excess is tolerated.
  expect_null(viol(max_obs_leverage = 1 + 1e-9))
  # Real violations are described, with the offending number in the message.
  expect_match(
    viol(max_obs_leverage = 22.6),
    "per-observation leverage is 22.6"
  )
  expect_match(viol(max_leverage = 8.18e84), "hat eigenvalue")
  expect_match(viol(min_block_eig_rel = -0.3), "positive semi-definite")
  # A non-finite monitor must not trip the check.
  expect_null(viol(max_obs_leverage = Inf, max_leverage = NaN))
})

test_that("a benign Poisson CL2 fit reports leverages inside the bounds", {
  skip_on_cran()
  fit <- fit_lb5_fixture(make_lb5_fixture(amp = 1, n_grid = 20L, k = 8L, 21L))
  V_exact <- expect_no_warning(
    refund:::pffr_vcov(fit, sandwich = "cl2", cl2_adjustment = "exact")
  )
  V_short <- expect_no_warning(
    refund:::pffr_vcov(fit, sandwich = "cl2", cl2_adjustment = "shortcut")
  )
  expect_lte(attr(V_exact, "max_obs_leverage"), 1)
  expect_lte(attr(V_short, "max_leverage"), 1)
  expect_null(attr(V_exact, "hat_invariant_violation"))
  expect_null(attr(V_short, "hat_invariant_violation"))
  # exact and shortcut agree closely in the well-behaved regime
  ratio <- sqrt(max(abs(diag(V_exact)))) / sqrt(max(abs(diag(V_short))))
  expect_gt(ratio, 0.5)
  expect_lt(ratio, 2)
})

test_that("a degenerate Poisson fit warns instead of returning a silent 1e30", {
  skip_on_cran()
  # amp = 10 drives individual fitted means many orders of magnitude above the
  # nominal marginal mean, the study-LB mechanism; the penalized hat then
  # breaks its own bounds by dozens of orders of magnitude.
  fit <- fit_lb5_fixture(make_lb5_fixture(amp = 10, n_grid = 30L, k = 12L, 21L))

  expect_warning(
    V_exact <- refund:::pffr_vcov(
      fit,
      sandwich = "cl2",
      cl2_adjustment = "exact"
    ),
    "NOT trustworthy"
  )
  expect_warning(
    V_short <- refund:::pffr_vcov(
      fit,
      sandwich = "cl2",
      cl2_adjustment = "shortcut"
    ),
    "NOT trustworthy"
  )
  expect_gt(attr(V_exact, "max_obs_leverage"), 1)
  expect_type(attr(V_exact, "hat_invariant_violation"), "character")

  se_exact <- sqrt(max(abs(diag(V_exact))))
  se_short <- sqrt(max(abs(diag(V_short))))
  expect_true(is.finite(se_exact) && is.finite(se_short))
  # The exact block is a principal block of a squared symmetric matrix and so
  # stays positive semi-definite; the shortcut's I - H_gg turns indefinite and
  # is floored into the largest legitimate variance inflation. The exact path
  # is therefore far less explosive here, though neither is meaningful.
  expect_lt(se_exact, se_short)
})
