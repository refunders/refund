# A18 exact-CL2 hardening gate. The standalone script writes the same table
# for release records; this fast test keeps its family/design coverage in CI.

source(testthat::test_path("..", "..", "hardening-exactcl2.R"), local = TRUE)

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
  expect_identical(resolve("auto", G = 20, maxDg = 500, p = 300), "shortcut")
  expect_identical(resolve("exact", G = 1000, maxDg = 5000, p = 500), "exact")
  expect_identical(resolve("shortcut", G = 2, maxDg = 2, p = 2), "shortcut")
})
