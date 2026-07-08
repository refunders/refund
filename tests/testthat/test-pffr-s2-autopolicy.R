#--------------------------------------
# Tests for sandwich = "auto" and the S2 resolution policy
#
# pffr_sandwich_auto_policy() / family_has_exact_score(): the family-score
# dispatch and the (G, max D_g) thresholds that promote a fit to CL2, the
# options(refund.pffr.autopolicy=) override, and the fit-time resolution +
# one-line message in pffr(). The factory DEFAULT is unchanged (still
# "cluster"); "auto" is opt-in pending the PI decision.
#--------------------------------------

pol <- function(...) refund:::pffr_sandwich_auto_policy(...)
fhes <- function(...) refund:::family_has_exact_score(...)

test_that("family_has_exact_score follows the actual sandwich dispatch", {
  # standard exponential-dispersion GLM families: exact score
  expect_true(fhes(gaussian()))
  expect_true(fhes(poisson()))
  expect_true(fhes(binomial()))
  expect_true(fhes(Gamma()))
  # gaulss: two-block Fisher-whitened exact score (special-cased)
  expect_true(fhes(mgcv::gaulss()))
  # extended families reuse the working-residual APPROXIMATION -> not exact
  expect_false(fhes(mgcv::scat()))
  expect_false(fhes(mgcv::nb()))
  # a family defining $sandwich (other than gaulss) -> HC fallback, no score path
  fake <- gaussian()
  fake$family <- "multinom"
  fake$sandwich <- function() NULL
  expect_false(fhes(fake))
  # NULL family
  expect_false(fhes(NULL))
})

test_that("policy resolves each branch as specified", {
  # eligible: exact-score family, moderate G, small D_g -> cl2
  expect_identical(pol(92, 55, gaussian()), "cl2")
  expect_identical(pol(80, 120, poisson()), "cl2")
  expect_identical(pol(150, 500, gaussian()), "cl2") # boundaries inclusive
  # G above threshold -> cluster
  expect_identical(pol(151, 55, gaussian()), "cluster")
  expect_identical(pol(200, 55, gaussian()), "cluster")
  # max D_g above threshold -> cluster
  expect_identical(pol(92, 501, gaussian()), "cluster")
  # family without an exact/two-block score path -> cluster
  expect_identical(pol(92, 55, mgcv::scat()), "cluster")
  expect_identical(pol(92, 55, mgcv::nb()), "cluster")
  # non-finite guards -> cluster
  expect_identical(pol(NA_real_, 55, gaussian()), "cluster")
  expect_identical(pol(92, NA_real_, gaussian()), "cluster")
})

test_that("DTI-shaped case (G=92, max D_g=55) resolves to cl2", {
  # Makes the paper's application consistent with its own rule (review C2).
  expect_identical(pol(92, 55, gaussian()), "cl2")
})

test_that("options(refund.pffr.autopolicy=) overrides the policy", {
  # named-list threshold override
  withr::with_options(
    list(refund.pffr.autopolicy = list(G_max = 80, Dg_max = 300)),
    {
      expect_identical(pol(92, 55, gaussian()), "cluster") # G>80 now
      expect_identical(pol(80, 55, gaussian()), "cl2")
      expect_identical(pol(50, 400, gaussian()), "cluster") # Dg>300 now
    }
  )
  # full-function override
  withr::with_options(
    list(refund.pffr.autopolicy = function(G, maxDg, family) "cl2"),
    expect_identical(pol(999, 999, mgcv::scat()), "cl2")
  )
  withr::with_options(
    list(refund.pffr.autopolicy = function(G, maxDg, family) "cluster"),
    expect_identical(pol(5, 5, gaussian()), "cluster")
  )
})

# One shared small fit reused across the fit-time tests.
make_auto_data <- local({
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
      cache <<- list(dat = dat, yind = attr(dat, "yindex"))
    }
    cache
  }
})

test_that("sandwich='auto' resolves at fit time and messages once", {
  skip_on_cran()
  d <- make_auto_data()
  expect_message(
    fit <- pffr(Y ~ ff(X1), data = d$dat, yind = d$yind, sandwich = "auto"),
    "sandwich='auto' resolved to 'cl2' \\(G=25, max D_g=20\\)"
  )
  # the resolved type is stored and served
  expect_identical(fit$pffr$sandwich_info$type, "cl2")
  expect_identical(fit$pffr$sandwich, "cl2")
})

test_that("auto-resolved cl2 equals an explicit cl2 fit; explicit values intact", {
  skip_on_cran()
  d <- make_auto_data()
  fit_auto <- suppressMessages(
    pffr(Y ~ ff(X1), data = d$dat, yind = d$yind, sandwich = "auto")
  )
  fit_cl2 <- suppressMessages(
    pffr(Y ~ ff(X1), data = d$dat, yind = d$yind, sandwich = "cl2")
  )
  expect_equal(
    fit_auto$pffr$Vsandwich,
    fit_cl2$pffr$Vsandwich,
    tolerance = 1e-12
  )
  # explicit choices still work and are stored verbatim
  fit_cluster <- suppressMessages(
    pffr(Y ~ ff(X1), data = d$dat, yind = d$yind, sandwich = "cluster")
  )
  expect_identical(fit_cluster$pffr$sandwich_info$type, "cluster")
  fit_none <- suppressMessages(
    pffr(Y ~ ff(X1), data = d$dat, yind = d$yind, sandwich = "none")
  )
  expect_null(fit_none$pffr$sandwich_info)
  # $Vp is model-based and identical across sandwich choices (S1 contract)
  expect_identical(fit_auto$Vp, fit_none$Vp)
})

test_that("an autopolicy that forces cluster is honored at fit time", {
  skip_on_cran()
  d <- make_auto_data()
  withr::with_options(
    list(refund.pffr.autopolicy = list(G_max = 10)),
    {
      fit <- suppressMessages(
        pffr(Y ~ ff(X1), data = d$dat, yind = d$yind, sandwich = "auto")
      )
      expect_identical(fit$pffr$sandwich_info$type, "cluster")
    }
  )
})
