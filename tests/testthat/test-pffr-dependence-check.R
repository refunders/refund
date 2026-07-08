#--------------------------------------
# S5: within-curve dependence diagnostic
#--------------------------------------
#
# pffr_dependence_check() is a DESCRIPTIVE flag (not a test/estimator) telling
# users which inference regime they are in. These tests check that AR1(0.9)
# within-curve errors produce a large design effect (DE >> 1.5, "detected"),
# that iid errors give DE ~ 1 ("weak"), and that an irregular-grid fit does not
# error.

# AR1 noise matrix (autocorrelation along the functional index, i.e. columns)
ar1_noise_matrix <- function(n, D, rho, sd = 1) {
  noise <- matrix(0, n, D)
  innov_sd <- sd * sqrt(1 - rho^2)
  for (i in seq_len(n)) {
    e <- numeric(D)
    e[1] <- rnorm(1, sd = sd)
    for (j in 2:D) {
      e[j] <- rho * e[j - 1] + rnorm(1, sd = innov_sd)
    }
    noise[i, ] <- e
  }
  noise
}


test_that("AR1(0.9) within-curve errors give a large design effect", {
  skip_on_cran()
  set.seed(30101)
  dat <- pffr_simulate(
    Y ~ xlin,
    n = 40,
    nygrid = 30,
    effects = list(xlin = "dnorm"),
    SNR = 50 # tiny iid noise, so the injected AR1 dominates the residuals
  )
  t <- attr(dat, "yindex")
  dat$Y <- dat$Y +
    ar1_noise_matrix(nrow(dat$Y), ncol(dat$Y), rho = 0.9, sd = 1)

  m <- pffr(
    Y ~ xlin,
    yind = t,
    data = dat,
    bs.yindex = list(bs = "ps", k = 12, m = c(2, 1)),
    sandwich = "none"
  )
  d <- pffr_dependence_check(m)

  expect_s3_class(d, "pffr_dependence_check")
  expect_gt(d$DE, 2.5) # well above the 1.5 threshold
  expect_identical(d$regime, "detected")
  # Wrong-order regression guard: pffr stacks the response curve-major
  # (as.vector(t(Y)), see pffr:::pffr_setup_response). A t-major misread would
  # scramble curves and collapse the per-curve lag-1 acf to ~0. Correct
  # alignment recovers a strongly positive rho1 (~0.6-0.8 here).
  expect_gt(d$rho1_mean, 0.5)
  expect_match(d$advisory, "anti-conservative")
  # N_eff = N / DE and stays below the raw sample size
  expect_equal(d$N_eff, d$N / d$DE, tolerance = 1e-10)
  expect_lt(d$N_eff, d$N)
})


test_that("iid within-curve errors give a design effect near 1", {
  skip_on_cran()
  set.seed(30202)
  dat <- pffr_simulate(
    Y ~ xlin,
    n = 40,
    nygrid = 30,
    effects = list(xlin = "dnorm"),
    SNR = 4
  )
  t <- attr(dat, "yindex")
  m <- pffr(
    Y ~ xlin,
    yind = t,
    data = dat,
    bs.yindex = list(bs = "ps", k = 12, m = c(2, 1)),
    sandwich = "none"
  )
  d <- pffr_dependence_check(m)

  expect_lt(d$DE, 1.5) # generous slack around DE ~ 1
  expect_gte(d$DE, 1) # DE is bounded below by 1 by construction
  expect_identical(d$regime, "weak")
  expect_lt(abs(d$rho1_mean), 0.15) # near-zero lag-1 autocorrelation
  expect_match(d$advisory, "weak")
})


test_that("iid errors on a dense grid (Dbar ~ 80) keep DE ~ 1 (noise-floor guard)", {
  skip_on_cran()
  # The literal "mean |acf|" design effect has a positive noise floor
  # ~sqrt(2/(pi D)) per lag that GROWS with the grid: at Dbar ~ 80 it would give
  # DE ~ 8 on independent data. Averaging SIGNED autocorrelations across curves
  # cancels that floor, so DE must stay near 1 even on a fine grid.
  set.seed(30205)
  dat <- pffr_simulate(
    Y ~ xlin,
    n = 40,
    nygrid = 80,
    effects = list(xlin = "dnorm"),
    SNR = 6
  )
  t <- attr(dat, "yindex")
  m <- pffr(
    Y ~ xlin,
    yind = t,
    data = dat,
    bs.yindex = list(bs = "ps", k = 12, m = c(2, 1)),
    sandwich = "none"
  )
  d <- pffr_dependence_check(m)
  expect_gt(d$Dbar, 70) # confirm the fine grid
  expect_lt(d$DE, 1.5) # noise floor did not inflate DE
  expect_identical(d$regime, "weak")
})


test_that("irregular-grid (sparse) fit does not error", {
  skip_on_cran()
  set.seed(30303)
  sp <- pffr_simulate(
    Y ~ s(xsmoo),
    n = 30,
    nygrid = 25,
    effects = list(xsmoo = "sine"),
    propmissing = 0.4,
    SNR = 8
  )
  ts <- attr(sp, "yindex")
  m <- pffr(
    Y ~ s(xsmoo),
    data = sp$data,
    ydata = sp$ydata,
    yind = ts,
    sandwich = "none"
  )

  expect_no_error(d <- pffr_dependence_check(m))
  expect_s3_class(d, "pffr_dependence_check")
  expect_true(is.finite(d$DE))
  expect_true(is.finite(d$N_eff))
  expect_gte(d$DE, 1)
  expect_true(d$G >= 1)
})


test_that("print method and summary() one-liner show the dependence flag", {
  skip_on_cran()
  set.seed(30404)
  dat <- pffr_simulate(
    Y ~ xlin,
    n = 30,
    nygrid = 25,
    effects = list(xlin = "dnorm"),
    SNR = 6
  )
  t <- attr(dat, "yindex")
  m <- pffr(
    Y ~ xlin,
    yind = t,
    data = dat,
    bs.yindex = list(bs = "ps", k = 12, m = c(2, 1)),
    sandwich = "none"
  )

  d <- pffr_dependence_check(m)
  out <- capture.output(print(d))
  expect_true(any(grepl("Within-curve dependence", out)))
  expect_true(any(grepl("design effect DE", out)))
  expect_true(any(grepl("vs G =", out)))

  # summary.pffr carries the flag and its printer emits the same headline
  sm <- summary(m)
  expect_s3_class(sm$dependence, "pffr_dependence_check")
  sout <- capture.output(print(sm))
  expect_true(any(grepl("Within-curve dependence", sout)))
})
