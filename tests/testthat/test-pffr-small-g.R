#--------------------------------------
# Small-cluster-count warning (plan S-C, amended 2026-09-08: threshold is
# G < 40, the paper's own recommendation boundary for pffr_coefboot(), not
# the earlier G < 20 draft). pffr() warns once, at fit time, when the
# resolved sandwich is "cluster" or "cl2" and the number of clusters G is
# below 40. The warning carries class "pffr_small_G_warning" so it can be
# muffled by class.
#--------------------------------------

make_smallg_dat <- function(n, seed = 5150) {
  pffr_simulate(
    Y ~ xlin,
    n = n,
    nygrid = 25,
    SNR = 5,
    effects = list(xlin = "dnorm"),
    seed = seed
  )
}

test_that("cl2 at G < 40 emits exactly one pffr_small_G_warning", {
  skip_on_cran()
  dat <- make_smallg_dat(n = 15)
  yind <- attr(dat, "yindex")

  warnings_seen <- character(0)
  fit <- withCallingHandlers(
    pffr(Y ~ xlin, data = dat, yind = yind, sandwich = "cl2"),
    pffr_small_G_warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_length(warnings_seen, 1)
  expect_match(warnings_seen, "Only G = 15 clusters")
  expect_match(warnings_seen, "pffr_coefboot")
  expect_identical(fit$pffr$sandwich_info$G, 15L)
})

test_that("cluster (CR1) at G < 40 emits exactly one pffr_small_G_warning", {
  skip_on_cran()
  dat <- make_smallg_dat(n = 12)
  yind <- attr(dat, "yindex")

  expect_warning(
    fit <- pffr(Y ~ xlin, data = dat, yind = yind, sandwich = "cluster"),
    class = "pffr_small_G_warning"
  )
  expect_identical(fit$pffr$sandwich_info$G, 12L)
})

test_that("cl2 at G >= 40 emits no pffr_small_G_warning", {
  skip_on_cran()
  dat <- make_smallg_dat(n = 45)
  yind <- attr(dat, "yindex")

  fired <- FALSE
  withCallingHandlers(
    pffr(Y ~ xlin, data = dat, yind = yind, sandwich = "cl2"),
    pffr_small_G_warning = function(w) {
      fired <<- TRUE
      invokeRestart("muffleWarning")
    }
  )
  expect_false(fired)
})

test_that("sandwich = 'none' at small G never warns", {
  skip_on_cran()
  dat <- make_smallg_dat(n = 10)
  yind <- attr(dat, "yindex")

  fired <- FALSE
  withCallingHandlers(
    pffr(Y ~ xlin, data = dat, yind = yind, sandwich = "none"),
    pffr_small_G_warning = function(w) {
      fired <<- TRUE
      invokeRestart("muffleWarning")
    }
  )
  expect_false(fired)
})

test_that("coef() on an already-warned small-G fit does not repeat the warning", {
  skip_on_cran()
  dat <- make_smallg_dat(n = 15)
  yind <- attr(dat, "yindex")

  fit <- withCallingHandlers(
    pffr(Y ~ xlin, data = dat, yind = yind, sandwich = "cl2"),
    pffr_small_G_warning = function(w) invokeRestart("muffleWarning")
  )

  fired <- FALSE
  withCallingHandlers(
    coef(fit),
    pffr_small_G_warning = function(w) {
      fired <<- TRUE
      invokeRestart("muffleWarning")
    }
  )
  expect_false(fired)
})
