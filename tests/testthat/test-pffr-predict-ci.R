# Pointwise intervals for predictions: pffr_predict_ci().

predict_ci_data <- function(family = gaussian()) {
  dat <- if (family$family == "gaussian") {
    ncv_test_data()
  } else {
    ncv_test_glm_data(family)
  }
  dat$z <- rnorm(nrow(dat))
  dat
}

predict_ci_fit <- function(dat, family = gaussian(), ...) {
  quiet_pffr(
    Y ~
      ff(
        X,
        splinepars = list(bs = "ps", m = list(c(2, 1), c(2, 1)), k = c(5, 5))
      ) +
        c(z),
    data = dat,
    yind = seq(0, 1, length.out = ncol(dat$Y)),
    family = family,
    bs.yindex = list(bs = "ps", k = 6, m = c(2, 1)),
    bs.int = list(bs = "ps", k = 6, m = c(2, 1)),
    ...
  )
}

expect_direct_intervals <- function(family) {
  set.seed(4211)
  dat <- predict_ci_data(family)
  fit <- predict_ci_fit(dat, family)
  newdata <- list(X = I(dat$X[1:4, ]), z = dat$z[1:4])
  L <- predict(fit, newdata, type = "lpmatrix", reformat = FALSE)
  V <- fit$pffr$Vsandwich
  est <- as.vector(L %*% fit$coefficients)
  se <- unname(sqrt(rowSums((L %*% V) * L)))
  df <- pffr_influence_df(pffr_influence(fit), L)$df
  half <- qt(0.975, df) * se

  link <- pffr_predict_ci(fit, newdata = newdata, type = "link")
  response <- pffr_predict_ci(fit, newdata = newdata, type = "response")
  expect_equal(link$fit, est, tolerance = 1e-10)
  expect_equal(link$se_link, se, tolerance = 1e-10)
  expect_equal(link$df, df, tolerance = 1e-10)
  expect_equal(link$lower, est - half, tolerance = 1e-10)
  expect_equal(link$upper, est + half, tolerance = 1e-10)
  linkinv <- fit$family$linkinv
  expect_equal(response$fit, linkinv(est), tolerance = 1e-10)
  expect_equal(response$lower, linkinv(est - half), tolerance = 1e-10)
  expect_equal(response$upper, linkinv(est + half), tolerance = 1e-10)
  expect_equal(response$se_link, link$se_link)
  expect_equal(nrow(link), 4 * ncol(dat$Y))
  expect_equal(link$.obs, rep(1:4, each = ncol(dat$Y)))
  expect_equal(link$.index, rep(fit$pffr$yind, 4))
  expect_identical(attr(link, "crit_used"), "satterthwaite")
}

test_that("pffr_predict_ci reproduces CL2 + Satterthwaite (Gaussian)", {
  skip_on_cran()
  expect_direct_intervals(gaussian())
})

test_that("pffr_predict_ci reproduces CL2 + Satterthwaite (Poisson)", {
  skip_on_cran()
  expect_direct_intervals(poisson())
})

test_that("pffr_predict_ci model-based intervals use Vc (REML) or Vp (NCV)", {
  skip_on_cran()
  skip_if_not_installed("mgcv", "1.9.0")
  set.seed(4214)
  dat <- predict_ci_data()
  fit <- predict_ci_fit(dat)
  L <- predict(fit, type = "lpmatrix", reformat = FALSE)
  ci <- pffr_predict_ci(fit, sandwich = FALSE, level = 0.9)
  se <- unname(sqrt(rowSums((L %*% fit$Vc) * L)))
  expect_equal(ci$se_link, se, tolerance = 1e-10)
  expect_equal(ci$upper - ci$fit, qnorm(0.95) * se, tolerance = 1e-10)
  expect_equal(ci$fit, as.vector(fit$linear.predictors))
  expect_true(all(is.infinite(ci$df)))

  fit_ncv <- predict_ci_fit(dat, method = "NCV", sandwich = FALSE)
  rm(list = intersect("ncv_intervals", ls(.pffr_state)), envir = .pffr_state)
  expect_message(
    ci_ncv <- pffr_predict_ci(fit_ncv),
    "Intervals around NCV estimates undercover"
  )
  # once per session
  expect_no_message(pffr_predict_ci(fit_ncv))
  expect_equal(ci_ncv$se_link, unname(sqrt(rowSums((L %*% fit_ncv$Vp) * L))))
})

test_that("pffr_predict_ci at zero covariates gives the full intercept alpha(t)", {
  skip_on_cran()
  set.seed(4219)
  dat <- predict_ci_data()
  fit <- predict_ci_fit(dat)
  zero <- list(X = I(matrix(0, 1, ncol(dat$X))), z = 0)
  alpha <- pffr_predict_ci(fit, newdata = zero)
  # The full intercept's rows: the scalar intercept plus the centred
  # Intercept(yindex) basis at the response grid.
  sm <- fit$smooth[["s(yindex.vec)"]]
  yind <- fit$pffr$yind
  L <- matrix(0, length(yind), length(fit$coefficients))
  L[, sm$first.para:sm$last.para] <- mgcv::PredictMat(
    sm,
    data.frame(yindex.vec = yind)
  )
  L[, names(fit$coefficients) == "(Intercept)"] <- 1
  se <- unname(sqrt(rowSums((L %*% fit$pffr$Vsandwich) * L)))
  df <- pffr_influence_df(pffr_influence(fit), L)$df
  expect_equal(alpha$.index, yind)
  expect_equal(alpha$fit, as.vector(L %*% fit$coefficients), tolerance = 1e-10)
  expect_equal(alpha$se_link, se, tolerance = 1e-10)
  expect_equal(alpha$df, df, tolerance = 1e-10)
  expect_equal(alpha$upper - alpha$fit, qt(0.975, df) * se, tolerance = 1e-10)
  # Its estimate is the centred term plus the scalar intercept from coef().
  cf <- suppressMessages(coef(
    fit,
    eval_grid = list(`Intercept(yindex)` = data.frame(yindex.vec = yind))
  ))
  expect_equal(
    alpha$fit,
    cf$smterms[["Intercept(yindex)"]]$coef$value +
      cf$pterms["(Intercept)", "value"],
    tolerance = 1e-10
  )
})

test_that("pffr_predict_ci aligns fitted points with missing responses", {
  skip_on_cran()
  set.seed(4215)
  dat <- predict_ci_data()
  dat$Y[cbind(c(2, 5, 5), c(3, 1, 16))] <- NA
  fit <- predict_ci_fit(dat)
  ci <- pffr_predict_ci(fit)
  n_grid <- ncol(dat$Y)
  expect_equal(nrow(ci), length(dat$Y) - 3)
  observed <- which(!is.na(t(dat$Y)))
  expect_equal(ci$.obs, (observed - 1) %/% n_grid + 1)
  expect_equal(ci$.index, fit$pffr$yind[(observed - 1) %% n_grid + 1])
  expect_equal(ci$fit, as.vector(fit$linear.predictors))
  expect_true(all(is.finite(ci$df)))
})

test_that("pffr_predict_ci aligns sparse responses, skipping missing values", {
  skip_on_cran()
  set.seed(4217)
  n <- 16
  dat <- data.frame(z = rnorm(n))
  yd <- data.frame(
    .obs = rep(seq_len(n), times = rep(c(6, 8, 10, 12), n / 4)),
    .index = runif(n * 9)
  )
  yd$.value <- dat$z[yd$.obs] * sin(2 * pi * yd$.index) + rnorm(nrow(yd))
  yd <- yd[sample(nrow(yd)), ]
  yd$.value[c(3, 40)] <- NA
  # The default CL2 sandwich clusters the retained (non-missing) rows by curve.
  fit <- quiet_pffr(
    Y ~ z,
    data = dat,
    ydata = yd,
    bs.yindex = list(bs = "ps", k = 5, m = c(2, 1)),
    bs.int = list(bs = "ps", k = 5, m = c(2, 1))
  )
  expect_identical(fit$pffr$sandwich, "cl2")
  ci <- pffr_predict_ci(fit)
  expect_true(all(is.finite(ci$df)))
  kept <- yd[!is.na(yd$.value), ]
  expect_equal(nrow(ci), nrow(kept))
  expect_equal(ci$.obs, kept$.obs)
  expect_equal(ci$.index, kept$.index)
  expect_equal(ci$fit, as.vector(fit$linear.predictors))
})

test_that("pffr_predict_ci rejects multi-linear-predictor families", {
  skip_on_cran()
  fit <- get_gaulss_model()
  expect_error(pffr_predict_ci(fit), "single linear-predictor")
})
