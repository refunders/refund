# Bias-aware pointwise intervals: coef.pffr(bias_ref = ) and pffr_predict_ci().

bias_aware_data <- function(family = gaussian()) {
  dat <- if (family$family == "gaussian") {
    ncv_test_data()
  } else {
    ncv_test_glm_data(family)
  }
  dat$z <- rnorm(nrow(dat))
  dat
}

bias_aware_fit <- function(dat, method, family = gaussian()) {
  pffr(
    Y ~
      ff(
        X,
        splinepars = list(bs = "ps", m = list(c(2, 1), c(2, 1)), k = c(5, 5))
      ) +
        c(z),
    data = dat,
    yind = seq(0, 1, length.out = ncol(dat$Y)),
    family = family,
    method = method,
    bs.yindex = list(bs = "ps", k = 6, m = c(2, 1)),
    bs.int = list(bs = "ps", k = 6, m = c(2, 1)),
    sandwich = "none"
  )
}

# Exact CL2 covariance in Bayesian form, as scored by the study.
bias_aware_cl2 <- function(fit) {
  V <- pffr_vcov(fit, sandwich = "cl2", cl2_adjustment = "exact", freq = FALSE)
  expect_identical(attr(V, "cl2_adjustment"), "exact")
  V
}

# The study's interval from rows L of an estimand: NCV estimate, CL2 SE of the
# NCV fit and the NCV - REML difference in quadrature, z critical value.
bias_aware_direct <- function(L, fit_ncv, fit_reml, V, level = 0.95) {
  est <- as.vector(L %*% fit_ncv$coefficients)
  se <- unname(sqrt(rowSums((L %*% V) * L)))
  delta <- as.vector(L %*% (fit_ncv$coefficients - fit_reml$coefficients))
  half <- qnorm((1 + level) / 2) * sqrt(se^2 + delta^2)
  list(
    est = est,
    se = se,
    delta = delta,
    lower = est - half,
    upper = est + half
  )
}

# Rows of the ff coefficient surface on the grid coef.pffr() used (by = 1, no
# intercept columns: the seWithMean = FALSE convention).
bias_aware_ff_rows <- function(fit, grid) {
  sm <- fit$smooth[[grep("X.smat", names(fit$smooth), fixed = TRUE)]]
  L <- matrix(0, nrow(grid), length(fit$coefficients))
  L[, sm$first.para:sm$last.para] <- mgcv::PredictMat(sm, grid)
  L
}

bias_aware_newdata <- function(dat, rows = 1:4) {
  list(X = I(dat$X[rows, ]), z = dat$z[rows])
}

expect_study_intervals <- function(family) {
  set.seed(4211)
  dat <- bias_aware_data(family)
  fit_ncv <- bias_aware_fit(dat, "NCV", family)
  fit_reml <- bias_aware_fit(dat, "REML", family)
  V <- bias_aware_cl2(fit_ncv)

  # Coefficient surface of the ff term.
  cf <- coef(
    fit_ncv,
    sandwich = "cl2",
    cl2_adjustment = "exact",
    ci = "pointwise",
    seWithMean = FALSE,
    bias_ref = fit_reml
  )
  ff_coef <- cf$smterms[["ff(X)"]]$coef
  ref <- bias_aware_direct(
    bias_aware_ff_rows(fit_ncv, ff_coef[, c("X.smat", "X.tmat", "L.X")]),
    fit_ncv,
    fit_reml,
    V
  )
  expect_gt(max(abs(ref$delta)), 1e-4)
  expect_equal(ff_coef$value, ref$est, tolerance = 1e-10)
  expect_equal(ff_coef$se, ref$se, tolerance = 1e-10)
  expect_equal(ff_coef$delta, ref$delta, tolerance = 1e-10)
  expect_equal(ff_coef$lower, ref$lower, tolerance = 1e-10)
  expect_equal(ff_coef$upper, ref$upper, tolerance = 1e-10)
  expect_identical(cf$ci_meta$bias_ref_method, "REML")

  # Parametric coefficient of z.
  j <- which(names(fit_ncv$coefficients) == "z")
  ref_z <- bias_aware_direct(
    diag(length(fit_ncv$coefficients))[j, , drop = FALSE],
    fit_ncv,
    fit_reml,
    V
  )
  expect_equal(
    unname(cf$pterms["z", c("value", "se", "delta", "lower", "upper")]),
    unlist(ref_z, use.names = FALSE),
    tolerance = 1e-10
  )

  # Conditional mean for new curves: link scale, then transformed endpoints.
  newdata <- bias_aware_newdata(dat)
  L <- predict(fit_ncv, newdata, type = "lpmatrix", reformat = FALSE)
  ref_mean <- bias_aware_direct(L, fit_ncv, fit_reml, V)
  args <- list(
    fit_ncv,
    newdata = newdata,
    sandwich = "cl2",
    cl2_adjustment = "exact",
    bias_ref = fit_reml
  )
  link <- do.call(pffr_predict_ci, c(args, type = "link"))
  response <- do.call(pffr_predict_ci, c(args, type = "response"))
  expect_equal(link$fit, ref_mean$est, tolerance = 1e-10)
  expect_equal(link$se_link, ref_mean$se, tolerance = 1e-10)
  expect_equal(link$delta_link, ref_mean$delta, tolerance = 1e-10)
  expect_equal(link$lower, ref_mean$lower, tolerance = 1e-10)
  expect_equal(link$upper, ref_mean$upper, tolerance = 1e-10)
  linkinv <- fit_ncv$family$linkinv
  expect_equal(response$fit, linkinv(ref_mean$est), tolerance = 1e-10)
  expect_equal(response$lower, linkinv(ref_mean$lower), tolerance = 1e-10)
  expect_equal(response$upper, linkinv(ref_mean$upper), tolerance = 1e-10)
  expect_equal(response$se_link, link$se_link)
  expect_equal(nrow(link), 4 * ncol(dat$Y))
  expect_equal(link$.obs, rep(1:4, each = ncol(dat$Y)))
  expect_equal(link$.index, rep(fit_ncv$pffr$yind, 4))
  expect_identical(attr(link, "bias_ref_method"), "REML")
}

test_that("bias-aware intervals reproduce the study's computation (Gaussian)", {
  skip_if_not_installed("mgcv", "1.9.0")
  expect_study_intervals(gaussian())
})

test_that("bias-aware intervals reproduce the study's computation (Poisson)", {
  skip_if_not_installed("mgcv", "1.9.0")
  expect_study_intervals(poisson())
})

test_that("a fit as its own bias reference gives delta = 0 and plain intervals", {
  skip_if_not_installed("mgcv", "1.9.0")
  set.seed(4212)
  dat <- bias_aware_data()
  fit <- bias_aware_fit(dat, "NCV")
  args <- list(
    fit,
    sandwich = "cl2",
    cl2_adjustment = "exact",
    ci = "pointwise"
  )
  plain <- suppressMessages(do.call(coef, args))
  self <- suppressMessages(do.call(coef, c(args, list(bias_ref = fit))))
  for (term in names(plain$smterms)) {
    expect_true(all(self$smterms[[term]]$coef$delta == 0))
    expect_equal(
      self$smterms[[term]]$coef[, c("value", "se", "lower", "upper")],
      plain$smterms[[term]]$coef[, c("value", "se", "lower", "upper")]
    )
  }
  expect_true(all(self$pterms[, "delta"] == 0))
  expect_equal(
    self$pterms[, c("lower", "upper")],
    plain$pterms[, c("lower", "upper")]
  )

  pred_plain <- pffr_predict_ci(fit, sandwich = "cl2", cl2_adjustment = "exact")
  pred_self <- pffr_predict_ci(
    fit,
    sandwich = "cl2",
    cl2_adjustment = "exact",
    bias_ref = fit
  )
  expect_true(all(pred_self$delta_link == 0))
  for (column in names(pred_plain)) {
    expect_equal(pred_self[[column]], pred_plain[[column]])
  }
  expect_false("delta_link" %in% names(pred_plain))
})

test_that("with seWithMean = TRUE delta uses the same linear map as the SE", {
  skip_if_not_installed("mgcv", "1.9.0")
  set.seed(4213)
  # Poisson: for a Gaussian fit the mean level (mean fitted value) is the same
  # for both fits, so the seWithMean shift would be zero.
  dat <- bias_aware_data(poisson())
  fit_ncv <- bias_aware_fit(dat, "NCV", poisson())
  fit_reml <- bias_aware_fit(dat, "REML", poisson())
  get_intercept <- function(seWithMean) {
    cf <- suppressMessages(coef(
      fit_ncv,
      sandwich = "cl2",
      cl2_adjustment = "exact",
      ci = "pointwise",
      seWithMean = seWithMean,
      bias_ref = fit_reml
    ))
    cf$smterms[["Intercept(yindex)"]]$coef
  }
  with_mean <- get_intercept(TRUE)
  without <- get_intercept(FALSE)
  sm <- fit_ncv$smooth[["s(yindex.vec)"]]
  expect_gt(attr(sm, "nCons"), 0)
  # The mean-level contrast: column means of the other model-matrix columns
  # (incl. the scalar intercept) times the coefficient difference.
  others <- setdiff(seq_along(fit_ncv$coefficients), sm$first.para:sm$last.para)
  dtheta <- fit_ncv$coefficients - fit_reml$coefficients
  shift <- sum(fit_ncv$cmX[others] * dtheta[others]) / (sm$meanL1 %||% 1)
  expect_gt(abs(shift), 1e-6)
  expect_equal(with_mean$delta - without$delta, rep(shift, nrow(with_mean)))
  expect_equal(
    with_mean$upper - with_mean$value,
    qnorm(0.975) * sqrt(with_mean$se^2 + with_mean$delta^2)
  )
})

test_that("pffr_predict_ci passes covariance options through", {
  skip_if_not_installed("mgcv", "1.9.0")
  set.seed(4214)
  dat <- bias_aware_data()
  fit <- bias_aware_fit(dat, "NCV")
  L <- predict(fit, type = "lpmatrix", reformat = FALSE)
  for (freq in c(FALSE, TRUE)) {
    V <- pffr_vcov(fit, sandwich = "cl2", cl2_adjustment = "exact", freq = freq)
    ci <- pffr_predict_ci(
      fit,
      sandwich = "cl2",
      cl2_adjustment = "exact",
      freq = freq,
      level = 0.9
    )
    se <- unname(sqrt(rowSums((L %*% V) * L)))
    expect_equal(ci$se_link, se, tolerance = 1e-10)
    expect_equal(ci$upper - ci$fit, qnorm(0.95) * se, tolerance = 1e-10)
    expect_equal(ci$fit, as.vector(fit$linear.predictors))
  }
  # Model-based covariance of an NCV fit is Vp.
  ci_model <- pffr_predict_ci(fit, sandwich = "none")
  expect_equal(ci_model$se_link, unname(sqrt(rowSums((L %*% fit$Vp) * L))))
})

test_that("bias_ref defaults to exact Bayesian CL2 unless overridden", {
  skip_if_not_installed("mgcv", "1.9.0")
  set.seed(4218)
  dat <- bias_aware_data()
  fit_ncv <- bias_aware_fit(dat, "NCV")
  fit_reml <- bias_aware_fit(dat, "REML")
  explicit <- list(sandwich = "cl2", cl2_adjustment = "exact", freq = FALSE)
  ff_se <- function(...) {
    cf <- suppressMessages(coef(
      fit_ncv,
      ci = "pointwise",
      seWithMean = FALSE,
      bias_ref = fit_reml,
      ...
    ))
    cf$smterms[["ff(X)"]]$coef$se
  }
  expect_equal(ff_se(), do.call(ff_se, explicit))
  expect_equal(
    pffr_predict_ci(fit_ncv, bias_ref = fit_reml),
    do.call(pffr_predict_ci, c(list(fit_ncv, bias_ref = fit_reml), explicit))
  )
  # Explicit choices are respected: the fit-time model-based covariance, the
  # CL2 shortcut and the frequentist form each change the SE.
  for (override in list(
    list(sandwich = "none"),
    list(sandwich = "cl2", cl2_adjustment = "shortcut"),
    list(freq = TRUE)
  )) {
    expect_false(isTRUE(all.equal(do.call(ff_se, override), ff_se())))
  }
  L <- predict(fit_ncv, type = "lpmatrix", reformat = FALSE)
  model_based <- pffr_predict_ci(
    fit_ncv,
    bias_ref = fit_reml,
    sandwich = "none"
  )
  expect_equal(
    model_based$se_link,
    unname(sqrt(rowSums((L %*% fit_ncv$Vp) * L)))
  )
  # Without bias_ref, NULL still inherits the fit-time (here model-based) choice.
  expect_equal(pffr_predict_ci(fit_ncv)$se_link, model_based$se_link)
})

test_that("pffr_predict_ci at zero covariates gives the full intercept alpha(t)", {
  skip_if_not_installed("mgcv", "1.9.0")
  set.seed(4219)
  dat <- bias_aware_data()
  fit_ncv <- bias_aware_fit(dat, "NCV")
  fit_reml <- bias_aware_fit(dat, "REML")
  zero <- list(X = I(matrix(0, 1, ncol(dat$X))), z = 0)
  alpha <- pffr_predict_ci(fit_ncv, newdata = zero, bias_ref = fit_reml)
  # The full intercept's rows: the scalar intercept plus the centred
  # Intercept(yindex) basis at the response grid (the study's alpha rows).
  sm <- fit_ncv$smooth[["s(yindex.vec)"]]
  yind <- fit_ncv$pffr$yind
  L <- matrix(0, length(yind), length(fit_ncv$coefficients))
  L[, sm$first.para:sm$last.para] <- mgcv::PredictMat(
    sm,
    data.frame(yindex.vec = yind)
  )
  L[, names(fit_ncv$coefficients) == "(Intercept)"] <- 1
  ref <- bias_aware_direct(L, fit_ncv, fit_reml, bias_aware_cl2(fit_ncv))
  expect_equal(alpha$.index, yind)
  expect_equal(alpha$fit, ref$est, tolerance = 1e-10)
  expect_equal(alpha$se_link, ref$se, tolerance = 1e-10)
  expect_equal(alpha$delta_link, ref$delta, tolerance = 1e-10)
  expect_equal(alpha$lower, ref$lower, tolerance = 1e-10)
  expect_equal(alpha$upper, ref$upper, tolerance = 1e-10)
  # Its estimate is the centred term plus the scalar intercept from coef().
  cf <- suppressMessages(coef(
    fit_ncv,
    bias_ref = fit_reml,
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
  skip_if_not_installed("mgcv", "1.9.0")
  set.seed(4215)
  dat <- bias_aware_data()
  dat$Y[cbind(c(2, 5, 5), c(3, 1, 16))] <- NA
  fit_ncv <- bias_aware_fit(dat, "NCV")
  fit_reml <- bias_aware_fit(dat, "REML")
  ci <- pffr_predict_ci(fit_ncv, bias_ref = fit_reml, sandwich = "cl2")
  n_grid <- ncol(dat$Y)
  expect_equal(nrow(ci), length(dat$Y) - 3)
  observed <- which(!is.na(t(dat$Y)))
  expect_equal(ci$.obs, (observed - 1) %/% n_grid + 1)
  expect_equal(ci$.index, fit_ncv$pffr$yind[(observed - 1) %% n_grid + 1])
  expect_equal(ci$fit, as.vector(fit_ncv$linear.predictors))
  expect_equal(
    ci$delta_link,
    as.vector(fit_ncv$linear.predictors - fit_reml$linear.predictors),
    tolerance = 1e-8
  )
})

test_that("pffr_predict_ci aligns sparse responses, skipping missing values", {
  skip_if_not_installed("mgcv", "1.9.0")
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
  fit_sparse <- function(method) {
    pffr(
      Y ~ z,
      data = dat,
      ydata = yd,
      method = method,
      bs.yindex = list(bs = "ps", k = 5, m = c(2, 1)),
      bs.int = list(bs = "ps", k = 5, m = c(2, 1)),
      sandwich = "none"
    )
  }
  # pffr's NCV does not accept missing sparse responses, so an ML fit stands in
  # for the second fit here: the alignment does not depend on the method.
  fit_ncv <- fit_sparse("REML")
  fit_reml <- fit_sparse("ML")
  # Model-based SEs: pffr_vcov() cannot yet build cluster ids for sparse fits
  # with missing .value rows.
  ci <- pffr_predict_ci(fit_ncv, bias_ref = fit_reml, sandwich = "none")
  kept <- yd[!is.na(yd$.value), ]
  expect_equal(nrow(ci), nrow(kept))
  expect_equal(ci$.obs, kept$.obs)
  expect_equal(ci$.index, kept$.index)
  expect_equal(ci$fit, as.vector(fit_ncv$linear.predictors))
  expect_gt(max(abs(ci$delta_link)), 1e-6)
  expect_equal(
    ci$delta_link,
    as.vector(fit_ncv$linear.predictors - fit_reml$linear.predictors),
    tolerance = 1e-8
  )
  expect_identical(attr(ci, "bias_ref_method"), "ML")
})

test_that("bias references of a different model or data are rejected", {
  skip_if_not_installed("mgcv", "1.9.0")
  set.seed(4216)
  dat <- bias_aware_data()
  fit <- bias_aware_fit(dat, "NCV")
  other <- dat
  other$X[1, ] <- other$X[1, ] + 1
  fit_other <- bias_aware_fit(other, "REML")
  expect_error(coef(fit, bias_ref = fit_other), "same model and data")
  expect_error(
    pffr_predict_ci(fit, bias_ref = fit_other),
    "same model and data"
  )
  fit_k <- pffr(
    Y ~
      ff(
        X,
        splinepars = list(bs = "ps", m = list(c(2, 1), c(2, 1)), k = c(5, 5))
      ) +
        c(z),
    data = dat,
    yind = seq(0, 1, length.out = ncol(dat$Y)),
    method = "REML",
    bs.yindex = list(bs = "ps", k = 6, m = c(2, 1)),
    bs.int = list(bs = "ps", k = 6, m = c(2, 1)),
    knots = list(yindex.vec = seq(-1.2, 2.2, length.out = 10)),
    sandwich = "none"
  )
  expect_error(coef(fit, bias_ref = fit_k), "same model and data")
  expect_error(coef(fit, bias_ref = fit$coefficients), "fitted pffr model")
  expect_error(
    coef(fit, ci = "simultaneous", bias_ref = fit),
    "pointwise only"
  )
  expect_error(
    coef(fit, ci = "pointwise", crit = "tG1", bias_ref = fit),
    "crit = \"z\""
  )
  expect_error(coef(fit, crit = "tG1", bias_ref = fit), "crit = \"z\"")
  expect_error(coef(fit, raw = TRUE, bias_ref = fit), "raw = TRUE")
})

test_that("bias references are rejected for multi-linear-predictor families", {
  skip_if_not_installed("mgcv", "1.9.0")
  fit <- get_gaulss_model()
  expect_error(coef(fit, bias_ref = fit), "single linear-predictor")
  expect_error(pffr_predict_ci(fit, bias_ref = fit), "single linear-predictor")
})
