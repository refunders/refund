test_that("NCV neighbourhoods partition shuffled, unequal blocks with both APIs", {
  ids <- c("b", "a", "b", "c", "a", "b", "d")
  nei <- pffr_ncv_nei(ids)
  expect_identical(nei$a, nei$k)
  expect_identical(nei$ma, nei$m)
  expect_identical(nei$d, nei$i)
  expect_identical(nei$md, nei$mi)
  expect_identical(nei$a, nei$d)
  expect_identical(nei$ma, nei$md)
  expect_false(nei$jackknife)
  expect_equal(sort(nei$a), seq_along(ids))
  groups <- ncv_test_groups(nei)
  expect_equal(unname(sort(lengths(groups))), c(1L, 1L, 2L, 3L))
  for (rows in groups) expect_equal(rows, which(ids == ids[rows[1]]))
  expect_error(pffr_ncv_nei(c(1, NA)), "nonmissing")
  expect_error(pffr_ncv_nei(rep(1, 3)), "at least two blocks")
})

test_that("installed mgcv honours blocks and the cached guard preserves RNG", {
  skip_if_not_installed("mgcv", "1.9.0")
  local_mocked_bindings(.pffr_ncv_cache = new.env(parent = emptyenv()))
  set.seed(719)
  seed <- .Random.seed
  expect_true(pffr_ncv_check_blocks())
  expect_identical(.Random.seed, seed)
  local_mocked_bindings(pffr_ncv_probe = function() stop("must use cache"))
  expect_true(pffr_ncv_check_blocks())
})

test_that("guard fails loudly if mgcv silently ignores blocks", {
  skip_if_not_installed("mgcv", "1.9.0")
  local_mocked_bindings(
    .pffr_ncv_cache = new.env(parent = emptyenv()),
    pffr_ncv_probe = function() list(blocked = c(1, 2), point = c(1, 2))
  )
  expect_error(pffr_ncv_check_blocks(), "ignored the NCV block structure")
  expect_false(isTRUE(.pffr_ncv_cache$checked))
})

test_that("dense NCV aligns retained model-frame rows and curve clusters", {
  skip_if_not_installed("mgcv", "1.9.0")
  set.seed(19)
  dat <- data.frame(
    Y = I(matrix(rnorm(120), 12, 10)),
    z = 1:12,
    subject = rep(c("b", "a", "c", "d"), 3)
  )
  fit <- pffr(
    Y ~ c(z),
    data = dat,
    method = "NCV",
    bs.int = list(k = 4, bs = "ps"),
    sandwich = "none"
  )
  ncv_expect_alignment(fit, fit$model$z)
  expect_identical(fit$method, "NCV")
  dat$Y[cbind(c(1, 4, 8), c(3, 1, 10))] <- NA
  missing <- pffr(
    Y ~ c(z),
    data = dat,
    method = "NCV",
    bs.int = list(k = 4, bs = "ps"),
    sandwich = "none"
  )
  ncv_expect_alignment(missing, missing$model$z)
  expect_equal(nrow(missing$model), 117L)
  grouped <- pffr(
    Y ~ c(z),
    data = dat,
    cluster = subject,
    method = "NCV",
    bs.int = list(k = 4, bs = "ps"),
    sandwich = "none"
  )
  ncv_expect_alignment(grouped, dat$subject[grouped$model$z])
  expect_equal(grouped$pffr$ncv$n_blocks, 4L)
  # Explicit nei uses the same retained-row contract even with missing y.
  expect_message(
    manual <- pffr(
      Y ~ c(z),
      data = dat,
      method = "NCV",
      nei = missing$pffr$ncv$nei,
      bs.int = list(k = 4, bs = "ps"),
      sandwich = "none"
    ),
    "user-supplied"
  )
  expect_equal(manual$coefficients, missing$coefficients, tolerance = 1e-12)
  expect_equal(manual$sp, missing$sp, tolerance = 1e-12)
})

test_that("irregular NCV follows supplied ydata order, not curve sorting", {
  skip_if_not_installed("mgcv", "1.9.0")
  set.seed(21)
  dat <- data.frame(z = 1:12)
  yd <- data.frame(
    .obs = rep(1:12, times = rep(c(5, 7, 9), 4)),
    .index = runif(84),
    .value = rnorm(84)
  )
  yd <- yd[sample(nrow(yd)), ]
  fit <- pffr(
    Y ~ c(z),
    data = dat,
    ydata = yd,
    method = "NCV",
    bs.int = list(k = 4, bs = "ps"),
    sandwich = "none"
  )
  ncv_expect_alignment(fit, fit$model$z)
  expect_equal(as.numeric(fit$model$Y), yd$.value)
})

test_that("NCV rejects unsafe backends and extra row omissions", {
  set.seed(2)
  dat <- data.frame(Y = I(matrix(rnorm(120), 12, 10)), z = 1:12)
  for (backend in c("bam", "gamm", "gamm4", "jagam")) {
    expect_error(
      pffr(
        Y ~ c(z),
        data = dat,
        method = "NCV",
        algorithm = backend,
        sandwich = "none"
      ),
      "supports only algorithm"
    )
  }
  expect_error(
    pffr(
      Y ~ c(z),
      data = dat,
      method = "NCV",
      discrete = TRUE,
      sandwich = "none"
    ),
    "supports only algorithm"
  )
  expect_error(
    pffr(
      Y ~ c(z),
      data = dat,
      method = "NCV",
      discrete = 10,
      sandwich = "none"
    ),
    "supports only algorithm"
  )
  skip_if_not_installed("mgcv", "1.9.0")
  expect_error(
    pffr(
      Y ~ c(z),
      data = dat,
      method = "NCV",
      subset = 1:5,
      bs.int = list(k = 4),
      sandwich = "none"
    ),
    "subset.*not supported"
  )
  dat$z[2] <- NA
  expect_error(
    pffr(
      Y ~ c(z),
      data = dat,
      method = "NCV",
      bs.int = list(k = 4, bs = "ps"),
      sandwich = "none"
    ),
    "cannot align the model frame"
  )
})

test_that("curve NCV smooths dependent ff data more and explicit nei reproduces it", {
  skip_on_cran()
  skip_if_not_installed("mgcv", "1.9.0")
  set.seed(1)
  dat <- ncv_test_data()
  blocked <- ncv_test_fit(dat)
  point <- ncv_test_fit(dat, ncv_blocks = "point")
  expect_true(all(blocked$sp[-1] > point$sp[-1]))
  expect_false(isTRUE(all.equal(blocked$coefficients, point$coefficients)))
  expect_message(
    manual <- ncv_test_fit(
      dat,
      nei = blocked$pffr$ncv$nei,
      ncv_blocks = "point"
    ),
    "user-supplied"
  )
  expect_identical(manual$pffr$ncv$blocks, "user")
  expect_identical(manual$pffr$ncv$nei, blocked$pffr$ncv$nei)
  expect_equal(manual$sp, blocked$sp, tolerance = 1e-12)
  expect_equal(manual$coefficients, blocked$coefficients, tolerance = 1e-12)
})

test_that("Gaussian NCV equals brute-force fixed-penalty curve deletion loss", {
  skip_on_cran()
  skip_if_not_installed("mgcv", "1.9.0")
  set.seed(102)
  dat <- ncv_test_data()
  fit <- ncv_test_fit(dat)
  X <- predict(fit, type = "lpmatrix", reformat = FALSE)
  S <- matrix(0, ncol(X), ncol(X))
  j <- 0L
  for (sm in fit$smooth) {
    ii <- seq.int(sm$first.para, sm$last.para)
    for (penalty in sm$S) {
      j <- j + 1L
      S[ii, ii] <- S[ii, ii] + fit$sp[j] * penalty
    }
  }
  expect_equal(j, length(fit$sp))
  expect_equal(
    as.numeric(solve(crossprod(X) + S, crossprod(X, fit$y))),
    as.numeric(fit$coefficients),
    tolerance = 1e-7
  )
  # Brute-force penalized least-squares refits, preserving the full-data basis,
  # constraints and penalty scaling. Reconstructing smooths after each deletion
  # could change those, and would no longer check the same fixed-sp criterion.
  predictions <- numeric(nrow(X))
  for (rows in ncv_test_groups(fit$pffr$ncv$nei)) {
    beta <- solve(
      crossprod(X[-rows, ]) + S,
      crossprod(X[-rows, ], fit$y[-rows])
    )
    predictions[rows] <- as.numeric(X[rows, ] %*% beta)
  }
  # ?mgcv::NCV defines V as a SUM of losses. gam.fit3's NCV branch uses
  # sum(family$dev.resids(...)) for gamma=1. Gaussian unit-weight deviance
  # is squared error: no division by N, number of curves, scale, or two.
  expect_equal(
    sum((fit$y - predictions)^2),
    as.numeric(fit$gcv.ubre),
    tolerance = 1e-7
  )
  implied <- attr(fit$gcv.ubre, "eta.cv")
  if (!is.null(implied))
    expect_equal(predictions, as.numeric(implied), tolerance = 1e-7)
})

test_that("NCV supports CL2 coefficients and model covariance fallback", {
  skip_on_cran()
  skip_if_not_installed("mgcv", "1.9.0")
  set.seed(102)
  fit <- ncv_test_fit(ncv_test_data())
  robust <- coef(fit, sandwich = "cl2", n1 = 4, n2 = 4)
  model <- coef(fit, sandwich = "none", n1 = 4, n2 = 4)
  rse <- unlist(lapply(robust$smterms, function(x) x$se))
  mse <- unlist(lapply(model$smterms, function(x) x$se))
  expect_true(length(rse) > 0L && all(is.finite(rse)))
  expect_gt(max(abs(rse - mse)), 1e-4)
  fit$Vc <- NULL
  expect_equal(pffr_vcov(fit, sandwich = "none"), fit$Vp)
  fit$edf2 <- NULL
  expect_warning(
    v <- pffr_vcov(
      fit,
      sandwich = "cluster",
      dof_correction = "edf",
      edf_type = "edf2"
    ),
    "edf2 is unavailable for NCV"
  )
  expect_true(all(is.finite(v)))
})


test_that("NCV rejects NA padding and missing weights and supports fit-time CL2", {
  skip_if_not_installed("mgcv", "1.9.0")
  set.seed(51)
  dat <- data.frame(Y = I(matrix(rnorm(120), 12, 10)), z = 1:12)
  dat$Y[2, 3] <- NA
  expect_error(
    pffr(
      Y ~ c(z),
      data = dat,
      method = "NCV",
      bs.int = list(bs = "ps", k = 4),
      sandwich = "none",
      na.action = na.exclude
    ),
    "requires na.action = na.omit"
  )
  weights <- rep(1, 12)
  weights[5] <- NA
  expect_error(
    pffr(
      Y ~ c(z),
      data = dat,
      method = "NCV",
      bs.int = list(bs = "ps", k = 4),
      sandwich = "none",
      weights = weights
    ),
    "cannot align the model frame"
  )
  expect_warning(
    fit <- pffr(
      Y ~ c(z),
      data = dat,
      method = "NCV",
      bs.int = list(bs = "ps", k = 4),
      sandwich = "cl2"
    ),
    class = "pffr_small_G_warning"
  )
  expect_true(all(is.finite(fit$pffr$Vsandwich)))
  fixed <- pffr(
    Y ~ c(z),
    data = dat,
    method = "NCV",
    bs.int = list(bs = "ps", k = 4),
    sandwich = "none",
    sp = fit$sp
  )
  expect_equal(fixed$coefficients, fit$coefficients, tolerance = 1e-7)
})
