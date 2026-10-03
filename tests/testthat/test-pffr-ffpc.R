# Tests for ffpc() terms: implied coefficient surface (ffpcplot, coef.pffr),
# its standard errors, and FPC scores of new data in predict.pffr.

# Covariate curves of rank 5 on [0, xmax] plus white noise, and a response
# from int X(s) beta(s,t) ds with trapezoidal weights.
sim_ffpc_data <- function(
  n = 60,
  S = 40,
  nt = 30,
  xmax = 10,
  noise = 0.01
) {
  s <- seq(0, xmax, length = S)
  t <- seq(0, 1, length = nt)
  basis <- cbind(1, stats::poly(s, degree = 4))
  X <- sapply(5:1, \(l) stats::rnorm(n, sd = sqrt(l))) %*%
    t(basis) +
    noise * matrix(stats::rnorm(n * S), n, S)
  beta <- outer(s, t, \(s, t) cos(2 * pi * s / xmax * t)) / xmax
  w <- quadWeights(s, method = "trapezoidal")
  y <- X %*% (w * beta) + 0.05 * matrix(stats::rnorm(n * nt), n, nt)
  list(data = list(y = y, X = X), s = s, t = t, beta = beta, w = w)
}

test_that("ffpcplot works and its surface equals the coef.pffr surface", {
  skip_on_cran()
  set.seed(4401)
  sim <- sim_ffpc_data()
  s <- sim$s
  m <- pffr(
    y ~ c(1) + 0 + ffpc(X, xind = s, decomppars = list(npc = 5)),
    data = sim$data,
    yind = sim$t,
    sandwich = "none"
  )

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_no_error(ffpcplot(m))
  fp <- ffpcplot(m, plot = FALSE)
  expect_equal(dim(fp$phibeta[[1]]), c(length(s), length(sim$t)))
  expect_equal(dim(fp$betatilde), c(length(sim$t), 5))

  cf <- coef(m, n2 = length(sim$t))
  expect_named(cf$smterms, "ffpc(X)")
  surf <- cf$smterms[["ffpc(X)"]]
  expect_equal(surf$dim, 2)
  expect_equal(surf$x, s)
  expect_equal(surf$y, sim$t)
  expect_equal(surf$coef[, 1], rep(s, length(sim$t)))
  expect_equal(surf$coef[, 2], rep(sim$t, each = length(s)))
  expect_equal(
    surf$coef[, "value"],
    as.vector(fp$phibeta[[1]]),
    tolerance = 1e-10
  )
  # and it estimates the true surface, not a multiple of it
  rel_rmse <- sqrt(mean((fp$phibeta[[1]] - sim$beta)^2) / mean(sim$beta^2))
  expect_lt(rel_rmse, 0.5)
})

test_that("the integrated ffpc surface reproduces the fitted ffpc term", {
  skip_on_cran()
  set.seed(4402)
  sim <- sim_ffpc_data()
  s <- sim$s
  m <- pffr(
    y ~ c(1) + 0 + ffpc(X, xind = s, decomppars = list(npc = 5)),
    data = sim$data,
    yind = sim$t,
    sandwich = "none"
  )
  trm <- m$pffr$ffpc[[1]]
  fp <- ffpcplot(m, plot = FALSE)
  beta_hat <- fp$phibeta[[1]]

  # the fitted term sum_k xi_k b_k(t), with the scores used in the fit
  # (predict() returns one term per FPC)
  xi <- ffpc(sim$data$X, xind = s, decomppars = list(npc = 5))$data
  fitted_term <- xi %*% t(fp$betatilde)
  expect_equal(
    as.vector(fitted_term),
    as.vector(Reduce(`+`, lapply(predict(m, type = "terms"), unclass))),
    tolerance = 1e-8
  )

  # exactly, for the FPCA reconstruction of the centred curves ...
  X_rec <- xi %*% t(trm$PCMat)
  expect_equal(X_rec %*% (sim$w * beta_hat), fitted_term, tolerance = 1e-10)
  # ... and up to score shrinkage and truncation, for the centred curves
  X_c <- sweep(sim$data$X, 2, trm$meanX)
  expect_equal(X_c %*% (sim$w * beta_hat), fitted_term, tolerance = 0.01)
})

test_that("the ffpc surface does not depend on the units of xind", {
  skip_on_cran()
  set.seed(4403)
  sim <- sim_ffpc_data(xmax = 1)
  s1 <- sim$s
  s10 <- 10 * sim$s
  m1 <- pffr(
    y ~ ffpc(X, xind = s1, decomppars = list(npc = 5)),
    data = sim$data,
    yind = sim$t,
    sandwich = "none"
  )
  m10 <- pffr(
    y ~ ffpc(X, xind = s10, decomppars = list(npc = 5)),
    data = sim$data,
    yind = sim$t,
    sandwich = "none"
  )
  expect_equal(fitted(m1), fitted(m10), tolerance = 1e-5)
  b1 <- coef(m1, se = FALSE)$smterms[["ffpc(X)"]]$coef[, "value"]
  b10 <- coef(m10, se = FALSE)$smterms[["ffpc(X)"]]$coef[, "value"]
  # beta(s, t) ds is invariant: the surface on [0, 10] is 1/10 of that on [0, 1]
  expect_equal(b10, b1 / 10, tolerance = 1e-4)

  # terms from older versions (FPCA on seq(0, 1)) are rescaled to xind
  expect_equal(ffpc_beta_scale(list(xind = s10)), 1 / 10)
  expect_equal(ffpc_beta_scale(m10$pffr$ffpc[[1]]), 1)
})

test_that("coef.pffr SEs of the ffpc surface are sqrt(diag(L V L'))", {
  skip_on_cran()
  set.seed(4404)
  sim <- sim_ffpc_data(n = 50)
  s <- sim$s
  m <- pffr(
    y ~ ffpc(X, xind = s, decomppars = list(npc = 4)),
    data = sim$data,
    yind = sim$t,
    sandwich = "none"
  )
  trm <- m$pffr$ffpc[[1]]
  n2 <- 12
  tg <- seq(min(sim$t), max(sim$t), length = n2)
  ns <- length(s)

  # L built independently: row (s_i, t_j) = sum_k psi_k(s_i) B_k(t_j)
  L <- matrix(0, ns * n2, length(m$coefficients))
  for (k in 1:4) {
    sm <- m$smooth[[which(vapply(
      m$smooth,
      \(sm) sm$by == paste0("X.PC", k),
      logical(1)
    ))]]
    nd <- data.frame(tg, 1)
    names(nd) <- c(sm$term, sm$by)
    B <- mgcv::PredictMat(sm, nd)
    cols <- sm$first.para:sm$last.para
    L[, cols] <- B[rep(seq_len(n2), each = ns), ] *
      trm$PCMat[rep(seq_len(ns), times = n2), k]
  }

  for (sw in c("none", "cl2")) {
    V <- pffr_vcov(m, sandwich = sw, freq = FALSE)
    cf <- coef(m, sandwich = sw, n2 = n2, ci = "pointwise")
    surf <- cf$smterms[["ffpc(X)"]]
    expect_equal(surf$coef[, "value"], drop(L %*% m$coefficients))
    expect_equal(surf$coef[, "se"], sqrt(rowSums((L %*% V) * L)))
    expect_equal(
      surf$coef[, "upper"] - surf$coef[, "value"],
      stats::qnorm(0.975) * surf$coef[, "se"]
    )
  }

  cf_satt <- coef(
    m,
    sandwich = "cl2",
    n2 = n2,
    ci = "pointwise",
    crit = "satterthwaite"
  )
  df_satt <- cf_satt$smterms[["ffpc(X)"]]$coef[, "df"]
  expect_true(all(is.finite(df_satt) & df_satt > 0 & df_satt < 50))

  cf_sim <- coef(m, n2 = n2, ci = "simultaneous", n_sim = 200, sim_seed = 1)
  surf <- cf_sim$smterms[["ffpc(X)"]]
  expect_gt(surf$crit, stats::qnorm(0.975))
  expect_true(all(surf$coef[, "lower"] <= surf$coef[, "upper"]))
})

test_that("predict.pffr uses the fitted FPC scores for new curves", {
  skip_on_cran()
  set.seed(4405)
  # noisy curves (score shrinkage matters) and npc.max below the FPCA's npc
  sim <- sim_ffpc_data(n = 50, noise = 1)
  m <- pffr(
    y ~ ffpc(X, decomppars = list(npc = 5), npc.max = 3),
    data = sim$data,
    yind = sim$t,
    sandwich = "none"
  )
  trm <- m$pffr$ffpc[[1]]
  expect_gt(ncol(trm$score_pars$efunctions), 3)
  expect_equal(
    ffpc_scores(trm, sim$data$X),
    unname(ffpc(sim$data$X, decomppars = list(npc = 5), npc.max = 3)$data),
    tolerance = 1e-10
  )
  expect_equal(
    predict(m, newdata = sim$data),
    fitted(m),
    tolerance = 1e-10,
    ignore_attr = TRUE
  )
  # curves with missing values get BLUPs over their observed points
  X_na <- sim$data$X[1:5, ]
  X_na[, 1:5] <- NA
  sc <- ffpc_scores(trm, X_na)
  expect_equal(dim(sc), c(5, 3))
  expect_false(anyNA(sc))
})
