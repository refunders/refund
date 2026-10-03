# The default intervals of pffr() fits: curve-clustered CL2 sandwich with the
# full-block ("exact") leverage adjustment in Bayesian form, and Satterthwaite
# critical values, for every estimand type. The reference below is an
# independent dense implementation of
#   A_g = {(I - H)^2}_gg^(-1/2),  u_g = Xw_g' A_g z_g,
#   V = G/(G-1) Vp (sum_g u_g u_g') Vp + (Vp - Ve),
# and of the moment df tr(Gamma)^2 / tr(Gamma^2) with
#   Gamma = P'P,  P[, g] = (I - H)[g, ]' A_g Xw_g Vp a.

defaults_data <- function(n = 40, nt = 12, ns = 15) {
  set.seed(6021)
  t <- seq(0, 1, length.out = nt)
  s <- seq(0, 1, length.out = ns)
  X1 <- matrix(rnorm(n * ns), n) + outer(rnorm(n), sin(2 * pi * s))
  # low-rank covariate for ffpc(): three components
  X2 <- outer(rnorm(n, sd = 2), sin(pi * s)) +
    outer(rnorm(n), cos(pi * s)) +
    outer(rnorm(n, sd = 0.5), sin(2 * pi * s))
  zs <- runif(n)
  zl <- rnorm(n)
  zc <- rnorm(n)
  L1 <- X1 %*% outer(s, t, \(s, t) s * t) / ns
  L2 <- X2 %*% outer(s, t, \(s, t) cos(pi * s) * t) / ns
  E <- t(replicate(n, as.numeric(arima.sim(list(ar = 0.7), n = nt))))
  Y <- L1 + L2 + outer(sin(2 * pi * zs), t) + outer(zl, t) + 0.5 * zc + E
  list(
    data = data.frame(
      Y = I(Y),
      X1 = I(X1),
      X2 = I(X2),
      zs = zs,
      zl = zl,
      zc = zc
    ),
    t = t,
    s = s
  )
}

defaults_fit <- function(d) {
  pffr(
    Y ~
      ff(
        X1,
        xind = d$s,
        splinepars = list(bs = "ps", m = list(c(2, 1), c(2, 1)), k = c(5, 5))
      ) +
        ffpc(X2, xind = d$s, splinepars = list(bs = "ps", m = c(2, 1), k = 5)) +
        s(zs, k = 5) +
        zl +
        c(zc),
    yind = d$t,
    data = d$data,
    bs.yindex = list(bs = "ps", k = 6, m = c(2, 1)),
    bs.int = list(bs = "ps", k = 6, m = c(2, 1))
  )
}

# Dense exact-CL2 reference for single-linear-predictor exponential families.
dense_cl2 <- function(fit, cap = 0.999) {
  b <- fit
  class(b) <- setdiff(class(b), "pffr")
  X <- model.matrix(b)
  mu <- b$fitted.values
  mu_eta <- b$family$mu.eta(b$linear.predictors)
  pw <- b$prior.weights
  common <- pw / (b$sig2 * b$family$variance(mu))
  Xw <- X * (mu_eta * sqrt(common))
  z <- (b$y - mu) * sqrt(common)
  cid <- rep(seq_len(fit$pffr$nobs), each = fit$pffr$nyindex)
  Vp <- b$Vp
  M <- diag(nrow(X)) - Xw %*% Vp %*% t(Xw)
  ids <- unique(cid)
  G <- length(ids)
  A <- U <- vector("list", G)
  meat <- matrix(0, ncol(X), ncol(X))
  for (g in seq_len(G)) {
    ii <- which(cid == ids[g])
    Bg <- tcrossprod(M[ii, , drop = FALSE])
    ee <- eigen((Bg + t(Bg)) / 2, symmetric = TRUE)
    A[[g]] <- ee$vectors %*%
      diag(1 / sqrt(pmax(ee$values, (1 - cap)^2))) %*%
      t(ee$vectors)
    u <- crossprod(Xw[ii, , drop = FALSE], A[[g]] %*% z[ii])
    meat <- meat + tcrossprod(u)
  }
  V <- G / (G - 1) * Vp %*% meat %*% Vp + (Vp - b$Ve)
  df <- function(Lrows) {
    apply(Lrows, 1, \(a) {
      P <- vapply(
        seq_len(G),
        \(g) {
          ii <- which(cid == ids[g])
          q <- A[[g]] %*% Xw[ii, , drop = FALSE] %*% Vp %*% a
          as.vector(crossprod(M[ii, , drop = FALSE], q))
        },
        numeric(nrow(X))
      )
      Gamma <- crossprod(P)
      sum(diag(Gamma))^2 / sum(Gamma^2)
    })
  }
  list(V = V, df = df)
}

# Rows of a smooth term on coef()'s grid, with mgcv's seWithMean convention.
term_rows <- function(fit, label, grid) {
  sm <- fit$smooth[[label]]
  X <- mgcv::PredictMat(sm, grid)
  p <- length(fit$coefficients)
  ii <- sm$first.para:sm$last.para
  if (attr(sm, "nCons") > 0) {
    L <- matrix(fit$cmX, nrow(X), p, byrow = TRUE)
    if (!is.null(sm$meanL1)) L <- L / sm$meanL1
  } else {
    L <- matrix(0, nrow(X), p)
  }
  L[, ii] <- X
  L
}

expect_cl2_satt <- function(coef_table, L, ref, tolerance = 1e-7) {
  se <- sqrt(rowSums((L %*% ref$V) * L))
  df <- ref$df(L)
  expect_equal(unname(coef_table[, "se"]), unname(se), tolerance = tolerance)
  expect_equal(unname(coef_table[, "df"]), unname(df), tolerance = 1e-6)
  expect_true(all(df < Inf))
  expect_equal(
    unname(coef_table[, "upper"] - coef_table[, "value"]),
    unname(qt(0.975, df) * se),
    tolerance = tolerance
  )
}

test_that("default intervals are exact CL2 + Satterthwaite for every estimand", {
  skip_on_cran()
  d <- defaults_data()
  fit <- defaults_fit(d)
  expect_identical(fit$pffr$sandwich, "cl2")
  expect_identical(fit$pffr$sandwich_info$cl2_adjustment, "exact")
  ref <- dense_cl2(fit)
  expect_equal(
    unname(matrix(vcov(fit), nrow(ref$V))),
    unname(ref$V),
    tolerance = 1e-7
  )

  cf <- suppressMessages(coef(fit, ci = "pointwise", n1 = 15, n2 = 8))
  expect_identical(cf$ci_meta$crit_used, "satterthwaite")
  expect_true(cf$ci_meta$sandwich)
  # scalar effect c(zc) and the scalar intercept: parametric coefficients
  pnames <- rownames(cf$pterms)
  Lp <- diag(length(fit$coefficients))[
    match(pnames, names(fit$coefficients)),
    ,
    drop = FALSE
  ]
  expect_setequal(pnames, c("(Intercept)", "zc"))
  expect_cl2_satt(cf$pterms, Lp, ref)

  # functional intercept (centred, seWithMean = TRUE), coefficient function of
  # a scalar covariate zl(t), smooth effect f(zs, t), ff surface
  for (term in c("Intercept(yindex)", "zl(yindex)", "s(zs", "ff(X1")) {
    i <- grep(term, names(cf$smterms), fixed = TRUE)
    expect_length(i, 1)
    tab <- cf$smterms[[i]]$coef
    label <- names(fit$pffr$short_labels)[
      fit$pffr$short_labels == names(cf$smterms)[i]
    ]
    expect_length(label, 1)
    L <- term_rows(fit, label, tab)
    expect_cl2_satt(
      as.matrix(tab[, c("value", "se", "lower", "upper", "df")]),
      L,
      ref
    )
  }

  # ffpc surface: a fixed linear map of the coefficients
  ffpc_tab <- cf$smterms[[grep("ffpc", names(cf$smterms))]]$coef
  map <- ffpc_surface_map(fit, 1, seq(min(d$t), max(d$t), length.out = 8))
  expect_cl2_satt(
    as.matrix(ffpc_tab[, c("value", "se", "lower", "upper", "df")]),
    map$L,
    ref
  )

  # fitted means and new predictions
  nd <- list(
    X1 = unclass(d$data$X1)[1:3, ],
    X2 = unclass(d$data$X2)[1:3, ],
    zs = d$data$zs[1:3],
    zl = d$data$zl[1:3],
    zc = d$data$zc[1:3]
  )
  Lnew <- predict(fit, newdata = nd, type = "lpmatrix", reformat = FALSE)
  p <- predict(fit, newdata = nd, se.fit = TRUE, reformat = FALSE)
  expect_equal(
    as.vector(p$se.fit),
    unname(sqrt(rowSums((Lnew %*% ref$V) * Lnew))),
    tolerance = 1e-7
  )
  expect_equal(as.vector(p$df), unname(ref$df(Lnew)), tolerance = 1e-6)
  expect_equal(as.vector(p$crit), qt(0.975, as.vector(p$df)))
  ci <- pffr_predict_ci(fit, newdata = nd)
  expect_equal(ci$df, as.vector(p$df))
  expect_equal(ci$upper - ci$fit, as.vector(p$crit * p$se.fit))
})

test_that("sandwich = FALSE gives model-based intervals with Gaussian values", {
  skip_on_cran()
  d <- defaults_data()
  fit <- pffr(
    Y ~ zl + c(zc),
    yind = d$t,
    data = d$data,
    bs.yindex = list(bs = "ps", k = 6, m = c(2, 1)),
    bs.int = list(bs = "ps", k = 6, m = c(2, 1)),
    sandwich = FALSE
  )
  expect_identical(fit$pffr$sandwich, "none")
  expect_null(fit$pffr[["Vsandwich"]])
  expect_identical(vcov(fit), fit$Vc)
  cf <- coef(fit, ci = "pointwise", n1 = 10)
  expect_identical(cf$ci_meta$crit_used, "z")
  expect_equal(
    unname(cf$pterms[, "se"]),
    unname(sqrt(diag(fit$Vc))[c(1, which(names(fit$coefficients) == "zc"))])
  )
  # the CL2 sandwich can still be requested afterwards
  cf2 <- suppressMessages(coef(fit, ci = "pointwise", sandwich = TRUE, n1 = 10))
  expect_identical(cf2$ci_meta$crit_used, "satterthwaite")
})

test_that("deprecated character values of sandwich map to TRUE/FALSE", {
  skip_on_cran()
  d <- defaults_data()
  f <- function(sandwich) {
    pffr(
      Y ~ c(zc),
      yind = d$t,
      data = d$data,
      bs.yindex = list(bs = "ps", k = 5, m = c(2, 1)),
      bs.int = list(bs = "ps", k = 5, m = c(2, 1)),
      sandwich = sandwich
    )
  }
  expect_warning(m <- f("none"), "deprecated; use sandwich = FALSE")
  expect_identical(m$pffr$sandwich, "none")
  expect_warning(m <- f("cl2"), "deprecated; use sandwich = TRUE")
  expect_identical(m$pffr$sandwich, "cl2")
  for (old in c("cluster", "hc")) {
    expect_warning(m <- f(old), "no longer available")
    expect_identical(m$pffr$sandwich_info$cl2_adjustment, "exact")
  }
  expect_error(f("auto"), "must be TRUE")
  expect_warning(
    coef(m, sandwich = "none", n1 = 5),
    "deprecated; use sandwich = FALSE"
  )
  expect_warning(coef(m, freq = TRUE, n1 = 5), "`freq` is deprecated")
})

test_that("binary responses get a one-time message about intercepts", {
  skip_on_cran()
  d <- defaults_data()
  set.seed(6022)
  d$data$Y <- I(matrix(
    rbinom(length(d$data$Y), 1, plogis(scale(d$data$Y))),
    nrow(d$data$Y)
  ))
  rm(list = intersect("binary_cl2", ls(.pffr_state)), envir = .pffr_state)
  f <- function() {
    pffr(
      Y ~ zl,
      yind = d$t,
      data = d$data,
      family = binomial(),
      bs.yindex = list(bs = "ps", k = 5, m = c(2, 1)),
      bs.int = list(bs = "ps", k = 5, m = c(2, 1))
    )
  }
  expect_message(f(), "Binary response: CL2 intervals for the functional")
  expect_no_message(f())
})

test_that("fits without a cluster score path fall back to model-based", {
  skip_on_cran()
  d <- defaults_data()
  expect_message(
    m <- pffr(
      Y ~ c(zc),
      yind = d$t,
      data = d$data,
      algorithm = "gamm",
      bs.yindex = list(bs = "ps", k = 5, m = c(2, 1)),
      bs.int = list(bs = "ps", k = 5, m = c(2, 1))
    ),
    "model-based intervals \\(sandwich = FALSE\\) because the sandwich is not available"
  )
  expect_identical(m$gam$pffr$sandwich, "none")
  expect_warning(
    pffr(
      Y ~ c(zc),
      yind = d$t,
      data = d$data,
      algorithm = "gamm",
      bs.yindex = list(bs = "ps", k = 5, m = c(2, 1)),
      bs.int = list(bs = "ps", k = 5, m = c(2, 1)),
      sandwich = TRUE
    ),
    "not available for algorithm"
  )
  expect_match(
    pffr_cl2_unavailable(list(family = mgcv::multinom(K = 2)), "gam", 1:3),
    "no cluster-robust covariance"
  )
  expect_match(
    pffr_cl2_unavailable(list(family = gaussian()), "gam", rep(1, 3)),
    "at least two"
  )
  expect_null(pffr_cl2_unavailable(list(family = gaussian()), "gam", 1:3))
})
