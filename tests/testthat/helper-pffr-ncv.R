ncv_test_data <- function() {
  n <- 24L
  tt <- seq(0, 1, length.out = 16)
  X <- matrix(rnorm(n * length(tt)), n)
  X <- X + outer(rnorm(n), sin(2 * pi * tt))
  E <- t(replicate(n, as.numeric(arima.sim(list(ar = .95), n = length(tt)))))
  Y <- outer(rowMeans(X), cos(pi * tt)) + E
  data.frame(Y = I(Y), X = I(X), z = seq_len(n))
}

ncv_test_fit <- function(dat, ...) {
  pffr(
    Y ~
      ff(
        X,
        splinepars = list(bs = "ps", m = list(c(2, 1), c(2, 1)), k = c(4, 4))
      ),
    data = dat,
    yind = seq(0, 1, length.out = ncol(dat$Y)),
    method = "NCV",
    bs.int = list(bs = "ps", k = 4),
    sandwich = "none",
    ...
  )
}

# Response on the scale of `family`, driven by the dependent Gaussian data.
ncv_test_glm_data <- function(family) {
  dat <- ncv_test_data()
  eta <- 0.5 * scale(dat$Y)
  y <- if (family$family == "poisson") {
    rpois(length(eta), 3 * exp(eta))
  } else {
    rbinom(length(eta), 1, plogis(eta - 0.5))
  }
  dat$Y <- I(matrix(y, nrow(eta)))
  dat
}

# Total penalty matrix sum_j sp_j S_j of a fitted model, in coefficient order.
ncv_test_penalty <- function(fit) {
  p <- length(fit$coefficients)
  S <- matrix(0, p, p)
  j <- 0L
  for (sm in fit$smooth) {
    ii <- seq.int(sm$first.para, sm$last.para)
    for (penalty in sm$S) {
      j <- j + 1L
      S[ii, ii] <- S[ii, ii] + fit$sp[j] * penalty
    }
  }
  stopifnot(j == length(fit$sp))
  S
}

# Penalized IRLS at fixed penalty S (minimizes deviance + beta' S beta).
ncv_test_pirls <- function(X, y, S, family, beta, tol = 1e-12) {
  for (iter in 1:100) {
    eta <- as.numeric(X %*% beta)
    mu <- family$linkinv(eta)
    deta <- family$mu.eta(eta)
    w <- deta^2 / family$variance(mu)
    z <- eta + (y - mu) / deta
    new <- solve(crossprod(X, w * X) + S, crossprod(X, w * z))
    if (max(abs(new - beta)) < tol) break
    beta <- new
  }
  as.numeric(new)
}

ncv_test_groups <- function(nei) {
  split(nei$a, rep(seq_along(nei$ma), diff(c(0L, nei$ma))))
}

ncv_expect_alignment <- function(fit, ids) {
  # ids are independently recovered from actual fitted model-frame values.
  groups <- ncv_test_groups(fit$pffr$ncv$nei)
  expect_equal(
    sort(unlist(groups, use.names = FALSE)),
    seq_len(nrow(fit$model))
  )
  expect_length(groups, length(unique(ids)))
  for (rows in groups) {
    expect_length(unique(ids[rows]), 1L)
    expect_equal(sort(rows), which(ids == ids[rows[1]]))
  }
}
