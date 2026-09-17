# Independent dense reference: full N-dimensional residualization, no compression.
dense_reference <- function(
  Z,
  B,
  cid,
  z,
  contrasts,
  adjustment = "exact",
  cap = .999
) {
  H <- Z %*% B %*% t(Z)
  L <- diag(nrow(Z)) - H
  ids <- unique(cid)
  K <- matrix(0, ncol(Z), length(ids))
  qs <- vector("list", length(ids))
  for (g in seq_along(ids)) {
    ii <- which(cid == ids[g])
    Zg <- Z[ii, , drop = FALSE]
    Rg <- if (adjustment == "exact") tcrossprod(L[ii, , drop = FALSE]) else
      diag(length(ii)) - H[ii, ii, drop = FALSE]
    ee <- eigen((Rg + t(Rg)) / 2, symmetric = TRUE)
    floor <- if (adjustment == "exact") (1 - cap)^2 else 1 - cap
    A <- if (adjustment == "none") diag(length(ii)) else
      tcrossprod(
        sweep(ee$vectors, 2, sqrt(pmax(ee$values, floor)), "/"),
        ee$vectors
      )
    K[, g] <- B %*% crossprod(Zg, A %*% z[ii])
    qs[[g]] <- A %*% Zg %*% B %*% t(contrasts)
  }
  df <- vapply(
    seq_len(nrow(contrasts)),
    function(j) {
      P <- matrix(0, nrow(Z), length(ids))
      for (g in seq_along(ids))
        P[, g] <- t(L[which(cid == ids[g]), , drop = FALSE]) %*% qs[[g]][, j]
      Gamma <- crossprod(P)
      sum(diag(Gamma))^2 / sum(Gamma^2)
    },
    0.
  )
  list(K = K, df = df)
}

testthat::test_that("compressed scores and full residualization moments match dense reference", {
  set.seed(84001)
  G <- 7L
  p <- 12L
  D <- 19L
  r <- 4L
  Z <- do.call(
    rbind,
    lapply(
      seq_len(G),
      function(g) matrix(rnorm(D * r), D) %*% matrix(rnorm(r * p), r)
    )
  )
  cid <- rep(seq_len(G), each = D)
  B <- solve(crossprod(Z) + diag(seq_len(p)))
  z <- rnorm(nrow(Z))
  Xp <- matrix(rnorm(5 * p), 5)
  for (adjustment in c("exact", "shortcut", "none")) {
    core <- pffr_influence_core(Z, B, cid, z, adjustment)
    ref <- dense_reference(Z, B, cid, z, Xp, adjustment)
    testthat::expect_equal(core$K, ref$K, tolerance = 1e-11)
    testthat::expect_equal(
      pffr_influence_df(core, Xp)$df,
      ref$df,
      tolerance = 1e-10
    )
    testthat::expect_equal(
      pffr_influence_df(core, Xp, 1)$df,
      pffr_influence_df(core, Xp)$df
    )
    testthat::expect_equal(core$diagnostics$rank, rep(r, G))
    ord <- sample(nrow(Z))
    other <- pffr_influence_core(Z[ord, ], B, cid[ord], z[ord], adjustment)
    testthat::expect_equal(
      other$K[, match(core$groups, other$groups)],
      core$K,
      tolerance = 1e-10
    )
    testthat::expect_equal(
      pffr_influence_df(other, Xp)$df,
      ref$df,
      tolerance = 1e-10
    )
  }
})

testthat::test_that("intercept-only OLS recovers G minus one, not G", {
  for (G in c(3L, 7L, 20L)) {
    core <- pffr_influence_core(matrix(1, G, 1), matrix(1 / G), seq_len(G))
    testthat::expect_equal(
      pffr_influence_df(core, matrix(1))$df,
      G - 1,
      tolerance = 1e-10
    )
    testthat::expect_true(is.na(pffr_influence_df(core, matrix(0))$df))
  }
})

testthat::test_that("OLS moment df agrees with clubSandwich CR2", {
  testthat::skip_if_not_installed("clubSandwich")
  set.seed(84002)
  cid <- rep(1:13, times = 2:14)
  n <- length(cid)
  dat <- data.frame(y = rnorm(n), x = rnorm(n), v = rnorm(n))
  fit <- lm(y ~ x + v, data = dat)
  Z <- model.matrix(fit)
  core <- pffr_influence_core(Z, solve(crossprod(Z)), cid, residuals(fit))
  ref <- clubSandwich::coef_test(
    fit,
    vcov = "CR2",
    cluster = cid,
    test = "Satterthwaite"
  )
  testthat::expect_equal(
    pffr_influence_df(core, diag(ncol(Z)))$df,
    ref$df_Satt,
    tolerance = 1e-9
  )
})

testthat::test_that("zero, full rank and saturated blocks are explicit", {
  Z <- rbind(matrix(0, 3, 3), diag(3))
  cid <- rep(1:2, each = 3)
  B <- diag(3)
  core <- pffr_influence_core(Z, B, cid, 1:6)
  testthat::expect_equal(core$diagnostics$rank, c(0L, 3L))
  testthat::expect_equal(core$K[, 1], rep(0, 3))
  testthat::expect_equal(core$diagnostics$n_floored, c(0L, 3L))
  testthat::expect_equal(
    core$diagnostics$min_block_eig,
    c(1, 0),
    tolerance = 1e-14
  )
  testthat::expect_true(all(is.na(pffr_influence_df(core, diag(3))$df)))
  testthat::expect_error(pffr_influence_core(Z, B, rep(1, 6)), "at least two")
  testthat::expect_error(pffr_influence_core(Z, B, c(NA, cid[-1])), "missing")
  testthat::expect_error(
    pffr_influence_core(Z, B, cid, rep(Inf, 6)),
    "finite residual"
  )
  Z <- rbind(diag(3), diag(3))
  B <- diag(.4, 3)
  cid <- rep(1:2, each = 3)
  z <- 1:6
  core <- pffr_influence_core(Z, B, cid, z)
  ref <- dense_reference(Z, B, cid, z, diag(3))
  testthat::expect_equal(core$K, ref$K, tolerance = 1e-10)
  testthat::expect_equal(
    pffr_influence_df(core, diag(3))$df,
    ref$df,
    tolerance = 1e-10
  )
})

testthat::test_that("sampling, B2, finite factor and centering stay separate", {
  set.seed(84003)
  Z <- matrix(rnorm(120), 30)
  B <- solve(crossprod(Z) + diag(4))
  core <- pffr_influence_core(Z, B, rep(1:10, each = 3), rnorm(30))
  core$B2 <- B %*% B
  plain <- function(V) matrix(V, nrow(V))
  testthat::expect_equal(
    plain(pffr_influence_vcov(core)),
    core$correction * tcrossprod(core$K) + core$B2
  )
  testthat::expect_equal(
    plain(pffr_influence_vcov(core, freq = TRUE)),
    core$correction * tcrossprod(core$K)
  )
  testthat::expect_equal(
    pffr_influence_vcov(core, b2 = FALSE),
    pffr_influence_vcov(core, freq = TRUE)
  )
  testthat::expect_equal(
    plain(pffr_influence_vcov(core, center_scores = TRUE)),
    core$correction * tcrossprod(core$K - rowMeans(core$K)) + core$B2
  )
})

# The diagonal moment df of THIS influence object's geometry: same formula as
# the historical shortcut, but evaluated with whatever leverage weight the
# covariance resolved to. Used to check the df_gram = "diagonal" branch against
# an independent dense implementation for every adjustment; for backwards
# compatibility with the pre-2026-09 numbers see historical_satterthwaite_df()
# (helper-historical-df.R), which always uses the shortcut weight.
diagonal_df_reference <- function(
  Z,
  B,
  cid,
  contrasts,
  adjustment = "exact",
  cap = .999
) {
  H <- Z %*% B %*% t(Z)
  L <- diag(nrow(Z)) - H
  ids <- unique(cid)
  s2 <- s4 <- numeric(nrow(contrasts))
  for (g in seq_along(ids)) {
    ii <- which(cid == ids[g])
    Zg <- Z[ii, , drop = FALSE]
    Rg <- if (adjustment == "exact") tcrossprod(L[ii, , drop = FALSE]) else
      diag(length(ii)) - H[ii, ii, drop = FALSE]
    ee <- eigen((Rg + t(Rg)) / 2, symmetric = TRUE)
    floor <- if (adjustment == "exact") (1 - cap)^2 else 1 - cap
    A <- if (adjustment == "none") diag(length(ii)) else
      tcrossprod(
        sweep(ee$vectors, 2, sqrt(pmax(ee$values, floor)), "/"),
        ee$vectors
      )
    cn2 <- colSums((A %*% Zg %*% B %*% t(contrasts))^2)
    s2 <- s2 + cn2
    s4 <- s4 + cn2^2
  }
  df <- s2^2 / s4
  df[!is.finite(df)] <- NA_real_
  ok <- is.finite(df)
  df[ok] <- pmin(pmax(df[ok], 1), length(ids))
  df
}

testthat::test_that("df_gram = 'diagonal' matches the same-adjustment diagonal df", {
  set.seed(84001)
  G <- 7L
  p <- 12L
  D <- 19L
  r <- 4L
  Z <- do.call(
    rbind,
    lapply(
      seq_len(G),
      function(g) matrix(rnorm(D * r), D) %*% matrix(rnorm(r * p), r)
    )
  )
  cid <- rep(seq_len(G), each = D)
  B <- solve(crossprod(Z) + diag(seq_len(p)))
  z <- rnorm(nrow(Z))
  Xp <- matrix(rnorm(5 * p), 5)
  for (adjustment in c("exact", "shortcut", "none")) {
    core <- pffr_influence_core(Z, B, cid, z, adjustment)
    testthat::expect_equal(
      pffr_influence_df(core, Xp, df_gram = "diagonal")$df,
      diagonal_df_reference(Z, B, cid, Xp, adjustment),
      tolerance = 1e-10
    )
    # The default is unchanged by the new argument.
    testthat::expect_identical(
      pffr_influence_df(core, Xp, df_gram = "full")$df,
      pffr_influence_df(core, Xp)$df
    )
    testthat::expect_equal(
      pffr_influence_df(core, Xp, 1L, df_gram = "diagonal")$df,
      pffr_influence_df(core, Xp, df_gram = "diagonal")$df
    )
  }
})

testthat::test_that("diagonal and residualized df differ by the dropped cross terms", {
  set.seed(84002)
  cid <- rep(1:13, times = 2:14)
  n <- length(cid)
  dat <- data.frame(y = rnorm(n), x = rnorm(n), v = rnorm(n))
  fit <- lm(y ~ x + v, data = dat)
  Z <- model.matrix(fit)
  B <- solve(crossprod(Z))
  core <- pffr_influence_core(Z, B, cid, residuals(fit))
  Xp <- diag(ncol(Z))
  testthat::expect_equal(
    pffr_influence_df(core, Xp, df_gram = "diagonal")$df,
    diagonal_df_reference(Z, B, cid, Xp),
    tolerance = 1e-10
  )
  testthat::expect_false(isTRUE(all.equal(
    pffr_influence_df(core, Xp, df_gram = "diagonal")$df,
    pffr_influence_df(core, Xp)$df
  )))
  # Balanced intercept-only OLS: the residualized moment df is the known
  # G - 1, the historical diagonal shortcut returns G.
  for (G in c(3L, 7L, 20L)) {
    core <- pffr_influence_core(matrix(1, G, 1), matrix(1 / G), seq_len(G))
    testthat::expect_equal(
      pffr_influence_df(core, matrix(1), df_gram = "diagonal")$df,
      as.numeric(G),
      tolerance = 1e-10
    )
    testthat::expect_equal(
      pffr_influence_df(core, matrix(1))$df,
      G - 1,
      tolerance = 1e-10
    )
    testthat::expect_true(is.na(
      pffr_influence_df(core, matrix(0), df_gram = "diagonal")$df
    ))
  }
})

testthat::test_that("df_gram = 'diagonal' is historical only with the shortcut weight", {
  set.seed(84003)
  G <- 9L
  p <- 8L
  D <- 11L
  r <- 5L
  Z <- do.call(
    rbind,
    lapply(
      seq_len(G),
      function(g) matrix(rnorm(D * r), D) %*% matrix(rnorm(r * p), r)
    )
  )
  cid <- rep(seq_len(G), each = D)
  B <- solve(crossprod(Z) + diag(seq_len(p)) / 10)
  Xp <- matrix(rnorm(6 * p), 6)
  historical <- historical_satterthwaite_df(Z, cid, B, Xp, use_cl2 = TRUE)$df
  short <- pffr_influence_core(Z, B, cid, adjustment = "shortcut")
  testthat::expect_equal(
    pffr_influence_df(short, Xp, df_gram = "diagonal")$df,
    historical,
    tolerance = 1e-10
  )
  # The exact block is a different leverage weight, so its diagonal df is NOT
  # the historical number (the documented caveat on df_gram = "diagonal").
  exact <- pffr_influence_core(Z, B, cid, adjustment = "exact")
  testthat::expect_false(isTRUE(all.equal(
    pffr_influence_df(exact, Xp, df_gram = "diagonal")$df,
    historical
  )))
  # The CR1 path (no adjustment) matches the historical use_cl2 = FALSE kernel.
  testthat::expect_equal(
    pffr_influence_df(
      pffr_influence_core(Z, B, cid, adjustment = "none"),
      Xp,
      df_gram = "diagonal"
    )$df,
    historical_satterthwaite_df(Z, cid, B, Xp, use_cl2 = FALSE)$df,
    tolerance = 1e-10
  )
})

testthat::test_that("an indefinite bread trips the hat lower-bound monitors", {
  # Study-LB P-LB5 lower bounds: a Vp with a negative eigenvalue makes the
  # penalized hat indefinite. The upper-bound monitors (max h_ii, max
  # eigen(H_gg)) never see it, so only min_obs_leverage / min_hat_eig can.
  G <- 4L
  Z <- do.call(rbind, rep(list(diag(2)), G))
  cid <- rep(seq_len(G), each = 2L)
  B <- diag(c(0.4, -0.5))
  z <- as.numeric(seq_len(nrow(Z)))
  for (adjustment in c("exact", "shortcut", "none")) {
    core <- pffr_influence_core(Z, B, cid, z, adjustment)
    V <- pffr_influence_vcov(core)
    testthat::expect_equal(attr(V, "min_hat_eig"), -0.5)
    testthat::expect_equal(attr(V, "min_obs_leverage"), -0.5)
    testthat::expect_lte(attr(V, "max_obs_leverage"), 1)
    violation <- attr(V, "hat_invariant_violation")
    testthat::expect_type(violation, "character")
    testthat::expect_match(violation, "below the bound 0")
    testthat::expect_match(violation, "positive semi-definite by construction")
    # The covariance builder turns it into exactly one warning.
    w <- testthat::capture_warnings(
      gam_sandwich_cluster_cl2(NULL, NULL, influence = core)
    )
    testthat::expect_length(w, 1L)
    testthat::expect_match(w, "NOT trustworthy")
  }
  # A well-behaved bread leaves every monitor inside its bounds.
  ok <- pffr_influence_core(Z, diag(c(0.4, 0.5)), cid, z, "exact")
  testthat::expect_null(attr(
    pffr_influence_vcov(ok),
    "hat_invariant_violation"
  ))
})

testthat::test_that("hat-invariant lower bounds are reported and tolerant", {
  viol <- pffr_hat_invariant_violation
  testthat::expect_null(viol(min_obs_leverage = 0, min_hat_eig = 0))
  testthat::expect_null(viol(min_obs_leverage = -1e-9, min_hat_eig = -1e-9))
  testthat::expect_null(viol(min_obs_leverage = NA_real_, min_hat_eig = NaN))
  testthat::expect_match(
    viol(min_obs_leverage = -0.25),
    "smallest per-observation leverage is -0.25"
  )
  testthat::expect_match(
    viol(min_hat_eig = -3.5),
    "smallest per-cluster hat eigenvalue is -3.5"
  )
})

testthat::test_that("max_block_kappa uses |max eig| / |min eig|", {
  # Upstream's definition: the ratio is built from the residual block's
  # largest and smallest eigenvalues BY VALUE, each taken in absolute value,
  # so max_block_kappa itself is always non-negative -- including on the
  # indefinite block below, where it is |1.5| / |-0.5| = 3. What records the
  # indefiniteness is the signed min_block_eig (-0.5), asserted alongside it.
  G <- 4L
  Z <- do.call(rbind, rep(list(diag(2)), G))
  cid <- rep(seq_len(G), each = 2L)
  core <- pffr_influence_core(
    Z,
    diag(c(0.4, -0.5)),
    cid,
    adjustment = "shortcut"
  )
  d <- core$diagnostics
  # residual block I - H_gg = diag(0.6, 1.5): both positive, kappa = 1.5/0.6
  testthat::expect_equal(unique(d$max_block_kappa), 1.5 / 0.6)
  # An indefinite residual block: I - H_gg = diag(-0.5, 1.5) from h = (1.5, -0.5)
  core2 <- pffr_influence_core(
    Z,
    diag(c(1.5, -0.5)),
    cid,
    adjustment = "shortcut"
  )
  testthat::expect_equal(unique(core2$diagnostics$max_block_kappa), 1.5 / 0.5)
  testthat::expect_equal(unique(core2$diagnostics$min_block_eig), -0.5)
})

testthat::test_that("undefined moment df yields missing limits with one warning", {
  # A zero contrast has zero sampling variance, so the moment df is undefined.
  # compute_pointwise_ci() must report NA limits, not silently substitute the
  # Gaussian quantile.
  set.seed(84004)
  G <- 6L
  Z <- matrix(rnorm(G * 3L * 4L), G * 3L)
  cid <- rep(seq_len(G), each = 3L)
  B <- solve(crossprod(Z) + diag(4) / 10)
  ctx <- list(
    ok = TRUE,
    type = "cl2",
    core = pffr_influence_core(Z, B, cid, rnorm(nrow(Z)), "exact"),
    Vp = B,
    G = G,
    df_gram = "full"
  )
  linear_map <- list(X = matrix(0, 3L, 4L), trmind = seq_len(4L))
  w <- testthat::capture_warnings(
    pw <- compute_pointwise_ci("satterthwaite", 0.95, linear_map, ctx)
  )
  testthat::expect_length(w, 1L)
  testthat::expect_match(w, "Undefined working-model moment df")
  testthat::expect_true(all(is.na(pw$crit)))
  testthat::expect_true(all(is.na(pw$df)))
  # A nonzero contrast is unaffected.
  linear_map$X <- matrix(1, 3L, 4L)
  pw_ok <- compute_pointwise_ci("satterthwaite", 0.95, linear_map, ctx)
  testthat::expect_true(all(is.finite(pw_ok$crit)))
})
