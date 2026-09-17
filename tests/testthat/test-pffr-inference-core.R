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

# Historical working-iid diagonal df: the deleted satterthwaite_df_kernel()
# formula, kept available through df_gram = "diagonal" for re-scoring
# comparisons only. Re-implemented here directly from that formula.
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

testthat::test_that("df_gram = 'diagonal' reproduces the historical shortcut", {
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
