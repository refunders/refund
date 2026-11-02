# Independent dense reference: full N-dimensional residualization, no
# compression, exact CL2 block ((I - H)^2)_gg.
dense_reference <- function(Z, B, cid, z, contrasts, cap = .999) {
  H <- Z %*% B %*% t(Z)
  L <- diag(nrow(Z)) - H
  ids <- unique(cid)
  K <- matrix(0, ncol(Z), length(ids))
  qs <- vector("list", length(ids))
  for (g in seq_along(ids)) {
    ii <- which(cid == ids[g])
    Zg <- Z[ii, , drop = FALSE]
    Rg <- tcrossprod(L[ii, , drop = FALSE])
    ee <- eigen((Rg + t(Rg)) / 2, symmetric = TRUE)
    A <- tcrossprod(
      sweep(ee$vectors, 2, sqrt(pmax(ee$values, (1 - cap)^2)), "/"),
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
  core <- pffr_influence_core(Z, B, cid, z)
  ref <- dense_reference(Z, B, cid, z, Xp)
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
  other <- pffr_influence_core(Z[ord, ], B, cid[ord], z[ord])
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

testthat::test_that("the covariance is the scaled sampling part plus B2", {
  set.seed(84003)
  Z <- matrix(rnorm(120), 30)
  B <- solve(crossprod(Z) + diag(4))
  core <- pffr_influence_core(Z, B, rep(1:10, each = 3), rnorm(30))
  core$B2 <- B %*% B
  V <- pffr_influence_vcov(core)
  testthat::expect_equal(
    matrix(V, nrow(V)),
    core$correction * tcrossprod(core$K) + core$B2
  )
  testthat::expect_equal(core$correction, 10 / 9)
  testthat::expect_identical(attr(V, "cl2_adjustment"), "exact")
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
  core <- pffr_influence_core(Z, B, cid, z)
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
  # A well-behaved bread leaves every monitor inside its bounds.
  ok <- pffr_influence_core(Z, diag(c(0.4, 0.5)), cid, z)
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

testthat::test_that("block diagnostics describe the exact residual block", {
  # Identity designs per cluster: H_gg = B and the exact block is
  # I - (2B - B Z'Z B) = I - 2B + G B^2 (here G = 4), i.e. diag(0.84, 3) for
  # B = diag(0.4, -0.5). max_block_kappa is |max eig| / |min eig|.
  G <- 4L
  Z <- do.call(rbind, rep(list(diag(2)), G))
  cid <- rep(seq_len(G), each = 2L)
  d <- pffr_influence_core(Z, diag(c(0.4, -0.5)), cid)$diagnostics
  testthat::expect_equal(unique(d$min_block_eig), 0.84)
  testthat::expect_equal(unique(d$max_block_kappa), 3 / 0.84)
  testthat::expect_equal(unique(d$min_block_eig_rel), 0.84 / 3)
})

testthat::test_that("undefined moment df falls back to z with one message", {
  # A zero contrast has zero sampling variance, so the moment df is undefined;
  # the Gaussian critical value is used there, with a message.
  set.seed(84004)
  G <- 6L
  Z <- matrix(rnorm(G * 3L * 4L), G * 3L)
  cid <- rep(seq_len(G), each = 3L)
  B <- solve(crossprod(Z) + diag(4) / 10)
  setup <- list(
    mode = "satterthwaite",
    core = pffr_influence_core(Z, B, cid, rnorm(nrow(Z)))
  )
  Xp <- rbind(matrix(0, 2L, 4L), 1)
  testthat::expect_message(
    pw <- pffr_pointwise_crit_values(setup, 0.95, Xp),
    "Satterthwaite df are undefined"
  )
  testthat::expect_equal(pw$df[1:2], c(Inf, Inf))
  testthat::expect_equal(pw$crit[1:2], rep(qnorm(0.975), 2))
  testthat::expect_true(is.finite(pw$df[3]) && pw$df[3] < G)
  testthat::expect_equal(pw$crit[3], qt(0.975, pw$df[3]))
  # Inside a collection window the message is emitted once, at the end.
  msgs <- testthat::capture_messages({
    opened <- pffr_begin_undefined_df()
    pffr_pointwise_crit_values(setup, 0.95, Xp)
    pffr_pointwise_crit_values(setup, 0.95, Xp)
    pffr_end_undefined_df(opened)
  })
  testthat::expect_length(msgs, 1L)
})

testthat::test_that("moment df with many unequal-rank clusters matches the dense reference", {
  # Covers the two cases the other fixtures do not: G above the default chunk
  # size (so the df is assembled from several chunks) and unequal per-cluster
  # ranks, including a rank-deficient cluster. Both are what the cached
  # per-cluster residualization blocks R T_g' have to get right.
  set.seed(84005)
  G <- 45L
  p <- 9L
  sizes <- rep(c(3L, 5L, 8L, 11L, 6L), length.out = G)
  cid <- rep(seq_len(G), times = sizes)
  Z <- matrix(rnorm(length(cid) * p), ncol = p)
  # one exactly rank-deficient cluster: every row repeated (row-wise; a plain
  # block assignment would recycle the source row column-major instead)
  dup <- which(cid == 4L)
  Z[dup, ] <- matrix(Z[dup[1L], ], length(dup), p, byrow = TRUE)
  B <- solve(crossprod(Z) + diag(seq_len(p)) / 10)
  z <- rnorm(nrow(Z))
  Xp <- matrix(rnorm(37L * p), 37L)
  core <- pffr_influence_core(Z, B, cid, z)
  testthat::expect_true(core$df_precompute)
  testthat::expect_gt(diff(range(core$diagnostics$rank)), 0)
  testthat::expect_equal(core$diagnostics$rank[4L], 1L)
  ref <- dense_reference(Z, B, cid, z, Xp)
  # G > p here, so the Frobenius norm of the Gram comes from the p x p side.
  testthat::expect_equal(
    pffr_influence_df(core, Xp)$df,
    ref$df,
    tolerance = 1e-10
  )
  # chunking stays order invariant across and beyond the chunk boundary
  for (chunk in c(1L, 7L, 32L, 10000L))
    testthat::expect_equal(
      pffr_influence_df(core, Xp, chunk)$df,
      ref$df,
      tolerance = 1e-10
    )
  # the fallback path (no cached blocks) returns the same numbers
  plain <- pffr_influence_core(Z, B, cid, z, df_precompute_bytes = 0)
  testthat::expect_false(plain$df_precompute)
  testthat::expect_null(plain$blocks[[1L]]$RT)
  testthat::expect_equal(
    pffr_influence_df(plain, Xp)$df,
    pffr_influence_df(core, Xp)$df,
    tolerance = 1e-10
  )
  # A bread that makes C indefinite has no usable Cholesky factor: the
  # precompute switches itself off instead of erroring, and the df still comes
  # out of the general path.
  bad <- pffr_influence_core(Z, -B, cid, z)
  testthat::expect_false(bad$df_precompute)
  testthat::expect_length(pffr_influence_df(bad, Xp)$df, nrow(Xp))
})

testthat::test_that("a zero df cache budget skips chol(C) entirely (Copilot review of #126)", {
  # The budget test `8 * p * sum(rank) <= df_precompute_bytes` must run before
  # C is factored, not after: with df_precompute_bytes = 0 the fallback is
  # unreachable-by-budget, so chol(C) must never be attempted (it would be a
  # wasted O(p^3) cost on a large fit, and defeats the documented opt-out).
  set.seed(84006)
  G <- 6L
  p <- 8L
  D <- 5L
  r <- 3L
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

  default_budget <- pffr_influence_core(Z, B, cid, z)
  testthat::expect_true(default_budget$df_precompute)

  testthat::local_mocked_bindings(
    chol = function(...) stop("chol must not be called"),
    .package = "base"
  )
  zero_budget <- pffr_influence_core(Z, B, cid, z, df_precompute_bytes = 0)
  testthat::expect_false(zero_budget$df_precompute)
  testthat::expect_null(zero_budget$blocks[[1L]]$RT)

  df_default <- pffr_influence_df(default_budget, Xp)
  df_zero <- pffr_influence_df(zero_budget, Xp)
  testthat::expect_equal(df_zero$df, df_default$df, tolerance = 1e-12)
  testthat::expect_equal(
    df_zero$expected_sampling_variance,
    df_default$expected_sampling_variance,
    tolerance = 1e-12
  )
})
