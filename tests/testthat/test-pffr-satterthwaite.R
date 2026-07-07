#--------------------------------------
# Tests for the Satterthwaite degrees-of-freedom machinery (S3)
#
# Per-point Bell-McCaffrey df for cluster-robust (CR1) and leverage-adjusted
# (CL2) pointwise intervals: kernel correctness, monotonicity, bounds, and the
# coef.pffr `crit` API (auto / z / tG1 / satterthwaite).
#--------------------------------------

test_that("balanced unpenalized OLS with identical clusters gives df ~ G-1", {
  # Identical intercept-only clusters, unpenalized OLS bread (Vp = (X'X)^{-1}).
  # The working-iid shortcut returns exactly G for perfectly identical clusters
  # (vs. the exact Bell-McCaffrey G-1, which needs the dropped cross-cluster
  # residual terms); at large G, G is within 1% of G-1.
  G <- 200
  Xw <- matrix(1, nrow = G, ncol = 1)
  cid <- seq_len(G)
  Vp <- solve(crossprod(Xw))
  Xp <- matrix(1, nrow = 1, ncol = 1)

  k_cr1 <- satterthwaite_df_kernel(Xw, cid, Vp, Xp, use_cl2 = FALSE)
  k_cl2 <- satterthwaite_df_kernel(Xw, cid, Vp, Xp, use_cl2 = TRUE)

  expect_equal(k_cr1$G, G)
  # within 1% of G-1 (the classic balanced-design target)
  expect_equal(k_cr1$df, G - 1, tolerance = 0.01)
  expect_equal(k_cl2$df, G - 1, tolerance = 0.01)
  # honesty: the shortcut is exactly G for identical clusters
  expect_equal(k_cr1$df, G, tolerance = 1e-8)
  expect_equal(k_cl2$df, G, tolerance = 1e-8)

  # identical multi-observation blocks (D=5, p=3) also give exactly G
  set.seed(11)
  Gb <- 40
  Z <- matrix(rnorm(5 * 3), 5, 3)
  Xwb <- do.call(rbind, replicate(Gb, Z, simplify = FALSE))
  cidb <- rep(seq_len(Gb), each = 5)
  Vpb <- solve(crossprod(Xwb))
  Xpb <- matrix(rnorm(4 * 3), 4, 3)
  kb <- satterthwaite_df_kernel(Xwb, cidb, Vpb, Xpb, use_cl2 = TRUE)
  expect_equal(kb$df, rep(Gb, 4), tolerance = 1e-8)
})

test_that("df decreases monotonically as one cluster's leverage is inflated", {
  set.seed(101)
  G <- 30
  D <- 4
  p <- 2
  Xw0 <- matrix(rnorm(G * D * p), G * D, p)
  cid <- rep(seq_len(G), each = D)
  Xp <- matrix(rnorm(p), 1, p)

  df_at_scale <- function(scale) {
    Xw <- Xw0
    Xw[cid == 1, ] <- Xw[cid == 1, ] * scale
    Vp <- solve(crossprod(Xw))
    satterthwaite_df_kernel(Xw, cid, Vp, Xp, use_cl2 = TRUE)$df
  }
  # large-weight inflation of cluster 1 (the brief's construction): df falls
  # monotonically once the inflated cluster dominates the leverage.
  dfs <- vapply(c(1, 4, 8, 16, 32), df_at_scale, numeric(1))

  # strictly decreasing as cluster 1's leverage grows
  expect_true(all(diff(dfs) < 0))
  # and it collapses toward 1 (that one cluster dominates the meat)
  expect_lt(dfs[length(dfs)], 2)
  expect_gte(min(dfs), 1)
})

test_that("per-point df stays within [1, G] on assorted pffr fits", {
  skip_on_cran()

  check_bounds <- function(fit, type) {
    G <- length(unique(build_cluster_id(fit$pffr)))
    co <- suppressMessages(coef(
      fit,
      ci = "pointwise",
      crit = "satterthwaite",
      sandwich = type,
      n1 = 25,
      n2 = 12
    ))
    dfv <- unlist(lapply(co$smterms, function(x) x$coef$df))
    dfv <- dfv[is.finite(dfv)]
    expect_gt(length(dfv), 0)
    expect_true(all(dfv >= 1 - 1e-8))
    expect_true(all(dfv <= G + 1e-6))
    dfv
  }

  # cluster (CR1) and cl2 (leverage-adjusted) on an ff + xlin fit
  m <- get_basic_pffr_model()
  d_cr1 <- check_bounds(m, "cluster")
  d_cl2 <- check_bounds(m, "cl2")

  # gaulss two-block whitened path
  mg <- get_gaulss_model()
  d_g <- check_bounds(mg, "cluster")

  # medians land in a sensible interior range for these G ~ 30 fits
  expect_gt(stats::median(d_cl2), 2)
})

test_that("crit = 'z' reproduces the Gaussian pointwise intervals exactly", {
  skip_on_cran()
  m <- get_basic_pffr_model()

  co_z <- coef(m, ci = "pointwise", crit = "z", sandwich = "cl2", n1 = 30)
  sm <- co_z$smterms[[1]]$coef
  z <- stats::qnorm(0.975)

  expect_equal(sm$upper - sm$value, z * sm$se)
  expect_equal(sm$value - sm$lower, z * sm$se)
  expect_true(all(is.infinite(sm$df)))
  expect_identical(co_z$ci_meta$crit_used, "z")

  # parametric terms too
  expect_true("df" %in% colnames(co_z$pterms))
  expect_true(all(is.infinite(co_z$pterms[, "df"])))
})

test_that("crit = 'tG1' uses a constant t_{G-1} reference", {
  skip_on_cran()
  m <- get_basic_pffr_model()
  G <- length(unique(build_cluster_id(m$pffr)))

  co <- coef(m, ci = "pointwise", crit = "tG1", sandwich = "cl2", n1 = 30)
  sm <- co$smterms[[1]]$coef

  expect_true(all(sm$df == G - 1))
  expect_equal(sm$upper - sm$value, stats::qt(0.975, G - 1) * sm$se)
  # t_{G-1} is wider than the Gaussian reference
  expect_gt(stats::qt(0.975, G - 1), stats::qnorm(0.975))
})

test_that("crit = 'auto' resolves per sandwich path and G", {
  skip_on_cran()
  m <- get_basic_pffr_model()

  # cl2 fit, G ~ 30 < 150  ->  satterthwaite
  co_auto <- coef(m, ci = "pointwise", sandwich = "cl2", n1 = 30)
  co_sat <- coef(
    m,
    ci = "pointwise",
    crit = "satterthwaite",
    sandwich = "cl2",
    n1 = 30
  )
  expect_identical(co_auto$ci_meta$crit_used, "satterthwaite")
  expect_equal(co_auto$smterms[[1]]$coef$upper, co_sat$smterms[[1]]$coef$upper)
  expect_equal(co_auto$smterms[[1]]$coef$df, co_sat$smterms[[1]]$coef$df)

  # sandwich = "none"  ->  z (no cluster structure)
  co_none <- coef(m, ci = "pointwise", sandwich = "none", n1 = 30)
  expect_identical(co_none$ci_meta$crit_used, "z")
  expect_true(all(is.infinite(co_none$smterms[[1]]$coef$df)))
})

test_that("crit = 'satterthwaite' on a non-cluster covariance degrades to z", {
  skip_on_cran()
  m <- get_basic_pffr_model()

  expect_warning(
    co <- coef(
      m,
      ci = "pointwise",
      crit = "satterthwaite",
      sandwich = "none",
      n1 = 20
    ),
    "requires a cluster-robust covariance"
  )
  expect_identical(co$ci_meta$crit_used, "z")
  expect_true(all(is.infinite(co$smterms[[1]]$coef$df)))
})

test_that("CL2 df is never above the CR1 df (leverage adjustment shrinks df)", {
  skip_on_cran()
  # A_g = (I - H_gg)^{-1/2} inflates ||q_g|| for high-leverage clusters, which
  # can only lower (or leave equal) the Satterthwaite df relative to A_g = I.
  m <- get_basic_pffr_model()
  co_cr1 <- suppressMessages(coef(
    m,
    ci = "pointwise",
    crit = "satterthwaite",
    sandwich = "cluster",
    n1 = 25
  ))
  co_cl2 <- suppressMessages(coef(
    m,
    ci = "pointwise",
    crit = "satterthwaite",
    sandwich = "cl2",
    n1 = 25
  ))
  df_cr1 <- co_cr1$smterms[[1]]$coef$df
  df_cl2 <- co_cl2$smterms[[1]]$coef$df
  ok <- is.finite(df_cr1) & is.finite(df_cl2)
  expect_true(all(df_cl2[ok] <= df_cr1[ok] + 1e-6))
})

test_that("summary() reports median/min Satterthwaite df for cluster-robust fits", {
  skip_on_cran()
  dat <- get_basic_pffr_data()
  s <- attr(dat, "xindex")
  t <- attr(dat, "yindex")
  m_cl2 <- pffr(
    Y ~ ff(X1, xind = s) + xlin,
    yind = t,
    data = dat,
    sandwich = "cl2"
  )
  sm <- summary(m_cl2)
  expect_false(is.null(sm$satterthwaite_df))
  expect_identical(sm$satterthwaite_df$type, "cl2")
  expect_true(sm$satterthwaite_df$median >= sm$satterthwaite_df$min)
  expect_true(sm$satterthwaite_df$min >= 1)
  out <- capture.output(print(sm))
  expect_true(any(grepl("Satterthwaite df", out)))

  # sandwich = "none" fits print no df line
  m_none <- get_basic_pffr_model()
  expect_null(summary(m_none)$satterthwaite_df)
})
