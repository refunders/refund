#--------------------------------------
# Tests for the Satterthwaite critical values of the default CL2 intervals:
# kernel correctness, monotonicity, bounds, and their use in coef(),
# predict() and pffr_predict_ci().
#--------------------------------------

test_that("balanced unpenalized OLS with identical clusters gives df ~ G-1", {
  # Identical intercept-only clusters, unpenalized OLS bread (Vp = (X'X)^{-1}).
  # Full residualization recovers G-1 even at small G.
  G <- 7
  Xw <- matrix(1, nrow = G, ncol = 1)
  core <- pffr_influence_core(Xw, solve(crossprod(Xw)), seq_len(G))
  k <- pffr_influence_df(core, matrix(1, nrow = 1, ncol = 1))
  expect_equal(k$G, G)
  expect_equal(k$df, G - 1, tolerance = 1e-10)

  # identical multi-observation blocks (D=5, p=3) also give exactly G-1
  set.seed(11)
  Gb <- 40
  Z <- matrix(rnorm(5 * 3), 5, 3)
  Xwb <- do.call(rbind, replicate(Gb, Z, simplify = FALSE))
  coreb <- pffr_influence_core(
    Xwb,
    solve(crossprod(Xwb)),
    rep(seq_len(Gb), each = 5)
  )
  Xpb <- matrix(rnorm(4 * 3), 4, 3)
  expect_equal(
    pffr_influence_df(coreb, Xpb)$df,
    rep(Gb - 1, 4),
    tolerance = 1e-8
  )
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
    core <- pffr_influence_core(Xw, solve(crossprod(Xw)), cid)
    pffr_influence_df(core, Xp)$df
  }
  dfs <- vapply(c(1, 4, 8, 16, 32), df_at_scale, numeric(1))

  # strictly decreasing as cluster 1's leverage grows
  expect_true(all(diff(dfs) < 0))
  # and it collapses toward 1 (that one cluster dominates the meat)
  expect_lt(dfs[length(dfs)], 2)
  expect_gte(min(dfs), 1)
})

test_that("per-point df stays within [1, G] on assorted pffr fits", {
  skip_on_cran()

  check_bounds <- function(fit) {
    G <- length(unique(build_cluster_id(fit$pffr)))
    co <- suppressMessages(coef(
      fit,
      ci = "pointwise",
      sandwich = TRUE,
      n1 = 25,
      n2 = 12
    ))
    expect_identical(co$ci_meta$crit_used, "satterthwaite")
    dfv <- unlist(lapply(co$smterms, function(x) x$coef$df))
    dfv <- dfv[is.finite(dfv)]
    expect_gt(length(dfv), 0)
    expect_true(all(dfv >= 1 - 1e-8))
    expect_true(all(dfv <= G + 1e-6))
    dfv
  }

  d_cl2 <- check_bounds(get_basic_pffr_model())
  # gaulss two-block whitened path
  check_bounds(get_gaulss_model())
  # medians land in a sensible interior range for these G ~ 30 fits
  expect_gt(stats::median(d_cl2), 2)
})

test_that("CL2 intervals use Satterthwaite, model-based ones Gaussian values", {
  skip_on_cran()
  m <- get_basic_pffr_model()

  co <- suppressMessages(coef(m, ci = "pointwise", sandwich = TRUE, n1 = 30))
  sm <- co$smterms[[1]]$coef
  expect_identical(co$ci_meta$crit_used, "satterthwaite")
  expect_true(all(is.finite(sm$df)))
  expect_equal(sm$upper - sm$value, stats::qt(0.975, sm$df) * sm$se)
  expect_equal(
    unname(co$pterms[, "upper"] - co$pterms[, "value"]),
    unname(stats::qt(0.975, co$pterms[, "df"]) * co$pterms[, "se"])
  )

  # The fixture is model-based (sandwich = FALSE): Gaussian critical values.
  co_mb <- suppressMessages(coef(m, ci = "pointwise", n1 = 30))
  sm_mb <- co_mb$smterms[[1]]$coef
  expect_identical(co_mb$ci_meta$crit_used, "z")
  expect_false(co_mb$ci_meta$sandwich)
  expect_true(all(is.infinite(sm_mb$df)))
  expect_equal(sm_mb$upper - sm_mb$value, stats::qnorm(0.975) * sm_mb$se)
  expect_true(all(is.infinite(co_mb$pterms[, "df"])))
})

test_that("predict() and pffr_predict_ci() return per-point crit and df", {
  skip_on_cran()
  m <- get_basic_pffr_model()
  dat <- get_basic_pffr_data()
  s <- attr(dat, "xindex")
  fit_cl2 <- quiet_pffr(
    Y ~ ff(X1, xind = s) + xlin,
    yind = attr(dat, "yindex"),
    data = dat
  )
  p <- predict(fit_cl2, se.fit = TRUE, reformat = FALSE)
  expect_named(p, c("fit", "se.fit", "crit", "df"))
  expect_true(all(is.finite(p$df) & p$df >= 1))
  expect_equal(p$crit, stats::qt(0.975, p$df))

  p_90 <- predict(fit_cl2, se.fit = TRUE, reformat = FALSE, level = 0.9)
  expect_equal(p_90$crit, stats::qt(0.95, p$df))
  # reformatted output: matrices like fit
  p_mat <- predict(fit_cl2, se.fit = TRUE)
  expect_equal(dim(p_mat$df), dim(p_mat$fit))

  ci <- pffr_predict_ci(fit_cl2)
  expect_identical(attr(ci, "crit_used"), "satterthwaite")
  expect_equal(ci$df, unname(p$df))
  expect_equal(ci$upper - ci$fit, stats::qt(0.975, ci$df) * ci$se_link)

  # model-based covariances use Gaussian critical values
  expect_identical(
    attr(pffr_predict_ci(fit_cl2, sandwich = FALSE), "crit_used"),
    "z"
  )
  expect_true(all(is.infinite(predict(m, se.fit = TRUE, reformat = FALSE)$df)))
})

test_that("Satterthwaite df in predictions are capped by prediction points", {
  skip_on_cran()
  fit_cl2 <- quiet_pffr(
    Y ~ xlin,
    yind = attr(get_basic_pffr_data(), "yindex"),
    data = get_basic_pffr_data()
  )
  old <- options(refund.pffr_satterthwaite_max_points = 10)
  on.exit(options(old), add = TRUE)
  expect_message(
    p <- predict(fit_cl2, se.fit = TRUE, reformat = FALSE),
    "Using Gaussian critical values"
  )
  expect_true(all(is.infinite(p$df)))
  expect_equal(p$crit, rep(stats::qnorm(0.975), length(p$fit)))
})

test_that("coef() and predict() give the same df for the same contrasts", {
  skip_on_cran()
  dat <- get_basic_pffr_data()
  t <- attr(dat, "yindex")
  # Without an intercept, the prediction rows at xlin = 1 are exactly the
  # coefficient-function rows of xlin(t) that coef() evaluates.
  fit <- quiet_pffr(Y ~ 0 + xlin, yind = t, data = dat)
  co <- coef(fit, ci = "pointwise", n1 = length(t))
  cf <- co$smterms[[1]]$coef
  expect_equal(cf[[1]], t)

  nd <- list(xlin = 1)
  p <- predict(fit, newdata = nd, se.fit = TRUE, reformat = FALSE)
  ci <- pffr_predict_ci(fit, newdata = nd)
  expect_equal(as.vector(p$se.fit), cf$se, tolerance = 1e-8)
  expect_equal(as.vector(p$df), cf$df, tolerance = 1e-8)
  expect_equal(ci$df, cf$df, tolerance = 1e-8)
  expect_equal(ci$lower, cf$lower, tolerance = 1e-8)
})

test_that("removed critical-value and covariance options error informatively", {
  skip_on_cran()
  m <- get_basic_pffr_model()
  for (arg in c("crit", "cl2_adjustment", "df_gram", "bias_ref", "ci_ref")) {
    args <- list(m, ci = "pointwise")
    args[[arg]] <- "z"
    expect_error(do.call(coef, args), "no longer exist")
  }
  expect_error(predict(m, se.fit = TRUE, crit = "z"), "no longer exist")
  expect_error(
    predict(m, se.fit = TRUE, se_method = "jackknife"),
    "no longer exist"
  )
  expect_false(any(
    c("crit", "cl2_adjustment", "bias_ref", "dof_correction") %in%
      c(
        names(formals(coef.pffr)),
        names(formals(predict.pffr)),
        names(formals(pffr_predict_ci)),
        names(formals(pffr))
      )
  ))
})
