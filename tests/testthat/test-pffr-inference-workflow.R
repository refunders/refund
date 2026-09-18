testthat::test_that("subject clustering and covariance survive coef predict plot", {
  set.seed(84101)
  G <- 20L
  D <- 18L
  subject <- rep(seq_len(G), times = rep(c(1L, 2L), length.out = G))
  n <- length(subject)
  tt <- seq(0, 1, length.out = D)
  dat <- list(Y = matrix(rnorm(n * D), n, D), x = rnorm(n), subject = subject)
  dat$Y <- dat$Y + outer(dat$x, sin(2 * pi * tt))
  fit <- suppressMessages(quiet_pffr(
    Y ~ x,
    yind = tt,
    data = dat,
    bs.yindex = list(bs = "ps", k = 6, m = c(2, 1)),
    bs.int = list(bs = "ps", k = 6, m = c(2, 1)),
    cluster = subject
  ))
  testthat::expect_identical(fit$pffr$sandwich_info$type, "cl2")
  testthat::expect_identical(fit$pffr$sandwich_info$G, G)
  testthat::expect_equal(fit$pffr$sandwich_info$cluster_var, subject)
  testthat::expect_equal(build_cluster_id(fit$pffr), rep(subject, each = D))
  co <- coef(fit, ci = "pointwise", n1 = 18)
  explicit <- coef(
    fit,
    sandwich = "cl2",
    cluster = subject,
    crit = "z",
    ci = "pointwise",
    n1 = 18
  )
  testthat::expect_identical(co$ci_meta$crit_used, "z")
  testthat::expect_equal(co$smterms, explicit$smterms)
  V <- pffr_vcov(fit)
  B <- fit$Vp
  X <- predict(fit, type = "lpmatrix", reformat = FALSE)
  pr <- predict(fit, se.fit = TRUE, type = "link")
  testthat::expect_equal(
    as.vector(t(pr$se.fit)),
    unname(sqrt(rowSums((X %*% V) * X))),
    tolerance = 1e-9
  )
  new_pr <- predict(fit, newdata = dat, se.fit = TRUE, type = "link")
  testthat::expect_equal(new_pr$se.fit, pr$se.fit, tolerance = 1e-9)
  pdf_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(pdf_file)
  on.exit(
    {
      grDevices::dev.off()
      unlink(pdf_file)
    },
    add = TRUE
  )
  plt <- plot(fit, n = 18, se = TRUE, pages = 0)
  ref <- pffr_model_based_gam(fit)
  ref$Vp <- V
  ref$Vc <- V
  expected <- mgcv::plot.gam(ref, n = 18, se = TRUE, pages = 0)
  testthat::expect_equal(
    lapply(plt, function(x) x$se),
    lapply(expected, function(x) x$se),
    tolerance = 1e-9
  )
  testthat::expect_equal(fit$Vp, B)
  core <- pffr_influence(fit)
  testthat::expect_equal(
    pffr_influence_df(core, diag(ncol(B)))$df,
    satterthwaite_df_kernel(
      model.matrix(pffr_model_based_gam(fit)) * sqrt(1 / fit$sig2),
      build_cluster_id(fit$pffr),
      B,
      diag(ncol(B)),
      TRUE
    )$df,
    tolerance = 1e-8
  )
  testthat::expect_true(any(grepl("^influence", ls(fit$pffr$Vsandwich_cache))))
  testthat::expect_error(
    refund::pffr(
      Y ~ x,
      yind = tt,
      data = dat,
      cluster = replace(subject, 1, NA)
    ),
    "missing"
  )
  testthat::expect_error(
    refund::pffr(Y ~ x, yind = tt, data = dat, cluster = subject[-1]),
    "one entry per curve"
  )
  for (a in c("exact", "shortcut")) {
    V <- pffr_vcov(fit, cl2_adjustment = a)
    ctx <- pffr_df_context(fit, "cl2", cl2_adjustment = a)
    testthat::expect_identical(ctx$core$adjustment, attr(V, "cl2_adjustment"))
    co <- coef(
      fit,
      ci = "pointwise",
      crit = "satterthwaite",
      cl2_adjustment = a,
      n1 = 12
    )
    testthat::expect_true(all(is.finite(co$smterms[[1]]$coef$df)))
  }
  # --- df_gram end to end (coef.pffr -> pffr_df_context -> pffr_df_from_context)
  co_full <- coef(fit, ci = "pointwise", crit = "satterthwaite", n1 = 12)
  co_diag <- coef(
    fit,
    ci = "pointwise",
    crit = "satterthwaite",
    df_gram = "diagonal",
    n1 = 12
  )
  # "full" is the default.
  testthat::expect_equal(
    coef(
      fit,
      ci = "pointwise",
      crit = "satterthwaite",
      df_gram = "full",
      n1 = 12
    )$smterms[[1]]$coef$df,
    co_full$smterms[[1]]$coef$df
  )
  width <- function(co) {
    with(co$smterms[[1]]$coef, upper - lower)
  }
  for (tm in seq_along(co_full$smterms)) {
    testthat::expect_false(isTRUE(all.equal(
      co_full$smterms[[tm]]$coef$df,
      co_diag$smterms[[tm]]$coef$df
    )))
  }
  # The diagonal df is the larger one (it drops the residualization), so its
  # t quantile - and therefore every interval - is strictly narrower.
  testthat::expect_true(all(
    co_diag$smterms[[1]]$coef$df > co_full$smterms[[1]]$coef$df
  ))
  testthat::expect_true(all(width(co_diag) < width(co_full)))
  testthat::expect_false(isTRUE(all.equal(
    co_full$pterms[, "df"],
    co_diag$pterms[, "df"]
  )))
  # df_gram must survive the whole call chain, not be swallowed on the way:
  # the parametric contrasts are unit vectors, so they can be rebuilt exactly.
  smind <- unlist(lapply(
    fit$smooth,
    function(s) seq(s$first.para, s$last.para)
  ))
  pind <- seq_along(fit$coefficients)[-smind]
  Xp_p <- matrix(0, length(pind), length(fit$coefficients))
  Xp_p[cbind(seq_along(pind), pind)] <- 1
  ctx_default <- pffr_df_context(fit, "cl2")
  testthat::expect_identical(ctx_default$df_gram, "full")
  testthat::expect_equal(
    unname(co_diag$pterms[, "df"]),
    pffr_df_from_context(ctx_default, Xp_p, df_gram = "diagonal"),
    tolerance = 1e-10
  )
  testthat::expect_equal(
    unname(co_full$pterms[, "df"]),
    pffr_df_from_context(ctx_default, Xp_p),
    tolerance = 1e-10
  )
  # --- df_gram = "diagonal" is backwards compatible only with the shortcut
  b <- pffr_model_based_gam(fit)
  Xw <- model.matrix(b) * sqrt(1 / fit$sig2)
  cid <- build_cluster_id(fit$pffr)
  historical <- historical_satterthwaite_df(
    Xw,
    cid,
    b$Vp,
    Xp_p,
    use_cl2 = TRUE
  )$df
  co_hist <- coef(
    fit,
    sandwich = "cl2",
    cl2_adjustment = "shortcut",
    crit = "satterthwaite",
    df_gram = "diagonal",
    ci = "pointwise",
    n1 = 12
  )
  testthat::expect_equal(
    unname(co_hist$pterms[, "df"]),
    historical,
    tolerance = 1e-10
  )
  co_exact <- coef(
    fit,
    sandwich = "cl2",
    cl2_adjustment = "exact",
    crit = "satterthwaite",
    df_gram = "diagonal",
    ci = "pointwise",
    n1 = 12
  )
  # ... but the exact geometry gives a different diagonal df. (The intercept
  # contrasts happen to agree to ~1e-10 here, so assert on the x(yindex)
  # coefficient surface and on a random contrast set, where it is ~5e-3 / 2e-2.)
  testthat::expect_false(isTRUE(all.equal(
    co_hist$smterms[["x(yindex)"]]$coef$df,
    co_exact$smterms[["x(yindex)"]]$coef$df
  )))
  set.seed(84105)
  Xp_r <- matrix(rnorm(8 * ncol(b$Vp)), 8)
  testthat::expect_equal(
    pffr_df_from_context(
      pffr_df_context(
        fit,
        "cl2",
        cl2_adjustment = "shortcut",
        df_gram = "diagonal"
      ),
      Xp_r
    ),
    historical_satterthwaite_df(Xw, cid, b$Vp, Xp_r, use_cl2 = TRUE)$df,
    tolerance = 1e-10
  )
  testthat::expect_false(isTRUE(all.equal(
    pffr_df_from_context(
      pffr_df_context(
        fit,
        "cl2",
        cl2_adjustment = "exact",
        df_gram = "diagonal"
      ),
      Xp_r
    ),
    historical_satterthwaite_df(Xw, cid, b$Vp, Xp_r, use_cl2 = TRUE)$df
  )))
})

testthat::test_that("an explicit dof_correction override is not served from cache", {
  set.seed(84103)
  G <- 12L
  D <- 10L
  tt <- seq(0, 1, length.out = D)
  dat <- list(Y = matrix(rnorm(G * D), G, D), x = rnorm(G))
  dat$Y <- dat$Y + outer(dat$x, sin(2 * pi * tt))
  fit <- suppressMessages(quiet_pffr(
    Y ~ x,
    yind = tt,
    data = dat,
    bs.yindex = list(bs = "ps", k = 5, m = c(2, 1)),
    bs.int = list(bs = "ps", k = 5, m = c(2, 1)),
    sandwich = "cluster"
  ))
  Xp <- diag(length(fit$coefficients))
  none <- pffr_influence(fit, "cluster", dof_correction = "none")
  edf <- pffr_influence(fit, "cluster", dof_correction = "edf")
  testthat::expect_gt(edf$correction, none$correction)
  testthat::expect_false(isTRUE(all.equal(
    pffr_influence_df(none, Xp)$expected_sampling_variance,
    pffr_influence_df(edf, Xp)$expected_sampling_variance
  )))
  # The df itself is scale free, so only the expected sampling variance moves.
  testthat::expect_equal(
    pffr_influence_df(none, Xp)$df,
    pffr_influence_df(edf, Xp)$df
  )
  # Both live in the cache under distinct keys, and re-reading returns the
  # matching object rather than whichever was computed first.
  keys <- ls(fit$pffr$Vsandwich_cache)
  testthat::expect_true(any(grepl("\\|none\\|", keys)))
  testthat::expect_true(any(grepl("\\|edf\\|", keys)))
  testthat::expect_equal(
    pffr_influence(fit, "cluster", dof_correction = "none")$correction,
    none$correction
  )
  testthat::expect_equal(
    pffr_influence(fit, "cluster", dof_correction = "edf")$correction,
    edf$correction
  )
})

testthat::test_that("missing interval limits do not break plot or summary", {
  set.seed(84104)
  G <- 10L
  D <- 8L
  tt <- seq(0, 1, length.out = D)
  dat <- list(Y = matrix(rnorm(G * D), G, D), x = rnorm(G))
  dat$Y <- dat$Y + outer(dat$x, sin(2 * pi * tt))
  fit <- suppressMessages(quiet_pffr(
    Y ~ x,
    yind = tt,
    data = dat,
    bs.yindex = list(bs = "ps", k = 5, m = c(2, 1)),
    bs.int = list(bs = "ps", k = 5, m = c(2, 1)),
    sandwich = "cl2"
  ))
  ctx <- pffr_df_context(fit, "cl2")
  # No fitted-model contrast in coef() is exactly zero, so the undefined-df
  # branch is reached through compute_pointwise_ci() with a zero contrast.
  linear_map <- list(
    X = matrix(0, 4L, ncol(ctx$Vp)),
    trmind = seq_len(ncol(ctx$Vp))
  )
  w <- testthat::capture_warnings(
    pw <- compute_pointwise_ci("satterthwaite", 0.95, linear_map, ctx)
  )
  testthat::expect_length(w, 1L)
  testthat::expect_match(w, "interval limits are missing")
  testthat::expect_true(all(is.na(pw$crit * 1) & is.na(pw$df)))
  # plot.pffr() draws standard-error bands and never consumes coef()'s
  # interval limits; summary() drops non-finite df. Both must stay clean.
  pdf_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(pdf_file)
  on.exit(
    {
      grDevices::dev.off()
      unlink(pdf_file)
    },
    add = TRUE
  )
  testthat::expect_no_error(plot(fit, pages = 1))
  testthat::expect_no_error(print(summary(fit)))
  testthat::expect_null(pffr_summary_df(
    structure(
      list(pffr = list(sandwich_info = list(type = "none"))),
      class = "pffr"
    )
  ))
})

testthat::test_that("missing-response bookkeeping preserves group alignment", {
  meta <- list(
    nobs = 3L,
    nyindex = 4L,
    cluster = c("b", "a", "b"),
    missing_indices = c(2L, 9L),
    is_sparse = FALSE
  )
  testthat::expect_equal(
    build_cluster_id(meta),
    rep(meta$cluster, each = 4)[-c(2, 9)]
  )
  meta$missing_indices <- integer(0)
  testthat::expect_length(build_cluster_id(meta), 12L)
})

testthat::test_that("unsupported families cannot silently substitute HC", {
  set.seed(84102)
  b <- mgcv::gam(y ~ x, data = data.frame(y = rnorm(30), x = rnorm(30)))
  b$family$sandwich <- function(...) NULL
  testthat::expect_error(
    gam_sandwich_cluster(b, rep(1:10, each = 3)),
    "No cluster-robust"
  )
  testthat::expect_error(
    gam_sandwich_cluster_cl2(b, rep(1:10, each = 3)),
    "No cluster-robust"
  )
})

#--------------------------------------
# Round-2 review: the leverage/invariant diagnostics through the public API
#--------------------------------------
#
# The monitors were previously asserted only through refund:::pffr_vcov() and
# refund:::gam_sandwich_cluster_cl2(). These exercise the same behaviour the
# way a user meets it: pffr(sandwich = "cl2") at fit time and coef() after.

testthat::test_that("the shortcut cap warns once at fit time, the exact path not at all", {
  testthat::skip_on_cran()
  # The "influential" design saturates one cluster's leverage (one covariate
  # value is three orders of magnitude off), so the shortcut's H_gg hits the
  # cap while the exact path only floors residual-block eigenvalues.
  fixture <- make_exactcl2_fixture("poisson", 4L, "influential")
  fit_cl2 <- function(adjustment) {
    suppressMessages(refund::pffr(
      Y ~ xlin,
      data = fixture$data,
      yind = fixture$yind,
      family = fixture$family,
      bs.yindex = list(bs = "ps", k = fixture$k, m = c(2, 1)),
      sandwich = "cl2",
      cl2_adjustment = adjustment,
      cluster = fixture$cluster
    ))
  }
  # Count only cap warnings: the Poisson fit itself may warn about other
  # things, and this assertion is about the cap warning not repeating.
  w_short <- testthat::capture_warnings(fit_short <- fit_cl2("shortcut"))
  testthat::expect_length(grep("hit the leverage cap", w_short), 1L)
  testthat::expect_gt(fit_short$pffr$sandwich_info$n_capped, 0)

  w_exact <- testthat::capture_warnings(fit_exact <- fit_cl2("exact"))
  testthat::expect_length(grep("hit the leverage cap", w_exact), 0L)
  testthat::expect_null(fit_exact$pffr$sandwich_info$hat_invariant_violation)

  # The covariance is computed once, at fit time; reading it back through
  # coef() must not re-run the adjustment and warn a second time.
  testthat::expect_no_warning(
    coef(fit_short, ci = "pointwise", crit = "satterthwaite", n1 = 12)
  )
  testthat::expect_no_warning(
    coef(fit_exact, ci = "pointwise", crit = "satterthwaite", n1 = 12)
  )
})

testthat::test_that("sandwich_info of a benign cl2 fit carries the hat monitors", {
  testthat::skip_on_cran()
  fixture <- make_lb5_fixture(amp = 1, n_grid = 20L, k = 8L, 21L)
  fit <- suppressWarnings(suppressMessages(refund::pffr(
    Y ~ xlin,
    data = fixture$data,
    yind = fixture$yind,
    family = stats::poisson(),
    bs.yindex = list(bs = "ps", k = fixture$k, m = c(2, 1)),
    sandwich = "cl2"
  )))
  info <- fit$pffr$sandwich_info
  testthat::expect_identical(info$type, "cl2")
  for (slot in c("max_obs_leverage", "min_obs_leverage", "min_hat_eig")) {
    testthat::expect_true(
      is.finite(info[[slot]]),
      info = paste("sandwich_info slot", slot)
    )
  }
  # A benign fit respects the bounds, so nothing is flagged.
  testthat::expect_lte(info$max_obs_leverage, 1)
  testthat::expect_null(info$hat_invariant_violation)
})

testthat::test_that("a degenerate fit surfaces one invariant warning through pffr and coef", {
  testthat::skip_on_cran()
  fixture <- make_lb5_fixture(amp = 10, n_grid = 30L, k = 12L, 21L)
  w <- testthat::capture_warnings(
    fit <- suppressMessages(refund::pffr(
      Y ~ xlin,
      data = fixture$data,
      yind = fixture$yind,
      family = stats::poisson(),
      bs.yindex = list(bs = "ps", k = fixture$k, m = c(2, 1)),
      sandwich = "cl2"
    ))
  )
  testthat::expect_length(grep("NOT trustworthy", w), 1L)
  testthat::expect_type(
    fit$pffr$sandwich_info$hat_invariant_violation,
    "character"
  )
  testthat::expect_gt(fit$pffr$sandwich_info$max_obs_leverage, 1)
  # Reading the stored covariance back does not re-raise it.
  testthat::expect_length(
    grep(
      "NOT trustworthy",
      testthat::capture_warnings(coef(fit, ci = "none", n1 = 12))
    ),
    0L
  )
})

testthat::test_that("one coef() call warns at most once about undefined df", {
  # Both the smooth-term block (compute_pointwise_ci) and the parametric
  # block detect an undefined moment df independently; before round 2 each
  # warned on its own. Mock the df kernel so BOTH blocks see a non-finite df.
  set.seed(84105)
  G <- 10L
  D <- 8L
  tt <- seq(0, 1, length.out = D)
  dat <- list(Y = matrix(rnorm(G * D), G, D), x = rnorm(G))
  dat$Y <- dat$Y + outer(dat$x, sin(2 * pi * tt))
  fit <- suppressMessages(quiet_pffr(
    Y ~ x,
    yind = tt,
    data = dat,
    bs.yindex = list(bs = "ps", k = 5, m = c(2, 1)),
    bs.int = list(bs = "ps", k = 5, m = c(2, 1)),
    sandwich = "cl2"
  ))
  testthat::local_mocked_bindings(
    # One undefined contrast per block, the rest finite, so each block takes
    # the NA branch and every other interval limit stays usable.
    pffr_df_from_context = function(ctx, Xp, df_gram = NULL) {
      n <- nrow(Xp)
      c(NA_real_, rep(8, max(n - 1L, 0L)))[seq_len(n)]
    }
  )
  w <- testthat::capture_warnings(
    cf <- coef(fit, ci = "pointwise", crit = "satterthwaite", n1 = 12)
  )
  testthat::expect_length(grep("Undefined working-model moment df", w), 1L)

  # ... and both blocks really did produce missing limits.
  smooth_lims <- cf$smterms[[1]]$coef
  testthat::expect_true(any(is.na(smooth_lims[, "lower"])))
  testthat::expect_true(any(is.finite(smooth_lims[, "lower"])))
  testthat::expect_true(any(is.na(cf$pterms[, "lower"])))
})
