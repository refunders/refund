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
    sandwich = TRUE,
    cluster = subject,
    ci = "pointwise",
    n1 = 18
  )
  testthat::expect_identical(co$ci_meta$crit_used, "satterthwaite")
  testthat::expect_equal(co$smterms, explicit$smterms)
  # By-curve clustering of the same fit differs.
  by_curve <- coef(fit, cluster = seq_len(n), ci = "pointwise", n1 = 18)
  testthat::expect_false(isTRUE(all.equal(
    by_curve$smterms[[1]]$coef$se,
    co$smterms[[1]]$coef$se
  )))
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
  # The parametric df are the moment df of unit contrasts at subject level.
  smind <- unlist(lapply(fit$smooth, \(s) seq(s$first.para, s$last.para)))
  pind <- seq_along(fit$coefficients)[-smind]
  Xp_p <- diag(length(fit$coefficients))[pind, , drop = FALSE]
  testthat::expect_equal(
    unname(co$pterms[, "df"]),
    pffr_influence_df(pffr_influence(fit), Xp_p)$df,
    tolerance = 1e-10
  )
  testthat::expect_true(all(co$pterms[, "df"] <= G))
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
})

testthat::test_that("plot and summary of a CL2 fit stay clean", {
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
    bs.int = list(bs = "ps", k = 5, m = c(2, 1))
  ))
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
  out <- utils::capture.output(print(summary(fit)))
  testthat::expect_true(any(grepl("model-based, which is invalid", out)))
  testthat::expect_true(summary(fit)$sandwich)
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
  # Sparse responses: rows with a missing value are not in the model frame.
  sparse <- list(
    is_sparse = TRUE,
    ydata = data.frame(.obs = c(1, 1, 2, 2, 3), .value = c(1, NA, 2, 3, NA))
  )
  testthat::expect_equal(build_cluster_id(sparse), c(1, 2, 2))
})

testthat::test_that("unsupported families cannot silently substitute HC", {
  set.seed(84102)
  b <- mgcv::gam(y ~ x, data = data.frame(y = rnorm(30), x = rnorm(30)))
  b$family$sandwich <- function(...) NULL
  testthat::expect_error(
    gam_sandwich_cluster_cl2(b, rep(1:10, each = 3)),
    "No cluster-robust"
  )
})

#--------------------------------------
# The leverage/invariant diagnostics through the public API
#--------------------------------------

testthat::test_that("the exact leverage floor does not warn at fit time or in coef", {
  testthat::skip_on_cran()
  # The "influential" design saturates one cluster's leverage (one covariate
  # value is three orders of magnitude off); the exact path only floors
  # residual-block eigenvalues.
  fixture <- make_exactcl2_fixture("poisson", 4L, "influential")
  w <- testthat::capture_warnings(
    fit <- suppressMessages(refund::pffr(
      Y ~ xlin,
      data = fixture$data,
      yind = fixture$yind,
      family = fixture$family,
      bs.yindex = list(bs = "ps", k = fixture$k, m = c(2, 1)),
      cluster = fixture$cluster
    ))
  )
  testthat::expect_length(grep("leverage|trustworthy", w), 0L)
  testthat::expect_gt(fit$pffr$sandwich_info$n_adjusted, 0)
  testthat::expect_null(fit$pffr$sandwich_info$hat_invariant_violation)
  testthat::expect_no_warning(suppressMessages(coef(
    fit,
    ci = "pointwise",
    n1 = 12
  )))
})

testthat::test_that("sandwich_info of a benign cl2 fit carries the hat monitors", {
  testthat::skip_on_cran()
  fixture <- make_lb5_fixture(amp = 1, n_grid = 20L, k = 8L, 21L)
  fit <- suppressWarnings(suppressMessages(refund::pffr(
    Y ~ xlin,
    data = fixture$data,
    yind = fixture$yind,
    family = stats::poisson(),
    bs.yindex = list(bs = "ps", k = fixture$k, m = c(2, 1))
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
      bs.yindex = list(bs = "ps", k = fixture$k, m = c(2, 1))
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

testthat::test_that("one coef() call notes at most once about undefined df", {
  # Both the smooth-term block and the parametric block can meet an undefined
  # moment df. Mock the df kernel so BOTH blocks see a non-finite df.
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
    bs.int = list(bs = "ps", k = 5, m = c(2, 1))
  ))
  testthat::local_mocked_bindings(
    # One undefined contrast per block, the rest finite.
    pffr_influence_df = function(core, Xp, chunk_size = 32L) {
      n <- nrow(Xp)
      list(df = c(NA_real_, rep(8, max(n - 1L, 0L)))[seq_len(n)], G = G)
    }
  )
  msgs <- testthat::capture_messages(
    cf <- coef(fit, ci = "pointwise", n1 = 12)
  )
  testthat::expect_length(grep("Satterthwaite df are undefined", msgs), 1L)
  # ... and both blocks fell back to the Gaussian critical value there.
  smooth_tab <- cf$smterms[[1]]$coef
  testthat::expect_identical(smooth_tab$df[1], Inf)
  testthat::expect_equal(
    smooth_tab$upper[1] - smooth_tab$value[1],
    stats::qnorm(0.975) * smooth_tab$se[1]
  )
  testthat::expect_true(all(smooth_tab$df[-1] == 8))
  testthat::expect_identical(unname(cf$pterms[1, "df"]), Inf)
})
