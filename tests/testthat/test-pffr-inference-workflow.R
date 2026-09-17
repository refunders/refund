testthat::test_that("subject clustering and covariance survive coef predict plot", {
  set.seed(84101)
  G <- 20L
  D <- 18L
  subject <- rep(seq_len(G), times = rep(c(1L, 2L), length.out = G))
  n <- length(subject)
  tt <- seq(0, 1, length.out = D)
  dat <- list(Y = matrix(rnorm(n * D), n, D), x = rnorm(n), subject = subject)
  dat$Y <- dat$Y + outer(dat$x, sin(2 * pi * tt))
  fit <- suppressMessages(refund::pffr(
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
