# S1 storage contract: $Vp/$Vc/$Ve are model-based ALWAYS; the robust
# covariance lives in $pffr$Vsandwich (+ $pffr$sandwich_info metadata) and all
# refund consumers resolve covariances through pffr_vcov(). These tests are the
# brief's (a)-(e) acceptance tests for the redesign that eliminated the
# double-sandwich footgun (recomputing a sandwich on a fit whose $Vp had been
# overwritten with the robust covariance double-applied the correction, up to
# 284x SE inflation; see notes/dti-discrepancy-reconciliation.md in pffr-ci).

# Build the shared fixtures once for this file.
make_storage_fits <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      dat <- pffr_simulate(
        Y ~ ff(X1),
        n = 25,
        nxgrid = 20,
        nygrid = 20,
        SNR = 5,
        effects = list(X1 = "random"),
        intercept = "random",
        seed = 5150
      )
      yind <- attr(dat, "yindex")
      cache <<- list(
        fit_none = pffr(
          Y ~ ff(X1),
          data = dat,
          yind = yind,
          sandwich = FALSE
        ),
        fit_cluster = quiet_pffr(
          Y ~ ff(X1),
          data = dat,
          yind = yind
        )
      )
    }
    cache
  }
})

# Devolve a current-format fit to the old (format-1) layout: robust matrices
# in $Vp/$Vc/$Ve, model-based stashed in $pffr$model_cov.
make_old_format <- function(fit) {
  stopifnot(!is.null(fit$pffr$Vsandwich))
  old <- fit
  old$pffr$model_cov <- list(Vp = fit$Vp, Ve = fit$Ve, Vc = fit$Vc)
  old$Vp <- fit$pffr$Vsandwich
  old$Vc <- fit$pffr$Vsandwich
  old$Ve <- fit$pffr$Vsandwich
  old$pffr$Vsandwich <- NULL
  old$pffr$sandwich_info <- NULL
  old$pffr$Vsandwich_cache <- NULL
  old$pffr$cov_format <- NULL
  old
}

reset_pffr_warn_once <- function() {
  st <- refund:::.pffr_state
  rm(list = ls(envir = st), envir = st)
}

ff_se <- function(fit, ...) {
  coef(fit, n1 = 20, n2 = 8, ...)$smterms[[2]]$coef$se
}

test_that("(a) $Vp/$Vc/$Ve on a sandwich fit are bitwise the model-based ones", {
  skip_on_cran()

  fits <- make_storage_fits()
  fit_none <- fits$fit_none
  fit_cluster <- fits$fit_cluster

  # bitwise identical to the sandwich = FALSE fit on the same data/seed
  expect_identical(fit_cluster$Vp, fit_none$Vp)
  expect_identical(fit_cluster$Vc, fit_none$Vc)
  expect_identical(fit_cluster$Ve, fit_none$Ve)

  # the robust covariance is stored in its own slot, with metadata
  expect_true(is.matrix(fit_cluster$pffr$Vsandwich))
  expect_gt(max(abs(fit_cluster$pffr$Vsandwich - fit_cluster$Vp)), 0)
  info <- fit_cluster$pffr$sandwich_info
  expect_identical(info$type, "cl2")
  expect_identical(info$cl2_adjustment, "exact")
  expect_identical(info$G, 25L)
  expect_true(is.numeric(info$n_adjusted))
  expect_true(!is.null(info$version))
  # no legacy stash on new fits
  expect_null(fit_cluster$pffr$model_cov)
})

test_that("(b) coef(fit, sandwich = FALSE) equals the sandwich-free fit", {
  skip_on_cran()

  fits <- make_storage_fits()
  # this exact call silently returned robust SEs pre-fix
  se_robustfit_none <- ff_se(fits$fit_cluster, sandwich = FALSE)
  se_cleanfit_none <- ff_se(fits$fit_none, sandwich = FALSE)
  expect_equal(se_robustfit_none, se_cleanfit_none, tolerance = 1e-12)

  # and they must differ from the fit-time robust SEs
  se_robust <- ff_se(fits$fit_cluster)
  expect_false(isTRUE(all.equal(se_robustfit_none, se_robust)))
})

test_that("(c) recompute-on-robust-fit equals compute-on-clean-fit", {
  skip_on_cran()

  fits <- make_storage_fits()
  # kills the double-sandwich class of bug forever
  se_cl2_on_robust <- ff_se(fits$fit_cluster, sandwich = TRUE)
  se_cl2_on_clean <- ff_se(fits$fit_none, sandwich = TRUE)
  expect_equal(se_cl2_on_robust, se_cl2_on_clean, tolerance = 1e-10)

  # same for the pterm SEs
  co_robust <- coef(fits$fit_cluster, sandwich = TRUE, n1 = 20, n2 = 8)
  co_clean <- coef(fits$fit_none, sandwich = TRUE, n1 = 20, n2 = 8)
  expect_equal(
    co_robust$pterms[, "se"],
    co_clean$pterms[, "se"],
    tolerance = 1e-10
  )

  # recomputation caches: repeated identical requests reuse the cached matrix
  V1 <- refund:::pffr_vcov(fits$fit_none, sandwich = TRUE)
  expect_true(!is.null(fits$fit_none$pffr$Vsandwich_cache[["cl2"]]))
  V2 <- refund:::pffr_vcov(fits$fit_none, sandwich = TRUE)
  expect_identical(V1, V2)
})

test_that("(d) predict and plot honor the fit-time sandwich type", {
  skip_on_cran()

  fits <- make_storage_fits()
  fit_none <- fits$fit_none
  fit_cluster <- fits$fit_cluster

  # predict(se.fit = TRUE): robust fit serves SEs from $pffr$Vsandwich ...
  p_cl <- predict(fit_cluster, se.fit = TRUE, reformat = FALSE)
  X <- predict(fit_cluster, type = "lpmatrix", reformat = FALSE)
  se_manual <- sqrt(rowSums((X %*% fit_cluster$pffr$Vsandwich) * X))
  expect_equal(as.numeric(p_cl$se.fit), as.numeric(se_manual), tolerance = 1e-8)

  # ... while the clean fit serves model-based SEs (mgcv's own Vp), unchanged
  p_none <- predict(fit_none, se.fit = TRUE, reformat = FALSE)
  se_manual_none <- sqrt(rowSums((X %*% fit_none$Vp) * X))
  expect_equal(
    as.numeric(p_none$se.fit),
    as.numeric(se_manual_none),
    tolerance = 1e-8
  )
  expect_gt(max(abs(p_cl$se.fit - p_none$se.fit)), 0)

  # the caller's fit object is not modified by predict's injection
  expect_identical(fit_cluster$Vp, fit_none$Vp)

  # plot.pffr: standard-error bands come from the resolved covariance
  # (compare the 1-D intercept smooth s(yindex.vec), plot data element 1)
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  pl_cl <- plot(fit_cluster, pages = 1, se = TRUE)
  pl_none <- plot(fit_none, pages = 1, se = TRUE)
  se_pl_cl <- pl_cl[[1]]$se
  se_pl_none <- pl_none[[1]]$se
  expect_true(length(se_pl_cl) > 0)
  expect_gt(max(abs(se_pl_cl - se_pl_none)), 0)
})

test_that("(e) old-format objects load with a warning and correct numbers", {
  skip_on_cran()

  fits <- make_storage_fits()
  fit_none <- fits$fit_none
  fit_cluster <- fits$fit_cluster
  fit_old <- make_old_format(fit_cluster)

  # one-time warning on first read; correct model-based numbers from the stash
  reset_pffr_warn_once()
  expect_warning(
    se_none_old <- ff_se(fit_old, sandwich = FALSE),
    "older refund version"
  )
  expect_equal(
    se_none_old,
    ff_se(fit_none, sandwich = FALSE),
    tolerance = 1e-10
  )

  # second read in the same session: no repeated warning; the CL2 sandwich
  # is recomputed from the stashed model-based bread
  expect_no_warning(se_cl_old <- ff_se(fit_old))
  expect_equal(se_cl_old, ff_se(fit_cluster), tolerance = 1e-10)
  co_old <- coef(fit_old, ci = "pointwise", n1 = 20, n2 = 8)
  expect_identical(co_old$ci_meta$crit_used, "satterthwaite")

  # objects of refund <= 0.1-40 (no stash at all) warn that recovery is
  # impossible and use their stored covariance with Gaussian critical values
  fit_ancient <- fit_old
  fit_ancient$pffr$model_cov <- NULL
  reset_pffr_warn_once()
  expect_warning(
    co_anc <- coef(fit_ancient, ci = "pointwise", n1 = 20, n2 = 8),
    "refund 0.1-40 or earlier"
  )
  expect_identical(co_anc$ci_meta$crit_used, "z")
  expect_equal(
    co_anc$smterms[[2]]$coef$se,
    ff_se(fit_cluster),
    tolerance = 1e-10
  )
  reset_pffr_warn_once()
})
