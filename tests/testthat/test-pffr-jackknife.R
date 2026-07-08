#--------------------------------------
# Tests for the leave-one-cluster-out (LOCO) jackknife SE (X3-SHIP)
#
# pffr_jackknife_se() / pffr_jackknife_core(): exact SMW downdate of the
# penalized normal equations for Gaussian-identity fits (1e-10 vs a direct
# deleted-cluster solve); a one-step approximation at fixed (lambda, weights)
# for GLM families (validated against a direct fixed-sp deleted-cluster IRLS
# refit); the crit / se_scale API; and the predict.pffr opt-in.
#
# Provenance of the estimator + its coverage: paper sec-eymean and
# notes/X3-phase2-findings.md (jackknife 0.856 (z) / 0.864 (t_{G-1}) at
# AR(1), G=20 vs 0.231 for the plug-in cluster sandwich).
#--------------------------------------

# Shared small fixtures (built once for this file).
make_jack_fits <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      dat <- pffr_simulate(
        Y ~ ff(X1),
        n = 25,
        nxgrid = 15,
        nygrid = 15,
        SNR = 5,
        effects = list(X1 = "random"),
        intercept = "random",
        seed = 5150
      )
      yind <- attr(dat, "yindex")
      fit_gauss <- suppressWarnings(pffr(
        Y ~ ff(X1),
        data = dat,
        yind = yind,
        sandwich = "cluster"
      ))
      cache <<- list(dat = dat, yind = yind, fit_gauss = fit_gauss)
    }
    cache
  }
})

make_poisson_fit <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      set.seed(202)
      dat <- pffr_simulate(
        Y ~ ff(X1),
        n = 30,
        nxgrid = 12,
        nygrid = 12,
        SNR = 8,
        effects = list(X1 = "random"),
        intercept = "random"
      )
      lam <- exp(0.4 * scale(dat$Y))
      dat$Ycount <- matrix(
        rpois(length(lam), pmax(lam, 0.05)),
        nrow = nrow(dat$Y)
      )
      yind <- attr(dat, "yindex")
      fit <- suppressWarnings(pffr(
        Ycount ~ ff(X1),
        data = dat,
        yind = yind,
        family = poisson(),
        sandwich = "cluster"
      ))
      cache <<- list(dat = dat, yind = yind, fit = fit)
    }
    cache
  }
})

test_that("SMW downdate matches the direct deleted-cluster solve to 1e-10", {
  skip_on_cran()
  fit <- make_jack_fits()$fit_gauss
  core <- refund:::pffr_jackknife_core(fit, smw_check = TRUE)
  expect_true(is.finite(core$smw_max_abs_err))
  expect_lt(core$smw_max_abs_err, 1e-10)
})

test_that("Gaussian jackknife SEs are finite, positive and correctly sized", {
  skip_on_cran()
  f <- make_jack_fits()
  jk <- pffr_jackknife_se(f$fit_gauss)
  # one row per fitted evaluation point (curve-major, index-fastest)
  expect_equal(nrow(jk), 25L * 15L)
  expect_true(all(is.finite(jk$se)))
  expect_true(all(jk$se > 0))
  expect_true(all(c("fit", "se", "lower", "upper") %in% names(jk)))
  expect_true(all(jk$lower < jk$upper))
})

test_that("jackknife diagnostics attributes are present and sane", {
  skip_on_cran()
  jk <- pffr_jackknife_se(make_jack_fits()$fit_gauss)
  expect_identical(attr(jk, "G"), 25L)
  expect_length(attr(jk, "max_eig_Hgg"), 25L)
  expect_true(is.numeric(attr(jk, "max_eig_Hgg")))
  expect_true(all(attr(jk, "max_eig_Hgg") <= 1 + 1e-8))
  expect_true(
    is.numeric(attr(jk, "n_floored")) ||
      is.integer(attr(jk, "n_floored"))
  )
  expect_identical(attr(jk, "crit"), "tG1")
  expect_identical(attr(jk, "df"), 24)
  expect_identical(attr(jk, "se_scale"), "link")
})

test_that("crit = 'tG1' gives wider intervals than crit = 'z'", {
  skip_on_cran()
  fit <- make_jack_fits()$fit_gauss
  jk_t <- pffr_jackknife_se(fit, crit = "tG1")
  jk_z <- pffr_jackknife_se(fit, crit = "z")
  # same SEs, wider critical value under t_{G-1}
  expect_equal(jk_t$se, jk_z$se, tolerance = 1e-12)
  expect_gt(attr(jk_t, "crit_value"), attr(jk_z, "crit_value"))
  expect_equal(
    attr(jk_t, "crit_value"),
    qt(0.975, df = 24),
    tolerance = 1e-12
  )
  expect_equal(attr(jk_z, "crit_value"), qnorm(0.975), tolerance = 1e-12)
  # t intervals strictly wider
  expect_true(all((jk_t$upper - jk_t$lower) > (jk_z$upper - jk_z$lower)))
})

test_that("se_scale = 'response' applies the delta method via mu.eta", {
  skip_on_cran()
  fit <- make_poisson_fit()$fit
  jk_link <- pffr_jackknife_se(fit, se_scale = "link")
  jk_resp <- pffr_jackknife_se(fit, se_scale = "response")
  expect_true(all(is.finite(jk_resp$se)))
  expect_true(all(jk_resp$se > 0))
  # response-scale fit is on the mean (mu) scale, link-scale fit is eta
  expect_equal(
    jk_resp$fit,
    as.numeric(fit$family$linkinv(jk_link$fit)),
    tolerance = 1e-8
  )
  # response SE = |mu.eta(eta)| * link SE (delta method)
  expect_equal(
    jk_resp$se,
    jk_link$se * abs(as.numeric(fit$family$mu.eta(jk_link$fit))),
    tolerance = 1e-8
  )
  expect_identical(attr(jk_resp, "se_scale"), "response")
})

test_that("GLM one-step matches a direct fixed-sp deleted-cluster refit", {
  skip_on_cran()
  fit <- make_poisson_fit()$fit
  fam <- fit$family

  cluster_id <- refund:::build_cluster_id(fit$pffr)
  Vp <- refund:::pffr_canonicalize_cov(fit)$model$Vp
  sig2 <- fit$sig2
  if (is.null(sig2) || !is.finite(sig2) || sig2 <= 0) sig2 <- 1
  A_inv <- Vp / sig2
  theta <- fit$coefficients
  Xfull <- predict(fit, type = "lpmatrix", reformat = FALSE)
  mi <- fit$pffr$missing_indices
  Xtr <- if (!is.null(mi)) Xfull[-mi, , drop = FALSE] else Xfull
  y <- as.vector(fit$y)
  mu <- as.vector(fit$fitted.values)
  eta <- as.vector(fit$linear.predictors)
  pw <- fit$prior.weights
  if (is.null(pw)) pw <- rep(1, length(y))
  mu_eta <- fam$mu.eta(eta)
  var_mu <- fam$variance(mu)
  s <- sqrt(pw * mu_eta^2 / var_mu)
  s[!is.finite(s)] <- 0
  Xt <- Xtr * s
  rt <- sign(mu_eta) * sqrt(pw / var_mu) * (y - mu)
  A <- solve(A_inv)
  S_lam <- A - crossprod(Xt) # fixed penalty (constrained coefficient space)

  groups <- unique(cluster_id)
  test_groups <- groups[seq_len(min(4L, length(groups)))]
  gaps <- vapply(
    test_groups,
    function(g) {
      idx <- which(cluster_id == g)
      Xtg <- Xt[idx, , drop = FALSE]
      rtg <- rt[idx]
      Hgg <- Xtg %*% A_inv %*% t(Xtg)
      Hgg <- 0.5 * (Hgg + t(Hgg))
      ee <- eigen(diag(length(idx)) - Hgg, symmetric = TRUE)
      vals <- pmax(ee$values, 1e-8)
      u <- ee$vectors %*% (crossprod(ee$vectors, rtg) / vals)
      delta <- as.vector(A_inv %*% crossprod(Xtg, u))
      theta_onestep <- theta - delta

      keep <- setdiff(seq_along(y), idx)
      Xk <- Xtr[keep, , drop = FALSE]
      yk <- y[keep]
      pwk <- pw[keep]
      th <- theta
      for (it in 1:200) {
        eta_it <- as.vector(Xk %*% th)
        mu_it <- fam$linkinv(eta_it)
        mue_it <- fam$mu.eta(eta_it)
        var_it <- fam$variance(mu_it)
        w_it <- pwk * mue_it^2 / var_it
        z_it <- eta_it + (yk - mu_it) / mue_it
        XtW <- t(Xk * w_it)
        th_new <- as.vector(solve(XtW %*% Xk + S_lam, XtW %*% z_it))
        if (max(abs(th_new - th)) < 1e-13) {
          th <- th_new
          break
        }
        th <- th_new
      }
      # RELATIVE gap: one-step error relative to the deleted-cluster
      # coefficient scale, so the criterion documents the one-step
      # approximation error itself rather than an absolute magic constant
      # vulnerable to RNG/BLAS drift.
      max(abs(theta_onestep - th)) / max(abs(th))
    },
    numeric(1)
  )
  # Measured relative one-step gap on this fit is ~2e-3 (0.2% of the
  # coefficient scale; see notes). Assert 1% relative -- holds with a ~5x
  # margin. The gap is the one-step (fixed W, lambda) approximation error, NOT
  # a bug: it vanishes for Gaussian-identity.
  expect_lt(max(gaps), 1e-2)
})

test_that("predict.pffr opt-in overrides SE but not fit; default unchanged", {
  skip_on_cran()
  fit <- make_jack_fits()$fit_gauss

  p_norm <- predict(fit, se.fit = TRUE, reformat = FALSE)
  p_jack <- predict(
    fit,
    se.fit = TRUE,
    reformat = FALSE,
    se_method = "jackknife"
  )

  # point predictions are exactly predict.gam's, unchanged
  expect_equal(
    as.numeric(p_jack$fit),
    as.numeric(p_norm$fit),
    tolerance = 1e-12
  )
  # SEs are replaced (differ from the sandwich SEs)
  expect_gt(max(abs(p_jack$se.fit - p_norm$se.fit)), 0)
  # and equal the standalone jackknife SEs
  jk <- pffr_jackknife_se(fit)
  expect_equal(as.numeric(p_jack$se.fit), jk$se, tolerance = 1e-10)

  # DEFAULT path (se_method = "normal") is byte-for-byte the pre-feature result:
  # se.fit comes from the fit-time sandwich covariance via pffr_vcov()
  X <- predict(fit, type = "lpmatrix", reformat = FALSE)
  se_sandwich <- sqrt(rowSums((X %*% fit$pffr$Vsandwich) * X))
  expect_equal(
    as.numeric(p_norm$se.fit),
    as.numeric(se_sandwich),
    tolerance = 1e-8
  )

  # reformat = TRUE (matrix output) also honors the opt-in
  p_jack_mat <- predict(fit, se.fit = TRUE, se_method = "jackknife")
  expect_true(is.matrix(p_jack_mat$se.fit))
  expect_equal(dim(p_jack_mat$se.fit), c(25L, 15L))
})

test_that("predict jackknife rejects type = 'terms'", {
  skip_on_cran()
  fit <- make_jack_fits()$fit_gauss
  expect_error(
    predict(fit, se.fit = TRUE, type = "terms", se_method = "jackknife"),
    "only supports type"
  )
})

test_that("unsupported families and missing cluster structure error clearly", {
  skip_on_cran()
  # gaulss (two linear predictors, no top-level mu.eta) is unsupported
  gm <- get_gaulss_model()
  expect_error(
    pffr_jackknife_se(gm),
    "not supported"
  )
  # a single-cluster resolution errors informatively
  fit <- make_jack_fits()$fit_gauss
  expect_error(
    pffr_jackknife_se(fit, cluster = rep(1L, fit$pffr$nobs)),
    "at least two independent"
  )
})

test_that("dense fits with missing responses work (no double row removal)", {
  skip_on_cran()
  # An NA response makes predict(type = "lpmatrix") return the FITTED
  # (missing-removed) rows already: nrow == length(fit$y) == 159 here, not the
  # full 20*8 = 160 grid. Subsetting missing_indices out again removed a fitted
  # row and made every NA-response fit error (council review, MUST-FIX 1).
  # (n = 20, not 10: the n = 10 version of this fixture has an ill-conditioned
  # penalized bread, cond(A_inv) ~ 1e8, which pushes the double-inversion
  # direct check's roundoff to ~1.5e-10; at n = 20 cond ~ 2e3 and the SMW
  # check sits at ~2e-14.)
  set.seed(11)
  dat <- pffr_simulate(
    Y ~ ff(X1),
    n = 20,
    nxgrid = 12,
    nygrid = 8,
    SNR = 5,
    effects = list(X1 = "random"),
    intercept = "random"
  )
  yind <- attr(dat, "yindex")
  dat$Y[3, 5] <- NA # missing point 21 in curve-major order: (3-1)*8 + 5
  fit <- suppressWarnings(suppressMessages(
    pffr(Y ~ ff(X1), data = dat, yind = yind, sandwich = "cluster")
  ))
  expect_identical(fit$pffr$missing_indices, 21L)
  expect_length(fit$y, 159L)

  # standalone path: one row per FITTED point, SMW still exact
  jk <- pffr_jackknife_se(fit)
  expect_equal(nrow(jk), 159L)
  expect_true(all(is.finite(jk$se)))
  expect_true(all(jk$se > 0))
  core <- refund:::pffr_jackknife_core(fit, smw_check = TRUE)
  expect_lt(core$smw_max_abs_err, 1e-10)

  # predict opt-in: reformatting pads the missing point with NA, and the
  # non-missing entries match the standalone jackknife SEs
  p <- predict(fit, se.fit = TRUE, se_method = "jackknife")
  expect_equal(dim(p$se.fit), c(20L, 8L))
  expect_true(is.na(p$se.fit[3, 5]))
  expect_identical(sum(is.na(p$se.fit)), 1L)
  se_vec <- as.vector(t(p$se.fit))
  expect_equal(se_vec[-21], jk$se, tolerance = 1e-10)
})

test_that("offset fits: exact at fitted points, newdata rejected", {
  skip_on_cran()
  set.seed(12)
  dat <- pffr_simulate(
    Y ~ ff(X1),
    n = 15,
    nxgrid = 12,
    nygrid = 10,
    SNR = 5,
    effects = list(X1 = "random"),
    intercept = "random"
  )
  yind <- attr(dat, "yindex")
  off_mat <- matrix(
    0.3 * sin(seq(0, 2, length.out = 10)),
    nrow = 15,
    ncol = 10,
    byrow = TRUE
  )
  fit_off <- suppressWarnings(suppressMessages(
    pffr(
      Y ~ ff(X1),
      data = dat,
      yind = yind,
      sandwich = "cluster",
      offset = off_mat
    )
  ))
  expect_true(any(fit_off$offset != 0))

  # returned fit is the offset-INCLUSIVE linear predictor (council MUST-FIX 2)
  jk <- pffr_jackknife_se(fit_off)
  expect_equal(
    jk$fit,
    as.vector(fit_off$linear.predictors),
    tolerance = 1e-10
  )
  # ... which differs from the offset-free X %*% theta by exactly the offset
  X <- predict(fit_off, type = "lpmatrix", reformat = FALSE)
  expect_equal(
    jk$fit - as.vector(X %*% fit_off$coefficients),
    as.vector(fit_off$offset),
    tolerance = 1e-10
  )
  # response-scale delta method evaluates linkinv/mu.eta at the CORRECT eta
  jk_resp <- pffr_jackknife_se(fit_off, se_scale = "response")
  expect_equal(
    jk_resp$fit,
    as.numeric(fit_off$family$linkinv(jk$fit)),
    tolerance = 1e-10
  )
  expect_equal(
    jk_resp$se,
    jk$se * abs(as.numeric(fit_off$family$mu.eta(jk$fit))),
    tolerance = 1e-10
  )
  # newdata + offset fit: rejected with a clear message (option: reject)
  expect_error(
    pffr_jackknife_se(fit_off, newdata = dat[1:2, ]),
    "model offset"
  )
})

test_that("non-unit prior weights: whitening keeps the SMW downdate exact", {
  skip_on_cran()
  set.seed(13)
  dat <- pffr_simulate(
    Y ~ ff(X1),
    n = 15,
    nxgrid = 12,
    nygrid = 10,
    SNR = 5,
    effects = list(X1 = "random"),
    intercept = "random"
  )
  yind <- attr(dat, "yindex")
  w_curve <- runif(15, 0.5, 2)
  W <- matrix(w_curve, nrow = 15, ncol = 10) # constant per curve
  fit_w <- suppressWarnings(suppressMessages(
    pffr(Y ~ ff(X1), data = dat, yind = yind, sandwich = "cluster", weights = W)
  ))
  expect_true(any(fit_w$prior.weights != 1))
  # Gaussian-identity with fixed prior weights: the weighted SMW downdate is
  # still EXACT (W does not depend on mu), so 1e-10 must hold
  core <- refund:::pffr_jackknife_core(fit_w, smw_check = TRUE)
  expect_lt(core$smw_max_abs_err, 1e-10)
  jk <- pffr_jackknife_se(fit_w)
  expect_true(all(is.finite(jk$se)))
  expect_true(all(jk$se > 0))
  # weights change the answer relative to the unweighted fit machinery
  fit_u <- suppressWarnings(suppressMessages(
    pffr(Y ~ ff(X1), data = dat, yind = yind, sandwich = "cluster")
  ))
  expect_false(isTRUE(all.equal(jk$se, pffr_jackknife_se(fit_u)$se)))
})

test_that("coarser cluster= grouping: variance matches direct LOO solves", {
  skip_on_cran()
  fit <- make_jack_fits()$fit_gauss # 25 curves x 15 points, gaussian identity
  grp <- rep(1:5, each = 5) # 5 superclusters of 5 curves

  jk <- pffr_jackknife_se(fit, cluster = grp, crit = "z")
  expect_identical(attr(jk, "G"), 5L)
  expect_identical(attr(jk, "df"), NA_real_)

  # Independent re-derivation: direct deleted-SUPERCLUSTER penalized solves
  # (Gaussian identity: (A - Xg'Xg) theta_(-g) = A theta - Xg' y_g), then the
  # mean-centered CV3 variance -- no SMW, no shared code path beyond inputs.
  Vp <- fit$Vp
  sig2 <- fit$sig2
  A <- solve(Vp / sig2)
  theta <- fit$coefficients
  X <- predict(fit, type = "lpmatrix", reformat = FALSE)
  y <- as.vector(fit$y)
  cid <- refund:::build_cluster_id(fit$pffr, cluster = grp)
  G <- length(unique(cid))
  eta_loo <- matrix(NA_real_, nrow = G, ncol = nrow(X))
  for (g in seq_len(G)) {
    idx <- which(cid == g)
    Xg <- X[idx, , drop = FALSE]
    theta_g <- solve(
      A - crossprod(Xg),
      as.vector(A %*% theta) - as.vector(crossprod(Xg, y[idx]))
    )
    eta_loo[g, ] <- as.vector(X %*% theta_g)
  }
  var_direct <- (G - 1) /
    G *
    colSums(sweep(eta_loo, 2L, colMeans(eta_loo))^2)
  expect_equal(jk$se, sqrt(var_direct), tolerance = 1e-8)

  # the coarser grouping gives genuinely different SEs than by-curve
  expect_false(isTRUE(all.equal(jk$se, pffr_jackknife_se(fit, crit = "z")$se)))
})

test_that("predict jackknife rejects the (parallel) cluster dot", {
  skip_on_cran()
  fit <- make_jack_fits()$fit_gauss
  expect_error(
    predict(fit, se.fit = TRUE, se_method = "jackknife", cluster = 1),
    "does not accept"
  )
})
