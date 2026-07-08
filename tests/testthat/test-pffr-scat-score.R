#--------------------------------------
# S4: exact scaled-t (scat) sandwich score
#--------------------------------------
#
# scat() is a scaled-t LOCATION family, not an exponential family, so the
# generic working-residual score is only an approximation. These tests validate
# the exact score dl/dmu = (nu+1) r / (nu sig^2 + r^2), its Gaussian limit, the
# EDF-consistency of the Fisher-whitened hat, and the disclosure warning that
# fires for families with neither an exact nor a two-block score.

# scaled-t log-density (constant terms retained; irrelevant for d/dmu)
scat_logdens <- function(y, mu, sig, nu) {
  lgamma((nu + 1) / 2) -
    lgamma(nu / 2) -
    log(sig * sqrt(pi * nu)) -
    ((nu + 1) / 2) * log1p(((y - mu) / sig)^2 / nu)
}
scat_score <- function(r, nu, sig) (nu + 1) * r / (nu * sig^2 + r^2)

# inline central difference (fallback when numDeriv is unavailable)
central_diff <- function(f, x, h = 1e-5) (f(x + h) - f(x - h)) / (2 * h)

# capture warning messages emitted while evaluating `expr`
captured_warnings <- function(expr) {
  msgs <- character(0)
  withCallingHandlers(
    force(expr),
    warning = function(cnd) {
      msgs <<- c(msgs, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    }
  )
  msgs
}

approx_msg <- "use the exponential-family working-residual approximation"

# Reset ONLY the approx-score warn-once keys (not the whole .pffr_state env,
# which also carries unrelated session flags such as the old-format back-compat
# notice used by other test files in the same process).
clear_approx_score_warnings <- function() {
  st <- refund:::.pffr_state
  keys <- grep("^approx_score_", ls(envir = st), value = TRUE)
  if (length(keys) > 0) {
    rm(list = keys, envir = st)
  }
  invisible(NULL)
}

# Shared scat fixture for the fit-based tests (heavy-tailed contamination so
# scat estimates a genuine finite nu); computed once per test session.
scat_test_env <- new.env(parent = emptyenv())
get_scat_test_model <- function() {
  if (is.null(scat_test_env$model)) {
    set.seed(21001)
    dat <- pffr_simulate(
      Y ~ xlin,
      n = 45,
      nygrid = 20,
      effects = list(xlin = "dnorm"),
      SNR = 6
    )
    t <- attr(dat, "yindex")
    dat$Y <- dat$Y +
      matrix(rt(length(dat$Y), df = 3) * 0.15, nrow = nrow(dat$Y))
    scat_test_env$model <- pffr(
      Y ~ xlin,
      yind = t,
      data = dat,
      family = mgcv::scat(),
      bs.yindex = list(bs = "ps", k = 12, m = c(2, 1)),
      sandwich = "none"
    )
  }
  scat_test_env$model
}


test_that("(i) scat score dl/dmu matches the numerical derivative to 1e-6", {
  have_nd <- requireNamespace("numDeriv", quietly = TRUE)
  max_err <- 0
  for (nu in c(3.5, 6, 15, 50, 300)) {
    for (sig in c(0.4, 1, 3)) {
      for (r in seq(-8, 8, by = 0.5)) {
        y <- 0.9
        mu <- y - r
        num <- if (have_nd) {
          numDeriv::grad(function(mm) scat_logdens(y, mm, sig, nu), mu)
        } else {
          central_diff(function(mm) scat_logdens(y, mm, sig, nu), mu)
        }
        max_err <- max(max_err, abs(num - scat_score(r, nu, sig)))
      }
    }
  }
  expect_lt(max_err, 1e-6)
})


test_that("(ii) large nu recovers the Gaussian score", {
  r <- c(-3, -1, 0.5, 2, 5)
  sig <- 1.4
  gaussian_score <- r / sig^2 # Gaussian location score for known sigma

  expect_equal(scat_score(r, 1e8, sig), gaussian_score, tolerance = 1e-6)

  # convergence is monotone in nu
  e_moderate <- max(abs(scat_score(r, 10, sig) - gaussian_score))
  e_large <- max(abs(scat_score(r, 1e5, sig) - gaussian_score))
  expect_lt(e_large, e_moderate)
})


test_that("(iii) scat CL2 whitened hat trace equals model EDF (within 1%)", {
  skip_on_cran()

  m <- get_scat_test_model()
  b <- m
  class(b) <- setdiff(class(b), "pffr")
  cluster_id <- build_cluster_id(m$pffr)

  # the fit must actually be a scaled-t fit
  expect_match(b$family$family, "^Scaled t")

  # (a) exact-score reconstruction: Xw^T z reproduces the scat score
  X <- model.matrix(b)
  work <- refund:::build_cl2_working_scat(b, cluster_id)
  score_total <- colSums(refund:::compute_scat_scores(b, X))
  reconstructed <- as.vector(crossprod(work$Xw, work$z))
  expect_equal(reconstructed, unname(score_total), tolerance = 1e-8)

  # (b) EDF-consistency: trace of the Fisher-whitened per-cluster hat = sum(edf)
  Vp <- b$Vp
  trace_H <- 0
  for (g in unique(work$cluster_id)) {
    idx <- which(work$cluster_id == g)
    Xwg <- work$Xw[idx, , drop = FALSE]
    trace_H <- trace_H + sum(diag(Xwg %*% Vp %*% t(Xwg)))
  }
  edf_total <- sum(b$edf)
  expect_equal(trace_H, edf_total, tolerance = 0.01 * edf_total)

  # (c) resulting CL2 covariance is finite and symmetric
  V_cl2 <- gam_sandwich_cluster_cl2(b, cluster_id, freq = FALSE)
  expect_equal(V_cl2, t(V_cl2), tolerance = 1e-10)
  expect_true(all(is.finite(V_cl2)))
  expect_true(all(diag(V_cl2) >= 0))

  # (d) coef.pffr() with the scat robust path returns finite SEs, no approx warn
  warn <- captured_warnings(
    co <- coef(m, se = TRUE, sandwich = "cl2", n1 = 20)
  )
  expect_false(any(grepl(approx_msg, warn)))
})


test_that("(iv) approx-score warning fires only for families without exact/two-block scores", {
  # classifier routes each family correctly
  expect_identical(refund:::pffr_score_kind(gaussian()), "exact")
  expect_identical(refund:::pffr_score_kind(poisson()), "exact")
  expect_identical(refund:::pffr_score_kind(Gamma()), "exact")
  expect_identical(refund:::pffr_score_kind(mgcv::gaulss()), "gaulss")
  expect_identical(refund:::pffr_score_kind(mgcv::scat()), "scat")
  expect_identical(refund:::pffr_score_kind(mgcv::nb()), "approx")
  expect_identical(refund:::pffr_score_kind(mgcv::tw()), "approx")

  # the warn-once helper emits the exact review-mandated message for approx
  clear_approx_score_warnings()
  expect_warning(
    refund:::pffr_warn_approx_score(mgcv::nb()),
    approx_msg
  )
  # once per session per family: the second call is silent
  expect_false(
    any(grepl(
      approx_msg,
      captured_warnings(
        refund:::pffr_warn_approx_score(mgcv::nb())
      )
    ))
  )

  # fitted extended families embed estimated parameters in the family string
  # ("Negative Binomial(2.403)"); the warn-once key must strip them so a refit
  # with a different theta does NOT re-warn.
  clear_approx_score_warnings()
  f_th1 <- mgcv::nb()
  f_th1$family <- "Negative Binomial(2.403)"
  f_th2 <- mgcv::nb()
  f_th2$family <- "Negative Binomial(5.1)"
  expect_warning(refund:::pffr_warn_approx_score(f_th1), approx_msg)
  expect_false(
    any(grepl(
      approx_msg,
      captured_warnings(refund:::pffr_warn_approx_score(f_th2))
    ))
  )
  # ... and the unfitted spelling ("negative binomial") shares the same key
  expect_false(
    any(grepl(
      approx_msg,
      captured_warnings(refund:::pffr_warn_approx_score(mgcv::nb()))
    ))
  )
  # gaussian / gaulss / scat never route to the approx warning
  for (fam in list(gaussian(), mgcv::gaulss(), mgcv::scat())) {
    expect_false(refund:::pffr_score_kind(fam) == "approx")
  }

  # real builder path: the warning must actually fire when the sandwich is
  # computed for an approx family, and stay silent for an exact one. Exercise
  # this by relabelling a real gaussian fit's family as an extended family: that
  # flips pffr_score_kind() to "approx" while leaving mu.eta/variance intact, so
  # the generic working-residual path runs and (because gaussian's working
  # residual IS the exact score) the covariance is unchanged.
  skip_on_cran()
  clear_approx_score_warnings()
  m0 <- get_basic_pffr_model()
  b0 <- m0
  class(b0) <- setdiff(class(b0), "pffr")
  cid0 <- build_cluster_id(m0$pffr)

  # exact gaussian path: silent
  warn_exact <- captured_warnings(V_exact <- gam_sandwich_cluster(b0, cid0))
  expect_false(any(grepl(approx_msg, warn_exact)))

  # relabelled as extended.family -> "approx" -> warns; covariance identical
  b_ext <- b0
  class(b_ext$family) <- c("extended.family", class(b_ext$family))
  expect_identical(refund:::pffr_score_kind(b_ext$family), "approx")
  warn_ext <- captured_warnings(V_ext <- gam_sandwich_cluster(b_ext, cid0))
  expect_true(any(grepl(approx_msg, warn_ext)))
  expect_equal(V_ext, V_exact, tolerance = 1e-12)

  # EVERY approximate-score consumer discloses, not just the sandwich
  # builders: the Satterthwaite df context and the X5/X6 shares diagnostic.
  m_ext <- m0
  class(m_ext$family) <- c("extended.family", class(m_ext$family))
  clear_approx_score_warnings()
  warn_ctx <- captured_warnings(
    ctx <- refund:::pffr_df_context(m_ext, "cluster")
  )
  expect_true(any(grepl(approx_msg, warn_ctx)))
  expect_true(ctx$ok)
  clear_approx_score_warnings()
  warn_sh <- captured_warnings(sh <- refund:::pffr_sandwich_shares(m_ext))
  expect_true(any(grepl(approx_msg, warn_sh)))
  expect_true(is.finite(sh$pen_share))
})


test_that("(v) scat with log link: exact score and EDF-consistent whitened hat", {
  skip_on_cran()

  set.seed(50004)
  dat <- pffr_simulate(
    Y ~ xlin,
    n = 40,
    nygrid = 20,
    effects = list(xlin = "dnorm"),
    SNR = 8
  )
  t <- attr(dat, "yindex")
  dat$Y <- exp(0.25 * dat$Y) +
    matrix(rt(length(dat$Y), df = 4) * 0.05, nrow = nrow(dat$Y))
  m <- pffr(
    Y ~ xlin,
    yind = t,
    data = dat,
    family = mgcv::scat(link = "log"),
    bs.yindex = list(bs = "ps", k = 12, m = c(2, 1)),
    sandwich = "none"
  )
  b <- m
  class(b) <- setdiff(class(b), "pffr")
  cluster_id <- build_cluster_id(m$pffr)
  X <- model.matrix(b)

  # the dmu/deta chain-rule factor is genuinely non-constant here
  mu_eta <- as.vector(b$family$mu.eta(b$linear.predictors))
  expect_gt(diff(range(mu_eta)), 0.5)

  # exact score matches mgcv's own Dd()$Dmu on the eta scale
  scores <- refund:::compute_scat_scores(b, X)
  th_raw <- b$family$getTheta(FALSE)
  Dd <- b$family$Dd(
    as.vector(b$y),
    as.vector(b$fitted.values),
    th_raw,
    wt = rep(1, length(b$y))
  )
  scores_mgcv <- (-0.5 * Dd$Dmu) * mu_eta * X
  expect_lt(max(abs(scores - scores_mgcv)), 1e-10)

  # factorization reconstructs the score, and the whitened hat trace = EDF
  work <- refund:::build_cl2_working_scat(b, cluster_id)
  reconstructed <- as.vector(crossprod(work$Xw, work$z))
  expect_equal(reconstructed, unname(colSums(scores)), tolerance = 1e-8)
  Vp <- b$Vp
  trace_H <- 0
  for (g in unique(work$cluster_id)) {
    Xwg <- work$Xw[work$cluster_id == g, , drop = FALSE]
    trace_H <- trace_H + sum(diag(Xwg %*% Vp %*% t(Xwg)))
  }
  edf_total <- sum(b$edf)
  expect_equal(trace_H, edf_total, tolerance = 0.01 * edf_total)
})


test_that("(vi) scat CR1 end-to-end via fit-time sandwich='cluster'", {
  skip_on_cran()

  set.seed(50006)
  dat <- pffr_simulate(
    Y ~ xlin,
    n = 40,
    nygrid = 20,
    effects = list(xlin = "dnorm"),
    SNR = 6
  )
  t <- attr(dat, "yindex")
  dat$Y <- dat$Y +
    matrix(rt(length(dat$Y), df = 3) * 0.2, nrow = nrow(dat$Y))
  m <- pffr(
    Y ~ xlin,
    yind = t,
    data = dat,
    family = mgcv::scat(),
    bs.yindex = list(bs = "ps", k = 12, m = c(2, 1)),
    sandwich = "cluster"
  )

  # fit-time CR1 went through the exact scat branch: robust SEs are usable and
  # no approximation disclosure fired
  warn <- captured_warnings(
    co <- coef(m, se = TRUE, n1 = 20)
  )
  expect_false(any(grepl(approx_msg, warn)))
  ses <- unlist(lapply(co$smterms, function(tm) tm$coef$se))
  expect_true(length(ses) > 0)
  expect_true(all(is.finite(ses)))
  expect_true(all(ses > 0))

  # smoke: the Satterthwaite context and the shares diagnostic run on scat
  ctx <- refund:::pffr_df_context(m, "cl2")
  expect_true(ctx$ok)
  expect_identical(ctx$G, 40L)
  sh <- refund:::pffr_sandwich_shares(m)
  expect_true(is.finite(sh$pen_share))
  expect_gte(sh$pen_share, 0)
  expect_true(is.finite(sh$fro_ratio))
})


test_that("(vii) scat with non-unit prior weights: score exact, hat EDF-consistent", {
  skip_on_cran()

  set.seed(50005)
  n <- 35
  D <- 20
  dat <- pffr_simulate(
    Y ~ xlin,
    n = n,
    nygrid = D,
    effects = list(xlin = "dnorm"),
    SNR = 6
  )
  t <- attr(dat, "yindex")
  dat$Y <- dat$Y +
    matrix(rt(length(dat$Y), df = 3) * 0.2, nrow = nrow(dat$Y))
  W <- matrix(runif(n * D, 0.5, 2), n, D)
  m <- pffr(
    Y ~ xlin,
    yind = t,
    data = dat,
    family = mgcv::scat(),
    weights = W,
    bs.yindex = list(bs = "ps", k = 12, m = c(2, 1)),
    sandwich = "none"
  )
  b <- m
  class(b) <- setdiff(class(b), "pffr")
  cluster_id <- build_cluster_id(m$pffr)
  X <- model.matrix(b)
  pw <- as.vector(b$prior.weights)
  expect_true(any(pw != 1)) # weights actually reached the fit

  # prior-weighted exact score matches mgcv's Dd()$Dmu (which includes wt)
  scores <- refund:::compute_scat_scores(b, X)
  th_raw <- b$family$getTheta(FALSE)
  Dd <- b$family$Dd(
    as.vector(b$y),
    as.vector(b$fitted.values),
    th_raw,
    wt = pw
  )
  mu_eta <- as.vector(b$family$mu.eta(b$linear.predictors))
  scores_mgcv <- (-0.5 * Dd$Dmu) * mu_eta * X
  expect_lt(max(abs(scores - scores_mgcv)), 1e-10)

  # factorization: pw sits in z, NOT in the whitening weight, so the score is
  # still reconstructed exactly ...
  work <- refund:::build_cl2_working_scat(b, cluster_id)
  reconstructed <- as.vector(crossprod(work$Xw, work$z))
  expect_equal(reconstructed, unname(colSums(scores)), tolerance = 1e-8)

  # ... and the whitened hat matches the fit's EDF (mgcv's scat Fisher/EDF
  # weights exclude prior weights; including pw here measured a +21% trace
  # error on exactly this configuration)
  Vp <- b$Vp
  trace_H <- 0
  for (g in unique(work$cluster_id)) {
    Xwg <- work$Xw[work$cluster_id == g, , drop = FALSE]
    trace_H <- trace_H + sum(diag(Xwg %*% Vp %*% t(Xwg)))
  }
  edf_total <- sum(b$edf)
  expect_equal(trace_H, edf_total, tolerance = 0.01 * edf_total)

  # whitening weight equals mgcv's stored Fisher weights (fit$weights)
  expect_equal(
    as.vector(colSums(work$Xw^2)),
    as.vector(colSums((X * sqrt(as.vector(b$weights)))^2)),
    tolerance = 1e-8
  )
})
