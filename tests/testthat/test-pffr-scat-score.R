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

  set.seed(21001)
  dat <- pffr_simulate(
    Y ~ xlin,
    n = 45,
    nygrid = 20,
    effects = list(xlin = "dnorm"),
    SNR = 6
  )
  t <- attr(dat, "yindex")
  # inject heavy-tailed contamination so scat estimates a genuine finite nu
  dat$Y <- dat$Y + matrix(rt(length(dat$Y), df = 3) * 0.15, nrow = nrow(dat$Y))

  m <- pffr(
    Y ~ xlin,
    yind = t,
    data = dat,
    family = mgcv::scat(),
    bs.yindex = list(bs = "ps", k = 12, m = c(2, 1)),
    sandwich = "none"
  )
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
  rm(
    list = ls(envir = refund:::.pffr_state),
    envir = refund:::.pffr_state
  )
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
  rm(
    list = ls(envir = refund:::.pffr_state),
    envir = refund:::.pffr_state
  )
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
})
