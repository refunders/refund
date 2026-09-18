#--------------------------------------
# S-F: quasi-likelihood families through the exact GLM score path
#--------------------------------------
#
# quasipoisson()/quasibinomial() carry an estimated dispersion phi-hat
# (`fit$sig2`). The cluster-robust sandwich must use it exactly once, and these
# tests pin the resulting known answer.
#
# Write W = diag(omega_i (dmu/deta)_i^2 / V(mu_i)) (dispersion-free) and
# S = sum_j lambda_j S_j. mgcv's bread is
#     V_p = phi (X'WX + S)^{-1} = (X'WX/phi + S/phi)^{-1},
# i.e. the Fisher information X'WX/phi penalized by S/phi, so
#     V_p(quasi) = phi-hat * V_p(fixed-dispersion fit)     [at the same lambda]
# and likewise for V_e. The per-observation score used by both cluster paths
# (gam_sandwich_cluster() and build_cl2_working_standard()) is
#     s_i = x_i omega_i (y_i - mu_i) (dmu/deta)_i / (phi V(mu_i)),
# so the meat M = sum_g U_g U_g' scales as phi^-2. Hence the sampling core
#     c V_p M V_p  ~  phi^2 * phi^-2 = phi^0
# is DISPERSION-FREE: the quasi fit's CR1/CL2 sampling core is *identical* to
# the fixed-dispersion fit's at the same smoothing parameters. The additive
# Bayesian smoothing-bias allowance B_2 = V_p - V_e is a model-based quantity
# and scales with phi-hat, consistently with the penalty convention S/phi
# above. So the full Bayesian sandwich satisfies
#     V(quasi) - V(fixed) = (phi-hat - 1) * B_2(fixed).
# Per-cluster leverage H_gg = Xw_g V_p Xw_g' is dispersion-free too
# (Xw ~ phi^-1/2, V_p ~ phi), so CL2's leverage adjustment is unchanged.
#
# Both fits below are refitted at the SAME fixed `sp`, so lambda is shared and
# the identities hold to machine precision (a fresh IRLS run in each case; a
# warm-started refit only converges to the fit's own tolerance).

quasi_env <- new.env(parent = emptyenv())

# Overdispersed counts on a functional-covariate design, small enough to fit
# fast and with effective rank(X1) above the ff() identifiability guard.
make_quasi_count_fixture <- function() {
  if (!is.null(quasi_env$count)) {
    return(quasi_env$count)
  }
  set.seed(7)
  G <- 12
  ns <- 10
  ny <- 12
  sgrid <- seq(0, 1, length.out = ns)
  tgrid <- seq(0, 1, length.out = ny)
  X1 <- matrix(rnorm(G * ns), nrow = G, ncol = ns)
  beta <- outer(sgrid, tgrid, function(s, t) 0.6 * cos(pi * s) * (1 + t))
  mu <- exp(1.2 + (X1 %*% beta) / ns)
  # negative-binomial draws => genuine overdispersion => phi-hat clearly > 1
  Y <- matrix(
    rnbinom(length(mu), mu = as.vector(mu), size = as.vector(mu) / 1.5),
    nrow = G,
    ncol = ny
  )
  quasi_env$count <- list(
    data = list(Y = Y, X1 = X1),
    sgrid = sgrid,
    tgrid = tgrid,
    G = G
  )
  quasi_env$count
}

make_quasi_binary_fixture <- function() {
  if (!is.null(quasi_env$binary)) {
    return(quasi_env$binary)
  }
  set.seed(11)
  G <- 12
  ns <- 10
  ny <- 12
  sgrid <- seq(0, 1, length.out = ns)
  tgrid <- seq(0, 1, length.out = ny)
  X1 <- matrix(rnorm(G * ns), nrow = G, ncol = ns)
  beta <- outer(sgrid, tgrid, function(s, t) 1.5 * cos(pi * s) * (1 + t))
  p <- plogis(0.2 + (X1 %*% beta) / ns)
  Y <- matrix(rbinom(length(p), 1, as.vector(p)), nrow = G, ncol = ny)
  quasi_env$binary <- list(
    data = list(Y = Y, X1 = X1),
    sgrid = sgrid,
    tgrid = tgrid,
    G = G
  )
  quasi_env$binary
}

fit_quasi_fixture <- function(fx, family, sp = NULL) {
  args <- list(
    formula = Y ~
      ff(
        X1,
        xind = fx$sgrid,
        splinepars = list(
          bs = "ps",
          k = c(5, 5),
          m = list(c(2, 1), c(2, 1))
        )
      ),
    data = fx$data,
    yind = fx$tgrid,
    family = family,
    sandwich = "none",
    bs.yindex = list(bs = "ps", k = 5, m = c(2, 1))
  )
  if (!is.null(sp)) {
    args$sp <- sp
  }
  suppressWarnings(do.call(pffr, args))
}

# One fixed-dispersion fit selects lambda; BOTH arms are then refitted at that
# sp, so the two fits differ only in the dispersion treatment.
get_quasi_pair <- function(which = c("count", "binary")) {
  which <- match.arg(which)
  if (!is.null(quasi_env[[paste0("pair_", which)]])) {
    return(quasi_env[[paste0("pair_", which)]])
  }
  if (which == "count") {
    fx <- make_quasi_count_fixture()
    fixed_family <- poisson()
    quasi_family <- quasipoisson()
  } else {
    fx <- make_quasi_binary_fixture()
    fixed_family <- binomial()
    quasi_family <- quasibinomial()
  }
  sp <- fit_quasi_fixture(fx, fixed_family)$sp
  pair <- list(
    fx = fx,
    fixed = fit_quasi_fixture(fx, fixed_family, sp = sp),
    quasi = fit_quasi_fixture(fx, quasi_family, sp = sp)
  )
  quasi_env[[paste0("pair_", which)]] <- pair
  pair
}

rel_diff <- function(a, b) max(abs(a - b)) / max(abs(b))

test_that("quasi families are classified as an exact GLM score path", {
  expect_identical(refund:::pffr_score_kind(quasipoisson()), "exact")
  expect_identical(refund:::pffr_score_kind(quasibinomial()), "exact")
  expect_true(refund:::family_has_exact_score(quasipoisson()))
  expect_true(refund:::family_has_exact_score(quasibinomial()))
})

test_that("sandwich = 'auto' routes quasi families through the exact path", {
  # (c): quasi families are promoted to CL2 exactly like their fixed-dispersion
  # counterparts, and are NOT lumped in with the "approx" extended families.
  expect_identical(
    refund:::pffr_sandwich_auto_policy(12, 12, quasipoisson()),
    refund:::pffr_sandwich_auto_policy(12, 12, poisson())
  )
  expect_identical(
    refund:::pffr_sandwich_auto_policy(12, 12, quasipoisson()),
    "cl2"
  )
  expect_identical(
    refund:::pffr_sandwich_auto_policy(12, 12, quasibinomial()),
    "cl2"
  )
  # contrast: an extended family with only the working-residual APPROXIMATION
  expect_identical(
    refund:::pffr_sandwich_auto_policy(12, 12, mgcv::nb()),
    "cluster"
  )
})

test_that("the quasi score carries phi-hat exactly once", {
  pair <- get_quasi_pair("count")
  phi <- pair$quasi$sig2
  expect_gt(phi, 1.2)
  expect_equal(pair$fixed$sig2, 1)
  # same lambda, same IRLS solution
  expect_equal(
    unname(pair$quasi$coefficients),
    unname(pair$fixed$coefficients),
    tolerance = 1e-10
  )

  strip <- function(b) {
    class(b) <- setdiff(class(b), "pffr")
    b
  }
  cid <- refund:::build_cluster_id(pair$fixed$pffr)
  # the exact score path is taken: no working-residual disclosure fires
  expect_no_warning(
    refund:::gam_sandwich_cluster(strip(pair$quasi), cid)
  )
  wq <- refund:::build_cl2_working_standard(strip(pair$quasi), cid)
  wf <- refund:::build_cl2_working_standard(strip(pair$fixed), cid)
  # per-cluster score sums U_g = sum_i s_i: phi enters the denominator once
  Uq <- rowsum(wq$Xw * wq$z, cid)
  Uf <- rowsum(wf$Xw * wf$z, cid)
  expect_lt(rel_diff(Uq, Uf / phi), 1e-8)
  # ... and the bread carries it once, in the numerator
  expect_lt(rel_diff(pair$quasi$Vp, phi * pair$fixed$Vp), 1e-8)
  expect_lt(rel_diff(pair$quasi$Ve, phi * pair$fixed$Ve), 1e-8)
})

test_that("quasipoisson CR1 sampling core equals the Poisson CR1 core", {
  # (a): B M B with B ~ phi and M ~ phi^-2 => the core is phi^0, i.e. the two
  # covariances are IDENTICAL, not rescaled.
  pair <- get_quasi_pair("count")
  phi <- pair$quasi$sig2
  Vq <- refund:::pffr_vcov(pair$quasi, sandwich = "cluster", freq = TRUE)
  Vf <- refund:::pffr_vcov(pair$fixed, sandwich = "cluster", freq = TRUE)
  expect_lt(rel_diff(Vq, Vf), 1e-8)

  # the full Bayesian sandwich differs only through the phi-scaled B2 term
  Bq <- refund:::pffr_vcov(pair$quasi, sandwich = "cluster")
  Bf <- refund:::pffr_vcov(pair$fixed, sandwich = "cluster")
  B2_fixed <- pair$fixed$Vp - pair$fixed$Ve
  expect_lt(rel_diff(Bq - Bf, (phi - 1) * B2_fixed), 1e-8)
  expect_lt(rel_diff(Bq, Vf + phi * B2_fixed), 1e-8)
})

test_that("quasibinomial CR1 sampling core equals the binomial CR1 core", {
  # (b)
  pair <- get_quasi_pair("binary")
  phi <- pair$quasi$sig2
  expect_false(isTRUE(all.equal(phi, 1, tolerance = 1e-3)))
  Vq <- refund:::pffr_vcov(pair$quasi, sandwich = "cluster", freq = TRUE)
  Vf <- refund:::pffr_vcov(pair$fixed, sandwich = "cluster", freq = TRUE)
  expect_lt(rel_diff(Vq, Vf), 1e-8)

  Bq <- refund:::pffr_vcov(pair$quasi, sandwich = "cluster")
  Bf <- refund:::pffr_vcov(pair$fixed, sandwich = "cluster")
  B2_fixed <- pair$fixed$Vp - pair$fixed$Ve
  expect_lt(rel_diff(Bq - Bf, (phi - 1) * B2_fixed), 1e-8)
})

test_that("CL2 works for quasi fits and its leverage is dispersion-free", {
  # (d)
  pair <- get_quasi_pair("count")
  cl_q <- refund:::pffr_vcov(pair$quasi, sandwich = "cl2", freq = TRUE)
  cl_f <- refund:::pffr_vcov(pair$fixed, sandwich = "cl2", freq = TRUE)
  cr_q <- refund:::pffr_vcov(pair$quasi, sandwich = "cluster", freq = TRUE)
  cr_f <- refund:::pffr_vcov(pair$fixed, sandwich = "cluster", freq = TRUE)
  # leverage H_gg = Xw_g Vp Xw_g' is phi-free, so the CL2 inflation matches
  expect_lt(rel_diff(cl_q, cl_f), 1e-8)
  expect_lt(
    max(abs(sqrt(diag(cl_q) / diag(cr_q)) - sqrt(diag(cl_f) / diag(cr_f)))),
    1e-8
  )
  expect_equal(
    attr(cl_q, "max_obs_leverage"),
    attr(cl_f, "max_obs_leverage"),
    tolerance = 1e-8
  )

  # coef(fit, sandwich = "cl2") returns finite SEs on a quasipoisson fit
  cf <- suppressWarnings(coef(pair$quasi, sandwich = "cl2"))
  ses <- unlist(lapply(cf$smterms, function(s) s$coef$se))
  expect_gt(length(ses), 0)
  expect_true(all(is.finite(ses)))
  expect_true(all(ses > 0))
})
