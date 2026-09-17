# Fixtures for the exact-CL2 hardening tests (DESIGN amendments A18 and the
# 2026-08-04 study-LB P-LB5 follow-up). The full standalone gates, which also
# write the archived result tables, live in the pffr-ci repo
# (diagnostics/hardening-exactcl2.R); these helpers keep the same coverage
# available to CI without a non-standard top-level file in the package.

make_exactcl2_fixture <- function(family_name, k, design) {
  set.seed(
    match(family_name, c("gaussian", "poisson", "binomial", "Gamma")) *
      100 +
      k +
      match(design, c("unbalanced", "influential"))
  )
  n_curve <- 12L
  n_grid <- 10L
  grid <- seq(0, 1, length.out = n_grid)
  xlin <- stats::rnorm(n_curve)
  eta <- outer(xlin, rep(0.35, n_grid)) +
    matrix(rep(0.2 * sin(2 * pi * grid), each = n_curve), n_curve, n_grid) -
    0.4
  Y <- switch(
    family_name,
    gaussian = eta + matrix(stats::rnorm(n_curve * n_grid, sd = 0.35), n_curve),
    poisson = matrix(
      stats::rpois(n_curve * n_grid, lambda = exp(pmin(eta, 2))),
      n_curve
    ),
    binomial = matrix(
      stats::rbinom(n_curve * n_grid, 1, stats::plogis(eta)),
      n_curve
    ),
    Gamma = matrix(
      stats::rgamma(n_curve * n_grid, shape = 5, scale = exp(pmin(eta, 2)) / 5),
      n_curve
    )
  )
  # Keep the response in its ordinary range, then make one covariate curve
  # extreme. This isolates high leverage from a separable/outcome-saturation
  # failure and exercises the shipped eigenvalue floor.
  if (design == "influential") xlin[1] <- 3000
  # Five deliberately unbalanced independent clusters, with cluster 1 made
  # influential above. Their curve counts 1, 1, 2, 3, 5 yield unequal blocks.
  cluster <- rep(seq_len(5), times = c(1, 1, 2, 3, 5))
  list(
    data = data.frame(Y = I(Y), xlin = xlin),
    yind = grid,
    cluster = cluster,
    family = switch(
      family_name,
      gaussian = stats::gaussian(),
      poisson = stats::poisson(),
      binomial = stats::binomial(),
      Gamma = stats::Gamma(link = "log")
    ),
    k = k
  )
}

run_exactcl2_hardening <- function() {
  cases <- expand.grid(
    family = c("gaussian", "poisson", "binomial", "Gamma"),
    k = c(4L, 7L),
    design = c("unbalanced", "influential"),
    stringsAsFactors = FALSE
  )
  rows <- vector("list", nrow(cases))

  for (i in seq_len(nrow(cases))) {
    case <- cases[i, ]
    fixture <- make_exactcl2_fixture(case$family, case$k, case$design)
    fit <- suppressWarnings(suppressMessages(refund::pffr(
      Y ~ xlin,
      data = fixture$data,
      yind = fixture$yind,
      family = fixture$family,
      bs.yindex = list(bs = "ps", k = fixture$k, m = c(2, 1)),
      sandwich = "none"
    )))

    # `1 - sqrt(tol)` makes the exact floor equal tol: Study EX's
    # cl2_exact_nocap comparator without duplicating the implementation.
    V_exact <- suppressWarnings(refund:::pffr_vcov(
      fit,
      sandwich = "cl2",
      cluster = fixture$cluster,
      cl2_adjustment = "exact"
    ))
    V_nocap <- suppressWarnings(refund:::gam_sandwich_cluster_cl2(
      refund:::pffr_model_based_gam(fit),
      refund:::build_cluster_id(fit$pffr, cluster = fixture$cluster),
      cl2_adjustment = "exact",
      leverage_cap = 1 - sqrt(1e-8)
    ))
    V_shortcut <- suppressWarnings(refund:::pffr_vcov(
      fit,
      sandwich = "cl2",
      cluster = fixture$cluster,
      cl2_adjustment = "shortcut"
    ))

    se_exact <- sqrt(pmax(diag(V_exact), 0))
    se_shortcut <- sqrt(pmax(diag(V_shortcut), 0))
    positive <- is.finite(se_exact) &
      is.finite(se_shortcut) &
      se_exact > 0 &
      se_shortcut > 0
    ratio <- se_exact[positive] / se_shortcut[positive]
    rows[[i]] <- data.frame(
      family = case$family,
      k = case$k,
      design = case$design,
      G = length(unique(fixture$cluster)),
      max_cluster_size = max(table(fixture$cluster)) * length(fixture$yind),
      n_adjusted = attr(V_exact, "n_adjusted"),
      min_block_eig = attr(V_exact, "min_block_eig"),
      max_block_kappa = attr(V_exact, "max_block_kappa"),
      max_abs_exact_nocap = max(abs(V_exact - V_nocap)),
      median_se_ratio_to_shortcut = stats::median(ratio),
      se_finite_positive = all(positive),
      se_ratio_sane = length(ratio) > 0 &&
        min(ratio) > 1e-3 &&
        max(ratio) < 1e3,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

#' Study-LB P-LB5 fixture: a Poisson design whose latent roughness amplitude
#' `amp` pushes individual fitted means far above the nominal marginal mean of
#' 4, the mechanism behind the exploded interval widths in the locked
#' benchmark. `amp = 1` is the benign control; `amp = 10` is degenerate.
make_lb5_fixture <- function(amp, n_grid, k, seed) {
  set.seed(seed)
  n_curve <- 8L
  grid <- seq(0, 1, length.out = n_grid)
  xlin <- stats::rnorm(n_curve)
  nfpc <- 6L
  rough <- matrix(stats::rnorm(n_curve * nfpc), n_curve, nfpc) %*%
    t(sapply(seq_len(nfpc), function(j) sin(j * pi * grid)))
  eta <- outer(xlin, rep(0.35, n_grid)) + amp * rough + log(4)
  Y <- matrix(
    stats::rpois(n_curve * n_grid, lambda = exp(pmin(eta, 600))),
    n_curve
  )
  list(data = data.frame(Y = I(Y), xlin = xlin), yind = grid, k = k)
}

fit_lb5_fixture <- function(fixture) {
  suppressWarnings(suppressMessages(refund::pffr(
    Y ~ xlin,
    data = fixture$data,
    yind = fixture$yind,
    family = stats::poisson(),
    bs.yindex = list(bs = "ps", k = fixture$k, m = c(2, 1)),
    sandwich = "none"
  )))
}
