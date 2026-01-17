# confint-benchmark.R
#
# Concept + scaffold: benchmark pointwise CI coverage/width from `coef.pffr()`
# under different residual covariance structures and error distributions, for:
# - pffr (default covariance)
# - pffr_gls (oracle hatSigma)
# - pffr with AR(1) residuals (oracle rho; bam only)
# - pffr with sandwich correction (post-hoc `coef(..., sandwich = TRUE)`)
#
# Two likelihoods for fitting:
# - gaussian()
# - mgcv::gaulss() (location-scale)
#
# Data generation uses `pffr_simulate()` to generate covariates and baseline signal,
# then overwrites the response with structured errors and adds a concurrent
# functional effect. This keeps the workflow close to pffr/pffr_simulate while giving
# explicit control over (hetero-)skedasticity, autocorrelation, SNR, and tails.
#
# Design sketch (factorial; start small, then expand):
# - term types: ff (beta(s,t)), concurrent functional (beta_c(t)),
#   linear varying coefficient of scalar (beta_z(t)), smooth varying effect
#   of scalar (f(z,t))
# - true error structure: iid, heteroskedastic only, autocorrelated only,
#   and combined heteroskedastic + autocorrelation; multiple shapes/strengths
# - SNR: e.g. 5, 10, 25 (global var(signal)/var(error))
# - error distribution: gaussian vs t(df=6)
# - fit likelihood: gaussian vs gaulss
# - replicates: e.g. 200+ (more for stable tail scenarios)
#
# Evaluation sketch:
# - use `coef(fit, sandwich = ...)` to extract pointwise estimates + SEs
# - build 1 - alpha CIs: estimate +/- qnorm(1 - alpha/2) * SE
# - coverage: mean(truth within CI) pointwise, aggregated over grid
# - width: mean(2 * qnorm(...) * SE) pointwise, aggregated over grid
# - report both integrated metrics (averaged over grid) and domain-resolved
#   patterns (coverage/width vs t, and heatmaps over (s,t) for ff and (z,t) for
#   s(z)-type terms)
#
# VERIFIED FINDINGS on undercoverage (2024-01 investigation):
#
# ROOT CAUSE: Smoothing bias, NOT concurvity.
# - At extreme SNR (1e6), ff term z-score SD ~3-4, causing ~47% coverage
# - This occurs even in single-term (ff-only) models with simple truth
# - The penalized spline shrinks estimates toward smoother functions, creating
#   systematic bias that the Bayesian posterior SEs don't account for
#
# SOLUTION: Use "spline-friendly" truth functions + moderate SNR
# - Truth functions should lie in the low-rank spline space (low-frequency
#   sinusoids, polynomials) to avoid smoothing bias
# - SNR=25 with simple truth achieves near-nominal coverage across all terms:
#   * ff: 97% (Z-SD=0.82)
#   * linear: 97% (Z-SD=0.96)
#   * concurrent: 93-96%
#   * smooth: 95-97%
# - SNR=100 also works well with simple truth
# - SNR>=1000 causes undercoverage even with simple truth
#
# IMPLICATION: coef.pffr() CIs are valid for functions that splines can
# represent without aggressive penalization. For "wiggly" truth, coverage
# degrades due to smoothing bias - this is expected behavior, not a bug.
#
# NONLINEARITY TESTS (2024-01):
# Coverage holds for moderate nonlinearities:
#   - linear: 95.6% (Z-SD=1.02)
#   - quadratic (monotonic): 94.2% (Z-SD=1.04)
#   - quadratic (U-shape): 92.8% (Z-SD=1.09)
#   - cos(pi*t) (half period): 93.9% (Z-SD=1.09)
# But fails for high-frequency/full-period functions:
#   - sin(2*pi*t): 12% coverage - identifiability issues at boundaries
#
# T-DISTRIBUTED ERRORS:
# The scat() family (scaled t) is supported for fitting models when errors
# are t-distributed. Coverage is similar to Gaussian family in tested cases.
#
# Visualization sketch (ggplot2 recommended):
# - coverage vs method: faceted by term type x (error dist / covariance / SNR)
# - width vs method: same faceting
# - 1D term functions (beta(t)): lines over t with ribbons for widths
# - 2D surfaces: heatmaps of coverage (or coverage - nominal) over (s,t) and
#   (z,t), faceted by method

# Setup ----

load_refund_dev <- function(pkg_dir = ".") {
  if (
    requireNamespace("devtools", quietly = TRUE) &&
      file.exists(file.path(pkg_dir, "DESCRIPTION"))
  ) {
    devtools::load_all(pkg_dir, quiet = TRUE)
  } else {
    suppressPackageStartupMessages(library(refund))
  }
  invisible(TRUE)
}

ensure_pffr_simulate_bench <- function(pkg_dir = ".") {
  if (
    exists("pffr_simulate_bench", envir = .GlobalEnv) &&
      is.function(get("pffr_simulate_bench", envir = .GlobalEnv)) &&
      "effects" %in%
        names(formals(get("pffr_simulate_bench", envir = .GlobalEnv)))
  ) {
    return(invisible(TRUE))
  }

  util_path <- file.path(pkg_dir, "R", "pffr-utilities.R")
  if (file.exists(util_path)) {
    # Prefer the repo-local simulation helper (avoids installed-version drift
    # and avoids calling deprecated wrappers like pffrSim()).
    util_env <- new.env(parent = .GlobalEnv)
    sys.source(util_path, envir = util_env)

    sim_fun <- NULL
    if ("pffr_simulate" %in% names(util_env)) {
      sim_fun <- util_env$pffr_simulate
    } else if ("pffrSim" %in% names(util_env)) {
      sim_fun <- util_env$pffrSim
    }

    if (is.null(sim_fun) || !"effects" %in% names(formals(sim_fun))) {
      stop(
        "Local `R/pffr-utilities.R` did not provide formula-interface simulation."
      )
    }

    assign("pffr_simulate_bench", sim_fun, envir = .GlobalEnv)
    return(invisible(TRUE))
  }

  # Fallback to search path.
  if (exists("pffr_simulate") && "effects" %in% names(formals(pffr_simulate))) {
    assign("pffr_simulate_bench", pffr_simulate, envir = .GlobalEnv)
    return(invisible(TRUE))
  }
  if (exists("pffrSim") && "effects" %in% names(formals(pffrSim))) {
    assign("pffr_simulate_bench", pffrSim, envir = .GlobalEnv)
    return(invisible(TRUE))
  }
  stop("Could not find a formula-interface simulation function to use.")
}

require_pkgs <- function() {
  load_refund_dev(".")
  stopifnot(
    requireNamespace("mgcv", quietly = TRUE),
    requireNamespace("mvtnorm", quietly = TRUE)
  )

  ensure_pffr_simulate_bench(".")
  invisible(TRUE)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

# Truth specification ----
#
# Per-replicate truth functions are drawn from low-rank B-spline bases.
# The wiggliness parameter scales coefficient variance so truth varies
# across replicates and across low/high wiggliness conditions.

make_bs_spec <- function(x, df) {
  basis <- splines::bs(x, df = df, intercept = TRUE)
  list(
    knots = attr(basis, "knots"),
    boundary_knots = attr(basis, "Boundary.knots"),
    degree = attr(basis, "degree"),
    intercept = attr(basis, "intercept"),
    df = df
  )
}

eval_bs <- function(x, spec) {
  splines::bs(
    x,
    df = spec$df,
    knots = spec$knots,
    Boundary.knots = spec$boundary_knots,
    degree = spec$degree,
    intercept = spec$intercept
  )
}

make_truth_functions <- function(
  s_grid,
  t_grid,
  z_ref,
  wiggliness,
  seed = NULL,
  k_te = 5L,
  k_uni = 7L
) {
  if (!is.null(seed)) set.seed(seed)

  rw_coef <- function(k, sd) {
    coef <- cumsum(rnorm(k, sd = sd))
    coef - mean(coef)
  }

  rw_coef_matrix <- function(k, sd) {
    base <- matrix(rnorm(k * k, sd = sd), nrow = k)
    mat <- apply(base, 2, cumsum)
    mat <- t(apply(mat, 1, cumsum))
    mat - mean(mat)
  }

  s_spec <- make_bs_spec(s_grid, df = k_te)
  t_te_spec <- make_bs_spec(t_grid, df = k_te)
  t_uni_spec <- make_bs_spec(t_grid, df = k_uni)
  z_ref <- as.numeric(z_ref %||% 0)
  if (length(z_ref) < 2) z_ref <- c(-0.5, 0.5)
  z_te_spec <- make_bs_spec(z_ref, df = k_te)

  coef_ff <- rw_coef_matrix(k_te, wiggliness)
  coef_con <- rw_coef(k_uni, wiggliness)
  coef_lin <- rw_coef(k_uni, wiggliness)
  coef_mu <- rw_coef(k_uni, wiggliness)
  coef_smooth <- rw_coef_matrix(k_te, wiggliness)

  beta_ff <- function(s, t) {
    Bs <- eval_bs(s, s_spec)
    Bt <- eval_bs(t, t_te_spec)
    drop(Bs %*% coef_ff %*% t(Bt))
  }

  beta_concurrent <- function(t) {
    Bt <- eval_bs(t, t_uni_spec)
    drop(Bt %*% coef_con)
  }

  beta_linear <- function(t) {
    Bt <- eval_bs(t, t_uni_spec)
    drop(Bt %*% coef_lin)
  }

  mu_t <- function(t) {
    Bt <- eval_bs(t, t_uni_spec)
    drop(Bt %*% coef_mu)
  }

  f_smooth_scalar <- function(z, t) {
    Bz <- eval_bs(z, z_te_spec)
    Bt <- eval_bs(t, t_te_spec)
    drop(Bz %*% coef_smooth %*% t(Bt))
  }

  f_smooth_matrix <- function(z_vec, t_vec) {
    f_smooth_scalar(z_vec, t_vec)
  }

  list(
    beta_ff = beta_ff,
    beta_concurrent = beta_concurrent,
    beta_linear = beta_linear,
    mu_t = mu_t,
    f_smooth_scalar = f_smooth_scalar,
    f_smooth_matrix = f_smooth_matrix,
    coef = list(
      ff = coef_ff,
      concurrent = coef_con,
      linear = coef_lin,
      mu_t = coef_mu,
      smooth = coef_smooth
    ),
    spec = list(
      s = s_spec,
      t_te = t_te_spec,
      t_uni = t_uni_spec,
      z_te = z_te_spec
    )
  )
}

# Effect constructors (DGP side) ----

ff_weights <- function(xind, method = c("simpson", "trapezoidal", "riemann")) {
  method <- match.arg(method)
  nx <- length(xind)
  if (nx < 2) return(1)

  if (method == "simpson") {
    w <- ((xind[nx] - xind[1]) / (nx - 1)) /
      3 *
      c(1, rep(c(4, 2), length = nx - 2), 1)
    return(w)
  }

  if (method == "trapezoidal") {
    diffs <- diff(xind)
    w <- c(
      diffs[1] / 2,
      (diffs[-1] + diffs[-length(diffs)]) / 2,
      diffs[length(diffs)] / 2
    )
    return(w)
  }

  diffs <- diff(xind)
  w <- c(mean(diffs), diffs)
  w
}

center_beta_ff <- function(beta_st, s_grid, weights) {
  w <- weights %||% ff_weights(s_grid, method = "simpson")
  w_sum <- sum(w)
  col_means <- colSums(beta_st * w) / w_sum
  sweep(beta_st, 2, col_means, "-")
}

compute_ff_effect <- function(X, beta_st, xind) {
  # Integration weights (Simpson's rule, consistent with ff() defaults)
  w <- ff_weights(xind, method = "simpson")
  L <- matrix(rep(w, times = nrow(X)), nrow = nrow(X), byrow = TRUE)
  (L * X) %*% beta_st
}

compute_linear_effect <- function(x, beta_t) {
  outer(x, beta_t)
}

# Covariate + response generation ----

random_function_matrix <- function(n, grid, bs_dim = 25L) {
  # Generate n random smooth-ish curves observed on `grid` using spline basis.
  # Returns an n x length(grid) matrix.
  stopifnot(is.numeric(n), length(n) == 1, n > 0)
  stopifnot(is.numeric(grid), length(grid) > 1)
  stopifnot(is.numeric(bs_dim), length(bs_dim) == 1, bs_dim >= 5)

  X <- splines::bs(grid, df = bs_dim, intercept = TRUE)

  coef_mat <- matrix(rnorm(n * bs_dim), nrow = bs_dim, ncol = n)
  t(X %*% coef_mat)
}

make_error_structure <- function(
  t_grid,
  corr_type = c("iid", "ar1", "gauss", "fourier"),
  corr_param = NULL,
  hetero_type = c("none", "linear", "u", "bump"),
  hetero_param = NULL
) {
  corr_type <- match.arg(corr_type)
  hetero_type <- match.arg(hetero_type)
  ny <- length(t_grid)
  dist_mat <- abs(outer(t_grid, t_grid, "-"))

  make_pd <- function(mat) {
    diag_mean <- mean(diag(mat))
    if (!is.finite(diag_mean) || diag_mean <= 0) diag_mean <- 1

    jitter <- 1e-10 * diag_mean
    for (k in seq_len(8)) {
      ok <- tryCatch(
        {
          chol(mat + diag(jitter, nrow(mat)))
          TRUE
        },
        error = function(e) FALSE
      )
      if (ok) return(mat + diag(jitter, nrow(mat)))
      jitter <- jitter * 10
    }
    mat + diag(jitter, nrow(mat))
  }

  corr_mat <- switch(
    corr_type,
    iid = diag(1, ny),
    ar1 = {
      rho <- corr_param %||% 0.4
      rho^abs(outer(seq_len(ny), seq_len(ny), "-"))
    },
    gauss = {
      phi <- corr_param %||% 0.15
      exp(-(dist_mat^2) / (2 * phi^2))
    },
    fourier = {
      # Non-banded, non-monotone correlation constructed via latent Fourier
      # factors; PD by construction.
      period1 <- corr_param %||% 0.33
      period2 <- 0.17
      B <- cbind(
        sin(2 * pi * t_grid / period1),
        cos(2 * pi * t_grid / period1),
        sin(2 * pi * t_grid / period2),
        cos(2 * pi * t_grid / period2)
      )
      B <- scale(B, center = TRUE, scale = FALSE)
      base_cov <- B %*% t(B)
      base_cov <- base_cov + diag(0.25, ny)
      stats::cov2cor(base_cov)
    }
  )
  corr_mat <- make_pd(corr_mat)

  sigma_t <- switch(
    hetero_type,
    none = rep(1, ny),
    linear = {
      amp <- hetero_param %||% 0.6
      1 + amp * (t_grid - mean(t_grid))
    },
    u = {
      amp <- hetero_param %||% 0.9
      shape <- (t_grid - 0.5)^2
      shape <- shape / max(shape)
      1 + amp * shape
    },
    bump = {
      amp <- hetero_param %||% 1.0
      center <- 0.7
      width <- 0.10
      1 + amp * exp(-0.5 * ((t_grid - center) / width)^2)
    }
  )
  sigma_t <- pmax(sigma_t, 1e-6)

  cov_mat <- diag(sigma_t, ny) %*% corr_mat %*% diag(sigma_t, ny)
  cov_mat <- make_pd(cov_mat)

  list(
    t_grid = t_grid,
    corr_type = corr_type,
    corr_param = corr_param,
    hetero_type = hetero_type,
    hetero_param = hetero_param,
    sigma_t = sigma_t,
    corr_mat = corr_mat,
    cov_mat = cov_mat,
    rho = if (corr_type == "ar1") (corr_param %||% 0.4) else NA_real_
  )
}

sample_errors <- function(n, cov_mat, error_dist = c("gaussian", "t6")) {
  error_dist <- match.arg(error_dist)
  ny <- nrow(cov_mat)
  stopifnot(ncol(cov_mat) == ny)

  if (error_dist == "gaussian") {
    return(mvtnorm::rmvnorm(n = n, mean = rep(0, ny), sigma = cov_mat))
  }

  df <- 6
  # rmvt() parameterizes via scale matrix; Cov = sigma * df/(df-2)
  sigma_scale <- cov_mat * (df - 2) / df
  mvtnorm::rmvt(n = n, sigma = sigma_scale, df = df, delta = rep(0, ny))
}

compute_snr_scale_factor <- function(eps, eta, snr_target) {
  stopifnot(is.matrix(eps), is.matrix(eta))
  stopifnot(all(dim(eps) == dim(eta)))
  stopifnot(is.numeric(snr_target), length(snr_target) == 1, snr_target > 0)

  var_eta <- stats::var(as.vector(eta))
  var_eps <- stats::var(as.vector(eps))
  if (!is.finite(var_eta) || var_eta <= 0) var_eta <- 1
  if (!is.finite(var_eps) || var_eps <= 0) var_eps <- 1

  sqrt(var_eta / (snr_target * var_eps))
}

rescale_errors_to_snr <- function(eps, eta, snr_target) {
  scale_fac <- compute_snr_scale_factor(eps, eta, snr_target)
  eps * scale_fac
}

simulate_dataset <- function(
  n,
  nxgrid,
  nygrid,
  snr,
  wiggliness,
  error_dist = c("gaussian", "t6"),
  corr_type = c("iid", "ar1", "gauss", "fourier"),
  corr_param = NULL,
  hetero_type = c("none", "linear", "u", "bump"),
  hetero_param = NULL,
  seed = NULL,
  terms = c("ff", "concurrent", "linear", "smooth")
) {
  covar_seed <- if (is.null(seed)) NULL else as.integer(seed)
  truth_seed <- if (is.null(seed)) NULL else as.integer(seed) + 100000L
  error_seed <- if (is.null(seed)) NULL else as.integer(seed) + 200000L
  error_dist <- match.arg(error_dist)
  corr_type <- match.arg(corr_type)
  hetero_type <- match.arg(hetero_type)
  terms <- unique(terms)

  base_terms <- c(
    if ("ff" %in% terms) "ff(X1, xind = s)" else NULL,
    if ("linear" %in% terms) "zlin" else NULL,
    if ("smooth" %in% terms) "s(zsmoo)" else NULL
  )
  base_formula <- if (length(base_terms)) {
    as.formula(paste("Y ~", paste(base_terms, collapse = " + ")))
  } else {
    Y ~ 1
  }

  if (!is.null(covar_seed)) set.seed(covar_seed)
  base <- pffr_simulate_bench(
    base_formula,
    n = n,
    nxgrid = nxgrid,
    nygrid = nygrid,
    SNR = snr, # used only for base (ignored once we overwrite Y)
    family = gaussian(),
    effects = list(
      X1 = 0,
      zlin = 0,
      zsmoo = 0
    ),
    intercept = "beta"
  )
  s_grid <- attr(base, "xindex")
  t_grid <- attr(base, "yindex")
  base_df <- base

  truth <- make_truth_functions(
    s_grid = s_grid,
    t_grid = t_grid,
    z_ref = base_df$zsmoo,
    wiggliness = wiggliness,
    seed = truth_seed
  )

  # Add concurrent functional covariate Xc(t) and its signal contribution.
  Xc <- if ("concurrent" %in% terms) {
    random_function_matrix(n = n, grid = t_grid, bs_dim = 9L)
  } else {
    NULL
  }

  # Enforce pffr-friendly centering so the functional intercept equals mu(t),
  # and term contributions have mean 0 (across i) for each t.
  X1 <- if ("ff" %in% terms) {
    X1 <- as.matrix(base_df$X1)
    X1 <- sweep(X1, 2, colMeans(X1), "-")
    base_df$X1 <- I(X1)
    X1
  } else {
    NULL
  }

  zlin <- if ("linear" %in% terms) {
    zlin <- drop(base_df$zlin)
    zlin <- zlin - mean(zlin)
    base_df$zlin <- zlin
    zlin
  } else {
    NULL
  }

  zsmoo <- if ("smooth" %in% terms) {
    zsmoo <- drop(base_df$zsmoo)
    zsmoo <- zsmoo - mean(zsmoo)
    base_df$zsmoo <- zsmoo
    zsmoo
  } else {
    NULL
  }

  if (!is.null(Xc)) {
    Xc <- sweep(Xc, 2, colMeans(Xc), "-")
  }

  # Recompute the true eta from centered covariates to match pffr constraints.
  mu_t <- truth$mu_t(t_grid)
  eta_int <- matrix(rep(mu_t, each = n), nrow = n, ncol = length(t_grid))

  beta_st <- if ("ff" %in% terms) truth$beta_ff(s_grid, t_grid) else NULL
  w_s <- if ("ff" %in% terms) ff_weights(s_grid, method = "simpson") else NULL
  beta_st <- if ("ff" %in% terms) center_beta_ff(beta_st, s_grid, w_s) else NULL
  eta_ff <- if ("ff" %in% terms)
    compute_ff_effect(
      X = X1,
      beta_st = beta_st,
      xind = s_grid
    ) else NULL

  beta_lin <- if ("linear" %in% terms) truth$beta_linear(t_grid) else NULL
  eta_lin <- if ("linear" %in% terms) compute_linear_effect(zlin, beta_lin) else
    NULL

  beta_con <- if ("concurrent" %in% terms) truth$beta_concurrent(t_grid) else
    NULL
  eta_con <- if ("concurrent" %in% terms) Xc * rep(beta_con, each = n) else NULL

  eta_smooth <- NULL
  if ("smooth" %in% terms) {
    eta_smooth_raw <- truth$f_smooth_matrix(zsmoo, t_grid)
    # Center smooth contribution per t (pffr uses sum_i f(z_i, t)=0 for all t).
    eta_smooth <- sweep(eta_smooth_raw, 2, colMeans(eta_smooth_raw), "-")
  }

  eta_total <- eta_int
  if (!is.null(eta_ff)) eta_total <- eta_total + eta_ff
  if (!is.null(eta_con)) eta_total <- eta_total + eta_con
  if (!is.null(eta_lin)) eta_total <- eta_total + eta_lin
  if (!is.null(eta_smooth)) eta_total <- eta_total + eta_smooth

  # Structured errors (common cov across i), rescaled to target SNR.
  err <- make_error_structure(
    t_grid = t_grid,
    corr_type = corr_type,
    corr_param = corr_param,
    hetero_type = hetero_type,
    hetero_param = hetero_param
  )
  if (!is.null(error_seed)) set.seed(error_seed)
  eps <- sample_errors(n = n, cov_mat = err$cov_mat, error_dist = error_dist)
  snr_scale <- compute_snr_scale_factor(
    eps = eps,
    eta = eta_total,
    snr_target = snr
  )
  eps <- eps * snr_scale
  err$cov_mat <- err$cov_mat * snr_scale^2

  Y <- eta_total + eps

  dat <- base_df
  dat$Y <- I(Y)
  if (!is.null(Xc)) dat$Xc <- I(Xc)
  attr(dat, "truth") <- list(
    eta = eta_total,
    etaTerms = list(
      intercept = eta_int,
      X1 = eta_ff,
      Xc = eta_con,
      zlin = eta_lin,
      zsmoo = eta_smooth
    ),
    beta = list(
      intercept = mu_t,
      X1 = beta_st,
      zlin = beta_lin,
      zsmoo = truth$f_smooth_scalar,
      Xc = beta_con
    ),
    epsilon = eps,
    truth_funs = list(
      beta_ff = truth$beta_ff,
      beta_concurrent = truth$beta_concurrent,
      beta_linear = truth$beta_linear,
      f_smooth_scalar = truth$f_smooth_scalar,
      mu_t = truth$mu_t
    ),
    wiggliness = wiggliness,
    center_ref = list(
      zsmoo = zsmoo,
      ff_weights = w_s
    ),
    term_set = terms,
    grids = list(s = s_grid, t = t_grid),
    dgp = list(
      snr = snr,
      wiggliness = wiggliness,
      error_dist = error_dist,
      corr_type = corr_type,
      corr_param = corr_param,
      hetero_type = hetero_type,
      hetero_param = hetero_param,
      cov_mat = err$cov_mat
    )
  )

  list(data = dat, s_grid = s_grid, t_grid = t_grid, err = err)
}

# Fitting methods ----

# Default basis dimensions: k=8 per margin for tensor products (2D), k=12 for univariate (1D)
# These larger values prevent smoothing bias for moderately complex truth functions
default_fit_args <- function() {
  list(
    # Tensor product bases for response direction (yindex) and intercept
    bs.yindex = list(bs = "ps", k = 12, m = c(2, 1)),
    bs.int = list(bs = "ps", k = 12, m = c(2, 1)),
    # For ff() term: k per margin in the (s,t) tensor product
    # Passed via formula: ff(X1, xind=s_grid, basistype="s", k=8)
    # For s(zsmoo): univariate smooth in z direction, tensor with t
    # These are set in the formula construction
    k_ff = 8, # k per margin for ff() tensor product
    k_smooth = 12 # k for univariate smooths like s(zsmoo)
  )
}

fit_family_object <- function(fit_family = c("gaussian", "gaulss", "scat")) {
  fit_family <- match.arg(fit_family)
  if (fit_family == "gaussian") return(stats::gaussian())
  if (fit_family == "gaulss") return(mgcv::gaulss())
  if (fit_family == "scat") return(mgcv::scat()) # scaled t-distribution
  stop("Unknown fit_family: ", fit_family)
}

fitted_mean_matrix <- function(fit) {
  fm <- tryCatch(fitted(fit, which = "mean"), error = function(e) fitted(fit))
  if (is.list(fm) && !is.null(fm$mean)) return(fm$mean)
  fm
}

estimate_hatSigma_from_fit <- function(fit, dat) {
  y_mat <- as.matrix(dat$Y)
  fit_mat <- as.matrix(fitted_mean_matrix(fit))
  stopifnot(all(dim(y_mat) == dim(fit_mat)))

  resid_mat <- y_mat - fit_mat
  hatSigma <- stats::cov(resid_mat)
  hatSigma <- 0.5 * (hatSigma + t(hatSigma))

  diag_mean <- mean(diag(hatSigma))
  if (!is.finite(diag_mean) || diag_mean <= 0) diag_mean <- 1

  ridge <- 1e-8 * diag_mean
  hatSigma <- hatSigma + diag(ridge, nrow(hatSigma))
  hatSigma
}

estimate_rho_from_fit <- function(fit, dat, clamp = 0.99) {
  y_mat <- as.matrix(dat$Y)
  fit_mat <- as.matrix(fitted_mean_matrix(fit))
  stopifnot(all(dim(y_mat) == dim(fit_mat)))

  resid_mat <- y_mat - fit_mat
  ny <- ncol(resid_mat)
  if (ny < 2) return(0)

  r1 <- as.vector(resid_mat[, 2:ny, drop = FALSE])
  r0 <- as.vector(resid_mat[, 1:(ny - 1), drop = FALSE])
  rho <- suppressWarnings(stats::cor(r0, r1))
  if (!is.finite(rho)) rho <- 0
  rho <- max(-clamp, min(clamp, rho))
  rho
}

# CI extraction + scoring ----

coef_ci <- function(est, se, alpha = 0.10) {
  z <- stats::qnorm(1 - alpha / 2)
  list(
    lower = est - z * se,
    upper = est + z * se,
    width = 2 * z * se
  )
}

select_lp1_smooths <- function(fit, coef_obj, lp = 1L) {
  sm_list <- unname(coef_obj$smterms)
  sm_names <- names(fit$smooth)
  stopifnot(length(sm_list) == length(sm_names))

  is_multi_lp <- !is.null(fit$family$nlp) && fit$family$nlp > 1
  if (!is_multi_lp) {
    return(list(sm_list = sm_list, sm_names = sm_names))
  }

  keep <- grepl(paste0("\\.", lp, "(\\(|$)"), sm_names)
  list(sm_list = sm_list[keep], sm_names = sm_names[keep])
}

select_target_terms <- function(fit, coef_obj, term_types = NULL) {
  lp1 <- select_lp1_smooths(fit, coef_obj, lp = 1L)
  sm_list <- lp1$sm_list
  sm_names <- lp1$sm_names

  is_ff <- vapply(
    seq_along(sm_list),
    \(i) isTRUE(sm_list[[i]]$dim == 2) && grepl("X1", sm_names[[i]]),
    logical(1)
  )
  is_smooth_zt <- vapply(
    seq_along(sm_list),
    \(i) isTRUE(sm_list[[i]]$dim == 2) && grepl("zsmoo", sm_names[[i]]),
    logical(1)
  )
  is_concurrent <- vapply(
    seq_along(sm_list),
    \(i) isTRUE(sm_list[[i]]$dim == 1) && grepl("Xc", sm_names[[i]]),
    logical(1)
  )
  is_linear <- vapply(
    seq_along(sm_list),
    \(i) isTRUE(sm_list[[i]]$dim == 1) && grepl("zlin", sm_names[[i]]),
    logical(1)
  )

  pick_one <- function(idx, label) {
    if (sum(idx) != 1) {
      stop("Could not uniquely identify term for: ", label)
    }
    which(idx)
  }

  out <- list()
  if (is.null(term_types) || "ff" %in% term_types) {
    if (sum(is_ff) == 1) {
      out$ff <- list(
        sm = sm_list[[pick_one(is_ff, "ff")]],
        smooth_name = sm_names[[pick_one(is_ff, "ff")]]
      )
    }
  }
  if (is.null(term_types) || "concurrent" %in% term_types) {
    if (sum(is_concurrent) == 1) {
      out$concurrent <- list(
        sm = sm_list[[pick_one(is_concurrent, "concurrent Xc")]],
        smooth_name = sm_names[[pick_one(is_concurrent, "concurrent Xc")]]
      )
    }
  }
  if (is.null(term_types) || "linear" %in% term_types) {
    if (sum(is_linear) == 1) {
      out$linear <- list(
        sm = sm_list[[pick_one(is_linear, "linear zlin")]],
        smooth_name = sm_names[[pick_one(is_linear, "linear zlin")]]
      )
    }
  }
  if (is.null(term_types) || "smooth" %in% term_types) {
    if (sum(is_smooth_zt) == 1) {
      out$smooth <- list(
        sm = sm_list[[pick_one(is_smooth_zt, "smooth s(zsmoo)")]],
        smooth_name = sm_names[[pick_one(is_smooth_zt, "smooth s(zsmoo)")]]
      )
    }
  }
  out
}

score_one_term <- function(term, truth_vec, alpha = 0.10) {
  est <- term$sm$coef[, "value"]
  se <- term$sm$coef[, "se"]
  ci <- coef_ci(est = est, se = se, alpha = alpha)
  sq_err <- (est - truth_vec)^2

  data.frame(
    smooth_name = term$smooth_name,
    x = term$sm$coef[, 1],
    y = if (term$sm$dim == 2) term$sm$coef[, 2] else NA_real_,
    est = est,
    se = se,
    truth = truth_vec,
    lower = ci$lower,
    upper = ci$upper,
    width = ci$width,
    sq_err = sq_err,
    truth_sq = truth_vec^2,
    covered = (truth_vec >= ci$lower) & (truth_vec <= ci$upper)
  )
}

integrate_over_grid <- function(df, value_col) {
  stopifnot(is.data.frame(df), value_col %in% names(df))
  vals <- df[[value_col]]
  if (all(is.na(df$y))) {
    xg <- sort(unique(df$x))
    dx <- if (length(xg) > 1) mean(diff(xg)) else 1
    return(sum(vals) * dx)
  }
  xg <- sort(unique(df$x))
  yg <- sort(unique(df$y))
  dx <- if (length(xg) > 1) mean(diff(xg)) else 1
  dy <- if (length(yg) > 1) mean(diff(yg)) else 1
  sum(vals) * dx * dy
}

evaluate_fit <- function(
  fit,
  truth,
  alpha = 0.10,
  use_sandwich = FALSE,
  coef_grid = list(n1 = 80, n2 = 30, n3 = 20),
  quiet = TRUE,
  term_types = NULL
) {
  coefs <- if (quiet) {
    tmp <- NULL
    invisible(utils::capture.output(
      tmp <- coef(
        fit,
        sandwich = use_sandwich,
        n1 = coef_grid$n1,
        n2 = coef_grid$n2,
        n3 = coef_grid$n3
      )
    ))
    tmp
  } else {
    coef(
      fit,
      sandwich = use_sandwich,
      n1 = coef_grid$n1,
      n2 = coef_grid$n2,
      n3 = coef_grid$n3
    )
  }
  terms <- select_target_terms(fit, coefs, term_types = term_types)

  truth_funs <- truth$truth_funs
  if (is.null(truth_funs)) {
    stop("Truth functions missing from dataset attributes.")
  }

  out_list <- list()
  if (!is.null(terms$ff)) {
    w_s <- truth$center_ref$ff_weights %||% ff_weights(terms$ff$sm$x)
    ff_vals <- truth_funs$beta_ff(terms$ff$sm$x, terms$ff$sm$y)
    ff_vals <- center_beta_ff(ff_vals, terms$ff$sm$x, w_s)
    ff_truth <- as.vector(ff_vals)
    out_list$ff <- score_one_term(terms$ff, ff_truth, alpha = alpha)
    out_list$ff$term_type <- "ff"
  }
  if (!is.null(terms$concurrent)) {
    con_truth <- truth_funs$beta_concurrent(terms$concurrent$sm$x)
    out_list$concurrent <- score_one_term(
      terms$concurrent,
      con_truth,
      alpha = alpha
    )
    out_list$concurrent$term_type <- "concurrent"
  }
  if (!is.null(terms$linear)) {
    lin_truth <- truth_funs$beta_linear(terms$linear$sm$x)
    out_list$linear <- score_one_term(terms$linear, lin_truth, alpha = alpha)
    out_list$linear$term_type <- "linear"
  }
  if (!is.null(terms$smooth)) {
    sm_truth_raw <- as.vector(
      truth_funs$f_smooth_scalar(terms$smooth$sm$x, terms$smooth$sm$y)
    )

    # Center smooth truth per t to match pffr identifiability constraints:
    # sum_i f(z_i, t) = 0 for all t.
    zsmoo_obs <- truth$center_ref$zsmoo %||% numeric(0)
    if (length(zsmoo_obs)) {
      m_t <- colMeans(
        truth_funs$f_smooth_scalar(zsmoo_obs, terms$smooth$sm$y)
      )
      n_x <- length(unique(terms$smooth$sm$x))
      m_t_vec <- rep(m_t, each = n_x)
      sm_truth <- sm_truth_raw - m_t_vec
    } else {
      sm_truth <- sm_truth_raw
    }

    out_list$smooth <- score_one_term(terms$smooth, sm_truth, alpha = alpha)
    out_list$smooth$term_type <- "smooth"
  }
  do.call(rbind, out_list)
}

# Benchmark runner ----

make_dgp_grid <- function() {
  # Full factorial design for benchmark
  expand.grid(
    n = 80L,
    nxgrid = 35L,
    nygrid = 45L,
    snr = c(2, 8, 25),
    wiggliness = c(2, 20),
    error_dist = c("gaussian", "t6"),
    corr_type = c("iid", "ar1", "gauss", "fourier"),
    corr_param = NA_real_,
    hetero_type = c("none", "bump"),
    hetero_param = NA_real_,
    stringsAsFactors = FALSE
  )
}

make_dgp_grid_tiny <- function() {
  # Minimal grid for quick smoke tests.
  expand.grid(
    n = c(30L),
    nxgrid = c(25L),
    nygrid = c(30L),
    snr = c(10),
    wiggliness = c(2),
    error_dist = c("gaussian"),
    corr_type = c("iid", "ar1"),
    corr_param = c(NA_real_),
    hetero_type = c("none"),
    hetero_param = c(NA_real_),
    stringsAsFactors = FALSE
  )
}

make_dgp_grid_single_terms <- function() {
  # Test each term type separately at SNR=25 (near-nominal coverage expected)
  data.frame(
    n = 200L,
    nxgrid = 35L,
    nygrid = 45L,
    snr = 25,
    wiggliness = 2,
    error_dist = "gaussian",
    corr_type = "iid",
    corr_param = NA_real_,
    hetero_type = "none",
    hetero_param = NA_real_,
    term_set = c("ff", "concurrent", "linear", "smooth"),
    stringsAsFactors = FALSE
  )
}

make_benchmark_grid <- function(
  dgp_grid,
  n_rep,
  fit_families,
  methods,
  seed = 1L
) {
  stopifnot(is.data.frame(dgp_grid), nrow(dgp_grid) > 0)
  stopifnot(is.numeric(n_rep), length(n_rep) == 1, n_rep >= 1)

  dgp_grid$dgp_id <- seq_len(nrow(dgp_grid))
  rep_grid <- data.frame(rep_id = seq_len(n_rep))
  fit_grid <- expand.grid(
    method = methods,
    fit_family = fit_families,
    stringsAsFactors = FALSE
  )

  bench_grid <- merge(dgp_grid, rep_grid, by = NULL)
  bench_grid <- merge(bench_grid, fit_grid, by = NULL)
  bench_grid$sim_seed <- as.integer(seed) +
    1000L * bench_grid$dgp_id +
    bench_grid$rep_id

  bench_grid <- subset(
    bench_grid,
    !(fit_family == "gaulss" & hetero_type == "none") &
      !(fit_family == "scat" & error_dist != "t6") &
      !(method == "pffr_ar" & fit_family %in% c("gaulss", "scat")) &
      !(method %in% c("pffr_gls", "pffr_gls_est") & fit_family != "gaussian")
  )

  bench_grid$bench_id <- seq_len(nrow(bench_grid))
  rownames(bench_grid) <- NULL
  bench_grid
}

# Helper to load and combine saved benchmark results
load_benchmark_results <- function(output_dir, pattern = "row.*\\.rds$") {
  files <- list.files(output_dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) {
    stop("No result files found in: ", output_dir)
  }

  results <- lapply(files, readRDS)
  combined <- do.call(rbind, results)
  cat(
    "Loaded",
    length(files),
    "result files with",
    nrow(combined),
    "total rows\n"
  )
  combined
}

# Resume a benchmark from where it left off
resume_benchmark <- function(
  dgp_grid,
  output_dir,
  save_prefix = "benchmark",
  n_rep = 50L,
  fit_families = c("gaussian", "gaulss", "scat"),
  methods = c("pffr", "pffr_gls", "pffr_gls_est", "pffr_ar", "pffr_sandwich"),
  seed = 1L,
  ...
) {
  # Find which DGPs have already been completed
  existing_files <- list.files(
    output_dir,
    pattern = paste0(save_prefix, "_row.*\\.rds$")
  )
  completed_rows <- as.integer(gsub(
    ".*row([0-9]+)\\.rds$",
    "\\1",
    existing_files
  ))

  bench_grid <- make_benchmark_grid(
    dgp_grid = dgp_grid,
    n_rep = n_rep,
    fit_families = fit_families,
    methods = methods,
    seed = seed
  )

  if (length(completed_rows) > 0) {
    cat(
      "Found",
      length(completed_rows),
      "completed rows; resuming remaining settings.\n"
    )
    bench_grid <- bench_grid[
      !(bench_grid$bench_id %in% completed_rows),
      ,
      drop = FALSE
    ]
  }

  if (nrow(bench_grid) == 0) {
    cat("All settings already completed. Loading existing results.\n")
    return(load_benchmark_results(output_dir, paste0(save_prefix, "_row")))
  }

  run_benchmark(
    dgp_grid,
    n_rep = n_rep,
    fit_families = fit_families,
    methods = methods,
    seed = seed,
    bench_grid = bench_grid,
    output_dir = output_dir,
    save_prefix = save_prefix,
    ...
  )
}

run_benchmark <- function(
  dgp_grid,
  n_rep = 50L,
  fit_families = c("gaussian", "gaulss", "scat"),
  methods = c("pffr", "pffr_gls", "pffr_gls_est", "pffr_ar", "pffr_sandwich"),
  alpha = 0.10, # 90% CI
  seed = 1L,
  coef_grid = list(n1 = 80, n2 = 30, n3 = 20),
  parallel = TRUE,
  n_workers = NULL,
  fit_args = NULL, # NULL means use default_fit_args()
  bench_grid = NULL,
  output_dir = NULL,
  save_prefix = "benchmark"
) {
  require_pkgs()
  if (is.null(bench_grid)) {
    stopifnot(is.data.frame(dgp_grid), nrow(dgp_grid) > 0)
    stopifnot(is.numeric(n_rep), length(n_rep) == 1, n_rep >= 1)
  } else {
    stopifnot(is.data.frame(bench_grid), nrow(bench_grid) > 0)
  }

  # Use default fit_args if not specified
  if (is.null(fit_args)) {
    fit_args <- default_fit_args()
  } else {
    # Merge with defaults (user args take precedence)
    fit_args <- utils::modifyList(default_fit_args(), fit_args)
  }

  if (is.null(bench_grid)) {
    bench_grid <- make_benchmark_grid(
      dgp_grid = dgp_grid,
      n_rep = n_rep,
      fit_families = fit_families,
      methods = methods,
      seed = seed
    )
  }

  # Setup output directory for incremental saving
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    cat("Results will be saved incrementally to:", output_dir, "\n")
  }

  merge_args <- function(base_args, extra_args) {
    # `base_args` should win over `extra_args`
    # Remove k_ff and k_smooth which are only for formula construction
    extra_args <- extra_args[!names(extra_args) %in% c("k_ff", "k_smooth")]
    utils::modifyList(extra_args, base_args)
  }

  add_meta <- function(scored, row) {
    scored$bench_id <- row$bench_id
    scored$dgp_id <- row$dgp_id
    scored$rep_id <- row$rep_id
    scored$sim_seed <- row$sim_seed
    scored$method <- row$method
    scored$fit_family <- row$fit_family
    scored$snr <- row$snr
    scored$wiggliness <- row$wiggliness
    scored$error_dist <- row$error_dist
    scored$corr_type <- row$corr_type
    scored$hetero_type <- row$hetero_type
    scored
  }

  run_one_setting <- function(row) {
    finalize_scored <- function(scored) {
      log_line <- paste(
        utils::capture.output(print(utils::tail(scored, 1))),
        collapse = "\n"
      )
      list(result = scored, log = log_line)
    }

    term_set <- if ("term_set" %in% names(row)) {
      strsplit(as.character(row$term_set), "\\s*,\\s*")[[1]]
    } else {
      c("ff", "concurrent", "linear", "smooth")
    }

    sim <- simulate_dataset(
      n = row$n,
      nxgrid = row$nxgrid,
      nygrid = row$nygrid,
      snr = row$snr,
      wiggliness = row$wiggliness,
      error_dist = row$error_dist,
      corr_type = row$corr_type,
      corr_param = if (is.na(row$corr_param)) NULL else row$corr_param,
      hetero_type = row$hetero_type,
      hetero_param = if (is.na(row$hetero_param)) NULL else row$hetero_param,
      seed = row$sim_seed,
      terms = term_set
    )
    dat <- sim$data
    s_grid <- sim$s_grid
    t_grid <- sim$t_grid
    hatSigma_true <- sim$err$cov_mat
    truth <- attr(dat, "truth")

    k_ff <- fit_args$k_ff %||% 8
    k_smooth <- fit_args$k_smooth %||% 12
    ff_splinepars <- sprintf(
      "list(bs = 'ps', m = list(c(2, 1), c(2, 1)), k = c(%d, %d))",
      k_ff,
      k_ff
    )
    formula_terms <- c(
      if ("ff" %in% term_set)
        sprintf("ff(X1, xind = s_grid, splinepars = %s)", ff_splinepars) else
        NULL,
      if ("concurrent" %in% term_set) "Xc" else NULL,
      if ("linear" %in% term_set) "zlin" else NULL,
      if ("smooth" %in% term_set) sprintf("s(zsmoo, k = %d)", k_smooth) else
        NULL
    )
    formula <- if (length(formula_terms)) {
      as.formula(paste("Y ~", paste(formula_terms, collapse = " + ")))
    } else {
      Y ~ 1
    }

    family <- fit_family_object(row$fit_family)
    method <- row$method

    need_std <- method %in%
      c("pffr", "pffr_sandwich", "pffr_gls_est", "pffr_ar")
    fit_std <- NULL
    if (need_std) {
      fit_std <- do.call(
        "pffr",
        merge_args(
          list(
            formula = formula,
            yind = t_grid,
            data = dat,
            family = family
          ),
          fit_args
        )
      )
    }

    if (method == "pffr") {
      scored <- evaluate_fit(
        fit = fit_std,
        truth = truth,
        alpha = alpha,
        use_sandwich = FALSE,
        coef_grid = coef_grid,
        term_types = term_set
      )
      return(finalize_scored(add_meta(scored, row)))
    }

    if (method == "pffr_sandwich") {
      scored <- evaluate_fit(
        fit = fit_std,
        truth = truth,
        alpha = alpha,
        use_sandwich = TRUE,
        coef_grid = coef_grid,
        term_types = term_set
      )
      return(finalize_scored(add_meta(scored, row)))
    }

    if (method == "pffr_gls_est") {
      hatSigma_est <- estimate_hatSigma_from_fit(fit_std, dat)
      fit_gls_est <- do.call(
        "pffr_gls",
        merge_args(
          list(
            formula = formula,
            yind = t_grid,
            data = dat,
            hatSigma = hatSigma_est
          ),
          fit_args
        )
      )
      scored <- evaluate_fit(
        fit = fit_gls_est,
        truth = truth,
        alpha = alpha,
        use_sandwich = FALSE,
        coef_grid = coef_grid,
        term_types = term_set
      )
      return(finalize_scored(add_meta(scored, row)))
    }

    if (method == "pffr_gls") {
      fit_gls <- do.call(
        "pffr_gls",
        merge_args(
          list(
            formula = formula,
            yind = t_grid,
            data = dat,
            hatSigma = hatSigma_true
          ),
          fit_args
        )
      )
      scored <- evaluate_fit(
        fit = fit_gls,
        truth = truth,
        alpha = alpha,
        use_sandwich = FALSE,
        coef_grid = coef_grid,
        term_types = term_set
      )
      return(finalize_scored(add_meta(scored, row)))
    }

    if (method == "pffr_ar") {
      rho_est <- estimate_rho_from_fit(fit_std, dat)
      fit_ar <- do.call(
        "pffr",
        merge_args(
          list(
            formula = formula,
            yind = t_grid,
            data = dat,
            family = family,
            algorithm = "bam",
            rho = rho_est
          ),
          fit_args
        )
      )
      scored <- evaluate_fit(
        fit = fit_ar,
        truth = truth,
        alpha = alpha,
        use_sandwich = FALSE,
        coef_grid = coef_grid,
        term_types = term_set
      )
      return(finalize_scored(add_meta(scored, row)))
    }

    stop("Unknown method: ", method)
  }

  if (
    parallel &&
      requireNamespace("foreach", quietly = TRUE) &&
      requireNamespace("doRNG", quietly = TRUE)
  ) {
    suppressPackageStartupMessages(library(foreach))
    suppressPackageStartupMessages(library(doRNG))

    if (!requireNamespace("doParallel", quietly = TRUE)) {
      foreach::registerDoSEQ()
    } else {
      if (is.null(n_workers)) {
        n_workers <- max(1L, parallel::detectCores() - 1L)
      }
      if (n_workers <= 1L) {
        foreach::registerDoSEQ()
      } else {
        doParallel::registerDoParallel(cores = n_workers)
      }
    }

    # Reproducible parallel RNG streams.
    doRNG::registerDoRNG(seed)

    combine_with_log <- function(acc, x) {
      if (!is.null(x$log)) {
        cat(x$log, "\n")
        flush.console()
      }
      if (is.null(acc)) return(x$result)
      rbind(acc, x$result)
    }

    combined_results <- foreach::foreach(
      row = seq_len(nrow(bench_grid)),
      .combine = combine_with_log,
      .init = NULL,
      .inorder = FALSE,
      .packages = c("mgcv", "mvtnorm", "splines"),
      .export = c(
        "run_one_setting",
        "bench_grid",
        "fit_args",
        "alpha",
        "coef_grid",
        "require_pkgs",
        "load_refund_dev",
        "ensure_pffr_simulate_bench",
        "simulate_dataset",
        "pffr_simulate_bench",
        "make_truth_functions",
        "make_error_structure",
        "sample_errors",
        "compute_snr_scale_factor",
        "rescale_errors_to_snr",
        "random_function_matrix",
        "ff_weights",
        "center_beta_ff",
        "compute_ff_effect",
        "compute_linear_effect",
        "fit_family_object",
        "estimate_hatSigma_from_fit",
        "estimate_rho_from_fit",
        "fitted_mean_matrix",
        "evaluate_fit",
        "select_lp1_smooths",
        "select_target_terms",
        "score_one_term",
        "coef_ci",
        "integrate_over_grid",
        "make_bs_spec",
        "eval_bs",
        "%||%"
      )
    ) %dorng%
      {
        require_pkgs() # loads refund (pffr, pffr_gls, etc.)
        run_one_setting(bench_grid[row, , drop = FALSE])
      }
    if (!is.null(output_dir) && !is.null(combined_results)) {
      combined_path <- file.path(
        output_dir,
        paste0(save_prefix, "_combined.rds")
      )
      saveRDS(combined_results, combined_path)
      cat(sprintf(
        "  Saved combined results: %s (%d rows)\n",
        combined_path,
        nrow(combined_results)
      ))
    }
    combined_results
  } else {
    set.seed(seed)
    all_results <- list()

    for (row in seq_len(nrow(bench_grid))) {
      dgp <- bench_grid[row, , drop = FALSE]
      cat(sprintf(
        "\n[%d/%d] Running DGP: SNR=%s, wiggle=%s, corr=%s, hetero=%s, error=%s\n",
        row,
        nrow(bench_grid),
        dgp$snr,
        dgp$wiggliness,
        dgp$corr_type,
        dgp$hetero_type,
        dgp$error_dist
      ))

      result_obj <- tryCatch(
        run_one_setting(dgp),
        error = function(e) {
          warning("Error in DGP row ", row, ": ", e$message)
          NULL
        }
      )

      if (!is.null(result_obj)) {
        if (!is.null(result_obj$log)) {
          cat(result_obj$log, "\n")
          flush.console()
        }
        result <- result_obj$result
        all_results[[row]] <- result

        # Save incrementally if output_dir specified
        if (!is.null(output_dir)) {
          dgp_label <- sprintf("%s_row%05d", save_prefix, dgp$bench_id)
          save_path <- file.path(output_dir, paste0(dgp_label, ".rds"))
          saveRDS(result, save_path)

          # Also save combined results so far
          combined <- do.call(rbind, all_results[!sapply(all_results, is.null)])
          combined_path <- file.path(
            output_dir,
            paste0(save_prefix, "_combined.rds")
          )
          saveRDS(combined, combined_path)

          cat(sprintf(
            "  Saved: %s (combined: %d rows)\n",
            save_path,
            nrow(combined)
          ))
        }
      }
    }

    do.call(rbind, all_results[!sapply(all_results, is.null)])
  }
}

summarize_results <- function(scored) {
  stopifnot(is.data.frame(scored), nrow(scored) > 0)
  agg_vars <- c(
    "method",
    "fit_family",
    "term_type",
    "snr",
    "wiggliness",
    "error_dist",
    "corr_type",
    "hetero_type"
  )

  # First summarize per replicate (so integrated metrics are well-defined),
  # then average across replicates.
  agg_rep_vars <- c(agg_vars, "rep_id")
  key_rep <- scored[, agg_rep_vars, drop = FALSE]

  summarize_one_rep <- function(df) {
    int_sq_err <- integrate_over_grid(df, "sq_err")
    int_truth_sq <- integrate_over_grid(df, "truth_sq")
    int_rel_rmse <- sqrt(int_sq_err) / (sqrt(int_truth_sq) + 1e-12)

    data.frame(
      coverage = mean(df$covered),
      width = mean(df$width),
      int_rel_rmse = int_rel_rmse,
      n_grid = nrow(df)
    )
  }

  split_idx_rep <- interaction(key_rep, drop = TRUE)
  parts_rep <- split(scored, split_idx_rep)
  sum_rep_list <- lapply(parts_rep, summarize_one_rep)
  sum_rep_df <- do.call(rbind, sum_rep_list)
  keys_rep_unique <- key_rep[
    match(names(parts_rep), split_idx_rep),
    ,
    drop = FALSE
  ]
  rownames(keys_rep_unique) <- NULL
  per_rep <- cbind(keys_rep_unique, sum_rep_df)

  # Aggregate across replicates
  split_idx <- interaction(per_rep[, agg_vars, drop = FALSE], drop = TRUE)
  parts <- split(per_rep, split_idx)
  sum_list <- lapply(
    parts,
    \(df)
      data.frame(
        coverage = mean(df$coverage),
        width = mean(df$width),
        int_rel_rmse = mean(df$int_rel_rmse),
        n_rep = nrow(df),
        n_grid = df$n_grid[1]
      )
  )
  sum_df <- do.call(rbind, sum_list)
  keys_unique <- per_rep[match(names(parts), split_idx), agg_vars, drop = FALSE]
  rownames(keys_unique) <- NULL

  cbind(keys_unique, sum_df)
}

# Visualization functions ----

plot_summary <- function(summary_df, nominal = 0.90) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("ggplot2 not installed; skipping plots.")
    return(NULL)
  }
  ggplot2::ggplot(
    summary_df,
    ggplot2::aes(x = method, y = coverage, color = fit_family)
  ) +
    ggplot2::geom_hline(yintercept = nominal, linetype = 2, color = "gray40") +
    ggplot2::geom_hline(yintercept = 0.90, linetype = 3, color = "gray60") +
    ggplot2::geom_point(
      position = ggplot2::position_dodge(width = 0.4),
      size = 2
    ) +
    ggplot2::facet_grid(
      term_type ~ error_dist + corr_type + hetero_type + snr,
      labeller = ggplot2::label_both
    ) +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
    ggplot2::labs(
      title = "CI Coverage by Term Type and DGP Settings",
      y = "Empirical Coverage",
      x = "Method"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      strip.text = ggplot2::element_text(size = 8)
    )
}

plot_metric <- function(
  summary_df,
  metric = c("coverage", "width", "int_rel_rmse"),
  nominal = 0.90
) {
  metric <- match.arg(metric)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("ggplot2 not installed; skipping plots.")
    return(NULL)
  }

  p <- ggplot2::ggplot(
    summary_df,
    ggplot2::aes(x = method, y = .data[[metric]], color = fit_family)
  ) +
    ggplot2::geom_point(
      position = ggplot2::position_dodge(width = 0.4),
      size = 2
    ) +
    ggplot2::facet_grid(
      term_type ~ error_dist + corr_type + hetero_type + snr,
      labeller = ggplot2::label_both
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      strip.text = ggplot2::element_text(size = 8)
    )

  if (metric == "coverage") {
    p <- p +
      ggplot2::geom_hline(
        yintercept = nominal,
        linetype = 2,
        color = "gray40"
      ) +
      ggplot2::geom_hline(yintercept = 0.90, linetype = 3, color = "gray60") +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::labs(y = "Coverage", title = "CI Coverage")
  } else if (metric == "width") {
    p <- p + ggplot2::labs(y = "Mean CI Width", title = "CI Width")
  } else {
    p <- p +
      ggplot2::labs(
        y = "Integrated Relative RMSE",
        title = "Estimation Accuracy"
      )
  }
  p
}

plot_surface_metric <- function(
  scored,
  term_type = c("ff", "smooth"),
  metric = "covered",
  show_delta = FALSE,
  nominal = 0.90
) {
  term_type <- match.arg(term_type)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("ggplot2 not installed; skipping plots.")
    return(NULL)
  }

  df <- scored[scored$term_type == term_type, , drop = FALSE]
  if (nrow(df) == 0) {
    message("No data for term_type: ", term_type)
    return(NULL)
  }

  agg <- aggregate(
    df[[metric]],
    by = list(x = df$x, y = df$y, method = df$method),
    FUN = mean
  )
  names(agg)[4] <- metric

  if (show_delta && metric == "covered") {
    agg$fill_val <- agg[[metric]] - nominal
    fill_label <- paste0("Coverage - ", nominal)
    fill_limits <- c(-0.5, 0.5)
    fill_colors <- c("darkred", "white", "darkblue")
  } else {
    agg$fill_val <- agg[[metric]]
    fill_label <- if (metric == "covered") "Coverage" else metric
    fill_limits <- if (metric == "covered") c(0, 1) else NULL
    fill_colors <- NULL
  }

  p <- ggplot2::ggplot(agg, ggplot2::aes(x = x, y = y, fill = fill_val)) +
    ggplot2::geom_tile() +
    ggplot2::facet_wrap(~method) +
    ggplot2::labs(
      title = sprintf("%s: %s by Grid Position", term_type, fill_label),
      x = if (term_type == "ff") "s" else "z",
      y = "t",
      fill = fill_label
    ) +
    ggplot2::theme_bw() +
    ggplot2::coord_fixed()

  if (!is.null(fill_colors)) {
    p <- p +
      ggplot2::scale_fill_gradient2(
        low = fill_colors[1],
        mid = fill_colors[2],
        high = fill_colors[3],
        midpoint = 0,
        limits = fill_limits
      )
  } else if (!is.null(fill_limits)) {
    p <- p + ggplot2::scale_fill_viridis_c(limits = fill_limits)
  } else {
    p <- p + ggplot2::scale_fill_viridis_c()
  }
  p
}

plot_curve_metric <- function(
  scored,
  term_type = c("concurrent", "linear"),
  metric = "covered",
  show_ribbon = TRUE,
  nominal = 0.90
) {
  term_type <- match.arg(term_type)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("ggplot2 not installed; skipping plots.")
    return(NULL)
  }

  df <- scored[scored$term_type == term_type, , drop = FALSE]
  if (nrow(df) == 0) {
    message("No data for term_type: ", term_type)
    return(NULL)
  }

  # Aggregate mean and SE
  agg <- aggregate(
    df[[metric]],
    by = list(x = df$x, method = df$method),
    FUN = function(v) c(mean = mean(v), se = sd(v) / sqrt(length(v)))
  )
  agg <- cbind(agg[, 1:2], as.data.frame(agg[[3]]))
  names(agg)[3:4] <- c("mean_val", "se_val")

  p <- ggplot2::ggplot(agg, ggplot2::aes(x = x, y = mean_val, color = method))

  if (show_ribbon) {
    p <- p +
      ggplot2::geom_ribbon(
        ggplot2::aes(
          ymin = mean_val - 1.96 * se_val,
          ymax = mean_val + 1.96 * se_val,
          fill = method
        ),
        alpha = 0.2,
        color = NA
      )
  }

  p <- p +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::labs(
      title = sprintf("%s: %s vs t", term_type, metric),
      x = "t",
      y = if (metric == "covered") "Coverage" else metric
    ) +
    ggplot2::theme_bw()

  if (metric == "covered") {
    p <- p +
      ggplot2::geom_hline(
        yintercept = nominal,
        linetype = 2,
        color = "gray40"
      ) +
      ggplot2::scale_y_continuous(limits = c(0, 1))
  }
  p
}

# Z-score distribution plot
plot_zscore_dist <- function(scored, facet_by = "term_type") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("ggplot2 not installed; skipping plots.")
    return(NULL)
  }

  scored$z_score <- (scored$est - scored$truth) / scored$se

  p <- ggplot2::ggplot(scored, ggplot2::aes(x = z_score)) +
    ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(density)),
      bins = 50,
      fill = "steelblue",
      alpha = 0.7
    ) +
    ggplot2::stat_function(fun = dnorm, color = "red", linewidth = 1) +
    ggplot2::facet_wrap(as.formula(paste("~", facet_by))) +
    ggplot2::labs(
      title = "Z-score Distribution vs Standard Normal",
      x = "Z-score = (estimate - truth) / SE",
      y = "Density"
    ) +
    ggplot2::theme_bw() +
    ggplot2::xlim(-5, 5)
  p
}

# Estimate vs Truth comparison plot
plot_estimate_vs_truth <- function(
  scored,
  term_type,
  method = NULL,
  rep_id = 1
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("ggplot2 not installed; skipping plots.")
    return(NULL)
  }

  df <- scored[scored$term_type == term_type, , drop = FALSE]
  if (!is.null(method)) df <- df[df$method == method, , drop = FALSE]
  df <- df[df$rep_id == rep_id, , drop = FALSE]

  if (nrow(df) == 0) {
    message("No data for specified filters")
    return(NULL)
  }

  is_2d <- !all(is.na(df$y))

  if (is_2d) {
    # For 2D terms (ff, smooth), show heatmaps side by side
    df_long <- rbind(
      data.frame(x = df$x, y = df$y, value = df$truth, type = "Truth"),
      data.frame(x = df$x, y = df$y, value = df$est, type = "Estimate")
    )
    p <- ggplot2::ggplot(df_long, ggplot2::aes(x = x, y = y, fill = value)) +
      ggplot2::geom_tile() +
      ggplot2::facet_wrap(~type) +
      ggplot2::scale_fill_viridis_c() +
      ggplot2::coord_fixed() +
      ggplot2::labs(title = sprintf("%s: Estimate vs Truth", term_type)) +
      ggplot2::theme_bw()
  } else {
    # For 1D terms, show lines with CI ribbon
    p <- ggplot2::ggplot(df, ggplot2::aes(x = x)) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = lower, ymax = upper),
        fill = "steelblue",
        alpha = 0.3
      ) +
      ggplot2::geom_line(
        ggplot2::aes(y = est, color = "Estimate"),
        linewidth = 1
      ) +
      ggplot2::geom_line(
        ggplot2::aes(y = truth, color = "Truth"),
        linewidth = 1,
        linetype = 2
      ) +
      ggplot2::scale_color_manual(
        values = c("Estimate" = "steelblue", "Truth" = "red")
      ) +
      ggplot2::labs(
        title = sprintf("%s: Estimate vs Truth (rep %d)", term_type, rep_id),
        x = "t",
        y = "Coefficient",
        color = ""
      ) +
      ggplot2::theme_bw()
  }
  p
}

# Coverage vs SNR plot (useful for showing degradation at high SNR)
plot_coverage_vs_snr <- function(summary_df, nominal = 0.90) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("ggplot2 not installed; skipping plots.")
    return(NULL)
  }

  ggplot2::ggplot(
    summary_df,
    ggplot2::aes(
      x = factor(snr),
      y = coverage,
      color = method,
      shape = fit_family
    )
  ) +
    ggplot2::geom_hline(yintercept = nominal, linetype = 2, color = "gray40") +
    ggplot2::geom_hline(yintercept = 0.90, linetype = 3, color = "gray60") +
    ggplot2::geom_point(
      size = 3,
      position = ggplot2::position_dodge(width = 0.3)
    ) +
    ggplot2::facet_wrap(~term_type) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      title = "Coverage vs SNR by Term Type",
      x = "SNR",
      y = "Coverage"
    ) +
    ggplot2::theme_bw()
}

# Multi-panel summary for a single DGP
plot_dgp_summary <- function(scored, summary_df, dgp_label = "") {
  if (
    !requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("patchwork", quietly = TRUE)
  ) {
    message("ggplot2 and/or patchwork not installed; skipping plots.")
    return(NULL)
  }

  p1 <- plot_metric(summary_df, "coverage") +
    ggplot2::ggtitle(paste("Coverage", dgp_label))
  p2 <- plot_zscore_dist(scored) +
    ggplot2::ggtitle("Z-score Distribution")

  patchwork::wrap_plots(p1, p2, ncol = 1)
}

# Main benchmark run ----
# Runs automatically when sourced. Results saved to ./ci-benchmark/

require_pkgs()

output_dir <- file.path(getwd(), "ci-benchmark")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

dgp_grid <- make_dgp_grid()
bench_grid <- make_benchmark_grid(
  dgp_grid = dgp_grid,
  n_rep = 10L,
  fit_families = c("gaussian", "gaulss", "scat"),
  methods = c("pffr", "pffr_gls", "pffr_gls_est", "pffr_ar", "pffr_sandwich"),
  seed = 2024L
)

cat("\n")
cat("==============================================================\n")
cat("  CI COVERAGE BENCHMARK\n")
cat("==============================================================\n")
cat("  DGP grid:   ", nrow(dgp_grid), " settings (make_dgp_grid)\n", sep = "")
cat("  Reps/DGP:    10\n")
cat("  Benchmark:  ", nrow(bench_grid), " fits\n", sep = "")
cat("  Methods:     pffr, pffr_gls, pffr_gls_est, pffr_ar, pffr_sandwich\n")
cat("  Families:    gaussian (+gaulss if hetero, +scat if t6)\n")
cat("  Parallel:    3 cores\n")
cat("  Output:     ", output_dir, "\n")
cat("==============================================================\n\n")

start_time <- Sys.time()
cat("Started at:", format(start_time), "\n\n")

scored <- run_benchmark(
  dgp_grid,
  n_rep = 10L,
  fit_families = c("gaussian", "gaulss", "scat"),
  methods = c("pffr", "pffr_gls", "pffr_gls_est", "pffr_ar", "pffr_sandwich"),
  alpha = 0.10,
  seed = 2024L,
  parallel = TRUE,
  n_workers = 3L,
  bench_grid = bench_grid,
  output_dir = output_dir,
  save_prefix = "ci_bench"
)

end_time <- Sys.time()
elapsed <- difftime(end_time, start_time, units = "mins")

cat("\n")
cat("==============================================================\n")
cat("  BENCHMARK COMPLETE\n")
cat("  Elapsed time:", round(as.numeric(elapsed), 1), "minutes\n")
cat("==============================================================\n\n")

# Summarize and save
summary_df <- summarize_results(scored)

saveRDS(scored, file.path(output_dir, "scored_full.rds"))
saveRDS(summary_df, file.path(output_dir, "summary.rds"))
write.csv(summary_df, file.path(output_dir, "summary.csv"), row.names = FALSE)

cat("Results saved to:\n")
cat("  ", file.path(output_dir, "scored_full.rds"), "\n")
cat("  ", file.path(output_dir, "summary.rds"), "\n")
cat("  ", file.path(output_dir, "summary.csv"), "\n\n")

# Print summary
cat("Coverage summary by term type:\n")
print(aggregate(coverage ~ term_type + method, data = summary_df, FUN = mean))

cat("\nDone.\n")
