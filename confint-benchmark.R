# confint-benchmark.R
#
# Benchmark CI coverage/width from coef.pffr() under different error structures.
#
# Methods tested:
#   - pffr: standard penalized functional regression
#   - pffr_sandwich: pffr with sandwich SEs
#   - pffr_gls: oracle covariance (upper bound)
#   - pffr_gls_est: estimated covariance from residuals
#   - pffr_ar: AR(1) residuals (bam only)
#
# Key design: parallelize over (dgp x rep), not over methods.
# For each dataset, fit pffr once and reuse for methods that need it.

# Setup ----

`%||%` <- function(x, y) if (is.null(x)) y else x

load_refund_dev <- function(pkg_dir = ".") {
  if (requireNamespace("devtools", quietly = TRUE) &&
      file.exists(file.path(pkg_dir, "DESCRIPTION"))) {
    devtools::load_all(pkg_dir, quiet = TRUE)
  } else {
    suppressPackageStartupMessages(library(refund))
  }
  invisible(TRUE)
}

require_pkgs <- function() {
  load_refund_dev(".")
  stopifnot(
    requireNamespace("mgcv", quietly = TRUE),
    requireNamespace("mvtnorm", quietly = TRUE)
  )
  invisible(TRUE)
}

# Truth functions ----

make_bs_spec <- function(x, df) {
  basis <- splines::bs(x, df = df, intercept = TRUE)
  list(
    knots = attr(basis, "knots"),
    boundary_knots = attr(basis, "Boundary.knots"),
    degree = attr(basis, "degree"),
    df = df
  )
}

eval_bs <- function(x, spec) {
  splines::bs(x, df = spec$df, knots = spec$knots,
              Boundary.knots = spec$boundary_knots,
              degree = spec$degree, intercept = TRUE)
}

make_truth_functions <- function(s_grid, t_grid, z_ref, wiggliness, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  rw_coef <- function(k, sd) {
    coef <- cumsum(rnorm(k, sd = sd))
    coef - mean(coef)
  }

  rw_coef_matrix <- function(k, sd) {
    mat <- matrix(rnorm(k * k, sd = sd), nrow = k)
    mat <- apply(mat, 2, cumsum)
    mat <- t(apply(mat, 1, cumsum))
    mat - mean(mat)
  }

  k_te <- 5L
  k_uni <- 7L

  s_spec <- make_bs_spec(s_grid, df = k_te)
  t_te_spec <- make_bs_spec(t_grid, df = k_te)
  t_uni_spec <- make_bs_spec(t_grid, df = k_uni)
  z_ref <- if (length(z_ref) < 2) c(-0.5, 0.5) else as.numeric(z_ref)
  z_te_spec <- make_bs_spec(z_ref, df = k_te)

  coef_ff <- rw_coef_matrix(k_te, wiggliness)
  coef_con <- rw_coef(k_uni, wiggliness)
  coef_lin <- rw_coef(k_uni, wiggliness)
  coef_mu <- rw_coef(k_uni, wiggliness)
  coef_smooth <- rw_coef_matrix(k_te, wiggliness)

  list(
    beta_ff = function(s, t) drop(eval_bs(s, s_spec) %*% coef_ff %*% t(eval_bs(t, t_te_spec))),
    beta_concurrent = function(t) drop(eval_bs(t, t_uni_spec) %*% coef_con),
    beta_linear = function(t) drop(eval_bs(t, t_uni_spec) %*% coef_lin),
    mu_t = function(t) drop(eval_bs(t, t_uni_spec) %*% coef_mu),
    f_smooth_scalar = function(z, t) drop(eval_bs(z, z_te_spec) %*% coef_smooth %*% t(eval_bs(t, t_te_spec)))
  )
}

# Error structure ----

make_error_structure <- function(t_grid, corr_type = "iid", hetero_type = "none") {
  ny <- length(t_grid)

  # Correlation matrix
  corr_mat <- switch(corr_type,
    iid = diag(1, ny),
    ar1 = {
      rho <- 0.4
      rho^abs(outer(seq_len(ny), seq_len(ny), "-"))
    },
    gauss = {
      dist_mat <- abs(outer(t_grid, t_grid, "-"))
      exp(-(dist_mat^2) / (2 * 0.15^2))
    },
    fourier = {
      B <- cbind(sin(2 * pi * t_grid / 0.33), cos(2 * pi * t_grid / 0.33),
                 sin(2 * pi * t_grid / 0.17), cos(2 * pi * t_grid / 0.17))
      B <- scale(B, center = TRUE, scale = FALSE)
      base_cov <- B %*% t(B) + diag(0.25, ny)
      stats::cov2cor(base_cov)
    }
  )

  # Heteroskedasticity
  sigma_t <- switch(hetero_type,
    none = rep(1, ny),
    bump = 1 + exp(-0.5 * ((t_grid - 0.7) / 0.10)^2)
  )

  cov_mat <- diag(sigma_t) %*% corr_mat %*% diag(sigma_t)
  # Ensure positive definiteness
  cov_mat <- cov_mat + diag(1e-8 * mean(diag(cov_mat)), ny)

  list(cov_mat = cov_mat, rho = if (corr_type == "ar1") 0.4 else NA_real_)
}

sample_errors <- function(n, cov_mat, error_dist = "gaussian") {
  ny <- nrow(cov_mat)
  if (error_dist == "gaussian") {
    return(mvtnorm::rmvnorm(n, mean = rep(0, ny), sigma = cov_mat))
  }
  df <- 6
  sigma_scale <- cov_mat * (df - 2) / df
  mvtnorm::rmvt(n, sigma = sigma_scale, df = df, delta = rep(0, ny))
}

# Data simulation ----

ff_weights <- function(xind) {
  nx <- length(xind)
  if (nx < 2) return(1)
  ((xind[nx] - xind[1]) / (nx - 1)) / 3 * c(1, rep(c(4, 2), length = nx - 2), 1)
}

center_beta_ff <- function(beta_st, s_grid, weights = NULL) {
  w <- weights %||% ff_weights(s_grid)
  col_means <- colSums(beta_st * w) / sum(w)
  sweep(beta_st, 2, col_means, "-")
}

random_curves <- function(n, grid, bs_dim = 25L) {
  X <- splines::bs(grid, df = bs_dim, intercept = TRUE)
  coef_mat <- matrix(rnorm(n * bs_dim), nrow = bs_dim)
  t(X %*% coef_mat)
}

simulate_dataset <- function(dgp, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n <- dgp$n
  s_grid <- seq(0, 1, length.out = dgp$nxgrid)
  t_grid <- seq(0, 1, length.out = dgp$nygrid)

  # Generate covariates (centered)
  X1 <- random_curves(n, s_grid, bs_dim = 20L)
  X1 <- sweep(X1, 2, colMeans(X1), "-")

  Xc <- random_curves(n, t_grid, bs_dim = 9L)
  Xc <- sweep(Xc, 2, colMeans(Xc), "-")

  zlin <- rnorm(n)
  zlin <- zlin - mean(zlin)

  zsmoo <- rnorm(n)
  zsmoo <- zsmoo - mean(zsmoo)

  # Truth functions
  truth_funs <- make_truth_functions(s_grid, t_grid, zsmoo, dgp$wiggliness,
                                     seed = if (!is.null(seed)) seed + 100000L else NULL)

  # Compute true signal components
  w_s <- ff_weights(s_grid)
  beta_st <- center_beta_ff(truth_funs$beta_ff(s_grid, t_grid), s_grid, w_s)
  beta_lin <- truth_funs$beta_linear(t_grid)
  beta_con <- truth_funs$beta_concurrent(t_grid)
  mu_t <- truth_funs$mu_t(t_grid)

  # eta components
  eta_int <- matrix(rep(mu_t, each = n), nrow = n)
  eta_ff <- (sweep(X1, 2, w_s, "*")) %*% beta_st
  eta_lin <- outer(zlin, beta_lin)
  eta_con <- Xc * rep(beta_con, each = n)
  eta_smooth_raw <- truth_funs$f_smooth_scalar(zsmoo, t_grid)
  eta_smooth <- sweep(eta_smooth_raw, 2, colMeans(eta_smooth_raw), "-")

  eta_total <- eta_int + eta_ff + eta_lin + eta_con + eta_smooth

  # Structured errors, rescaled to target SNR
  err <- make_error_structure(t_grid, dgp$corr_type, dgp$hetero_type)
  if (!is.null(seed)) set.seed(seed + 200000L)
  eps <- sample_errors(n, err$cov_mat, dgp$error_dist)

  var_eta <- var(as.vector(eta_total))
  var_eps <- var(as.vector(eps))
  snr_scale <- sqrt(var_eta / (dgp$snr * var_eps))
  eps <- eps * snr_scale
  err$cov_mat <- err$cov_mat * snr_scale^2

  Y <- eta_total + eps

  # Build data frame
  dat <- data.frame(zlin = zlin, zsmoo = zsmoo)
  dat$X1 <- I(X1)
  dat$Xc <- I(Xc)
  dat$Y <- I(Y)

  list(
    data = dat,
    s_grid = s_grid,
    t_grid = t_grid,
    cov_mat = err$cov_mat,
    rho = err$rho,
    truth = list(
      truth_funs = truth_funs,
      beta_st = beta_st,
      beta_lin = beta_lin,
      beta_con = beta_con,
      mu_t = mu_t,
      ff_weights = w_s,
      zsmoo = zsmoo
    )
  )
}

# Fitting ----

fit_pffr <- function(dat, s_grid, t_grid, family = gaussian()) {
  formula <- Y ~ ff(X1, xind = s_grid,
                    splinepars = list(bs = "ps", m = list(c(2, 1), c(2, 1)), k = c(8, 8))) +
                 Xc + zlin + s(zsmoo, k = 12)
  pffr(formula, yind = t_grid, data = dat, family = family,
       bs.yindex = list(bs = "ps", k = 12, m = c(2, 1)),
       bs.int = list(bs = "ps", k = 12, m = c(2, 1)))
}

fit_pffr_ar <- function(dat, s_grid, t_grid, rho, family = gaussian()) {
  formula <- Y ~ ff(X1, xind = s_grid,
                    splinepars = list(bs = "ps", m = list(c(2, 1), c(2, 1)), k = c(8, 8))) +
                 Xc + zlin + s(zsmoo, k = 12)
  pffr(formula, yind = t_grid, data = dat, family = family, algorithm = "bam", rho = rho,
       bs.yindex = list(bs = "ps", k = 12, m = c(2, 1)),
       bs.int = list(bs = "ps", k = 12, m = c(2, 1)))
}

fit_pffr_gls <- function(dat, s_grid, t_grid, hatSigma) {
  formula <- Y ~ ff(X1, xind = s_grid,
                    splinepars = list(bs = "ps", m = list(c(2, 1), c(2, 1)), k = c(8, 8))) +
                 Xc + zlin + s(zsmoo, k = 12)
  pffr_gls(formula, yind = t_grid, data = dat, hatSigma = hatSigma,
           bs.yindex = list(bs = "ps", k = 12, m = c(2, 1)),
           bs.int = list(bs = "ps", k = 12, m = c(2, 1)))
}

estimate_hatSigma <- function(fit, dat) {
  y_mat <- as.matrix(dat$Y)
  fit_mat <- tryCatch(fitted(fit, which = "mean"), error = function(e) fitted(fit))
  if (is.list(fit_mat) && !is.null(fit_mat$mean)) fit_mat <- fit_mat$mean
  resid_mat <- y_mat - fit_mat
  hatSigma <- cov(resid_mat)
  hatSigma <- 0.5 * (hatSigma + t(hatSigma))
  hatSigma + diag(1e-8 * mean(diag(hatSigma)), nrow(hatSigma))
}

estimate_rho <- function(fit, dat) {
  y_mat <- as.matrix(dat$Y)
  fit_mat <- tryCatch(fitted(fit, which = "mean"), error = function(e) fitted(fit))
  if (is.list(fit_mat) && !is.null(fit_mat$mean)) fit_mat <- fit_mat$mean
  resid_mat <- y_mat - fit_mat
  ny <- ncol(resid_mat)
  if (ny < 2) return(0)
  r1 <- as.vector(resid_mat[, 2:ny])
  r0 <- as.vector(resid_mat[, 1:(ny - 1)])
  rho <- suppressWarnings(cor(r0, r1))
  if (!is.finite(rho)) rho <- 0
  max(-0.99, min(0.99, rho))
}

# Scoring ----

score_fit <- function(fit, truth, alpha = 0.10, use_sandwich = FALSE) {
  # Extract coefficients on evaluation grid
  coefs <- tryCatch({
    tmp <- NULL
    invisible(capture.output(tmp <- coef(fit, sandwich = use_sandwich, n1 = 50, n2 = 25, n3 = 15)))
    tmp
  }, error = function(e) NULL)
  if (is.null(coefs)) return(NULL)

  z_crit <- qnorm(1 - alpha / 2)
  results <- list()

  # Get smooth terms from LP1 only (for gaulss/scat families)
  sm_list <- unname(coefs$smterms)
  sm_names <- names(fit$smooth)
  is_multi_lp <- !is.null(fit$family$nlp) && fit$family$nlp > 1
  if (is_multi_lp) {
    keep <- grepl("\\.1(\\(|$)", sm_names)
    sm_list <- sm_list[keep]
    sm_names <- sm_names[keep]
  }

  for (i in seq_along(sm_list)) {
    sm <- sm_list[[i]]
    if (is.null(sm)) next
    nm <- sm_names[[i]]

    # Identify term type
    term_type <- if (grepl("X1", nm) && sm$dim == 2) "ff"
                 else if (grepl("zsmoo", nm) && sm$dim == 2) "smooth"
                 else if (grepl("Xc", nm) && sm$dim == 1) "concurrent"
                 else if (grepl("zlin", nm) && sm$dim == 1) "linear"
                 else next

    est <- sm$coef[, "value"]
    se <- sm$coef[, "se"]

    # Extract grid coordinates from coef dataframe
    # For 1D: first col is x; For 2D: first col is x, second is y
    grid_x <- sm$coef[, 1]
    grid_y <- if (sm$dim >= 2) sm$coef[, 2] else NULL

    # Evaluate truth at evaluation grid points
    # For 2D terms: sm$x and sm$y are unique values; sm$coef has expand.grid(x,y) ordering
    truth_vec <- switch(term_type,
      ff = {
        beta_mat <- truth$truth_funs$beta_ff(sm$x, sm$y)
        beta_mat_centered <- center_beta_ff(beta_mat, sm$x, truth$ff_weights)
        as.vector(beta_mat_centered)  # column-major matches expand.grid
      },
      concurrent = truth$truth_funs$beta_concurrent(grid_x),
      linear = truth$truth_funs$beta_linear(grid_x),
      smooth = {
        raw_mat <- truth$truth_funs$f_smooth_scalar(sm$x, sm$y)
        # Center per t using observed z values
        m_t <- colMeans(truth$truth_funs$f_smooth_scalar(truth$zsmoo, sm$y))
        as.vector(sweep(raw_mat, 2, m_t, "-"))
      }
    )

    # Compute pointwise metrics then aggregate
    lower <- est - z_crit * se
    upper <- est + z_crit * se
    covered <- (truth_vec >= lower) & (truth_vec <= upper)
    width <- 2 * z_crit * se
    sq_err <- (est - truth_vec)^2

    results[[term_type]] <- data.frame(
      term_type = term_type,
      coverage = mean(covered),
      width = mean(width),
      rmse = sqrt(mean(sq_err)),
      n_grid = length(est)
    )
  }

  do.call(rbind, results)
}

# Main benchmark runner ----

make_dgp_grid <- function(tiny = FALSE) {
  if (tiny) {
    # Minimal grid that still tests all method combinations:
    # - iid+none: basic pffr methods
    # - ar1+none: tests AR residual estimation
    # - iid+bump: tests gaulss (heteroskedasticity)
    # - t6 variant: tests scat family
    return(rbind(
      expand.grid(n = 40L, nxgrid = 25L, nygrid = 30L, snr = 10, wiggliness = 2,
                  error_dist = "gaussian", corr_type = c("iid", "ar1"),
                  hetero_type = "none", stringsAsFactors = FALSE),
      expand.grid(n = 40L, nxgrid = 25L, nygrid = 30L, snr = 10, wiggliness = 2,
                  error_dist = "gaussian", corr_type = "iid",
                  hetero_type = "bump", stringsAsFactors = FALSE),
      expand.grid(n = 40L, nxgrid = 25L, nygrid = 30L, snr = 10, wiggliness = 2,
                  error_dist = "t6", corr_type = "iid",
                  hetero_type = "none", stringsAsFactors = FALSE)
    ))
  }
  expand.grid(
    n = 80L, nxgrid = 35L, nygrid = 45L,
    snr = c(2, 8, 25),
    wiggliness = c(2, 20),
    error_dist = c("gaussian", "t6"),
    corr_type = c("iid", "ar1", "gauss", "fourier"),
    hetero_type = c("none", "bump"),
    stringsAsFactors = FALSE
  )
}

fit_all_methods <- function(dgp, dat, s_grid, t_grid, cov_mat) {

  # Fit all relevant methods for this DGP, reusing base fit where possible
  fits <- list()

  # Base gaussian pffr (needed for pffr, sandwich, gls_est, ar)
  fits$pffr_gaussian <- tryCatch(
    fit_pffr(dat, s_grid, t_grid, gaussian()),
    error = function(e) NULL
  )
  if (is.null(fits$pffr_gaussian)) return(NULL)

  # pffr_gls with oracle covariance
  fits$pffr_gls <- tryCatch(
    fit_pffr_gls(dat, s_grid, t_grid, cov_mat),
    error = function(e) NULL
  )

  # pffr_gls with estimated covariance (from base fit residuals)
  hatSigma_est <- estimate_hatSigma(fits$pffr_gaussian, dat)
  fits$pffr_gls_est <- tryCatch(
    fit_pffr_gls(dat, s_grid, t_grid, hatSigma_est),
    error = function(e) NULL
  )

  # pffr_ar with estimated rho (from base fit residuals)
  rho_est <- estimate_rho(fits$pffr_gaussian, dat)
  fits$pffr_ar <- tryCatch(
    fit_pffr_ar(dat, s_grid, t_grid, rho_est, gaussian()),
    error = function(e) NULL
  )

  # gaulss (if heteroskedastic DGP)
  if (dgp$hetero_type != "none") {
    fits$pffr_gaulss <- tryCatch(
      fit_pffr(dat, s_grid, t_grid, mgcv::gaulss()),
      error = function(e) NULL
    )
  }

  # scat (if t-distributed errors)
  if (dgp$error_dist == "t6") {
    fits$pffr_scat <- tryCatch(
      fit_pffr(dat, s_grid, t_grid, mgcv::scat()),
      error = function(e) NULL
    )
  }

  fits
}

score_all_methods <- function(fits, truth, dgp, alpha = 0.10) {
  # Score all fitted models and return aggregated results
  results <- list()

  add_result <- function(fit, method, family, sandwich = FALSE) {
    if (is.null(fit)) return()
    scored <- score_fit(fit, truth, alpha, use_sandwich = sandwich)
    if (!is.null(scored)) {
      scored$method <- method
      scored$fit_family <- family
      results[[paste0(method, "_", family, if (sandwich) "_sw" else "")]] <<- scored
    }
  }

  # Gaussian family methods
  add_result(fits$pffr_gaussian, "pffr", "gaussian", FALSE)
  add_result(fits$pffr_gaussian, "pffr_sandwich", "gaussian", TRUE)
  add_result(fits$pffr_gls, "pffr_gls", "gaussian", FALSE)
  add_result(fits$pffr_gls_est, "pffr_gls_est", "gaussian", FALSE)
  add_result(fits$pffr_ar, "pffr_ar", "gaussian", FALSE)

  # gaulss methods (if applicable)
  if (dgp$hetero_type != "none") {
    add_result(fits$pffr_gaulss, "pffr", "gaulss", FALSE)
    add_result(fits$pffr_gaulss, "pffr_sandwich", "gaulss", TRUE)
  }

  # scat methods (if applicable)
  if (dgp$error_dist == "t6") {
    add_result(fits$pffr_scat, "pffr", "scat", FALSE)
    add_result(fits$pffr_scat, "pffr_sandwich", "scat", TRUE)
  }

  if (length(results) == 0) return(NULL)
  do.call(rbind, results)
}

run_one_rep <- function(dgp, rep_id, seed, alpha = 0.10) {
  sim_seed <- seed + 1000L * dgp$dgp_id + rep_id

  # Generate dataset once
  sim <- simulate_dataset(dgp, seed = sim_seed)

  # Fit all methods for this DGP (reusing fits where possible)
  fits <- fit_all_methods(dgp, sim$data, sim$s_grid, sim$t_grid, sim$cov_mat)
  if (is.null(fits)) {
    warning("Base pffr fit failed for dgp ", dgp$dgp_id, " rep ", rep_id)
    return(NULL)
  }

  # Score all fits
  out <- score_all_methods(fits, sim$truth, dgp, alpha)
  if (is.null(out)) return(NULL)

  # Add metadata
  out$dgp_id <- dgp$dgp_id
  out$rep_id <- rep_id
  out$snr <- dgp$snr
  out$wiggliness <- dgp$wiggliness
  out$error_dist <- dgp$error_dist
  out$corr_type <- dgp$corr_type
  out$hetero_type <- dgp$hetero_type
  rownames(out) <- NULL
  out
}

run_benchmark <- function(dgp_grid, n_rep = 10L, alpha = 0.10, seed = 2024L,
                          parallel = FALSE, n_workers = NULL, verbose = TRUE) {
  require_pkgs()

  dgp_grid$dgp_id <- seq_len(nrow(dgp_grid))
  n_dgp <- nrow(dgp_grid)
  total_tasks <- n_dgp * n_rep

  if (verbose) {
    cat("Running benchmark:\n")
    cat("  DGP settings:", n_dgp, "\n")
    cat("  Reps per DGP:", n_rep, "\n")
    cat("  Total tasks: ", total_tasks, "\n\n")
  }

  all_results <- list()
  task_idx <- 0

  for (i in seq_len(n_dgp)) {
    dgp <- dgp_grid[i, , drop = FALSE]
    for (r in seq_len(n_rep)) {
      task_idx <- task_idx + 1
      if (verbose && task_idx %% 10 == 0) {
        cat(sprintf("[%d/%d] DGP %d rep %d: snr=%s corr=%s hetero=%s\n",
                    task_idx, total_tasks, i, r, dgp$snr, dgp$corr_type, dgp$hetero_type))
      }

      result <- tryCatch(
        run_one_rep(dgp, r, seed, alpha),
        error = function(e) {
          if (verbose) warning("Error in DGP ", i, " rep ", r, ": ", e$message)
          NULL
        }
      )
      if (!is.null(result)) {
        all_results[[length(all_results) + 1]] <- result
      }
    }
  }

  if (length(all_results) == 0) {
    stop("No results collected")
  }

  do.call(rbind, all_results)
}

summarize_results <- function(scored) {
  # Aggregate coverage/width/rmse by method, family, term_type, and DGP settings
  agg_vars <- c("method", "fit_family", "term_type", "snr", "wiggliness",
                "error_dist", "corr_type", "hetero_type")
  scored_split <- split(scored, interaction(scored[, agg_vars], drop = TRUE))

  summary_list <- lapply(names(scored_split), function(nm) {
    df <- scored_split[[nm]]
    data.frame(
      method = df$method[1],
      fit_family = df$fit_family[1],
      term_type = df$term_type[1],
      snr = df$snr[1],
      wiggliness = df$wiggliness[1],
      error_dist = df$error_dist[1],
      corr_type = df$corr_type[1],
      hetero_type = df$hetero_type[1],
      coverage = mean(df$coverage),
      width = mean(df$width),
      rmse = mean(df$rmse),
      n_rep = nrow(df)
    )
  })

  do.call(rbind, summary_list)
}

# Main ----

if (sys.nframe() == 0 || identical(environment(), globalenv())) {
  require_pkgs()

  # Parse command line args or use tiny grid for testing
  args <- commandArgs(trailingOnly = TRUE)
  grid_size <- if (length(args) > 0) args[1] else "tiny"
  n_rep <- if (length(args) > 1) as.integer(args[2]) else 3L

  dgp_grid <- make_dgp_grid(tiny = (grid_size == "tiny"))

  cat("\n")
  cat("==============================================================\n")
  cat("  CI COVERAGE BENCHMARK\n")
  cat("==============================================================\n")
  cat("  Grid size:   ", grid_size, "\n")
  cat("  DGP settings:", nrow(dgp_grid), "\n")
  cat("  Reps per DGP:", n_rep, "\n")
  cat("==============================================================\n\n")

  start_time <- Sys.time()
  scored <- run_benchmark(dgp_grid, n_rep = n_rep, alpha = 0.10, seed = 2024L, verbose = TRUE)
  elapsed <- difftime(Sys.time(), start_time, units = "mins")

  cat("\n==============================================================\n")
  cat("  BENCHMARK COMPLETE\n")
  cat("  Elapsed:", round(as.numeric(elapsed), 1), "minutes\n")
  cat("  Total rows:", nrow(scored), "\n")
  cat("==============================================================\n\n")

  # Summarize
  summary_df <- summarize_results(scored)
  cat("Coverage summary by term type and method:\n")
  print(aggregate(coverage ~ term_type + method, data = summary_df, FUN = mean))
}
