# confint-benchmark.R
#
# Benchmark pointwise CI coverage/width from coef.pffr() under different
# residual covariance structures and error distributions.
#
# Methods: pffr, pffr_sandwich, pffr_gls_est, pffr_ar
# Families: gaussian, gaulss, scat
#
# Key design:
# - Loop over (dgp × rep), generate data ONCE per combo
# - Fit base pffr once, reuse for sandwich/estimation of rho/hatSigma
# - Return term-averaged metrics (not pointwise)
# - Parallel over (dgp × rep), incremental saving

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

require_pkgs <- function(pkg_dir = ".") {
  load_refund_dev(pkg_dir)
  stopifnot(
    requireNamespace("mgcv", quietly = TRUE),
    requireNamespace("mvtnorm", quietly = TRUE),
    requireNamespace("parallel", quietly = TRUE)
  )
  invisible(TRUE)
}

# Truth specification ----

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
  splines::bs(x, df = spec$df, knots = spec$knots,
              Boundary.knots = spec$boundary_knots,
              degree = spec$degree, intercept = spec$intercept)
}

make_truth_functions <- function(s_grid, t_grid, wiggliness, seed = NULL,
                                  k_te = 5L, k_uni = 7L) {
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

  coef_ff <- rw_coef_matrix(k_te, wiggliness)
  coef_lin <- rw_coef(k_uni, wiggliness)
  coef_mu <- rw_coef(k_uni, wiggliness)

  beta_ff <- function(s, t) {
    Bs <- eval_bs(s, s_spec)
    Bt <- eval_bs(t, t_te_spec)
    drop(Bs %*% coef_ff %*% t(Bt))
  }
  beta_linear <- function(t) drop(eval_bs(t, t_uni_spec) %*% coef_lin)
  mu_t <- function(t) drop(eval_bs(t, t_uni_spec) %*% coef_mu)

  list(beta_ff = beta_ff, beta_linear = beta_linear, mu_t = mu_t,
       spec = list(s = s_spec, t_te = t_te_spec, t_uni = t_uni_spec))
}

# Effect computation ----

ff_weights <- function(xind, method = "simpson") {
  nx <- length(xind)
  if (nx < 2) return(1)
  if (method == "simpson") {
    return(((xind[nx] - xind[1]) / (nx - 1)) / 3 *
             c(1, rep(c(4, 2), length = nx - 2), 1))
  }
  diffs <- diff(xind)
  c(diffs[1] / 2, (diffs[-1] + diffs[-length(diffs)]) / 2,
    diffs[length(diffs)] / 2)
}

center_beta_ff <- function(beta_st, s_grid, weights = NULL) {
  w <- weights %||% ff_weights(s_grid)
  sweep(beta_st, 2, colSums(beta_st * w) / sum(w), "-")
}

compute_ff_effect <- function(X, beta_st, xind) {
  w <- ff_weights(xind)
  L <- matrix(rep(w, each = nrow(X)), nrow = nrow(X), byrow = FALSE)
  (L * X) %*% beta_st
}

# Error structure ----

make_error_structure <- function(t_grid, corr_type = "iid", corr_param = NULL,
                                  hetero_type = "none", hetero_param = NULL) {
  ny <- length(t_grid)

  make_pd <- function(mat) {
    diag_mean <- mean(diag(mat))
    if (!is.finite(diag_mean) || diag_mean <= 0) diag_mean <- 1
    jitter <- 1e-10 * diag_mean
    for (k in 1:8) {
      ok <- tryCatch({ chol(mat + diag(jitter, nrow(mat))); TRUE },
                     error = function(e) FALSE)
      if (ok) return(mat + diag(jitter, nrow(mat)))
      jitter <- jitter * 10
    }
    mat + diag(jitter, nrow(mat))
  }

  corr_mat <- switch(corr_type,
    iid = diag(1, ny),
    ar1 = {
      rho <- corr_param %||% 0.4
      rho^abs(outer(seq_len(ny), seq_len(ny), "-"))
    },
    gauss = {
      phi <- corr_param %||% 0.15
      dist_mat <- abs(outer(t_grid, t_grid, "-"))
      exp(-(dist_mat^2) / (2 * phi^2))
    },
    fourier = {
      period1 <- corr_param %||% 0.33
      B <- cbind(sin(2 * pi * t_grid / period1), cos(2 * pi * t_grid / period1),
                 sin(2 * pi * t_grid / 0.17), cos(2 * pi * t_grid / 0.17))
      B <- scale(B, center = TRUE, scale = FALSE)
      stats::cov2cor(B %*% t(B) + diag(0.25, ny))
    }
  )
  corr_mat <- make_pd(corr_mat)

  sigma_t <- switch(hetero_type,
    none = rep(1, ny),
    linear = 1 + (hetero_param %||% 0.6) * (t_grid - mean(t_grid)),
    u = {
      shape <- (t_grid - 0.5)^2
      1 + (hetero_param %||% 0.9) * shape / max(shape)
    },
    bump = 1 + (hetero_param %||% 1.0) * exp(-0.5 * ((t_grid - 0.7) / 0.10)^2)
  )
  sigma_t <- pmax(sigma_t, 1e-6)

  cov_mat <- diag(sigma_t) %*% corr_mat %*% diag(sigma_t)
  cov_mat <- make_pd(cov_mat)

  list(corr_mat = corr_mat, cov_mat = cov_mat, sigma_t = sigma_t,
       rho = if (corr_type == "ar1") (corr_param %||% 0.4) else NA_real_)
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

# Dataset simulation ----

simulate_dataset <- function(n, nxgrid, nygrid, snr, wiggliness,
                              error_dist = "gaussian", corr_type = "iid",
                              corr_param = NULL, hetero_type = "none",
                              hetero_param = NULL, seed = NULL,
                              terms = c("ff", "linear")) {
  if (!is.null(seed)) set.seed(seed)

  s_grid <- seq(0, 1, length.out = nxgrid)
  t_grid <- seq(0, 1, length.out = nygrid)

  truth <- make_truth_functions(s_grid, t_grid, wiggliness,
                                 seed = if (!is.null(seed)) seed + 1e5L else NULL)

  # Covariates
  X1 <- if ("ff" %in% terms) {
    X <- splines::bs(s_grid, df = 15, intercept = TRUE)
    coef_mat <- matrix(rnorm(n * 15), nrow = 15)
    mat <- t(X %*% coef_mat)
    sweep(mat, 2, colMeans(mat), "-")  # center per s
  } else NULL

  zlin <- if ("linear" %in% terms) {
    z <- rnorm(n)
    z - mean(z)
  } else NULL

  # True effects
  mu_t <- truth$mu_t(t_grid)
  eta_int <- matrix(rep(mu_t, each = n), nrow = n)

  w_s <- if ("ff" %in% terms) ff_weights(s_grid) else NULL
  beta_st <- if ("ff" %in% terms) center_beta_ff(truth$beta_ff(s_grid, t_grid), s_grid, w_s) else NULL
  eta_ff <- if ("ff" %in% terms) compute_ff_effect(X1, beta_st, s_grid) else NULL

  beta_lin <- if ("linear" %in% terms) truth$beta_linear(t_grid) else NULL
  eta_lin <- if ("linear" %in% terms) outer(zlin, beta_lin) else NULL

  eta_total <- eta_int
  if (!is.null(eta_ff)) eta_total <- eta_total + eta_ff
  if (!is.null(eta_lin)) eta_total <- eta_total + eta_lin

  # Errors
  err <- make_error_structure(t_grid, corr_type, corr_param, hetero_type, hetero_param)
  eps <- sample_errors(n, err$cov_mat, error_dist)

  # Scale to target SNR
  var_eta <- var(as.vector(eta_total))
  var_eps <- var(as.vector(eps))
  scale_fac <- sqrt(var_eta / (snr * var_eps))
  eps <- eps * scale_fac
  err$cov_mat <- err$cov_mat * scale_fac^2

  Y <- eta_total + eps

  dat <- data.frame(id = seq_len(n))
  dat$Y <- I(Y)
  if (!is.null(X1)) dat$X1 <- I(X1)
  if (!is.null(zlin)) dat$zlin <- zlin

  truth_out <- list(
    eta = eta_total,
    beta = list(intercept = mu_t, X1 = beta_st, zlin = beta_lin),
    truth_funs = list(beta_ff = truth$beta_ff, beta_linear = truth$beta_linear, mu_t = truth$mu_t),
    ff_weights = w_s,
    terms = terms
  )

  list(data = dat, s_grid = s_grid, t_grid = t_grid, err = err, truth = truth_out)
}

# Fitting helpers ----

fit_family_object <- function(fit_family) {
  switch(fit_family,
    gaussian = stats::gaussian(),
    gaulss = mgcv::gaulss(),
    scat = mgcv::scat(),
    stop("Unknown family: ", fit_family)
  )
}

estimate_hatSigma <- function(fit, dat) {
  y_mat <- as.matrix(dat$Y)
  fit_mat <- tryCatch(fitted(fit, which = "mean"), error = function(e) fitted(fit))
  if (is.list(fit_mat) && !is.null(fit_mat$mean)) fit_mat <- fit_mat$mean
  resid_mat <- y_mat - as.matrix(fit_mat)
  hatSigma <- cov(resid_mat)
  hatSigma <- 0.5 * (hatSigma + t(hatSigma))
  hatSigma + diag(1e-8 * mean(diag(hatSigma)), nrow(hatSigma))
}

estimate_rho <- function(fit, dat, clamp = 0.99) {
  y_mat <- as.matrix(dat$Y)
  fit_mat <- tryCatch(fitted(fit, which = "mean"), error = function(e) fitted(fit))
  if (is.list(fit_mat) && !is.null(fit_mat$mean)) fit_mat <- fit_mat$mean
  resid_mat <- y_mat - as.matrix(fit_mat)
  ny <- ncol(resid_mat)
  if (ny < 2) return(0)
  r1 <- as.vector(resid_mat[, 2:ny])
  r0 <- as.vector(resid_mat[, 1:(ny - 1)])
  rho <- suppressWarnings(cor(r0, r1))
  if (!is.finite(rho)) rho <- 0
  max(-clamp, min(clamp, rho))
}

# Scoring ----

score_fit <- function(fit, truth, s_grid, t_grid, use_sandwich = FALSE,
                      alpha = 0.10, terms = c("ff", "linear")) {
  # Extract coefficients with SEs
  coefs <- suppressMessages(coef(fit, sandwich = use_sandwich,
                                  n1 = length(t_grid), n2 = length(s_grid), n3 = 20))

  z_crit <- qnorm(1 - alpha / 2)
  results <- list()

  # Score each term
  for (term in terms) {
    sm_name <- if (term == "ff") "X1" else if (term == "linear") "zlin" else term
    sm_idx <- grep(sm_name, names(coefs$smterms), fixed = TRUE)
    if (length(sm_idx) == 0) next

    sm <- coefs$smterms[[sm_idx[1]]]
    est <- sm$coef[, "value"]
    se <- sm$coef[, "se"]

    # Get truth on same grid
    if (term == "ff") {
      truth_raw <- truth$truth_funs$beta_ff(sm$coef[, 1], sm$coef[, 2])
      truth_vec <- as.vector(center_beta_ff(truth_raw, sm$coef[, 1], truth$ff_weights))
    } else if (term == "linear") {
      truth_vec <- truth$truth_funs$beta_linear(sm$coef[, 1])
    } else next

    lower <- est - z_crit * se
    upper <- est + z_crit * se
    covered <- (truth_vec >= lower) & (truth_vec <= upper)

    results[[term]] <- data.frame(
      term_type = term,
      coverage = mean(covered),
      width = mean(2 * z_crit * se),
      rmse = sqrt(mean((est - truth_vec)^2)),
      n_grid = length(est)
    )
  }

  # Score E(y)
  pred <- predict(fit, se.fit = TRUE)
  if (is.list(pred) && !is.null(pred$fit)) {
    fit_vals <- as.vector(pred$fit)
    se_vals <- as.vector(pred$se.fit)
  } else {
    fit_vals <- as.vector(pred)
    se_vals <- rep(NA_real_, length(fit_vals))
  }

  truth_eta <- as.vector(truth$eta)
  if (length(fit_vals) == length(truth_eta) && all(!is.na(se_vals))) {
    lower_ey <- fit_vals - z_crit * se_vals
    upper_ey <- fit_vals + z_crit * se_vals
    covered_ey <- (truth_eta >= lower_ey) & (truth_eta <= upper_ey)

    results$Ey <- data.frame(
      term_type = "Ey",
      coverage = mean(covered_ey),
      width = mean(2 * z_crit * se_vals),
      rmse = sqrt(mean((fit_vals - truth_eta)^2)),
      n_grid = length(fit_vals)
    )
  }

  do.call(rbind, results)
}

# Fit all methods for one dataset ----

fit_all_methods <- function(dat, formula, t_grid, s_grid, truth, methods,
                            family, alpha = 0.10, fit_args = list()) {
  results <- list()
  terms <- truth$terms

  # Always fit base pffr first (needed for most methods)
  t0 <- Sys.time()
  fit_base <- tryCatch(
    do.call(pffr, c(list(formula = formula, yind = t_grid, data = dat,
                         family = fit_family_object(family)), fit_args)),
    error = function(e) { message("pffr failed: ", e$message); NULL }
  )
  time_base <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  if (is.null(fit_base)) return(NULL)

  # Estimate rho and hatSigma from base fit
  rho_est <- estimate_rho(fit_base, dat)
  hatSigma_est <- estimate_hatSigma(fit_base, dat)

  for (method in methods) {
    t0 <- Sys.time()

    fit <- NULL
    use_sandwich <- FALSE

    if (method == "pffr") {
      fit <- fit_base
      time_fit <- time_base
    } else if (method == "pffr_sandwich") {
      fit <- fit_base
      use_sandwich <- TRUE
      time_fit <- time_base
    } else if (method == "pffr_gls_est") {
      if (family != "gaussian") next  # GLS only for gaussian
      fit <- tryCatch(
        do.call(pffr_gls, c(list(formula = formula, yind = t_grid, data = dat,
                                  hatSigma = hatSigma_est), fit_args)),
        error = function(e) { message("pffr_gls_est failed: ", e$message); NULL }
      )
      time_fit <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    } else if (method == "pffr_ar") {
      if (family != "gaussian") next  # AR only for gaussian with bam
      fit <- tryCatch(
        do.call(pffr, c(list(formula = formula, yind = t_grid, data = dat,
                             family = fit_family_object(family),
                             algorithm = "bam", rho = rho_est), fit_args)),
        error = function(e) { message("pffr_ar failed: ", e$message); NULL }
      )
      time_fit <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    } else next

    if (is.null(fit)) next

    t0_score <- Sys.time()
    scored <- tryCatch(
      score_fit(fit, truth, s_grid, t_grid, use_sandwich, alpha, terms),
      error = function(e) { message("score_fit failed: ", e$message); NULL }
    )
    time_score <- as.numeric(difftime(Sys.time(), t0_score, units = "secs"))

    if (!is.null(scored)) {
      scored$method <- method
      scored$family <- family
      scored$time_fit_sec <- time_fit
      scored$time_score_sec <- time_score
      results[[paste(method, family, sep = "_")]] <- scored
    }
  }

  if (length(results) == 0) return(NULL)
  do.call(rbind, results)
}

# Main benchmark ----

make_dgp_grid <- function() {
  expand.grid(
    n = 80L, nxgrid = 35L, nygrid = 45L,
    snr = c(2, 8, 25), wiggliness = c(2, 20),
    error_dist = c("gaussian", "t6"),
    corr_type = c("iid", "ar1", "gauss", "fourier"),
    hetero_type = c("none", "bump"),
    stringsAsFactors = FALSE
  )
}

make_test_grid <- function() {
  # Smart edge-case coverage, not full cartesian
  rbind(
    data.frame(n=30L, nxgrid=15L, nygrid=20L, snr=10, wiggliness=2,
               error_dist="gaussian", corr_type="iid", hetero_type="none"),
    data.frame(n=30L, nxgrid=15L, nygrid=20L, snr=10, wiggliness=2,
               error_dist="gaussian", corr_type="ar1", hetero_type="none"),
    data.frame(n=30L, nxgrid=15L, nygrid=20L, snr=10, wiggliness=2,
               error_dist="gaussian", corr_type="iid", hetero_type="bump"),
    data.frame(n=30L, nxgrid=15L, nygrid=20L, snr=10, wiggliness=2,
               error_dist="t6", corr_type="iid", hetero_type="none"),
    data.frame(n=30L, nxgrid=15L, nygrid=20L, snr=10, wiggliness=2,
               error_dist="t6", corr_type="ar1", hetero_type="bump")
  )
}

get_families_for_dgp <- function(dgp) {
  fams <- "gaussian"
  if (dgp$hetero_type != "none") fams <- c(fams, "gaulss", "scat")
  else if (dgp$error_dist == "t6") fams <- c(fams, "scat")
  fams
}

run_one_dgp_rep <- function(dgp, rep_id, seed, methods, alpha, terms, fit_args) {
  sim <- simulate_dataset(
    n = dgp$n, nxgrid = dgp$nxgrid, nygrid = dgp$nygrid,
    snr = dgp$snr, wiggliness = dgp$wiggliness,
    error_dist = dgp$error_dist, corr_type = dgp$corr_type,
    hetero_type = dgp$hetero_type, seed = seed, terms = terms
  )

  # Build formula
  formula_terms <- c(
    if ("ff" %in% terms) "ff(X1, xind = s_grid)" else NULL,
    if ("linear" %in% terms) "zlin" else NULL
  )
  formula <- as.formula(paste("Y ~", paste(formula_terms, collapse = " + ")))

  families <- get_families_for_dgp(dgp)
  results <- list()

  for (fam in families) {
    res <- fit_all_methods(
      dat = sim$data, formula = formula,
      t_grid = sim$t_grid, s_grid = sim$s_grid,
      truth = sim$truth, methods = methods,
      family = fam, alpha = alpha, fit_args = fit_args
    )
    if (!is.null(res)) results[[fam]] <- res
  }

  if (length(results) == 0) return(NULL)

  out <- do.call(rbind, results)
  out$dgp_id <- dgp$dgp_id
  out$rep_id <- rep_id
  out$snr <- dgp$snr
  out$wiggliness <- dgp$wiggliness
  out$error_dist <- dgp$error_dist
  out$corr_type <- dgp$corr_type
  out$hetero_type <- dgp$hetero_type
  out
}

run_benchmark <- function(dgp_grid, n_rep = 10L,
                          methods = c("pffr", "pffr_sandwich", "pffr_gls_est", "pffr_ar"),
                          alpha = 0.10, seed = 2024L, n_workers = 1L,
                          output_dir = NULL, terms = c("ff", "linear"),
                          fit_args = list()) {
  require_pkgs()

  dgp_grid$dgp_id <- seq_len(nrow(dgp_grid))

  # Build task list: all (dgp_id, rep_id) combinations
  tasks <- expand.grid(dgp_idx = seq_len(nrow(dgp_grid)), rep_id = seq_len(n_rep))
  tasks$task_id <- seq_len(nrow(tasks))
  tasks$seed <- seed + 1000L * tasks$dgp_idx + tasks$rep_id

  n_tasks <- nrow(tasks)
  cat(sprintf("\n=== BENCHMARK: %d DGP settings x %d reps = %d tasks ===\n",
              nrow(dgp_grid), n_rep, n_tasks))
  cat(sprintf("Methods: %s\n", paste(methods, collapse = ", ")))
  cat(sprintf("Terms: %s\n", paste(terms, collapse = ", ")))
  cat(sprintf("Workers: %d\n\n", n_workers))

  if (!is.null(output_dir) && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  run_task <- function(task_row) {
    dgp <- dgp_grid[task_row$dgp_idx, , drop = FALSE]
    run_one_dgp_rep(dgp, task_row$rep_id, task_row$seed, methods, alpha, terms, fit_args)
  }

  start_time <- Sys.time()
  all_results <- list()

  if (n_workers > 1 && .Platform$OS.type != "windows") {
    # Parallel execution with mclapply
    chunk_size <- max(1, n_workers * 2)
    task_chunks <- split(seq_len(n_tasks), ceiling(seq_len(n_tasks) / chunk_size))

    for (chunk_idx in seq_along(task_chunks)) {
      chunk_tasks <- task_chunks[[chunk_idx]]
      chunk_results <- parallel::mclapply(chunk_tasks, function(i) {
        run_task(tasks[i, , drop = FALSE])
      }, mc.cores = n_workers)

      for (j in seq_along(chunk_tasks)) {
        task_i <- chunk_tasks[j]
        res <- chunk_results[[j]]
        if (!is.null(res)) {
          all_results[[task_i]] <- res
          # Log summary
          dgp <- dgp_grid[tasks$dgp_idx[task_i], ]
          summary_str <- paste(sapply(split(res, res$method), function(m) {
            sprintf("%s:%.2f", m$method[1], mean(m$coverage))
          }), collapse = " | ")
          cat(sprintf("[%d/%d] dgp=%d rep=%d %s/%s/%s | %s\n",
                      task_i, n_tasks, dgp$dgp_id, tasks$rep_id[task_i],
                      dgp$corr_type, dgp$hetero_type, dgp$error_dist, summary_str))
        }
      }

      # Incremental save after each chunk
      if (!is.null(output_dir)) {
        combined <- do.call(rbind, all_results[!sapply(all_results, is.null)])
        if (!is.null(combined) && nrow(combined) > 0) {
          saveRDS(combined, file.path(output_dir, "results_partial.rds"))
        }
      }
    }
  } else {
    # Sequential execution
    for (i in seq_len(n_tasks)) {
      res <- tryCatch(run_task(tasks[i, , drop = FALSE]),
                      error = function(e) { message("Task ", i, " failed: ", e$message); NULL })
      if (!is.null(res)) {
        all_results[[i]] <- res
        dgp <- dgp_grid[tasks$dgp_idx[i], ]
        summary_str <- paste(sapply(split(res, res$method), function(m) {
          sprintf("%s:%.2f", m$method[1], mean(m$coverage))
        }), collapse = " | ")
        cat(sprintf("[%d/%d] dgp=%d rep=%d %s/%s/%s | %s\n",
                    i, n_tasks, dgp$dgp_id, tasks$rep_id[i],
                    dgp$corr_type, dgp$hetero_type, dgp$error_dist, summary_str))

        # Incremental save
        if (!is.null(output_dir) && i %% 5 == 0) {
          combined <- do.call(rbind, all_results[!sapply(all_results, is.null)])
          if (!is.null(combined)) saveRDS(combined, file.path(output_dir, "results_partial.rds"))
        }
      }
    }
  }

  elapsed <- difftime(Sys.time(), start_time, units = "mins")
  cat(sprintf("\n=== COMPLETE: %.1f minutes ===\n", as.numeric(elapsed)))

  combined <- do.call(rbind, all_results[!sapply(all_results, is.null)])
  rownames(combined) <- NULL

  if (!is.null(output_dir)) {
    saveRDS(combined, file.path(output_dir, "results_final.rds"))
    write.csv(combined, file.path(output_dir, "results_final.csv"), row.names = FALSE)
    cat(sprintf("Saved to: %s\n", output_dir))
  }

  combined
}

# Summary ----

summarize_results <- function(results) {
  agg_vars <- c("method", "family", "term_type", "snr", "wiggliness",
                "error_dist", "corr_type", "hetero_type")

  agg <- aggregate(
    cbind(coverage, width, rmse, time_fit_sec) ~ method + family + term_type +
      snr + wiggliness + error_dist + corr_type + hetero_type,
    data = results, FUN = mean
  )
  agg$n_rep <- as.vector(table(interaction(results[, agg_vars], drop = TRUE)))
  agg
}

# Run if sourced directly ----

if (sys.nframe() == 0) {
  require_pkgs()

  output_dir <- file.path(getwd(), "ci-benchmark")
  dgp_grid <- make_test_grid()

  results <- run_benchmark(
    dgp_grid, n_rep = 1L,
    methods = c("pffr", "pffr_sandwich", "pffr_gls_est", "pffr_ar"),
    alpha = 0.10, seed = 2024L, n_workers = 1L,
    output_dir = output_dir, terms = c("ff", "linear")
  )

  if (!is.null(results) && nrow(results) > 0) {
    cat("\n=== SUMMARY ===\n")
    summary_df <- summarize_results(results)
    print(summary_df[, c("method", "family", "term_type", "corr_type",
                         "hetero_type", "coverage", "width", "rmse")])
  }
}
