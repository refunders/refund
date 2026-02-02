# ===========================================================================
# Sandwich Benchmark Part 2: Heteroskedasticity, Complex Correlation, ff Terms
# ===========================================================================
#
# Three scenarios:
#   A) IID errors + strong heteroskedasticity (x5, x10, x20 variance ratio),
#      coverage split by high-variance vs low-variance regions
#   B) Homoskedastic errors + complex non-monotonic autocorrelation (fourier),
#      compared to AR(1) of similar strength
#   C) Coverage for ff(X1) term under correlated errors
#
# Usage:
#   Rscript ci-benchmark/sandwich-benchmark-2.R
# ===========================================================================

devtools::load_all(".")
library(mgcv)
library(dplyr)
library(tidyr)
library(ggplot2)

source("ci-benchmark/benchmark-utils.R")

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

#' Fit pffr and extract coefficient CIs for all 4 sandwich methods
#'
#' @param m A fitted pffr model (no sandwich).
#' @param term_name Grep pattern for coef.pffr smterms name.
#' @param truth_vec Numeric vector of true values on the eval grid.
#' @param n2 Grid size for 2D terms.
#' @returns Data frame with method, t (or s,t index), est, se, covered, width.
extract_all_methods <- function(m, term_name, truth_vec, n2 = 25) {
  methods <- list(
    nosandwich = list(sandwich = FALSE, freq = FALSE),
    hc = list(sandwich = "hc", freq = FALSE),
    cluster_freq = list(sandwich = TRUE, freq = TRUE),
    cluster_bayes = list(sandwich = TRUE, freq = FALSE)
  )

  results <- list()
  for (mname in names(methods)) {
    args <- methods[[mname]]
    cc <- tryCatch(
      coef(
        m,
        sandwich = args$sandwich,
        freq = args$freq,
        seWithMean = FALSE,
        n1 = 50,
        n2 = n2
      ),
      error = function(e) NULL
    )
    if (is.null(cc)) next

    sm_names <- names(cc$smterms)
    idx <- grep(term_name, sm_names, ignore.case = TRUE)
    if (length(idx) == 0) next
    ti <- cc$smterms[[idx[1]]]
    if (is.null(ti) || is.null(ti$coef)) next

    est <- ti$coef$value
    se <- ti$coef$se
    if (is.null(se) || length(est) != length(truth_vec)) next

    results[[mname]] <- data.frame(
      method = mname,
      idx = seq_along(est),
      est = est,
      se = se,
      truth = truth_vec,
      lower = est - 1.96 * se,
      upper = est + 1.96 * se,
      stringsAsFactors = FALSE
    ) |>
      mutate(
        covered = truth >= lower & truth <= upper,
        width = upper - lower
      )
  }
  bind_rows(results)
}

# =========================================================================
# SCENARIO A: Strong heteroskedasticity, IID errors
# =========================================================================

cat("=== SCENARIO A: Strong Heteroskedasticity ===\n\n")

#' Simulate FOS data with strong heteroskedasticity
#'
#' @param n Number of curves.
#' @param T_grid Number of time points.
#' @param var_ratio Ratio of max to min variance (e.g., 5, 10, 20).
#' @param hetero_shape "linear", "bump", "step" — shape of variance function.
#' @param snr Signal-to-noise ratio.
simulate_strong_hetero <- function(
  n,
  T_grid = 40,
  var_ratio = 10,
  hetero_shape = "linear",
  snr = 5
) {
  t_grid <- seq(0, 1, length.out = T_grid)

  beta0 <- sin(2 * pi * t_grid)
  beta1 <- cos(2 * pi * t_grid)

  x <- rnorm(n)
  x <- x - mean(x)

  signal <- outer(rep(1, n), beta0) + outer(x, beta1)

  # Build variance pattern: min SD = 1, max SD = sqrt(var_ratio)
  # so variance ranges from 1 to var_ratio
  sd_ratio <- sqrt(var_ratio)
  sigma_t <- switch(
    hetero_shape,
    linear = {
      # Linearly increasing: left = 1, right = sd_ratio
      1 + (sd_ratio - 1) * t_grid
    },
    bump = {
      # Gaussian bump centered at t=0.7: baseline SD=1, peak SD=sd_ratio
      1 + (sd_ratio - 1) * exp(-0.5 * ((t_grid - 0.7) / 0.1)^2)
    },
    step = {
      # Step function: low variance for t<0.5, high for t>=0.5
      ifelse(t_grid < 0.5, 1, sd_ratio)
    }
  )

  # Scale to achieve target SNR
  signal_var <- var(as.vector(signal))
  base_var <- mean(sigma_t^2)
  error_scale <- sqrt(signal_var / (snr * base_var))
  sigma_t_scaled <- sigma_t * error_scale

  # IID errors (no autocorrelation) — each time point independent
  epsilon <- matrix(rnorm(n * T_grid), nrow = n)
  epsilon <- sweep(epsilon, 2, sigma_t_scaled, "*")

  Y <- signal + epsilon

  list(
    data = data.frame(Y = I(Y), x = x),
    yind = t_grid,
    truth = list(beta0 = beta0, beta1 = beta1),
    sigma_t = sigma_t_scaled
  )
}

run_scenario_a <- function(
  n_values = c(50, 100),
  var_ratios = c(5, 10, 20),
  hetero_shapes = c("linear", "bump"),
  B_rep = 200,
  seed = 123
) {
  set.seed(seed)

  settings <- expand.grid(
    n = n_values,
    var_ratio = var_ratios,
    hetero_shape = hetero_shapes,
    stringsAsFactors = FALSE
  )

  cat("Scenario A:", nrow(settings), "settings x", B_rep, "reps\n")

  all_results <- list()

  for (s in seq_len(nrow(settings))) {
    n <- settings$n[s]
    vr <- settings$var_ratio[s]
    hs <- settings$hetero_shape[s]

    cat(sprintf(
      "[%d/%d] n=%d, var_ratio=%d, shape=%s\n",
      s,
      nrow(settings),
      n,
      vr,
      hs
    ))

    for (b in seq_len(B_rep)) {
      sim <- simulate_strong_hetero(n = n, var_ratio = vr, hetero_shape = hs)

      m <- tryCatch(
        pffr(Y ~ x, yind = sim$yind, data = sim$data),
        error = function(e) NULL
      )
      if (is.null(m)) next

      # Extract beta1 coefficients
      c_base <- tryCatch(
        coef(m, sandwich = FALSE, seWithMean = FALSE, n1 = 50),
        error = function(e) NULL
      )
      if (is.null(c_base)) next

      sm_names <- names(c_base$smterms)
      x_idx <- grep("^x\\(", sm_names)
      if (length(x_idx) == 0) next
      ti <- c_base$smterms[[x_idx[1]]]
      t_vals <- ti$x

      beta1_true <- approx(
        seq(0, 1, length.out = length(sim$truth$beta1)),
        sim$truth$beta1,
        xout = t_vals
      )$y

      # Get all methods
      res <- extract_all_methods(m, "^x\\(", beta1_true)
      if (nrow(res) == 0) next

      # Tag with setting info
      res$n <- n
      res$var_ratio <- vr
      res$hetero_shape <- hs
      res$rep <- b
      res$t <- t_vals[res$idx]

      # Classify region: interpolate sigma_t to eval grid
      sigma_at_t <- approx(
        seq(0, 1, length.out = length(sim$sigma_t)),
        sim$sigma_t,
        xout = t_vals
      )$y
      median_sigma <- median(sigma_at_t)
      res$region <- ifelse(
        sigma_at_t[res$idx] > median_sigma,
        "high_var",
        "low_var"
      )

      all_results[[length(all_results) + 1]] <- res
    }
  }

  raw <- bind_rows(all_results)

  # Summary: overall and by region
  summary_overall <- raw |>
    filter(!is.na(covered)) |>
    group_by(n, var_ratio, hetero_shape, method) |>
    summarise(
      coverage = mean(covered),
      mean_width = mean(width),
      se_ratio = mean(se, na.rm = TRUE) / sd(est - truth, na.rm = TRUE),
      .groups = "drop"
    )

  summary_by_region <- raw |>
    filter(!is.na(covered)) |>
    group_by(n, var_ratio, hetero_shape, method, region) |>
    summarise(
      coverage = mean(covered),
      mean_width = mean(width),
      .groups = "drop"
    )

  list(
    raw = raw,
    summary_overall = summary_overall,
    summary_by_region = summary_by_region
  )
}

# =========================================================================
# SCENARIO B: Complex non-monotonic autocorrelation
# =========================================================================

cat("=== SCENARIO B: Complex Autocorrelation ===\n\n")

#' Simulate FOS data with arbitrary correlation structure
simulate_complex_corr <- function(
  n,
  T_grid = 40,
  corr_type = "fourier",
  corr_param = NULL,
  snr = 5
) {
  t_grid <- seq(0, 1, length.out = T_grid)

  beta0 <- sin(2 * pi * t_grid)
  beta1 <- cos(2 * pi * t_grid)

  x <- rnorm(n)
  x <- x - mean(x)

  signal <- outer(rep(1, n), beta0) + outer(x, beta1)

  # Use benchmark-utils error structure
  err_struct <- make_error_structure(
    t_grid,
    corr_type = corr_type,
    corr_param = corr_param,
    hetero_type = "none"
  )

  # Scale for SNR
  signal_var <- var(as.vector(signal))
  base_var <- mean(diag(err_struct$cov_mat))
  error_scale <- sqrt(signal_var / (snr * base_var))
  cov_scaled <- err_struct$cov_mat * error_scale^2

  # Generate correlated errors
  L <- t(chol(cov_scaled))
  epsilon <- matrix(rnorm(n * T_grid), nrow = n) %*% t(L)

  Y <- signal + epsilon

  list(
    data = data.frame(Y = I(Y), x = x),
    yind = t_grid,
    truth = list(beta0 = beta0, beta1 = beta1),
    corr_mat = err_struct$corr_mat,
    corr_type = corr_type
  )
}

run_scenario_b <- function(
  n_values = c(50, 100),
  B_rep = 200,
  seed = 456
) {
  set.seed(seed)

  # Compare: fourier (non-monotonic), gauss (smooth monotonic), ar1 (standard)
  # Use default corr_param for each — benchmark-utils sets them to comparable
  # RMS off-diagonal strength (~0.4)
  corr_types <- c("ar1", "gauss", "fourier")

  settings <- expand.grid(
    n = n_values,
    corr_type = corr_types,
    stringsAsFactors = FALSE
  )

  cat("Scenario B:", nrow(settings), "settings x", B_rep, "reps\n")

  all_results <- list()

  for (s in seq_len(nrow(settings))) {
    n <- settings$n[s]
    ct <- settings$corr_type[s]

    cat(sprintf("[%d/%d] n=%d, corr=%s\n", s, nrow(settings), n, ct))

    for (b in seq_len(B_rep)) {
      sim <- simulate_complex_corr(n = n, corr_type = ct)

      m <- tryCatch(
        pffr(Y ~ x, yind = sim$yind, data = sim$data),
        error = function(e) NULL
      )
      if (is.null(m)) next

      c_base <- tryCatch(
        coef(m, sandwich = FALSE, seWithMean = FALSE, n1 = 50),
        error = function(e) NULL
      )
      if (is.null(c_base)) next

      sm_names <- names(c_base$smterms)
      x_idx <- grep("^x\\(", sm_names)
      if (length(x_idx) == 0) next
      ti <- c_base$smterms[[x_idx[1]]]
      t_vals <- ti$x

      beta1_true <- approx(
        seq(0, 1, length.out = length(sim$truth$beta1)),
        sim$truth$beta1,
        xout = t_vals
      )$y

      res <- extract_all_methods(m, "^x\\(", beta1_true)
      if (nrow(res) == 0) next

      res$n <- n
      res$corr_type <- ct
      res$rep <- b
      res$t <- t_vals[res$idx]

      all_results[[length(all_results) + 1]] <- res
    }
  }

  raw <- bind_rows(all_results)

  summary_b <- raw |>
    filter(!is.na(covered)) |>
    group_by(n, corr_type, method) |>
    summarise(
      coverage = mean(covered),
      mean_width = mean(width),
      se_ratio = mean(se, na.rm = TRUE) / sd(est - truth, na.rm = TRUE),
      .groups = "drop"
    )

  list(raw = raw, summary = summary_b)
}

# =========================================================================
# SCENARIO C: ff term coverage under correlated errors
# =========================================================================

cat("=== SCENARIO C: ff Term Coverage ===\n\n")

#' Simulate data with ff(X1) term and correlated errors
simulate_ff_data <- function(
  n,
  nxgrid = 25,
  nygrid = 30,
  rho = 0.6,
  snr = 5,
  seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  # Use pffr_simulate for the ff term
  dat <- pffr_simulate(
    Y ~ ff(X1, xind = s),
    n = n,
    nxgrid = nxgrid,
    nygrid = nygrid,
    SNR = snr,
    effects = list(X1 = "cosine"),
    intercept = "zero",
    seed = seed
  )

  s_grid <- attr(dat, "xindex")
  t_grid <- attr(dat, "yindex")
  truth <- attr(dat, "truth")

  # Now add AR(1) correlated errors (replacing the iid errors from pffr_simulate)
  eta_mat <- truth$eta # true signal without error

  # Generate AR(1) errors
  ny <- length(t_grid)
  tt <- seq_len(ny)
  R <- rho^abs(outer(tt, tt, "-"))

  # Scale to maintain SNR
  signal_var <- var(as.vector(eta_mat))
  error_scale <- sqrt(signal_var / snr)

  L <- t(chol(R)) * error_scale
  epsilon <- matrix(rnorm(n * ny), nrow = n) %*% t(L)

  dat$Y <- eta_mat + epsilon

  # Extract the beta(s,t) matrix — key is "ff(X1, xind = s)"
  beta_key <- grep("^ff\\(", names(truth$beta), value = TRUE)
  beta_mat <- if (length(beta_key) > 0) truth$beta[[beta_key[1]]] else NULL

  list(
    data = dat,
    s_grid = s_grid,
    t_grid = t_grid,
    truth = truth,
    beta_mat = beta_mat
  )
}

#' Interpolate ff truth beta(s,t) to coef evaluation grid and center
interpolate_ff_truth <- function(
  beta_mat,
  s_grid,
  t_grid,
  x_eval,
  y_eval,
  fit
) {
  if (is.null(beta_mat) || !is.matrix(beta_mat)) return(NULL)

  # Bilinear interpolation
  beta_eval <- outer(
    x_eval,
    y_eval,
    Vectorize(function(s, t) {
      si <- findInterval(s, s_grid, all.inside = TRUE)
      ti <- findInterval(t, t_grid, all.inside = TRUE)
      s_frac <- (s - s_grid[si]) / (s_grid[si + 1] - s_grid[si])
      t_frac <- (t - t_grid[ti]) / (t_grid[ti + 1] - t_grid[ti])
      (1 - s_frac) *
        (1 - t_frac) *
        beta_mat[si, ti] +
        s_frac * (1 - t_frac) * beta_mat[si + 1, ti] +
        (1 - s_frac) * t_frac * beta_mat[si, ti + 1] +
        s_frac * t_frac * beta_mat[si + 1, ti + 1]
    })
  )

  # Center to match pffr identifiability constraints
  # sum_s w(s) * beta(s, t) = 0 for all t
  w_s <- ff_weights(x_eval, method = "simpson")
  beta_eval <- center_beta_ff(beta_eval, x_eval, w_s)

  as.vector(beta_eval)
}

#' Extract ff term estimates with all 4 sandwich methods
extract_ff_all_methods <- function(m, beta_mat, s_grid, t_grid) {
  methods <- list(
    nosandwich = list(sandwich = FALSE, freq = FALSE),
    hc = list(sandwich = "hc", freq = FALSE),
    cluster_freq = list(sandwich = TRUE, freq = TRUE),
    cluster_bayes = list(sandwich = TRUE, freq = FALSE)
  )

  # Get baseline coefs first to determine eval grid and truth
  cc0 <- tryCatch(
    coef(m, sandwich = FALSE, seWithMean = FALSE, n2 = 25),
    error = function(e) NULL
  )
  if (is.null(cc0)) return(NULL)

  sm_names <- names(cc0$smterms)
  ff_idx <- grep("ff", sm_names, ignore.case = TRUE)
  if (length(ff_idx) == 0) return(NULL)

  ti0 <- cc0$smterms[[ff_idx[1]]]
  if (is.null(ti0) || is.null(ti0$coef) || is.null(ti0$x) || is.null(ti0$y)) {
    return(NULL)
  }

  # Compute truth on eval grid
  truth_vec <- interpolate_ff_truth(
    beta_mat,
    s_grid,
    t_grid,
    ti0$x,
    ti0$y,
    m
  )
  if (is.null(truth_vec) || length(truth_vec) != length(ti0$coef$value)) {
    return(NULL)
  }

  results <- list()
  for (mname in names(methods)) {
    args <- methods[[mname]]
    cc <- tryCatch(
      coef(
        m,
        sandwich = args$sandwich,
        freq = args$freq,
        seWithMean = FALSE,
        n2 = 25
      ),
      error = function(e) NULL
    )
    if (is.null(cc)) next

    ti <- cc$smterms[[ff_idx[1]]]
    if (is.null(ti) || is.null(ti$coef)) next

    est <- ti$coef$value
    se <- ti$coef$se
    if (is.null(se) || length(est) != length(truth_vec)) next

    results[[mname]] <- data.frame(
      method = mname,
      idx = seq_along(est),
      est = est,
      se = se,
      truth = truth_vec,
      lower = est - 1.96 * se,
      upper = est + 1.96 * se,
      stringsAsFactors = FALSE
    ) |>
      mutate(
        covered = truth >= lower & truth <= upper,
        width = upper - lower
      )
  }
  bind_rows(results)
}

run_scenario_c <- function(
  n_values = c(50, 100),
  rho_values = c(0, 0.3, 0.6, 0.9),
  B_rep = 100,
  seed = 789
) {
  settings <- expand.grid(
    n = n_values,
    rho = rho_values,
    stringsAsFactors = FALSE
  )

  cat("Scenario C:", nrow(settings), "settings x", B_rep, "reps\n")

  all_results <- list()

  for (s in seq_len(nrow(settings))) {
    n <- settings$n[s]
    rho <- settings$rho[s]

    cat(sprintf("[%d/%d] n=%d, rho=%.1f\n", s, nrow(settings), n, rho))

    for (b in seq_len(B_rep)) {
      sim <- tryCatch(
        simulate_ff_data(
          n = n,
          rho = rho,
          seed = seed * 1000 + (s - 1) * B_rep + b
        ),
        error = function(e) NULL
      )
      if (is.null(sim) || is.null(sim$beta_mat)) next

      # s_grid must be visible in the formula environment
      s_grid <- sim$s_grid
      m <- tryCatch(
        pffr(
          Y ~ ff(X1, xind = s_grid),
          yind = sim$t_grid,
          data = sim$data
        ),
        error = function(e) NULL
      )
      if (is.null(m)) next

      res <- tryCatch(
        extract_ff_all_methods(m, sim$beta_mat, sim$s_grid, sim$t_grid),
        error = function(e) NULL
      )
      if (is.null(res) || nrow(res) == 0) next

      res$n <- n
      res$rho <- rho
      res$rep <- b

      all_results[[length(all_results) + 1]] <- res
    }
  }

  raw <- bind_rows(all_results)

  if (nrow(raw) == 0) {
    cat("WARNING: No ff term results collected!\n")
    return(list(raw = raw, summary = tibble()))
  }

  summary_c <- raw |>
    filter(!is.na(covered)) |>
    group_by(n, rho, method) |>
    summarise(
      coverage = mean(covered),
      mean_width = mean(width),
      se_ratio = if (sd(est - truth, na.rm = TRUE) > 0) {
        mean(se, na.rm = TRUE) / sd(est - truth, na.rm = TRUE)
      } else {
        NA_real_
      },
      .groups = "drop"
    )

  list(raw = raw, summary = summary_c)
}

# =========================================================================
# Main entry point
# =========================================================================

if (!interactive()) {
  outdir <- "ci-benchmark/results"
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  # --- Scenario A ---
  res_a <- run_scenario_a(
    n_values = c(50, 100),
    var_ratios = c(5, 10, 20),
    hetero_shapes = c("linear", "bump"),
    B_rep = 200,
    seed = 123
  )

  cat("\n========== SCENARIO A: COVERAGE OVERALL ==========\n")
  print(
    res_a$summary_overall |>
      arrange(n, var_ratio, hetero_shape, method),
    n = 100
  )

  cat("\n========== SCENARIO A: COVERAGE BY REGION ==========\n")
  print(
    res_a$summary_by_region |>
      arrange(n, var_ratio, hetero_shape, method, region),
    n = 200
  )

  saveRDS(res_a, file.path(outdir, "sandwich-scenarioA.rds"))

  # --- Scenario B ---
  res_b <- run_scenario_b(
    n_values = c(50, 100),
    B_rep = 200,
    seed = 456
  )

  cat("\n========== SCENARIO B: COVERAGE (complex correlation) ==========\n")
  print(
    res_b$summary |>
      arrange(n, corr_type, method),
    n = 50
  )

  saveRDS(res_b, file.path(outdir, "sandwich-scenarioB.rds"))

  # --- Scenario C ---
  res_c <- run_scenario_c(
    n_values = c(50, 100),
    rho_values = c(0, 0.3, 0.6, 0.9),
    B_rep = 100,
    seed = 789
  )

  cat("\n========== SCENARIO C: ff TERM COVERAGE ==========\n")
  print(
    res_c$summary |>
      arrange(n, rho, method),
    n = 50
  )

  saveRDS(res_c, file.path(outdir, "sandwich-scenarioC.rds"))

  cat("\nAll results saved to", outdir, "\n")
}
