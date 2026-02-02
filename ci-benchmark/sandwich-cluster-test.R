# ===========================================================================
# Cluster-Robust Sandwich Estimator — Monte Carlo Benchmark
# ===========================================================================
#
# Compares coverage and covariance calibration of:
#   1. pffr (no sandwich)          — baseline
#   2. pffr(sandwich = "hc")       — mgcv observation-level HC sandwich
#   3. pffr(sandwich = TRUE)       — cluster-robust sandwich (freq, B2=0)
#   4. pffr(sandwich = TRUE) Bayes — cluster-robust sandwich (B2 = Vp-Ve)
#
# DGP: Y_i(t) = beta_0(t) + x_i * beta_1(t) + epsilon_i(t)
# with AR(1) and heteroskedastic within-curve errors.
#
# Usage:
#   Rscript ci-benchmark/sandwich-cluster-test.R
#   # or source interactively and call run_sandwich_benchmark()
# ===========================================================================

library(refund)
library(mgcv)
library(dplyr)
library(tidyr)
library(ggplot2)

# ---------------------------------------------------------------------------
# 1. Error generation
# ---------------------------------------------------------------------------

#' Generate AR(1) covariance matrix with heteroskedasticity
#'
#' @param T_grid Integer, number of time points.
#' @param rho AR(1) correlation parameter.
#' @param sigma_t Numeric vector of length T_grid giving marginal SDs at each t.
#' @returns T_grid x T_grid covariance matrix.
make_ar1_cov <- function(T_grid, rho, sigma_t = rep(1, T_grid)) {
  tt <- seq_len(T_grid)
  R <- rho^abs(outer(tt, tt, "-"))
  D <- diag(sigma_t)
  D %*% R %*% D
}

#' Generate heteroskedasticity patterns
#'
#' @param T_grid Number of time points.
#' @param type One of "none", "linear", "ushape".
#' @returns Numeric vector of marginal SDs.
make_hetero_pattern <- function(T_grid, type = "none") {
  t_grid <- seq(0, 1, length.out = T_grid)
  switch(
    type,
    none = rep(1, T_grid),
    linear = 0.5 + 1.5 * t_grid,
    ushape = 0.5 + 2 * (t_grid - 0.5)^2,
    stop("Unknown hetero type: ", type)
  )
}

# ---------------------------------------------------------------------------
# 2. DGP
# ---------------------------------------------------------------------------

#' Simulate one dataset from the function-on-scalar DGP
#'
#' @param n Number of curves.
#' @param T_grid Number of time points per curve.
#' @param rho AR(1) correlation.
#' @param hetero_type Heteroskedasticity type.
#' @param snr Signal-to-noise ratio.
#' @returns List with data, yind, truth (beta0, beta1).
simulate_fos_data <- function(
  n,
  T_grid = 40,
  rho = 0,
  hetero_type = "none",
  snr = 5
) {
  t_grid <- seq(0, 1, length.out = T_grid)

  # True coefficients
  beta0 <- sin(2 * pi * t_grid)
  beta1 <- cos(2 * pi * t_grid)

  # Scalar covariate
  x <- rnorm(n)
  x <- x - mean(x) # center for identifiability

  # Signal
  signal <- outer(rep(1, n), beta0) + outer(x, beta1)

  # Error covariance
  sigma_t <- make_hetero_pattern(T_grid, hetero_type)
  # Scale errors so that Var(signal) / Var(error) ≈ snr
  signal_var <- var(as.vector(signal))
  base_error_var <- mean(sigma_t^2) # average marginal variance
  error_scale <- sqrt(signal_var / (snr * base_error_var))
  sigma_t_scaled <- sigma_t * error_scale

  Sigma <- make_ar1_cov(T_grid, rho, sigma_t_scaled)

  # Generate correlated errors via Cholesky
  L <- t(chol(Sigma))
  epsilon <- matrix(rnorm(n * T_grid), nrow = n) %*% t(L)

  Y <- signal + epsilon

  list(
    data = data.frame(Y = I(Y), x = x),
    yind = t_grid,
    truth = list(beta0 = beta0, beta1 = beta1)
  )
}

# ---------------------------------------------------------------------------
# 3. Fit and extract
# ---------------------------------------------------------------------------

#' Fit pffr and extract coefficient estimates + SEs for beta_1(t)
#'
#' @param sim_data Output from simulate_fos_data().
#' @returns Data frame with columns: method, t, beta1_hat, se, lower, upper,
#'   beta1_true.
fit_and_extract <- function(sim_data) {
  d <- sim_data$data
  yind <- sim_data$yind
  truth <- sim_data$truth

  # Fit base model (no sandwich)
  m <- tryCatch(
    pffr(Y ~ x, yind = yind, data = d),
    error = function(e) NULL
  )
  if (is.null(m)) return(NULL)

  results <- list()

  # Method 1: no sandwich (default Vc/Vp)
  c_nosw <- tryCatch(coef(m, sandwich = FALSE), error = function(e) NULL)
  if (!is.null(c_nosw)) {
    results$nosandwich <- extract_beta1_from_coef(c_nosw, truth, "nosandwich")
  }

  # Method 2: HC sandwich (observation-level)
  c_hc <- tryCatch(coef(m, sandwich = "hc"), error = function(e) NULL)
  if (!is.null(c_hc)) {
    results$hc <- extract_beta1_from_coef(c_hc, truth, "hc")
  }

  # Method 3: cluster-robust sandwich (frequentist, B2=0)
  c_cl_freq <- tryCatch(
    coef(m, sandwich = TRUE, freq = TRUE),
    error = function(e) NULL
  )
  if (!is.null(c_cl_freq)) {
    results$cluster_freq <- extract_beta1_from_coef(
      c_cl_freq,
      truth,
      "cluster_freq"
    )
  }

  # Method 4: cluster-robust sandwich (Bayesian, B2=Vp-Ve)
  c_cl_bayes <- tryCatch(
    coef(m, sandwich = TRUE, freq = FALSE),
    error = function(e) NULL
  )
  if (!is.null(c_cl_bayes)) {
    results$cluster_bayes <- extract_beta1_from_coef(
      c_cl_bayes,
      truth,
      "cluster_bayes"
    )
  }

  # Also extract beta0 (intercept)
  c_nosw_int <- extract_beta0_from_coef(c_nosw, truth, "nosandwich")
  c_cl_bayes_int <- extract_beta0_from_coef(c_cl_bayes, truth, "cluster_bayes")

  # Combine beta1 results
  beta1_df <- bind_rows(results)

  # Extract Vp diagonal and beta1_hat for covariance comparison
  beta1_hat_vec <- if (!is.null(c_nosw)) {
    smterms <- c_nosw$smterms
    x_term <- smterms[[grep("^x\\(", names(smterms))[1]]]
    if (!is.null(x_term)) x_term$coef[, "value"] else NULL
  } else {
    NULL
  }

  list(
    beta1 = beta1_df,
    beta0_nosandwich = c_nosw_int,
    beta0_cluster = c_cl_bayes_int,
    beta1_hat = beta1_hat_vec
  )
}


#' Extract beta_1(t) estimates from coef.pffr output
#'
#' @param coef_obj Output from coef.pffr().
#' @param truth List with beta0, beta1 vectors.
#' @param method_name Character label for this method.
#' @returns Data frame with t, beta1_hat, se, lower, upper, beta1_true, method.
extract_beta1_from_coef <- function(coef_obj, truth, method_name) {
  smterms <- coef_obj$smterms

  # Find the x(yindex) term — it's the varying coefficient for x
  x_idx <- grep("^x\\(", names(smterms))
  if (length(x_idx) == 0) return(NULL)
  x_term <- smterms[[x_idx[1]]]

  coef_mat <- x_term$coef
  t_vals <- x_term$x # evaluation grid for t

  # Interpolate truth to evaluation grid
  beta1_true <- approx(
    seq(0, 1, length.out = length(truth$beta1)),
    truth$beta1,
    xout = t_vals
  )$y

  data.frame(
    t = t_vals,
    beta1_hat = coef_mat[, "value"],
    se = if ("se" %in% colnames(coef_mat)) coef_mat[, "se"] else NA_real_,
    beta1_true = beta1_true,
    method = method_name,
    stringsAsFactors = FALSE
  ) |>
    mutate(
      lower = beta1_hat - 1.96 * se,
      upper = beta1_hat + 1.96 * se,
      covered = beta1_true >= lower & beta1_true <= upper,
      width = upper - lower
    )
}

#' Extract beta_0(t) (intercept) estimates from coef.pffr output
extract_beta0_from_coef <- function(coef_obj, truth, method_name) {
  if (is.null(coef_obj)) return(NULL)
  smterms <- coef_obj$smterms

  int_idx <- grep("Intercept", names(smterms))
  if (length(int_idx) == 0) return(NULL)
  int_term <- smterms[[int_idx[1]]]

  coef_mat <- int_term$coef
  t_vals <- int_term$x

  beta0_true <- approx(
    seq(0, 1, length.out = length(truth$beta0)),
    truth$beta0,
    xout = t_vals
  )$y

  data.frame(
    t = t_vals,
    beta0_hat = coef_mat[, "value"],
    se = if ("se" %in% colnames(coef_mat)) coef_mat[, "se"] else NA_real_,
    beta0_true = beta0_true,
    method = method_name,
    stringsAsFactors = FALSE
  ) |>
    mutate(
      lower = beta0_hat - 1.96 * se,
      upper = beta0_hat + 1.96 * se,
      covered = beta0_true >= lower & beta0_true <= upper,
      width = upper - lower
    )
}

# ---------------------------------------------------------------------------
# 4. Monte Carlo experiment
# ---------------------------------------------------------------------------

#' Run full Monte Carlo benchmark
#'
#' @param n_values Sample sizes to test.
#' @param rho_values AR(1) correlations to test.
#' @param hetero_types Heteroskedasticity patterns to test.
#' @param B_rep Number of MC replications per setting.
#' @param T_grid Number of time points.
#' @param snr Signal-to-noise ratio.
#' @param seed Random seed.
#' @returns List with raw results and summary tables.
#' @export
run_sandwich_benchmark <- function(
  n_values = c(20, 50, 100, 200),
  rho_values = c(0, 0.3, 0.6, 0.9),
  hetero_types = c("none", "linear", "ushape"),
  B_rep = 500,
  T_grid = 40,
  snr = 5,
  seed = 42
) {
  set.seed(seed)

  # Build experiment grid
  settings <- expand.grid(
    n = n_values,
    rho = rho_values,
    hetero = hetero_types,
    stringsAsFactors = FALSE
  )

  cat(
    "Running",
    nrow(settings),
    "settings x",
    B_rep,
    "replications =",
    nrow(settings) * B_rep,
    "total fits\n"
  )

  all_beta1 <- list()
  all_beta0 <- list()
  all_beta1_hats <- list() # for covariance comparison

  for (s in seq_len(nrow(settings))) {
    n <- settings$n[s]
    rho <- settings$rho[s]
    hetero <- settings$hetero[s]

    cat(sprintf(
      "[%d/%d] n=%d, rho=%.1f, hetero=%s ...\n",
      s,
      nrow(settings),
      n,
      rho,
      hetero
    ))

    setting_beta1_hats <- list()

    for (b in seq_len(B_rep)) {
      sim <- simulate_fos_data(
        n = n,
        T_grid = T_grid,
        rho = rho,
        hetero_type = hetero,
        snr = snr
      )

      res <- tryCatch(fit_and_extract(sim), error = function(e) NULL)
      if (is.null(res)) next

      if (!is.null(res$beta1) && nrow(res$beta1) > 0) {
        res$beta1$n <- n
        res$beta1$rho <- rho
        res$beta1$hetero <- hetero
        res$beta1$rep <- b
        all_beta1[[length(all_beta1) + 1]] <- res$beta1
      }

      # Collect beta0
      for (nm in c("beta0_nosandwich", "beta0_cluster")) {
        if (!is.null(res[[nm]]) && nrow(res[[nm]]) > 0) {
          res[[nm]]$n <- n
          res[[nm]]$rho <- rho
          res[[nm]]$hetero <- hetero
          res[[nm]]$rep <- b
          all_beta0[[length(all_beta0) + 1]] <- res[[nm]]
        }
      }

      # Collect beta1_hat for empirical covariance
      if (!is.null(res$beta1_hat)) {
        setting_beta1_hats[[b]] <- res$beta1_hat
      }
    }

    all_beta1_hats[[paste(n, rho, hetero, sep = "_")]] <- setting_beta1_hats
  }

  raw_beta1 <- bind_rows(all_beta1)
  raw_beta0 <- bind_rows(all_beta0)

  # ---------------------------------------------------------------------------
  # 5. Compute summaries
  # ---------------------------------------------------------------------------

  # Level 1: Bias diagnostics
  bias_summary <- raw_beta1 |>
    group_by(n, rho, hetero, method) |>
    summarise(
      mean_bias = mean(beta1_hat - beta1_true, na.rm = TRUE),
      mean_abs_bias = mean(abs(beta1_hat - beta1_true), na.rm = TRUE),
      .groups = "drop"
    )

  # Level 3: Coverage and width
  coverage_summary <- raw_beta1 |>
    filter(!is.na(covered)) |>
    group_by(n, rho, hetero, method) |>
    summarise(
      coverage = mean(covered, na.rm = TRUE),
      mean_width = mean(width, na.rm = TRUE),
      mc_se = sqrt(coverage * (1 - coverage) / n()),
      .groups = "drop"
    )

  # Beta0 coverage
  coverage_beta0 <- raw_beta0 |>
    filter(!is.na(covered)) |>
    group_by(n, rho, hetero, method) |>
    summarise(
      coverage = mean(covered, na.rm = TRUE),
      mean_width = mean(width, na.rm = TRUE),
      .groups = "drop"
    )

  # Level 2: Covariance calibration — empirical SD vs mean estimated SE
  se_calibration <- raw_beta1 |>
    group_by(n, rho, hetero, method, t) |>
    summarise(
      emp_sd = sd(beta1_hat, na.rm = TRUE),
      mean_se = mean(se, na.rm = TRUE),
      .groups = "drop"
    )

  list(
    raw_beta1 = raw_beta1,
    raw_beta0 = raw_beta0,
    bias_summary = bias_summary,
    coverage_summary = coverage_summary,
    coverage_beta0 = coverage_beta0,
    se_calibration = se_calibration,
    beta1_hats = all_beta1_hats,
    settings = settings
  )
}


# ---------------------------------------------------------------------------
# 6. Plotting
# ---------------------------------------------------------------------------

#' Plot 1: Bias check
plot_bias <- function(results) {
  results$bias_summary |>
    ggplot(aes(
      x = factor(rho),
      y = mean_bias,
      colour = method,
      group = method
    )) +
    geom_point(position = position_dodge(0.3)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_grid(hetero ~ n, labeller = label_both) +
    labs(
      title = "Pointwise bias of beta_1(t)",
      x = "rho",
      y = "Mean bias"
    ) +
    theme_bw()
}

#' Plot 2: SE calibration (estimated SE vs empirical SD)
plot_se_calibration <- function(results) {
  results$se_calibration |>
    ggplot(aes(x = emp_sd, y = mean_se, colour = method)) +
    geom_point(alpha = 0.3, size = 0.5) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    facet_grid(hetero ~ n, labeller = label_both) +
    labs(
      title = "SE calibration: estimated SE vs empirical SD",
      subtitle = "45-degree line = perfect calibration",
      x = "Empirical SD of beta1_hat(t)",
      y = "Mean estimated SE"
    ) +
    theme_bw()
}

#' Plot 3: Coverage vs rho
plot_coverage <- function(results) {
  results$coverage_summary |>
    ggplot(aes(x = rho, y = coverage, colour = method, group = method)) +
    geom_line() +
    geom_point() +
    geom_hline(yintercept = 0.95, linetype = "dashed") +
    geom_hline(yintercept = 0.90, linetype = "dotted", alpha = 0.5) +
    facet_grid(hetero ~ n, labeller = label_both) +
    labs(
      title = "Pointwise coverage of 95% CI for beta_1(t)",
      x = "rho",
      y = "Coverage"
    ) +
    coord_cartesian(ylim = c(0.5, 1)) +
    theme_bw()
}

#' Plot 4: CI width vs rho
plot_width <- function(results) {
  results$coverage_summary |>
    ggplot(aes(x = rho, y = mean_width, colour = method, group = method)) +
    geom_line() +
    geom_point() +
    facet_grid(hetero ~ n, labeller = label_both, scales = "free_y") +
    labs(
      title = "Mean CI width for beta_1(t)",
      x = "rho",
      y = "Mean width"
    ) +
    theme_bw()
}


# ---------------------------------------------------------------------------
# 7. Main entry point
# ---------------------------------------------------------------------------

if (!interactive()) {
  cat("=== Sandwich Cluster Benchmark ===\n\n")

  # Quick run for CI: fewer settings
  results <- run_sandwich_benchmark(
    n_values = c(50, 100),
    rho_values = c(0, 0.3, 0.6, 0.9),
    hetero_types = c("none", "linear"),
    B_rep = 200,
    seed = 42
  )

  cat("\n=== Coverage Summary ===\n")
  print(
    results$coverage_summary |>
      select(n, rho, hetero, method, coverage, mean_width) |>
      arrange(n, rho, hetero, method),
    n = 100
  )

  # Save results
  outdir <- "ci-benchmark/results"
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  saveRDS(results, file.path(outdir, "sandwich-cluster-results.rds"))

  # Save plots
  ggsave(
    file.path(outdir, "sandwich-bias.pdf"),
    plot_bias(results),
    width = 10,
    height = 6
  )
  ggsave(
    file.path(outdir, "sandwich-se-calibration.pdf"),
    plot_se_calibration(results),
    width = 10,
    height = 6
  )
  ggsave(
    file.path(outdir, "sandwich-coverage.pdf"),
    plot_coverage(results),
    width = 10,
    height = 6
  )
  ggsave(
    file.path(outdir, "sandwich-width.pdf"),
    plot_width(results),
    width = 10,
    height = 6
  )

  cat("\nResults saved to", outdir, "\n")
}
