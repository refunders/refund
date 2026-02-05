# benchmark-utils.R
#
# Utility functions for pffr CI coverage benchmark:
# - Data generation helpers (error structures, covariate generation)
# - Visualization functions for benchmark results
#
# Source this file from confint-benchmark.R

# Null coalescing -------------------------------------------------------------

`%||%` <- function(x, y) if (is.null(x)) y else x

# Sandwich argument normalization ---------------------------------------------

#' Normalize sandwich argument to expected pffr conventions
#'
#' Several scripts historically used `TRUE/FALSE` for sandwich SEs. The main CI
#' benchmark uses explicit types: `"cluster"` and `"hc"`. This helper converts:
#' - `FALSE`/`NULL` -> `FALSE`
#' - `TRUE` -> `"cluster"`
#' - `"cluster"`/`"hc"` -> unchanged
#'
#' @param sandwich One of FALSE/NULL, TRUE, "cluster", "hc".
#' @returns `FALSE` or a character scalar ("cluster" or "hc").
normalize_sandwich_type <- function(sandwich) {
  if (is.null(sandwich) || isFALSE(sandwich)) return(FALSE)
  if (isTRUE(sandwich)) return("cluster")
  if (is.character(sandwich) && length(sandwich) == 1) {
    if (sandwich %in% c("cluster", "hc")) return(sandwich)
  }
  stop(
    "Invalid sandwich type. Use FALSE/NULL, TRUE, 'cluster', or 'hc'.",
    call. = FALSE
  )
}

# Error Structure Generation --------------------------------------------------

#' Create error covariance structure
#'
#' Generates correlation matrix and heteroskedasticity pattern for structured
#' errors in simulation studies.
#'
#' @param t_grid Numeric vector of response evaluation points.
#' @param corr_type One of "iid", "ar1", "gauss", "fourier", "fourier_pos".
#' @param corr_param Correlation parameter (e.g., rho for AR1).
#' @param hetero_type One of "none", "linear", "u", "bump".
#' @param hetero_param Heteroskedasticity amplitude parameter.
#' @returns List with t_grid, sigma_t, corr_mat, cov_mat, rho.
make_error_structure <- function(
  t_grid,
  corr_type = c("iid", "ar1", "gauss", "fourier", "fourier_pos"),
  corr_param = NULL,
  hetero_type = c("none", "linear", "u", "bump"),
  hetero_param = NULL
) {
  corr_type <- match.arg(corr_type)
  hetero_type <- match.arg(hetero_type)
  ny <- length(t_grid)
  dist_mat <- abs(outer(t_grid, t_grid, "-"))
  corr_param_used <- corr_param

  make_pd <- function(mat) {
    # Check if already PD (fast path)
    ok <- tryCatch(
      {
        chol(mat)
        TRUE
      },
      error = function(e) FALSE
    )
    if (ok) return(mat)

    # Project to nearest PD matrix (preserves structure better than jitter)
    pd <- as.matrix(Matrix::nearPD(mat, corr = FALSE, keepDiag = TRUE)$mat)

    # Add tiny jitter for numerical safety
    pd + diag(1e-10 * mean(diag(pd)), nrow(pd))
  }

  corr_mat <- switch(
    corr_type,
    iid = diag(1, ny),
    ar1 = {
      rho <- if (!is.null(corr_param)) {
        corr_param
      } else {
        # Choose rho so that Corr(t, t + 0.5 * range(t_grid)) ≈ 0.4.
        # For an equally-spaced grid this corresponds to k_half ≈ 0.5 / dt steps.
        dt <- mean(diff(t_grid))
        k_half <- max(1L, as.integer(round(0.5 / dt)))
        0.4^(1 / k_half)
      }
      corr_param_used <- rho
      rho^abs(outer(seq_len(ny), seq_len(ny), "-"))
    },
    gauss = {
      phi <- if (!is.null(corr_param)) {
        corr_param
      } else {
        # Choose phi so that Corr(t, t + 0.5 * range(t_grid)) ≈ 0.4 under the
        # squared-exponential kernel exp(-d^2/(2*phi^2)).
        d_half <- 0.5 * diff(range(t_grid))
        if (!is.finite(d_half) || d_half <= 0) {
          0.15
        } else {
          d_half / sqrt(2 * log(1 / 0.4))
        }
      }
      corr_param_used <- phi
      exp(-(dist_mat^2) / (2 * phi^2))
    },
    fourier = {
      # Low-rank periodic correlation:
      #   S_ij = cos(2*pi*(t_i - t_j)/period)
      # obtained via 2D Fourier features [sin, cos]. We add a diagonal "nugget"
      # and convert to correlation. The nugget shrinks all off-diagonals:
      #   Corr_ij = S_ij / (1 + nugget).
      #
      # Default nugget is chosen so that the RMS off-diagonal correlation is
      # ~0.4 (comparable strength to other defaults), while retaining a
      # sinusoidal shape.
      period <- 0.4

      B <- cbind(
        sin(2 * pi * t_grid / period),
        cos(2 * pi * t_grid / period)
      )
      S <- B %*% t(B)

      nugget <- if (!is.null(corr_param)) {
        corr_param
      } else if (ny < 2) {
        0
      } else {
        r_target <- 0.4
        off <- S[upper.tri(S, diag = FALSE)]
        rms0 <- sqrt(mean(off^2))
        if (!is.finite(rms0) || rms0 <= 0) {
          0
        } else {
          max(0, rms0 / r_target - 1)
        }
      }
      corr_param_used <- nugget

      base_cov <- S + diag(nugget, ny)
      stats::cov2cor(base_cov)
    },
    fourier_pos = {
      # Non-monotonic correlation clamped to be non-negative:
      #   raw(d) = cos(2*pi*d/period), corr(d) = max(0, raw(d))
      # This gives a periodic correlation that never goes negative,
      # creating a non-monotonic but all-positive dependence structure.
      period <- if (!is.null(corr_param)) corr_param else 0.4
      corr_param_used <- period
      raw_corr <- cos(2 * pi * dist_mat / period)
      mat <- pmax(raw_corr, 0)
      diag(mat) <- 1
      mat
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
    corr_param = corr_param_used,
    hetero_type = hetero_type,
    hetero_param = hetero_param,
    sigma_t = sigma_t,
    corr_mat = corr_mat,
    cov_mat = cov_mat,
    rho = if (corr_type == "ar1") corr_param_used else NA_real_
  )
}

#' Sample multivariate errors
#'
#' @param n Number of observations.
#' @param cov_mat Covariance matrix (ny x ny).
#' @param error_dist One of "gaussian", "t6".
#' @returns Matrix (n x ny) of error samples.
sample_errors <- function(n, cov_mat, error_dist = c("gaussian", "t6")) {
  error_dist <- match.arg(error_dist)
  ny <- nrow(cov_mat)
  stopifnot(ncol(cov_mat) == ny)

  if (error_dist == "gaussian") {
    return(mvtnorm::rmvnorm(n = n, mean = rep(0, ny), sigma = cov_mat))
  }

  df <- 6
  sigma_scale <- cov_mat * (df - 2) / df

  mvtnorm::rmvt(n = n, sigma = sigma_scale, df = df, delta = rep(0, ny))
}

#' Compute scale factor to achieve target SNR
#'
#' @param eps Error matrix (n x ny).
#' @param eta Signal matrix (n x ny).
#' @param snr_target Target signal-to-noise ratio.
#' @returns Numeric scale factor for errors.
compute_snr_scale_factor <- function(eps, eta, snr_target) {
  var_eta <- stats::var(as.vector(eta))
  var_eps <- stats::var(as.vector(eps))
  if (!is.finite(var_eta) || var_eta <= 0) var_eta <- 1
  if (!is.finite(var_eps) || var_eps <= 0) var_eps <- 1
  sqrt(var_eta / (snr_target * var_eps))
}

# Covariate Generation --------------------------------------------------------

#' Generate random smooth functions
#'
#' @param n Number of curves to generate.
#' @param grid Evaluation points.
#' @param bs_dim B-spline basis dimension.
#' @returns Matrix (n x length(grid)).
random_function_matrix <- function(n, grid, bs_dim = 25L) {
  stopifnot(is.numeric(n), length(n) == 1, n > 0)
  stopifnot(is.numeric(grid), length(grid) > 1)
  stopifnot(is.numeric(bs_dim), length(bs_dim) == 1, bs_dim >= 5)

  X <- splines::bs(grid, df = bs_dim, intercept = TRUE)
  coef_mat <- matrix(rnorm(n * bs_dim), nrow = bs_dim, ncol = n)
  t(X %*% coef_mat)
}

# Integration Weights ---------------------------------------------------------

#' Compute integration weights for ff terms
#'
#' @param xind Covariate evaluation points.
#' @param method One of "simpson", "trapezoidal", "riemann".
#' @returns Numeric vector of weights.
ff_weights <- function(xind, method = c("simpson", "trapezoidal", "riemann")) {
  method <- match.arg(method)
  nx <- length(xind)
  if (nx < 2) return(1)

  if (method == "simpson") {
    # Match pffr's formula: ((b-a)/n)/3 * [1, 4, 2, ..., 1] (divides by n, not n-1)
    w <- ((xind[nx] - xind[1]) / nx) /
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

  # riemann
  diffs <- diff(xind)
  c(mean(diffs), diffs)
}

#' Center beta(s,t) surface for ff terms
#'
#' Enforces sum_s w(s) * beta(s,t) = 0 for all t.
#'
#' @param beta_st Matrix (ns x nt).
#' @param s_grid s evaluation points.
#' @param weights Integration weights.
#' @returns Centered beta matrix.
center_beta_ff <- function(beta_st, s_grid, weights = NULL) {
  w <- weights %||% ff_weights(s_grid, method = "simpson")
  w_sum <- sum(w)
  col_means <- colSums(beta_st * w) / w_sum
  sweep(beta_st, 2, col_means, "-")
}

# Covariance Estimation from Fit ----------------------------------------------

#' Extract fitted mean matrix
fitted_mean_matrix <- function(fit) {
  fm <- tryCatch(fitted(fit, which = "mean"), error = function(e) fitted(fit))
  if (is.list(fm) && !is.null(fm$mean)) return(fm$mean)
  fm
}

#' Estimate residual covariance matrix from pffr fit
#'
#' @param fit A pffr fit object.
#' @param dat The data used for fitting.
#' @returns Covariance matrix (nyindex x nyindex).
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
  hatSigma + diag(ridge, nrow(hatSigma))
}

#' Estimate AR(1) correlation from pffr fit
#'
#' @param fit A pffr fit object.
#' @param dat The data used for fitting.
#' @param clamp Maximum absolute value for rho.
#' @returns Estimated rho.
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
  max(-clamp, min(clamp, rho))
}

# Visualization Functions -----------------------------------------------------

#' Plot coverage summary
#'
#' @param summary_df Summary data frame from summarize_benchmark().
#' @param nominal Nominal coverage level.
#' @returns ggplot object.
plot_coverage_summary <- function(summary_df, nominal = 0.90) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("ggplot2 not installed; skipping plots.")
    return(NULL)
  }

  ggplot2::ggplot(
    summary_df,
    ggplot2::aes(x = method, y = mean_coverage, color = term_type)
  ) +
    ggplot2::geom_hline(yintercept = nominal, linetype = 2, color = "gray40") +
    ggplot2::geom_point(
      position = ggplot2::position_dodge(width = 0.4),
      size = 2.5
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = mean_coverage - 1.96 * se_coverage,
        ymax = mean_coverage + 1.96 * se_coverage
      ),
      position = ggplot2::position_dodge(width = 0.4),
      width = 0.2
    ) +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
    ggplot2::labs(
      title = "CI Coverage by Method and Term Type",
      y = "Empirical Coverage",
      x = "Method",
      color = "Term Type"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' Plot coverage vs SNR
#'
#' @param summary_df Summary data frame with snr column.
#' @param nominal Nominal coverage level.
#' @returns ggplot object.
plot_coverage_vs_snr <- function(summary_df, nominal = 0.90) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("ggplot2 not installed; skipping plots.")
    return(NULL)
  }

  ggplot2::ggplot(
    summary_df,
    ggplot2::aes(
      x = factor(snr),
      y = mean_coverage,
      color = method,
      shape = term_type
    )
  ) +
    ggplot2::geom_hline(yintercept = nominal, linetype = 2, color = "gray40") +
    ggplot2::geom_point(
      size = 3,
      position = ggplot2::position_dodge(width = 0.3)
    ) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      title = "Coverage vs SNR",
      x = "SNR",
      y = "Coverage"
    ) +
    ggplot2::theme_bw()
}

#' Plot CI width comparison
#'
#' @param summary_df Summary data frame.
#' @returns ggplot object.
plot_width_comparison <- function(summary_df) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("ggplot2 not installed; skipping plots.")
    return(NULL)
  }

  ggplot2::ggplot(
    summary_df,
    ggplot2::aes(x = method, y = mean_width, fill = term_type)
  ) +
    ggplot2::geom_col(position = ggplot2::position_dodge()) +
    ggplot2::labs(
      title = "Mean CI Width by Method and Term",
      x = "Method",
      y = "Mean CI Width",
      fill = "Term Type"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' Plot timing comparison
#'
#' @param timing_df Timing summary data frame.
#' @returns ggplot object.
plot_timing_comparison <- function(timing_df) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("ggplot2 not installed; skipping plots.")
    return(NULL)
  }

  ggplot2::ggplot(
    timing_df,
    ggplot2::aes(x = method, y = mean_fit_time_total, fill = method)
  ) +
    ggplot2::geom_col() +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = mean_fit_time_total - sd_fit_time_total,
        ymax = mean_fit_time_total + sd_fit_time_total
      ),
      width = 0.2
    ) +
    ggplot2::labs(
      title = "Fit Time by Method (includes pilot for GLS/AR)",
      x = "Method",
      y = "Time (seconds)"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none")
}

#' Plot Z-score distribution
#'
#' @param results Raw results with est, truth, se columns.
#' @returns ggplot object.
plot_zscore_dist <- function(results) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("ggplot2 not installed; skipping plots.")
    return(NULL)
  }

  if (!"z_score" %in% names(results)) {
    results$z_score <- (results$est - results$truth) / results$se
  }

  ggplot2::ggplot(results, ggplot2::aes(x = z_score)) +
    ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(density)),
      bins = 50,
      fill = "steelblue",
      alpha = 0.7
    ) +
    ggplot2::stat_function(fun = dnorm, color = "red", linewidth = 1) +
    ggplot2::facet_wrap(~term_type) +
    ggplot2::labs(
      title = "Z-score Distribution vs Standard Normal",
      x = "Z-score = (estimate - truth) / SE",
      y = "Density"
    ) +
    ggplot2::theme_bw() +
    ggplot2::xlim(-5, 5)
}

# Family helpers --------------------------------------------------------------

#' Get family object from string
#'
#' @param fit_family One of "gaussian", "gaulss", "scat".
#' @returns Family object.
fit_family_object <- function(fit_family = c("gaussian", "gaulss", "scat")) {
  fit_family <- match.arg(fit_family)
  switch(
    fit_family,
    gaussian = stats::gaussian(),
    gaulss = mgcv::gaulss(),
    scat = mgcv::scat()
  )
}
