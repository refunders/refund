# ===========================================================================
# Sandwich Benchmark: Non-Gaussian Families (Poisson, Gamma)
# ===========================================================================
#
# Two scenarios:
#   D) Poisson with log link + correlated latent process on linear predictor
#   E) Gamma with log link + correlated latent process on linear predictor
#
# Key insight: for log-link families, a latent Gaussian process z_i(t) on the
# linear predictor scale preserves beta_1(t) in the marginal model (the
# intercept shifts by sigma_z^2 / 2, but the slope is unchanged). This means
# the working model Y ~ x, family=poisson/Gamma is correctly specified for
# the slope -- the sandwich corrects for overdispersion and autocorrelation.
#
# The score formula in gam_sandwich_cluster() is the general GLM score:
#   w = (d mu/d eta) * (y - mu) / (phi * V(mu))
# which is correct for all standard exponential families. No code changes
# are needed to support non-Gaussian pffr models.
#
# Usage:
#   Rscript ci-benchmark/sandwich-benchmark-nongaussian.R
# ===========================================================================

devtools::load_all(".")
library(mgcv)
library(dplyr)
library(tidyr)

source("ci-benchmark/benchmark-utils.R")

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

#' Fit pffr and extract coefficient CIs for all 4 sandwich methods
#'
#' @param m A fitted pffr model (no sandwich).
#' @param term_name Grep pattern for coef.pffr smterms name.
#' @param truth_vec Numeric vector of true values on the eval grid.
#' @returns Data frame with method, idx, est, se, covered, width.
extract_all_methods <- function(m, term_name, truth_vec) {
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
        n1 = 50
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
      dplyr::mutate(
        covered = truth >= lower & truth <= upper,
        width = upper - lower
      )
  }
  dplyr::bind_rows(results)
}

# =========================================================================
# SCENARIO D: Poisson with correlated latent process
# =========================================================================

cat("=== SCENARIO D: Poisson with Correlated Latent Process ===\n\n")

#' Simulate Poisson FOS data with correlated latent process
#'
#' DGP:
#'   eta_i(t) = beta_0(t) + x_i * beta_1(t) + z_i(t)
#'   Y_i(t) | eta_i(t) ~ Poisson(exp(eta_i(t)))
#'   z_i(t) ~ GP(0, sigma_z^2 * R(rho))
#'
#' The marginal model preserves beta_1(t) exactly (log link is collapsible).
#' Only beta_0 shifts by sigma_z^2/2.
#'
#' @param n Number of curves.
#' @param T_grid Number of time points.
#' @param rho AR(1) correlation parameter for latent process.
#' @param sigma_z SD of latent process (controls overdispersion).
#' @returns List with data, yind, truth.
simulate_poisson_corr <- function(n, T_grid = 40, rho = 0.6, sigma_z = 0.5) {
  t_grid <- seq(0, 1, length.out = T_grid)

  # Modest coefficients to avoid huge counts
  beta0 <- 0.5 * sin(2 * pi * t_grid)
  beta1 <- 0.3 * cos(2 * pi * t_grid)

  x <- rnorm(n)
  x <- x - mean(x)

  signal <- outer(rep(1, n), beta0) + outer(x, beta1)

  # Correlated latent process (AR(1))
  R <- rho^abs(outer(seq_len(T_grid), seq_len(T_grid), "-"))
  L <- t(chol(R)) * sigma_z
  z <- matrix(rnorm(n * T_grid), nrow = n) %*% t(L)

  eta <- signal + z
  mu <- exp(eta)
  Y <- matrix(
    rpois(n * T_grid, lambda = as.vector(mu)),
    nrow = n,
    ncol = T_grid
  )

  # True beta1 for the marginal model is the same beta1
  list(
    data = data.frame(Y = I(Y), x = x),
    yind = t_grid,
    truth = list(beta0 = beta0, beta1 = beta1)
  )
}

run_scenario_d <- function(
  n_values = c(50, 100, 200),
  rho_values = c(0, 0.3, 0.6, 0.9),
  sigma_z_values = c(0.3, 0.7),
  B_rep = 200,
  seed = 2026
) {
  set.seed(seed)

  settings <- expand.grid(
    n = n_values,
    rho = rho_values,
    sigma_z = sigma_z_values,
    stringsAsFactors = FALSE
  )

  cat("Scenario D:", nrow(settings), "settings x", B_rep, "reps\n")

  all_results <- list()

  for (s in seq_len(nrow(settings))) {
    n <- settings$n[s]
    rho <- settings$rho[s]
    sigma_z <- settings$sigma_z[s]

    cat(sprintf(
      "[%d/%d] n=%d, rho=%.1f, sigma_z=%.1f\n",
      s,
      nrow(settings),
      n,
      rho,
      sigma_z
    ))

    for (b in seq_len(B_rep)) {
      sim <- tryCatch(
        simulate_poisson_corr(n = n, rho = rho, sigma_z = sigma_z),
        error = function(e) NULL
      )
      if (is.null(sim)) next

      m <- tryCatch(
        pffr(Y ~ x, yind = sim$yind, data = sim$data, family = poisson()),
        error = function(e) NULL
      )
      if (is.null(m)) next

      # Get eval grid from baseline coefs
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
      res$rho <- rho
      res$sigma_z <- sigma_z
      res$rep <- b
      res$t <- t_vals[res$idx]

      all_results[[length(all_results) + 1]] <- res
    }
  }

  raw <- dplyr::bind_rows(all_results)

  if (nrow(raw) == 0) {
    cat("WARNING: No Poisson results collected!\n")
    return(list(raw = raw, summary = tibble()))
  }

  summary_d <- raw |>
    dplyr::filter(!is.na(covered)) |>
    dplyr::group_by(n, rho, sigma_z, method) |>
    dplyr::summarise(
      coverage = mean(covered),
      mean_width = mean(width),
      se_ratio = if (sd(est - truth, na.rm = TRUE) > 0) {
        mean(se, na.rm = TRUE) / sd(est - truth, na.rm = TRUE)
      } else {
        NA_real_
      },
      .groups = "drop"
    )

  list(raw = raw, summary = summary_d)
}

# =========================================================================
# SCENARIO E: Gamma with correlated latent process
# =========================================================================

cat("=== SCENARIO E: Gamma with Correlated Latent Process ===\n\n")

#' Simulate Gamma FOS data with correlated latent process
#'
#' DGP:
#'   eta_i(t) = beta_0(t) + x_i * beta_1(t) + z_i(t)
#'   Y_i(t) | eta_i(t) ~ Gamma(shape, rate = shape / exp(eta_i(t)))
#'   z_i(t) ~ GP(0, sigma_z^2 * R(rho))
#'
#' Gamma with log link is collapsible: beta_1(t) is preserved in the
#' marginal model.
#'
#' @param n Number of curves.
#' @param T_grid Number of time points.
#' @param rho AR(1) correlation parameter for latent process.
#' @param sigma_z SD of latent process.
#' @param shape Gamma shape parameter (lower = more dispersed).
#' @returns List with data, yind, truth.
simulate_gamma_corr <- function(
  n,
  T_grid = 40,
  rho = 0.6,
  sigma_z = 0.5,
  shape = 5
) {
  t_grid <- seq(0, 1, length.out = T_grid)

  # Positive intercept to keep mu well above 0
  beta0 <- 1 + 0.5 * sin(2 * pi * t_grid)
  beta1 <- 0.3 * cos(2 * pi * t_grid)

  x <- rnorm(n)
  x <- x - mean(x)

  signal <- outer(rep(1, n), beta0) + outer(x, beta1)

  # Correlated latent process (AR(1))
  R <- rho^abs(outer(seq_len(T_grid), seq_len(T_grid), "-"))
  L <- t(chol(R)) * sigma_z
  z <- matrix(rnorm(n * T_grid), nrow = n) %*% t(L)

  eta <- signal + z
  mu <- exp(eta)

  # Gamma: mean = mu, Var(Y) = mu^2 / shape
  Y <- matrix(
    rgamma(n * T_grid, shape = shape, rate = shape / as.vector(mu)),
    nrow = n,
    ncol = T_grid
  )

  list(
    data = data.frame(Y = I(Y), x = x),
    yind = t_grid,
    truth = list(beta0 = beta0, beta1 = beta1)
  )
}

run_scenario_e <- function(
  n_values = c(50, 100, 200),
  rho_values = c(0, 0.3, 0.6, 0.9),
  sigma_z_values = c(0.3, 0.7),
  shape = 5,
  B_rep = 200,
  seed = 2027
) {
  set.seed(seed)

  settings <- expand.grid(
    n = n_values,
    rho = rho_values,
    sigma_z = sigma_z_values,
    stringsAsFactors = FALSE
  )

  cat("Scenario E:", nrow(settings), "settings x", B_rep, "reps\n")

  all_results <- list()

  for (s in seq_len(nrow(settings))) {
    n <- settings$n[s]
    rho <- settings$rho[s]
    sigma_z <- settings$sigma_z[s]

    cat(sprintf(
      "[%d/%d] n=%d, rho=%.1f, sigma_z=%.1f\n",
      s,
      nrow(settings),
      n,
      rho,
      sigma_z
    ))

    for (b in seq_len(B_rep)) {
      sim <- tryCatch(
        simulate_gamma_corr(
          n = n,
          rho = rho,
          sigma_z = sigma_z,
          shape = shape
        ),
        error = function(e) NULL
      )
      if (is.null(sim)) next

      m <- tryCatch(
        pffr(
          Y ~ x,
          yind = sim$yind,
          data = sim$data,
          family = Gamma(link = "log")
        ),
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
      res$rho <- rho
      res$sigma_z <- sigma_z
      res$shape <- shape
      res$rep <- b
      res$t <- t_vals[res$idx]

      all_results[[length(all_results) + 1]] <- res
    }
  }

  raw <- dplyr::bind_rows(all_results)

  if (nrow(raw) == 0) {
    cat("WARNING: No Gamma results collected!\n")
    return(list(raw = raw, summary = tibble()))
  }

  summary_e <- raw |>
    dplyr::filter(!is.na(covered)) |>
    dplyr::group_by(n, rho, sigma_z, method) |>
    dplyr::summarise(
      coverage = mean(covered),
      mean_width = mean(width),
      se_ratio = if (sd(est - truth, na.rm = TRUE) > 0) {
        mean(se, na.rm = TRUE) / sd(est - truth, na.rm = TRUE)
      } else {
        NA_real_
      },
      .groups = "drop"
    )

  list(raw = raw, summary = summary_e)
}

# =========================================================================
# Main entry point
# =========================================================================

if (!interactive()) {
  outdir <- "ci-benchmark/results"
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  # --- Scenario D: Poisson ---
  res_d <- run_scenario_d(
    n_values = c(50, 100, 200),
    rho_values = c(0, 0.3, 0.6, 0.9),
    sigma_z_values = c(0.3, 0.7),
    B_rep = 200,
    seed = 2026
  )

  cat("\n========== SCENARIO D: POISSON COVERAGE ==========\n")
  print(
    res_d$summary |>
      dplyr::arrange(sigma_z, n, rho, method),
    n = 200
  )

  saveRDS(res_d, file.path(outdir, "sandwich-scenarioD.rds"))

  # --- Scenario E: Gamma ---
  res_e <- run_scenario_e(
    n_values = c(50, 100, 200),
    rho_values = c(0, 0.3, 0.6, 0.9),
    sigma_z_values = c(0.3, 0.7),
    shape = 5,
    B_rep = 200,
    seed = 2027
  )

  cat("\n========== SCENARIO E: GAMMA COVERAGE ==========\n")
  print(
    res_e$summary |>
      dplyr::arrange(sigma_z, n, rho, method),
    n = 200
  )

  saveRDS(res_e, file.path(outdir, "sandwich-scenarioE.rds"))

  cat("\nAll results saved to", outdir, "\n")
}
