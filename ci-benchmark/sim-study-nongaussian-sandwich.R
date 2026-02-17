# ===========================================================================
# Study 1: Non-Gaussian Sandwich with ff + linear Model
# ===========================================================================
#
# Evaluates cluster-robust sandwich CI coverage for pffr under:
#   - Poisson(log) and Binomial(logit) families
#   - Model: Y ~ ff(X1) + zlin
#   - Three correlation structures: IID, AR1(0.9), fourier_pos(0.3)
#   - Two sample sizes: n = 200, 400
#
# Key design decisions:
#   Poisson: eps additive on linear predictor, log link is collapsible
#     → conditional β = marginal β. Intercept shift by σ²_eps/2 corrected.
#   Binomial: Gaussian copula preserves marginal β (logit NOT collapsible).
#     → P(Y=1|X) = p exactly; within-curve correlation via copula.
#
# Methods: default (no sandwich), HC sandwich, cluster-robust sandwich.
#
# Usage:
#   Rscript ci-benchmark/sim-study-nongaussian-sandwich.R [mode]
#   mode: "smoke" (2 reps), "pilot" (10 reps), "full" (150 reps)
# ===========================================================================

# Setup -----------------------------------------------------------------------

library(tidyverse)

if (file.exists("DESCRIPTION")) {
  devtools::load_all(".", quiet = TRUE)
} else {
  library(refund)
}

library(mvtnorm)

source("ci-benchmark/benchmark-utils.R")

# Source confint-benchmark.R for shared metric functions
# (find_term_index, evaluate_truth_on_grid, compute_term_metrics, etc.)
# Guarded main block won't execute when sourced.
source("ci-benchmark/confint-benchmark.R")

# Constants -------------------------------------------------------------------

STUDY1_BASE_SEED <- 3001L
STUDY1_NXGRID <- 30L
STUDY1_NYGRID <- 40L
STUDY1_SIGMA_EPS <- 0.5 # Poisson latent error SD (on log scale)
STUDY1_MU_TARGET <- 4 # Poisson target mean count for average subject
STUDY1_P_TARGET <- 0.5 # Binomial target probability for average subject
STUDY1_K_FF <- c(12L, 12L) # Basis dimensions for ff term
STUDY1_K_YINDEX <- 12L # Basis dimension for varying coefficient terms

# Truth Generation Helpers ----------------------------------------------------

#' Generate random smooth 2D surface using tensor B-splines
#'
#' @param s_grid Evaluation points for first dimension.
#' @param t_grid Evaluation points for second dimension.
#' @param k_s Number of B-spline basis functions for s.
#' @param k_t Number of B-spline basis functions for t.
#' @param amplitude Target marginal SD of the surface values.
#' @returns Matrix (length(s_grid) x length(t_grid)).
generate_smooth_surface <- function(
  s_grid,
  t_grid,
  k_s = 6,
  k_t = 6,
  amplitude = 0.3
) {
  Bs <- splines::bs(s_grid, df = k_s, intercept = TRUE)
  Bt <- splines::bs(t_grid, df = k_t, intercept = TRUE)
  coef_mat <- matrix(rnorm(k_s * k_t), nrow = k_s, ncol = k_t)
  surface <- Bs %*% coef_mat %*% t(Bt)
  sd_surface <- sd(as.vector(surface))
  if (sd_surface > 0) surface <- surface / sd_surface * amplitude
  surface
}

#' Generate random smooth 1D function using B-splines
#'
#' @param grid Evaluation points.
#' @param k Number of B-spline basis functions.
#' @param amplitude Target SD of the function values.
#' @returns Numeric vector of length(grid).
generate_smooth_function <- function(grid, k = 6, amplitude = 0.3) {
  B <- splines::bs(grid, df = k, intercept = TRUE)
  coefs <- rnorm(k)
  f <- as.vector(B %*% coefs)
  sd_f <- sd(f)
  if (sd_f > 0) f <- f / sd_f * amplitude
  f
}

# DGP Simulation Functions ---------------------------------------------------

#' Simulate Poisson ff + linear data with correlated latent process
#'
#' DGP:
#'   eta_i(t) = beta0(t) + integral beta(s,t)*X1_i(s)ds + z_lin_i*beta_lin(t) + eps_i(t)
#'   Y_i(t) | eta ~ Poisson(exp(eta_i(t)))
#'   eps_i(t) ~ N(0, sigma_eps^2 * R)
#'
#' Log link is collapsible: conditional beta = marginal beta.
#' Intercept shifts by sigma_eps^2/2 in marginal model.
#'
#' @param n Number of observations (curves).
#' @param nxgrid Grid size for functional covariate.
#' @param nygrid Grid size for response.
#' @param corr_type Correlation type: "iid", "ar1", "fourier_pos".
#' @param corr_param Correlation parameter (rho for AR1, period for fourier_pos).
#' @param sigma_eps SD of latent Gaussian process on log scale.
#' @param seed RNG seed.
#' @returns List with data, s_grid, t_grid, truth, err_struct, terms.
simulate_poisson_ff_linear <- function(
  n,
  nxgrid = STUDY1_NXGRID,
  nygrid = STUDY1_NYGRID,
  corr_type = "iid",
  corr_param = NULL,
  sigma_eps = STUDY1_SIGMA_EPS,
  seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  s_grid <- seq(0, 1, length.out = nxgrid)
  t_grid <- seq(0, 1, length.out = nygrid)

  # Generate truth functions
  beta_ff <- generate_smooth_surface(
    s_grid,
    t_grid,
    k_s = 6,
    k_t = 6,
    amplitude = 0.3
  )
  # Center ff truth: sum_s w(s) * beta(s,t) = 0 for all t
  w_s <- ff_weights(s_grid)
  beta_ff <- center_beta_ff(beta_ff, s_grid, w_s)

  beta_lin <- generate_smooth_function(t_grid, k = 6, amplitude = 0.3)

  # Generate covariates
  X1 <- random_function_matrix(n, s_grid, bs_dim = 15L)
  z_lin <- rnorm(n)
  z_lin <- z_lin - mean(z_lin)

  # Compute ff integral: sum_s w(s) * beta(s,t) * X1_i(s)
  X1_w <- sweep(X1, 2, w_s, "*") # n x nxgrid
  ff_contrib <- X1_w %*% beta_ff # n x nygrid

  # Linear contribution
  lin_contrib <- outer(z_lin, beta_lin) # n x nygrid

  # Signal without intercept
  eta_signal <- ff_contrib + lin_contrib

  # Generate correlated errors on log scale
  err_struct <- make_error_structure(
    t_grid = t_grid,
    corr_type = corr_type,
    corr_param = corr_param,
    hetero_type = "none"
  )
  eps <- sample_errors(n, err_struct$cov_mat) * sigma_eps

  # Calibrate intercept for marginal model:
  # E[Y_i(t)|X_i,z_i] = exp(beta0(t) + ff_i(t) + lin_i(t) + sigma^2/2)
  # Set beta0 so average E[Y] ≈ mu_target
  # NOTE (R1-M2): This targets E[eta] on log scale. Due to Jensen's inequality,

  # E[exp(eta)] > exp(E[eta]) when signal variance is high. The sigma^2/2 term
  # corrects for eps variance but not for covariate signal variance. In practice,
  # achieved counts may drift slightly from mu_target when signal variance is large.
  # The DGP validation check (validate_poisson_ranges) confirms acceptable range.
  beta0 <- log(STUDY1_MU_TARGET) - colMeans(eta_signal) - sigma_eps^2 / 2

  # Full linear predictor (conditional, includes eps)
  eta_full <- sweep(eta_signal + eps, 2, beta0, "+")

  # Generate Poisson responses
  mu <- exp(eta_full)
  Y <- matrix(
    rpois(n * nygrid, lambda = as.vector(mu)),
    nrow = n,
    ncol = nygrid
  )

  # Package data for pffr
  dat <- data.frame(zlin = z_lin)
  dat$Y <- I(Y)
  dat$X1 <- I(X1)

  # Truth: marginal model coefficients
  # For Poisson log-link: slope coefficients are same conditional = marginal
  # Marginal intercept = beta0 + sigma_eps^2/2
  beta0_marginal <- beta0 + sigma_eps^2 / 2

  truth <- list(
    beta = list(
      `ff(X1)` = beta_ff, # nxgrid x nygrid matrix
      zlin = beta_lin, # nygrid vector
      intercept = beta0_marginal # nygrid vector (marginal)
    ),
    # eta for E(Y) comparison: marginal signal per subject (on link scale)
    eta = sweep(eta_signal, 2, beta0_marginal, "+"),
    sigma_eps = sigma_eps
  )

  list(
    data = dat,
    s_grid = s_grid,
    t_grid = t_grid,
    truth = truth,
    err_struct = err_struct,
    terms = c("ff", "linear")
  )
}

#' Simulate Binomial ff + linear data using Gaussian copula
#'
#' DGP (Gaussian copula preserves marginal β):
#'   p_i(t) = logit^{-1}(beta0(t) + integral beta(s,t)*X1_i(s)ds + z_lin_i*beta_lin(t))
#'   Z_i(t) ~ N(0, R) where R is the target correlation matrix
#'   Y_i(t) = 1{Phi(Z_i(t)) < p_i(t)}
#'
#' This preserves P(Y_i(t)=1 | X_i, z_lin_i) = p_i(t) exactly.
#' Within-curve correlation is controlled by R.
#'
#' @inheritParams simulate_poisson_ff_linear
#' @returns List with data, s_grid, t_grid, truth, err_struct, terms.
simulate_binomial_ff_linear <- function(
  n,
  nxgrid = STUDY1_NXGRID,
  nygrid = STUDY1_NYGRID,
  corr_type = "iid",
  corr_param = NULL,
  seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  s_grid <- seq(0, 1, length.out = nxgrid)
  t_grid <- seq(0, 1, length.out = nygrid)

  # Generate truth functions (moderate amplitude for logit scale)
  beta_ff <- generate_smooth_surface(
    s_grid,
    t_grid,
    k_s = 6,
    k_t = 6,
    amplitude = 0.4
  )
  w_s <- ff_weights(s_grid)
  beta_ff <- center_beta_ff(beta_ff, s_grid, w_s)

  beta_lin <- generate_smooth_function(t_grid, k = 6, amplitude = 0.3)

  # Generate covariates
  X1 <- random_function_matrix(n, s_grid, bs_dim = 15L)
  z_lin <- rnorm(n)
  z_lin <- z_lin - mean(z_lin)

  # Compute ff integral
  X1_w <- sweep(X1, 2, w_s, "*")
  ff_contrib <- X1_w %*% beta_ff

  # Linear contribution
  lin_contrib <- outer(z_lin, beta_lin)

  # Signal without intercept
  eta_signal <- ff_contrib + lin_contrib

  # Calibrate intercept: p_target = expit(beta0 + mean_signal)
  # logit(0.5) = 0, so beta0 = -colMeans(eta_signal)
  beta0 <- stats::qlogis(STUDY1_P_TARGET) - colMeans(eta_signal)

  # Full linear predictor (on logit scale)
  eta <- sweep(eta_signal, 2, beta0, "+")

  # Probability matrix
  p_mat <- stats::plogis(eta) # n x nygrid

  # Generate within-curve correlation via Gaussian copula
  err_struct <- make_error_structure(
    t_grid = t_grid,
    corr_type = corr_type,
    corr_param = corr_param,
    hetero_type = "none"
  )
  R <- err_struct$corr_mat # nygrid x nygrid

  Z <- mvtnorm::rmvnorm(n, mean = rep(0, nygrid), sigma = R)
  U <- pnorm(Z) # Marginally Uniform(0,1), correlated via R

  # Binary responses: Y_i(t) = 1{U_i(t) < p_i(t)}
  Y <- matrix(as.integer(U < p_mat), nrow = n, ncol = nygrid)

  # Package data
  dat <- data.frame(zlin = z_lin)
  dat$Y <- I(Y)
  dat$X1 <- I(X1)

  # Truth: marginal model coefficients (preserved by copula construction)
  truth <- list(
    beta = list(
      `ff(X1)` = beta_ff,
      zlin = beta_lin,
      intercept = beta0
    ),
    eta = eta # signal on logit scale per subject
  )

  list(
    data = dat,
    s_grid = s_grid,
    t_grid = t_grid,
    truth = truth,
    err_struct = err_struct,
    terms = c("ff", "linear")
  )
}

# DGP Settings ----------------------------------------------------------------

#' Create Study 1 DGP settings grid
#'
#' 2 families x 3 correlation structures x 2 sample sizes = 12 DGP cells.
#'
#' @returns Tibble with DGP configurations.
make_study1_settings <- function() {
  tidyr::crossing(
    family = c("poisson", "binomial"),
    corr_type = c("iid", "ar1", "fourier_pos"),
    n = c(200L, 400L)
  ) |>
    dplyr::mutate(
      corr_param = dplyr::case_when(
        corr_type == "ar1" ~ 0.9,
        corr_type == "fourier_pos" ~ 0.3,
        TRUE ~ NA_real_
      ),
      nxgrid = STUDY1_NXGRID,
      nygrid = STUDY1_NYGRID,
      dgp_id = dplyr::row_number()
    )
}

# Formula Building ------------------------------------------------------------

#' Build pffr formula for Study 1 (ff + linear)
#'
#' @param s_grid Covariate evaluation grid.
#' @param k_ff Basis dimensions for ff term.
#' @returns Formula object with environment containing s_grid.
build_study1_formula <- function(s_grid, k_ff = STUDY1_K_FF) {
  frml <- as.formula(sprintf(
    "Y ~ ff(X1, xind = s_grid, splinepars = list(bs = 'ps', k = c(%d, %d))) + zlin",
    k_ff[1],
    k_ff[2]
  ))
  e <- new.env(parent = globalenv())
  e$s_grid <- s_grid
  environment(frml) <- e
  frml
}

# Metric Extraction -----------------------------------------------------------

#' Extract metrics for all three methods from one pffr fit
#'
#' Fits the model once and extracts coefficients with different sandwich types.
#'
#' @param fit A fitted pffr model.
#' @param sim Simulation result from simulate_*_ff_linear().
#' @param alpha Significance level for CI.
#' @returns Tibble with metrics for each (method, term) combination.
extract_study1_metrics <- function(fit, sim, alpha = 0.10) {
  sandwich_types <- list(
    default = FALSE,
    hc = "hc",
    cluster = "cluster"
  )

  results <- list()
  for (method_name in names(sandwich_types)) {
    sw_type <- sandwich_types[[method_name]]

    for (term_type in c("ff", "linear")) {
      metrics <- tryCatch(
        compute_term_metrics(
          fit = fit,
          truth = sim$truth,
          term_type = term_type,
          alpha = alpha,
          use_sandwich = sw_type,
          s_grid = sim$s_grid,
          t_grid = sim$t_grid,
          data = sim$data,
          err_struct = sim$err_struct
        ),
        error = function(e) {
          warning(sprintf(
            "compute_term_metrics failed for method=%s term=%s: %s",
            method_name,
            term_type,
            conditionMessage(e)
          ))
          NULL
        }
      )

      if (!is.null(metrics)) {
        metrics$method <- method_name
        results[[length(results) + 1]] <- metrics
      }
    }

    # E(Y) metrics
    ey_metrics <- tryCatch(
      compute_ey_metrics(
        fit = fit,
        truth = sim$truth,
        alpha = alpha,
        use_sandwich = sw_type,
        err_struct = sim$err_struct,
        t_grid = sim$t_grid
      ),
      error = function(e) {
        warning(sprintf(
          "compute_ey_metrics failed for method=%s: %s",
          method_name,
          conditionMessage(e)
        ))
        NULL
      }
    )

    if (!is.null(ey_metrics)) {
      ey_metrics$method <- method_name
      results[[length(results) + 1]] <- ey_metrics
    }
  }

  dplyr::bind_rows(results)
}

# Runner Functions ------------------------------------------------------------

#' Run one Study 1 replicate
#'
#' @param row One-row tibble/list with DGP settings.
#' @param rep_id Replicate number.
#' @param alpha Significance level.
#' @returns Tibble with metrics, or NULL on failure.
run_study1_rep <- function(row, rep_id, alpha = 0.10) {
  if (inherits(row, "data.frame")) row <- as.list(row)

  seed <- STUDY1_BASE_SEED + 1000L * row$dgp_id + rep_id

  # Simulate data
  sim_fn <- if (row$family == "poisson") {
    simulate_poisson_ff_linear
  } else {
    simulate_binomial_ff_linear
  }

  sim <- sim_fn(
    n = row$n,
    nxgrid = row$nxgrid,
    nygrid = row$nygrid,
    corr_type = row$corr_type,
    corr_param = if (is.na(row$corr_param)) NULL else row$corr_param,
    seed = seed
  )

  # Build formula and fit
  frml <- build_study1_formula(sim$s_grid)
  bs_yindex <- list(bs = "ps", k = STUDY1_K_YINDEX, m = c(2, 1))

  fam <- if (row$family == "poisson") poisson() else binomial()

  t0 <- Sys.time()
  fit <- pffr(
    frml,
    yind = sim$t_grid,
    data = sim$data,
    family = fam,
    bs.yindex = bs_yindex
  )
  fit_time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  # Extract metrics
  metrics <- extract_study1_metrics(fit, sim, alpha = alpha)

  if (nrow(metrics) == 0) {
    # Return explicit failure row instead of silently dropping (R1-H1)
    return(tibble(
      term_type = NA_character_,
      coverage = NA_real_,
      coverage_high_var = NA_real_,
      coverage_low_var = NA_real_,
      mean_width = NA_real_,
      rmse = NA_real_,
      bias = NA_real_,
      mean_abs_error = NA_real_,
      mean_se = NA_real_,
      median_se = NA_real_,
      z_mean = NA_real_,
      z_sd = NA_real_,
      z_kurtosis = NA_real_,
      z2_mean = NA_real_,
      n_grid = NA_integer_,
      method = NA_character_,
      dgp_id = row$dgp_id,
      rep_id = rep_id,
      seed = seed,
      family = row$family,
      n = row$n,
      nxgrid = row$nxgrid,
      nygrid = row$nygrid,
      corr_type = row$corr_type,
      corr_param = row$corr_param,
      fit_time = fit_time,
      converged = FALSE,
      error_msg = "metrics extraction returned 0 rows"
    ))
  }

  # Annotate with settings
  metrics |>
    dplyr::mutate(
      dgp_id = row$dgp_id,
      rep_id = rep_id,
      seed = seed,
      family = row$family,
      n = row$n,
      nxgrid = row$nxgrid,
      nygrid = row$nygrid,
      corr_type = row$corr_type,
      corr_param = row$corr_param,
      fit_time = fit_time,
      converged = TRUE,
      error_msg = NA_character_
    )
}

#' Atomic save: write to tempfile then rename
#'
#' @param obj Object to save.
#' @param path Target file path.
atomic_saveRDS <- function(obj, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  tmp <- tempfile(tmpdir = dirname(path), fileext = ".rds.tmp")
  saveRDS(obj, tmp)
  if (!file.rename(tmp, path)) {
    warning("atomic_saveRDS: file.rename failed for ", path)
    # Fallback: try direct copy
    file.copy(tmp, path, overwrite = TRUE)
    unlink(tmp)
  }
}

#' Run Study 1 benchmark
#'
#' @param n_rep Number of replications per DGP cell.
#' @param parallel Use parallel processing?
#' @param n_workers Number of workers.
#' @param output_dir Output directory.
#' @param alpha Significance level.
#' @returns Combined results tibble.
run_study1 <- function(
  n_rep = 150L,
  parallel = TRUE,
  n_workers = max(1L, min(6L, parallel::detectCores() - 1L)),
  output_dir = "ci-benchmark/study1-nongaussian",
  alpha = 0.10
) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  settings <- make_study1_settings()

  # Expand to full grid
  grid <- settings |>
    tidyr::crossing(rep_id = seq_len(n_rep))

  cat("Study 1: Non-Gaussian Sandwich with ff + linear\n")
  cat("================================================\n")
  cat("DGP cells:", nrow(settings), "\n")
  cat("Replications:", n_rep, "\n")
  cat("Total fits:", nrow(grid), "\n")
  cat("Output dir:", output_dir, "\n\n")

  # Skip already-completed settings
  existing_files <- list.files(
    output_dir,
    pattern = "^dgp\\d+_rep\\d+\\.rds$",
    full.names = TRUE
  )
  if (length(existing_files) > 0) {
    existing_keys <- sub("\\.rds$", "", basename(existing_files))
    grid_keys <- sprintf("dgp%03d_rep%03d", grid$dgp_id, grid$rep_id)
    already_done <- grid_keys %in% existing_keys
    n_skip <- sum(already_done)
    if (n_skip > 0) {
      cat("Skipping", n_skip, "already-completed settings\n")
      grid <- grid[!already_done, , drop = FALSE]
    }
  }

  if (nrow(grid) == 0) {
    cat("All settings completed — loading from disk.\n")
    return(load_study1_results(output_dir))
  }

  # Worker function
  run_one <- function(row_df) {
    row <- as.list(row_df)
    result <- tryCatch(
      run_study1_rep(row, row$rep_id, alpha = alpha),
      error = function(e) {
        warning(sprintf(
          "run_study1_rep failed for dgp=%d rep=%d: %s",
          row$dgp_id,
          row$rep_id,
          conditionMessage(e)
        ))
        tibble(
          term_type = NA_character_,
          coverage = NA_real_,
          coverage_high_var = NA_real_,
          coverage_low_var = NA_real_,
          mean_width = NA_real_,
          rmse = NA_real_,
          bias = NA_real_,
          mean_abs_error = NA_real_,
          mean_se = NA_real_,
          median_se = NA_real_,
          z_mean = NA_real_,
          z_sd = NA_real_,
          z_kurtosis = NA_real_,
          z2_mean = NA_real_,
          n_grid = NA_integer_,
          method = NA_character_,
          dgp_id = row$dgp_id,
          rep_id = row$rep_id,
          seed = NA_integer_,
          family = row$family,
          n = row$n,
          nxgrid = row$nxgrid,
          nygrid = row$nygrid,
          corr_type = row$corr_type,
          corr_param = row$corr_param,
          fit_time = NA_real_,
          converged = FALSE,
          error_msg = conditionMessage(e)
        )
      }
    )

    # Always save, including failure rows (R1-H1)
    if (!is.null(result) && nrow(result) > 0) {
      save_path <- file.path(
        output_dir,
        sprintf("dgp%03d_rep%03d.rds", row$dgp_id, row$rep_id)
      )
      atomic_saveRDS(result, save_path)
    }
    gc()
    result
  }

  if (
    parallel &&
      requireNamespace("furrr", quietly = TRUE) &&
      requireNamespace("future", quietly = TRUE)
  ) {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multicore, workers = n_workers)

    is_dev_mode <- file.exists("DESCRIPTION")
    pkg_dir <- normalizePath(".")
    bench_utils_path <- normalizePath("ci-benchmark/benchmark-utils.R")
    confint_path <- normalizePath("ci-benchmark/confint-benchmark.R")
    study1_path <- normalizePath(
      "ci-benchmark/sim-study-nongaussian-sandwich.R"
    )

    rows <- split(grid, seq_len(nrow(grid)))

    results <- furrr::future_map_dfr(
      rows,
      function(row_df) {
        if (is_dev_mode) {
          devtools::load_all(pkg_dir, quiet = TRUE)
        } else {
          library(refund)
        }
        source(bench_utils_path, local = TRUE)
        source(confint_path, local = TRUE)
        source(study1_path, local = TRUE)
        run_one(row_df)
      },
      .options = furrr::furrr_options(seed = TRUE),
      .progress = TRUE
    )
  } else {
    results <- list()
    for (i in seq_len(nrow(grid))) {
      row_df <- grid[i, ]
      cat(sprintf(
        "\r[%d/%d] dgp=%d (family=%s, corr=%s, n=%d) rep=%d",
        i,
        nrow(grid),
        row_df$dgp_id,
        row_df$family,
        row_df$corr_type,
        row_df$n,
        row_df$rep_id
      ))
      results[[i]] <- run_one(row_df)
    }
    cat("\n")
    results <- dplyr::bind_rows(results)
  }

  # Load all results including previously completed
  all_results <- load_study1_results(output_dir)

  # Save combined file here only, not in load function (R1-M5)
  atomic_saveRDS(all_results, file.path(output_dir, "results_combined.rds"))
  cat("Total results:", nrow(all_results), "rows\n")

  all_results
}

#' Load Study 1 results from incremental saves
#'
#' @param output_dir Directory with per-rep RDS files.
#' @returns Combined tibble.
load_study1_results <- function(
  output_dir = "ci-benchmark/study1-nongaussian"
) {
  all_files <- list.files(
    output_dir,
    pattern = "^dgp\\d+_rep\\d+\\.rds$",
    full.names = TRUE
  )
  if (length(all_files) == 0) {
    warning("No result files found in ", output_dir)
    return(tibble())
  }
  results <- dplyr::bind_rows(lapply(all_files, readRDS))
  # Combined file saved only by run_study1(), not here (R1-M5: avoid race)
  results
}

# Summarization ---------------------------------------------------------------

#' Summarize Study 1 results
#'
#' @param results Results tibble from run_study1().
#' @returns Summary tibble.
summarize_study1 <- function(results) {
  # Compute failure counts BEFORE filtering NAs (R1-M1)
  failure_counts <- results |>
    dplyr::group_by(family, corr_type, corr_param, n, method, term_type) |>
    dplyr::summarise(
      n_total = dplyr::n(),
      n_failed = sum(!converged | is.na(coverage), na.rm = TRUE),
      .groups = "drop"
    )

  results |>
    dplyr::filter(!is.na(coverage)) |>
    dplyr::group_by(family, corr_type, corr_param, n, method, term_type) |>
    dplyr::summarise(
      mean_coverage = mean(coverage, na.rm = TRUE),
      se_coverage = sd(coverage, na.rm = TRUE) / sqrt(dplyr::n()),
      mean_width = mean(mean_width, na.rm = TRUE),
      mean_rmse = mean(rmse, na.rm = TRUE),
      mean_bias = mean(bias, na.rm = TRUE),
      mean_z_sd = mean(z_sd, na.rm = TRUE),
      mean_fit_time = mean(fit_time, na.rm = TRUE),
      n_reps = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::left_join(
      failure_counts,
      by = c("family", "corr_type", "corr_param", "n", "method", "term_type")
    )
}

# DGP Validation Checks -------------------------------------------------------

#' Validate Poisson DGP mean ranges
#'
#' Check that >95% of E[Y_i(t)] values fall in [0.5, 15].
#'
#' @param n_curves Number of curves to sample.
#' @param corr_type Correlation type.
#' @param corr_param Correlation parameter.
#' @param seed RNG seed.
#' @returns List with pass (logical), fraction_in_range, range.
validate_poisson_ranges <- function(
  n_curves = 1000,
  corr_type = "ar1",
  corr_param = 0.9,
  seed = 99
) {
  sim <- simulate_poisson_ff_linear(
    n = n_curves,
    corr_type = corr_type,
    corr_param = corr_param,
    seed = seed
  )
  # E[Y|X,z] = exp(eta_signal + sigma^2/2) = exp(truth$eta per subject)
  ey_mat <- exp(sim$truth$eta)
  ey_vals <- as.vector(ey_mat)
  frac_ok <- mean(ey_vals >= 0.5 & ey_vals <= 15)

  list(
    pass = frac_ok > 0.95,
    fraction_in_range = frac_ok,
    range = range(ey_vals),
    mean_ey = mean(ey_vals)
  )
}

#' Validate Binomial copula marginal probabilities
#'
#' Verify empirical P(Y=1) matches theoretical p_i(t) within MC tolerance.
#'
#' @param n Number of curves.
#' @param corr_type Correlation type.
#' @param corr_param Correlation parameter.
#' @param seed RNG seed.
#' @param n_mc Number of MC replicates to average over.
#' @returns List with pass (logical), max_deviation.
validate_binomial_copula <- function(
  n = 400,
  corr_type = "ar1",
  corr_param = 0.9,
  seed = 99,
  n_mc = 200
) {
  # Fix covariates and truth once, then re-draw Y multiple times
  # to verify E[Y|X] = p(X) under the copula construction.
  set.seed(seed)

  # Generate truth and covariates once (first draw sets the truth)
  sim0 <- simulate_binomial_ff_linear(
    n = n,
    corr_type = corr_type,
    corr_param = corr_param,
    seed = seed
  )
  p_mat <- stats::plogis(sim0$truth$eta) # n x nygrid (theoretical)

  # Now re-draw Y from the copula n_mc times, keeping X and truth fixed
  # We need to regenerate only the copula U and Y, not the covariates
  nygrid <- ncol(p_mat)
  err_struct <- sim0$err_struct
  R <- err_struct$corr_mat

  Y_sum <- matrix(0, nrow = n, ncol = nygrid)
  for (mc in seq_len(n_mc)) {
    set.seed(seed * 10000L + mc)
    Z <- mvtnorm::rmvnorm(n, mean = rep(0, nygrid), sigma = R)
    U <- pnorm(Z)
    Y <- matrix(as.integer(U < p_mat), nrow = n, ncol = nygrid)
    Y_sum <- Y_sum + Y
  }

  # Empirical rate across MC draws
  empirical_rate <- Y_sum / n_mc

  # Compare to theoretical probabilities
  deviation <- abs(empirical_rate - p_mat)

  # MC SE per cell: sqrt(p*(1-p)/n_mc)
  mc_se <- sqrt(p_mat * (1 - p_mat) / n_mc)

  # z-score of deviation
  z_dev <- deviation / pmax(mc_se, 1e-10)

  # Pass criteria: fraction of cells with |z| > 3 should be < 5%
  frac_large_z <- mean(z_dev > 3, na.rm = TRUE)

  list(
    pass = frac_large_z < 0.05,
    frac_large_z = frac_large_z,
    max_deviation = max(deviation),
    mean_deviation = mean(deviation),
    median_z = median(z_dev, na.rm = TRUE),
    max_z = max(z_dev, na.rm = TRUE)
  )
}

#' Validate Poisson collapsibility
#'
#' Verify coverage of conditional β ≈ coverage of marginal β.
#'
#' @param n_rep Number of reps.
#' @param n Sample size.
#' @param seed RNG seed.
#' @returns List with pass, coverage_diff.
validate_collapsibility <- function(n_rep = 20, n = 200, seed = 77) {
  # For Poisson log-link: conditional β = marginal β (collapsibility)
  # Coverage of both should be identical (same estimand)
  # This is inherent in the DGP construction — just verify coverage is similar
  # to what we'd expect under correctly specified model.

  coverages <- numeric(n_rep)
  for (rep in seq_len(n_rep)) {
    res <- tryCatch(
      {
        sim <- simulate_poisson_ff_linear(
          n = n,
          corr_type = "iid",
          seed = seed * 100 + rep
        )
        frml <- build_study1_formula(sim$s_grid)
        bs_yindex <- list(bs = "ps", k = STUDY1_K_YINDEX, m = c(2, 1))
        fit <- pffr(
          frml,
          yind = sim$t_grid,
          data = sim$data,
          family = poisson(),
          bs.yindex = bs_yindex
        )
        metrics <- compute_term_metrics(
          fit,
          sim$truth,
          "linear",
          alpha = 0.10,
          use_sandwich = FALSE,
          s_grid = sim$s_grid,
          t_grid = sim$t_grid,
          data = sim$data
        )
        if (!is.null(metrics)) metrics$coverage else NA_real_
      },
      error = function(e) NA_real_
    )
    coverages[rep] <- res
  }

  mean_cov <- mean(coverages, na.rm = TRUE)
  list(
    pass = abs(mean_cov - 0.90) < 0.15, # generous threshold for small n_rep
    mean_coverage = mean_cov,
    coverages = coverages
  )
}

# Phase 0: Known-Answer Validation --------------------------------------------

#' Phase 0: Fit OLS-equivalent pffr on trivial Gaussian DGP
#'
#' Verify ~90% coverage to confirm infrastructure works.
#'
#' @param n_rep Number of reps.
#' @param seed Base seed.
#' @returns List with pass and mean_coverage.
phase0_known_answer <- function(n_rep = 30, seed = 9999) {
  cat("Phase 0: Known-answer validation (Gaussian DGP, expect ~90% coverage)\n")

  coverages <- numeric(n_rep)
  for (rep in seq_len(n_rep)) {
    cov <- tryCatch(
      {
        # Use generate_benchmark_data from confint-benchmark.R for Gaussian DGP
        sim <- generate_benchmark_data(
          n = 100,
          nxgrid = 20,
          nygrid = 30,
          snr = 25,
          error_dist = "gaussian",
          corr_type = "iid",
          terms = c("ff", "linear"),
          wiggliness = 5,
          seed = seed + rep
        )

        frml <- build_pffr_formula(
          sim,
          sim$s_grid,
          k_smooth = 12,
          k_ff = c(12, 12)
        )
        bs_yindex <- list(bs = "ps", k = 12, m = c(2, 1))
        fit <- pffr(
          frml,
          yind = sim$t_grid,
          data = sim$data,
          bs.yindex = bs_yindex
        )

        metrics <- compute_term_metrics(
          fit,
          sim$truth,
          "linear",
          alpha = 0.10,
          use_sandwich = FALSE,
          s_grid = sim$s_grid,
          t_grid = sim$t_grid,
          data = sim$data
        )
        if (!is.null(metrics)) metrics$coverage else NA_real_
      },
      error = function(e) {
        message("Phase 0 error: ", e$message)
        NA_real_
      }
    )
    coverages[rep] <- cov
    cat(sprintf("\r  Rep %d/%d: coverage = %.2f", rep, n_rep, cov))
  }
  cat("\n")

  mean_cov <- mean(coverages, na.rm = TRUE)
  pass <- abs(mean_cov - 0.90) < 0.10
  cat(sprintf(
    "  Mean coverage: %.3f [%s]\n",
    mean_cov,
    if (pass) "PASS" else "FAIL"
  ))

  list(pass = pass, mean_coverage = mean_cov, coverages = coverages)
}

# Main Entry Point -------------------------------------------------------------

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  mode <- if (length(args) > 0) args[1] else "smoke"

  n_rep <- switch(
    mode,
    smoke = 2L,
    pilot = 10L,
    full = 150L,
    as.integer(mode) # allow numeric argument
  )
  if (is.na(n_rep)) n_rep <- 2L

  cat("Running Study 1 in", mode, "mode (", n_rep, "reps per cell)\n\n")

  # Phase 0: known-answer check
  if (mode != "full") {
    p0 <- phase0_known_answer(n_rep = 10)
    if (!p0$pass) {
      cat("WARNING: Phase 0 failed — infrastructure may have issues.\n")
    }
  }

  # Run study
  results <- run_study1(
    n_rep = n_rep,
    parallel = (n_rep > 5),
    n_workers = max(1L, min(6L, parallel::detectCores() - 1L))
  )

  # Summarize
  if (nrow(results) > 0) {
    summary_df <- summarize_study1(results)

    cat("\n========== STUDY 1 COVERAGE SUMMARY ==========\n")
    print(
      summary_df |>
        dplyr::select(
          family,
          corr_type,
          n,
          method,
          term_type,
          mean_coverage,
          se_coverage,
          mean_z_sd,
          n_reps
        ) |>
        dplyr::arrange(family, corr_type, n, method, term_type),
      n = 200
    )

    # Save summary
    output_dir <- "ci-benchmark/study1-nongaussian"
    write.csv(
      summary_df,
      file.path(output_dir, "summary.csv"),
      row.names = FALSE
    )
    cat("\nSummary saved to", file.path(output_dir, "summary.csv"), "\n")
  }

  cat("\nStudy 1 complete.\n")
}
