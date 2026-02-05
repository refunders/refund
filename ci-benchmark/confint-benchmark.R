# confint-benchmark.R
#
# Benchmark CI coverage/width for pffr coefficient estimates under different
# residual covariance structures and fitting methods (Round 2).
#
# Methods compared:
#   - pffr (default Bayesian Vc from mgcv)
#   - pffr_sandwich (cluster-robust sandwich SEs)
#   - pffr_hc (observation-level HC sandwich SEs -- included to show it doesn't
#     help under within-curve correlation)
#   - pffr_ar (AR(1) residuals via bam, gaussian + correlated only)
#   - pffr_gaulss (gaulss family, gaussian + heteroskedastic only)
#
# Metrics per term type (ff, concurrent, linear, smooth, intercept, E(Y)):
#   - coverage: fraction of grid points where CI contains truth
#   - coverage_high_var / coverage_low_var: region-specific coverage
#   - mean_width: average CI width across grid
#   - rmse, bias, mean_abs_error, mean_se, median_se
#   - z_mean, z_sd, z_kurtosis, z2_mean: z-score diagnostics
#   - fit_time_total: fitting time

# Setup -----------------------------------------------------------------------

library(tidyverse)

# Load refund (dev version if available)
if (file.exists("DESCRIPTION")) {
  devtools::load_all(".", quiet = TRUE)
} else {
  library(refund)
}

# Source utilities
source("ci-benchmark/benchmark-utils.R")

# DGP Settings ----------------------------------------------------------------

#' Create DGP settings grid
#'
#' Targeted (not full factorial) design with 4 arms:
#'   Arm 1 (correlation): vary rho with AR(1), no heteroskedasticity
#'   Arm 2 (fourier_pos): non-monotonic correlation at two strengths
#'   Arm 3 (heteroskedasticity): vary hetero under iid errors
#'   Arm 4 (interaction): correlation + heteroskedasticity at large n
#'
#' @param size One of "tiny", "small", "full".
#' @returns Tibble with DGP configurations.
make_dgp_settings <- function(size = c("tiny", "small", "full")) {
  size <- match.arg(size)

  # Fixed grid dimensions
  NXGRID <- 30L
  NYGRID <- 40L

  if (size == "tiny") {
    # Minimal for quick testing: 2 DGPs
    base <- tibble(
      n = 50L,
      nxgrid = NXGRID,
      nygrid = NYGRID,
      snr = 25,
      wiggliness = 5,
      error_dist = "gaussian",
      corr_type = c("iid", "ar1"),
      corr_param = c(NA_real_, 0.3),
      hetero_type = "none",
      hetero_param = NA_real_
    )
    return(
      base |>
        mutate(
          dgp_id = row_number(),
          terms = list(c("ff", "linear", "smooth", "concurrent"))
        )
    )
  }

  if (size == "small") {
    # Moderate for piloting: ~12 DGPs
    base <- bind_rows(
      # Correlation arm subset
      tidyr::crossing(
        n = c(50L, 200L),
        snr = 25,
        wiggliness = 5,
        corr_type = "ar1",
        corr_param = c(0, 0.3, 0.9),
        hetero_type = "none",
        hetero_param = NA_real_,
        error_dist = "gaussian"
      ),
      # Hetero arm subset
      tidyr::crossing(
        n = 200L,
        snr = 25,
        wiggliness = 5,
        corr_type = "iid",
        corr_param = NA_real_,
        hetero_type = c("bump"),
        hetero_param = 3.0,
        error_dist = "gaussian"
      )
    ) |>
      mutate(nxgrid = NXGRID, nygrid = NYGRID)

    return(
      base |>
        mutate(
          dgp_id = row_number(),
          terms = list(c("ff", "linear", "smooth", "concurrent"))
        )
    )
  }

  # Full design: ~72 DGPs

  # Arm 1: Correlation (48 DGPs)
  # vary n, snr, wiggliness, rho, error_dist; hetero=none
  # Map rho=0 to corr_type="iid", rho>0 to "ar1"
  arm1_raw <- tidyr::crossing(
    n = c(50L, 200L),
    snr = c(5, 25),
    wiggliness = c(0.01, 5),
    rho = c(0, 0.3, 0.9),
    error_dist = c("gaussian", "t6")
  ) |>
    mutate(
      corr_type = ifelse(rho == 0, "iid", "ar1"),
      corr_param = ifelse(rho == 0, NA_real_, rho),
      hetero_type = "none",
      hetero_param = NA_real_
    ) |>
    select(-rho)

  # Arm 2: Non-monotonic correlation (4 DGPs)
  # fourier_pos at two period values, n={50, 200}
  arm2 <- tidyr::crossing(
    n = c(50L, 200L),
    snr = 25,
    wiggliness = 5,
    corr_type = "fourier_pos",
    corr_param = c(0.3, 0.6),
    error_dist = "gaussian",
    hetero_type = "none",
    hetero_param = NA_real_
  )

  # Arm 3: Heteroskedasticity (16 DGPs)
  # vary n, snr, wiggliness, hetero level; iid errors, gaussian
  arm3 <- tidyr::crossing(
    n = c(50L, 200L),
    snr = c(5, 25),
    wiggliness = c(0.01, 5),
    corr_type = "iid",
    corr_param = NA_real_,
    error_dist = "gaussian",
    hetero_type = c("bump"),
    hetero_param = c(1.0, 3.0)
  )

  # Arm 4: Interaction (4 DGPs)
  # correlation + heteroskedasticity at n=200, high SNR
  arm4 <- tidyr::crossing(
    n = 200L,
    snr = 25,
    wiggliness = 5,
    corr_type = "ar1",
    corr_param = c(0.3, 0.9),
    error_dist = "gaussian",
    hetero_type = "bump",
    hetero_param = c(1.0, 3.0)
  )

  base <- bind_rows(arm1_raw, arm2, arm3, arm4) |>
    distinct() |>
    mutate(nxgrid = NXGRID, nygrid = NYGRID)

  base |>
    mutate(
      dgp_id = row_number(),
      terms = list(c("ff", "linear", "smooth", "concurrent"))
    )
}

#' Determine which methods are applicable for a DGP
#'
#' Methods:
#'   pffr: default Bayesian Vc from mgcv (no sandwich)
#'   pffr_sandwich: cluster-robust sandwich (sandwich = "cluster")
#'   pffr_hc: observation-level HC sandwich (sandwich = "hc")
#'   pffr_ar: bam with estimated AR(1) rho (gaussian + correlated only)
#'   pffr_gaulss: gaulss family (gaussian + hetero != "none" only)
#'
#' @param corr_type Correlation type.
#' @param hetero_type Heteroskedasticity type.
#' @param error_dist Error distribution.
#' @returns Character vector of method names.
get_applicable_methods <- function(corr_type, hetero_type, error_dist) {
  # Be robust to factors / accidental whitespace / case
  corr_type <- tolower(trimws(as.character(corr_type %||% "")))
  hetero_type <- tolower(trimws(as.character(hetero_type %||% "")))
  error_dist <- tolower(trimws(as.character(error_dist %||% "")))

  # Core trio: always run default, cluster sandwich, and HC sandwich
  # (HC is included to demonstrate it does NOT help under correlation)
  methods <- c("pffr", "pffr_sandwich", "pffr_hc")

  # pffr_ar: only meaningful under correlated errors, but can be run under any
  # error distribution (it is still a valid fitting strategy under
  # misspecification, e.g. t errors).
  if (corr_type %in% c("ar1", "fourier_pos")) {
    methods <- c(methods, "pffr_ar")
  }

  # pffr_gaulss: gaussian + heteroskedastic only
  if (hetero_type != "none" && error_dist == "gaussian") {
    methods <- c(methods, "pffr_gaulss")
  }

  methods
}

# Data Generation -------------------------------------------------------------

#' Generate benchmark data
#'
#' Creates simulated data with known ground truth for CI evaluation.
#'
#' @param n Number of observations.
#' @param nxgrid Grid size for functional covariate.
#' @param nygrid Grid size for response.
#' @param snr Signal-to-noise ratio.
#' @param error_dist Error distribution ("gaussian" or "t6").
#' @param corr_type Correlation type ("iid", "ar1").
#' @param corr_param Correlation parameter.
#' @param hetero_type Heteroskedasticity type.
#' @param hetero_param Heteroskedasticity parameter.
#' @param terms Character vector of terms to include.
#' @param wiggliness Wiggliness parameter for random truth generation (default 3).
#' @param seed RNG seed.
#' @returns List with data, grids, truth, and error structure.
generate_benchmark_data <- function(
  n,
  nxgrid,
  nygrid,
  snr,
  error_dist = "gaussian",
  corr_type = "iid",
  corr_param = NA,
  hetero_type = "none",
  hetero_param = NA,
  terms = c("ff", "linear", "smooth"),
  wiggliness = 3,
  seed = NULL,
  balance_terms = TRUE
) {
  if (!is.null(seed)) set.seed(seed)

  # Build formula from terms
  formula_parts <- c(
    if ("ff" %in% terms) "ff(X1)" else NULL,
    if ("linear" %in% terms) "zlin" else NULL,
    if ("smooth" %in% terms) "s(zsmoo)" else NULL,
    if ("concurrent" %in% terms) "Xconc" else NULL
  )
  formula <- as.formula(paste("Y ~", paste(formula_parts, collapse = " + ")))

  # Use "random" effects with specified wiggliness for reproducible random truth
  # Concurrent terms need list syntax: list(type = "concurrent", effect = "random")
  effects <- list(
    X1 = "random",
    zlin = "random",
    zsmoo = "random",
    Xconc = list(type = "concurrent", effect = "random")
  )

  # Generate base data with pffr_simulate (iid gaussian errors, SNR applied)
  # k_truth = 6 for all terms to ensure k_model/k_truth ratio of 2x per margin
  # (model uses k=12 for smooth/ff), reducing smoothing bias
  dat <- pffr_simulate(
    formula,
    n = n,
    nxgrid = nxgrid,
    nygrid = nygrid,
    SNR = snr,
    effects = effects,
    intercept = "random",
    wiggliness = wiggliness,
    k_truth = list(
      ff_s = 6,
      ff_t = 6,
      smooth_z = 6,
      smooth_t = 6,
      linear = 6,
      concurrent = 6,
      intercept = 6
    ),
    seed = seed
  )

  s_grid <- attr(dat, "xindex")
  t_grid <- attr(dat, "yindex")
  truth <- attr(dat, "truth")

  # Balance term contributions if requested (done before generating errors so the
  # requested SNR applies to the final linear predictor).
  balance_term_contributions <- function(
    truth,
    term_labels = names(truth$etaTerms),
    target = c("median", "mean"),
    measure = c("sd", "rms")
  ) {
    target <- match.arg(target)
    measure <- match.arg(measure)

    term_magnitude <- function(x) {
      v <- as.vector(x)
      v <- v[is.finite(v)]
      if (!length(v)) return(NA_real_)
      if (measure == "sd") return(stats::sd(v))
      sqrt(mean(v^2))
    }

    eta_terms <- truth$etaTerms
    term_labels <- intersect(term_labels, names(eta_terms))
    term_labels <- term_labels[term_labels != "epsilon"]

    mags <- vapply(
      term_labels,
      function(lbl) term_magnitude(eta_terms[[lbl]]),
      numeric(1)
    )
    mags_ok <- mags[is.finite(mags) & mags > 0]
    if (!length(mags_ok)) return(truth)

    target_mag <- if (target == "median") stats::median(mags_ok) else
      mean(mags_ok)
    scales <- target_mag / mags
    scales[!is.finite(scales)] <- 1

    for (lbl in term_labels) {
      sc <- unname(scales[[lbl]])
      if (!is.finite(sc) || sc == 1) next

      eta_terms[[lbl]] <- eta_terms[[lbl]] * sc

      b <- truth$beta[[lbl]]
      if (is.function(b)) {
        truth$beta[[lbl]] <- local({
          b0 <- b
          sc0 <- sc
          function(...) sc0 * b0(...)
        })
      } else {
        truth$beta[[lbl]] <- b * sc
      }
    }

    truth$etaTerms <- eta_terms
    truth$eta <- Reduce(`+`, eta_terms[term_labels])
    truth$term_scales <- as.list(scales)
    truth
  }

  if (isTRUE(balance_terms)) {
    truth <- balance_term_contributions(
      truth,
      term_labels = c("intercept", "ff(X1)", "zlin", "s(zsmoo)", "Xconc"),
      target = "median",
      measure = "sd"
    )
  }

  # Get the signal (eta) before adding errors
  eta <- truth$eta

  # Generate errors with the requested structure, then scale to achieve target SNR.
  err_struct <- make_error_structure(
    t_grid = t_grid,
    corr_type = corr_type,
    corr_param = if (is.na(corr_param)) NULL else corr_param,
    hetero_type = hetero_type,
    hetero_param = if (is.na(hetero_param)) NULL else hetero_param
  )

  eps_new <- sample_errors(n, err_struct$cov_mat, error_dist)

  scale_fac <- compute_snr_scale_factor(eps_new, eta, snr)
  eps_new <- eps_new * scale_fac
  err_struct$cov_mat <- err_struct$cov_mat * scale_fac^2

  dat$Y <- I(eta + eps_new)
  truth$epsilon <- eps_new
  attr(dat, "truth") <- truth

  # Note: Do NOT center covariates after data generation - this would break
  # the relationship between the stored truth and actual effect values.
  # pffr handles identifiability constraints internally.

  list(
    data = dat,
    s_grid = s_grid,
    t_grid = t_grid,
    truth = truth,
    err_struct = err_struct,
    terms = terms
  )
}

# Model Fitting ---------------------------------------------------------------

#' Build pffr formula from simulation
#'
#' @param sim Simulation result from generate_benchmark_data().
#' @param s_grid s evaluation grid (for ff xind).
#' @param k_smooth Basis dimension for smooth terms (default 12).
#' @param k_ff Basis dimensions for ff term as c(k_s, k_t) (default c(12, 12)).
#' @returns Formula object with environment containing s_grid.
#' @details Default basis dimensions are set larger than truth generation
#'   defaults (k=8 for all terms) to avoid smoothing bias.
build_pffr_formula <- function(sim, s_grid, k_smooth = 12, k_ff = c(12, 12)) {
  terms <- sim$terms

  # xind for ff is passed via s_grid variable
  # Use larger k for ff term to reduce smoothing bias
  formula_parts <- c(
    if ("ff" %in% terms) {
      sprintf(
        "ff(X1, xind = s_grid, splinepars = list(bs = 'ps', k = c(%d, %d)))",
        k_ff[1],
        k_ff[2]
      )
    } else {
      NULL
    },
    if ("linear" %in% terms) "zlin" else NULL,
    if ("smooth" %in% terms) sprintf("s(zsmoo, k = %d)", k_smooth) else NULL,
    if ("concurrent" %in% terms) "Xconc" else NULL
  )

  frml <- as.formula(paste("Y ~", paste(formula_parts, collapse = " + ")))

  # Set formula environment to include s_grid
  e <- new.env(parent = globalenv())
  e$s_grid <- s_grid
  environment(frml) <- e

  frml
}

#' Fit all methods for one dataset
#'
#' @param sim Simulation result.
#' @param methods Character vector of methods to fit.
#' @param k_yindex Basis dimension for varying coefficient terms (default 12).
#' @returns List with fits and timings.
#' @details Default k_yindex must be larger than truth generation default
#'   (k=8 for all terms) to avoid smoothing bias.
fit_all_methods <- function(sim, methods, k_yindex = 12) {
  dat <- sim$data
  t_grid <- sim$t_grid
  s_grid <- sim$s_grid

  formula <- build_pffr_formula(sim, s_grid, k_smooth = 12, k_ff = c(12, 12))

  # Use larger basis for varying coefficient terms to reduce smoothing bias
  bs_yindex <- list(bs = "ps", k = k_yindex, m = c(2, 1))

  fits <- list()
  timings <- list()

  # Always fit pffr first (pilot for GLS/AR)
  t0 <- Sys.time()
  fits$pffr <- pffr(formula, yind = t_grid, data = dat, bs.yindex = bs_yindex)
  pffr_time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  timings$pffr <- pffr_time

  # pffr_sandwich: same fit, cluster-robust SE extraction
  if ("pffr_sandwich" %in% methods) {
    fits$pffr_sandwich <- fits$pffr
    timings$pffr_sandwich <- pffr_time
  }

  # pffr_hc: same fit, observation-level HC sandwich SE extraction
  if ("pffr_hc" %in% methods) {
    fits$pffr_hc <- fits$pffr
    timings$pffr_hc <- pffr_time
  }

  # pffr_gaulss: fit with gaulss family to model heteroskedasticity
  if ("pffr_gaulss" %in% methods) {
    t0 <- Sys.time()
    fits$pffr_gaulss <- pffr(
      formula,
      yind = t_grid,
      data = dat,
      family = mgcv::gaulss(),
      bs.yindex = bs_yindex
    )
    gaulss_time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    timings$pffr_gaulss <- gaulss_time
  }

  # pffr_ar: requires pilot fit for rho estimation
  if ("pffr_ar" %in% methods) {
    t0 <- Sys.time()
    rho <- estimate_rho_from_fit(fits$pffr, dat)
    fits$pffr_ar <- pffr(
      formula,
      yind = t_grid,
      data = dat,
      algorithm = "bam",
      rho = rho,
      bs.yindex = bs_yindex
    )
    ar_time <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    timings$pffr_ar <- pffr_time + ar_time
    timings$pffr_ar_only <- ar_time
  }

  list(fits = fits, timings = timings, pffr_pilot_time = pffr_time)
}

# Metric Computation ----------------------------------------------------------

#' Find term index in coef.pffr output by term type
#'
#' @param sm_names Names of smterms from coef().
#' @param term_type One of "ff", "linear", "smooth", "intercept", "concurrent".
#' @returns Integer index or NULL.
find_term_index <- function(sm_names, term_type) {
  for (i in seq_along(sm_names)) {
    nm <- tolower(sm_names[i])
    if (term_type == "ff" && grepl("ff\\(x1\\)", nm)) return(i)
    if (term_type == "linear" && grepl("zlin", nm)) return(i)
    if (term_type == "smooth" && grepl("zsmoo", nm)) return(i)
    if (term_type == "intercept" && grepl("^intercept\\(", nm)) return(i)
    if (term_type == "concurrent" && grepl("xconc", nm)) return(i)
  }
  NULL
}

ff_smooth_index <- function(fit) {
  if (is.null(fit) || is.null(fit$pffr)) return(NULL)
  if (is.null(fit$pffr$where$ff) || length(fit$pffr$where$ff) < 1) return(NULL)
  if (is.null(fit$pffr$labelmap) || is.null(fit$smooth)) return(NULL)

  ff_label_idx <- fit$pffr$where$ff[1]
  ff_smooth_name <- unname(fit$pffr$labelmap[ff_label_idx])
  idx <- match(ff_smooth_name, names(fit$smooth))
  if (is.na(idx)) return(NULL)
  idx
}

should_center_ff_truth <- function(fit) {
  idx <- ff_smooth_index(fit)
  if (is.null(idx)) return(FALSE)
  n_cons <- attr(fit$smooth[[idx]], "nCons")
  isTRUE(!is.null(n_cons) && n_cons > 0)
}

#' Evaluate truth on coefficient grid
#'
#' @param truth Truth list from simulation (with beta element).
#' @param term_type Term type.
#' @param term_info Coefficient term info from coef.pffr.
#' @param s_grid Original s grid from simulation.
#' @param t_grid Original t grid from simulation.
#' @param data Data frame used for fitting (needed for smooth term centering).
#' @param fit Optional pffr fit used to decide whether to apply centering
#'   constraints when evaluating truth.
#' @param center_ff One of "match_fit", "always", "never". If "match_fit",
#'   center ff truth only when the fitted ff term has constraints.
#' @returns Numeric vector of truth values matching est length.
evaluate_truth_on_grid <- function(
  truth,
  term_type,
  term_info,
  s_grid,
  t_grid,
  data = NULL,
  fit = NULL,
  center_ff = c("match_fit", "always", "never")
) {
  beta <- truth$beta
  if (is.null(beta)) return(NULL)

  center_ff <- match.arg(center_ff)

  # Get evaluation grid from term_info
  x_eval <- term_info$x
  y_eval <- term_info$y

  if (term_type == "ff") {
    # beta(s, t) - stored as matrix ns x nt
    beta_mat <- beta[["ff(X1)"]]
    if (is.null(beta_mat)) return(NULL)

    # Interpolate to evaluation grid
    # beta_mat is on s_grid (rows) x t_grid (cols)
    # Need to interpolate to x_eval x y_eval
    if (!is.null(x_eval) && !is.null(y_eval)) {
      # Use bilinear interpolation
      interp_beta <- function(s_new, t_new) {
        # For each (s_new, t_new) point, interpolate from beta_mat
        outer(
          s_new,
          t_new,
          Vectorize(function(s, t) {
            # Find surrounding indices in original grids
            si <- findInterval(s, s_grid, all.inside = TRUE)
            ti <- findInterval(t, t_grid, all.inside = TRUE)
            # Bilinear weights
            s_frac <- (s - s_grid[si]) / (s_grid[si + 1] - s_grid[si])
            t_frac <- (t - t_grid[ti]) / (t_grid[ti + 1] - t_grid[ti])
            # Interpolate
            (1 - s_frac) *
              (1 - t_frac) *
              beta_mat[si, ti] +
              s_frac * (1 - t_frac) * beta_mat[si + 1, ti] +
              (1 - s_frac) * t_frac * beta_mat[si, ti + 1] +
              s_frac * t_frac * beta_mat[si + 1, ti + 1]
          })
        )
      }
      beta_eval <- interp_beta(x_eval, y_eval)
    } else {
      beta_eval <- beta_mat
    }

    do_center <- switch(
      center_ff,
      always = TRUE,
      never = FALSE,
      match_fit = should_center_ff_truth(fit)
    )
    if (do_center) {
      # Center to match pffr constraints: sum_s w(s) * beta(s, t) = 0 for all t
      w_s <- ff_weights(x_eval %||% s_grid, method = "simpson")
      beta_eval <- center_beta_ff(beta_eval, x_eval %||% s_grid, w_s)
    }

    return(as.vector(beta_eval))
  }

  if (term_type == "linear") {
    # beta(t) - stored as vector of length nt
    beta_t <- beta$zlin
    if (is.null(beta_t)) return(NULL)

    # Interpolate to x_eval (which is the t grid for 1D terms)
    if (!is.null(x_eval) && length(x_eval) != length(beta_t)) {
      beta_eval <- stats::approx(t_grid, beta_t, xout = x_eval, rule = 2)$y
    } else {
      beta_eval <- beta_t
    }

    return(beta_eval)
  }

  if (term_type == "smooth") {
    # f(z, t) - stored as function or matrix
    f_zt <- beta[["s(zsmoo)"]]
    if (is.null(f_zt)) return(NULL)

    if (is.function(f_zt)) {
      # Evaluate function on grid
      if (!is.null(x_eval) && !is.null(y_eval)) {
        f_eval <- f_zt(x_eval, y_eval)
      } else {
        return(NULL)
      }
    } else if (is.matrix(f_zt)) {
      f_eval <- f_zt
    } else {
      return(NULL)
    }

    # Center per t to match pffr constraints
    # Must use actual observed z values, not the evaluation grid
    if (is.matrix(f_eval) && is.function(f_zt)) {
      z_obs <- if (!is.null(data) && !is.null(data$zsmoo)) {
        as.vector(data$zsmoo)
      } else {
        warning(
          "No z observations available for smooth centering, using eval grid"
        )
        x_eval
      }
      f_on_obs <- f_zt(z_obs, y_eval)
      f_eval <- sweep(f_eval, 2, colMeans(f_on_obs), "-")
    }

    return(as.vector(f_eval))
  }

  if (term_type == "intercept") {
    # mu(t) - stored as 1D vector over t
    beta_t <- beta$intercept
    if (is.null(beta_t)) return(NULL)

    # Interpolate to evaluation grid if needed
    if (!is.null(x_eval) && length(x_eval) != length(beta_t)) {
      beta_eval <- stats::approx(t_grid, beta_t, xout = x_eval, rule = 2)$y
    } else {
      beta_eval <- beta_t
    }
    return(beta_eval)
  }

  if (term_type == "concurrent") {
    # beta(t) for concurrent term - stored as 1D vector
    beta_t <- beta$Xconc
    if (is.null(beta_t)) return(NULL)

    # Interpolate to evaluation grid if needed
    if (!is.null(x_eval) && length(x_eval) != length(beta_t)) {
      beta_eval <- stats::approx(t_grid, beta_t, xout = x_eval, rule = 2)$y
    } else {
      beta_eval <- beta_t
    }
    return(beta_eval)
  }

  NULL
}

#' Extract coefficient estimates + CIs on coef.pffr grid
#'
#' Returns a tidy data frame on the exact grid used by `coef.pffr`, including
#' point estimates, SEs, CIs, truth, and coverage indicators. This is useful
#' for interactive exploration/plotting while using the same code-path as the
#' benchmark metric computation.
#'
#' @param fit A pffr fit.
#' @param truth Truth list from simulation.
#' @param term_type One of "ff", "linear", "smooth", "intercept", "concurrent".
#' @param alpha Significance level for CI.
#' @param use_sandwich Sandwich type: FALSE, TRUE, "cluster", or "hc".
#' @param s_grid s evaluation points (for ff).
#' @param t_grid t evaluation points.
#' @param data Data frame used for fitting (needed for smooth term centering).
#' @param coefs Optional pre-computed coefficients from coef.pffr.
#' @param center_ff One of "match_fit", "always", "never".
#' @returns Tibble with columns x, y (if 2D), estimate, se, lower, upper, truth,
#'   covered, error, z, term_type.
extract_term_ci_df <- function(
  fit,
  truth,
  term_type,
  alpha = 0.10,
  use_sandwich = FALSE,
  s_grid = NULL,
  t_grid = NULL,
  data = NULL,
  coefs = NULL,
  center_ff = c("match_fit", "always", "never")
) {
  use_sandwich <- normalize_sandwich_type(use_sandwich)
  center_ff <- match.arg(center_ff)

  if (is.null(coefs)) {
    coefs <- tryCatch(
      coef(
        fit,
        sandwich = use_sandwich,
        seWithMean = FALSE,
        n1 = 50,
        n2 = 25,
        n3 = 15
      ),
      error = function(e) NULL
    )
  }
  if (is.null(coefs) || is.null(coefs$smterms)) return(NULL)

  sm_names <- names(coefs$smterms)
  term_idx <- find_term_index(sm_names, term_type)
  if (is.null(term_idx)) return(NULL)

  term_info <- coefs$smterms[[term_idx]]
  if (is.null(term_info) || is.null(term_info$coef)) return(NULL)

  est <- term_info$coef$value
  se <- term_info$coef$se
  if (is.null(est) || is.null(se) || !length(se)) return(NULL)

  truth_vals <- evaluate_truth_on_grid(
    truth,
    term_type,
    term_info,
    s_grid,
    t_grid,
    data = data,
    fit = fit,
    center_ff = center_ff
  )
  if (is.null(truth_vals)) return(NULL)

  if (length(truth_vals) != length(est)) {
    warning(sprintf(
      "Length mismatch for %s: est=%d, truth=%d",
      term_type,
      length(est),
      length(truth_vals)
    ))
    return(NULL)
  }

  z_crit <- qnorm(1 - alpha / 2)
  lower <- est - z_crit * se
  upper <- est + z_crit * se
  err <- est - truth_vals
  z_score <- ifelse(is.finite(se) & se > 0, err / se, NA_real_)
  covered <- (truth_vals >= lower) & (truth_vals <= upper)

  is_2d <- !is.null(term_info$y)
  if (is_2d) {
    x_eval <- term_info$x
    y_eval <- term_info$y
    if (length(est) != length(x_eval) * length(y_eval)) return(NULL)
    grid_df <- expand.grid(x = x_eval, y = y_eval)
    return(tibble(
      term_type = term_type,
      x = grid_df$x,
      y = grid_df$y,
      estimate = est,
      se = se,
      lower = lower,
      upper = upper,
      truth = truth_vals,
      covered = covered,
      error = err,
      z = z_score
    ))
  }

  x_eval <- term_info$x
  if (length(est) != length(x_eval)) return(NULL)

  tibble(
    term_type = term_type,
    x = x_eval,
    estimate = est,
    se = se,
    lower = lower,
    upper = upper,
    truth = truth_vals,
    covered = covered,
    error = err,
    z = z_score
  )
}

#' Compute metrics for one term
#'
#' @param fit A pffr fit.
#' @param truth Truth list from simulation.
#' @param term_type One of "ff", "linear", "smooth", "intercept", "concurrent".
#' @param alpha Significance level for CI.
#' @param use_sandwich Sandwich type: FALSE, "cluster", or "hc".
#' @param s_grid s evaluation points (for ff).
#' @param t_grid t evaluation points.
#' @param data Data frame used for fitting (needed for smooth term centering).
#' @param coefs Pre-computed coefficients from coef.pffr.
#' @param err_struct Error structure from simulation (for region-specific coverage).
#' @returns Tibble with one row of metrics, or NULL if term not found.
compute_term_metrics <- function(
  fit,
  truth,
  term_type,
  alpha = 0.10,
  use_sandwich = FALSE,
  s_grid = NULL,
  t_grid = NULL,
  data = NULL,
  coefs = NULL,
  err_struct = NULL
) {
  use_sandwich <- normalize_sandwich_type(use_sandwich)

  df <- extract_term_ci_df(
    fit,
    truth,
    term_type,
    alpha = alpha,
    use_sandwich = use_sandwich,
    s_grid = s_grid,
    t_grid = t_grid,
    data = data,
    coefs = coefs
  )
  if (is.null(df) || nrow(df) == 0) return(NULL)

  z_finite <- df$z[is.finite(df$z)]

  z_kurtosis <- if (length(z_finite) > 3) {
    m4 <- mean((z_finite - mean(z_finite))^4)
    s2 <- stats::var(z_finite)
    if (s2 > 0) m4 / s2^2 - 3 else NA_real_
  } else {
    NA_real_
  }

  coverage_high_var <- NA_real_
  coverage_low_var <- NA_real_

  if (!is.null(err_struct) && err_struct$hetero_type != "none") {
    sigma_t <- err_struct$sigma_t
    sigma_median <- stats::median(sigma_t)
    if (!is.null(df$y)) {
      t_to_sigma <- stats::approx(t_grid, sigma_t, xout = df$y, rule = 2)$y
      hi <- t_to_sigma >= sigma_median
      lo <- t_to_sigma < sigma_median
      if (any(hi, na.rm = TRUE)) {
        coverage_high_var <- mean(df$covered[hi], na.rm = TRUE)
      }
      if (any(lo, na.rm = TRUE)) {
        coverage_low_var <- mean(df$covered[lo], na.rm = TRUE)
      }
    } else if (!is.null(df$x)) {
      t_to_sigma <- stats::approx(t_grid, sigma_t, xout = df$x, rule = 2)$y
      hi <- t_to_sigma >= sigma_median
      lo <- t_to_sigma < sigma_median
      if (any(hi, na.rm = TRUE)) {
        coverage_high_var <- mean(df$covered[hi], na.rm = TRUE)
      }
      if (any(lo, na.rm = TRUE)) {
        coverage_low_var <- mean(df$covered[lo], na.rm = TRUE)
      }
    }
  }

  tibble(
    term_type = term_type,
    coverage = mean(df$covered, na.rm = TRUE),
    coverage_high_var = coverage_high_var,
    coverage_low_var = coverage_low_var,
    mean_width = mean(df$upper - df$lower, na.rm = TRUE),
    rmse = sqrt(mean(df$error^2, na.rm = TRUE)),
    bias = mean(df$error, na.rm = TRUE),
    mean_abs_error = mean(abs(df$error), na.rm = TRUE),
    mean_se = mean(df$se, na.rm = TRUE),
    median_se = stats::median(df$se, na.rm = TRUE),
    z_mean = mean(z_finite),
    z_sd = stats::sd(z_finite),
    z_kurtosis = z_kurtosis,
    z2_mean = mean(z_finite^2),
    n_grid = sum(!is.na(df$covered))
  )
}

#' Compute E(Y) metrics (fitted values vs true eta)
#'
#' @param fit A pffr fit.
#' @param truth Truth list from simulation.
#' @param alpha Significance level for CI.
#' @param use_sandwich Sandwich type: FALSE, "cluster", or "hc". When specified,
#'   applies sandwich correction to the fit's covariance matrices before
#'   computing prediction SEs.
#' @param err_struct Error structure from simulation (for region-specific coverage).
#' @param t_grid t evaluation points (needed to map sigma_t to columns).
#' @returns Tibble with one row of E(Y) metrics.
compute_ey_metrics <- function(
  fit,
  truth,
  alpha = 0.10,
  use_sandwich = FALSE,
  err_struct = NULL,
  t_grid = NULL
) {
  use_sandwich <- normalize_sandwich_type(use_sandwich)

  # Apply sandwich correction if requested
  pred_fit <- fit
  if (!isFALSE(use_sandwich)) {
    algorithm <- fit$pffr$algorithm %||% "gam"
    sandwich_fn <- if (exists("apply_sandwich_correction", mode = "function")) {
      apply_sandwich_correction
    } else {
      refund:::apply_sandwich_correction
    }
    pred_fit <- sandwich_fn(fit, algorithm, type = use_sandwich)
  }

  # Get fitted values with SE (suppress gaulss reshape warnings)
  pred <- tryCatch(
    suppressWarnings(predict(pred_fit, se.fit = TRUE, type = "link")),
    error = function(e) NULL
  )
  if (is.null(pred)) return(NULL)

  # Truth is the full eta matrix (n x nyindex)
  eta_true <- truth$eta
  if (is.null(eta_true)) return(NULL)

  n <- nrow(eta_true)
  nygrid <- ncol(eta_true)

  # predict.pffr returns n×nygrid matrices (even for gaulss, which internally
  # truncates to the mu linear predictor). Use as-is if already a matrix with
  # correct dimensions, otherwise reshape.
  if (is.matrix(pred$fit) && nrow(pred$fit) == n && ncol(pred$fit) == nygrid) {
    fit_mat <- pred$fit
    se_mat <- pred$se.fit
  } else {
    # Fallback: vector output or unexpected shape — take first n*nygrid elements
    fit_vec <- as.vector(pred$fit)[seq_len(n * nygrid)]
    se_vec <- as.vector(pred$se.fit)[seq_len(n * nygrid)]
    fit_mat <- matrix(fit_vec, nrow = n, ncol = nygrid)
    se_mat <- matrix(se_vec, nrow = n, ncol = nygrid)
  }

  # Compute pointwise CIs
  z <- qnorm(1 - alpha / 2)
  lower <- fit_mat - z * se_mat
  upper <- fit_mat + z * se_mat

  # Coverage, width, RMSE
  covered <- (eta_true >= lower) & (eta_true <= upper)

  err <- fit_mat - eta_true
  z_score <- ifelse(is.finite(se_mat) & se_mat > 0, err / se_mat, NA_real_)
  z_finite <- as.vector(z_score[is.finite(z_score)])

  z_kurtosis <- if (length(z_finite) > 3) {
    m4 <- mean((z_finite - mean(z_finite))^4)
    s2 <- stats::var(z_finite)
    if (s2 > 0) m4 / s2^2 - 3 else NA_real_
  } else {
    NA_real_
  }

  # Region-specific coverage: split columns by sigma_t pattern
  coverage_high_var <- NA_real_
  coverage_low_var <- NA_real_

  if (
    !is.null(err_struct) && err_struct$hetero_type != "none" && !is.null(t_grid)
  ) {
    sigma_t <- err_struct$sigma_t
    sigma_median <- stats::median(sigma_t)
    # Map sigma_t (on original t_grid) to the nygrid columns of the covered matrix
    # The response grid used by pffr may differ from t_grid, so use the fit's yindex
    y_eval <- fit$pffr$yind %||% t_grid
    if (length(y_eval) == nygrid) {
      t_to_sigma <- stats::approx(t_grid, sigma_t, xout = y_eval, rule = 2)$y
      high_cols <- which(t_to_sigma >= sigma_median)
      low_cols <- which(t_to_sigma < sigma_median)
      if (length(high_cols) > 0) {
        coverage_high_var <- mean(covered[, high_cols], na.rm = TRUE)
      }
      if (length(low_cols) > 0) {
        coverage_low_var <- mean(covered[, low_cols], na.rm = TRUE)
      }
    }
  }

  tibble(
    term_type = "E(Y)",
    coverage = mean(covered, na.rm = TRUE),
    coverage_high_var = coverage_high_var,
    coverage_low_var = coverage_low_var,
    mean_width = mean(2 * z * se_mat, na.rm = TRUE),
    rmse = sqrt(mean(err^2, na.rm = TRUE)),
    bias = mean(err, na.rm = TRUE),
    mean_abs_error = mean(abs(err), na.rm = TRUE),
    mean_se = mean(se_mat, na.rm = TRUE),
    median_se = stats::median(se_mat, na.rm = TRUE),
    z_mean = mean(z_finite),
    z_sd = stats::sd(z_finite),
    z_kurtosis = z_kurtosis,
    z2_mean = mean(z_finite^2),
    n_grid = sum(!is.na(covered))
  )
}

#' Compute all metrics for one fit
#'
#' @param fit A pffr fit.
#' @param sim Simulation result.
#' @param method Method name.
#' @param use_sandwich Sandwich type: FALSE, "cluster", or "hc".
#' @param alpha Significance level.
#' @param err_struct Error structure from simulation (for region-specific coverage).
#' @returns Tibble with one row per term.
compute_all_metrics <- function(
  fit,
  sim,
  method,
  use_sandwich = FALSE,
  alpha = 0.10,
  err_struct = NULL
) {
  use_sandwich <- normalize_sandwich_type(use_sandwich)
  terms <- sim$terms
  truth <- sim$truth

  # Extract coefficients once per fit (coef.pffr can be expensive)
  coefs <- tryCatch(
    coef(
      fit,
      sandwich = use_sandwich,
      seWithMean = FALSE,
      n1 = 50,
      n2 = 25,
      n3 = 15
    ),
    error = function(e) {
      message("coef error: ", e$message)
      NULL
    }
  )
  if (is.null(coefs) || is.null(coefs$smterms)) return(tibble())

  results <- map_dfr(terms, function(term_type) {
    metrics <- compute_term_metrics(
      fit,
      truth,
      term_type,
      alpha = alpha,
      use_sandwich = use_sandwich,
      s_grid = sim$s_grid,
      t_grid = sim$t_grid,
      data = sim$data,
      coefs = coefs,
      err_struct = err_struct
    )
    if (!is.null(metrics)) {
      metrics$method <- method
    }
    metrics
  })

  # Add intercept metrics (always present in pffr models)
  intercept_metrics <- compute_term_metrics(
    fit,
    truth,
    "intercept",
    alpha = alpha,
    use_sandwich = use_sandwich,
    s_grid = sim$s_grid,
    t_grid = sim$t_grid,
    data = sim$data,
    coefs = coefs,
    err_struct = err_struct
  )
  if (!is.null(intercept_metrics)) {
    intercept_metrics$method <- method
    results <- bind_rows(results, intercept_metrics)
  }

  # Add E(Y) metrics (fitted values vs true eta)
  ey_metrics <- compute_ey_metrics(
    fit,
    truth,
    alpha,
    use_sandwich = use_sandwich,
    err_struct = err_struct,
    t_grid = sim$t_grid
  )
  if (!is.null(ey_metrics)) {
    ey_metrics$method <- method
    results <- bind_rows(results, ey_metrics)
  }

  results
}

# Benchmark Runner ------------------------------------------------------------

#' Create a standardized failure metrics row
#'
#' @param row One-row tibble/list of benchmark settings (may include dgp_id/rep_id).
#' @param error_msg Error message string.
#' @returns Tibble with the same columns as successful benchmark output.
make_failure_metrics_row <- function(row, error_msg) {
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
    dgp_id = row$dgp_id %||% NA_integer_,
    rep_id = row$rep_id %||% NA_integer_,
    seed = row$seed %||% NA_integer_,
    n = row$n %||% NA_integer_,
    nxgrid = row$nxgrid %||% NA_integer_,
    nygrid = row$nygrid %||% NA_integer_,
    snr = row$snr %||% NA_real_,
    wiggliness = row$wiggliness %||% NA_real_,
    corr_type = row$corr_type %||% NA_character_,
    corr_param = row$corr_param %||% NA_real_,
    hetero_type = row$hetero_type %||% NA_character_,
    hetero_param = row$hetero_param %||% NA_real_,
    error_dist = row$error_dist %||% NA_character_,
    fit_time_total = NA_real_,
    fit_time_step_only = NA_real_,
    pffr_pilot_time = NA_real_,
    converged = FALSE,
    error_msg = error_msg
  )
}

#' Run one benchmark setting (data + fits + metrics)
#'
#' This is the shared core used by both the full benchmark runner and
#' interactive exploration tools.
#'
#' @param row One-row tibble/list with benchmark settings. Must include at least
#'   `n`, `nxgrid`, `nygrid`, `snr`, `wiggliness`, `corr_type`, `corr_param`,
#'   `hetero_type`, `hetero_param`, `error_dist`, and `seed`. If `dgp_id` and
#'   `rep_id` are present they are included in returned metrics.
#' @param alpha Significance level.
#' @param methods Optional character vector of methods to run. When provided,
#'   it is intersected with `get_applicable_methods()` for the DGP.
#' @returns List with `sim`, `fits`, `timings`, `pffr_pilot_time`, and `results`
#'   (standardized metrics tibble). On failure, returns `results` containing
#'   a single failure row and an `error` element.
run_one_setting <- function(row, alpha = 0.10, methods = NULL) {
  # Ensure row is list-like with $ access
  if (inherits(row, "data.frame")) row <- as.list(row)

  terms <- if (!is.null(row$terms)) {
    if (is.list(row$terms)) row$terms[[1]] else row$terms
  } else {
    c("ff", "linear", "smooth", "concurrent")
  }

  tryCatch(
    {
      sim <- generate_benchmark_data(
        n = row$n,
        nxgrid = row$nxgrid,
        nygrid = row$nygrid,
        snr = row$snr,
        error_dist = row$error_dist,
        corr_type = row$corr_type,
        corr_param = row$corr_param,
        hetero_type = row$hetero_type,
        hetero_param = row$hetero_param,
        terms = terms,
        wiggliness = row$wiggliness,
        seed = row$seed
      )

      applicable <- get_applicable_methods(
        row$corr_type,
        row$hetero_type,
        row$error_dist
      )
      if (is.null(methods)) {
        methods <- applicable
      } else {
        methods <- intersect(unique(methods), applicable)
        if (!length(methods)) {
          stop(
            "No requested methods are applicable for this DGP.",
            call. = FALSE
          )
        }
      }

      # Fit
      fit_result <- fit_all_methods(sim, methods)
      fits <- fit_result$fits
      timings <- fit_result$timings

      # Compute metrics for each fit/method
      results <- map_dfr(names(fits), function(method_name) {
        sandwich_type <- switch(
          method_name,
          pffr_sandwich = "cluster",
          pffr_hc = "hc",
          FALSE
        )
        fit <- if (method_name %in% c("pffr_sandwich", "pffr_hc")) {
          fits$pffr
        } else {
          fits[[method_name]]
        }

        metrics <- compute_all_metrics(
          fit,
          sim,
          method_name,
          use_sandwich = sandwich_type,
          alpha = alpha,
          err_struct = sim$err_struct
        )
        if (nrow(metrics) == 0) return(metrics)

        fit_time <- timings[[method_name]]
        step_time <- timings[[paste0(method_name, "_only")]] %||% fit_time

        metrics |>
          mutate(
            dgp_id = row$dgp_id %||% NA_integer_,
            rep_id = row$rep_id %||% NA_integer_,
            seed = row$seed,
            n = row$n,
            nxgrid = row$nxgrid,
            nygrid = row$nygrid,
            snr = row$snr,
            wiggliness = row$wiggliness,
            corr_type = row$corr_type,
            corr_param = row$corr_param,
            hetero_type = row$hetero_type,
            hetero_param = row$hetero_param,
            error_dist = row$error_dist,
            fit_time_total = fit_time,
            fit_time_step_only = step_time,
            pffr_pilot_time = fit_result$pffr_pilot_time,
            converged = TRUE,
            error_msg = NA_character_
          )
      })

      list(
        sim = sim,
        fits = fits,
        timings = timings,
        pffr_pilot_time = fit_result$pffr_pilot_time,
        results = results,
        error = NULL
      )
    },
    error = function(e) {
      msg <- conditionMessage(e)
      list(
        sim = NULL,
        fits = list(),
        timings = list(),
        pffr_pilot_time = NA_real_,
        results = make_failure_metrics_row(row, msg),
        error = msg
      )
    }
  )
}

#' Run benchmark for one (dgp, rep) combination
#'
#' @param row One row of the benchmark grid.
#' @param alpha Significance level.
#' @returns Tibble with metrics for all methods and terms.
run_one_dgp_rep <- function(row, alpha = 0.10) {
  out <- run_one_setting(row, alpha = alpha)
  if (!is.null(out$error)) {
    warning(
      "Failure in dgp_id=",
      row$dgp_id,
      " rep=",
      row$rep_id,
      ": ",
      out$error
    )
  }
  out$results
}

#' Run full benchmark
#'
#' @param dgp_settings DGP settings tibble.
#' @param n_rep Number of replications per DGP.
#' @param seed Base RNG seed.
#' @param parallel Use parallel processing?
#' @param n_workers Number of workers for parallel.
#' @param output_dir Directory for incremental saves.
#' @param alpha Significance level.
#' @returns Combined results tibble.
run_benchmark <- function(
  dgp_settings = make_dgp_settings("tiny"),
  n_rep = 10L,
  seed = 2024L,
  parallel = FALSE,
  n_workers = parallel::detectCores() - 1,
  output_dir = "ci-benchmark/results",
  alpha = 0.10
) {
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Expand to full grid with replications
  grid <- dgp_settings |>
    tidyr::crossing(rep_id = seq_len(n_rep)) |>
    mutate(seed = seed + 1000L * dgp_id + rep_id)

  cat("Benchmark grid:", nrow(grid), "total (dgp x rep) combinations\n")
  cat("DGP settings:", nrow(dgp_settings), "\n")
  cat("Replications:", n_rep, "\n")
  cat("Output dir:", output_dir, "\n\n")

  # Run sequentially or in parallel
  if (
    parallel &&
      requireNamespace("furrr", quietly = TRUE) &&
      requireNamespace("future", quietly = TRUE)
  ) {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)

    future::plan(future::multisession, workers = n_workers)

    # Split to preserve list-columns (e.g., terms) as 1-row tibbles
    rows <- split(grid, seq_len(nrow(grid)))

    # Detect if running in dev mode (via devtools::load_all)
    # Capture absolute paths so workers don't depend on working directory
    is_dev_mode <- file.exists("DESCRIPTION")
    pkg_dir <- normalizePath(".")
    bench_utils_path <- normalizePath("ci-benchmark/benchmark-utils.R")

    results <- furrr::future_map_dfr(
      rows,
      function(row) {
        # Load package in worker using absolute paths
        if (is_dev_mode) {
          devtools::load_all(pkg_dir, quiet = TRUE)
        } else {
          library(refund)
        }
        source(bench_utils_path, local = TRUE)

        res <- tryCatch(
          run_one_dgp_rep(row, alpha),
          error = function(e) {
            warning(
              "Error in dgp_id=",
              row$dgp_id,
              " rep=",
              row$rep_id,
              ": ",
              e$message
            )
            NULL
          }
        )
        if (!is.null(res) && nrow(res) > 0) {
          print(res)
          # Incremental save
          save_path <- file.path(
            output_dir,
            sprintf("dgp%03d_rep%03d.rds", row$dgp_id, row$rep_id)
          )
          saveRDS(res, save_path)
        }
        res
      },
      .options = furrr::furrr_options(seed = TRUE),
      .progress = TRUE
    )
  } else {
    # Sequential with progress
    results <- list()

    for (i in seq_len(nrow(grid))) {
      row <- grid[i, ]
      cat(sprintf(
        "\r[%d/%d] dgp=%d rep=%d",
        i,
        nrow(grid),
        row$dgp_id,
        row$rep_id
      ))

      res <- tryCatch(
        run_one_dgp_rep(row, alpha),
        error = function(e) {
          warning(
            "\nError in dgp_id=",
            row$dgp_id,
            " rep=",
            row$rep_id,
            ": ",
            e$message
          )
          NULL
        }
      )

      if (!is.null(res) && nrow(res) > 0) {
        results[[i]] <- res
        print(res)
        # Incremental save
        save_path <- file.path(
          output_dir,
          sprintf("dgp%03d_rep%03d.rds", row$dgp_id, row$rep_id)
        )
        saveRDS(res, save_path)
      }
    }
    cat("\n")

    results <- bind_rows(results)
  }

  # Save combined results
  saveRDS(results, file.path(output_dir, "results_combined.rds"))
  cat("Saved combined results:", nrow(results), "rows\n")

  results
}

# Summarization ---------------------------------------------------------------

#' Summarize benchmark results
#'
#' @param results Results tibble from run_benchmark().
#' @returns Summary tibble.
summarize_benchmark <- function(results) {
  # Add converged column if missing (backward compatibility)
  if (!"converged" %in% names(results)) {
    results$converged <- TRUE
  }

  n_failures_total <- sum(!results$converged, na.rm = TRUE)
  if (n_failures_total > 0) {
    message(
      "Note: ",
      n_failures_total,
      " failed row(s) excluded from summary statistics."
    )
  }

  # Filter to converged rows for metric computation
  results |>
    dplyr::filter(.data$converged) |>
    group_by(
      method,
      term_type,
      dgp_id,
      n,
      nxgrid,
      nygrid,
      snr,
      wiggliness,
      corr_type,
      corr_param,
      hetero_type,
      hetero_param,
      error_dist
    ) |>
    summarize(
      mean_coverage = mean(coverage, na.rm = TRUE),
      se_coverage = sd(coverage, na.rm = TRUE) / sqrt(n()),
      mean_coverage_high_var = mean(coverage_high_var, na.rm = TRUE),
      mean_coverage_low_var = mean(coverage_low_var, na.rm = TRUE),
      mean_width = mean(mean_width, na.rm = TRUE),
      mean_rmse = mean(rmse, na.rm = TRUE),
      mean_bias = mean(bias, na.rm = TRUE),
      mean_mean_se = mean(mean_se, na.rm = TRUE),
      mean_z_mean = mean(z_mean, na.rm = TRUE),
      # SE calibration: mean z_sd across reps (should be ~1 if SEs are correct)
      mean_z_sd = mean(z_sd, na.rm = TRUE),
      mean_z_kurtosis = mean(z_kurtosis, na.rm = TRUE),
      mean_fit_time_total = mean(fit_time_total, na.rm = TRUE),
      mean_fit_time_step = mean(fit_time_step_only, na.rm = TRUE),
      mean_pffr_pilot_time = mean(pffr_pilot_time, na.rm = TRUE),
      n_converged = n(),
      .groups = "drop"
    )
}

#' Summarize timings
#'
#' @param results Results tibble.
#' @returns Timing summary tibble.
summarize_timings <- function(results) {
  results |>
    distinct(dgp_id, rep_id, method, .keep_all = TRUE) |>
    group_by(
      method,
      dgp_id,
      n,
      nxgrid,
      nygrid,
      snr,
      wiggliness,
      corr_type,
      corr_param,
      hetero_type,
      hetero_param,
      error_dist
    ) |>
    summarize(
      mean_fit_time_total = mean(fit_time_total, na.rm = TRUE),
      sd_fit_time_total = sd(fit_time_total, na.rm = TRUE),
      mean_fit_time_step = mean(fit_time_step_only, na.rm = TRUE),
      mean_pffr_pilot_time = mean(pffr_pilot_time, na.rm = TRUE),
      n_rep = n(),
      .groups = "drop"
    )
}

# Main Execution (when sourced) -----------------------------------------------

if (FALSE) {
  # Production run (Round 2)
  results <- run_benchmark(
    dgp_settings = make_dgp_settings("full"),
    n_rep = 200,
    seed = 2026,
    parallel = TRUE,
    n_workers = 3
  )
}

if (sys.nframe() == 0) {
  # Running as script
  cat("CI Coverage Benchmark\n")
  cat("=====================\n\n")

  results <- run_benchmark(
    dgp_settings = make_dgp_settings("tiny"),
    n_rep = 2,
    seed = 2024,
    parallel = TRUE
  )

  summary_df <- summarize_benchmark(results)

  cat("\nCoverage summary:\n")
  print(
    summary_df |>
      select(method, term_type, mean_coverage, se_coverage, n_converged)
  )

  cat("\nTiming summary:\n")
  print(summarize_timings(results))

  cat("\nDone.\n")
}
